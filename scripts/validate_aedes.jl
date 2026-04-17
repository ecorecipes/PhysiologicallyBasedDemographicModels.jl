#!/usr/bin/env julia
# Validation script for Ae. albopictus thermal biology
# matching the vignette (Pasquali et al. 2020 parameters).

using PhysiologicallyBasedDemographicModels
using CairoMakie

figdir = joinpath(@__DIR__, "figures", "aedes")
mkpath(figdir)

# ══════════════════════════════════════════════════════════════════════
# 1. Development rate parameters (Pasquali et al. 2020, Table 1)
# ══════════════════════════════════════════════════════════════════════

# Egg development: cubic r(T) = a·T²·(T_sup − T)  (Eq. 9)
egg_devrate(T) = (T <= 0 || T >= 37.3253) ? 0.0 : 2.948e-5 * T^2 * (37.3253 - T)

# Larval, pupal, immature-adult: Brière function (Eq. 10)
larva_dev  = BriereDevelopmentRate(8.604e-5, 8.2934, 36.0729)
pupa_dev   = BriereDevelopmentRate(3.102e-4, 11.9433, 40.0)
immadult_dev = BriereDevelopmentRate(1.812e-4, 7.7804, 35.2937)

# Reproductive adult: constant (Section 2.1.3)
const ADULT_DEV_RATE = 0.015

# Print development-rate table
println("Development rates (1/day) at selected temperatures:")
println("T(°C) |   Egg    |  Larva   |  Pupa    | Imm.Adult")
println("-"^60)
for T in [10.0, 15.0, 20.0, 25.0, 30.0, 35.0]
    re = egg_devrate(T)
    rl = development_rate(larva_dev, T)
    rp = development_rate(pupa_dev, T)
    ra = development_rate(immadult_dev, T)
    println("  $(lpad(T, 4))  | $(lpad(round(re, digits=5), 7)) | ",
            "$(lpad(round(rl, digits=5), 7)) | ",
            "$(lpad(round(rp, digits=5), 7)) | ",
            "$(lpad(round(ra, digits=5), 7))")
end

# ══════════════════════════════════════════════════════════════════════
# 2. Temperature-dependent mortality (Tables 2–3, Pasquali et al. 2020)
# ══════════════════════════════════════════════════════════════════════

# Proportional mortality m(T) = a·T² + b·T + c (Eq. 11, Table 2)
const EGG_MORT_A   =  0.002869;  const EGG_MORT_B   = -0.1417;  const EGG_MORT_C   = 2.1673
const LARVA_MORT_A =  0.002793;  const LARVA_MORT_B = -0.1255;  const LARVA_MORT_C = 1.5768
const PUPA_MORT_A  =  0.003289;  const PUPA_MORT_B  = -0.1437;  const PUPA_MORT_C  = 1.6197

# Extension coefficients (Eq. 12, Table 3)
const EGG_ALOW    = 0.05;   const EGG_AHIGH    = 0.018
const EGG_TLOW    = 12.0;   const EGG_THIGH    = 30.5
const LARVA_ALOW  = 0.2;    const LARVA_AHIGH  = 0.15
const LARVA_TLOW  = 17.0;   const LARVA_THIGH  = 30.5
const PUPA_ALOW   = 0.1;    const PUPA_AHIGH   = 0.2
const PUPA_TLOW   = 17.0;   const PUPA_THIGH   = 37.5

const ADULT_MORTALITY = 0.067   # Section 2.1.3

proportional_mortality(T, a, b, c) = clamp(a * T^2 + b * T + c, 0.0, 0.9)

function _central_mortality(T, dev_fn, a, b, c)
    m = proportional_mortality(T, a, b, c)
    v = isa(dev_fn, Function) ? dev_fn(T) : development_rate(dev_fn, T)
    (m <= 0.0 || m >= 1.0 || v <= 0.0) && return 0.0
    return -v * log(1.0 - m)
end

function stage_mortality_rate(T, dev_fn, ma, mb, mc, alow, ahigh, Tlow, Thigh)
    if T < Tlow
        μ_ref = _central_mortality(Tlow, dev_fn, ma, mb, mc)
        return max(0.0, μ_ref + alow * (Tlow - T))
    elseif T > Thigh
        μ_ref = _central_mortality(Thigh, dev_fn, ma, mb, mc)
        return max(0.0, μ_ref + ahigh * (T - Thigh))
    else
        return _central_mortality(T, dev_fn, ma, mb, mc)
    end
end

egg_mortality(T)   = stage_mortality_rate(T, egg_devrate,
    EGG_MORT_A, EGG_MORT_B, EGG_MORT_C, EGG_ALOW, EGG_AHIGH, EGG_TLOW, EGG_THIGH)
larva_mortality(T) = stage_mortality_rate(T, larva_dev,
    LARVA_MORT_A, LARVA_MORT_B, LARVA_MORT_C, LARVA_ALOW, LARVA_AHIGH, LARVA_TLOW, LARVA_THIGH)
pupa_mortality(T)  = stage_mortality_rate(T, pupa_dev,
    PUPA_MORT_A, PUPA_MORT_B, PUPA_MORT_C, PUPA_ALOW, PUPA_AHIGH, PUPA_TLOW, PUPA_THIGH)

println("\nMortality rates (1/day) at selected temperatures:")
println("T(°C) |   Egg    |  Larva   |  Pupa    | Adults")
println("-"^60)
for T in [10.0, 15.0, 20.0, 25.0, 30.0, 35.0]
    println("  $(lpad(T, 4))  | $(lpad(round(egg_mortality(T), digits=4), 7)) | ",
            "$(lpad(round(larva_mortality(T), digits=4), 7)) | ",
            "$(lpad(round(pupa_mortality(T), digits=4), 7)) | ",
            "$(lpad(round(ADULT_MORTALITY, digits=4), 7))")
end

# ══════════════════════════════════════════════════════════════════════
# 3. Figure 1 — All 4 stage development rates on one panel
# ══════════════════════════════════════════════════════════════════════

temps = 0.0:0.5:42.0
tv = collect(temps)

fig1 = Figure(size=(800, 500))
ax1 = Axis(fig1[1,1],
    xlabel="Temperature (°C)",
    ylabel="Development rate (1/day)",
    title="Ae. albopictus stage-specific development rates\n(Pasquali et al. 2020)")

lines!(ax1, tv, [egg_devrate(T) for T in temps],
    label="Eggs (cubic)", linewidth=2)
lines!(ax1, tv, [development_rate(larva_dev, T) for T in temps],
    label="Larvae (Brière)", linewidth=2)
lines!(ax1, tv, [development_rate(pupa_dev, T) for T in temps],
    label="Pupae (Brière)", linewidth=2)
lines!(ax1, tv, [development_rate(immadult_dev, T) for T in temps],
    label="Immature adults (Brière)", linewidth=2)

axislegend(ax1, position=:lt)
save(joinpath(figdir, "devrates.png"), fig1, px_per_unit=2)
println("\nSaved devrates.png")
println("  Egg peak: ", round(maximum(egg_devrate.(tv)), digits=4),
        " at T≈", tv[argmax(egg_devrate.(tv))], "°C")
println("  Larva peak: ", round(maximum([development_rate(larva_dev, T) for T in tv]), digits=4))
println("  Pupa peak: ", round(maximum([development_rate(pupa_dev, T) for T in tv]), digits=4))
println("  Imm.adult peak: ", round(maximum([development_rate(immadult_dev, T) for T in tv]), digits=4))

# ══════════════════════════════════════════════════════════════════════
# 4. Figure 2 — Temperature-dependent mortality rates
# ══════════════════════════════════════════════════════════════════════

mort_temps = 0.0:0.5:42.0
mtv = collect(mort_temps)

fig2 = Figure(size=(800, 500))
ax2 = Axis(fig2[1,1],
    xlabel="Temperature (°C)",
    ylabel="Mortality rate (1/day)",
    title="Ae. albopictus stage-specific mortality rates\n(Pasquali et al. 2020, Tables 2–3)")

lines!(ax2, mtv, [egg_mortality(T) for T in mort_temps],
    label="Eggs", linewidth=2)
lines!(ax2, mtv, [larva_mortality(T) for T in mort_temps],
    label="Larvae", linewidth=2)
lines!(ax2, mtv, [pupa_mortality(T) for T in mort_temps],
    label="Pupae", linewidth=2)
hlines!(ax2, [ADULT_MORTALITY], linestyle=:dash, color=:gray,
    label="Adults (constant)")

axislegend(ax2, position=:lt)
save(joinpath(figdir, "mortality.png"), fig2, px_per_unit=2)
println("\nSaved mortality.png")

# ══════════════════════════════════════════════════════════════════════
# 5. Custom types needed for LifeStage construction
# ══════════════════════════════════════════════════════════════════════

# Constant development rate for reproductive adults
struct ConstantDevelopmentRate{T<:Real} <: AbstractDevelopmentRate
    rate::T
end

function PhysiologicallyBasedDemographicModels.development_rate(
        m::ConstantDevelopmentRate, T::Real)
    return m.rate
end

function PhysiologicallyBasedDemographicModels.degree_days(
        m::ConstantDevelopmentRate, T::Real)
    return m.rate
end

# Brière approximation for egg cubic (similar peak location and magnitude)
# Cubic peak: T≈24.9°C, value≈0.227
# Brière fitted: a·T·(T − T_lower)·√(T_upper − T) ≈ 0.227 at T=24.9
const EGG_BRIERE_A     = 1.06e-4
const EGG_BRIERE_TINF  = 0.5
const EGG_BRIERE_TSUP  = 37.3253
egg_dev_approx = BriereDevelopmentRate(EGG_BRIERE_A, EGG_BRIERE_TINF, EGG_BRIERE_TSUP)

println("\nEgg dev rate comparison (cubic vs Brière approx):")
for T in [15.0, 20.0, 25.0, 30.0]
    rc = egg_devrate(T)
    rb = development_rate(egg_dev_approx, T)
    println("  T=$(T)°C: cubic=$(round(rc, digits=4)), brière=$(round(rb, digits=4))")
end

# ══════════════════════════════════════════════════════════════════════
# 6. Build Ae. albopictus population model
# ══════════════════════════════════════════════════════════════════════

# Fecundity: Brière function (Eq. 14, Table 5)
const FECUND_A    = 0.0032
const FECUND_TINF = 16.24
const FECUND_TSUP = 35.02

function fecundity(T::Real)
    (T <= FECUND_TINF || T >= FECUND_TSUP) && return 0.0
    return FECUND_A * T * (T - FECUND_TINF) * sqrt(FECUND_TSUP - T)
end

adult_dev = ConstantDevelopmentRate(ADULT_DEV_RATE)

function build_aedes_population(; N0_adults=100.0, ref_temp=25.0)
    egg_μ   = egg_mortality(ref_temp)
    larva_μ = larva_mortality(ref_temp)
    pupa_μ  = pupa_mortality(ref_temp)

    stages = [
        LifeStage(:eggs,
            DistributedDelay(15, 45.0; W0=0.0),
            egg_dev_approx, egg_μ),
        LifeStage(:larvae,
            DistributedDelay(20, 120.0; W0=0.0),
            larva_dev, larva_μ),
        LifeStage(:pupae,
            DistributedDelay(10, 30.0; W0=0.0),
            pupa_dev, pupa_μ),
        LifeStage(:immature_adults,
            DistributedDelay(10, 50.0; W0=0.0),
            immadult_dev, ADULT_MORTALITY),
        LifeStage(:adults,
            DistributedDelay(25, 67.0; W0=N0_adults),
            adult_dev, ADULT_MORTALITY),
    ]
    Population(:aedes_albopictus, stages)
end

# ══════════════════════════════════════════════════════════════════════
# 7. Simulation helper
# ══════════════════════════════════════════════════════════════════════

function simulate_aedes(weather, tspan; N0_adults=100.0)
    pop = build_aedes_population(N0_adults=N0_adults)

    reproduction_fn = (p, w, params, day) -> begin
        T = w.T_mean
        N_adults = delay_total(p.stages[5].delay)
        return max(0.0, fecundity(T) * N_adults)
    end

    prob = PBDMProblem(DensityDependent(), pop, weather, tspan)
    sol = solve(prob, DirectIteration(); reproduction_fn=reproduction_fn)
    return sol
end

# ══════════════════════════════════════════════════════════════════════
# 8. Figure 3 — Constant 25°C simulation (365 days)
# ══════════════════════════════════════════════════════════════════════

println("\n── Constant 25°C simulation (365 days) ──")

const_weather = WeatherSeries([25.0 for _ in 1:365])
sol_const = simulate_aedes(const_weather, (1, 365))

pop_total = total_population(sol_const)
println("  Return code: ", sol_const.retcode)
println("  Peak total population: ", round(maximum(pop_total), digits=1))
println("  Final total population: ", round(pop_total[end], digits=1))

stage_names = ["Eggs", "Larvae", "Pupae", "Imm. adults", "Adults"]
stage_colors = [:goldenrod, :forestgreen, :darkorchid, :coral, :navy]

fig3 = Figure(size=(900, 600))
ax3a = Axis(fig3[1,1],
    xlabel="Day",
    ylabel="Population",
    title="Ae. albopictus at constant 25°C — stage dynamics")

for (j, (sname, scol)) in enumerate(zip(stage_names, stage_colors))
    traj = stage_trajectory(sol_const, j)
    lines!(ax3a, sol_const.t, traj, label=sname, color=scol, linewidth=1.5)
end
axislegend(ax3a, position=:rt)

ax3b = Axis(fig3[2,1],
    xlabel="Day",
    ylabel="Total population",
    title="Total population")
lines!(ax3b, sol_const.t, pop_total, color=:black, linewidth=2)

save(joinpath(figdir, "sim_constant_25C.png"), fig3, px_per_unit=2)
println("Saved sim_constant_25C.png")

# ══════════════════════════════════════════════════════════════════════
# 9. Figure 4 — Seasonal simulation (sinusoidal, min=5°C, max=30°C)
# ══════════════════════════════════════════════════════════════════════

println("\n── Seasonal simulation (T_min=5°C, T_max=30°C, 365 days) ──")

# T(d) = T_mean + amplitude·sin(2π(d−phase)/365)
# min = T_mean − amplitude = 5  → T_mean = 17.5, amplitude = 12.5
# phase = 200 puts peak in late July (day ~200)
seasonal_weather = SinusoidalWeather(17.5, 12.5; phase=200.0)
sol_seasonal = simulate_aedes(seasonal_weather, (1, 365))

pop_seasonal = total_population(sol_seasonal)
println("  Return code: ", sol_seasonal.retcode)
println("  Peak total population: ", round(maximum(pop_seasonal), digits=1))
println("  Final total population: ", round(pop_seasonal[end], digits=1))

# Temperature profile for annotation
temp_profile = [17.5 + 12.5 * sin(2π * (d - 200) / 365) for d in 1:365]
println("  Temperature range: ", round(minimum(temp_profile), digits=1),
        "°C to ", round(maximum(temp_profile), digits=1), "°C")

fig4 = Figure(size=(900, 750))

# Top: temperature forcing
ax4t = Axis(fig4[1,1],
    ylabel="Temperature (°C)",
    title="Seasonal forcing: sinusoidal T (5–30°C)")
lines!(ax4t, 1:365, temp_profile, color=:red, linewidth=1.5)
hlines!(ax4t, [17.5], linestyle=:dash, color=:gray60)
hidexdecorations!(ax4t, grid=false)

# Middle: per-stage dynamics
ax4m = Axis(fig4[2,1],
    ylabel="Population",
    title="Ae. albopictus stage dynamics")

for (j, (sname, scol)) in enumerate(zip(stage_names, stage_colors))
    traj = stage_trajectory(sol_seasonal, j)
    lines!(ax4m, sol_seasonal.t, traj, label=sname, color=scol, linewidth=1.5)
end
axislegend(ax4m, position=:rt)
hidexdecorations!(ax4m, grid=false)

# Bottom: total population
ax4b = Axis(fig4[3,1],
    xlabel="Day of year",
    ylabel="Total population",
    title="Total population")
lines!(ax4b, sol_seasonal.t, pop_seasonal, color=:black, linewidth=2)

save(joinpath(figdir, "sim_seasonal.png"), fig4, px_per_unit=2)
println("Saved sim_seasonal.png")

# ══════════════════════════════════════════════════════════════════════
# Summary
# ══════════════════════════════════════════════════════════════════════

println("\n══ All Aedes validation figures saved to: ", figdir, " ══")
println("Files:")
for f in readdir(figdir)
    fpath = joinpath(figdir, f)
    sz = round(filesize(fpath) / 1024, digits=1)
    println("  $f  ($(sz) KB)")
end
