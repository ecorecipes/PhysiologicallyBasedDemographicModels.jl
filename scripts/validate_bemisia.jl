#!/usr/bin/env julia
# Validation script for Bemisia tabaci MEAM1 (B biotype) whitefly PBDM
# Parameters from Bonato et al. (2007) and Muñiz & Nombela (2001).
#
# Generates 5 figures in scripts/figures/bemisia/:
#   1. devrate_curves.png     — Brière development rate curves + literature data
#   2. mortality_curves.png   — U-shaped temperature-dependent mortality
#   3. sim_constant_25C.png   — Cohort simulation at constant 25°C
#   4. sim_rome_london.png    — Seasonal Rome vs London comparison
#   5. warming_comparison.png — Rome baseline / +1°C / +2°C invasion risk

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using PhysiologicallyBasedDemographicModels
using CairoMakie

figdir = joinpath(@__DIR__, "figures", "bemisia")
mkpath(figdir)

# ══════════════════════════════════════════════════════════════════════
# 1. Brière development rate models for B. tabaci MEAM1
#    r(T) = a · T · (T − T_lower) · √(T_upper − T)
#
#    Literature thermal thresholds:
#      Lower: 11.3°C (Muñiz & Nombela 2001)
#      Upper: 33–35°C (Bonato et al. 2007)
#      Optimal: 26–27°C
#      Total DD egg–adult: 333–385 DD
# ══════════════════════════════════════════════════════════════════════

const T_LOWER = 11.3   # common lower threshold (°C)

# Stage-specific Brière parameters calibrated to literature durations
egg_dev   = BriereDevelopmentRate(1.40e-4, T_LOWER, 35.5)
nymph_dev = BriereDevelopmentRate(5.00e-5, T_LOWER, 35.0)
adult_dev = BriereDevelopmentRate(3.50e-5, T_LOWER, 36.0)

println("Bemisia tabaci MEAM1 — Development rates (1/day):")
println("T(°C) |    Egg   |  Nymph   |  Adult")
println("-"^50)
for T in [15.0, 20.0, 25.0, 27.0, 30.0, 33.0]
    re = development_rate(egg_dev, T)
    rn = development_rate(nymph_dev, T)
    ra = development_rate(adult_dev, T)
    println("  $(lpad(T, 4))  | $(lpad(round(re, digits=5), 7)) | ",
            "$(lpad(round(rn, digits=5), 7)) | ",
            "$(lpad(round(ra, digits=5), 7))")
end
println()

# ══════════════════════════════════════════════════════════════════════
# 2. Temperature-dependent mortality — quadratic U-shaped
#    μ(T) = a·T² − b·T + c   (≥ 0)
#    High mortality below ~12°C and above ~33°C
# ══════════════════════════════════════════════════════════════════════

# Egg:   min ≈ 0.014 at ~26°C, high at 10°C (0.18) and 34°C (0.12)
const EGG_MORT_A  = 0.00110;  const EGG_MORT_B  = 0.0513;  const EGG_MORT_C  = 0.609
# Nymph: min ≈ 0.008 at ~26°C
const NYM_MORT_A  = 0.00087;  const NYM_MORT_B  = 0.0411;  const NYM_MORT_C  = 0.488
# Adult: min ≈ 0.015 at ~26°C
const ADU_MORT_A  = 0.00100;  const ADU_MORT_B  = 0.0442;  const ADU_MORT_C  = 0.486

egg_mortality(T)   = max(0.0, EGG_MORT_A * T^2 - EGG_MORT_B * T + EGG_MORT_C)
nymph_mortality(T) = max(0.0, NYM_MORT_A * T^2 - NYM_MORT_B * T + NYM_MORT_C)
adult_mortality(T) = max(0.0, ADU_MORT_A * T^2 - ADU_MORT_B * T + ADU_MORT_C)

println("Mortality rates (per day) at selected temperatures:")
println("T(°C) |   Egg    |  Nymph   |  Adult")
println("-"^50)
for T in [10.0, 15.0, 20.0, 25.0, 30.0, 35.0]
    println("  $(lpad(T, 4))  | $(lpad(round(egg_mortality(T), digits=4), 7)) | ",
            "$(lpad(round(nymph_mortality(T), digits=4), 7)) | ",
            "$(lpad(round(adult_mortality(T), digits=4), 7))")
end

# Optimal temperatures (vertex of parabola)
T_opt_egg   = EGG_MORT_B / (2 * EGG_MORT_A)
T_opt_nymph = NYM_MORT_B / (2 * NYM_MORT_A)
T_opt_adult = ADU_MORT_B / (2 * ADU_MORT_A)
println("\nMortality minima: egg=$(round(T_opt_egg, digits=1))°C, ",
        "nymph=$(round(T_opt_nymph, digits=1))°C, ",
        "adult=$(round(T_opt_adult, digits=1))°C")

# ══════════════════════════════════════════════════════════════════════
# 3. Fecundity — Brière-like temperature dependence
#    Peak ~2 eggs/female/day at ~27°C
# ══════════════════════════════════════════════════════════════════════

const FECUND_A    = 0.002014
const FECUND_TINF = 14.0
const FECUND_TSUP = 35.0

function fecundity(T::Real)
    (T <= FECUND_TINF || T >= FECUND_TSUP) && return 0.0
    return FECUND_A * T * (T - FECUND_TINF) * sqrt(FECUND_TSUP - T)
end

println("\nFecundity (eggs/female/day):")
for T in [20.0, 25.0, 27.0, 30.0]
    println("  T=$(T)°C: $(round(fecundity(T), digits=2))")
end

# ══════════════════════════════════════════════════════════════════════
# Figure 1 — Development rate curves with published data points
# ══════════════════════════════════════════════════════════════════════

temps = 0.0:0.5:40.0
tv = collect(temps)

# Published data points (Bonato et al. 2007; Muñiz & Nombela 2001)
egg_data_T   = [20.0, 25.0, 30.0]
egg_data_r   = [0.08, 0.14, 0.20]
nymph_data_T = [20.0, 25.0, 30.0]
nymph_data_r = [0.03, 0.06, 0.08]

fig1 = Figure(size=(900, 600))
ax1 = Axis(fig1[1,1],
    xlabel="Temperature (°C)",
    ylabel="Development rate (1/day)",
    title="B. tabaci MEAM1 — Brière development rate curves\n(Bonato et al. 2007; Muñiz & Nombela 2001)")

lines!(ax1, tv, [development_rate(egg_dev, T) for T in temps],
    label="Egg (Brière)", linewidth=2.5, color=:goldenrod)
lines!(ax1, tv, [development_rate(nymph_dev, T) for T in temps],
    label="Nymph (Brière)", linewidth=2.5, color=:forestgreen)
lines!(ax1, tv, [development_rate(adult_dev, T) for T in temps],
    label="Adult (Brière)", linewidth=2.5, color=:navy)

scatter!(ax1, egg_data_T, egg_data_r,
    marker=:circle, markersize=12, color=:goldenrod, strokewidth=1.5,
    strokecolor=:black, label="Egg data (Bonato)")
scatter!(ax1, nymph_data_T, nymph_data_r,
    marker=:utriangle, markersize=12, color=:forestgreen, strokewidth=1.5,
    strokecolor=:black, label="Nymph data (Muñiz)")

xlims!(ax1, 0, 40)
ylims!(ax1, 0, nothing)
axislegend(ax1, position=:lt)

save(joinpath(figdir, "devrate_curves.png"), fig1, px_per_unit=2)
println("\nSaved devrate_curves.png")
println("  Egg peak: $(round(maximum(development_rate.(Ref(egg_dev), tv)), digits=4))")
println("  Nymph peak: $(round(maximum(development_rate.(Ref(nymph_dev), tv)), digits=4))")
println("  Adult peak: $(round(maximum(development_rate.(Ref(adult_dev), tv)), digits=4))")

# ══════════════════════════════════════════════════════════════════════
# Figure 2 — U-shaped mortality curves
# ══════════════════════════════════════════════════════════════════════

mort_temps = 0.0:0.5:42.0
mtv = collect(mort_temps)

fig2 = Figure(size=(900, 600))
ax2 = Axis(fig2[1,1],
    xlabel="Temperature (°C)",
    ylabel="Mortality rate (per day)",
    title="B. tabaci MEAM1 — Temperature-dependent mortality\n(quadratic U-shaped model)")

lines!(ax2, mtv, [egg_mortality(T) for T in mort_temps],
    label="Egg", linewidth=2.5, color=:goldenrod)
lines!(ax2, mtv, [nymph_mortality(T) for T in mort_temps],
    label="Nymph", linewidth=2.5, color=:forestgreen)
lines!(ax2, mtv, [adult_mortality(T) for T in mort_temps],
    label="Adult", linewidth=2.5, color=:navy)

# Mark the optimal zone
vspan!(ax2, 22.0, 30.0, color=(:green, 0.08))
text!(ax2, 26.0, 0.005,
    text="Optimal zone", align=(:center, :bottom), fontsize=10, color=:gray40)

xlims!(ax2, 0, 42)
ylims!(ax2, 0, nothing)
axislegend(ax2, position=:rt)

save(joinpath(figdir, "mortality_curves.png"), fig2, px_per_unit=2)
println("Saved mortality_curves.png")

# ══════════════════════════════════════════════════════════════════════
# Population builder — creates a fresh population for each simulation
# ══════════════════════════════════════════════════════════════════════

function build_bemisia(; N0_egg=0.0, N0_nymph=0.0, N0_adult=0.0,
                        egg_μ=0.02, nymph_μ=0.02, adult_μ=0.03)
    stages = [
        LifeStage(:egg,
            DistributedDelay(6, 1.0; W0=N0_egg / 6),
            egg_dev, egg_μ),
        LifeStage(:nymph,
            DistributedDelay(10, 1.0; W0=N0_nymph / 10),
            nymph_dev, nymph_μ),
        LifeStage(:adult,
            DistributedDelay(8, 1.0; W0=N0_adult / 8),
            adult_dev, adult_μ),
    ]
    Population(:bemisia_tabaci, stages)
end

# ══════════════════════════════════════════════════════════════════════
# Figure 3 — Cohort simulation at constant 25°C (180 days)
#   Start with 500 eggs, no reproduction (density-independent)
# ══════════════════════════════════════════════════════════════════════

println("\n── Cohort simulation: 500 eggs at constant 25°C, 180 days ──")

pop_cohort = build_bemisia(N0_egg=500.0, egg_μ=0.02, nymph_μ=0.02, adult_μ=0.03)
weather_25C = WeatherSeries([DailyWeather(25.0) for _ in 1:180])
prob_cohort = PBDMProblem(pop_cohort, weather_25C, (1, 180))
sol_cohort = solve(prob_cohort, DirectIteration())

println("  Return code: $(sol_cohort.retcode)")
for (i, sname) in enumerate([:egg, :nymph, :adult])
    traj = sol_cohort.stage_totals[i, 2:end]
    peak = maximum(traj)
    peak_day = argmax(traj)
    println("  $(sname): peak=$(round(peak, digits=1)) at day $(peak_day)")
end
pop_total_cohort = total_population(sol_cohort)
println("  Peak total: $(round(maximum(pop_total_cohort), digits=1))")

stage_names = ["Eggs", "Nymphs", "Adults"]
stage_colors = [:goldenrod, :forestgreen, :navy]

fig3 = Figure(size=(900, 600))
ax3a = Axis(fig3[1,1],
    xlabel="Day",
    ylabel="Population",
    title="B. tabaci cohort at constant 25°C — 500 initial eggs")

for (j, (sname, scol)) in enumerate(zip(stage_names, stage_colors))
    traj = sol_cohort.stage_totals[j, 2:end]
    lines!(ax3a, sol_cohort.t[2:end], traj, label=sname, color=scol, linewidth=2)
end
axislegend(ax3a, position=:rt)

ax3b = Axis(fig3[2,1],
    xlabel="Day",
    ylabel="Total population")
lines!(ax3b, sol_cohort.t, pop_total_cohort, color=:black, linewidth=2)

save(joinpath(figdir, "sim_constant_25C.png"), fig3, px_per_unit=2)
println("Saved sim_constant_25C.png")

# ══════════════════════════════════════════════════════════════════════
# Seasonal simulation helper (density-dependent with reproduction)
# ══════════════════════════════════════════════════════════════════════

const CARRYING_CAPACITY = 100000.0

function simulate_bemisia_seasonal(weather, tspan; N0_adult=100.0)
    pop = build_bemisia(N0_adult=N0_adult,
                        egg_μ=0.30, nymph_μ=0.40, adult_μ=0.20)

    reproduction_fn = (p, w, params, day) -> begin
        T = w.T_mean
        N_adults = delay_total(p.stages[3].delay)
        N_total = total_population(p)
        density_factor = max(0.0, 1.0 - N_total / CARRYING_CAPACITY)
        return max(0.0, fecundity(T) * N_adults * density_factor)
    end

    prob = PBDMProblem(DensityDependent(), pop, weather, tspan)
    sol = solve(prob, DirectIteration(); reproduction_fn=reproduction_fn)
    return sol
end

# ══════════════════════════════════════════════════════════════════════
# Figure 4 — Rome vs London seasonal comparison (365 days)
#   Rome:   mean 16°C, amplitude 8°C  → range 8–24°C
#   London: mean 11°C, amplitude 7°C  → range 4–18°C
# ══════════════════════════════════════════════════════════════════════

println("\n── Seasonal simulations: Rome vs London, 365 days ──")

weather_rome   = SinusoidalWeather(16.0, 8.0)
weather_london = SinusoidalWeather(11.0, 7.0)

sol_rome   = simulate_bemisia_seasonal(weather_rome,   (1, 365))
sol_london = simulate_bemisia_seasonal(weather_london, (1, 365))

pop_rome   = total_population(sol_rome)
pop_london = total_population(sol_london)

println("  Rome:   peak=$(round(maximum(pop_rome), digits=1)), " *
        "final=$(round(pop_rome[end], digits=1))")
println("  London: peak=$(round(maximum(pop_london), digits=1)), " *
        "final=$(round(pop_london[end], digits=1))")

# Temperature profiles for annotation
temp_rome   = [16.0 + 8.0 * sin(2π * (d - 200) / 365) for d in 1:365]
temp_london = [11.0 + 7.0 * sin(2π * (d - 200) / 365) for d in 1:365]

fig4 = Figure(size=(1200, 700))

# Rome panel
ax4a_t = Axis(fig4[1,1], ylabel="°C", title="Rome (mean 16°C, amp 8°C)")
lines!(ax4a_t, 1:365, temp_rome, color=:red, linewidth=1.2)
hlines!(ax4a_t, [T_LOWER], linestyle=:dash, color=:gray60, linewidth=0.8)
hidexdecorations!(ax4a_t, grid=false)

ax4a_p = Axis(fig4[2,1], xlabel="Day of year", ylabel="Total population")
lines!(ax4a_p, sol_rome.t, pop_rome, color=:firebrick, linewidth=2)

# London panel
ax4b_t = Axis(fig4[1,2], ylabel="°C", title="London (mean 11°C, amp 7°C)")
lines!(ax4b_t, 1:365, temp_london, color=:blue, linewidth=1.2)
hlines!(ax4b_t, [T_LOWER], linestyle=:dash, color=:gray60, linewidth=0.8)
hidexdecorations!(ax4b_t, grid=false)

ax4b_p = Axis(fig4[2,2], xlabel="Day of year", ylabel="Total population")
lines!(ax4b_p, sol_london.t, pop_london, color=:steelblue, linewidth=2)

# Link y-axes for comparison
linkyaxes!(ax4a_p, ax4b_p)
linkyaxes!(ax4a_t, ax4b_t)

Label(fig4[0, :], "B. tabaci seasonal dynamics — Rome vs London (100 initial adults)",
    fontsize=16, font=:bold)

save(joinpath(figdir, "sim_rome_london.png"), fig4, px_per_unit=2)
println("Saved sim_rome_london.png")

# ══════════════════════════════════════════════════════════════════════
# Figure 5 — Warming comparison: Rome at baseline, +1°C, +2°C
# ══════════════════════════════════════════════════════════════════════

println("\n── Warming comparison: Rome baseline / +1°C / +2°C ──")

# Use uncapped growth (no K) to reveal temperature-driven divergence
function simulate_warming(weather, tspan; N0_adult=10.0)
    pop = build_bemisia(N0_adult=N0_adult,
                        egg_μ=0.80, nymph_μ=1.00, adult_μ=0.40)
    reproduction_fn = (p, w, params, day) -> begin
        T = w.T_mean
        N_adults = delay_total(p.stages[3].delay)
        return max(0.0, fecundity(T) * N_adults)
    end
    prob = PBDMProblem(DensityDependent(), pop, weather, tspan)
    sol = solve(prob, DirectIteration(); reproduction_fn=reproduction_fn)
    return sol
end

weather_rome_0 = SinusoidalWeather(16.0, 8.0)
weather_rome_1 = SinusoidalWeather(17.0, 8.0)
weather_rome_2 = SinusoidalWeather(18.0, 8.0)

sol_w0 = simulate_warming(weather_rome_0, (1, 365))
sol_w1 = simulate_warming(weather_rome_1, (1, 365))
sol_w2 = simulate_warming(weather_rome_2, (1, 365))

pop_w0 = total_population(sol_w0)
pop_w1 = total_population(sol_w1)
pop_w2 = total_population(sol_w2)

println("  Baseline: peak=$(round(maximum(pop_w0), digits=1))")
println("  +1°C:     peak=$(round(maximum(pop_w1), digits=1))")
println("  +2°C:     peak=$(round(maximum(pop_w2), digits=1))")

# Temperature profiles
temp_w0 = [16.0 + 8.0 * sin(2π * (d - 200) / 365) for d in 1:365]
temp_w1 = [17.0 + 8.0 * sin(2π * (d - 200) / 365) for d in 1:365]
temp_w2 = [18.0 + 8.0 * sin(2π * (d - 200) / 365) for d in 1:365]

fig5 = Figure(size=(1000, 700))

# Top: temperature profiles
ax5t = Axis(fig5[1,1], ylabel="Temperature (°C)",
    title="B. tabaci invasion risk under warming — Rome climate")
lines!(ax5t, 1:365, temp_w0, color=:blue,   linewidth=1.2, label="Baseline (16°C)")
lines!(ax5t, 1:365, temp_w1, color=:orange,  linewidth=1.2, label="+1°C (17°C)")
lines!(ax5t, 1:365, temp_w2, color=:red,     linewidth=1.2, label="+2°C (18°C)")
hlines!(ax5t, [T_LOWER], linestyle=:dash, color=:gray60, linewidth=0.8)
axislegend(ax5t, position=:rb)
hidexdecorations!(ax5t, grid=false)

# Bottom: total population (log scale to show divergent growth rates)
ax5b = Axis(fig5[2,1], xlabel="Day of year", ylabel="Total population (log₁₀)",
    yscale=log10)
lines!(ax5b, sol_w0.t, max.(pop_w0, 1.0), color=:blue,   linewidth=2.5, label="Baseline")
lines!(ax5b, sol_w1.t, max.(pop_w1, 1.0), color=:orange,  linewidth=2.5, label="+1°C")
lines!(ax5b, sol_w2.t, max.(pop_w2, 1.0), color=:red,     linewidth=2.5, label="+2°C")
axislegend(ax5b, position=:lt)

save(joinpath(figdir, "warming_comparison.png"), fig5, px_per_unit=2)
println("Saved warming_comparison.png")

# ══════════════════════════════════════════════════════════════════════
# Summary
# ══════════════════════════════════════════════════════════════════════

println("\n══ All Bemisia validation figures saved to: $(figdir) ══")
println("Files:")
for f in sort(readdir(figdir))
    fpath = joinpath(figdir, f)
    sz = round(filesize(fpath) / 1024, digits=1)
    println("  $f  ($(sz) KB)")
end
