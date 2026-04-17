#!/usr/bin/env julia
# Validation script for the olive–climate vignette (11_olive_climate.qmd).
#
# Compares PBDM outputs against published parameters from:
#   Ponti et al. (2014) PNAS 111(15):5598–5603
#   Ponti, Cossu & Gutierrez (2009) Glob. Change Biol. 15:2874–2884
#   Rochat & Gutierrez (2001) J. Anim. Ecol. 70:476–490
#
# Generates five PNG figures in scripts/figures/olive_climate/.

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using PhysiologicallyBasedDemographicModels
using CairoMakie
using Statistics

figdir = joinpath(@__DIR__, "figures", "olive_climate")
mkpath(figdir)

# ============================================================
# 1. Shared parameters
# ============================================================

# Olive tree phenology (degree-day model, base 9.1 °C)
const OLIVE_T_BASE  = 9.1     # °C — budburst/bloom accumulation threshold
const OLIVE_T_UPPER = 40.0    # °C — upper developmental ceiling
const DD_BLOOM      = 400.0   # DD>9.1 from Jan 1 to flowering (Ponti 2009)
const DD_FRUIT_SET  = 500.0   # DD to fruit set
const DD_PIT_HARDEN = 900.0   # DD to pit hardening (onset of fly susceptibility)
const DD_OIL_ACCUM  = 1400.0  # DD to oil accumulation phase
const DD_MATURITY   = 1800.0  # DD to fruit maturity / harvest

# Olive fly (Bactrocera oleae) parameters — Gutierrez et al. (2009)
const FLY_T_BASE    = 10.0    # °C — base development temperature
const FLY_T_UPPER   = 35.0    # °C — upper lethal / activity threshold
const FLY_BRIERE_A  = 2.0e-4  # Brière rate coefficient (for plotting)
# Stage-specific thresholds from validate_olive_fly.jl
const FLY_T_EGG     = 9.2
const FLY_T_LARVA   = 13.9
const FLY_T_PUPA    = 12.4
const FLY_T_ADULT   = 10.0
# Stage DD requirements (linear model)
const FLY_DD_EGG    = 65.0
const FLY_DD_LARVA  = 220.0
const FLY_DD_PUPA   = 220.0
const FLY_DD_ADULT  = 600.0

# Ponti et al. (2009) yield regression (Eqn 3, Sardinia):
#   Yield(g DM/tree) = 3363.7 + 1.6·DDa − 84.2·DDb − 0.9·mm   R²=0.47
olive_yield_model = WeatherYieldModel(1.6, -0.9, -0.002, 3363.7)

# Ponti et al. (2009) fly regression (Eqn 4):
#   log₁₀(pupae) = 0.00730·DDa − 0.034·DDb + 0.0087·Blmday − 0.00036·mm
# Paper reports avg effects dda→1.8, bloom→1.2, ddb→−0.23, mm→−0.21
# giving average log₁₀ ≈ 2.56 (i.e., ~360 pupae/tree at typical Sardinia DDa ≈ 250).
# The coefficient for DDa in the OCR is 10× too large (0.00730 vs 0.000730);
# we rescale to match the paper's stated average contributions.
const FLY_β_DDa   =  0.000730   # rescaled: avg effect ≈ 1.8 at DDa=2466
const FLY_β_DDb   = -0.000383   # rescaled: avg effect ≈ −0.23 at DDb=600
const FLY_β_BLOOM =  0.0092     # avg effect ≈ 1.2 at bloom=130
const FLY_β_RAIN  = -0.000356   # avg effect ≈ −0.21 at mm=589

# Fly damage (calibrated: ~5000 pupae/tree ≈ 60% fruit attacked)
fly_damage = ExponentialDamageFunction(0.00018)

# Typical Mediterranean rainfall Apr–Sep (mm)
const RAIN_MM = 250.0

# Sinusoidal weather phase: peak at day ~191 (mid-July) for N. hemisphere
const PHASE = 100.0

# Published reference values for annotation
const PAPER_BLOOM_SHIFT = -9.4   # days per °C (Ponti et al. 2009)
const PAPER_FRUIT_ATTACK_RANGE = (38.0, 87.0)  # % fruit attacked, Sardinia obs.
const PAPER_YIELD_CHANGE_1_8 = 4.1   # % basin-wide yield change at +1.8 °C
const PAPER_FLY_CHANGE_1_8 = -8.0    # % basin-wide fly change at +1.8 °C

# Development rate models
olive_dev = LinearDevelopmentRate(OLIVE_T_BASE, OLIVE_T_UPPER)
fly_briere = BriereDevelopmentRate(FLY_BRIERE_A, FLY_T_BASE, FLY_T_UPPER)
fly_linear = LinearDevelopmentRate(FLY_T_BASE, FLY_T_UPPER)
# Stage-specific linear rates for PBDM simulation
egg_dev   = LinearDevelopmentRate(FLY_T_EGG,   FLY_T_UPPER)
larva_dev = LinearDevelopmentRate(FLY_T_LARVA, FLY_T_UPPER)
pupa_dev  = LinearDevelopmentRate(FLY_T_PUPA,  FLY_T_UPPER)
adult_dev = LinearDevelopmentRate(FLY_T_ADULT, FLY_T_UPPER)

# ============================================================
# Helper functions
# ============================================================

"""Make sinusoidal weather with correct Mediterranean summer peak."""
sw_make(T_mean, amplitude) = SinusoidalWeather(T_mean, amplitude; phase=PHASE, radiation=22.0)

"""Compute DDa (Apr–Sep), DDb (bloom–Aug), and bloom day."""
function compute_climate_indices(T_mean, amplitude; rain_mm=RAIN_MM)
    sw = sw_make(T_mean, amplitude)
    dda = 0.0
    for d in 91:273
        w = get_weather(sw, d)
        dda += degree_days(olive_dev, w.T_mean)
    end
    cdd_est = 0.0; bloom_day = 120
    for d in 1:365
        w = get_weather(sw, d)
        cdd_est += degree_days(olive_dev, w.T_mean)
        if cdd_est >= DD_BLOOM
            bloom_day = d; break
        end
    end
    ddb = 0.0
    for d in bloom_day:243
        w = get_weather(sw, d)
        ddb += degree_days(olive_dev, w.T_mean)
    end
    return (; dda, ddb, bloom_day)
end

"""Predict log₁₀(pupae/tree) from regression (Ponti et al. 2009 Eqn 4)."""
function predict_fly_log10(dda, ddb, bloom_day, rain_mm)
    return FLY_β_DDa * dda + FLY_β_DDb * ddb +
           FLY_β_BLOOM * bloom_day + FLY_β_RAIN * rain_mm
end

"""Full regional simulation: yield, fly density, actual yield."""
function simulate_region(T_mean, amplitude, rain_mm)
    idx = compute_climate_indices(T_mean, amplitude; rain_mm)
    pot_y = max(0.0, predict_yield(olive_yield_model, idx.dda, rain_mm))
    log10_fly = predict_fly_log10(idx.dda, idx.ddb, idx.bloom_day, rain_mm)
    pupae = 10.0^log10_fly   # pupae per tree
    act_y = actual_yield(fly_damage, pupae, pot_y)
    return (; yield=act_y, fly=pupae, pot_yield=pot_y,
              bloom=idx.bloom_day, dda=idx.dda, ddb=idx.ddb, log10_fly)
end

"""Build olive fly population with linear DD rates (for PBDM simulation)."""
function make_olive_fly(; n_adults=0.0)
    k_e, k_l, k_p, k_a = 2, 8, 8, 15
    [
        LifeStage(:egg,   DistributedDelay(k_e, FLY_DD_EGG;   W0=0.0),         egg_dev,   0.004),
        LifeStage(:larva, DistributedDelay(k_l, FLY_DD_LARVA; W0=0.0),         larva_dev, 0.003),
        LifeStage(:pupa,  DistributedDelay(k_p, FLY_DD_PUPA;  W0=0.0),         pupa_dev,  0.002),
        LifeStage(:adult, DistributedDelay(k_a, FLY_DD_ADULT; W0=n_adults/k_a), adult_dev, 0.002),
    ]
end

# ============================================================
# Figure 1: Development rate — olive phenology vs temperature
# ============================================================

println("── Figure 1: Development rate (olive phenology vs temperature) ──")

Ts = collect(range(0.0, 45.0, length=500))
olive_dd = [degree_days(olive_dev, T) for T in Ts]
fly_briere_r = [development_rate(fly_briere, T) for T in Ts]
fly_linear_r = [degree_days(fly_linear, T) for T in Ts]

# Phenological milestones: bloom and maturity day vs mean annual temperature
mean_temps = collect(range(10.0, 32.0, length=200))
bloom_days_by_temp = Float64[]
maturity_days_by_temp = Float64[]
for Tm in mean_temps
    sw = sw_make(Tm, 8.0)
    cdd = 0.0
    bd, md = 365.0, 365.0
    for d in 1:365
        w = get_weather(sw, d)
        cdd += degree_days(olive_dev, w.T_mean)
        if bd == 365.0 && cdd >= DD_BLOOM
            bd = Float64(d)
        end
        if md == 365.0 && cdd >= DD_MATURITY
            md = Float64(d)
        end
    end
    push!(bloom_days_by_temp, bd)
    push!(maturity_days_by_temp, md)
end

fig1 = Figure(size=(1000, 700))

# Panel A: degree-day accumulation rates
ax1a = Axis(fig1[1, 1],
    title="A: Daily Degree-Day Accumulation",
    xlabel="Temperature (°C)", ylabel="DD / day",
    xlabelsize=13, ylabelsize=13)
lines!(ax1a, Ts, olive_dd, linewidth=2.5, color=:forestgreen, label="Olive (base 9.1°C)")
lines!(ax1a, Ts, fly_linear_r, linewidth=2.5, color=:firebrick, label="Fly linear (base 10°C)")
# Scale Brière rate for visibility (×100)
lines!(ax1a, Ts, fly_briere_r .* 100.0, linewidth=2, color=:darkorange, linestyle=:dash,
       label="Fly Brière ×100")
vlines!(ax1a, [OLIVE_T_BASE], color=(:forestgreen, 0.4), linestyle=:dot, linewidth=1)
vlines!(ax1a, [FLY_T_BASE], color=(:firebrick, 0.4), linestyle=:dot, linewidth=1)
vlines!(ax1a, [FLY_T_UPPER], color=(:firebrick, 0.4), linestyle=:dot, linewidth=1)
text!(ax1a, OLIVE_T_BASE + 0.3, 2.0, text="9.1°C", fontsize=9, color=:forestgreen)
text!(ax1a, FLY_T_BASE + 0.3, 1.0, text="10°C", fontsize=9, color=:firebrick)
text!(ax1a, FLY_T_UPPER - 3.5, 1.0, text="35°C", fontsize=9, color=:firebrick)
xlims!(ax1a, 0, 45)
axislegend(ax1a, position=:lt, framevisible=false)

# Panel B: phenological calendar vs mean temperature
ax1b = Axis(fig1[1, 2],
    title="B: Phenological Timing vs Mean Temperature",
    xlabel="Annual mean temperature (°C)", ylabel="Day of year",
    xlabelsize=13, ylabelsize=13)
lines!(ax1b, mean_temps, bloom_days_by_temp, linewidth=2.5, color=:goldenrod,
       label="Bloom (400 DD)")
lines!(ax1b, mean_temps, maturity_days_by_temp, linewidth=2.5, color=:purple,
       label="Maturity (1800 DD)")
# Paper reference: bloom shifts −9.4 d/°C
ref_idx = findfirst(t -> t >= 18.0, mean_temps)
if ref_idx !== nothing
    ref_bloom = bloom_days_by_temp[ref_idx]
    ref_xs = [14.0, 22.0]
    ref_ys = [ref_bloom + PAPER_BLOOM_SHIFT * (14.0 - 18.0),
              ref_bloom + PAPER_BLOOM_SHIFT * (22.0 - 18.0)]
    lines!(ax1b, ref_xs, ref_ys, linewidth=1.5, color=:gray40, linestyle=:dash,
           label="Paper: −9.4 d/°C")
end
hlines!(ax1b, [120.0, 150.0], color=(:olivedrab, 0.15), linewidth=8)
text!(ax1b, 10.5, 123.0, text="Typical bloom range", fontsize=9, color=:olivedrab)
axislegend(ax1b, position=:rt, framevisible=false)

save(joinpath(figdir, "development_rate.png"), fig1, px_per_unit=2)
println("  Saved development_rate.png")
println("  Olive DD/day at 20°C: $(round(degree_days(olive_dev, 20.0), digits=2))")
println("  Fly Brière peak: $(round(maximum(fly_briere_r), digits=5)) at $(round(Ts[argmax(fly_briere_r)], digits=1))°C")

# ============================================================
# Figure 2: Olive fly dynamics in Mediterranean climate
# ============================================================

println("\n── Figure 2: Olive fly dynamics (Mediterranean) ──")

# Reproduction function: adults produce eggs proportional to population & temperature
const FECUNDITY_MAX = 0.3
function olive_fly_reproduction(pop, w, p, day)
    adult_stage = pop.stages[end]
    n_adult = delay_total(adult_stage.delay)
    T = w.T_mean
    (T < 15.0 || T > 33.0) && return 0.0
    ovi_scale = max(0.0, 1.0 - ((T - 25.0) / 10.0)^2)
    return FECUNDITY_MAX * ovi_scale * n_adult
end

t0_fly, tf_fly = 120, 319   # late April to mid-November
baseline_weather = sw_make(18.0, 8.0)

fly_stages_sim = make_olive_fly(; n_adults=30.0)
fly_pop = Population(:olive_fly, fly_stages_sim)
temps_fly = [get_weather(baseline_weather, d).T_mean for d in t0_fly:tf_fly]
weather_fly = WeatherSeries([DailyWeather(T) for T in temps_fly]; day_offset=t0_fly)
prob_fly = PBDMProblem(DensityDependent(), fly_pop, weather_fly, (t0_fly, tf_fly))
sol_fly = solve(prob_fly, DirectIteration(); reproduction_fn=olive_fly_reproduction)

stage_names = [:egg, :larva, :pupa, :adult]
stage_colors = [:goldenrod, :forestgreen, :steelblue, :firebrick]

# Fruit susceptibility gate overlay (pit hardening → maturity)
dd_cum_gate = Float64[]
let cdd_g = 0.0
    for d in 1:tf_fly
        w = get_weather(baseline_weather, d)
        cdd_g += degree_days(olive_dev, w.T_mean)
        if d >= t0_fly
            push!(dd_cum_gate, cdd_g)
        end
    end
end
fruit_gate = [clamp((cdd - DD_PIT_HARDEN) / (DD_MATURITY - DD_PIT_HARDEN), 0.0, 1.0)
              for cdd in dd_cum_gate]

fig2 = Figure(size=(1000, 650))
ax2a = Axis(fig2[1, 1],
    title="B. oleae Population Dynamics — Baseline Mediterranean (18°C mean)",
    subtitle="cf. Ponti et al. (2009) Fig. 3: Villacidro, Sardinia",
    xlabel="Day of year", ylabel="Population",
    xlabelsize=13, ylabelsize=13)
ax2b = Axis(fig2[1, 1],
    ylabel="Fruit susceptibility",
    yaxisposition=:right, ylabelsize=12,
    yticklabelcolor=(:olivedrab, 0.7))
hidespines!(ax2b)
hidexdecorations!(ax2b)

for (i, sname) in enumerate(stage_names)
    traj = sol_fly.stage_totals[i, :]
    lines!(ax2a, sol_fly.t, traj, linewidth=2.5, color=stage_colors[i],
           label=String(sname))
end

band!(ax2b, collect(t0_fly:tf_fly), zeros(length(fruit_gate)), fruit_gate,
      color=(:olivedrab, 0.12))
lines!(ax2b, collect(t0_fly:tf_fly), fruit_gate,
       color=(:olivedrab, 0.5), linewidth=1.5, linestyle=:dash, label="Fruit gate")

xlims!(ax2a, t0_fly, tf_fly)
xlims!(ax2b, t0_fly, tf_fly)
ylims!(ax2b, -0.05, 1.15)
axislegend(ax2a, position=:lt, framevisible=false)
axislegend(ax2b, position=:rt, framevisible=false)

text!(ax2a, tf_fly - 80, maximum(sol_fly.stage_totals) * 0.85,
      text="Paper: pupae peak in\ncoastal Sardinia ≈ day 220–260",
      fontsize=9, color=:gray40)

save(joinpath(figdir, "fly_dynamics.png"), fig2, px_per_unit=2)
println("  Saved fly_dynamics.png")
for (i, sname) in enumerate(stage_names)
    traj = sol_fly.stage_totals[i, :]
    pk = maximum(traj)
    pk_day = sol_fly.t[argmax(traj)]
    println("  $(sname): peak=$(round(pk, digits=1)) at day $(pk_day)")
end

# ============================================================
# Figure 3: Climate scenarios — yield under warming
# ============================================================

println("\n── Figure 3: Climate scenarios (yield projections) ──")

warming_levels = [0.0, 1.0, 1.8, 2.5, 3.6]
warming_labels = ["Baseline", "+1°C", "+1.8°C", "+2.5°C", "+3.6°C"]

# Regional baselines from Ponti et al. (2014) supplementary
regions = [
    (name="S. Europe",      T_mean=18.0, amplitude=8.0,  rain=250.0),
    (name="Greece/Turkey",  T_mean=19.5, amplitude=9.0,  rain=180.0),
    (name="North Africa",   T_mean=21.0, amplitude=7.0,  rain=120.0),
    (name="Middle East",    T_mean=23.0, amplitude=10.0, rain=80.0),
    (name="Iberia (north)", T_mean=16.5, amplitude=7.5,  rain=320.0),
]

# Compute yields and fly pressure for all regions × warming
yield_matrix = zeros(length(regions), length(warming_levels))
fly_matrix   = zeros(length(regions), length(warming_levels))
bloom_matrix = zeros(length(regions), length(warming_levels))

for (ri, reg) in enumerate(regions)
    for (wi, ΔT) in enumerate(warming_levels)
        res = simulate_region(reg.T_mean + ΔT, reg.amplitude, reg.rain)
        yield_matrix[ri, wi] = res.yield
        fly_matrix[ri, wi]   = res.fly
        bloom_matrix[ri, wi] = res.bloom
    end
end

fig3 = Figure(size=(1100, 750))
region_colors = [:forestgreen, :steelblue, :goldenrod, :firebrick, :purple]

# Panel A: Yield by region across warming levels
ax3a = Axis(fig3[1, 1],
    title="A: Predicted Yield Under Climate Warming",
    subtitle="Ponti et al. (2014): basin-wide +$(PAPER_YIELD_CHANGE_1_8)% at +1.8°C",
    xlabel="Warming (°C)", ylabel="Actual yield (kg/ha, after fly damage)",
    xlabelsize=13, ylabelsize=13,
    xticks=(1:length(warming_levels), warming_labels))
for (ri, reg) in enumerate(regions)
    lines!(ax3a, 1:length(warming_levels), yield_matrix[ri, :],
           linewidth=2.5, color=region_colors[ri], label=reg.name)
    scatter!(ax3a, 1:length(warming_levels), yield_matrix[ri, :],
             color=region_colors[ri], markersize=8)
end
baseline_avg = mean(yield_matrix[:, 1])
warm18_avg = mean(yield_matrix[:, 3])
model_yield_change = baseline_avg > 0 ? (warm18_avg - baseline_avg) / baseline_avg * 100 : 0.0
text!(ax3a, 3.3, minimum(yield_matrix) + 200,
      text="Basin avg Δ: $(round(model_yield_change, digits=1))%\nPaper: +$(PAPER_YIELD_CHANGE_1_8)%",
      fontsize=10, color=:gray30)
axislegend(ax3a, position=:lt, framevisible=false)

# Panel B: Fly density by region
ax3b = Axis(fig3[1, 2],
    title="B: Predicted Fly Density (pupae/tree)",
    subtitle="Ponti et al. (2014): basin-wide $(PAPER_FLY_CHANGE_1_8)% at +1.8°C",
    xlabel="Warming (°C)", ylabel="Fly density (pupae/tree)",
    xlabelsize=13, ylabelsize=13,
    xticks=(1:length(warming_levels), warming_labels))
for (ri, reg) in enumerate(regions)
    lines!(ax3b, 1:length(warming_levels), fly_matrix[ri, :],
           linewidth=2.5, color=region_colors[ri], label=reg.name)
    scatter!(ax3b, 1:length(warming_levels), fly_matrix[ri, :],
             color=region_colors[ri], markersize=8)
end
fly_base_avg = mean(fly_matrix[:, 1])
fly_warm_avg = mean(fly_matrix[:, 3])
model_fly_change = fly_base_avg > 0 ? (fly_warm_avg - fly_base_avg) / fly_base_avg * 100 : 0.0
text!(ax3b, 3.3, maximum(fly_matrix) * 0.85,
      text="Basin avg Δ: $(round(model_fly_change, digits=1))%\nPaper: $(PAPER_FLY_CHANGE_1_8)%",
      fontsize=10, color=:gray30)
axislegend(ax3b, position=:rt, framevisible=false)

# Panel C: Bloom date shift
ax3c = Axis(fig3[2, 1:2],
    title="C: Bloom Date Shift vs Warming",
    subtitle="Paper reference: −9.4 days/°C (Ponti et al. 2009)",
    xlabel="Warming (°C)", ylabel="Bloom day of year",
    xlabelsize=13, ylabelsize=13)
for (ri, reg) in enumerate(regions)
    lines!(ax3c, warming_levels, bloom_matrix[ri, :],
           linewidth=2.5, color=region_colors[ri], label=reg.name)
    scatter!(ax3c, warming_levels, bloom_matrix[ri, :],
             color=region_colors[ri], markersize=8)
end
ref_bloom_base = bloom_matrix[1, 1]
lines!(ax3c, [0.0, 3.6], [ref_bloom_base, ref_bloom_base + PAPER_BLOOM_SHIFT * 3.6],
       linewidth=2, color=:gray40, linestyle=:dash, label="Paper: −9.4 d/°C")
axislegend(ax3c, position=:rt, framevisible=false, nbanks=2)

save(joinpath(figdir, "climate_scenarios.png"), fig3, px_per_unit=2)
println("  Saved climate_scenarios.png")
println("  Basin-wide yield change at +1.8°C: $(round(model_yield_change, digits=1))% (paper: +$(PAPER_YIELD_CHANGE_1_8)%)")
println("  Basin-wide fly change at +1.8°C:   $(round(model_fly_change, digits=1))% (paper: $(PAPER_FLY_CHANGE_1_8)%)")

# ============================================================
# Figure 4: Fruit infestation at different temperatures
# ============================================================

println("\n── Figure 4: Fly damage (infestation vs temperature) ──")

# Panel A: PBDM fly simulation at a range of constant temperatures
test_temps = collect(range(12.0, 36.0, length=25))
peak_adults_by_temp = Float64[]
cumul_pupae_by_temp = Float64[]

for Tc in test_temps
    n_days_c = 180
    weather_c = WeatherSeries([DailyWeather(Tc) for _ in 1:n_days_c])
    fly_c = Population(:olive_fly, make_olive_fly(; n_adults=30.0))
    prob_c = PBDMProblem(fly_c, weather_c, (1, n_days_c))
    sol_c = solve(prob_c, DirectIteration())
    push!(peak_adults_by_temp, maximum(sol_c.stage_totals[4, :]))
    push!(cumul_pupae_by_temp, sum(sol_c.stage_totals[3, :]))
end

# Panel B: Seasonal infestation — cumulative pupa production under warming
seasonal_scenarios = [("Baseline", 0.0, :black, :solid),
                      ("+1°C",     1.0, :dodgerblue, :dash),
                      ("+2°C",     2.0, :red, :dot)]

fig4 = Figure(size=(1000, 700))

ax4a = Axis(fig4[1, 1],
    title="A: Fly Peak Adults vs Constant Temperature",
    subtitle="180-day DI simulation, 30 initial adults",
    xlabel="Temperature (°C)", ylabel="Peak adult population",
    xlabelsize=13, ylabelsize=13)
lines!(ax4a, test_temps, peak_adults_by_temp, linewidth=2.5, color=:firebrick)
scatter!(ax4a, test_temps, peak_adults_by_temp, color=:firebrick, markersize=8)
vlines!(ax4a, [FLY_T_BASE], color=(:firebrick, 0.3), linestyle=:dot)
vlines!(ax4a, [FLY_T_UPPER], color=(:firebrick, 0.3), linestyle=:dot)
text!(ax4a, FLY_T_BASE + 0.3, maximum(peak_adults_by_temp) * 0.92,
      text="Base 10°C", fontsize=9, color=:firebrick)
text!(ax4a, FLY_T_UPPER - 3.0, maximum(peak_adults_by_temp) * 0.92,
      text="Upper 35°C", fontsize=9, color=:firebrick)
opt_idx = argmax(peak_adults_by_temp)
println("  Fly optimal T (peak adults): $(round(test_temps[opt_idx], digits=1))°C")

ax4b = Axis(fig4[1, 2],
    title="B: Cumulative Pupae vs Constant Temperature",
    subtitle="Higher total ≈ more seasonal generations",
    xlabel="Temperature (°C)", ylabel="Cumulative pupae (arb. units)",
    xlabelsize=13, ylabelsize=13)
lines!(ax4b, test_temps, cumul_pupae_by_temp, linewidth=2.5, color=:steelblue)
scatter!(ax4b, test_temps, cumul_pupae_by_temp, color=:steelblue, markersize=8)
# Paper: fruit attack range annotation
hspan!(ax4b, 0.0, maximum(cumul_pupae_by_temp),
       color=(:gray, 0.0))  # invisible, just for alignment
text!(ax4b, 30.0, maximum(cumul_pupae_by_temp) * 0.85,
      text="Paper obs. fruit attack:\n$(PAPER_FRUIT_ATTACK_RANGE[1])–$(PAPER_FRUIT_ATTACK_RANGE[2])%",
      fontsize=9, color=:gray40)

ax4c = Axis(fig4[2, 1:2],
    title="C: Seasonal Cumulative Pupa Production — Mediterranean Climate",
    xlabel="Day of year", ylabel="Cumulative pupae (arb. units)",
    xlabelsize=13, ylabelsize=13)
for (label, ΔT, col, sty) in seasonal_scenarios
    sw_t = sw_make(18.0 + ΔT, 8.0)
    fly_t = Population(:olive_fly, make_olive_fly(; n_adults=30.0))
    temps_t = [get_weather(sw_t, d).T_mean for d in t0_fly:tf_fly]
    weather_t = WeatherSeries([DailyWeather(T) for T in temps_t]; day_offset=t0_fly)
    prob_t = PBDMProblem(DensityDependent(), fly_t, weather_t, (t0_fly, tf_fly))
    sol_t = solve(prob_t, DirectIteration(); reproduction_fn=olive_fly_reproduction)
    cum_pupa = cumsum(sol_t.stage_totals[3, :])
    lines!(ax4c, sol_t.t, cum_pupa, linewidth=2.5, color=col, linestyle=sty, label=label)
    println("  $(label): final cumul. pupae = $(round(cum_pupa[end], digits=1))")
end
axislegend(ax4c, position=:lt, framevisible=false)

save(joinpath(figdir, "fly_damage.png"), fig4, px_per_unit=2)
println("  Saved fly_damage.png")

# ============================================================
# Figure 5: Geographic comparison — Sardinia vs Northern Italy
# ============================================================

println("\n── Figure 5: Geographic comparison (Sardinia vs N. Italy) ──")

# Sardinia (Villacidro): hot, dry, coastal (Ponti et al. 2009: lat 39.43°N)
sardinia = (name="Sardinia (Villacidro)", T_mean=17.5, amplitude=10.0, rain=200.0)
# Northern Italy (Liguria / Po margin): cooler, wetter
north_italy = (name="N. Italy (Liguria)", T_mean=14.5, amplitude=9.0, rain=450.0)

locations = [sardinia, north_italy]
loc_colors = [:darkorange, :steelblue]

fig5 = Figure(size=(1100, 750))

# Panel A: Fly population trajectories
ax5a = Axis(fig5[1, 1],
    title="A: Olive Fly Seasonal Dynamics",
    xlabel="Day of year", ylabel="Total fly population",
    xlabelsize=13, ylabelsize=13)
for (li, loc) in enumerate(locations)
    sw_l = sw_make(loc.T_mean, loc.amplitude)
    fly_l = Population(:olive_fly, make_olive_fly(; n_adults=30.0))
    temps_l = [get_weather(sw_l, d).T_mean for d in t0_fly:tf_fly]
    weather_l = WeatherSeries([DailyWeather(T) for T in temps_l]; day_offset=t0_fly)
    prob_l = PBDMProblem(DensityDependent(), fly_l, weather_l, (t0_fly, tf_fly))
    sol_l = solve(prob_l, DirectIteration(); reproduction_fn=olive_fly_reproduction)
    total_pop_l = [sum(sol_l.stage_totals[:, j]) for j in 1:size(sol_l.stage_totals, 2)]
    lines!(ax5a, sol_l.t, total_pop_l, linewidth=2.5, color=loc_colors[li],
           label=loc.name)
    pk = maximum(total_pop_l)
    pk_day = sol_l.t[argmax(total_pop_l)]
    println("  $(loc.name): peak total=$(round(pk, digits=1)) at day $(pk_day)")
end
axislegend(ax5a, position=:lt, framevisible=false)

# Panel B: Yield under warming at both locations
ax5b = Axis(fig5[1, 2],
    title="B: Yield Response to Warming",
    xlabel="Warming (°C)", ylabel="Actual yield (kg/ha)",
    xlabelsize=13, ylabelsize=13)
for (li, loc) in enumerate(locations)
    yields_loc = [simulate_region(loc.T_mean + ΔT, loc.amplitude, loc.rain).yield
                  for ΔT in warming_levels]
    lines!(ax5b, warming_levels, yields_loc, linewidth=2.5, color=loc_colors[li],
           label=loc.name)
    scatter!(ax5b, warming_levels, yields_loc, color=loc_colors[li], markersize=8)
end
axislegend(ax5b, position=:lt, framevisible=false)

# Panel C: Annual temperature profiles
ax5c = Axis(fig5[2, 1],
    title="C: Annual Temperature Profiles",
    xlabel="Day of year", ylabel="Temperature (°C)",
    xlabelsize=13, ylabelsize=13)
for (li, loc) in enumerate(locations)
    sw_l = sw_make(loc.T_mean, loc.amplitude)
    temps_profile = [get_weather(sw_l, d).T_mean for d in 1:365]
    lines!(ax5c, 1:365, temps_profile, linewidth=2.5, color=loc_colors[li],
           label=loc.name)
end
hlines!(ax5c, [FLY_T_UPPER], color=(:firebrick, 0.4), linestyle=:dot)
text!(ax5c, 10, FLY_T_UPPER + 0.5, text="Fly upper limit (35°C)", fontsize=9,
      color=:firebrick)
hlines!(ax5c, [OLIVE_T_BASE], color=(:forestgreen, 0.4), linestyle=:dot)
text!(ax5c, 10, OLIVE_T_BASE + 0.5, text="Olive base (9.1°C)", fontsize=9,
      color=:forestgreen)
axislegend(ax5c, position=:lt, framevisible=false)

# Panel D: Winners/losers bar chart at +1.8 °C
ax5d = Axis(fig5[2, 2],
    title="D: Winners and Losers Summary (+1.8°C)",
    xlabel="Location", ylabel="Change (%)",
    xlabelsize=13, ylabelsize=13,
    xticks=(1:2, [sardinia.name, north_italy.name]))
for (li, loc) in enumerate(locations)
    base_res = simulate_region(loc.T_mean, loc.amplitude, loc.rain)
    warm_res = simulate_region(loc.T_mean + 1.8, loc.amplitude, loc.rain)
    yield_ch = base_res.yield > 0 ? (warm_res.yield - base_res.yield) / base_res.yield * 100 : 0.0
    fly_ch = base_res.fly > 0 ? (warm_res.fly - base_res.fly) / base_res.fly * 100 : 0.0
    barplot!(ax5d, [li - 0.15], [yield_ch], width=0.25,
             color=:forestgreen, label=(li == 1 ? "Yield Δ%" : nothing))
    barplot!(ax5d, [li + 0.15], [fly_ch], width=0.25,
             color=:firebrick, label=(li == 1 ? "Fly Δ%" : nothing))
    println("  $(loc.name) at +1.8°C: yield Δ=$(round(yield_ch, digits=1))%, fly Δ=$(round(fly_ch, digits=1))%")
end
hlines!(ax5d, [0.0], color=:gray50, linewidth=1)
text!(ax5d, 0.55, -8.0,
      text="Paper: Sardinia fly ↓\nin hot lowlands;\nN. Italy yield ↑",
      fontsize=8, color=:gray40)
axislegend(ax5d, position=:lt, framevisible=false)

save(joinpath(figdir, "geographic_comparison.png"), fig5, px_per_unit=2)
println("  Saved geographic_comparison.png")

# ============================================================
# Summary validation table
# ============================================================

println("\n" * "="^70)
println("VALIDATION SUMMARY")
println("="^70)
println("  Parameter                  | Vignette/Model  | Paper reference")
println("-"^70)
println("  Olive base temp            | $(OLIVE_T_BASE)°C         | 9.1°C (Ponti 2009)")
println("  Bloom DD threshold         | $(DD_BLOOM) DD        | 400 DD (calibrated)")
println("  Fly base temp (egg)        | $(FLY_T_EGG)°C         | 9.2°C (Gutierrez 2009)")
println("  Fly base temp (adult)      | $(FLY_T_ADULT)°C        | 10°C (Gutierrez 2009)")
println("  Fly upper temp             | $(FLY_T_UPPER)°C        | 35°C (Gutierrez 2009)")
println("  Bloom shift                | simulated       | −9.4 d/°C (Ponti 2009)")
println("  Yield change +1.8°C        | $(round(model_yield_change, digits=1))%         | +$(PAPER_YIELD_CHANGE_1_8)% (Ponti 2014)")
println("  Fly change +1.8°C          | $(round(model_fly_change, digits=1))%        | $(PAPER_FLY_CHANGE_1_8)% (Ponti 2014)")
println("  Fruit attack (obs. range)  | model computed  | 38–87% (Ponti 2009)")
println("  Yield regression R²        | (from paper)    | 0.47 (Ponti 2009 Eqn 3)")
println("  Fly regression R²          | (from paper)    | 0.80 (Ponti 2009 Eqn 4)")
println("="^70)

println("\n✓ All olive climate validation figures saved to: $(figdir)")
