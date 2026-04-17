#!/usr/bin/env julia
# Validation script for the Olive Fruit Fly (Bactrocera oleae) PBDM
# on olive (Olea europaea) in Puglia, Italy.
#
# Literature references:
#   Gutierrez et al. (2009) Ecol. Model. 220:1579–1597
#   Ponti et al. (2014) Ecol. Model. 296:66–79
#
# Generates five PNG figures in scripts/figures/olive_fly/.

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using PhysiologicallyBasedDemographicModels
using CairoMakie

figdir = joinpath(@__DIR__, "figures", "olive_fly")
mkpath(figdir)

# ============================================================
# Parameters from Gutierrez et al. (2009), Ponti et al. (2014)
# ============================================================

# Stage-specific lower developmental thresholds (°C)
const T_LOWER_EGG   = 9.2
const T_LOWER_LARVA = 13.9
const T_LOWER_PUPA  = 12.4
const T_LOWER_ADULT = 10.0

# Common upper threshold (°C)
const T_UPPER = 35.0

# Degree-day requirements per stage
const DD_EGG   = 65.0
const DD_LARVA = 220.0
const DD_PUPA  = 220.0
const DD_ADULT = 600.0   # reproductive adult lifespan

# Brière development rate models (for plotting rate curves — Fig. 1)
# r(T) = a * T * (T - T_lower) * sqrt(T_upper - T)
egg_briere   = BriereDevelopmentRate(2.0e-4,  T_LOWER_EGG,   T_UPPER)
larva_briere = BriereDevelopmentRate(7.0e-5,  T_LOWER_LARVA, T_UPPER)
pupa_briere  = BriereDevelopmentRate(7.0e-5,  T_LOWER_PUPA,  T_UPPER)

# Linear development rate models (for simulation — degree-day accumulation)
egg_dev   = LinearDevelopmentRate(T_LOWER_EGG,   T_UPPER)
larva_dev = LinearDevelopmentRate(T_LOWER_LARVA, T_UPPER)
pupa_dev  = LinearDevelopmentRate(T_LOWER_PUPA,  T_UPPER)
adult_dev = LinearDevelopmentRate(T_LOWER_ADULT, T_UPPER)

# Background mortality rates (per degree-day)
const μ_EGG   = 0.004
const μ_LARVA = 0.003
const μ_PUPA  = 0.002
const μ_ADULT = 0.002

# Distributed delay substage counts (chosen for numerical stability:
# need (k/τ)*dd_max < 1 where dd_max = T_UPPER - T_lower)
const K_EGG   = 2    # τ/dd_max = 65/25.8 = 2.5
const K_LARVA = 8    # τ/dd_max = 220/21.1 = 10.4
const K_PUPA  = 8    # τ/dd_max = 220/22.6 = 9.7
const K_ADULT = 15   # τ/dd_max = 600/25.0 = 24

# Puglia seasonal temperature model:
#   T(d) = T_mean + amplitude * sin(2π(d - phase_offset) / 365)
const PUGLIA_T_MEAN    = 18.0   # annual mean (°C)
const PUGLIA_AMPLITUDE = 11.0   # annual half-range (°C)
const PUGLIA_PHASE     = 100.0  # day offset (peak ≈ day 191)

puglia_temp(d) = PUGLIA_T_MEAN + PUGLIA_AMPLITUDE * sin(2π * (d - PUGLIA_PHASE) / 365)

# ============================================================
# Helper: build olive fly population with a given initial state
# ============================================================

function make_olive_fly(; n_eggs=0.0, n_larvae=0.0, n_pupae=0.0, n_adults=0.0)
    stages = [
        LifeStage(:egg,   DistributedDelay(K_EGG,   DD_EGG;   W0=n_eggs   / K_EGG),
                  egg_dev,   μ_EGG),
        LifeStage(:larva, DistributedDelay(K_LARVA, DD_LARVA; W0=n_larvae / K_LARVA),
                  larva_dev, μ_LARVA),
        LifeStage(:pupa,  DistributedDelay(K_PUPA,  DD_PUPA;  W0=n_pupae  / K_PUPA),
                  pupa_dev,  μ_PUPA),
        LifeStage(:adult, DistributedDelay(K_ADULT, DD_ADULT; W0=n_adults / K_ADULT),
                  adult_dev, μ_ADULT),
    ]
    Population(:bactrocera_oleae, stages)
end

# ============================================================
# Figure 1: Stage-specific development rate curves
# ============================================================

println("── Figure 1: Development rate curves ──")

Ts = range(0.0, 40.0, length=400)
egg_rates   = [development_rate(egg_briere,   T) for T in Ts]
larva_rates = [development_rate(larva_briere, T) for T in Ts]
pupa_rates  = [development_rate(pupa_briere,  T) for T in Ts]

# Published data points (approximate from Gutierrez et al. 2009)
obs_egg_T   = [20.0, 25.0];  obs_egg_r   = [0.06, 0.12]
obs_larva_T = [20.0, 25.0];  obs_larva_r = [0.015, 0.035]
obs_pupa_T  = [20.0, 25.0];  obs_pupa_r  = [0.02, 0.04]

fig1 = Figure(size=(900, 600))
ax1 = Axis(fig1[1, 1],
    title="Bactrocera oleae — Stage-Specific Development Rates (Brière model)",
    xlabel="Temperature (°C)", ylabel="Development rate (1/day)",
    xlabelsize=14, ylabelsize=14)

lines!(ax1, collect(Ts), egg_rates,   linewidth=2.5, color=:goldenrod, label="Egg (T_L=9.2°C)")
lines!(ax1, collect(Ts), larva_rates, linewidth=2.5, color=:forestgreen, label="Larva (T_L=13.9°C)")
lines!(ax1, collect(Ts), pupa_rates,  linewidth=2.5, color=:steelblue, label="Pupa (T_L=12.4°C)")

scatter!(ax1, obs_egg_T, obs_egg_r, color=:goldenrod, markersize=12, marker=:circle,
         label="Egg obs.")
scatter!(ax1, obs_larva_T, obs_larva_r, color=:forestgreen, markersize=12, marker=:utriangle,
         label="Larva obs.")
scatter!(ax1, obs_pupa_T, obs_pupa_r, color=:steelblue, markersize=12, marker=:diamond,
         label="Pupa obs.")

xlims!(ax1, 0, 40)
ylims!(ax1, 0, nothing)
axislegend(ax1, position=:lt)

save(joinpath(figdir, "devrate_curves.png"), fig1, px_per_unit=2)
println("  Saved devrate_curves.png")
println("  Egg peak rate:   $(round(maximum(egg_rates),   digits=4)) at $(round(Ts[argmax(egg_rates)],   digits=1))°C")
println("  Larva peak rate: $(round(maximum(larva_rates), digits=4)) at $(round(Ts[argmax(larva_rates)], digits=1))°C")
println("  Pupa peak rate:  $(round(maximum(pupa_rates),  digits=4)) at $(round(Ts[argmax(pupa_rates)],  digits=1))°C")

# ============================================================
# Figure 2: Olive fruit susceptibility gate (Puglia)
# ============================================================

println("\n── Figure 2: Fruit susceptibility gate ──")

# Olive phenology for Puglia:
#   Fruit set: day 140
#   Gate opens at 400 DD above 9°C base, fully open at 800 DD
const FRUIT_SET_DAY = 140
const GATE_BASE_T   = 9.0    # base temperature for DD accumulation toward gate
const GATE_DD_OPEN  = 400.0  # DD after fruit set when gate starts opening
const GATE_DD_FULL  = 800.0  # DD after fruit set when gate is fully open

days_year = 1:365
temps_year = [puglia_temp(d) for d in days_year]

# Accumulate DD from fruit set day
dd_fruit = zeros(365)
let cumdd = 0.0
    for d in FRUIT_SET_DAY:365
        cumdd += max(0.0, puglia_temp(d) - GATE_BASE_T)
        dd_fruit[d] = cumdd
    end
end

# Fruit gate: ramp from 0 to 1 between GATE_DD_OPEN and GATE_DD_FULL
fruit_gate = zeros(365)
for d in 1:365
    if d < FRUIT_SET_DAY
        fruit_gate[d] = 0.0
    else
        dd = dd_fruit[d]
        if dd < GATE_DD_OPEN
            fruit_gate[d] = 0.0
        elseif dd < GATE_DD_FULL
            fruit_gate[d] = (dd - GATE_DD_OPEN) / (GATE_DD_FULL - GATE_DD_OPEN)
        else
            fruit_gate[d] = 1.0
        end
    end
end

gate_open_day = findfirst(x -> x > 0, fruit_gate)
gate_full_day = findfirst(x -> x >= 1.0, fruit_gate)

fig2 = Figure(size=(900, 550))
ax2a = Axis(fig2[1, 1],
    title="Olive Fruit Susceptibility Gate — Puglia, Italy",
    xlabel="Day of year", ylabel="Temperature (°C)",
    xlabelsize=14, ylabelsize=14,
    yticklabelcolor=:firebrick)
ax2b = Axis(fig2[1, 1],
    ylabel="Fruit gate (0–1)",
    yaxisposition=:right,
    ylabelsize=14,
    yticklabelcolor=:forestgreen)

hidespines!(ax2b)
hidexdecorations!(ax2b)

lines!(ax2a, collect(days_year), temps_year, color=:firebrick, linewidth=2, label="Temperature")
lines!(ax2b, collect(days_year), fruit_gate, color=:forestgreen, linewidth=2.5, label="Fruit gate")

# Annotate key events
vlines!(ax2a, [FRUIT_SET_DAY], color=:gray, linestyle=:dash, linewidth=1)
text!(ax2a, FRUIT_SET_DAY + 3, minimum(temps_year) + 1,
      text="Fruit set\n(day $FRUIT_SET_DAY)", fontsize=10, color=:gray)
if gate_open_day !== nothing
    vlines!(ax2a, [gate_open_day], color=(:forestgreen, 0.5), linestyle=:dot, linewidth=1)
end
if gate_full_day !== nothing
    vlines!(ax2a, [gate_full_day], color=(:forestgreen, 0.5), linestyle=:dot, linewidth=1)
end

xlims!(ax2a, 1, 365)
xlims!(ax2b, 1, 365)
ylims!(ax2b, -0.05, 1.15)

axislegend(ax2a, position=:lt)
axislegend(ax2b, position=:rt)

save(joinpath(figdir, "fruit_gate.png"), fig2, px_per_unit=2)
println("  Saved fruit_gate.png")
println("  Gate starts opening: day $(gate_open_day !== nothing ? gate_open_day : "N/A")")
println("  Gate fully open:     day $(gate_full_day !== nothing ? gate_full_day : "N/A")")

# ============================================================
# Figure 3: Cohort simulation at constant 25°C
# ============================================================

println("\n── Figure 3: Constant 25°C simulation ──")

n_days_sim = 200
const_T = 25.0

fly_25 = make_olive_fly(n_eggs=200.0)
weather_25 = WeatherSeries([DailyWeather(const_T) for _ in 1:n_days_sim])
prob_25 = PBDMProblem(fly_25, weather_25, (1, n_days_sim))
sol_25 = solve(prob_25, DirectIteration())

stage_names = [:egg, :larva, :pupa, :adult]
stage_colors = [:goldenrod, :forestgreen, :steelblue, :firebrick]

fig3 = Figure(size=(900, 600))
ax3 = Axis(fig3[1, 1],
    title="B. oleae Cohort at Constant 25°C — 200 Eggs Initial",
    xlabel="Day", ylabel="Population",
    xlabelsize=14, ylabelsize=14)

for (i, sname) in enumerate(stage_names)
    traj = sol_25.stage_totals[i, :]
    lines!(ax3, sol_25.t, traj, linewidth=2.5, color=stage_colors[i],
           label=String(sname))
    pk = maximum(traj)
    pk_day = sol_25.t[argmax(traj)]
    println("  $(sname): peak=$(round(pk, digits=1)) at day $(pk_day)")
end

xlims!(ax3, 1, n_days_sim)
axislegend(ax3, position=:rt)

save(joinpath(figdir, "sim_constant_25C.png"), fig3, px_per_unit=2)
println("  Saved sim_constant_25C.png")
println("  Return code: $(sol_25.retcode)")

# ============================================================
# Figure 4: Seasonal Puglia simulation (with reproduction)
# ============================================================

println("\n── Figure 4: Seasonal Puglia simulation ──")

t0_season = 120   # late April start
tf_season = 319   # mid-November end
n_season = tf_season - t0_season

# Reproduction: adults produce eggs proportional to adult population and temperature
# Fecundity scaled to produce realistic seasonal dynamics
const FECUNDITY_MAX = 0.3   # eggs per adult per day (moderate oviposition)

function olive_fly_reproduction(pop, w, p, day)
    adult_stage = pop.stages[end]
    n_adult = delay_total(adult_stage.delay)
    T = w.T_mean
    # Temperature-dependent oviposition: active only between 15–33°C
    if T < 15.0 || T > 33.0
        return 0.0
    end
    # Bell-shaped scaling peaking at ~25°C
    ovi_scale = max(0.0, 1.0 - ((T - 25.0) / 10.0)^2)
    return FECUNDITY_MAX * ovi_scale * n_adult
end

# Generate Puglia weather for the fly season
temps_season = [puglia_temp(d) for d in t0_season:tf_season]
weather_season = WeatherSeries([DailyWeather(T) for T in temps_season]; day_offset=t0_season)

fly_puglia = make_olive_fly(n_adults=30.0)
prob_puglia = PBDMProblem(DensityDependent(), fly_puglia, weather_season,
                          (t0_season, tf_season))
sol_puglia = solve(prob_puglia, DirectIteration();
                   reproduction_fn=olive_fly_reproduction)

# Compute fruit gate for the season window
gate_season = fruit_gate[t0_season:tf_season]

fig4 = Figure(size=(1000, 650))
ax4a = Axis(fig4[1, 1],
    title="B. oleae Seasonal Dynamics — Puglia, Italy (Day $(t0_season)–$(tf_season))",
    xlabel="Day of year", ylabel="Population",
    xlabelsize=14, ylabelsize=14)
ax4b = Axis(fig4[1, 1],
    ylabel="Fruit gate",
    yaxisposition=:right,
    ylabelsize=14,
    yticklabelcolor=(:forestgreen, 0.6))
hidespines!(ax4b)
hidexdecorations!(ax4b)

for (i, sname) in enumerate(stage_names)
    traj = sol_puglia.stage_totals[i, :]
    lines!(ax4a, sol_puglia.t, traj, linewidth=2.5, color=stage_colors[i],
           label=String(sname))
end

# Overlay fruit gate as shaded band
band!(ax4b, collect(t0_season:tf_season), zeros(length(gate_season)), gate_season,
      color=(:forestgreen, 0.15), label="Fruit gate")
lines!(ax4b, collect(t0_season:tf_season), gate_season,
       color=(:forestgreen, 0.5), linewidth=1.5, linestyle=:dash)

xlims!(ax4a, t0_season, tf_season)
xlims!(ax4b, t0_season, tf_season)
ylims!(ax4b, -0.05, 1.15)

axislegend(ax4a, position=:lt)
axislegend(ax4b, position=:rt)

save(joinpath(figdir, "seasonal_puglia.png"), fig4, px_per_unit=2)
println("  Saved seasonal_puglia.png")
for (i, sname) in enumerate(stage_names)
    traj = sol_puglia.stage_totals[i, :]
    pk = maximum(traj)
    pk_day = sol_puglia.t[argmax(traj)]
    println("  $(sname): peak=$(round(pk, digits=1)) at day $(pk_day)")
end

# ============================================================
# Figure 5: Climate warming comparison
# ============================================================

println("\n── Figure 5: Warming comparison ──")

warming_offsets = [0.0, 1.0, 2.0]
warming_labels  = ["Baseline", "+1°C", "+2°C"]
warming_colors  = [:black, :darkorange, :red]
warming_styles  = [:solid, :dash, :dot]

fig5 = Figure(size=(1000, 600))
ax5 = Axis(fig5[1, 1],
    title="B. oleae — Climate Warming Impact on Total Population (Puglia)",
    xlabel="Day of year", ylabel="Total fly population",
    xlabelsize=14, ylabelsize=14)

for (idx, offset) in enumerate(warming_offsets)
    temps_w = [puglia_temp(d) + offset for d in t0_season:tf_season]
    weather_w = WeatherSeries([DailyWeather(T) for T in temps_w]; day_offset=t0_season)
    fly_w = make_olive_fly(n_adults=30.0)
    prob_w = PBDMProblem(DensityDependent(), fly_w, weather_w, (t0_season, tf_season))
    sol_w = solve(prob_w, DirectIteration();
                  reproduction_fn=olive_fly_reproduction)

    total_pop = [sum(sol_w.stage_totals[:, j]) for j in 1:size(sol_w.stage_totals, 2)]
    lines!(ax5, sol_w.t, total_pop, linewidth=2.5,
           color=warming_colors[idx], linestyle=warming_styles[idx],
           label=warming_labels[idx])

    pk = maximum(total_pop)
    pk_day = sol_w.t[argmax(total_pop)]
    println("  $(warming_labels[idx]): peak total=$(round(pk, digits=1)) at day $(pk_day)")
end

xlims!(ax5, t0_season, tf_season)
axislegend(ax5, position=:lt)

save(joinpath(figdir, "warming_comparison.png"), fig5, px_per_unit=2)
println("  Saved warming_comparison.png")

println("\n✓ All olive fly validation figures saved to: $(figdir)")
