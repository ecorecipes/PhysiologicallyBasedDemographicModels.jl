#!/usr/bin/env julia
# Validate Cowpea–Thrips PBDM for Megalurothrips sjostedti (bean flower thrips).
# Run: julia --project=. scripts/validate_cowpea_thrips.jl

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using PhysiologicallyBasedDemographicModels
using CairoMakie

figdir = joinpath(@__DIR__, "figures", "cowpea_thrips")
mkpath(figdir)

println("="^65)
println("Cowpea–Thrips (Megalurothrips sjostedti) PBDM Validation")
println("="^65)
println()

# ============================================================
# Parameters — estimated from related thrips species
# ============================================================

# Developmental thresholds (estimated from Frankliniella, Thrips tabaci)
const T_LOWER = 10.0   # lower developmental threshold (°C)
const T_UPPER = 35.0   # upper developmental threshold (°C)

# Degree-day requirements (estimated)
const DD_EGG   = 50.0   # ~3–4 days at 25°C
const DD_LARVA = 80.0   # ~5–6 days at 25°C
const DD_PUPA  = 40.0   # ~2–3 days at 25°C
const DD_ADULT = 80.0   # adult longevity ~5–6 days at 25°C
const DD_TOTAL = DD_EGG + DD_LARVA + DD_PUPA  # egg-to-adult ~170 DD

# Brière development rate coefficients
const egg_dev   = BriereDevelopmentRate(3.0e-4, T_LOWER, T_UPPER)
const larva_dev = BriereDevelopmentRate(1.8e-4, T_LOWER, T_UPPER)
const pupa_dev  = BriereDevelopmentRate(2.4e-4, T_LOWER, T_UPPER)

# Fecundity: eggs per adult per day (temperature-dependent)
const FECUNDITY_MAX = 4.0  # peak eggs/adult/day near optimum
function thrips_fecundity(T::Real)
    T <= T_LOWER && return 0.0
    T >= T_UPPER && return 0.0
    # Bell-shaped curve peaking around 28°C
    T_opt = 28.0
    σ = 6.0
    return FECUNDITY_MAX * exp(-0.5 * ((T - T_opt) / σ)^2)
end

# Reproduction function for DensityDependent solver
function thrips_reproduction(pop, w, p, day)
    n_adult = delay_total(pop.stages[4].delay)
    return max(0.0, n_adult * thrips_fecundity(w.T_mean))
end

# Cowpea crop parameters
const COWPEA_BASE_TEMP   = 10.0    # base temperature (°C)
const COWPEA_DD_FLOWER   = 600.0   # DD to onset of flowering
const COWPEA_DD_POD      = 800.0   # DD to onset of podding
const COWPEA_DD_MATURITY = 1200.0  # DD to crop maturity
const POTENTIAL_YIELD    = 1500.0  # kg/ha

# Sahel climate parameters
const SAHEL_T_MEAN  = 28.0   # mean annual temperature (°C)
const SAHEL_T_AMP   = 5.0    # annual amplitude (°C)
const SAHEL_T_PHASE = 210.0  # day of peak temperature

println("Parameters:")
println("  Thrips T_lower=$(T_LOWER)°C, T_upper=$(T_UPPER)°C")
println("  DD: egg=$DD_EGG, larva=$DD_LARVA, pupa=$DD_PUPA, adult=$DD_ADULT")
println("  DD egg-to-adult: $DD_TOTAL")
println("  Fecundity max: $(FECUNDITY_MAX) eggs/adult/day at 28°C")
println("  Cowpea base=$(COWPEA_BASE_TEMP)°C, DD to flower=$(COWPEA_DD_FLOWER)")
println("  Potential yield: $(POTENTIAL_YIELD) kg/ha")
println()

# ============================================================
# Helper: Sahel temperature for a given day of year
# ============================================================

function sahel_temperature(day::Int)
    return SAHEL_T_MEAN + SAHEL_T_AMP * sin(2π * (day - SAHEL_T_PHASE) / 365)
end

# ============================================================
# Helper: Build thrips population with given adult initial count
# ============================================================

function make_thrips_population(; adult_W0=5.0)
    # W0 is per-substage; total adults = k_adult * adult_W0
    stages = [
        LifeStage(:egg,   DistributedDelay(15, DD_EGG;   W0=0.0),      egg_dev,   0.01),
        LifeStage(:larva, DistributedDelay(20, DD_LARVA;  W0=0.0),      larva_dev, 0.008),
        LifeStage(:pupa,  DistributedDelay(10, DD_PUPA;   W0=0.0),      pupa_dev,  0.005),
        LifeStage(:adult, DistributedDelay(10, DD_ADULT;  W0=adult_W0), egg_dev,   0.015),
    ]
    return Population(:megalurothrips_sjostedti, stages)
end

# ============================================================
# Helper: Simple cowpea phenology (logistic biomass curves)
# ============================================================

function cowpea_phenology(planting_day::Int, n_days::Int)
    leaf_biomass   = zeros(n_days)
    flower_biomass = zeros(n_days)
    pod_biomass    = zeros(n_days)
    dd_cum = 0.0
    for d in 1:n_days
        day_of_year = planting_day + d - 1
        T = sahel_temperature(day_of_year)
        dd_cum += max(0.0, T - COWPEA_BASE_TEMP)
        if dd_cum < COWPEA_DD_MATURITY
            leaf_frac = min(1.0, dd_cum / COWPEA_DD_FLOWER)
            senesce   = dd_cum > COWPEA_DD_POD ? (dd_cum - COWPEA_DD_POD) / (COWPEA_DD_MATURITY - COWPEA_DD_POD) : 0.0
            leaf_biomass[d] = leaf_frac * (1.0 - 0.6 * senesce)
        end
        if COWPEA_DD_FLOWER <= dd_cum < COWPEA_DD_MATURITY
            flower_progress = (dd_cum - COWPEA_DD_FLOWER) / (COWPEA_DD_POD - COWPEA_DD_FLOWER)
            flower_biomass[d] = flower_progress < 1.0 ? sin(π * flower_progress) : max(0.0, 1.0 - (dd_cum - COWPEA_DD_POD) / 200.0)
        end
        if dd_cum >= COWPEA_DD_POD
            pod_progress = min(1.0, (dd_cum - COWPEA_DD_POD) / (COWPEA_DD_MATURITY - COWPEA_DD_POD))
            pod_biomass[d] = pod_progress
        end
    end
    return leaf_biomass, flower_biomass, pod_biomass
end

# ============================================================
# Figure 1: Development rate curves
# ============================================================

println("--- Figure 1: Development rate curves ---")

Ts = range(0.0, 40.0, length=300)
rates_egg   = [development_rate(egg_dev, T)   for T in Ts]
rates_larva = [development_rate(larva_dev, T) for T in Ts]
rates_pupa  = [development_rate(pupa_dev, T)  for T in Ts]

# F. occidentalis reference data (approximate from literature)
ref_temps_egg   = [20.0, 25.0]
ref_rates_egg   = [0.10, 0.17]
ref_temps_larva = [20.0, 25.0]
ref_rates_larva = [0.04, 0.08]

fig1 = Figure(size=(800, 600))
ax1 = Axis(fig1[1, 1],
           xlabel="Temperature (°C)", ylabel="Development rate (day⁻¹)",
           title="M. sjostedti — Stage-Specific Development Rates (Brière Model)",
           xlabelsize=14, ylabelsize=14)
lines!(ax1, collect(Ts), rates_egg,   color=:gold,  linewidth=2.5, label="Egg")
lines!(ax1, collect(Ts), rates_larva, color=:green, linewidth=2.5, label="Larva")
lines!(ax1, collect(Ts), rates_pupa,  color=:brown, linewidth=2.5, label="Pupa")
scatter!(ax1, ref_temps_egg, ref_rates_egg, color=:gold, marker=:diamond,
         markersize=12, label="F. occidentalis egg (ref)")
scatter!(ax1, ref_temps_larva, ref_rates_larva, color=:green, marker=:utriangle,
         markersize=12, label="F. occidentalis larva (ref)")
vlines!(ax1, [T_LOWER], color=:blue, linestyle=:dot, linewidth=1,
        label="T_lower = $(T_LOWER)°C")
vlines!(ax1, [T_UPPER], color=:red, linestyle=:dot, linewidth=1,
        label="T_upper = $(T_UPPER)°C")
xlims!(ax1, 0, 40)
ylims!(ax1, 0, maximum(rates_egg) * 1.2)
axislegend(ax1, position=:lt, framevisible=false)

save(joinpath(figdir, "devrate_curves.png"), fig1, px_per_unit=2)
println("  Saved devrate_curves.png")

println("  Development rate comparison:")
println("  T(°C)  | Egg r(T) | Larva r(T) | Pupa r(T)")
for T in [15.0, 20.0, 25.0, 28.0, 30.0, 33.0]
    r_e = development_rate(egg_dev, T)
    r_l = development_rate(larva_dev, T)
    r_p = development_rate(pupa_dev, T)
    println("  $(lpad(T, 5)) |  $(lpad(round(r_e, digits=4), 7)) |   $(lpad(round(r_l, digits=4), 7))  |  $(lpad(round(r_p, digits=4), 7))")
end
println()

# ============================================================
# Figure 2: Cowpea phenology for June planting
# ============================================================

println("--- Figure 2: Cowpea phenology (June planting) ---")

planting_june = 152  # ~June 1
n_grow = 120         # 120-day growing season
leaf_b, flower_b, pod_b = cowpea_phenology(planting_june, n_grow)
grow_days = planting_june .+ (0:n_grow-1)

fig2 = Figure(size=(800, 600))
ax2 = Axis(fig2[1, 1],
           xlabel="Day of year", ylabel="Relative biomass",
           title="Cowpea Phenology — June Planting (Day $planting_june), Sahel Climate",
           xlabelsize=14, ylabelsize=14)
lines!(ax2, grow_days, leaf_b,   color=:green,  linewidth=2.5, label="Leaf")
lines!(ax2, grow_days, flower_b, color=:orange, linewidth=2.5, label="Flower")
lines!(ax2, grow_days, pod_b,    color=:brown,  linewidth=2.5, label="Pod")
vlines!(ax2, [planting_june], color=:gray, linestyle=:dash, linewidth=1, label="Planting")
xlims!(ax2, planting_june - 5, planting_june + n_grow + 5)
ylims!(ax2, 0, 1.15)
axislegend(ax2, position=:rt, framevisible=false)

save(joinpath(figdir, "cowpea_phenology.png"), fig2, px_per_unit=2)
println("  Saved cowpea_phenology.png")
println()

# ============================================================
# Figure 3: Thrips cohort at constant 28°C with reproduction
# ============================================================

println("--- Figure 3: Thrips simulation at constant 28°C ---")

n_sim = 150
pop_28 = make_thrips_population(adult_W0=5.0)  # 10 substages × 5 = 50 adults
ws_28  = WeatherSeries([DailyWeather(28.0) for _ in 1:n_sim]; day_offset=1)
prob_28 = PBDMProblem(DensityDependent(), pop_28, ws_28, (1, n_sim))
sol_28  = solve(prob_28, DirectIteration(); reproduction_fn=thrips_reproduction)

egg_28   = stage_trajectory(sol_28, 1)
larva_28 = stage_trajectory(sol_28, 2)
pupa_28  = stage_trajectory(sol_28, 3)
adult_28 = stage_trajectory(sol_28, 4)
tp_28    = total_population(sol_28)
cdd_28   = cumulative_degree_days(sol_28)

println("  DD per day at 28°C (egg model): $(round(development_rate(egg_dev, 28.0), digits=4))")
println("  Fecundity at 28°C: $(round(thrips_fecundity(28.0), digits=2)) eggs/adult/day")
println("  Initial population:  $(round(tp_28[1], digits=1))")
println("  Peak population:     $(round(maximum(tp_28), digits=1)) at day $(sol_28.t[argmax(tp_28)])")
println("  Final population:    $(round(tp_28[end], digits=1))")
println("  Total DD accumulated: $(round(cdd_28[end], digits=0))")

fig3 = Figure(size=(900, 600))
ax3 = Axis(fig3[1, 1],
           xlabel="Day", ylabel="Population (individuals)",
           title="M. sjostedti at Constant 28°C — Stage Dynamics (50 adults initial)",
           xlabelsize=14, ylabelsize=14)
lines!(ax3, sol_28.t, egg_28,   label="Egg",   color=:gold,  linewidth=1.5)
lines!(ax3, sol_28.t, larva_28, label="Larva", color=:green, linewidth=1.5)
lines!(ax3, sol_28.t, pupa_28,  label="Pupa",  color=:brown, linewidth=1.5)
lines!(ax3, sol_28.t, adult_28, label="Adult", color=:red,   linewidth=2)
lines!(ax3, sol_28.t, tp_28,    label="Total", color=:black, linewidth=2, linestyle=:dash)
axislegend(ax3, position=:lt, framevisible=false)

save(joinpath(figdir, "sim_constant_28C.png"), fig3, px_per_unit=2)
println("  Saved sim_constant_28C.png")
println()

# ============================================================
# Figure 4: Planting date comparison — total thrips population
# ============================================================

println("--- Figure 4: Planting date comparison ---")

# Background thrips increase through the warm season: later planting
# encounters higher initial thrips populations.
planting_specs = [
    ("May 15",  135, 2.0),   # early season → low background thrips
    ("June 15", 166, 5.0),   # mid season → moderate
    ("July 15", 196, 12.0),  # peak season → high background
]
planting_order  = [s[1] for s in planting_specs]
planting_colors = [:blue, :orange, :red]
n_season = 120

sim_results = Dict{String, Any}()

for (label, pday, w0) in planting_specs
    temps = [sahel_temperature(pday + d - 1) for d in 1:n_season]
    ws = WeatherSeries([DailyWeather(T) for T in temps]; day_offset=pday)

    pop = make_thrips_population(adult_W0=w0)
    prob = PBDMProblem(DensityDependent(), pop, ws, (pday, pday + n_season))
    sol = solve(prob, DirectIteration(); reproduction_fn=thrips_reproduction)

    tp = total_population(sol)
    # Compute cumulative thrips-days during the flowering window (days 30–60)
    flower_start = 30
    flower_end   = min(60, length(tp) - 1)
    thrips_days_flower = sum(tp[flower_start+1:flower_end+1])
    peak = maximum(tp)
    peak_day = sol.t[argmax(tp)]

    sim_results[label] = (sol=sol, tp=tp, peak=peak, peak_day=peak_day,
                          thrips_days=thrips_days_flower, w0=w0)
    println("  $label (day $pday, W0=$w0): peak=$(round(peak, digits=1)) at day $peak_day, " *
            "flower thrips-days=$(round(thrips_days_flower, digits=0))")
end

fig4 = Figure(size=(900, 600))
ax4 = Axis(fig4[1, 1],
           xlabel="Day of year", ylabel="Total thrips population",
           title="M. sjostedti Population by Planting Date — Sahel Climate",
           xlabelsize=14, ylabelsize=14)
for (i, label) in enumerate(planting_order)
    res = sim_results[label]
    lines!(ax4, res.sol.t, res.tp, color=planting_colors[i], linewidth=2.5,
           label=label)
end
axislegend(ax4, position=:lt, framevisible=false)

save(joinpath(figdir, "planting_date_comparison.png"), fig4, px_per_unit=2)
println("  Saved planting_date_comparison.png")
println()

# ============================================================
# Figure 5: Yield loss bar chart
# ============================================================

println("--- Figure 5: Yield loss by planting date ---")

# Yield loss based on cumulative thrips-days during the flowering window.
# Exponential damage model: yield = potential × exp(-α × thrips_days_flower)
const DAMAGE_ALPHA = 5.0e-6

yields     = Float64[]
bar_labels = String[]

for label in planting_order
    res = sim_results[label]
    damage_fraction = 1.0 - exp(-DAMAGE_ALPHA * res.thrips_days)
    y = POTENTIAL_YIELD * (1.0 - damage_fraction)
    push!(yields, y)
    push!(bar_labels, label)
    println("  $label: thrips-days=$(round(res.thrips_days, digits=0)), " *
            "damage=$(round(damage_fraction*100, digits=1))%, yield=$(round(y, digits=0)) kg/ha")
end

fig5 = Figure(size=(700, 500))
ax5 = Axis(fig5[1, 1],
           xlabel="Planting date", ylabel="Cowpea yield (kg/ha)",
           title="Estimated Cowpea Yield by Planting Date",
           xlabelsize=14, ylabelsize=14,
           xticks=(1:3, bar_labels))
barplot!(ax5, 1:3, yields, color=[:blue, :orange, :red],
         strokewidth=1, strokecolor=:black)
hlines!(ax5, [POTENTIAL_YIELD], color=:gray, linestyle=:dash, linewidth=1,
        label="Potential yield ($POTENTIAL_YIELD kg/ha)")
ylims!(ax5, 0, POTENTIAL_YIELD * 1.15)
axislegend(ax5, position=:rt, framevisible=false)

save(joinpath(figdir, "yield_loss_bar.png"), fig5, px_per_unit=2)
println("  Saved yield_loss_bar.png")

println()
println("="^65)
println("All cowpea-thrips validation figures saved to: $figdir")
println("="^65)
