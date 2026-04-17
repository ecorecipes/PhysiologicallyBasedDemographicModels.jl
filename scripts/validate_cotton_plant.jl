#!/usr/bin/env julia
# Validation script for Cotton Plant (Gossypium hirsutum L., cv. IAC-17) PBDM
# matching key biological parameters from:
#   - Gutierrez et al. (1984) — cotton PBDM with distributed delay
#   - Vignette 02_cotton_plant.qmd — parameter choices for this implementation
#
# Generates 5 figures in scripts/figures/cotton_plant/:
#   1. cotton_dev_rate.png          — development rate vs temperature (base 12°C)
#   2. cotton_carbon_allocation.png — supply-demand C allocation to organs
#   3. cotton_constant_temp.png     — organ mass at 25°C constant temp
#   4. cotton_seasonal.png          — seasonal phenology in California climate
#   5. cotton_yield.png             — boll number and mass predictions

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using PhysiologicallyBasedDemographicModels
using CairoMakie
using Statistics

figdir = joinpath(@__DIR__, "figures", "cotton_plant")
mkpath(figdir)

# ============================================================
# Parameters from Gutierrez et al. (1984) Table I and vignette
# ============================================================

# Temperature threshold for cotton development (degree-days above 12°C)
const T_THRESHOLD = 12.0   # °C (Gutierrez et al. 1975, 1984)
const T_UPPER     = 40.0   # °C upper limit

# Development rate model
dev_rate = LinearDevelopmentRate(T_THRESHOLD, T_UPPER)

# Key phenological benchmarks (degree-days above 12°C)
# From paper: first fruiting branch at 415 D°, peak squaring at 940 D°,
# first open boll at 1200 D°
const DD_FFB   = 415.0    # First fruiting branch (Table I, item 2)
const DD_PEAK  = 940.0    # Peak squaring (from field data)
const DD_BOLL  = 1200.0   # First open boll (from field data)

# Organ developmental times from vignette
const DD_LEAF_SENESCE = 700.0  # Leaf senescence age (paper: leaves die at 700 D°)
const DD_STEM  = 2000.0   # Stems: long-lived structural tissue
const DD_ROOT  = 150.0    # Root: rapid turnover
const DD_FRUIT = 800.0    # Fruit (square → open boll)

# Distributed delay substages
const K_SUBSTAGES = 30

# Respiration parameters (Q₁₀ model, vignette values)
leaf_resp  = Q10Respiration(0.030, 2.3, 25.0)   # 3.0% at 25°C
stem_resp  = Q10Respiration(0.015, 2.3, 25.0)   # 1.5% at 25°C
root_resp  = Q10Respiration(0.010, 2.3, 25.0)   # 1.0% at 25°C
fruit_resp = Q10Respiration(0.010, 2.3, 25.0)   # 1.0% at 25°C

# Light interception (Frazer-Gilbert functional response)
light_response = FraserGilbertResponse(0.7)

# Approach-aware hybrid layer used by the migrated vignette examples.
canopy_resp = Q10Respiration(0.01625, 2.3, 25.0)
cotton_bdf = BiodemographicFunctions(dev_rate, light_response, canopy_resp;
                                     label=:cotton_bdf)
cotton_mp = MetabolicPool(1.0,
    [0.8, 0.6, 0.4, 1.2],
    [:leaf, :stem, :root, :fruit])
cotton_hybrid = CoupledPBDMModel(cotton_bdf, cotton_mp; label=:cotton_hybrid)

# ============================================================
# Helper: build a fresh cotton Population
# ============================================================

function build_cotton(; W0_leaf=0.5, W0_stem=0.3, W0_root=0.2, W0_fruit=0.0)
    leaf_d  = DistributedDelay(K_SUBSTAGES, DD_LEAF_SENESCE; W0=W0_leaf)
    stem_d  = DistributedDelay(K_SUBSTAGES, DD_STEM;  W0=W0_stem)
    root_d  = DistributedDelay(K_SUBSTAGES, DD_ROOT;  W0=W0_root)
    fruit_d = DistributedDelay(K_SUBSTAGES, DD_FRUIT; W0=W0_fruit)

    Population(:cotton_IAC17, [
        LifeStage(:leaf,  leaf_d,  dev_rate, 0.001),
        LifeStage(:stem,  stem_d,  dev_rate, 0.0005),
        LifeStage(:root,  root_d,  dev_rate, 0.002),
        LifeStage(:fruit, fruit_d, dev_rate, 0.001),
    ])
end

# ============================================================
# Figure 1: Development rate vs temperature
# Paper states degree-day calculations above 12°C threshold
# ============================================================

println("=== Figure 1: Development rate vs temperature ===")

Ts = range(0.0, 45.0, length=400)
dd_linear = [degree_days(dev_rate, T) for T in Ts]

# Literature reference: simple degree-day model above 12°C
# At 20°C: DD=8, at 25°C: DD=13, at 30°C: DD=18, at 35°C: DD=23
lit_temps = [15.0, 20.0, 25.0, 30.0, 35.0]
lit_dd    = [3.0,  8.0,  13.0, 18.0, 23.0]  # T - 12°C

# Days to reach phenological milestones at different constant temps
println("Days to reach key phenological events at constant temperatures:")
for T in [20.0, 25.0, 30.0]
    dd_per_day = max(0.0, T - T_THRESHOLD)
    days_ffb  = DD_FFB / dd_per_day
    days_peak = DD_PEAK / dd_per_day
    days_boll = DD_BOLL / dd_per_day
    println("  T=$(T)°C: DD/day=$(dd_per_day), FFB=$(round(days_ffb, digits=0))d, " *
            "peak_sq=$(round(days_peak, digits=0))d, open_boll=$(round(days_boll, digits=0))d")
end

fig1 = Figure(size=(950, 650))
ax1 = Axis(fig1[1, 1],
    title="Cotton (IAC-17) — Development Rate vs Temperature\n" *
          "(Linear degree-day model, base temp = $(T_THRESHOLD)°C; Gutierrez et al. 1984)",
    xlabel="Temperature (°C)",
    ylabel="Degree-days per day (°C·d)",
    xlabelsize=14, ylabelsize=14)

lines!(ax1, collect(Ts), dd_linear, linewidth=2.5, color=:forestgreen,
       label="DD/day = max(0, T − $(Int(T_THRESHOLD)))")

scatter!(ax1, lit_temps, lit_dd, color=:red, markersize=12, marker=:diamond,
         label="Reference (T − 12°C)")

# Mark base temperature
vlines!(ax1, [T_THRESHOLD], color=:gray40, linestyle=:dash, linewidth=1.5)
text!(ax1, T_THRESHOLD + 0.5, maximum(dd_linear) * 0.92,
    text="T_base = $(T_THRESHOLD)°C", fontsize=11, color=:gray40,
    align=(:left, :top))

# Mark upper threshold
vlines!(ax1, [T_UPPER], color=:gray40, linestyle=:dash, linewidth=1.5)
text!(ax1, T_UPPER - 0.5, maximum(dd_linear) * 0.92,
    text="T_upper = $(Int(T_UPPER))°C", fontsize=11, color=:gray40,
    align=(:right, :top))

# Annotate phenological durations at 25°C
dd25 = 25.0 - T_THRESHOLD
text!(ax1, 26.0, dd25 + 1.0,
    text="At 25°C: $(round(dd25, digits=0)) DD/day\n" *
         "FFB = $(round(DD_FFB/dd25, digits=0))d\n" *
         "Peak sq = $(round(DD_PEAK/dd25, digits=0))d\n" *
         "Open boll = $(round(DD_BOLL/dd25, digits=0))d",
    fontsize=10, color=:forestgreen, align=(:left, :bottom))

xlims!(ax1, 0, 45)
ylims!(ax1, 0, nothing)
axislegend(ax1, position=:lt)

save(joinpath(figdir, "cotton_dev_rate.png"), fig1, px_per_unit=2)
println("Saved cotton_dev_rate.png")

# ============================================================
# Figure 2: Carbon allocation — supply/demand at different φ
# Paper: priority order is respiration > fruit > vegetative > reserves
# ============================================================

println("\n=== Figure 2: Carbon allocation (supply/demand) ===")

# Simulate allocation at a range of supply levels to show priority scheme
# Fix organ masses at mid-season values
leaf_mass  = 50.0   # g
stem_mass  = 60.0
root_mass  = 30.0
fruit_mass = 40.0
T_alloc = 25.0

# Compute demands
resp_demand_leaf  = respiration_rate(leaf_resp, T_alloc) * leaf_mass
resp_demand_stem  = respiration_rate(stem_resp, T_alloc) * stem_mass
resp_demand_root  = respiration_rate(root_resp, T_alloc) * root_mass
resp_demand_fruit = respiration_rate(fruit_resp, T_alloc) * fruit_mass
total_resp = resp_demand_leaf + resp_demand_stem + resp_demand_root + resp_demand_fruit

growth_demand_veg   = 0.02 * (leaf_mass + stem_mass + root_mass)
growth_demand_fruit = 0.03 * fruit_mass

total_demand = total_resp + growth_demand_fruit + growth_demand_veg

println("At 25°C with L=$(leaf_mass), S=$(stem_mass), R=$(root_mass), F=$(fruit_mass) g:")
println("  Respiration demand: $(round(total_resp, digits=2))")
println("  Fruit growth demand: $(round(growth_demand_fruit, digits=2))")
println("  Vegetative growth demand: $(round(growth_demand_veg, digits=2))")
println("  Total demand: $(round(total_demand, digits=2))")

# Sweep supply from 0 to 2× total demand
supply_range = range(0.0, 2.0 * total_demand, length=200)
alloc_resp   = Float64[]
alloc_fruit  = Float64[]
alloc_veg    = Float64[]
alloc_reserve = Float64[]
phi_vals     = Float64[]

for supply in supply_range
    pool = MetabolicPool(supply,
        [total_resp, growth_demand_fruit, growth_demand_veg],
        [:respiration, :fruit_growth, :vegetative_growth])
    alloc = allocate(pool)
    φ = supply_demand_index(pool)
    push!(alloc_resp, alloc[1])
    push!(alloc_fruit, alloc[2])
    push!(alloc_veg, alloc[3])
    push!(alloc_reserve, max(0.0, supply - sum(alloc)))
    push!(phi_vals, φ)
end

fig2 = Figure(size=(1000, 700))
ax2 = Axis(fig2[1, 1],
    title="Cotton (IAC-17) — Photosynthate Allocation by Priority\n" *
          "(Gutierrez et al. 1984: respiration > fruit > vegetative > reserves)",
    xlabel="Photosynthate supply (g/day)",
    ylabel="Allocation (g/day)",
    xlabelsize=14, ylabelsize=14)

# Stacked area bands
band!(ax2, collect(supply_range), zeros(length(supply_range)), alloc_resp,
      color=(:firebrick, 0.5), label="Respiration")
band!(ax2, collect(supply_range), alloc_resp, alloc_resp .+ alloc_fruit,
      color=(:goldenrod, 0.5), label="Fruit growth")
band!(ax2, collect(supply_range), alloc_resp .+ alloc_fruit,
      alloc_resp .+ alloc_fruit .+ alloc_veg,
      color=(:forestgreen, 0.5), label="Vegetative growth")
band!(ax2, collect(supply_range), alloc_resp .+ alloc_fruit .+ alloc_veg,
      alloc_resp .+ alloc_fruit .+ alloc_veg .+ alloc_reserve,
      color=(:steelblue, 0.3), label="Reserves")

# 1:1 line (supply = demand)
lines!(ax2, collect(supply_range), collect(supply_range),
       color=:black, linestyle=:dash, linewidth=1.5, label="Supply = Allocation")

# Mark total demand
vlines!(ax2, [total_demand], color=:red, linestyle=:dot, linewidth=2)
text!(ax2, total_demand + 0.1, maximum(supply_range) * 0.85,
    text="φ = 1.0\n(total demand)", fontsize=10, color=:red,
    align=(:left, :top))

# Mark stress threshold (paper: shedding when φ < 1)
idx_phi_half = findfirst(p -> p >= 0.5, phi_vals)
if idx_phi_half !== nothing
    supply_half = supply_range[idx_phi_half]
    vlines!(ax2, [supply_half], color=:darkorange, linestyle=:dot, linewidth=1.5)
    text!(ax2, supply_half + 0.1, maximum(supply_range) * 0.7,
        text="φ ≈ 0.5\n(fruit shed)", fontsize=10, color=:darkorange,
        align=(:left, :top))
end

xlims!(ax2, 0, 2.0 * total_demand)
ylims!(ax2, 0, 2.0 * total_demand)
axislegend(ax2, position=:lt)

save(joinpath(figdir, "cotton_carbon_allocation.png"), fig2, px_per_unit=2)
println("Saved cotton_carbon_allocation.png")

# ============================================================
# Figure 3: Organ mass dynamics at 25°C constant temperature
# Paper: validates allocation partitioning and leaf senescence at 700 D°
# ============================================================

println("\n=== Figure 3: Organ mass at constant 25°C ===")

n_days_const = 180
cotton_const = build_cotton()
weather_const = WeatherSeries(fill(25.0, n_days_const))
prob_const = PBDMProblem(cotton_hybrid, cotton_const, weather_const, (1, n_days_const))
sol_const = solve(prob_const, DirectIteration())

stage_names = [:leaf, :stem, :root, :fruit]
stage_colors = [:forestgreen, :saddlebrown, :darkorange, :firebrick]

# Compute cumulative DD
cdd_const = cumulative_degree_days(sol_const)
dd_per_day_25 = 25.0 - T_THRESHOLD  # 13 DD/day

println("At constant 25°C ($(dd_per_day_25) DD/day):")
for (i, sname) in enumerate(stage_names)
    traj = stage_trajectory(sol_const, i)
    peak = maximum(traj)
    peak_day = argmax(traj)
    println("  $(sname): peak=$(round(peak, digits=2)) at day $(peak_day) " *
            "($(round(cdd_const[min(peak_day, length(cdd_const))], digits=0)) DD)")
end

# Literature reference points — expected DD milestones mapped to calendar days at 25°C
lit_days_ffb  = DD_FFB / dd_per_day_25    # ~32 days
lit_days_peak = DD_PEAK / dd_per_day_25   # ~72 days
lit_days_boll = DD_BOLL / dd_per_day_25   # ~92 days
lit_days_senesce = DD_LEAF_SENESCE / dd_per_day_25  # ~54 days

fig3 = Figure(size=(1000, 700))
ax3 = Axis(fig3[1, 1],
    title="Cotton (IAC-17) — Organ Mass Dynamics at Constant 25°C\n" *
          "($(dd_per_day_25) DD/day; dashed lines = phenological milestones from paper)",
    xlabel="Day",
    ylabel="Organ mass (g/plant)",
    xlabelsize=14, ylabelsize=14)

for (i, sname) in enumerate(stage_names)
    traj = stage_trajectory(sol_const, i)
    lines!(ax3, sol_const.t, traj, linewidth=2.5, color=stage_colors[i],
           label=String(sname))
end

# Overlay phenological milestone markers from literature
milestone_days = [lit_days_ffb, lit_days_peak, lit_days_boll, lit_days_senesce]
milestone_labels = ["FFB\n(415 DD)", "Peak sq.\n(940 DD)", "Open boll\n(1200 DD)",
                    "Leaf sen.\n(700 DD)"]
milestone_colors = [:purple, :blue, :red, :forestgreen]

for (md, ml, mc) in zip(milestone_days, milestone_labels, milestone_colors)
    if md <= n_days_const
        vlines!(ax3, [md], color=(mc, 0.5), linestyle=:dash, linewidth=1.5)
        text!(ax3, md + 1, maximum(stage_trajectory(sol_const, 1)) * 0.95,
            text=ml, fontsize=9, color=mc, align=(:left, :top))
    end
end

xlims!(ax3, 1, n_days_const)
ylims!(ax3, 0, nothing)
axislegend(ax3, position=:rt)

save(joinpath(figdir, "cotton_constant_temp.png"), fig3, px_per_unit=2)
println("Saved cotton_constant_temp.png")

# ============================================================
# Figure 4: Seasonal phenology in California San Joaquin climate
# Paper compares Londrina, Brazil; here we show a California analog
# as noted in paper: "very similar to Acala cotton in California"
# ============================================================

println("\n=== Figure 4: Seasonal phenology (California San Joaquin Valley) ===")

# San Joaquin Valley: planting ~April 15 (day 105), harvest ~October
# Approximate seasonal temperature: warm Mediterranean summer
n_days_season = 210
temps_sjv = Float64[]
rads_sjv  = Float64[]
for d in 1:n_days_season
    day_of_year = 104 + d  # Start April 15
    # San Joaquin Valley seasonal cycle: warm summers, mild winters
    T = 20.0 + 10.0 * sin(2π * (day_of_year - 80) / 365)
    push!(temps_sjv, T)
    # Solar radiation (MJ/m²/day)
    R = 20.0 + 8.0 * sin(2π * (day_of_year - 80) / 365)
    push!(rads_sjv, R)
end

weather_sjv_days = [DailyWeather(temps_sjv[d], temps_sjv[d] - 5, temps_sjv[d] + 5;
                                 radiation=rads_sjv[d]) for d in 1:n_days_season]
weather_sjv = WeatherSeries(weather_sjv_days; day_offset=1)

cotton_sjv = build_cotton()
prob_sjv = PBDMProblem(cotton_hybrid, cotton_sjv, weather_sjv, (1, n_days_season))
sol_sjv = solve(prob_sjv, DirectIteration())
cdd_sjv = cumulative_degree_days(sol_sjv)

println("San Joaquin Valley simulation ($(n_days_season) days):")
println("  Total DD accumulated: $(round(cdd_sjv[end], digits=0))")

# Find phenological milestone days
for (name, threshold) in [("First fruiting branch", DD_FFB),
                          ("Peak squaring", DD_PEAK),
                          ("First open boll", DD_BOLL)]
    idx = findfirst(c -> c >= threshold, cdd_sjv)
    if idx !== nothing
        println("  $name: day $(sol_sjv.t[idx]) ($(round(cdd_sjv[idx], digits=0)) DD)")
    else
        println("  $name: NOT reached (max DD = $(round(cdd_sjv[end], digits=0)))")
    end
end

# Also simulate Londrina, Brazil (the paper's field site)
n_days_lond = 210
temps_lond = Float64[]
rads_lond  = Float64[]
for d in 1:n_days_lond
    day_of_year = mod(297 + d, 365) + 1  # Planting Oct 25
    T = 22.0 + 5.0 * sin(2π * (day_of_year - 355) / 365)
    push!(temps_lond, T)
    R = 18.0 + 6.0 * sin(2π * (day_of_year - 355) / 365)
    push!(rads_lond, R)
end

weather_lond_days = [DailyWeather(temps_lond[d], temps_lond[d] - 4, temps_lond[d] + 4;
                                  radiation=rads_lond[d]) for d in 1:n_days_lond]
weather_lond = WeatherSeries(weather_lond_days; day_offset=1)

cotton_lond = build_cotton()
prob_lond = PBDMProblem(cotton_hybrid, cotton_lond, weather_lond, (1, n_days_lond))
sol_lond = solve(prob_lond, DirectIteration())
cdd_lond = cumulative_degree_days(sol_lond)

println("Londrina, Brazil simulation:")
println("  Total DD accumulated: $(round(cdd_lond[end], digits=0))")

fig4 = Figure(size=(1100, 800))

# Top: cumulative DD comparison
ax4a = Axis(fig4[1, 1],
    title="Cotton (IAC-17) — Seasonal Phenology Comparison\n" *
          "(Paper: Londrina, Brazil 1982–83; Analog: San Joaquin Valley, CA)",
    ylabel="Cumulative DD (°C·d)",
    xlabelsize=14, ylabelsize=14)

lines!(ax4a, sol_sjv.t, cdd_sjv, linewidth=2.5, color=:steelblue,
       label="San Joaquin Valley, CA")
lines!(ax4a, sol_lond.t, cdd_lond, linewidth=2.5, color=:firebrick,
       label="Londrina, Brazil")

# Horizontal lines for paper milestones
for (dd_val, label, col) in [(DD_FFB, "FFB (415 DD)", :purple),
                              (DD_PEAK, "Peak sq. (940 DD)", :blue),
                              (DD_BOLL, "Open boll (1200 DD)", :red)]
    hlines!(ax4a, [dd_val], color=(col, 0.4), linestyle=:dash, linewidth=1.5)
    text!(ax4a, 5, dd_val + 20, text=label, fontsize=9, color=col,
          align=(:left, :bottom))
end

axislegend(ax4a, position=:lt)
xlims!(ax4a, 1, n_days_season)
hidexdecorations!(ax4a, grid=false)

# Bottom: organ dynamics for SJV
ax4b = Axis(fig4[2, 1],
    xlabel="Day after planting",
    ylabel="Organ mass (g/plant)",
    xlabelsize=14, ylabelsize=14)

for (i, sname) in enumerate(stage_names)
    traj_sjv = stage_trajectory(sol_sjv, i)
    lines!(ax4b, sol_sjv.t, traj_sjv, linewidth=2, color=stage_colors[i],
           label="$(sname) (SJV)")
end
for (i, sname) in enumerate(stage_names)
    traj_lond = stage_trajectory(sol_lond, i)
    lines!(ax4b, sol_lond.t, traj_lond, linewidth=1.5, color=stage_colors[i],
           linestyle=:dash, label="$(sname) (Londrina)")
end

axislegend(ax4b, position=:rt, nbanks=2)
xlims!(ax4b, 1, n_days_season)
ylims!(ax4b, 0, nothing)
linkxaxes!(ax4a, ax4b)
rowsize!(fig4.layout, 1, Relative(0.4))

save(joinpath(figdir, "cotton_seasonal.png"), fig4, px_per_unit=2)
println("Saved cotton_seasonal.png")

# ============================================================
# Figure 5: Yield predictions — temperature sensitivity
# Paper: yields at low/standard/high T = 1.16, 1.04, 0.52 t/ha
# Paper: minimum solar radiation 650 cal/cm²/day needed for max yield
# ============================================================

println("\n=== Figure 5: Yield predictions (temperature sensitivity) ===")

# Paper simulation results (Table from p.244):
# Low temp (-10%): yield ~1.16 t/ha, bolls 10% higher
# Standard: yield ~1.04 t/ha
# High temp (+10%): yield ~0.52 t/ha, bolls 10-15% lower
lit_labels = ["-10% Temp", "Standard", "+10% Temp"]
lit_yields = [1.16, 1.04, 0.52]  # t/ha from paper

# Run three temperature scenarios using Londrina weather (as in paper)
offsets = [-2.2, 0.0, +2.2]  # Approximate ±10% of ~22°C mean
scenario_colors = [:steelblue, :black, :firebrick]

sim_results = []
for (label, offset, col) in zip(lit_labels, offsets, scenario_colors)
    mod_temps = temps_lond .+ offset
    mod_days = [DailyWeather(mod_temps[d], mod_temps[d] - 4, mod_temps[d] + 4;
                             radiation=rads_lond[d]) for d in 1:n_days_lond]
    mod_weather = WeatherSeries(mod_days; day_offset=1)

    cotton_mod = build_cotton()
    prob_mod = PBDMProblem(cotton_hybrid, cotton_mod, mod_weather, (1, n_days_lond))
    sol_mod = solve(prob_mod, DirectIteration())
    cdd_mod = cumulative_degree_days(sol_mod)
    fruit_traj = stage_trajectory(sol_mod, 4)

    push!(sim_results, (label=label, sol=sol_mod, cdd=cdd_mod,
                        fruit=fruit_traj, color=col))

    peak_fruit = maximum(fruit_traj)
    final_fruit = fruit_traj[end]
    total_dd = cdd_mod[end]
    println("  $label (offset=$(offset > 0 ? "+" : "")$(offset)°C): " *
            "total_DD=$(round(total_dd, digits=0)), peak_fruit=$(round(peak_fruit, digits=2)), " *
            "final_fruit=$(round(final_fruit, digits=2))")
end

fig5 = Figure(size=(1100, 850))

# Top-left: fruit mass trajectories
ax5a = Axis(fig5[1, 1],
    title="Cotton — Fruit Mass Dynamics\n(Londrina weather ± temperature)",
    xlabel="Day after planting",
    ylabel="Fruit mass (g/plant)",
    xlabelsize=13, ylabelsize=13)

for r in sim_results
    lines!(ax5a, r.sol.t, r.fruit, linewidth=2.5, color=r.color,
           label=r.label)
end

axislegend(ax5a, position=:lt)
xlims!(ax5a, 1, n_days_lond)
ylims!(ax5a, 0, nothing)

# Top-right: total vegetative mass
ax5b = Axis(fig5[1, 2],
    title="Cotton — Total Vegetative Mass\n(leaf + stem + root)",
    xlabel="Day after planting",
    ylabel="Vegetative mass (g/plant)",
    xlabelsize=13, ylabelsize=13)

for r in sim_results
    veg_traj = stage_trajectory(r.sol, 1) .+ stage_trajectory(r.sol, 2) .+
               stage_trajectory(r.sol, 3)
    lines!(ax5b, r.sol.t, veg_traj, linewidth=2.5, color=r.color,
           label=r.label)
end

axislegend(ax5b, position=:lt)
xlims!(ax5b, 1, n_days_lond)
ylims!(ax5b, 0, nothing)

# Bottom-left: simulated yield comparison bar chart
ax5c = Axis(fig5[2, 1],
    title="Simulated vs Literature Yields\n(Gutierrez et al. 1984: 1.16, 1.04, 0.52 t/ha)",
    xlabel="Temperature scenario",
    ylabel="Relative fruit mass (normalized)",
    xlabelsize=13, ylabelsize=13,
    xticks=(1:3, lit_labels))

# Normalize simulated final fruit mass to standard run for comparison
sim_final_fruits = [r.fruit[end] for r in sim_results]
std_fruit = sim_final_fruits[2]  # standard run
sim_normalized = std_fruit > 0 ? sim_final_fruits ./ std_fruit : sim_final_fruits
lit_normalized = lit_yields ./ lit_yields[2]

barplot!(ax5c, [0.8, 1.8, 2.8], sim_normalized,
         color=scenario_colors, width=0.35, label="Simulated")
barplot!(ax5c, [1.2, 2.2, 3.2], lit_normalized,
         color=[(:steelblue, 0.3), (:gray, 0.3), (:firebrick, 0.3)],
         width=0.35, label="Literature")

# Add text labels
for (i, (sn, ln)) in enumerate(zip(sim_normalized, lit_normalized))
    text!(ax5c, i - 0.2, sn + 0.03,
        text="$(round(sn, digits=2))", fontsize=9, color=:black,
        align=(:center, :bottom))
    text!(ax5c, i + 0.2, ln + 0.03,
        text="$(round(ln, digits=2))", fontsize=9, color=:gray50,
        align=(:center, :bottom))
end

axislegend(ax5c, position=:rt)

# Bottom-right: cumulative DD comparison
ax5d = Axis(fig5[2, 2],
    title="Cumulative Degree-Days\n(determines phenological timing)",
    xlabel="Day after planting",
    ylabel="Cumulative DD (°C·d)",
    xlabelsize=13, ylabelsize=13)

for r in sim_results
    lines!(ax5d, r.sol.t, r.cdd, linewidth=2.5, color=r.color,
           label=r.label)
end

# Reference lines from paper
hlines!(ax5d, [DD_FFB], color=(:purple, 0.3), linestyle=:dash, linewidth=1)
text!(ax5d, 5, DD_FFB + 15, text="FFB (415 DD)", fontsize=8, color=:purple,
      align=(:left, :bottom))
hlines!(ax5d, [DD_PEAK], color=(:blue, 0.3), linestyle=:dash, linewidth=1)
text!(ax5d, 5, DD_PEAK + 15, text="Peak sq. (940 DD)", fontsize=8, color=:blue,
      align=(:left, :bottom))
hlines!(ax5d, [DD_BOLL], color=(:red, 0.3), linestyle=:dash, linewidth=1)
text!(ax5d, 5, DD_BOLL + 15, text="Open boll (1200 DD)", fontsize=8, color=:red,
      align=(:left, :bottom))

axislegend(ax5d, position=:lt)
xlims!(ax5d, 1, n_days_lond)

save(joinpath(figdir, "cotton_yield.png"), fig5, px_per_unit=2)
println("Saved cotton_yield.png")

# ============================================================
# Summary
# ============================================================

println("\n=== Summary ===")
println("Paper reference values (Gutierrez et al. 1984):")
println("  Base temperature: 12°C ✓")
println("  First fruiting branch: 415 DD (8.5 nodes)")
println("  Peak squaring: ~940 DD")
println("  First open boll: ~1200 DD")
println("  Leaf senescence age: 700 DD")
println("  Yield sensitivity: -10%T → +12% yield, +10%T → -50% yield")
println("  Minimum radiation for max yield: 650 cal/cm²/day (~315 W/m²)")
println("\nAll cotton plant validation figures saved to: $(figdir)")
