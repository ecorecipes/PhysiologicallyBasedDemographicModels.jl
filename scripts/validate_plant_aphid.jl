#!/usr/bin/env julia
# Validation script for the Plant–Aphid–Parasitoid Tritrophic Dynamics vignette.
#
# Generates five figures in scripts/figures/plant_aphid/:
#   1. devrate_temperature.png     — Aphid and parasitoid development rate curves (0–35°C)
#   2. plant_phenology.png         — Thimbleberry leaf area index over the growing season
#   3. aphid_no_parasitoid.png     — Aphid population dynamics at constant 18°C
#   4. tritrophic_dynamics.png     — Aphid with parasitoid at Berkeley, CA climate
#   5. supply_demand_ratios.png    — Phi values at plant–aphid and aphid–parasitoid interfaces
#
# Literature: Gilbert & Gutierrez (1973) J. Anim. Ecol. 42:323–340

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using PhysiologicallyBasedDemographicModels
using CairoMakie

figdir = joinpath(@__DIR__, "figures", "plant_aphid")
mkpath(figdir)

# ============================================================
# Species parameters — from vignette & Gilbert & Gutierrez (1973)
# ============================================================

# --- Thermal thresholds ---
const APHID_T_BASE = 4.4    # °C — lower developmental threshold (temperate aphid)
const APHID_T_MAX  = 32.0   # °C — upper lethal threshold
const PARA_T_BASE  = 6.0    # °C — parasitoid lower threshold (slightly higher)
const PARA_T_MAX   = 33.0   # °C — parasitoid upper threshold

# --- Development rate models (Brière) ---
# Aphid nymph: ~110 DD above 4.4°C for development; r_max ≈ 0.15–0.30 /day at ~20°C
aphid_dev = BriereDevelopmentRate(0.00012, APHID_T_BASE, APHID_T_MAX)
# Parasitoid larva: ~160 DD above 6.0°C; slower development
para_dev  = BriereDevelopmentRate(0.00006, PARA_T_BASE, PARA_T_MAX)

# --- Degree-day durations ---
const DD_NYMPH = 3.0     # nymph τ (physiological time units)
const DD_ADULT = 2.0     # adult τ
const DD_PARA_LARVA = 4.0  # parasitoid larval τ
const DD_PARA_ADULT = 2.5  # parasitoid adult τ

# --- Mortality rates (per degree-day) ---
const μ_NYMPH = 0.03
const μ_ADULT = 0.02
const μ_PARA_LARVA = 0.04
const μ_PARA_ADULT = 0.03

# --- Thimbleberry leaf phenology ---
const LEAF_FLUSH_DAY = 60     # start of leaf expansion (~March 1)
const PEAK_LEAF_DAY  = 170    # peak leaf area (~mid-June)
const SENESCENCE_DAY = 220    # onset of senescence (~early August)
const DORMANCY_DAY   = 270    # full dormancy (~late September)
const K_MAX = 8000.0          # peak carrying capacity (aphids)

function leaf_capacity(day)
    if day < LEAF_FLUSH_DAY || day > DORMANCY_DAY
        return 0.0
    elseif day <= PEAK_LEAF_DAY
        t = (day - LEAF_FLUSH_DAY) / (PEAK_LEAF_DAY - LEAF_FLUSH_DAY)
        return K_MAX / (1.0 + exp(-10.0 * (t - 0.5)))
    elseif day <= SENESCENCE_DAY
        return K_MAX
    else
        t = (day - SENESCENCE_DAY) / (DORMANCY_DAY - SENESCENCE_DAY)
        return K_MAX * (1.0 - t)
    end
end

# Leaf area index (normalized 0–1 for phenology plot)
function leaf_area_index(day)
    return leaf_capacity(day) / K_MAX
end

# --- Berkeley, CA weather ---
const N_DAYS = 365
function make_berkeley_weather()
    days = DailyWeather{Float64}[]
    for d in 1:N_DAYS
        T_mean = 13.5 + 4.5 * sin(2π * (d - 100) / 365)
        T_min  = T_mean - 4.0
        T_max  = T_mean + 5.5
        rad = 15.0 + 8.0 * sin(2π * (d - 80) / 365)
        push!(days, DailyWeather(T_mean, T_min, T_max;
                                  radiation=rad, photoperiod=13.0))
    end
    return WeatherSeries(days; day_offset=1)
end

weather = make_berkeley_weather()

# Print diagnostics
println("=" ^ 60)
println("Plant–Aphid–Parasitoid (Gilbert & Gutierrez 1973)")
println("=" ^ 60)

println("\nAphid Brière: a=0.00012, T_base=$(APHID_T_BASE)°C, T_max=$(APHID_T_MAX)°C")
println("Parasitoid Brière: a=0.00006, T_base=$(PARA_T_BASE)°C, T_max=$(PARA_T_MAX)°C")
for T in [10.0, 15.0, 20.0, 25.0, 30.0]
    ra = development_rate(aphid_dev, T)
    rp = development_rate(para_dev, T)
    println("  T=$(T)°C: aphid r=$(round(ra, digits=4)), parasitoid r=$(round(rp, digits=4))")
end

println("\n--- Berkeley, CA Climate Snapshot ---")
for (label, d) in [("Jan", 15), ("Mar", 75), ("Jun", 165), ("Aug", 225), ("Oct", 290)]
    w = get_weather(weather, d)
    println("  $label: T_mean=$(round(w.T_mean, digits=1))°C, " *
            "range $(round(w.T_min, digits=1))–$(round(w.T_max, digits=1))°C")
end

println("\n--- Thimbleberry Leaf Capacity ---")
for (label, d) in [("Feb (dormant)", 45), ("Mar (flush)", 75),
                    ("May (expanding)", 130), ("Jun (peak)", 170),
                    ("Aug (plateau)", 210), ("Oct (senescent)", 260)]
    println("  $label (day $d): K = $(round(leaf_capacity(d), digits=0))")
end

# ============================================================
# Figure 1: Development rate vs temperature
# ============================================================

Ts = range(0.0, 35.0, length=300)
aphid_rates = [development_rate(aphid_dev, T) for T in Ts]
para_rates  = [development_rate(para_dev, T)  for T in Ts]

# Find optima
aphid_opt_idx = argmax(aphid_rates)
para_opt_idx  = argmax(para_rates)
aphid_T_opt = Ts[aphid_opt_idx]
para_T_opt  = Ts[para_opt_idx]

fig1 = Figure(size=(900, 550))
ax1 = Axis(fig1[1, 1],
    title="Development Rate vs Temperature — Aphid & Parasitoid",
    xlabel="Temperature (°C)",
    ylabel="Development rate (1/day)",
    xlabelsize=14, ylabelsize=14)

lines!(ax1, collect(Ts), aphid_rates, linewidth=2.5, color=:steelblue,
       label="Aphid (Masonaphis maxima)")
lines!(ax1, collect(Ts), para_rates, linewidth=2.5, color=:firebrick,
       label="Parasitoid (Praon sp.)")

# Annotate thermal optima
scatter!(ax1, [aphid_T_opt], [aphid_rates[aphid_opt_idx]], color=:steelblue, markersize=10)
scatter!(ax1, [para_T_opt], [para_rates[para_opt_idx]], color=:firebrick, markersize=10)
text!(ax1, aphid_T_opt + 0.5, aphid_rates[aphid_opt_idx],
      text="optimum $(round(aphid_T_opt, digits=1))°C",
      align=(:left, :center), fontsize=10, color=:steelblue)
text!(ax1, para_T_opt + 0.5, para_rates[para_opt_idx],
      text="optimum $(round(para_T_opt, digits=1))°C",
      align=(:left, :center), fontsize=10, color=:firebrick)

# Shade Berkeley growing season range
vspan!(ax1, 10.0, 22.0, color=(:green, 0.08))
text!(ax1, 16.0, maximum(aphid_rates) * 0.95,
      text="Berkeley\ngrowing season\nrange", align=(:center, :top),
      fontsize=9, color=(:darkgreen, 0.7))

# Base temperature annotations
vlines!(ax1, [APHID_T_BASE], color=:steelblue, linestyle=:dot, linewidth=0.8)
vlines!(ax1, [PARA_T_BASE], color=:firebrick, linestyle=:dot, linewidth=0.8)

xlims!(ax1, 0, 35)
ylims!(ax1, 0, nothing)
axislegend(ax1, position=:lt)

save(joinpath(figdir, "devrate_temperature.png"), fig1, px_per_unit=2)
println("\nSaved devrate_temperature.png — aphid optimum: $(round(aphid_T_opt, digits=1))°C, " *
        "parasitoid optimum: $(round(para_T_opt, digits=1))°C")

# ============================================================
# Figure 2: Plant resource phenology
# ============================================================

season_days = 1:365
lai_series = [leaf_area_index(d) for d in season_days]
K_series   = [leaf_capacity(d) for d in season_days]

# Monthly temperature for context
temps_monthly = [get_weather(weather, d).T_mean for d in season_days]

fig2 = Figure(size=(900, 550))
ax2 = Axis(fig2[1, 1],
    title="Thimbleberry (Rubus parviflorus) Leaf Phenology — Berkeley Hills",
    xlabel="Day of year",
    ylabel="Leaf area index (proportion of max)",
    xlabelsize=14, ylabelsize=14)

# LAI curve
band!(ax2, collect(season_days), zeros(365), lai_series,
      color=(:forestgreen, 0.2))
lines!(ax2, collect(season_days), lai_series, linewidth=2.5, color=:forestgreen,
       label="Leaf area index")

# Mark phenological phases
vlines!(ax2, [LEAF_FLUSH_DAY], color=:green, linestyle=:dash, linewidth=1)
vlines!(ax2, [PEAK_LEAF_DAY], color=:darkgreen, linestyle=:dash, linewidth=1)
vlines!(ax2, [SENESCENCE_DAY], color=:orange, linestyle=:dash, linewidth=1)
vlines!(ax2, [DORMANCY_DAY], color=:brown, linestyle=:dash, linewidth=1)

text!(ax2, LEAF_FLUSH_DAY + 2, 0.05, text="leaf\nflush", fontsize=9, color=:green,
      align=(:left, :bottom))
text!(ax2, PEAK_LEAF_DAY + 2, 0.95, text="peak\ncanopy", fontsize=9, color=:darkgreen,
      align=(:left, :top))
text!(ax2, SENESCENCE_DAY + 2, 0.85, text="senescence\nonset", fontsize=9, color=:orange,
      align=(:left, :top))
text!(ax2, DORMANCY_DAY - 2, 0.15, text="dormancy", fontsize=9, color=:brown,
      align=(:right, :bottom))

# Month labels
month_starts = [1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335]
month_labels = ["Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"]
ax2.xticks = (month_starts, month_labels)

# Temperature overlay
ax2b = Axis(fig2[1, 1],
    ylabel="Temperature (°C)",
    ylabelsize=12,
    yaxisposition=:right,
    yticklabelcolor=:darkorange,
    ylabelcolor=:darkorange)
hidexdecorations!(ax2b)
hidespines!(ax2b)
lines!(ax2b, collect(season_days), temps_monthly, linewidth=1.2,
       color=(:darkorange, 0.5), linestyle=:dot)
ylims!(ax2b, 5, 25)

xlims!(ax2, 1, 365)
ylims!(ax2, 0, 1.1)
axislegend(ax2, position=:lt)

save(joinpath(figdir, "plant_phenology.png"), fig2, px_per_unit=2)
println("Saved plant_phenology.png")

# ============================================================
# Figure 3: Aphid without parasitoid (constant 18°C)
# ============================================================

const T_CONST = 18.0
const N_SIM_CONST = 200

# Reproduction: adults produce nymphs, limited by leaf capacity
const FECUNDITY = 2.5  # nymphs per adult per DD unit
function reproduce_aphid_only(pop, w, p, day)
    T = w.T_mean
    adult_total = delay_total(pop.stages[2].delay)
    total_pop = delay_total(pop.stages[1].delay) + adult_total
    thermal_rate = development_rate(aphid_dev, T)
    K = leaf_capacity(day)
    density_effect = K > 0 ? max(0.0, 1.0 - total_pop / K) : 0.0
    return max(0.0, adult_total * thermal_rate * FECUNDITY * density_effect)
end

stages_const = [
    LifeStage(:nymph, DistributedDelay(15, DD_NYMPH; W0=10.0), aphid_dev, μ_NYMPH),
    LifeStage(:adult, DistributedDelay(10, DD_ADULT; W0=2.0),  aphid_dev, μ_ADULT),
]
pop_const = Population(:aphid, stages_const)

weather_const = WeatherSeries(fill(T_CONST, N_SIM_CONST); day_offset=LEAF_FLUSH_DAY)
prob_const = PBDMProblem(DensityDependent(), pop_const, weather_const,
                         (LEAF_FLUSH_DAY, LEAF_FLUSH_DAY + N_SIM_CONST);
                         p=nothing)
sol_const = solve(prob_const, DirectIteration(); reproduction_fn=reproduce_aphid_only)

nymph_const = stage_trajectory(sol_const, 1)
adult_const = stage_trajectory(sol_const, 2)
total_const = total_population(sol_const)

# Carrying capacity over same period
days_const = sol_const.t
K_const = [leaf_capacity(d) for d in days_const]

println("\n--- Aphid Without Parasitoid (constant $(T_CONST)°C) ---")
println("  Return code: $(sol_const.retcode)")
println("  Peak total: $(round(maximum(total_const), digits=1)) at day $(days_const[argmax(total_const)])")
println("  Final total: $(round(total_const[end], digits=1))")

fig3 = Figure(size=(900, 550))
ax3 = Axis(fig3[1, 1],
    title="Aphid Population Without Parasitoid (constant $(T_CONST)°C)",
    xlabel="Day of year",
    ylabel="Population",
    xlabelsize=14, ylabelsize=14)

# Carrying capacity envelope
band!(ax3, days_const, zeros(length(days_const)), K_const,
      color=(:green, 0.1))
lines!(ax3, days_const, K_const, linewidth=1, color=:green, linestyle=:dash,
       label="Carrying capacity K(t)")

lines!(ax3, days_const, nymph_const, linewidth=2, color=:steelblue, label="Nymphs")
lines!(ax3, days_const, adult_const, linewidth=2, color=:firebrick, label="Adults")
lines!(ax3, days_const, total_const, linewidth=2.5, color=:black, linestyle=:dash,
       label="Total")

# Annotate growth phase
peak_day = days_const[argmax(total_const)]
peak_val = maximum(total_const)
scatter!(ax3, [peak_day], [peak_val], color=:black, markersize=10)
text!(ax3, peak_day + 3, peak_val,
      text="peak day $(peak_day)\n$(round(Int, peak_val)) aphids",
      align=(:left, :center), fontsize=10, color=:black)

text!(ax3, LEAF_FLUSH_DAY + 15, K_MAX * 0.1,
      text="exponential\ngrowth phase", fontsize=10, color=(:steelblue, 0.7),
      align=(:left, :bottom))

xlims!(ax3, LEAF_FLUSH_DAY, LEAF_FLUSH_DAY + N_SIM_CONST)
ylims!(ax3, 0, nothing)
axislegend(ax3, position=:lt)

save(joinpath(figdir, "aphid_no_parasitoid.png"), fig3, px_per_unit=2)
println("Saved aphid_no_parasitoid.png")

# ============================================================
# Figure 4: Tritrophic dynamics (Berkeley climate)
# ============================================================

const PARA_ATTACK_RATE = 0.15
const PARA_INTRO_DAY = 100  # mid-April: parasitoid arrives

# --- Helper: run aphid-only simulation with Berkeley weather ---
function run_aphid_only_berkeley()
    stg = [
        LifeStage(:nymph, DistributedDelay(15, DD_NYMPH; W0=10.0), aphid_dev, μ_NYMPH),
        LifeStage(:adult, DistributedDelay(10, DD_ADULT; W0=2.0),  aphid_dev, μ_ADULT),
    ]
    pop = Population(:aphid, stg)
    t0, tf = 1, N_DAYS
    n = tf - t0

    nymph_traj = zeros(n + 1)
    adult_traj = zeros(n + 1)
    nymph_traj[1] = delay_total(pop.stages[1].delay)
    adult_traj[1] = delay_total(pop.stages[2].delay)

    for d in 1:n
        day = t0 + d - 1
        w = get_weather(weather, day)
        T = w.T_mean
        adult_total = delay_total(pop.stages[2].delay)
        total_pop = delay_total(pop.stages[1].delay) + adult_total
        thermal_rate = development_rate(aphid_dev, T)
        K = leaf_capacity(day)
        density_effect = K > 0 ? max(0.0, 1.0 - total_pop / K) : 0.0
        offspring = max(0.0, adult_total * thermal_rate * FECUNDITY * density_effect)
        pop.stages[1].delay.W[1] += offspring
        step_population!(pop, w)
        nymph_traj[d + 1] = delay_total(pop.stages[1].delay)
        adult_traj[d + 1] = delay_total(pop.stages[2].delay)
    end
    return nymph_traj .+ adult_traj
end

aphid_ref_traj = run_aphid_only_berkeley()

# --- Tritrophic simulation ---
stages_tri = [
    LifeStage(:nymph, DistributedDelay(15, DD_NYMPH; W0=10.0), aphid_dev, μ_NYMPH),
    LifeStage(:adult, DistributedDelay(10, DD_ADULT; W0=2.0),  aphid_dev, μ_ADULT),
]
pop_tri = Population(:aphid, stages_tri)

stages_para = [
    LifeStage(:larva, DistributedDelay(12, DD_PARA_LARVA; W0=0.0), para_dev, μ_PARA_LARVA),
    LifeStage(:adult, DistributedDelay(8, DD_PARA_ADULT; W0=0.0),  para_dev, μ_PARA_ADULT),
]
pop_para = Population(:parasitoid, stages_para)

t0, tf = 1, N_DAYS
n_days = tf - t0

aphid_nymph_traj = zeros(n_days + 1)
aphid_adult_traj = zeros(n_days + 1)
para_larva_traj  = zeros(n_days + 1)
para_adult_traj  = zeros(n_days + 1)
phi_plant_traj   = zeros(n_days)
phi_aphid_traj   = zeros(n_days)

aphid_nymph_traj[1] = delay_total(pop_tri.stages[1].delay)
aphid_adult_traj[1] = delay_total(pop_tri.stages[2].delay)
para_larva_traj[1]  = delay_total(pop_para.stages[1].delay)
para_adult_traj[1]  = delay_total(pop_para.stages[2].delay)

for d in 1:n_days
    day = t0 + d - 1
    w = get_weather(weather, day)
    T = w.T_mean

    # --- Plant-Aphid supply/demand ---
    aphid_total = delay_total(pop_tri.stages[1].delay) + delay_total(pop_tri.stages[2].delay)
    K = leaf_capacity(day)
    aphid_demand = max(aphid_total * 0.1, 1.0)
    phi_plant = K > 0 ? min(K / aphid_demand, 5.0) : 0.0
    phi_plant_traj[d] = phi_plant

    # --- Aphid reproduction ---
    adult_total = delay_total(pop_tri.stages[2].delay)
    thermal_rate = development_rate(aphid_dev, T)
    density_effect = K > 0 ? max(0.0, 1.0 - aphid_total / K) : 0.0

    # Parasitoid fecundity suppression (via supply-demand)
    para_total = delay_total(pop_para.stages[1].delay) + delay_total(pop_para.stages[2].delay)
    para_effect = 1.0
    if day >= PARA_INTRO_DAY && para_total > 0.5 && aphid_total > 0
        # Parasitoid demand on aphid supply — each parasitoid demands 20 hosts
        para_demand = para_total * 20.0
        phi_ap = supply_demand_ratio(FraserGilbertResponse(PARA_ATTACK_RATE),
                                     aphid_total, para_demand)
        phi_aphid_traj[d] = phi_ap
        # Dampen: parasitoid reduces fecundity partially (not full phi suppression)
        # Tuned to yield ~60% peak reduction consistent with Gilbert & Gutierrez (1973)
        para_effect = 1.0 - 0.65 * max(0.0, 1.0 - phi_ap)
        para_effect = max(0.15, para_effect)
    else
        phi_aphid_traj[d] = aphid_total > 0 ? 5.0 : 0.0
    end

    offspring = max(0.0, adult_total * thermal_rate * FECUNDITY * density_effect * para_effect)
    pop_tri.stages[1].delay.W[1] += offspring

    step_population!(pop_tri, w)

    # --- Parasitoid dynamics ---
    if day >= PARA_INTRO_DAY
        if day == PARA_INTRO_DAY
            pop_para.stages[1].delay.W[1] += 2.0
            pop_para.stages[2].delay.W[1] += 1.0
        end

        para_adults = delay_total(pop_para.stages[2].delay)
        aphid_nymphs = delay_total(pop_tri.stages[1].delay)
        if para_adults > 0.01 && aphid_nymphs > 0.01
            para_thermal = development_rate(para_dev, T)
            # Parasitoid reproduction: adults find and parasitize nymphs
            hosts_found = acquire(FraserGilbertResponse(PARA_ATTACK_RATE),
                                  aphid_nymphs, para_adults * 5.0)
            para_offspring = hosts_found * para_thermal * 1.5
            pop_para.stages[1].delay.W[1] += max(0.0, para_offspring)
        end

        step_population!(pop_para, w)
    end

    aphid_nymph_traj[d + 1] = delay_total(pop_tri.stages[1].delay)
    aphid_adult_traj[d + 1] = delay_total(pop_tri.stages[2].delay)
    para_larva_traj[d + 1]  = delay_total(pop_para.stages[1].delay)
    para_adult_traj[d + 1]  = delay_total(pop_para.stages[2].delay)
end

aphid_total_traj = aphid_nymph_traj .+ aphid_adult_traj
para_total_traj  = para_larva_traj .+ para_adult_traj
days_all = collect(t0:tf)

peak_ref = maximum(aphid_ref_traj)
peak_tri = maximum(aphid_total_traj)
reduction = peak_ref > 0 ? (1.0 - peak_tri / peak_ref) * 100 : 0.0

println("\n--- Tritrophic System (Berkeley Climate) ---")
println("  Aphid peak (no parasitoid): $(round(peak_ref, digits=0)) at day $(days_all[argmax(aphid_ref_traj)])")
println("  Aphid peak (with parasitoid): $(round(peak_tri, digits=0)) at day $(days_all[argmax(aphid_total_traj)])")
println("  Parasitoid peak: $(round(maximum(para_total_traj), digits=0)) at day $(days_all[argmax(para_total_traj)])")
println("  Peak reduction by parasitoid: $(round(reduction, digits=1))%")

fig4 = Figure(size=(1000, 650))

# Panel A: Population dynamics
ax4a = Axis(fig4[1, 1],
    title="(a) Tritrophic Dynamics — Berkeley, CA Climate",
    xlabel="Day of year",
    ylabel="Population",
    xlabelsize=14, ylabelsize=14)

# Carrying capacity envelope
K_all = [leaf_capacity(d) for d in days_all]
band!(ax4a, days_all, zeros(length(days_all)), K_all,
      color=(:green, 0.08))
lines!(ax4a, days_all, K_all, linewidth=0.8, color=:green, linestyle=:dash,
       label="K(t) — leaf capacity")

lines!(ax4a, days_all, aphid_ref_traj, linewidth=2, color=:red, linestyle=:dash,
       label="Aphid (no parasitoid)")
lines!(ax4a, days_all, aphid_total_traj, linewidth=2.5, color=:steelblue,
       label="Aphid (with Praon)")
lines!(ax4a, days_all, para_total_traj .* 20, linewidth=2, color=:firebrick,
       label="Parasitoid (×20)")

vlines!(ax4a, [PARA_INTRO_DAY], color=:red, linestyle=:dot, linewidth=0.8)
text!(ax4a, PARA_INTRO_DAY + 3, maximum(K_all) * 0.9,
      text="parasitoid\nintroduced", fontsize=9, color=:red)

month_starts_ax = [1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335]
month_labels_ax = ["Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"]
ax4a.xticks = (month_starts_ax, month_labels_ax)

xlims!(ax4a, 1, N_DAYS)
ylims!(ax4a, 0, nothing)
axislegend(ax4a, position=:rt)

# Panel B: Stage structure
ax4b = Axis(fig4[2, 1],
    title="(b) Stage-Structured Dynamics",
    xlabel="Day of year",
    ylabel="Population",
    xlabelsize=14, ylabelsize=14)

lines!(ax4b, days_all, aphid_nymph_traj, linewidth=2, color=:orange, label="Aphid nymphs")
lines!(ax4b, days_all, aphid_adult_traj, linewidth=2, color=:steelblue, label="Aphid adults")
lines!(ax4b, days_all, para_larva_traj, linewidth=1.5, color=:mediumpurple, label="Parasitoid larvae")
lines!(ax4b, days_all, para_adult_traj, linewidth=1.5, color=:firebrick, label="Parasitoid adults")

ax4b.xticks = (month_starts_ax, month_labels_ax)
xlims!(ax4b, 1, N_DAYS)
ylims!(ax4b, 0, nothing)
axislegend(ax4b, position=:rt)

# Temperature overlay on panel A
ax4t = Axis(fig4[1, 1],
    ylabel="Temperature (°C)",
    ylabelsize=11,
    yaxisposition=:right,
    yticklabelcolor=:darkorange,
    ylabelcolor=:darkorange)
hidexdecorations!(ax4t)
hidespines!(ax4t)
temps_plot = [get_weather(weather, d).T_mean for d in 1:N_DAYS]
lines!(ax4t, 1:N_DAYS, temps_plot, linewidth=1, color=(:darkorange, 0.4), linestyle=:dot)
ylims!(ax4t, 5, 25)

save(joinpath(figdir, "tritrophic_dynamics.png"), fig4, px_per_unit=2)
println("Saved tritrophic_dynamics.png")

# ============================================================
# Figure 5: Supply–Demand ratios
# ============================================================

fig5 = Figure(size=(900, 550))
ax5 = Axis(fig5[1, 1],
    title="Supply–Demand Ratios at Each Trophic Interface",
    xlabel="Day of year",
    ylabel="Supply / Demand ratio (ϕ)",
    xlabelsize=14, ylabelsize=14)

# Clip phi values for visualization
phi_plant_clipped = min.(phi_plant_traj, 5.0)
phi_aphid_clipped = min.(phi_aphid_traj, 5.0)
sim_days = collect(1:n_days)

lines!(ax5, sim_days, phi_plant_clipped, linewidth=2.5, color=:forestgreen,
       label="ϕ_plant (Plant → Aphid)")
lines!(ax5, sim_days, phi_aphid_clipped, linewidth=2.5, color=:mediumpurple,
       label="ϕ_aphid (Aphid → Parasitoid)")

# Critical threshold line
hlines!(ax5, [1.0], color=:black, linestyle=:dash, linewidth=1.5,
        label="ϕ = 1 (demand fully met)")

# Shade resource-limited zones
band!(ax5, sim_days,
      [min(p, 1.0) for p in phi_plant_clipped],
      ones(n_days),
      color=(:green, 0.05))

# Annotate key dynamics
text!(ax5, 140, 0.5,
      text="ϕ < 1 →\nresource limited", fontsize=10, color=(:black, 0.6),
      align=(:center, :center))
text!(ax5, 200, 3.5,
      text="ϕ > 1 →\ndemand satisfied", fontsize=10, color=(:black, 0.6),
      align=(:center, :center))

# Mark parasitoid introduction
vlines!(ax5, [PARA_INTRO_DAY], color=:red, linestyle=:dot, linewidth=0.8)
text!(ax5, PARA_INTRO_DAY + 3, 4.5,
      text="parasitoid\narrives", fontsize=9, color=:red)

ax5.xticks = (month_starts_ax, month_labels_ax)
xlims!(ax5, 1, N_DAYS)
ylims!(ax5, 0, 5.5)
axislegend(ax5, position=:rt)

save(joinpath(figdir, "supply_demand_ratios.png"), fig5, px_per_unit=2)
println("Saved supply_demand_ratios.png")

println("\n" * "=" ^ 60)
println("All plant–aphid–parasitoid validation figures saved to:\n  $(figdir)")
println("=" ^ 60)
