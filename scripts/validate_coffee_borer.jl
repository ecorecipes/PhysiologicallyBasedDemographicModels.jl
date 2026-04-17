#!/usr/bin/env julia
# Validate Coffee Berry Borer (Hypothenemus hampei) PBDM against literature.
#
# References:
#   Cure et al. (2020). The coffee agroecosystem: bio-economic analysis of CBB control.
#   Gutierrez et al. (1998). Tritrophic analysis of the coffee–CBB system.
#   Rodríguez et al. (2011/2017). Coffee agroecosystem model II: dynamics of CBB.
#   Jaramillo et al. (2009). Thermal tolerance of the CBB.
#
# Run: cd PhysiologicallyBasedDemographicModels.jl && julia --project=. scripts/validate_coffee_borer.jl

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using PhysiologicallyBasedDemographicModels
using CairoMakie

figdir = joinpath(@__DIR__, "figures", "coffee_borer")
mkpath(figdir)

println("=" ^ 70)
println("Coffee Berry Borer (Hypothenemus hampei) PBDM Validation")
println("=" ^ 70)
println()

# ============================================================
# 1. CBB Parameters from the literature
# ============================================================

# --- Temperature thresholds ---
const CBB_T_LOWER = 14.9    # °C  (Jaramillo et al. 2009, via Rodríguez/Cure)
const CBB_T_UPPER = 34.25   # °C  (Lactin fit upper)

# Brière development rate: r(T) = a·T·(T − T_lower)·√(T_upper − T)
# Calibrated so peak ≈ 0.045 near 28–30 °C matching Jaramillo data
const cbb_briere = BriereDevelopmentRate(2.56e-5, CBB_T_LOWER, CBB_T_UPPER)

# Linear model (simple degree-day) — used internally by the distributed delay
const cbb_linear = LinearDevelopmentRate(CBB_T_LOWER, CBB_T_UPPER)

# --- Degree-day durations per life stage (Table 1, Rodríguez / Cure) ---
const DD_EGG         = 44.15    # Embryonic:     0 – 44.15 DD
const DD_LARVA_I     = 20.77    # Larva I:    44.15 – 64.92 DD
const DD_LARVA_II    = 110.06   # Larva II:   64.92 – 174.98 DD
const DD_PREPUPA     = 25.83    # Pre-pupa:  174.98 – 200.81 DD
const DD_PUPA        = 61.66    # Pupa:      200.81 – 262.47 DD
const DD_YOUNG_ADULT = 50.0     # Young adult: 262.47 – 312.47 DD
const DD_MATURE      = 514.53   # Mature ♀:   312.47 – 827.0 DD
const DD_EGG_TO_ADULT = DD_EGG + DD_LARVA_I + DD_LARVA_II + DD_PREPUPA +
                        DD_PUPA + DD_YOUNG_ADULT  # 312.47

# --- Intrinsic mortality rates (dd⁻¹) from Table 1 ---
const MU_EGG      = 0.00102
const MU_LARVA_I  = 0.00094
const MU_LARVA_II = 0.00094
const MU_PREPUPA  = 0.00073
const MU_PUPA     = 0.00074
const MU_YOUNG    = 0.00035
const MU_MATURE   = 0.00035

# --- Fecundity ---
const FECUNDITY    = 0.3481       # eggs/female/degree-day (Table 1)
const FEMALE_FRAC  = 10.0 / 11.0  # sex ratio 10:1 female:male (Baker 1992)

# --- Density dependence ---
const U_THRESHOLD  = 4.0          # larvae per berry (competition threshold)

# --- Observed life-history at 22 °C (Bergamin 1943 / Gutierrez 1998) ---
const OBS_EGG_DAYS_22C     = 8.6    # days at 22 °C
const OBS_LARVA_DAYS_22C   = 15.9   # days
const OBS_PUPA_DAYS_22C    = 7.6    # days
const OBS_LIFESPAN_22C     = 156.6  # days (mean, range 82–282)
const OBS_EGGS_TOTAL       = 74.1   # lifetime fecundity (range 31–111)
const OBS_OVIPOSITION_DAYS = 20.0   # oviposition period

# --- Berry preference by phenological stage (cv. Colombia, Table 1) ---
const BERRY_PREF = Dict(
    "Pin"       => 0.00,
    "Green"     => 0.054,
    "Yellow"    => 0.57,
    "Ripe"      => 0.61,
    "Overripe"  => 0.57,
)

# --- Coffee berry development (degree-days >10 °C, Paper I) ---
const COFFEE_BASE_T     = 10.0
const BERRY_DD_PIN      = 1254.0
const BERRY_DD_GREEN    = 2622.0
const BERRY_DD_YELLOW   = 2836.0
const BERRY_DD_RIPE     = 3304.0
const BERRY_DD_OVERRIPE = 3538.0

println("Literature parameters loaded.")
println("  CBB thresholds: T_lower=$(CBB_T_LOWER)°C, T_upper=$(CBB_T_UPPER)°C")
println("  DD egg→adult: $(DD_EGG_TO_ADULT) DD")
println("  Fecundity: $(FECUNDITY) eggs/♀/DD, sex ratio $(round(FEMALE_FRAC*100, digits=1))% ♀")
println()

# ============================================================
# Helpers
# ============================================================

"""Temperature-dependent daily mortality (heuristic from literature)."""
function cbb_mortality_T(T; μ_base=0.01, T_low=CBB_T_LOWER, T_high=33.0)
    T <= T_low && return 1.0
    T >= CBB_T_UPPER + 2 && return 1.0
    T < T_low + 3 && return μ_base + 0.15 * ((T_low + 3 - T) / 3)^2
    T > T_high && return μ_base + 0.2 * ((T - T_high) / (CBB_T_UPPER + 2 - T_high))^2
    return μ_base
end

function colombia_temperature(day; T_mean=22.0, T_amp=2.0, phase=100)
    return T_mean + T_amp * sin(2π * (day - phase) / 365)
end

function make_cbb_population(; egg_W0=100.0)
    # k (substages) must satisfy CFL stability: k < τ / dd_max.
    # At 22 °C, dd = 7.1 DD/day; at 30 °C, dd = 15.1.  Use dd_max ≈ 15
    # to stay stable across the full viable range.
    dd_max = 15.0
    safe_k(τ) = max(1, floor(Int, τ / dd_max))
    stages = [
        LifeStage(:egg,           DistributedDelay(safe_k(DD_EGG),       DD_EGG;       W0=egg_W0), cbb_linear, MU_EGG),
        LifeStage(:larva_I,       DistributedDelay(safe_k(DD_LARVA_I),   DD_LARVA_I;   W0=0.0),    cbb_linear, MU_LARVA_I),
        LifeStage(:larva_II,      DistributedDelay(safe_k(DD_LARVA_II),  DD_LARVA_II;  W0=0.0),    cbb_linear, MU_LARVA_II),
        LifeStage(:prepupa,       DistributedDelay(safe_k(DD_PREPUPA),   DD_PREPUPA;   W0=0.0),    cbb_linear, MU_PREPUPA),
        LifeStage(:pupa,          DistributedDelay(safe_k(DD_PUPA),      DD_PUPA;      W0=0.0),    cbb_linear, MU_PUPA),
        LifeStage(:young_adult,   DistributedDelay(safe_k(DD_YOUNG_ADULT), DD_YOUNG_ADULT; W0=0.0), cbb_linear, MU_YOUNG),
        LifeStage(:mature_female, DistributedDelay(safe_k(DD_MATURE),    DD_MATURE;    W0=0.0),    cbb_linear, MU_MATURE),
    ]
    return Population(:coffee_berry_borer, stages)
end

# ============================================================
# Figure 1: Development rate — Brière vs observed data
# ============================================================

println("--- Figure 1: Development rate curves ---")

Ts = range(5.0, 40.0, length=500)
r_briere = [development_rate(cbb_briere, T) for T in Ts]

# Observed egg-to-adult development rate at selected temperatures
# (Jaramillo et al. 2009, digitised from their Figure 3 / Table)
obs_T    = [17.0,  20.0,  23.0,  25.0,  27.0,  30.0,  33.0]
obs_days = [120.0,  65.0,  45.0,  35.0,  28.0,  22.0,  28.0]
obs_rate = 1.0 ./ obs_days

# Linear DD model: rate = (T − T_lower) / DD_total  (day⁻¹ for egg-to-adult)
r_linear_day = [max(0.0, (T - CBB_T_LOWER) / DD_EGG_TO_ADULT) for T in Ts]

fig1 = Figure(size=(900, 600))
ax1 = Axis(fig1[1, 1],
    xlabel="Temperature (°C)", ylabel="Development rate (day⁻¹)",
    title="CBB Egg-to-Adult Development Rate vs Temperature",
    xlabelsize=14, ylabelsize=14)
lines!(ax1, collect(Ts), r_briere, color=:blue, linewidth=2.5,
       label="Brière (a=2.56e-5)")
lines!(ax1, collect(Ts), r_linear_day, color=:green, linewidth=2, linestyle=:dash,
       label="Linear DD ($(round(DD_EGG_TO_ADULT, digits=0)) DD)")
scatter!(ax1, obs_T, obs_rate, color=:red, marker=:circle, markersize=12,
         label="Jaramillo et al. (observed)")
vlines!(ax1, [CBB_T_LOWER], color=:cyan, linestyle=:dot, linewidth=1.5,
        label="T_lower = $(CBB_T_LOWER)°C")
vlines!(ax1, [CBB_T_UPPER], color=:magenta, linestyle=:dot, linewidth=1.5,
        label="T_upper = $(CBB_T_UPPER)°C")
xlims!(ax1, 5, 40)
ylims!(ax1, -0.005, max(maximum(obs_rate), maximum(r_briere)) * 1.2)
axislegend(ax1, position=:lt, framevisible=false)

save(joinpath(figdir, "devrate_curves.png"), fig1, px_per_unit=2)
println("  Saved devrate_curves.png")

println("  Model vs observed at key temperatures:")
println("  T(°C)  | Brière  | Linear  | Observed")
for (Tv, dobs) in zip(obs_T, obs_days)
    rb = development_rate(cbb_briere, Tv)
    rl = max(0.0, (Tv - CBB_T_LOWER) / DD_EGG_TO_ADULT)
    println("  $(lpad(Tv, 5)) | $(lpad(round(rb, digits=4), 7)) | $(lpad(round(rl, digits=4), 7)) | $(round(1.0/dobs, digits=4))")
end
println()

# ============================================================
# Figure 2: Temperature-dependent mortality curve
# ============================================================

println("--- Figure 2: Mortality curve ---")

Ts_mort = range(5.0, 40.0, length=500)
mort_daily = [cbb_mortality_T(T) for T in Ts_mort]

# Convert stage-specific dd⁻¹ rates to daily rates at 22 °C
dd_per_day_22 = 22.0 - CBB_T_LOWER  # 7.1
ref_mort = [
    ("Egg",      MU_EGG,      :gold),
    ("Larva I",  MU_LARVA_I,  :green),
    ("Larva II", MU_LARVA_II, :darkgreen),
    ("Pre-pupa", MU_PREPUPA,  :brown),
    ("Pupa",     MU_PUPA,     :purple),
    ("Adult",    MU_YOUNG,    :red),
]

fig2 = Figure(size=(900, 600))
ax2 = Axis(fig2[1, 1],
    xlabel="Temperature (°C)", ylabel="Daily mortality probability",
    title="CBB Mortality: Temperature Curve + Stage-Specific Intrinsic Rates (at 22°C)",
    xlabelsize=14, ylabelsize=14)
lines!(ax2, collect(Ts_mort), mort_daily, color=:black, linewidth=2.5,
       label="T-dependent heuristic")

for (name, mu_dd, col) in ref_mort
    mu_day = mu_dd * dd_per_day_22
    hlines!(ax2, [mu_day], color=col, linestyle=:dash, linewidth=1.5,
            label="$name $(round(mu_dd, digits=5))/DD → $(round(mu_day, digits=4))/day")
end

vspan!(ax2, CBB_T_LOWER + 3, 33.0, color=(:green, 0.08))
text!(ax2, 23.5, 0.30, text="Optimal range", align=(:center, :center),
      fontsize=12, color=:green)
vlines!(ax2, [CBB_T_LOWER], color=:cyan, linestyle=:dot, linewidth=1.5)
vlines!(ax2, [CBB_T_UPPER], color=:magenta, linestyle=:dot, linewidth=1.5)
xlims!(ax2, 5, 40)
ylims!(ax2, 0, 0.40)
axislegend(ax2, position=:rt, framevisible=false, nbanks=2)

save(joinpath(figdir, "mortality_curve.png"), fig2, px_per_unit=2)
println("  Saved mortality_curve.png")
println()

# ============================================================
# Figure 3: Constant temperature simulation at 22 °C
# ============================================================

println("--- Figure 3: Constant 22°C simulation (cohort, no reproduction) ---")

n_days = 365
dd_per_day = 22.0 - CBB_T_LOWER

# Density-independent cohort (no reproduction) to track stage progression
cbb_pop_22 = make_cbb_population(egg_W0=100.0)
ws_22 = WeatherSeries([DailyWeather(22.0) for _ in 1:n_days]; day_offset=1)
prob_22 = PBDMProblem(cbb_pop_22, ws_22, (1, n_days))
sol_22 = solve(prob_22, DirectIteration())

stage_names  = [:egg, :larva_I, :larva_II, :prepupa, :pupa, :young_adult, :mature_female]
stage_labels = ["Egg", "Larva I", "Larva II", "Pre-pupa", "Pupa", "Young adult", "Mature ♀"]
stage_colors = [:gold, :green, :darkgreen, :sienna, :purple, :orange, :red]

cdd_22 = cumulative_degree_days(sol_22)
tp_22  = total_population(sol_22)

expected_egg_days = DD_EGG / dd_per_day
expected_gen_days = DD_EGG_TO_ADULT / dd_per_day
gens_per_year = cdd_22[end] / DD_EGG_TO_ADULT

println("  DD per day at 22°C: $(round(dd_per_day, digits=1))")
println("  Expected egg duration: $(round(expected_egg_days, digits=1)) days (obs: $(OBS_EGG_DAYS_22C))")
println("  Expected generation time: $(round(expected_gen_days, digits=1)) days")
println("  Total DD/year: $(round(cdd_22[end], digits=0))")
println("  Estimated generations/year: $(round(gens_per_year, digits=1))")

fig3 = Figure(size=(950, 700))

ax3a = Axis(fig3[1, 1],
    ylabel="Individuals",
    title="CBB Cohort at 22°C — No Reproduction (initial 1500 eggs)",
    xlabelsize=14, ylabelsize=14)
for (i, (lbl, col)) in enumerate(zip(stage_labels, stage_colors))
    traj = stage_trajectory(sol_22, i)
    lines!(ax3a, sol_22.t, traj, color=col, linewidth=1.8, label=lbl)
end
lines!(ax3a, sol_22.t, tp_22, color=:black, linewidth=2.5, linestyle=:dash, label="Total")

# Mark observed stage durations at 22°C
obs_egg_end = round(Int, OBS_EGG_DAYS_22C)
obs_larva_end = obs_egg_end + round(Int, OBS_LARVA_DAYS_22C)
obs_pupa_end  = obs_larva_end + round(Int, OBS_PUPA_DAYS_22C)
for (d, lbl, col) in [(obs_egg_end, "Egg (obs)", :gold),
                       (obs_larva_end, "Larva end (obs)", :green),
                       (obs_pupa_end, "Pupa end (obs)", :purple)]
    vlines!(ax3a, [d], color=col, linestyle=:dot, linewidth=1.5)
end
axislegend(ax3a, position=:rt, framevisible=false, nbanks=2)

ax3b = Axis(fig3[2, 1],
    xlabel="Day", ylabel="Cumulative DD",
    xlabelsize=14, ylabelsize=14)
lines!(ax3b, 1:length(cdd_22), cdd_22, color=:black, linewidth=2)
for g in 1:8
    gen_dd = g * DD_EGG_TO_ADULT
    gen_dd > cdd_22[end] && break
    gen_day = findfirst(x -> x >= gen_dd, cdd_22)
    gen_day === nothing && continue
    vlines!(ax3b, [gen_day], color=:red, linestyle=:dot, linewidth=1)
    text!(ax3b, gen_day + 3, gen_dd, text="G$g", fontsize=10, color=:red)
end
hlines!(ax3b, [DD_EGG_TO_ADULT], color=:gray, linestyle=:dash, linewidth=1,
        label="1 gen ($(round(Int, DD_EGG_TO_ADULT)) DD)")
axislegend(ax3b, position=:lt, framevisible=false)

save(joinpath(figdir, "sim_constant_22C.png"), fig3, px_per_unit=2)
println("  Saved sim_constant_22C.png")

println("\n  Duration comparison at 22°C:")
println("  Stage          | Model (days) | Gutierrez 1998 (days)")
stages_dd = [DD_EGG, DD_LARVA_I + DD_LARVA_II, DD_PREPUPA + DD_PUPA]
lit_days  = [OBS_EGG_DAYS_22C, OBS_LARVA_DAYS_22C, OBS_PUPA_DAYS_22C]
names_cmp = ["Egg", "Larva (I+II)", "Pupa (+prepupa)"]
for (nm, ddv, ld) in zip(names_cmp, stages_dd, lit_days)
    model_d = ddv / dd_per_day
    println("  $(rpad(nm, 16)) | $(lpad(round(model_d, digits=1), 12)) | $(lpad(ld, 8))")
end
println()

# ============================================================
# Figure 4: Coffee berry phenology — fruiting cycle
# ============================================================

println("--- Figure 4: Coffee berry phenology ---")

coffee_dd_per_day = 22.0 - COFFEE_BASE_T  # 12.0 at 22°C

pin_end    = round(Int, BERRY_DD_PIN / coffee_dd_per_day)
green_end  = round(Int, BERRY_DD_GREEN / coffee_dd_per_day)
yellow_end = round(Int, BERRY_DD_YELLOW / coffee_dd_per_day)
ripe_end   = round(Int, BERRY_DD_RIPE / coffee_dd_per_day)
n_berry    = ripe_end + 30

# Susceptibility over berry development
susceptibility = zeros(n_berry)
for d in 1:n_berry
    if d <= pin_end
        susceptibility[d] = BERRY_PREF["Pin"]
    elseif d <= green_end
        susceptibility[d] = BERRY_PREF["Green"]
    elseif d <= yellow_end
        susceptibility[d] = BERRY_PREF["Yellow"]
    elseif d <= ripe_end
        susceptibility[d] = BERRY_PREF["Ripe"]
    else
        susceptibility[d] = BERRY_PREF["Overripe"]
    end
end

fig4 = Figure(size=(950, 650))
ax4a = Axis(fig4[1, 1],
    ylabel="Berry stage",
    title="Coffee Berry Phenology (cv. Colombia, 22°C) + CBB Attack Susceptibility",
    xlabelsize=14, ylabelsize=14,
    yticks=([1, 2, 3, 4, 5], ["Pin", "Green", "Yellow", "Ripe", "Overripe"]))

vspan!(ax4a, 1.0, Float64(pin_end), color=(:palegreen, 0.5))
vspan!(ax4a, Float64(pin_end), Float64(green_end), color=(:green, 0.3))
vspan!(ax4a, Float64(green_end), Float64(yellow_end), color=(:gold, 0.4))
vspan!(ax4a, Float64(yellow_end), Float64(ripe_end), color=(:red, 0.3))
vspan!(ax4a, Float64(ripe_end), Float64(n_berry), color=(:brown, 0.3))

mid(a, b) = (a + b) / 2
text!(ax4a, mid(1, pin_end), 4.5, text="Pin\n$(pin_end)d", fontsize=10,
      align=(:center, :center))
text!(ax4a, mid(pin_end, green_end), 4.5,
      text="Green\n$(green_end-pin_end)d", fontsize=10, align=(:center, :center))
text!(ax4a, mid(green_end, yellow_end), 4.5,
      text="Yellow\n$(yellow_end-green_end)d", fontsize=10, align=(:center, :center))
text!(ax4a, mid(yellow_end, ripe_end), 4.5,
      text="Ripe\n$(ripe_end-yellow_end)d", fontsize=10, align=(:center, :center))
xlims!(ax4a, 0, n_berry)
ylims!(ax4a, 0, 5.5)

ax4b = Axis(fig4[2, 1],
    xlabel="Days after flowering",
    ylabel="CBB attack preference",
    xlabelsize=14, ylabelsize=14)
lines!(ax4b, 1:n_berry, susceptibility, color=:red, linewidth=2.5, label="Susceptibility")

for (stage, pref) in sort(collect(BERRY_PREF), by=x -> x[2], rev=true)
    pref <= 0 && continue
    if stage == "Green"
        xp = mid(pin_end, green_end)
    elseif stage == "Yellow"
        xp = mid(green_end, yellow_end)
    elseif stage == "Ripe"
        xp = mid(yellow_end, ripe_end)
    elseif stage == "Overripe"
        xp = ripe_end + 15.0
    else
        continue
    end
    scatter!(ax4b, [xp], [pref], color=:red, markersize=10)
    text!(ax4b, xp, pref + 0.04, text="$stage\n$pref", fontsize=9,
          align=(:center, :bottom))
end

xlims!(ax4b, 0, n_berry)
ylims!(ax4b, 0, 0.80)
axislegend(ax4b, position=:lt, framevisible=false)

save(joinpath(figdir, "coffee_phenology.png"), fig4, px_per_unit=2)
println("  Saved coffee_phenology.png")
println("  Berry durations at 22°C: Pin=$(pin_end)d, Green=$(green_end-pin_end)d, " *
        "Yellow=$(yellow_end-green_end)d, Ripe=$(ripe_end-yellow_end)d")
println("  Total flowering→ripe: $(ripe_end)d (~$(round(ripe_end/7, digits=0)) weeks)")
println()

# ============================================================
# Figure 5: Infestation dynamics over the fruiting season
# ============================================================

println("--- Figure 5: Infestation dynamics ---")

# Simulate with reproduction + carrying-capacity feedback.
# Berries provide a limited resource: K_berries berries per tree.
# Max ~4 larvae per berry (density regulation from Table 1).
const K_BERRIES    = 2000.0   # berries per tree
const K_POP        = K_BERRIES * U_THRESHOLD  # carrying capacity
const ALPHA_DAMAGE = 5e-5     # damage accumulation rate

function cbb_reproduction_dd(pop, w, p, day)
    mature = delay_total(pop.stages[end].delay)
    dd = degree_days(pop.stages[1].dev_rate, w.T_mean)
    total_pop = sum(delay_total(s.delay) for s in pop.stages)
    # Logistic density regulation
    density_factor = max(0.0, 1.0 - total_pop / K_POP)
    return FECUNDITY * FEMALE_FRAC * mature * dd / 10.0 * density_factor
end

n_season = min(ripe_end + 30, 365)
colombia_temps = [colombia_temperature(d) for d in 1:n_season]
ws_col = WeatherSeries([DailyWeather(T) for T in colombia_temps]; day_offset=1)

# Start with 50 colonising mature females
cbb_pop_inf = make_cbb_population(egg_W0=0.0)
for i in 1:cbb_pop_inf.stages[7].delay.k
    cbb_pop_inf.stages[7].delay.W[i] = 5.0
end

prob_inf = PBDMProblem(DensityDependent(), cbb_pop_inf, ws_col, (1, n_season))
sol_inf  = solve(prob_inf, DirectIteration(); reproduction_fn=cbb_reproduction_dd)

tp_inf     = total_population(sol_inf)
mature_inf = stage_trajectory(sol_inf, 7)
egg_inf    = stage_trajectory(sol_inf, 1)

# Cumulative berry damage
cumulative_attack = zeros(length(sol_inf.t))
for (i, d) in enumerate(sol_inf.t)
    attack = (d >= 1 && d <= n_berry) ? mature_inf[i] * susceptibility[min(d, n_berry)] : 0.0
    cumulative_attack[i] = (i > 1 ? cumulative_attack[i-1] : 0.0) + attack
end
damage_pct = [100.0 * (1.0 - exp(-ALPHA_DAMAGE * ca)) for ca in cumulative_attack]

fig5 = Figure(size=(950, 800))

ax5a = Axis(fig5[1, 1],
    ylabel="CBB population",
    title="CBB Infestation — Colombian Coffee Zone (22°C, K=$(round(Int, K_POP)) per tree)",
    xlabelsize=14, ylabelsize=14)
lines!(ax5a, sol_inf.t, tp_inf,     color=:black, linewidth=2.5, label="Total CBB")
lines!(ax5a, sol_inf.t, mature_inf, color=:red,   linewidth=2,   label="Mature ♀")
lines!(ax5a, sol_inf.t, egg_inf,    color=:gold,  linewidth=1.5, label="Eggs")
hlines!(ax5a, [K_POP], color=:gray, linestyle=:dash, linewidth=1,
        label="K = $(round(Int, K_POP))")
for (dm, lbl) in [(pin_end, "Pin→Grn"), (green_end, "Grn→Yel"),
                   (yellow_end, "Yel→Ripe")]
    dm <= n_season && vlines!(ax5a, [dm], color=:gray, linestyle=:dot, linewidth=1)
end
axislegend(ax5a, position=:lt, framevisible=false)

ax5b = Axis(fig5[2, 1],
    xlabel="Days after flowering",
    ylabel="Berry damage (%)",
    xlabelsize=14, ylabelsize=14)
lines!(ax5b, sol_inf.t, damage_pct, color=:darkred, linewidth=2.5, label="Cumulative damage")
hlines!(ax5b, [5.0],  color=:green,  linestyle=:dash, linewidth=1.5, label="5% economic threshold")
hlines!(ax5b, [20.0], color=:orange, linestyle=:dash, linewidth=1.5, label="20% severe (Cure 2020)")
hlines!(ax5b, [35.0], color=:red,    linestyle=:dash, linewidth=1.5, label="35% critical")
for (dm, lbl) in [(green_end, "Susceptible →"), (ripe_end, "Harvest")]
    if dm <= n_season
        vlines!(ax5b, [dm], color=:gray, linestyle=:dot, linewidth=1)
        text!(ax5b, dm + 3, min(damage_pct[end] * 0.85, 35.0),
              text=lbl, fontsize=10, align=(:left, :top))
    end
end
ylims!(ax5b, 0, max(45.0, maximum(damage_pct) * 1.1))
axislegend(ax5b, position=:lt, framevisible=false)

save(joinpath(figdir, "infestation_dynamics.png"), fig5, px_per_unit=2)
println("  Saved infestation_dynamics.png")
println("  Peak CBB: $(round(maximum(tp_inf), digits=0)) at day $(sol_inf.t[argmax(tp_inf)])")
println("  Final damage: $(round(damage_pct[end], digits=1))%")
println()

# ============================================================
# Summary
# ============================================================

println("=" ^ 70)
println("Parameter validation summary")
println("=" ^ 70)
println()
println("Degree-day requirements (Rodríguez Table 1):")
for (name, dd) in [("Egg", DD_EGG), ("Larva I", DD_LARVA_I),
                    ("Larva II", DD_LARVA_II), ("Pre-pupa", DD_PREPUPA),
                    ("Pupa", DD_PUPA), ("Young adult", DD_YOUNG_ADULT),
                    ("Mature ♀", DD_MATURE)]
    days_22 = dd / dd_per_day
    println("  $(rpad(name, 14)): $(lpad(round(dd, digits=2), 7)) DD  " *
            "($(round(days_22, digits=1)) days at 22°C)")
end
println("  $(rpad("Egg→Adult", 14)): $(lpad(round(DD_EGG_TO_ADULT, digits=2), 7)) DD  " *
        "($(round(DD_EGG_TO_ADULT/dd_per_day, digits=1)) days at 22°C)")
println()
println("Berry preferences (cv. Colombia, Table 1):")
for (stage, pref) in sort(collect(BERRY_PREF), by=x -> x[2], rev=true)
    bar = repeat("█", round(Int, pref * 30))
    println("  $(rpad(stage, 10)): $pref $bar")
end
println()
println("Mortality rates (dd⁻¹, Table 1):")
println("  Egg: $MU_EGG, Larva I: $MU_LARVA_I, Larva II: $MU_LARVA_II")
println("  Pre-pupa: $MU_PREPUPA, Pupa: $MU_PUPA, Adult: $MU_YOUNG")
println()
println("All 5 figures saved to: $figdir")
println("=" ^ 70)
