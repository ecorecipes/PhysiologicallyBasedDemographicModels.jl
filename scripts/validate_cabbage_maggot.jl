#!/usr/bin/env julia
# Validation script for Cabbage Root Fly (Delia radicum) PBDM
# Based on Collier & Finch (1985), Johnsen et al. (1997)
#
# Generates figures in scripts/figures/cabbage_maggot/ for:
#   1. Temperature-dependent development rate curves (Brière)
#   2. Photoperiod response and diapause induction window
#   3. Cohort dynamics at constant 18°C
#   4. Seasonal UK simulation with photoperiod
#   5. Latitude comparison (UK, France, Scotland)

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using PhysiologicallyBasedDemographicModels
using CairoMakie

figdir = joinpath(@__DIR__, "figures", "cabbage_maggot")
mkpath(figdir)

# ============================================================
# Parameters from literature
# Collier & Finch (1985), Johnsen et al. (1997)
# ============================================================

const T_LOWER = 4.3     # Lower developmental threshold (°C)
const T_UPPER = 32.0    # Upper developmental threshold (°C)

# Degree-day requirements above T_LOWER
const DD_EGG   = 55.0   # Egg stage
const DD_LARVA = 250.0  # Larval stage (3 instars combined)
const DD_PUPA  = 200.0  # Pupal stage
const DD_ADULT = 400.0  # Adult reproductive lifespan
# Total egg-to-adult: ~505 DD

# Diapause parameters
const CRITICAL_PHOTOPERIOD = 14.5  # hours; short-day diapause induced below this

# Brière development rate coefficients (calibrated to literature data)
# r(T) = a · T · (T − T_lower) · √(T_upper − T)
const A_EGG   = 1.50e-4  # Egg: fitted to (10°C,0.05), (15°C,0.10), (20°C,0.17)
const A_LARVA = 3.10e-5  # Larva: fitted to (15°C,0.02), (20°C,0.035)
const A_PUPA  = 4.50e-5  # Pupa: interpolated from DD ratio

# Background mortality rates (per degree-day)
# Higher rates create faster cohort turnover → distinct generational peaks.
const MU_EGG   = 0.005
const MU_LARVA = 0.003
const MU_PUPA  = 0.002
const MU_ADULT = 0.005

# Substage counts (k) for distributed delay
# Stability requires (k/τ)*dd_max + μ*dd_max < 1 where dd_max ≈ 15.7 (France 20°C).
const K_EGG   = 3    # (3/55)*15.7 + 0.005*15.7 = 0.94
const K_LARVA = 14   # (14/250)*15.7 + 0.003*15.7 = 0.93
const K_PUPA  = 11   # (11/200)*15.7 + 0.002*15.7 = 0.90
const K_ADULT = 20   # (20/400)*15.7 + 0.005*15.7 = 0.86

println("=" ^ 60)
println("Delia radicum (Cabbage Root Fly) PBDM Validation")
println("=" ^ 60)
println("Lower threshold: $(T_LOWER)°C")
println("Upper threshold: $(T_UPPER)°C")
println("DD requirements: egg=$(DD_EGG), larva=$(DD_LARVA), pupa=$(DD_PUPA)")
println("Total egg-to-adult: $(DD_EGG + DD_LARVA + DD_PUPA) DD")
println("Critical photoperiod for diapause: $(CRITICAL_PHOTOPERIOD) h")

# ============================================================
# Helper: build Delia radicum population
# ============================================================

"""
Build a fresh Delia radicum population with specified initial numbers.
Initial individuals are distributed evenly across substages of the
specified life stage. All other stages start empty.
"""
function make_delia_population(;
    egg_init::Float64=0.0,
    larva_init::Float64=0.0,
    pupa_init::Float64=0.0,
    adult_init::Float64=0.0
)
    dev_rate = LinearDevelopmentRate(T_LOWER, T_UPPER)
    stages = [
        LifeStage(:egg,   DistributedDelay(K_EGG,   DD_EGG;   W0=egg_init   / K_EGG),   dev_rate, MU_EGG),
        LifeStage(:larva, DistributedDelay(K_LARVA,  DD_LARVA; W0=larva_init / K_LARVA),  dev_rate, MU_LARVA),
        LifeStage(:pupa,  DistributedDelay(K_PUPA,   DD_PUPA;  W0=pupa_init  / K_PUPA),   dev_rate, MU_PUPA),
        LifeStage(:adult, DistributedDelay(K_ADULT,  DD_ADULT; W0=adult_init / K_ADULT),  dev_rate, MU_ADULT),
    ]
    Population(:delia_radicum, stages)
end

# Reproduction: ~100 eggs per female over 400 DD lifespan, 50% female
# Tuned with higher mortality rates to produce visible multi-generational
# peaks in the density-dependent reproduction model.
const FECUNDITY = 0.05  # eggs per DD per individual (adjusted for model behavior)

"""
Reproduction function for density-dependent solver.
Adults produce eggs proportional to degree-day accumulation.
"""
function delia_reproduction(pop, w, p, day)
    adult_total = delay_total(pop.stages[4].delay)
    dd = degree_days(pop.stages[1].dev_rate, w.T_mean)
    return p.fecundity * dd * adult_total
end

const STAGE_NAMES  = [:egg, :larva, :pupa, :adult]
const STAGE_COLORS = [:goldenrod, :forestgreen, :steelblue, :firebrick]

# Month tick marks for seasonal plots
const MONTH_STARTS = [1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335]
const MONTH_LABELS = ["Jan","Feb","Mar","Apr","May","Jun",
                      "Jul","Aug","Sep","Oct","Nov","Dec"]

# ============================================================
# Figure 1: Development rate curves (Brière)
# ============================================================

egg_briere   = BriereDevelopmentRate(A_EGG,   T_LOWER, T_UPPER)
larva_briere = BriereDevelopmentRate(A_LARVA, T_LOWER, T_UPPER)
pupa_briere  = BriereDevelopmentRate(A_PUPA,  T_LOWER, T_UPPER)

Ts = range(0.0, 38.0, length=300)
egg_rates   = [development_rate(egg_briere,   T) for T in Ts]
larva_rates = [development_rate(larva_briere, T) for T in Ts]
pupa_rates  = [development_rate(pupa_briere,  T) for T in Ts]

# Observed data points from literature
obs_egg_T   = [10.0, 15.0, 20.0]
obs_egg_r   = [0.05, 0.10, 0.17]
obs_larva_T = [15.0, 20.0]
obs_larva_r = [0.02, 0.035]

fig1 = Figure(size=(900, 600))
ax1 = Axis(fig1[1, 1],
    title="Delia radicum Development Rate Curves (Brière model)\nCollier & Finch (1985), Johnsen et al. (1997)",
    xlabel="Temperature (°C)",
    ylabel="Development rate (1/day)",
    xlabelsize=14, ylabelsize=14)

lines!(ax1, collect(Ts), egg_rates,   linewidth=2.5, color=:goldenrod,
       label="Egg (a=$(A_EGG))")
lines!(ax1, collect(Ts), larva_rates, linewidth=2.5, color=:forestgreen,
       label="Larva (a=$(A_LARVA))")
lines!(ax1, collect(Ts), pupa_rates,  linewidth=2.5, color=:steelblue,
       label="Pupa (a=$(A_PUPA))")

scatter!(ax1, obs_egg_T, obs_egg_r, color=:goldenrod, markersize=12,
         marker=:circle, label="Egg obs.")
scatter!(ax1, obs_larva_T, obs_larva_r, color=:forestgreen, markersize=12,
         marker=:utriangle, label="Larva obs.")

vlines!(ax1, [T_LOWER], color=:gray, linestyle=:dash, linewidth=1)
vlines!(ax1, [T_UPPER], color=:gray, linestyle=:dash, linewidth=1)
text!(ax1, T_LOWER + 0.5, maximum(egg_rates) * 0.95,
    text="T_lower=$(T_LOWER)°C", fontsize=10, color=:gray)
text!(ax1, T_UPPER - 6.0, maximum(egg_rates) * 0.95,
    text="T_upper=$(T_UPPER)°C", fontsize=10, color=:gray)

xlims!(ax1, 0, 38)
ylims!(ax1, 0, nothing)
axislegend(ax1, position=:lt)

save(joinpath(figdir, "devrate_curves.png"), fig1, px_per_unit=2)
println("\nSaved devrate_curves.png")
println("  Egg peak rate:   $(round(maximum(egg_rates), digits=4)) at $(round(Ts[argmax(egg_rates)], digits=1))°C")
println("  Larva peak rate: $(round(maximum(larva_rates), digits=4)) at $(round(Ts[argmax(larva_rates)], digits=1))°C")
println("  Pupa peak rate:  $(round(maximum(pupa_rates), digits=4)) at $(round(Ts[argmax(pupa_rates)], digits=1))°C")

# ============================================================
# Figure 2: Photoperiod response and diapause induction
# ============================================================

days_year = 1:365
pp_uk     = [photoperiod(52.0, d) for d in days_year]
pp_france = [photoperiod(45.0, d) for d in days_year]

fig2 = Figure(size=(900, 600))
ax2 = Axis(fig2[1, 1],
    title="Photoperiod and Diapause Induction Window\nDelia radicum — Critical Photoperiod $(CRITICAL_PHOTOPERIOD) h",
    xlabel="Day of year",
    ylabel="Photoperiod (hours)",
    xlabelsize=14, ylabelsize=14,
    xticks=(MONTH_STARTS, MONTH_LABELS))

lines!(ax2, collect(days_year), pp_uk, linewidth=2.5, color=:steelblue,
       label="52°N (UK)")
lines!(ax2, collect(days_year), pp_france, linewidth=2.5, color=:firebrick,
       label="45°N (France)")

hlines!(ax2, [CRITICAL_PHOTOPERIOD], color=:black, linestyle=:dash, linewidth=2,
        label="Critical $(CRITICAL_PHOTOPERIOD) h")

# Diapause induction window: late-summer days when photoperiod drops below critical
uk_diapause_days = [d for d in 180:365 if pp_uk[d] < CRITICAL_PHOTOPERIOD]
if !isempty(uk_diapause_days)
    diapause_start_uk = first(uk_diapause_days)
    vspan!(ax2, diapause_start_uk, 365, color=(:steelblue, 0.1))
    text!(ax2, (diapause_start_uk + 365) / 2, CRITICAL_PHOTOPERIOD + 0.5,
        text="Diapause\ninduction (UK)", align=(:center, :bottom),
        fontsize=11, color=:steelblue)
end

fr_diapause_days = [d for d in 180:365 if pp_france[d] < CRITICAL_PHOTOPERIOD]
if !isempty(fr_diapause_days)
    diapause_start_fr = first(fr_diapause_days)
    vlines!(ax2, [diapause_start_fr], color=:firebrick, linestyle=:dot, linewidth=1.5)
    text!(ax2, diapause_start_fr - 2, CRITICAL_PHOTOPERIOD - 0.8,
        text="France\nonset", align=(:right, :top), fontsize=10, color=:firebrick)
end

xlims!(ax2, 1, 365)
ylims!(ax2, 6, 20)
axislegend(ax2, position=:cb)

save(joinpath(figdir, "photoperiod_response.png"), fig2, px_per_unit=2)
println("\nSaved photoperiod_response.png")
println("  UK photoperiod range:     $(round(minimum(pp_uk), digits=1))–$(round(maximum(pp_uk), digits=1)) h")
println("  France photoperiod range: $(round(minimum(pp_france), digits=1))–$(round(maximum(pp_france), digits=1)) h")
if !isempty(uk_diapause_days)
    println("  UK diapause onset:        day $(first(uk_diapause_days))")
end
if !isempty(fr_diapause_days)
    println("  France diapause onset:    day $(first(fr_diapause_days))")
end

# ============================================================
# Figure 3: Cohort dynamics at constant 18°C
# ============================================================

pop_18C = make_delia_population(egg_init=300.0)
n_days_3 = 200
weather_18C = WeatherSeries(fill(18.0, n_days_3); day_offset=1)

println("\n18°C cohort population:")
println("  Stages: $(n_stages(pop_18C))")
println("  Total substages: $(n_substages(pop_18C))")
println("  Initial population: $(total_population(pop_18C))")

prob_18C = PBDMProblem(pop_18C, weather_18C, (1, n_days_3))
sol_18C = solve(prob_18C, DirectIteration())

dd_per_day_18 = 18.0 - T_LOWER
println("\nConstant 18°C simulation ($(n_days_3) days, 300 eggs):")
println("  Return code: $(sol_18C.retcode)")
println("  DD/day at 18°C: $(round(dd_per_day_18, digits=1))")
for (i, sname) in enumerate(STAGE_NAMES)
    traj = stage_trajectory(sol_18C, i)
    peak = maximum(traj)
    peak_day = sol_18C.t[argmax(traj)]
    println("  $(sname): peak=$(round(peak, digits=1)) at day $(peak_day)")
end

# Expected stage transition days at 18°C
egg_end   = round(Int, DD_EGG / dd_per_day_18)
larva_end = round(Int, (DD_EGG + DD_LARVA) / dd_per_day_18)
pupa_end  = round(Int, (DD_EGG + DD_LARVA + DD_PUPA) / dd_per_day_18)
println("  Expected egg→larva:  day $(egg_end)")
println("  Expected larva→pupa: day $(larva_end)")
println("  Expected pupa→adult: day $(pupa_end)")

fig3 = Figure(size=(900, 600))
ax3 = Axis(fig3[1, 1],
    title="Delia radicum Cohort Dynamics at Constant 18°C\n300 initial eggs, $(n_days_3) days",
    xlabel="Day",
    ylabel="Population",
    xlabelsize=14, ylabelsize=14)

for (i, sname) in enumerate(STAGE_NAMES)
    traj = stage_trajectory(sol_18C, i)
    lines!(ax3, sol_18C.t, traj, linewidth=2.5, color=STAGE_COLORS[i],
           label=String(sname))
end

# Stage transition markers
vlines!(ax3, [egg_end, larva_end, pupa_end], color=:gray, linestyle=:dot, linewidth=1)

xlims!(ax3, 1, n_days_3)
ylims!(ax3, 0, nothing)
axislegend(ax3, position=:rt)

save(joinpath(figdir, "sim_constant_18C.png"), fig3, px_per_unit=2)
println("\nSaved sim_constant_18C.png")

# ============================================================
# Figure 4: Seasonal UK simulation with photoperiod
# ============================================================

n_days_4 = 365

uk_weather_days = [begin
    T  = 10.5 + 6.5 * sin(2π * (d - 80) / 365)
    pp = 12.0 + 4.0 * sin(2π * (d - 80) / 365)
    DailyWeather(T, T - 3.0, T + 3.0; photoperiod=pp)
end for d in 1:n_days_4]
weather_uk = WeatherSeries(uk_weather_days; day_offset=1)

pop_uk = make_delia_population(pupa_init=100.0)
prob_uk = PBDMProblem(DensityDependent(), pop_uk, weather_uk, (1, n_days_4);
                      p=(fecundity=FECUNDITY,))
sol_uk = solve(prob_uk, DirectIteration(); reproduction_fn=delia_reproduction)

println("\nSeasonal UK simulation (365 days, 100 overwintered pupae):")
println("  Return code: $(sol_uk.retcode)")
for (i, sname) in enumerate(STAGE_NAMES)
    traj = stage_trajectory(sol_uk, i)
    peak = maximum(traj)
    peak_day = sol_uk.t[argmax(traj)]
    println("  $(sname): peak=$(round(peak, digits=1)) at day $(peak_day)")
end
cum_dd_uk = cumulative_degree_days(sol_uk)
println("  Total DD accumulated: $(round(cum_dd_uk[end], digits=1))")
println("  Theoretical generations: $(round(cum_dd_uk[end] / (DD_EGG + DD_LARVA + DD_PUPA), digits=1))")

fig4 = Figure(size=(1000, 750))

# Top panel: stage dynamics
ax4a = Axis(fig4[1, 1],
    title="Delia radicum Seasonal Dynamics — UK (52°N)\n100 overwintered pupae, sinusoidal T (mean 10.5°C, amp 6.5°C)",
    ylabel="Population",
    xlabelsize=14, ylabelsize=14)

for (i, sname) in enumerate(STAGE_NAMES)
    traj = stage_trajectory(sol_uk, i)
    lines!(ax4a, sol_uk.t, traj, linewidth=2.5, color=STAGE_COLORS[i],
           label=String(sname))
end
axislegend(ax4a, position=:rt)
xlims!(ax4a, 1, n_days_4)
ylims!(ax4a, 0, nothing)
hidexdecorations!(ax4a)

# Bottom panel: temperature and photoperiod
ax4b = Axis(fig4[2, 1],
    xlabel="Day of year",
    ylabel="Temperature (°C)",
    xlabelsize=14, ylabelsize=14,
    yticklabelcolor=:red, ylabelcolor=:red,
    xticks=(MONTH_STARTS, MONTH_LABELS))

uk_temps = [10.5 + 6.5 * sin(2π * (d - 80) / 365) for d in 1:n_days_4]
uk_pp    = [12.0 + 4.0 * sin(2π * (d - 80) / 365) for d in 1:n_days_4]

lines!(ax4b, 1:n_days_4, uk_temps, linewidth=2, color=:red, label="Temperature")
hlines!(ax4b, [T_LOWER], color=:red, linestyle=:dash, linewidth=1)

ax4c = Axis(fig4[2, 1],
    ylabel="Photoperiod (h)",
    yaxisposition=:right,
    ylabelcolor=:purple, yticklabelcolor=:purple)
lines!(ax4c, 1:n_days_4, uk_pp, linewidth=2, color=:purple, linestyle=:dash)
hlines!(ax4c, [CRITICAL_PHOTOPERIOD], color=:purple, linestyle=:dot, linewidth=1)
hidespines!(ax4c)
hidexdecorations!(ax4c)

linkxaxes!(ax4a, ax4b)
xlims!(ax4b, 1, n_days_4)

rowsize!(fig4.layout, 1, Relative(0.65))

save(joinpath(figdir, "seasonal_uk.png"), fig4, px_per_unit=2)
println("\nSaved seasonal_uk.png")

# ============================================================
# Figure 5: Latitude comparison
# ============================================================

locations = [
    (name="UK (52°N)",       lat=52.0, T_mean=10.5, T_amp=6.5),
    (name="France (48°N)",   lat=48.0, T_mean=12.0, T_amp=8.0),
    (name="Scotland (55°N)", lat=55.0, T_mean=9.0,  T_amp=5.5),
]
loc_colors = [:steelblue, :firebrick, :forestgreen]

n_days_5 = 365

fig5 = Figure(size=(1000, 750))

ax5a = Axis(fig5[1, 1],
    title="Delia radicum — Latitude Comparison\n100 overwintered pupae, 365-day simulation",
    ylabel="Total population",
    xlabelsize=14, ylabelsize=14)

ax5b = Axis(fig5[2, 1],
    xlabel="Day of year",
    ylabel="Daily mean temperature (°C)",
    xlabelsize=14, ylabelsize=14,
    xticks=(MONTH_STARTS, MONTH_LABELS))

println("\nLatitude comparison:")
for (idx, loc) in enumerate(locations)
    loc_weather_days = [begin
        T  = loc.T_mean + loc.T_amp * sin(2π * (d - 80) / 365)
        pp = 12.0 + 4.0 * sin(2π * (d - 80) / 365)
        DailyWeather(T, T - 3.0, T + 3.0; photoperiod=pp)
    end for d in 1:n_days_5]
    loc_weather = WeatherSeries(loc_weather_days; day_offset=1)

    loc_pop = make_delia_population(pupa_init=100.0)
    loc_prob = PBDMProblem(DensityDependent(), loc_pop, loc_weather, (1, n_days_5);
                           p=(fecundity=FECUNDITY,))
    loc_sol = solve(loc_prob, DirectIteration(); reproduction_fn=delia_reproduction)

    total_pop = total_population(loc_sol)
    lines!(ax5a, loc_sol.t, total_pop, linewidth=2.5, color=loc_colors[idx],
           label=loc.name)

    loc_temps = [loc.T_mean + loc.T_amp * sin(2π * (d - 80) / 365) for d in 1:n_days_5]
    lines!(ax5b, 1:n_days_5, loc_temps, linewidth=1.5, color=loc_colors[idx],
           label=loc.name)

    cum_dd = cumulative_degree_days(loc_sol)
    peak_pop = maximum(total_pop)
    peak_day = loc_sol.t[argmax(total_pop)]
    println("  $(loc.name): peak=$(round(peak_pop, digits=1)) at day $(peak_day), " *
            "DD=$(round(cum_dd[end], digits=1)), " *
            "est. gen.=$(round(cum_dd[end] / (DD_EGG + DD_LARVA + DD_PUPA), digits=1))")
end

hlines!(ax5b, [T_LOWER], color=:gray, linestyle=:dash, linewidth=1, label="T_lower=$(T_LOWER)°C")
axislegend(ax5a, position=:rt)
axislegend(ax5b, position=:rt)
xlims!(ax5a, 1, n_days_5)
xlims!(ax5b, 1, n_days_5)
ylims!(ax5a, 0, nothing)
hidexdecorations!(ax5a)
linkxaxes!(ax5a, ax5b)

rowsize!(fig5.layout, 1, Relative(0.65))

save(joinpath(figdir, "latitude_comparison.png"), fig5, px_per_unit=2)
println("\nSaved latitude_comparison.png")

println("\n" * "=" ^ 60)
println("All Delia radicum validation figures saved to:")
println("  $(figdir)")
println("=" ^ 60)
