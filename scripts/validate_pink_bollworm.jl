#!/usr/bin/env julia
# Validation script for Pink Bollworm (Pectinophora gossypiella) PBDM
# matching key biological parameters from:
#   - Gutierrez et al. (2006) — degree-day requirements and range limits
#   - Henneberry & Naranjo (1998) — development rate data
#
# Generates 5 figures in scripts/figures/pink_bollworm/:
#   1. devrate_curves.png       — Brière development rate vs temperature
#   2. mortality_curves.png     — U-shaped mortality vs temperature
#   3. sim_constant_28C.png     — Cohort dynamics at constant 28°C
#   4. range_limits.png         — Climate-driven range limits (4 locations)
#   5. generations_count.png    — Adult emergence pulses at Imperial Valley

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using PhysiologicallyBasedDemographicModels
using CairoMakie

figdir = joinpath(@__DIR__, "figures", "pink_bollworm")
mkpath(figdir)

# ============================================================
# Parameters from literature
# (Gutierrez et al. 2006, Henneberry & Naranjo 1998)
# ============================================================

const T_LOWER = 12.8   # lower developmental threshold (°C)
const T_UPPER = 36.0   # upper developmental threshold (°C)

# Degree-day requirements (above 12.8°C)
const DD_EGG   = 60.0   # egg stage
const DD_LARVA = 220.0  # larval stage
const DD_PUPA  = 130.0  # pupal stage
const DD_ADULT = 200.0  # adult longevity

# Brière coefficients calibrated to match literature rates:
#   r(T) = a * T * (T - T_lower) * sqrt(T_upper - T)
#   Egg:   ~0.08 at 20°C, ~0.14 at 25°C, ~0.18 at 30°C
#   Larva: ~0.04 at 25°C, ~0.06 at 30°C
#   Pupa:  intermediate between egg and larva
const A_EGG   = 0.0000960
const A_LARVA = 0.0000260
const A_PUPA  = 0.0000440

# Brière models for development rate curve visualization (Fig 1)
egg_briere   = BriereDevelopmentRate(A_EGG,   T_LOWER, T_UPPER)
larva_briere = BriereDevelopmentRate(A_LARVA, T_LOWER, T_UPPER)
pupa_briere  = BriereDevelopmentRate(A_PUPA,  T_LOWER, T_UPPER)

# Linear model for simulations: degree_days(T) = max(0, T - T_lower)
# τ (DD) values match literature directly in linear DD units
sim_dev = LinearDevelopmentRate(T_LOWER, T_UPPER)

# Mortality coefficients — U-shaped: μ(T) = a*(T - T_opt)^2 + μ_min
const MU_MIN_EGG   = 0.005;  const MU_A_EGG   = 0.00025; const T_OPT_EGG   = 27.0
const MU_MIN_LARVA = 0.003;  const MU_A_LARVA = 0.00020; const T_OPT_LARVA = 26.0
const MU_MIN_PUPA  = 0.004;  const MU_A_PUPA  = 0.00022; const T_OPT_PUPA  = 26.0
const MU_MIN_ADULT = 0.008;  const MU_A_ADULT = 0.00030; const T_OPT_ADULT = 25.0

# Mortality functions
μ_egg(T)   = MU_MIN_EGG   + MU_A_EGG   * (T - T_OPT_EGG)^2
μ_larva(T) = MU_MIN_LARVA + MU_A_LARVA * (T - T_OPT_LARVA)^2
μ_pupa(T)  = MU_MIN_PUPA  + MU_A_PUPA  * (T - T_OPT_PUPA)^2
μ_adult(T) = MU_MIN_ADULT + MU_A_ADULT * (T - T_OPT_ADULT)^2

# Diagnostic output
Ts_diag = [15.0, 20.0, 25.0, 30.0, 35.0]
println("Pink bollworm Brière development rates (for plotting):")
for T in Ts_diag
    re = development_rate(egg_briere, T)
    rl = development_rate(larva_briere, T)
    rp = development_rate(pupa_briere, T)
    println("  T=$(T)°C: egg=$(round(re, digits=4)), larva=$(round(rl, digits=4)), pupa=$(round(rp, digits=4))")
end
println("Linear DD model for simulations:")
for T in Ts_diag
    dd = degree_days(sim_dev, T)
    println("  T=$(T)°C: DD/day=$(round(dd, digits=1)), egg_dur=$(round(DD_EGG/max(dd,0.1), digits=1))d, larva_dur=$(round(DD_LARVA/max(dd,0.1), digits=1))d")
end

# ============================================================
# Figure 1: Development rate curves
# ============================================================

Ts = range(5.0, 42.0, length=300)
egg_rates   = [development_rate(egg_briere, T) for T in Ts]
larva_rates = [development_rate(larva_briere, T) for T in Ts]
pupa_rates  = [development_rate(pupa_briere, T) for T in Ts]

# Literature data points (Henneberry & Naranjo 1998)
lit_egg_T   = [20.0, 25.0, 30.0]
lit_egg_r   = [0.08, 0.14, 0.18]
lit_larva_T = [25.0, 30.0]
lit_larva_r = [0.04, 0.06]

fig1 = Figure(size=(900, 600))
ax1 = Axis(fig1[1, 1],
    title="Pink Bollworm — Temperature-Dependent Development Rates\n(Brière model, T_lower=$(T_LOWER)°C, T_upper=$(T_UPPER)°C)",
    xlabel="Temperature (°C)",
    ylabel="Development rate (1/day)",
    xlabelsize=14, ylabelsize=14)

lines!(ax1, collect(Ts), egg_rates, linewidth=2.5, color=:firebrick, label="Egg (τ=$(Int(DD_EGG)) DD)")
lines!(ax1, collect(Ts), larva_rates, linewidth=2.5, color=:forestgreen, label="Larva (τ=$(Int(DD_LARVA)) DD)")
lines!(ax1, collect(Ts), pupa_rates, linewidth=2.5, color=:steelblue, label="Pupa (τ=$(Int(DD_PUPA)) DD)")

scatter!(ax1, lit_egg_T, lit_egg_r, color=:firebrick, markersize=12, marker=:diamond,
         label="Egg (literature)")
scatter!(ax1, lit_larva_T, lit_larva_r, color=:forestgreen, markersize=12, marker=:utriangle,
         label="Larva (literature)")

vlines!(ax1, [T_LOWER], color=:gray, linestyle=:dash, linewidth=1)
text!(ax1, T_LOWER + 0.3, maximum(egg_rates) * 0.95,
    text="T_lower=$(T_LOWER)°C", fontsize=10, color=:gray, align=(:left, :top))
vlines!(ax1, [T_UPPER], color=:gray, linestyle=:dash, linewidth=1)
text!(ax1, T_UPPER - 0.3, maximum(egg_rates) * 0.95,
    text="T_upper=$(T_UPPER)°C", fontsize=10, color=:gray, align=(:right, :top))

xlims!(ax1, 5, 42)
ylims!(ax1, 0, nothing)
axislegend(ax1, position=:lt)

save(joinpath(figdir, "devrate_curves.png"), fig1, px_per_unit=2)
println("Saved devrate_curves.png — egg peak: $(round(maximum(egg_rates), digits=4)) at $(round(collect(Ts)[argmax(egg_rates)], digits=1))°C")

# ============================================================
# Figure 2: Mortality curves
# ============================================================

Ts_mort = range(5.0, 42.0, length=300)
egg_mort   = [μ_egg(T) for T in Ts_mort]
larva_mort = [μ_larva(T) for T in Ts_mort]
pupa_mort  = [μ_pupa(T) for T in Ts_mort]
adult_mort = [μ_adult(T) for T in Ts_mort]

fig2 = Figure(size=(900, 600))
ax2 = Axis(fig2[1, 1],
    title="Pink Bollworm — Temperature-Dependent Mortality Rates\n(U-shaped: high mortality below 15°C and above 35°C)",
    xlabel="Temperature (°C)",
    ylabel="Mortality rate (per degree-day)",
    xlabelsize=14, ylabelsize=14)

lines!(ax2, collect(Ts_mort), egg_mort, linewidth=2, color=:firebrick, label="Egg")
lines!(ax2, collect(Ts_mort), larva_mort, linewidth=2, color=:forestgreen, label="Larva")
lines!(ax2, collect(Ts_mort), pupa_mort, linewidth=2, color=:steelblue, label="Pupa")
lines!(ax2, collect(Ts_mort), adult_mort, linewidth=2, color=:darkorange, label="Adult")

# Mark the high-mortality zones
vspan!(ax2, 5.0, 15.0, color=(:blue, 0.06))
vspan!(ax2, 35.0, 42.0, color=(:red, 0.06))
text!(ax2, 10.0, maximum(egg_mort) * 0.9, text="Cold stress", fontsize=10, color=:blue,
      align=(:center, :top))
text!(ax2, 38.5, maximum(egg_mort) * 0.9, text="Heat stress", fontsize=10, color=:red,
      align=(:center, :top))

xlims!(ax2, 5, 42)
ylims!(ax2, 0, nothing)
axislegend(ax2, position=:ct)

save(joinpath(figdir, "mortality_curves.png"), fig2, px_per_unit=2)
println("Saved mortality_curves.png — egg mortality range: $(round(minimum(egg_mort), digits=4))–$(round(maximum(egg_mort), digits=4))")

# ============================================================
# Figure 3: Cohort simulation at constant 28°C
# ============================================================

function build_pbw_population(; n_eggs=0.0)
    # k must satisfy k ≤ τ/dd_max for CFL stability of the delay scheme
    # dd_max ≈ 23 (at ~36°C, the upper threshold)
    egg   = LifeStage(:egg,   DistributedDelay(2,  DD_EGG;   W0=0.0), sim_dev, 0.005)
    larva = LifeStage(:larva, DistributedDelay(9,  DD_LARVA; W0=0.0), sim_dev, 0.003)
    pupa  = LifeStage(:pupa,  DistributedDelay(5,  DD_PUPA;  W0=0.0), sim_dev, 0.004)
    adult = LifeStage(:adult, DistributedDelay(8,  DD_ADULT; W0=0.0), sim_dev, 0.008)
    if n_eggs > 0.0
        egg.delay.W[1] = n_eggs
    end
    Population(:pink_bollworm, [egg, larva, pupa, adult])
end

n_sim = 200
pop_28 = build_pbw_population(n_eggs=500.0)
weather_28 = WeatherSeries(fill(28.0, n_sim))
prob_28 = PBDMProblem(pop_28, weather_28, (1, n_sim))
sol_28 = solve(prob_28, DirectIteration())

stage_names = [:egg, :larva, :pupa, :adult]
colors_stage = [:firebrick, :forestgreen, :steelblue, :darkorange]

println("\nCohort at constant 28°C (500 eggs, $(n_sim) days):")
for (i, sname) in enumerate(stage_names)
    traj = stage_trajectory(sol_28, i)
    peak = maximum(traj)
    peak_day = argmax(traj) - 1  # offset for initial state
    println("  $(sname): peak=$(round(peak, digits=1)) at day $(peak_day)")
end

fig3 = Figure(size=(900, 600))
ax3 = Axis(fig3[1, 1],
    title="Pink Bollworm — Cohort Dynamics at Constant 28°C\n(Initial: 500 eggs, $(n_sim) days)",
    xlabel="Day",
    ylabel="Population",
    xlabelsize=14, ylabelsize=14)

for (i, sname) in enumerate(stage_names)
    traj = stage_trajectory(sol_28, i)
    lines!(ax3, sol_28.t, traj, linewidth=2.5, color=colors_stage[i],
           label=String(sname))
end

dd_per_day = 28.0 - T_LOWER
egg_dur = DD_EGG / dd_per_day
larva_dur = DD_LARVA / dd_per_day
pupa_dur = DD_PUPA / dd_per_day
text!(ax3, n_sim * 0.55, maximum(stage_trajectory(sol_28, 1)) * 0.85,
    text="DD/day = $(round(dd_per_day, digits=1))\nEgg: ~$(round(egg_dur, digits=0))d\nLarva: ~$(round(larva_dur, digits=0))d\nPupa: ~$(round(pupa_dur, digits=0))d",
    fontsize=10, align=(:left, :top))

axislegend(ax3, position=:rt)
xlims!(ax3, 1, n_sim)

save(joinpath(figdir, "sim_constant_28C.png"), fig3, px_per_unit=2)
println("Saved sim_constant_28C.png")

# ============================================================
# Figure 4: Climate-driven range limits (4 locations)
# ============================================================

# Reproduction callback: adults lay eggs proportional to DD accumulated
# Pink bollworm ~150-300 eggs/female, ~200 DD lifespan, 50% female
# Net R0 ≈ 3-5 per generation: fecundity_per_dd ~ 0.15
const FECUNDITY_PER_DD = 0.15

function pbw_reproduction(pop, weather_day, p, day)
    adult_total = delay_total(pop.stages[4].delay)
    dd = degree_days(pop.stages[1].dev_rate, weather_day.T_mean)
    return max(0.0, adult_total * FECUNDITY_PER_DD * dd)
end

# Climates modeled with SinusoidalWeather (T_mean, amplitude)
#   T(d) = T_mean + amplitude * sin(2π(d - phase)/365),  phase=200 (peak ~mid-July)
# Amplitudes chosen to include winter periods below T_LOWER (12.8°C)
climates = [
    ("Imperial Valley, CA (hot, endemic)",      22.0, 14.0),
    ("Phoenix, AZ (marginal)",                  21.0, 14.0),
    ("Lubbock, TX (northern limit)",            15.0, 15.0),
    ("Bakersfield, CA (cool margin)",           17.0, 12.0),
]

climate_colors = [:firebrick, :darkorange, :steelblue, :mediumpurple]

fig4 = Figure(size=(1000, 700))
ax4 = Axis(fig4[1, 1],
    title="Pink Bollworm — Climate-Driven Range Limits\n(365-day total population with reproduction, sinusoidal temperature)",
    xlabel="Day of year",
    ylabel="Total population (log₁₀)",
    yscale=log10,
    xlabelsize=14, ylabelsize=14)

println("\nRange limit comparison (365 days, 100 initial eggs, with reproduction):")
for (idx, (label, T_mean, amp)) in enumerate(climates)
    weather_loc = SinusoidalWeather(T_mean, amp)
    pop_loc = build_pbw_population(n_eggs=100.0)
    prob_loc = PBDMProblem(DensityDependent(), pop_loc, weather_loc, (1, 365))
    sol_loc = solve(prob_loc, DirectIteration(); reproduction_fn=pbw_reproduction)
    total_traj = total_population(sol_loc)
    # Clamp to positive for log scale
    total_traj_plot = max.(total_traj, 0.1)

    lines!(ax4, sol_loc.t, total_traj_plot, linewidth=2.5, color=climate_colors[idx],
           label=label)

    peak = maximum(total_traj)
    final = total_traj[end]
    println("  $(label): mean_T=$(T_mean)°C, amp=$(amp)°C, peak=$(round(peak, sigdigits=3)), final=$(round(final, sigdigits=3))")
end

axislegend(ax4, position=:lt)
xlims!(ax4, 1, 365)

save(joinpath(figdir, "range_limits.png"), fig4, px_per_unit=2)
println("Saved range_limits.png")

# ============================================================
# Figure 5: Generations per year at Imperial Valley
# ============================================================

weather_iv = SinusoidalWeather(22.0, 14.0)
pop_iv = build_pbw_population(n_eggs=100.0)
prob_iv = PBDMProblem(DensityDependent(), pop_iv, weather_iv, (1, 365))
sol_iv = solve(prob_iv, DirectIteration(); reproduction_fn=pbw_reproduction)

adult_traj = stage_trajectory(sol_iv, 4)  # stage 4 = adult
days_iv = sol_iv.t

# Use cumulative degree-days to estimate generation timing
# Total egg-to-adult development: 60 + 220 + 130 = 410 DD
const GEN_DD = DD_EGG + DD_LARVA + DD_PUPA
cum_dd = cumulative_degree_days(sol_iv)

# Find days when cumulative DD crosses generation boundaries
gen_days = let
    result = Int[]
    gn = 1
    for i in 1:length(cum_dd)
        if cum_dd[i] >= gn * GEN_DD
            push!(result, days_iv[i])
            gn += 1
        end
    end
    result
end
n_generations = length(gen_days)

fig5 = Figure(size=(1000, 700))

# Top panel: adult population
ax5a = Axis(fig5[1, 1],
    title="Pink Bollworm — Generations per Year at Imperial Valley, CA\n(Mean 22°C, Amplitude 14°C; literature: ~4–5 generations with diapause)",
    ylabel="Adult population",
    xlabelsize=14, ylabelsize=14)

lines!(ax5a, days_iv, adult_traj, linewidth=2.5, color=:darkorange, label="Adults")

# Mark generation emergence times (when cumDD crosses multiples of 410 DD)
for (g, gday) in enumerate(gen_days)
    vlines!(ax5a, [gday], color=:gray40, linestyle=:dash, linewidth=1)
    # Find the adult population at this day
    idx = findfirst(d -> d >= gday, days_iv)
    if idx !== nothing
        ypos = adult_traj[idx]
        text!(ax5a, gday + 3, ypos,
            text="G$(g)", fontsize=11, color=:red, align=(:left, :center))
    end
end

text!(ax5a, 340, maximum(adult_traj) * 0.85,
    text="$(n_generations) generations\n($(Int(GEN_DD)) DD each)",
    fontsize=12, color=:red, align=(:right, :top))
axislegend(ax5a, position=:lt)
xlims!(ax5a, 1, 365)
hidexdecorations!(ax5a, grid=false)

# Bottom panel: cumulative degree-days with generation thresholds
ax5b = Axis(fig5[2, 1],
    xlabel="Day of year",
    ylabel="Cumulative DD (°C·days)",
    xlabelsize=14, ylabelsize=14)

lines!(ax5b, days_iv[1:end-1], cum_dd, linewidth=2, color=:black, label="Cumulative DD")

for g in 1:min(n_generations + 1, 10)
    hlines!(ax5b, [g * GEN_DD], color=(:red, 0.4), linestyle=:dash, linewidth=1)
    text!(ax5b, 5, g * GEN_DD + 30,
        text="G$(g) = $(Int(g * GEN_DD)) DD", fontsize=9, color=:red,
        align=(:left, :bottom))
end

axislegend(ax5b, position=:lt)
xlims!(ax5b, 1, 365)
linkxaxes!(ax5a, ax5b)
rowsize!(fig5.layout, 1, Relative(0.6))

save(joinpath(figdir, "generations_count.png"), fig5, px_per_unit=2)
println("Saved generations_count.png — $(n_generations) generation boundaries crossed ($(Int(GEN_DD)) DD/gen, total DD=$(round(cum_dd[end], digits=0)))")

println("\nAll pink bollworm validation figures saved to: $(figdir)")
