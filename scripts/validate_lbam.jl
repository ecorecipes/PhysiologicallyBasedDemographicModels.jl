#!/usr/bin/env julia
# Validation script for Light Brown Apple Moth (Epiphyas postvittana) PBDM
# matching literature values from Danthanarayana (1975) and Gutierrez et al.
#
# Generates 5 figures in scripts/figures/lbam/:
#   1. devrate_curves.png        — Brière development rate curves with observed data
#   2. mortality_curves.png      — U-shaped temperature-dependent mortality
#   3. sim_constant_20C.png      — Cohort dynamics at constant 20°C
#   4. three_climates.png        — Total population under three climatic regimes
#   5. generations_panel.png     — Adult emergence pulses (generation counting)

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using PhysiologicallyBasedDemographicModels
using CairoMakie

figdir = joinpath(@__DIR__, "figures", "lbam")
mkpath(figdir)

# ============================================================
# Species parameters — Epiphyas postvittana
# Literature: Danthanarayana (1975), Gutierrez et al. (2010)
# ============================================================

# Developmental thresholds (°C)
const T_LOWER = 7.5
const T_UPPER = 31.5

# Degree-day requirements per stage (DD above 7.5°C)
# At 20°C: DD/day = 12.5, so egg ~10d, larva ~28d, pupa ~10d
# Total immature ≈ 600 DD (literature: 526–570 DD)
const DD_EGG   = 125.0
const DD_LARVA = 350.0
const DD_PUPA  = 125.0
const DD_ADULT = 300.0   # adult lifespan in DD (~24 days at 20°C)

# Substage counts chosen for numerical stability: k/τ * DD_max < 1
# where DD_max ≈ 27.5 (at 35°C = max Phoenix temp)
const K_EGG   = 4
const K_LARVA = 12
const K_PUPA  = 4
const K_ADULT = 10

# Linear development rate for simulations (degree-day accumulation)
const lbam_lindev = LinearDevelopmentRate(T_LOWER, T_UPPER)

# Brière coefficients calibrated to Danthanarayana (1975) observed rates
# r(T) = a · T · (T − T_lower) · √(T_upper − T) for T_lower < T < T_upper
const A_EGG   = 1.15e-4
const A_LARVA = 3.50e-5
const A_PUPA  = 1.10e-4

egg_briere   = BriereDevelopmentRate(A_EGG,   T_LOWER, T_UPPER)
larva_briere = BriereDevelopmentRate(A_LARVA, T_LOWER, T_UPPER)
pupa_briere  = BriereDevelopmentRate(A_PUPA,  T_LOWER, T_UPPER)

# Mortality coefficients — U-shaped quadratic: μ(T) = max(0, a·T² − b·T + c)
const μ_EGG_A = 0.00035;  const μ_EGG_B = 0.0140;  const μ_EGG_C = 0.155
const μ_LAR_A = 0.00040;  const μ_LAR_B = 0.0165;  const μ_LAR_C = 0.185
const μ_PUP_A = 0.00030;  const μ_PUP_B = 0.0120;  const μ_PUP_C = 0.130
const μ_ADU_A = 0.00038;  const μ_ADU_B = 0.0155;  const μ_ADU_C = 0.170

μ_egg(T)   = max(0.0, μ_EGG_A * T^2 - μ_EGG_B * T + μ_EGG_C)
μ_larva(T) = max(0.0, μ_LAR_A * T^2 - μ_LAR_B * T + μ_LAR_C)
μ_pupa(T)  = max(0.0, μ_PUP_A * T^2 - μ_PUP_B * T + μ_PUP_C)
μ_adult(T) = max(0.0, μ_ADU_A * T^2 - μ_ADU_B * T + μ_ADU_C)

# Temperature range for curves
Ts = range(0.0, 40.0, length=300)

# Observed data from Danthanarayana (1975) — development rate (1/day)
const OBS_TEMPS = [15.0, 20.0, 25.0]
const OBS_EGG   = [0.05, 0.10, 0.16]
const OBS_LARVA = [0.015, 0.03, 0.05]
const OBS_PUPA  = [0.05, 0.10, 0.15]

# Shared plot styling
const STAGE_NAMES  = [:egg, :larva, :pupa, :adult]
const STAGE_COLORS = [:goldenrod, :forestgreen, :steelblue, :firebrick]

# ============================================================
# Figure 1: Brière development rate curves with observed data
# ============================================================

egg_rates   = [development_rate(egg_briere, T)   for T in Ts]
larva_rates = [development_rate(larva_briere, T) for T in Ts]
pupa_rates  = [development_rate(pupa_briere, T)  for T in Ts]

fig1 = Figure(size=(900, 600))
ax1 = Axis(fig1[1, 1],
    title="LBAM (E. postvittana) — Brière Development Rate Curves",
    xlabel="Temperature (°C)",
    ylabel="Development rate (1/day)",
    xlabelsize=14, ylabelsize=14)

lines!(ax1, collect(Ts), egg_rates, linewidth=2.5, color=:goldenrod, label="Egg")
lines!(ax1, collect(Ts), larva_rates, linewidth=2.5, color=:forestgreen, label="Larva")
lines!(ax1, collect(Ts), pupa_rates, linewidth=2.5, color=:steelblue, label="Pupa")

scatter!(ax1, OBS_TEMPS, OBS_EGG, color=:goldenrod, markersize=12, marker=:circle,
         strokewidth=1.5, strokecolor=:black, label="Egg (Danthanarayana 1975)")
scatter!(ax1, OBS_TEMPS, OBS_LARVA, color=:forestgreen, markersize=12, marker=:utriangle,
         strokewidth=1.5, strokecolor=:black, label="Larva (obs.)")
scatter!(ax1, OBS_TEMPS, OBS_PUPA, color=:steelblue, markersize=12, marker=:diamond,
         strokewidth=1.5, strokecolor=:black, label="Pupa (obs.)")

xlims!(ax1, 0, 40)
ylims!(ax1, 0, nothing)
axislegend(ax1, position=:lt)

save(joinpath(figdir, "devrate_curves.png"), fig1, px_per_unit=2)
println("Saved devrate_curves.png")
println("  Egg peak:   $(round(maximum(egg_rates), digits=4)) 1/day")
println("  Larva peak: $(round(maximum(larva_rates), digits=4)) 1/day")
println("  Pupa peak:  $(round(maximum(pupa_rates), digits=4)) 1/day")

# ============================================================
# Figure 2: Temperature-dependent mortality curves (U-shaped)
# ============================================================

egg_mort   = [μ_egg(T)   for T in Ts]
larva_mort = [μ_larva(T) for T in Ts]
pupa_mort  = [μ_pupa(T)  for T in Ts]
adult_mort = [μ_adult(T) for T in Ts]

fig2 = Figure(size=(900, 600))
ax2 = Axis(fig2[1, 1],
    title="LBAM — Temperature-Dependent Mortality Rates (U-shaped quadratic)",
    xlabel="Temperature (°C)",
    ylabel="Mortality rate (per DD)",
    xlabelsize=14, ylabelsize=14)

lines!(ax2, collect(Ts), egg_mort,   linewidth=2, color=:goldenrod,   label="Egg")
lines!(ax2, collect(Ts), larva_mort, linewidth=2, color=:forestgreen, label="Larva")
lines!(ax2, collect(Ts), pupa_mort,  linewidth=2, color=:steelblue,   label="Pupa")
lines!(ax2, collect(Ts), adult_mort, linewidth=2, color=:firebrick,   label="Adult")

T_opt_egg   = μ_EGG_B / (2 * μ_EGG_A)
T_opt_larva = μ_LAR_B / (2 * μ_LAR_A)
T_opt_pupa  = μ_PUP_B / (2 * μ_PUP_A)
T_opt_adult = μ_ADU_B / (2 * μ_ADU_A)

vlines!(ax2, [T_opt_egg], color=(:goldenrod, 0.4), linestyle=:dash, linewidth=1)
vlines!(ax2, [T_opt_larva], color=(:forestgreen, 0.4), linestyle=:dash, linewidth=1)

xlims!(ax2, 0, 40)
ylims!(ax2, 0, 0.30)
axislegend(ax2, position=:rt)

save(joinpath(figdir, "mortality_curves.png"), fig2, px_per_unit=2)
println("\nSaved mortality_curves.png")
println("  Optimal T (min mortality): egg=$(round(T_opt_egg, digits=1))°C, " *
        "larva=$(round(T_opt_larva, digits=1))°C, " *
        "pupa=$(round(T_opt_pupa, digits=1))°C, " *
        "adult=$(round(T_opt_adult, digits=1))°C")

# ============================================================
# Figure 3: Cohort at constant 20°C, 200 days — 500 eggs
# Uses LinearDevelopmentRate so DD units match DistributedDelay τ
# ============================================================

n_sim = 200

lbam_stages_20 = [
    LifeStage(:egg,   DistributedDelay(K_EGG,   DD_EGG;   W0=125.0), lbam_lindev, 0.001),
    LifeStage(:larva, DistributedDelay(K_LARVA, DD_LARVA; W0=0.0),   lbam_lindev, 0.0008),
    LifeStage(:pupa,  DistributedDelay(K_PUPA,  DD_PUPA;  W0=0.0),   lbam_lindev, 0.0005),
    LifeStage(:adult, DistributedDelay(K_ADULT, DD_ADULT; W0=0.0),   lbam_lindev, 0.0005),
]
lbam_20 = Population(:lbam, lbam_stages_20)

weather_20 = WeatherSeries([DailyWeather(20.0) for _ in 1:n_sim])
prob_20 = PBDMProblem(lbam_20, weather_20, (1, n_sim + 1))
sol_20 = solve(prob_20, DirectIteration())

dd_per_day = 20.0 - T_LOWER
println("\nConstant 20°C cohort ($(n_sim) days, 500 eggs):")
println("  DD/day: $(dd_per_day),  expected egg: $(round(DD_EGG/dd_per_day, digits=1))d, " *
        "larva: $(round(DD_LARVA/dd_per_day, digits=1))d, pupa: $(round(DD_PUPA/dd_per_day, digits=1))d")
for (i, sname) in enumerate(STAGE_NAMES)
    traj = sol_20.stage_totals[i, 2:end]
    peak_val = maximum(traj)
    peak_day = argmax(traj)
    println("  $(sname): peak=$(round(peak_val, digits=1)) at day $(peak_day)")
end

fig3 = Figure(size=(900, 600))
ax3 = Axis(fig3[1, 1],
    title="LBAM Cohort at Constant 20°C (500 eggs initial)",
    xlabel="Day", ylabel="Population",
    xlabelsize=14, ylabelsize=14)

for (i, sname) in enumerate(STAGE_NAMES)
    traj = sol_20.stage_totals[i, 2:end]
    lines!(ax3, 1:n_sim, traj, linewidth=2.5, color=STAGE_COLORS[i], label=String(sname))
end
xlims!(ax3, 1, n_sim)
ylims!(ax3, 0, nothing)
axislegend(ax3, position=:rt)

save(joinpath(figdir, "sim_constant_20C.png"), fig3, px_per_unit=2)
println("Saved sim_constant_20C.png")

# ============================================================
# Helpers for climate simulations
# ============================================================

function make_lbam_population(; n_adults=50.0, μ_bg=0.002)
    n_adult_per_sub = n_adults / K_ADULT
    stages = [
        LifeStage(:egg,   DistributedDelay(K_EGG,   DD_EGG;   W0=0.0),             lbam_lindev, μ_bg),
        LifeStage(:larva, DistributedDelay(K_LARVA, DD_LARVA; W0=0.0),             lbam_lindev, μ_bg),
        LifeStage(:pupa,  DistributedDelay(K_PUPA,  DD_PUPA;  W0=0.0),             lbam_lindev, μ_bg),
        LifeStage(:adult, DistributedDelay(K_ADULT, DD_ADULT; W0=n_adult_per_sub), lbam_lindev, μ_bg),
    ]
    Population(:lbam, stages)
end

function make_annual_weather(T_mean, amplitude; n_days=365, phase=200.0)
    temps = [T_mean + amplitude * sin(2π * (d - phase) / 365) for d in 1:n_days]
    WeatherSeries([DailyWeather(T) for T in temps])
end

# Reproduction callback: density-regulated with hard cap
function lbam_reproduction(pop, w, p, day)
    adult_total = delay_total(pop.stages[end].delay)
    T = w.T_mean
    (T < T_LOWER || T > T_UPPER) && return 0.0
    total_pop = total_population(pop)
    total_pop >= p.K && return 0.0
    fecundity = p.f * (1.0 - total_pop / p.K)
    return min(fecundity * adult_total, p.max_eggs)
end

const REPRO_PARAMS = (f=0.5, K=1000.0, max_eggs=50.0)

# ============================================================
# Figure 4: Three climates — total population comparison
# ============================================================

climates = [
    ("San Francisco", 14.0, 3.0),
    ("Sacramento",    16.0, 9.0),
    ("Phoenix",       23.0, 12.0),
]
climate_colors = [:steelblue, :forestgreen, :firebrick]
n_days_annual = 365

fig4 = Figure(size=(900, 600))
ax4 = Axis(fig4[1, 1],
    title="LBAM Total Population — Three Climatic Regimes (365 days)",
    xlabel="Day of year", ylabel="Total population",
    xlabelsize=14, ylabelsize=14)

for (idx, (name, T_mean, amp)) in enumerate(climates)
    pop = make_lbam_population(n_adults=50.0)
    weather = make_annual_weather(T_mean, amp; n_days=n_days_annual)
    prob = PBDMProblem(DensityDependent(), pop, weather, (1, n_days_annual + 1);
                       p=REPRO_PARAMS)
    sol = solve(prob, DirectIteration(); reproduction_fn=lbam_reproduction)

    total_pop = [sum(sol.stage_totals[:, d]) for d in 2:(n_days_annual + 1)]
    lines!(ax4, 1:n_days_annual, total_pop, linewidth=2.5, color=climate_colors[idx],
           label="$(name) ($(T_mean)°C ± $(amp)°C)")

    println("\n$(name) (mean=$(T_mean)°C, amp=$(amp)°C):")
    println("  Peak total pop: $(round(maximum(total_pop), digits=1)) at day $(argmax(total_pop))")
    println("  Final total pop: $(round(total_pop[end], digits=1))")
end

xlims!(ax4, 1, n_days_annual)
ylims!(ax4, 0, nothing)
axislegend(ax4, position=:lt)

save(joinpath(figdir, "three_climates.png"), fig4, px_per_unit=2)
println("\nSaved three_climates.png")

# ============================================================
# Figure 5: Generation counting — adult emergence pulses
# ============================================================

fig5 = Figure(size=(1000, 700))
panel_climates = [("San Francisco", 14.0, 3.0), ("Sacramento", 16.0, 9.0)]

for (row, (name, T_mean, amp)) in enumerate(panel_climates)
    pop = make_lbam_population(n_adults=50.0)
    weather = make_annual_weather(T_mean, amp; n_days=n_days_annual)
    prob = PBDMProblem(DensityDependent(), pop, weather, (1, n_days_annual + 1);
                       p=REPRO_PARAMS)
    sol = solve(prob, DirectIteration(); reproduction_fn=lbam_reproduction)

    adult_traj = sol.stage_totals[4, 2:end]  # stage 4 = adult

    ax = Axis(fig5[row, 1],
        title="$(name) — Adult Population (mean $(T_mean)°C ± $(amp)°C)",
        xlabel=row == 2 ? "Day of year" : "",
        ylabel="Adult population",
        xlabelsize=14, ylabelsize=14)

    lines!(ax, 1:n_days_annual, adult_traj, linewidth=2, color=climate_colors[row])

    # Detect generation peaks via local maxima with minimum separation
    peaks_days = Int[]
    for d in 3:(n_days_annual - 1)
        if adult_traj[d] > adult_traj[d-1] && adult_traj[d] > adult_traj[d+1] &&
           adult_traj[d] > 0.05 * maximum(adult_traj)
            if isempty(peaks_days) || (d - peaks_days[end]) > 30
                push!(peaks_days, d)
            elseif adult_traj[d] > adult_traj[peaks_days[end]]
                peaks_days[end] = d
            end
        end
    end

    if !isempty(peaks_days)
        scatter!(ax, peaks_days, adult_traj[peaks_days],
                 color=:red, markersize=10, marker=:star5)
        for (gi, pd) in enumerate(peaks_days)
            text!(ax, pd, adult_traj[pd],
                text="G$(gi)", align=(:center, :bottom), fontsize=11,
                color=:red, offset=(0, 5))
        end
    end

    xlims!(ax, 1, n_days_annual)
    ylims!(ax, 0, nothing)

    n_gen = length(peaks_days)
    println("\n$(name): detected $(n_gen) adult emergence peaks (generations)")
    for (gi, pd) in enumerate(peaks_days)
        println("  G$(gi): day $(pd), adults=$(round(adult_traj[pd], digits=1))")
    end
end

save(joinpath(figdir, "generations_panel.png"), fig5, px_per_unit=2)
println("\nSaved generations_panel.png")
println("\nAll LBAM validation figures saved to: $(figdir)")
