#!/usr/bin/env julia
# Validation script for the Spotted Alfalfa Aphid (Therioaphis maculata) vignette.
#
# Generates five figures in scripts/figures/alfalfa_aphid/:
#   1. devrate_curves.png        — Brière development rate curves (nymph & adult)
#   2. mortality_curves.png      — U-shaped temperature-dependent mortality
#   3. sim_constant_25C.png      — Cohort dynamics at constant 25°C
#   4. factorial_biocontrol.png  — Multi-factor suppression bar chart
#   5. seasonal_california.png   — Central Valley CA 365-day seasonal simulation

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using PhysiologicallyBasedDemographicModels
using CairoMakie

figdir = joinpath(@__DIR__, "figures", "alfalfa_aphid")
mkpath(figdir)

# ============================================================
# Species parameters — Therioaphis maculata
# Literature: Kindler & Spomer (1986), Berberet et al. (2009),
#             Summers et al. (1959), Gutierrez & Ponti (2013)
# ============================================================

# Thermal thresholds
const T_LOWER       = 6.0    # Lower developmental threshold (°C)
const T_UPPER_NYMPH = 32.0   # Upper thermal limit — nymphs
const T_UPPER_ADULT = 33.0   # Upper thermal limit — adults

# Brière development rate models
# r(T) = a · T · (T − T_lower) · √(T_upper − T)
nymph_dev = BriereDevelopmentRate(0.0001, T_LOWER, T_UPPER_NYMPH)
adult_dev = BriereDevelopmentRate(0.000040, T_LOWER, T_UPPER_ADULT)

# Degree-day durations (in Brière-DD units, where DD = development_rate(model, T))
# Cohort simulations (figs 3–4): larger τ so dynamics span ~60–100 days at 25°C,
# giving mortality time to accumulate visibly across biocontrol scenarios.
const DD_NYMPH_COHORT = 8.0   # nymph ≈ 64 days at 25°C
const DD_ADULT_COHORT = 5.0   # adult ≈ 93 days at 25°C

# Seasonal simulation (fig 5): smaller τ for fast generations enabling
# multiple population cycles within one year.
const DD_NYMPH_SEASONAL = 1.5  # nymph ≈ 12 days at 25°C
const DD_ADULT_SEASONAL = 1.0  # adult ≈ 19 days at 25°C

# Background mortality (per degree-day)
const μ_NYMPH = 0.04
const μ_ADULT = 0.03

# U-shaped mortality curves (quadratic, for plotting)
# μ(T) = a·T² − b·T + c   (minimum near optimal temperature)
const μ_N_A = 0.00045;  const μ_N_B = 0.0225;  const μ_N_C = 0.34
const μ_A_A = 0.00035;  const μ_A_B = 0.0175;  const μ_A_C = 0.26

μ_nymph_T(T) = max(0.0, μ_N_A * T^2 - μ_N_B * T + μ_N_C)
μ_adult_T(T) = max(0.0, μ_A_A * T^2 - μ_A_B * T + μ_A_C)

# Published nymph development rate data (T°C, rate 1/day)
const PUBLISHED_T    = [15.0, 20.0, 25.0, 30.0]
const PUBLISHED_RATE = [0.05, 0.09, 0.13, 0.10]

# Print parameter diagnostics
println("=" ^ 60)
println("Therioaphis maculata — Spotted Alfalfa Aphid")
println("=" ^ 60)
println("\nNymph Brière: a=0.0001, T_lower=$(T_LOWER)°C, T_upper=$(T_UPPER_NYMPH)°C")
println("Adult Brière: a=0.000040, T_lower=$(T_LOWER)°C, T_upper=$(T_UPPER_ADULT)°C")
for T in [15.0, 20.0, 25.0, 30.0]
    rn = development_rate(nymph_dev, T)
    ra = development_rate(adult_dev, T)
    println("  T=$(T)°C: nymph r=$(round(rn, digits=4)), adult r=$(round(ra, digits=4))")
end

dd25_nymph = development_rate(nymph_dev, 25.0)
dd25_adult = development_rate(adult_dev, 25.0)
println("\nAt 25°C (cohort τ):")
println("  Nymph DD/day = $(round(dd25_nymph, digits=4)), duration = $(round(DD_NYMPH_COHORT / dd25_nymph, digits=1)) days")
println("  Adult DD/day = $(round(dd25_adult, digits=4)), duration = $(round(DD_ADULT_COHORT / dd25_adult, digits=1)) days")
println("At 25°C (seasonal τ):")
println("  Nymph duration = $(round(DD_NYMPH_SEASONAL / dd25_nymph, digits=1)) days")
println("  Adult duration = $(round(DD_ADULT_SEASONAL / dd25_adult, digits=1)) days")

# ============================================================
# Figure 1: Development rate curves for nymph and adult stages
# ============================================================

Ts = range(0.0, 40.0, length=300)
nymph_rates = [development_rate(nymph_dev, T) for T in Ts]
adult_rates = [development_rate(adult_dev, T) for T in Ts]

fig1 = Figure(size=(900, 550))
ax1 = Axis(fig1[1, 1],
    title="Therioaphis maculata — Development Rate vs Temperature",
    xlabel="Temperature (°C)",
    ylabel="Development rate (1/day)",
    xlabelsize=14, ylabelsize=14)

lines!(ax1, collect(Ts), nymph_rates, linewidth=2.5, color=:steelblue,
       label="Nymph (a=0.0001, T_u=$(T_UPPER_NYMPH)°C)")
lines!(ax1, collect(Ts), adult_rates, linewidth=2.5, color=:firebrick,
       label="Adult (a=0.000040, T_u=$(T_UPPER_ADULT)°C)")

scatter!(ax1, PUBLISHED_T, PUBLISHED_RATE, color=:black, markersize=12,
         marker=:diamond, label="Published nymph data")

# Annotate optimal range
vspan!(ax1, 26.0, 29.0, color=(:green, 0.08))
text!(ax1, 27.5, maximum(nymph_rates) * 1.08,
      text="optimal\nrange", align=(:center, :bottom), fontsize=9, color=:green)

xlims!(ax1, 0, 40)
ylims!(ax1, 0, nothing)
axislegend(ax1, position=:lt)

save(joinpath(figdir, "devrate_curves.png"), fig1, px_per_unit=2)
println("\nSaved devrate_curves.png — nymph peak: $(round(maximum(nymph_rates), digits=4)), " *
        "adult peak: $(round(maximum(adult_rates), digits=4))")

# ============================================================
# Figure 2: Temperature-dependent mortality (U-shaped)
# ============================================================

Ts_mort = range(0.0, 42.0, length=300)
nymph_mort = [μ_nymph_T(T) for T in Ts_mort]
adult_mort = [μ_adult_T(T) for T in Ts_mort]

# Optimal mortality temperatures (minimum of U-curve)
T_opt_nymph = μ_N_B / (2 * μ_N_A)
T_opt_adult = μ_A_B / (2 * μ_A_A)

fig2 = Figure(size=(900, 550))
ax2 = Axis(fig2[1, 1],
    title="Therioaphis maculata — Temperature-Dependent Mortality",
    xlabel="Temperature (°C)",
    ylabel="Mortality rate (per day)",
    xlabelsize=14, ylabelsize=14)

lines!(ax2, collect(Ts_mort), nymph_mort, linewidth=2.5, color=:steelblue,
       label="Nymph")
lines!(ax2, collect(Ts_mort), adult_mort, linewidth=2.5, color=:firebrick,
       label="Adult")

# Mark the heat-sensitivity zone
vspan!(ax2, 32.0, 42.0, color=(:red, 0.10))
text!(ax2, 37.0, maximum(nymph_mort) * 0.7,
      text="Heat\nstress\nzone", align=(:center, :center), fontsize=10, color=:red)

# Annotate minima
scatter!(ax2, [T_opt_nymph], [μ_nymph_T(T_opt_nymph)], color=:steelblue, markersize=10)
scatter!(ax2, [T_opt_adult], [μ_adult_T(T_opt_adult)], color=:firebrick, markersize=10)
text!(ax2, T_opt_nymph + 1.0, μ_nymph_T(T_opt_nymph),
      text="min $(round(T_opt_nymph, digits=1))°C",
      align=(:left, :center), fontsize=9, color=:steelblue)
text!(ax2, T_opt_adult + 1.0, μ_adult_T(T_opt_adult),
      text="min $(round(T_opt_adult, digits=1))°C",
      align=(:left, :center), fontsize=9, color=:firebrick)

xlims!(ax2, 0, 42)
ylims!(ax2, 0, nothing)
axislegend(ax2, position=:rt)

save(joinpath(figdir, "mortality_curves.png"), fig2, px_per_unit=2)
println("Saved mortality_curves.png — nymph min mortality at $(round(T_opt_nymph, digits=1))°C, " *
        "adult min at $(round(T_opt_adult, digits=1))°C")

# ============================================================
# Figure 3: Cohort simulation at constant 25°C, 150 days
# ============================================================

n_sim = 150

stages_25 = [
    LifeStage(:nymph, DistributedDelay(6, DD_NYMPH_COHORT; W0=200.0), nymph_dev, μ_NYMPH),
    LifeStage(:adult, DistributedDelay(15, DD_ADULT_COHORT; W0=0.0),  adult_dev, μ_ADULT),
]
pop_25 = Population(:alfalfa_aphid, stages_25)

weather_25 = WeatherSeries(fill(25.0, n_sim); day_offset=1)
prob_25 = PBDMProblem(pop_25, weather_25, (1, n_sim))
sol_25 = solve(prob_25, DirectIteration())

nymph_traj = stage_trajectory(sol_25, 1)
adult_traj = stage_trajectory(sol_25, 2)
total_traj = total_population(sol_25)

println("\nConstant 25°C simulation ($(n_sim) days):")
println("  Return code: $(sol_25.retcode)")
println("  Nymph peak: $(round(maximum(nymph_traj), digits=1)) at day $(argmax(nymph_traj))")
println("  Adult peak: $(round(maximum(adult_traj), digits=1)) at day $(argmax(adult_traj))")

fig3 = Figure(size=(900, 550))
ax3 = Axis(fig3[1, 1],
    title="T. maculata Cohort at Constant 25°C ($(n_sim) days)",
    xlabel="Day",
    ylabel="Population",
    xlabelsize=14, ylabelsize=14)

lines!(ax3, sol_25.t, nymph_traj, linewidth=2.5, color=:steelblue, label="Nymphs")
lines!(ax3, sol_25.t, adult_traj, linewidth=2.5, color=:firebrick, label="Adults")
lines!(ax3, sol_25.t, total_traj, linewidth=1.5, color=:gray40, linestyle=:dash, label="Total")

# Mark adult peak
adult_peak_day = sol_25.t[argmax(adult_traj)]
adult_peak_val = maximum(adult_traj)
scatter!(ax3, [adult_peak_day], [adult_peak_val], color=:firebrick, markersize=10)
text!(ax3, adult_peak_day + 3, adult_peak_val,
      text="peak day $(adult_peak_day)",
      align=(:left, :center), fontsize=10, color=:firebrick)

xlims!(ax3, 1, n_sim)
ylims!(ax3, 0, nothing)
axislegend(ax3, position=:rt)

save(joinpath(figdir, "sim_constant_25C.png"), fig3, px_per_unit=2)
println("Saved sim_constant_25C.png")

# ============================================================
# Figure 4: Factorial biocontrol scenarios — bar chart
# ============================================================

# Helper: build population and run simulation, return peak adult count
function run_scenario(; W0_nymph=200.0, μ_n=μ_NYMPH, μ_a=μ_ADULT, n_days=150)
    stages = [
        LifeStage(:nymph, DistributedDelay(6, DD_NYMPH_COHORT; W0=W0_nymph), nymph_dev, μ_n),
        LifeStage(:adult, DistributedDelay(15, DD_ADULT_COHORT; W0=0.0),     adult_dev, μ_a),
    ]
    pop = Population(:alfalfa_aphid, stages)
    weather = WeatherSeries(fill(25.0, n_days); day_offset=1)
    prob = PBDMProblem(pop, weather, (1, n_days))
    sol = solve(prob, DirectIteration())
    peak_adult = maximum(stage_trajectory(sol, 2))
    return peak_adult
end

# Six scenarios
scenario_names = [
    "No control",
    "Parasitoid\nonly",
    "Predator\nonly",
    "Resistant\nvariety",
    "Parasitoid\n+ Predator",
    "All three\ncombined",
]

# (1) No control: baseline mortality, 200 initial nymphs
peak_1 = run_scenario(W0_nymph=200.0, μ_n=μ_NYMPH, μ_a=μ_ADULT)
# (2) Parasitoid only: increased aphid mortality
peak_2 = run_scenario(W0_nymph=200.0, μ_n=0.08, μ_a=0.08)
# (3) Predator only: moderate mortality increase
peak_3 = run_scenario(W0_nymph=200.0, μ_n=0.06, μ_a=0.06)
# (4) Resistant variety: reduced fecundity → fewer initial aphids
peak_4 = run_scenario(W0_nymph=50.0, μ_n=μ_NYMPH, μ_a=μ_ADULT)
# (5) Parasitoid + Predator: high combined mortality
peak_5 = run_scenario(W0_nymph=200.0, μ_n=0.10, μ_a=0.10)
# (6) All three: high mortality + resistant variety
peak_6 = run_scenario(W0_nymph=50.0, μ_n=0.10, μ_a=0.10)

peaks = [peak_1, peak_2, peak_3, peak_4, peak_5, peak_6]

println("\nBiocontrol factorial scenarios (peak adult population):")
for (i, (name, pk)) in enumerate(zip(scenario_names, peaks))
    label = replace(name, "\n" => " ")
    println("  ($i) $(label): $(round(pk, digits=1))")
end

bar_colors = [:gray60, :orange, :teal, :mediumpurple, :goldenrod, :forestgreen]

fig4 = Figure(size=(900, 550))
ax4 = Axis(fig4[1, 1],
    title="T. maculata — Multi-Factor Biocontrol Suppression at 25°C",
    xlabel="Scenario",
    ylabel="Peak adult population",
    xticks=(1:6, scenario_names),
    xticklabelsize=10,
    xlabelsize=14, ylabelsize=14)

barplot!(ax4, 1:6, peaks, color=bar_colors, strokewidth=1, strokecolor=:black)

# Annotate bars with values
for (i, pk) in enumerate(peaks)
    text!(ax4, Float64(i), pk + maximum(peaks) * 0.02,
          text=string(round(Int, pk)),
          align=(:center, :bottom), fontsize=10)
end

ylims!(ax4, 0, maximum(peaks) * 1.15)

save(joinpath(figdir, "factorial_biocontrol.png"), fig4, px_per_unit=2)
println("Saved factorial_biocontrol.png")

# ============================================================
# Figure 5: Seasonal California Central Valley simulation
# ============================================================

n_seasonal = 365

# Generate Central Valley CA temperatures: mean 17°C, amplitude 10°C
# Peak heat in mid-July (day ≈ 200), so sine phase offset = 200 − 365/4 ≈ 109
temps_ca = [17.0 + 10.0 * sin(2π * (d - 109) / 365) for d in 1:n_seasonal]
weather_ca = WeatherSeries(temps_ca; day_offset=1)

# Reproduction function: adults produce nymphs scaled by thermal suitability.
# A Gaussian heat penalty above 24°C captures reproductive suppression at
# high summer temperatures — the defining feature of T. maculata seasonality.
const FECUNDITY_SCALE = 1.5
const CARRYING_CAPACITY = 10_000.0
function reproduce_aphid(pop, w, p, day)
    T = w.T_mean
    adult_total = delay_total(pop.stages[2].delay)
    total_pop = delay_total(pop.stages[1].delay) + adult_total
    thermal_rate = development_rate(nymph_dev, T)
    # Wider Gaussian (σ≈2.2) centred at 25°C — suppresses but doesn't eliminate
    # reproduction during peak summer, leaving a residual population that
    # rebounds when autumn temperatures drop below 24°C.
    heat_penalty = T > 25.0 ? exp(-0.2 * (T - 25.0)^2) : 1.0
    # Logistic density dependence: finite alfalfa crop limits aphid population
    density_effect = max(0.0, 1.0 - total_pop / CARRYING_CAPACITY)
    return max(0.0, adult_total * thermal_rate * FECUNDITY_SCALE * heat_penalty * density_effect)
end

# Higher baseline mortality for seasonal simulation — faster turnover produces
# visible generation structure within one year.
const μ_NYMPH_SEASONAL = 0.08
const μ_ADULT_SEASONAL = 0.06

stages_ca = [
    LifeStage(:nymph, DistributedDelay(6, DD_NYMPH_SEASONAL; W0=100.0), nymph_dev, μ_NYMPH_SEASONAL),
    LifeStage(:adult, DistributedDelay(15, DD_ADULT_SEASONAL; W0=0.0),  adult_dev, μ_ADULT_SEASONAL),
]
pop_ca = Population(:alfalfa_aphid, stages_ca)

prob_ca = PBDMProblem(DensityDependent(), pop_ca, weather_ca, (1, n_seasonal))
sol_ca = solve(prob_ca, DirectIteration(); reproduction_fn=reproduce_aphid)

nymph_ca = stage_trajectory(sol_ca, 1)
adult_ca = stage_trajectory(sol_ca, 2)
total_ca = total_population(sol_ca)

println("\nSeasonal California simulation ($(n_seasonal) days):")
println("  Return code: $(sol_ca.retcode)")
println("  Nymph peak: $(round(maximum(nymph_ca), digits=1)) at day $(argmax(nymph_ca))")
println("  Adult peak: $(round(maximum(adult_ca), digits=1)) at day $(argmax(adult_ca))")
println("  Total peak: $(round(maximum(total_ca), digits=1)) at day $(argmax(total_ca))")

# Temperature for secondary axis
temps_plot = temps_ca

fig5 = Figure(size=(1000, 600))
ax5 = Axis(fig5[1, 1],
    title="T. maculata — Central Valley CA Seasonal Dynamics\n(mean 17°C, amplitude 10°C)",
    xlabel="Day of year",
    ylabel="Population",
    xlabelsize=14, ylabelsize=14)

lines!(ax5, sol_ca.t, nymph_ca, linewidth=2, color=:steelblue, label="Nymphs")
lines!(ax5, sol_ca.t, adult_ca, linewidth=2, color=:firebrick, label="Adults")
lines!(ax5, sol_ca.t, total_ca, linewidth=1.5, color=:gray40, linestyle=:dash, label="Total")

# Shade summer heat period (approximately June–August, days 152–243)
vspan!(ax5, 152, 243, color=(:red, 0.06))
text!(ax5, 197.0, maximum(total_ca) * 0.95,
      text="summer heat\nsuppression", align=(:center, :top),
      fontsize=9, color=:red)

xlims!(ax5, 1, n_seasonal)
ylims!(ax5, 0, nothing)
axislegend(ax5, position=:rt)

# Temperature overlay on secondary axis
ax5b = Axis(fig5[1, 1],
    ylabel="Temperature (°C)",
    ylabelsize=12,
    yaxisposition=:right,
    yticklabelcolor=:darkorange,
    ylabelcolor=:darkorange)
hidexdecorations!(ax5b)
hidespines!(ax5b)
lines!(ax5b, 1:n_seasonal, temps_plot, linewidth=1.2, color=(:darkorange, 0.6),
       linestyle=:dot)
ylims!(ax5b, 0, 35)

save(joinpath(figdir, "seasonal_california.png"), fig5, px_per_unit=2)
println("Saved seasonal_california.png")

println("\n" * "=" ^ 60)
println("All alfalfa aphid validation figures saved to:\n  $(figdir)")
println("=" ^ 60)
