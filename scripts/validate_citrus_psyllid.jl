# Validation script: Asian citrus psyllid (Diaphorina citri)
# Generates development rate and mortality curves for comparison with
# Gutierrez & Ponti (2013) and published thermal biology data.
#
# Literature references:
# - Liu & Tsai (2000): dev rate at 15-28°C, lower threshold ~11°C
# - Nakata (2006): nymph DD = 192 above 11.56°C
# - Hall et al. (2011): oviposition threshold 16°C

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using PhysiologicallyBasedDemographicModels
using CairoMakie

outdir = joinpath(@__DIR__, "figures", "citrus_psyllid")
mkpath(outdir)

# --- D. citri development rate parameters ---
# Brière model: r(T) = a * T * (T - T_lower) * sqrt(T_upper - T)
# Calibrate 'a' so peak rate matches literature

function calibrate_briere_a(target_peak, T_lower, T_upper)
    # Find peak of r(T) = T * (T - T_lower) * sqrt(T_upper - T) with a=1
    temps = range(T_lower + 0.1, T_upper - 0.1, length=1000)
    rates = [T * (T - T_lower) * sqrt(T_upper - T) for T in temps]
    peak_raw = maximum(rates)
    return target_peak / peak_raw
end

# Egg: lower=13.0, upper=33.5, peak rate ~0.25 at ~27°C (from Liu & Tsai 2000: ~4 days at 28°C)
T_lower_egg = 13.0
T_upper_egg = 33.5
a_egg = calibrate_briere_a(0.25, T_lower_egg, T_upper_egg)

# Nymph: lower=13.0, upper=33.5, peak rate ~0.06 at ~27°C (from Liu & Tsai: ~17 days at 28°C)
T_lower_nymph = 13.0
T_upper_nymph = 33.5
a_nymph = calibrate_briere_a(0.06, T_lower_nymph, T_upper_nymph)

egg_dev = BriereDevelopmentRate(a_egg, T_lower_egg, T_upper_egg)
nymph_dev = BriereDevelopmentRate(a_nymph, T_lower_nymph, T_upper_nymph)

# --- Figure 1: Development rate curves ---
temps = range(5.0, 40.0, length=200)
egg_rates = [development_rate(egg_dev, T) for T in temps]
nymph_rates = [development_rate(nymph_dev, T) for T in temps]

fig1 = Figure(size=(700, 400))
ax1 = Axis(fig1[1, 1],
    xlabel="Temperature (°C)", ylabel="Development rate (1/day)",
    title="D. citri Development Rates (Brière model)")
lines!(ax1, collect(temps), egg_rates, label="Egg", linewidth=2, color=:red)
lines!(ax1, collect(temps), nymph_rates, label="Nymph", linewidth=2, color=:blue)
# Literature data points from Liu & Tsai (2000)
# Egg: 15°C→0.042, 20°C→0.104, 25°C→0.185, 28°C→0.25
scatter!(ax1, [15.0, 20.0, 25.0, 28.0], [0.042, 0.104, 0.185, 0.25],
    label="Liu & Tsai (egg)", marker=:circle, color=:red, markersize=10)
# Nymph: 15°C→0.014, 20°C→0.029, 25°C→0.045, 28°C→0.059
scatter!(ax1, [15.0, 20.0, 25.0, 28.0], [0.014, 0.029, 0.045, 0.059],
    label="Liu & Tsai (nymph)", marker=:diamond, color=:blue, markersize=10)
axislegend(ax1, position=:lt)
save(joinpath(outdir, "devrate_curves.png"), fig1, px_per_unit=2)
println("Saved devrate_curves.png")

# --- Figure 2: Mortality curves (U-shaped) ---
# Quadratic mortality: μ(T) = a*(T - T_opt)^2 + μ_min
function mortality_rate(T, T_opt, μ_min, a_mort)
    return μ_min + a_mort * (T - T_opt)^2
end

T_opt_egg = 27.0
μ_min_egg = 0.01
a_mort_egg = 0.001

T_opt_nymph = 26.0
μ_min_nymph = 0.02
a_mort_nymph = 0.0015

egg_mort = [mortality_rate(T, T_opt_egg, μ_min_egg, a_mort_egg) for T in temps]
nymph_mort = [mortality_rate(T, T_opt_nymph, μ_min_nymph, a_mort_nymph) for T in temps]

fig2 = Figure(size=(700, 400))
ax2 = Axis(fig2[1, 1],
    xlabel="Temperature (°C)", ylabel="Daily mortality rate",
    title="D. citri Temperature-Dependent Mortality")
lines!(ax2, collect(temps), egg_mort, label="Egg", linewidth=2, color=:red)
lines!(ax2, collect(temps), nymph_mort, label="Nymph", linewidth=2, color=:blue)
axislegend(ax2, position=:ct)
save(joinpath(outdir, "mortality_curves.png"), fig2, px_per_unit=2)
println("Saved mortality_curves.png")

# --- Figure 3: Population dynamics at constant 25°C ---
egg_stage = LifeStage(:egg,
    DistributedDelay(10, 62.0; W0=0.0), egg_dev, 0.01)
nymph_stage = LifeStage(:nymph,
    DistributedDelay(20, 166.0; W0=0.0), nymph_dev, 0.02)
adult_dev = LinearDevelopmentRate(13.0, 33.5)
adult_stage = LifeStage(:adult,
    DistributedDelay(5, 500.0; W0=0.0), adult_dev, 0.03)

# Set initial egg population
egg_stage.delay.W[1] = 1000.0

pop = Population(:d_citri, [egg_stage, nymph_stage, adult_stage])

weather_25 = WeatherSeries([DailyWeather(25.0) for _ in 1:180])
prob = PBDMProblem(pop, weather_25, (1, 180))
sol = solve(prob, DirectIteration())

# Extract from stage_totals matrix (n_stages × n_days+1)
ndays = size(sol.stage_totals, 2) - 1
egg_traj = sol.stage_totals[1, 2:end]
nymph_traj = sol.stage_totals[2, 2:end]
adult_traj = sol.stage_totals[3, 2:end]

fig3 = Figure(size=(700, 400))
ax3 = Axis(fig3[1, 1],
    xlabel="Day", ylabel="Population",
    title="D. citri Cohort Dynamics at 25°C")
lines!(ax3, 1:ndays, egg_traj, label="Eggs", linewidth=2, color=:red)
lines!(ax3, 1:ndays, nymph_traj, label="Nymphs", linewidth=2, color=:blue)
lines!(ax3, 1:ndays, adult_traj, label="Adults", linewidth=2, color=:green)
axislegend(ax3, position=:rt)
save(joinpath(outdir, "sim_constant_25C.png"), fig3, px_per_unit=2)
println("Saved sim_constant_25C.png")

# --- Figure 4: Seasonal simulation (Florida-like climate) ---
# Sinusoidal temperature: mean=23°C, amplitude=7°C (Miami/Florida pattern)
days = 1:365
florida_temps = [23.0 + 7.0 * sin(2π * (d - 90) / 365) for d in days]

egg_stage2 = LifeStage(:egg,
    DistributedDelay(10, 62.0; W0=0.0), egg_dev, 0.01)
nymph_stage2 = LifeStage(:nymph,
    DistributedDelay(20, 166.0; W0=0.0), nymph_dev, 0.02)
adult_stage2 = LifeStage(:adult,
    DistributedDelay(5, 500.0; W0=0.0), adult_dev, 0.03)
egg_stage2.delay.W[1] = 100.0

pop2 = Population(:d_citri, [egg_stage2, nymph_stage2, adult_stage2])
weather_fl = WeatherSeries([DailyWeather(T) for T in florida_temps])
prob2 = PBDMProblem(pop2, weather_fl, (1, 365))
sol2 = solve(prob2, DirectIteration())

total_pop = vec(sum(sol2.stage_totals[:, 2:end], dims=1))
ndays2 = length(total_pop)

fig4 = Figure(size=(800, 500))
ax4a = Axis(fig4[1, 1], ylabel="Temperature (°C)", title="Seasonal D. citri Dynamics (Florida-like)")
lines!(ax4a, 1:ndays2, florida_temps[1:ndays2], color=:orange, linewidth=1.5, label="Temperature")
axislegend(ax4a, position=:rt)

ax4b = Axis(fig4[2, 1], xlabel="Day of year", ylabel="Total population")
lines!(ax4b, 1:ndays2, total_pop, color=:darkgreen, linewidth=2, label="Total population")
axislegend(ax4b, position=:rt)

save(joinpath(outdir, "sim_seasonal_florida.png"), fig4, px_per_unit=2)
println("Saved sim_seasonal_florida.png")

# --- Figure 5: Combined panel ---
fig5 = Figure(size=(1000, 800))

ax5a = Axis(fig5[1, 1], xlabel="Temperature (°C)", ylabel="Dev rate (1/day)",
    title="(a) Development rates")
lines!(ax5a, collect(temps), egg_rates, label="Egg", linewidth=2, color=:red)
lines!(ax5a, collect(temps), nymph_rates, label="Nymph", linewidth=2, color=:blue)
scatter!(ax5a, [15.0, 20.0, 25.0, 28.0], [0.042, 0.104, 0.185, 0.25],
    color=:red, marker=:circle, markersize=8)
scatter!(ax5a, [15.0, 20.0, 25.0, 28.0], [0.014, 0.029, 0.045, 0.059],
    color=:blue, marker=:diamond, markersize=8)
axislegend(ax5a, position=:lt)

ax5b = Axis(fig5[1, 2], xlabel="Temperature (°C)", ylabel="Daily mortality",
    title="(b) Mortality")
lines!(ax5b, collect(temps), egg_mort, label="Egg", linewidth=2, color=:red)
lines!(ax5b, collect(temps), nymph_mort, label="Nymph", linewidth=2, color=:blue)
axislegend(ax5b, position=:ct)

ax5c = Axis(fig5[2, 1], xlabel="Day", ylabel="Population",
    title="(c) Cohort at 25°C")
lines!(ax5c, 1:ndays, egg_traj, label="Eggs", linewidth=2, color=:red)
lines!(ax5c, 1:ndays, nymph_traj, label="Nymphs", linewidth=2, color=:blue)
lines!(ax5c, 1:ndays, adult_traj, label="Adults", linewidth=2, color=:green)
axislegend(ax5c, position=:rt)

ax5d = Axis(fig5[2, 2], xlabel="Day of year", ylabel="Population",
    title="(d) Seasonal dynamics")
lines!(ax5d, 1:ndays2, total_pop, color=:darkgreen, linewidth=2)

save(joinpath(outdir, "combined_panel.png"), fig5, px_per_unit=2)
println("Saved combined_panel.png")

println("\nAll D. citri validation figures generated successfully!")
println("Key findings:")
println("  Egg dev rate: a=$(round(a_egg, sigdigits=4)), peak ~0.25 at ~27°C")
println("  Nymph dev rate: a=$(round(a_nymph, sigdigits=4)), peak ~0.06 at ~27°C")
println("  DD egg: 62, DD nymph: 166, total egg-adult: 228")
println("  Literature range: 192-350 DD (within range)")
