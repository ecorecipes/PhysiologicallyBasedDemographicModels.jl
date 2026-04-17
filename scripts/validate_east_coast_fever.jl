# Validation script: East Coast Fever tick (Rhipicephalus appendiculatus)
# Generates development rate curves and population dynamics for comparison
# with Randolph (1993), Norval et al. (1992), Gitau et al. (2000)

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using PhysiologicallyBasedDemographicModels
using CairoMakie

outdir = joinpath(@__DIR__, "figures", "east_coast_fever")
mkpath(outdir)

# --- Tick development rates (Briere model) ---
egg_dev = BriereDevelopmentRate(5.0e-5, 12.0, 35.0)
larva_dev = BriereDevelopmentRate(8.0e-5, 12.0, 35.0)
nymph_dev = BriereDevelopmentRate(6.0e-5, 12.0, 35.0)
adult_dev = BriereDevelopmentRate(4.0e-5, 12.0, 35.0)

# --- Figure 1: Development rate curves ---
temps = range(5.0, 40.0, length=200)
egg_r = [development_rate(egg_dev, T) for T in temps]
larva_r = [development_rate(larva_dev, T) for T in temps]
nymph_r = [development_rate(nymph_dev, T) for T in temps]
adult_r = [development_rate(adult_dev, T) for T in temps]

fig1 = Figure(size=(700, 400))
ax1 = Axis(fig1[1, 1],
    xlabel="Temperature (C)", ylabel="Development rate (1/day)",
    title="R. appendiculatus Off-Host Development Rates")
lines!(ax1, collect(temps), egg_r, label="Egg", linewidth=2, color=:red)
lines!(ax1, collect(temps), larva_r, label="Larva", linewidth=2, color=:forestgreen)
lines!(ax1, collect(temps), nymph_r, label="Nymph", linewidth=2, color=:blue)
lines!(ax1, collect(temps), adult_r, label="Adult preoviposition", linewidth=2, color=:purple)
scatter!(ax1, [20.0, 25.0, 30.0], [0.008, 0.015, 0.020],
    label="Egg (Randolph)", marker=:circle, color=:red, markersize=10)
scatter!(ax1, [20.0, 25.0, 30.0], [0.012, 0.023, 0.030],
    label="Larva (Randolph)", marker=:diamond, color=:forestgreen, markersize=10)
axislegend(ax1, position=:lt)
save(joinpath(outdir, "devrate_curves.png"), fig1, px_per_unit=2)
println("Saved devrate_curves.png")

# --- Figure 2: Mortality curves ---
function tick_mortality(T, T_opt, mu_min, a_mort)
    return mu_min + a_mort * (T - T_opt)^2
end

mort_temps = collect(temps)
egg_mort = [tick_mortality(T, 25.0, 0.005, 0.0008) for T in mort_temps]
larva_mort = [tick_mortality(T, 25.0, 0.008, 0.001) for T in mort_temps]
nymph_mort = [tick_mortality(T, 24.0, 0.006, 0.0009) for T in mort_temps]

fig2 = Figure(size=(700, 400))
ax2 = Axis(fig2[1, 1],
    xlabel="Temperature (C)", ylabel="Daily mortality rate",
    title="R. appendiculatus Temperature-Dependent Mortality")
lines!(ax2, mort_temps, egg_mort, label="Egg", linewidth=2, color=:red)
lines!(ax2, mort_temps, larva_mort, label="Larva", linewidth=2, color=:forestgreen)
lines!(ax2, mort_temps, nymph_mort, label="Nymph", linewidth=2, color=:blue)
axislegend(ax2, position=:ct)
save(joinpath(outdir, "mortality_curves.png"), fig2, px_per_unit=2)
println("Saved mortality_curves.png")

# --- Figure 3: Cohort at 25C ---
e = LifeStage(:egg, DistributedDelay(8, 300.0; W0=0.0), egg_dev, 0.005)
l = LifeStage(:larva, DistributedDelay(8, 200.0; W0=0.0), larva_dev, 0.008)
n = LifeStage(:nymph, DistributedDelay(8, 250.0; W0=0.0), nymph_dev, 0.006)
a = LifeStage(:adult, DistributedDelay(5, 350.0; W0=0.0), adult_dev, 0.01)
e.delay.W[1] = 500.0

pop = Population(:r_appendiculatus, [e, l, n, a])
weather = WeatherSeries([DailyWeather(25.0) for _ in 1:365])
prob = PBDMProblem(pop, weather, (1, 365))
sol = solve(prob, DirectIteration())

nd = size(sol.stage_totals, 2) - 1
fig3 = Figure(size=(700, 400))
ax3 = Axis(fig3[1, 1], xlabel="Day", ylabel="Population",
    title="R. appendiculatus Cohort at 25C")
lines!(ax3, 1:nd, sol.stage_totals[1, 2:end], label="Egg", linewidth=2, color=:red)
lines!(ax3, 1:nd, sol.stage_totals[2, 2:end], label="Larva", linewidth=2, color=:forestgreen)
lines!(ax3, 1:nd, sol.stage_totals[3, 2:end], label="Nymph", linewidth=2, color=:blue)
lines!(ax3, 1:nd, sol.stage_totals[4, 2:end], label="Adult", linewidth=2, color=:purple)
axislegend(ax3, position=:rt)
save(joinpath(outdir, "sim_constant_25C.png"), fig3, px_per_unit=2)
println("Saved sim_constant_25C.png")

# --- Figure 4: Highland vs lowland vs mid-altitude ---
climates = [
    ("Highland (18C +/- 3)", 18.0, 3.0, :blue),
    ("Mid-altitude (22C +/- 4)", 22.0, 4.0, :orange),
    ("Lowland (28C +/- 4)", 28.0, 4.0, :red),
]

fig4 = Figure(size=(700, 400))
ax4 = Axis(fig4[1, 1], xlabel="Day", ylabel="Total tick population",
    title="Tick Population Across East African Altitude Zones")

for (label, tmean, tamp, col) in climates
    tvec = [tmean + tamp * sin(2 * pi * (d - 60) / 365) for d in 1:365]
    e4 = LifeStage(:egg, DistributedDelay(8, 300.0; W0=0.0), egg_dev, 0.005)
    l4 = LifeStage(:larva, DistributedDelay(8, 200.0; W0=0.0), larva_dev, 0.008)
    n4 = LifeStage(:nymph, DistributedDelay(8, 250.0; W0=0.0), nymph_dev, 0.006)
    a4 = LifeStage(:adult, DistributedDelay(5, 350.0; W0=0.0), adult_dev, 0.01)
    e4.delay.W[1] = 200.0
    pop4 = Population(:r_appendiculatus, [e4, l4, n4, a4])
    w4 = WeatherSeries([DailyWeather(T) for T in tvec])
    prob4 = PBDMProblem(pop4, w4, (1, 365))
    sol4 = solve(prob4, DirectIteration())
    nd4 = size(sol4.stage_totals, 2) - 1
    total4 = vec(sum(sol4.stage_totals[:, 2:end], dims=1))
    lines!(ax4, 1:nd4, total4, label=label, linewidth=2, color=col)
end
axislegend(ax4, position=:rt)
save(joinpath(outdir, "highland_lowland.png"), fig4, px_per_unit=2)
println("Saved highland_lowland.png")

# --- Figure 5: 2-year Kenya highland simulation ---
kenya_temps = [19.0 + 4.0 * (sin(2 * pi * (d - 60) / 365) + 0.3 * sin(4 * pi * (d - 60) / 365)) for d in 1:730]

e5 = LifeStage(:egg, DistributedDelay(8, 300.0; W0=0.0), egg_dev, 0.005)
l5 = LifeStage(:larva, DistributedDelay(8, 200.0; W0=0.0), larva_dev, 0.008)
n5 = LifeStage(:nymph, DistributedDelay(8, 250.0; W0=0.0), nymph_dev, 0.006)
a5 = LifeStage(:adult, DistributedDelay(5, 350.0; W0=0.0), adult_dev, 0.01)
e5.delay.W[1] = 200.0

pop5 = Population(:r_appendiculatus, [e5, l5, n5, a5])
w5 = WeatherSeries([DailyWeather(T) for T in kenya_temps])
prob5 = PBDMProblem(pop5, w5, (1, 730))
sol5 = solve(prob5, DirectIteration())

nd5 = size(sol5.stage_totals, 2) - 1
fig5 = Figure(size=(800, 500))
ax5a = Axis(fig5[1, 1], ylabel="Temperature (C)",
    title="2-Year Tick Dynamics - Kenya Highlands")
lines!(ax5a, 1:730, kenya_temps, color=:orange, linewidth=1.5)

ax5b = Axis(fig5[2, 1], xlabel="Day", ylabel="Population")
lines!(ax5b, 1:nd5, sol5.stage_totals[1, 2:end], label="Egg", color=:red)
lines!(ax5b, 1:nd5, sol5.stage_totals[2, 2:end], label="Larva", color=:forestgreen)
lines!(ax5b, 1:nd5, sol5.stage_totals[3, 2:end], label="Nymph", color=:blue)
lines!(ax5b, 1:nd5, sol5.stage_totals[4, 2:end], label="Adult", color=:purple)
axislegend(ax5b, position=:rt)
save(joinpath(outdir, "seasonal_kenya.png"), fig5, px_per_unit=2)
println("Saved seasonal_kenya.png")

println("\nAll 5 figures generated successfully!")
