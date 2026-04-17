# Validation script: Vine mealybug (Planococcus ficus)
# Generates development rate curves and population dynamics for comparison
# with Gutierrez et al. (2008), Daane et al. (2006)

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using PhysiologicallyBasedDemographicModels
using CairoMakie

outdir = joinpath(@__DIR__, "figures", "vine_mealybug")
mkpath(outdir)

# --- Vine mealybug development rates ---
# Literature: lower threshold 14-16°C, upper 35°C, optimal ~28°C
# DD egg-adult: ~400-500 above 14°C
egg_dev = BriereDevelopmentRate(2.5e-4, 14.0, 35.0)
nymph_dev = BriereDevelopmentRate(8.0e-5, 14.0, 35.0)
adult_dev = LinearDevelopmentRate(14.0, 35.0)

# Parasitoid Anagyrus pseudococci
para_dev = BriereDevelopmentRate(1.2e-4, 11.0, 35.0)

# --- Figure 1: Development rate curves ---
temps = range(5.0, 40.0, length=200)
egg_r = [development_rate(egg_dev, T) for T in temps]
nymph_r = [development_rate(nymph_dev, T) for T in temps]
para_r = [development_rate(para_dev, T) for T in temps]

fig1 = Figure(size=(700, 400))
ax1 = Axis(fig1[1, 1],
    xlabel="Temperature (°C)", ylabel="Development rate (1/day)",
    title="P. ficus and A. pseudococci Development Rates")
lines!(ax1, collect(temps), egg_r, label="Mealybug egg", linewidth=2, color=:red)
lines!(ax1, collect(temps), nymph_r, label="Mealybug nymph", linewidth=2, color=:blue)
lines!(ax1, collect(temps), para_r, label="Parasitoid", linewidth=2, color=:green, linestyle=:dash)
# Literature data (Daane et al. 2006)
scatter!(ax1, [20.0, 25.0, 30.0], [0.06, 0.12, 0.16],
    label="Egg (Daane)", marker=:circle, color=:red, markersize=10)
scatter!(ax1, [20.0, 25.0, 30.0], [0.02, 0.045, 0.06],
    label="Nymph (Daane)", marker=:diamond, color=:blue, markersize=10)
axislegend(ax1, position=:lt)
save(joinpath(outdir, "devrate_curves.png"), fig1, px_per_unit=2)
println("Saved devrate_curves.png")

# --- Figure 2: Mortality curves ---
function mb_mortality(T, T_opt, mu_min, a_mort)
    return mu_min + a_mort * (T - T_opt)^2
end

egg_mort = [mb_mortality(T, 27.0, 0.01, 0.001) for T in temps]
nymph_mort = [mb_mortality(T, 26.0, 0.015, 0.0012) for T in temps]

fig2 = Figure(size=(700, 400))
ax2 = Axis(fig2[1, 1],
    xlabel="Temperature (°C)", ylabel="Daily mortality rate",
    title="P. ficus Temperature-Dependent Mortality")
lines!(ax2, collect(temps), egg_mort, label="Egg", linewidth=2, color=:red)
lines!(ax2, collect(temps), nymph_mort, label="Nymph", linewidth=2, color=:blue)
axislegend(ax2, position=:ct)
save(joinpath(outdir, "mortality_curves.png"), fig2, px_per_unit=2)
println("Saved mortality_curves.png")

# --- Figure 3: Cohort at 25°C ---
e = LifeStage(:egg, DistributedDelay(8, 80.0; W0=0.0), egg_dev, 0.01)
n = LifeStage(:nymph, DistributedDelay(15, 350.0; W0=0.0), nymph_dev, 0.015)
a = LifeStage(:adult, DistributedDelay(5, 400.0; W0=0.0), adult_dev, 0.025)
e.delay.W[1] = 400.0

pop = Population(:p_ficus, [e, n, a])
weather = WeatherSeries([DailyWeather(25.0) for _ in 1:200])
prob = PBDMProblem(pop, weather, (1, 200))
sol = solve(prob, DirectIteration())

nd = size(sol.stage_totals, 2) - 1
fig3 = Figure(size=(700, 400))
ax3 = Axis(fig3[1, 1], xlabel="Day", ylabel="Population",
    title="P. ficus Cohort at 25°C")
lines!(ax3, 1:nd, sol.stage_totals[1, 2:end], label="Egg", linewidth=2, color=:red)
lines!(ax3, 1:nd, sol.stage_totals[2, 2:end], label="Nymph", linewidth=2, color=:blue)
lines!(ax3, 1:nd, sol.stage_totals[3, 2:end], label="Adult", linewidth=2, color=:green)
axislegend(ax3, position=:rt)
save(joinpath(outdir, "sim_constant_25C.png"), fig3, px_per_unit=2)
println("Saved sim_constant_25C.png")

# --- Figure 4: Refuge effect (bark vs leaves) ---
# Bark refuges protect mealybugs from parasitoids
# Model: no refuge (high parasitoid mortality=0.06) vs with refuge (low mortality=0.015)

scenarios = [
    ("No biocontrol", 0.015, :gray),
    ("Biocontrol, no refuge", 0.06, :blue),
    ("Biocontrol + bark refuge", 0.03, :orange),
]

fig4 = Figure(size=(700, 400))
ax4 = Axis(fig4[1, 1], xlabel="Day", ylabel="Total mealybug",
    title="Refuge Effect on Biocontrol Efficacy")

napa_temps = [15.0 + 8.0 * sin(2pi * (d - 100) / 365) for d in 1:365]

for (label, mort, col) in scenarios
    e4 = LifeStage(:egg, DistributedDelay(8, 80.0; W0=0.0), egg_dev, 0.01)
    n4 = LifeStage(:nymph, DistributedDelay(15, 350.0; W0=0.0), nymph_dev, mort)
    a4 = LifeStage(:adult, DistributedDelay(5, 400.0; W0=0.0), adult_dev, 0.025)
    e4.delay.W[1] = 100.0
    pop4 = Population(:p_ficus, [e4, n4, a4])
    w4 = WeatherSeries([DailyWeather(T) for T in napa_temps])
    prob4 = PBDMProblem(pop4, w4, (1, 365))
    sol4 = solve(prob4, DirectIteration())
    nd4 = size(sol4.stage_totals, 2) - 1
    total4 = vec(sum(sol4.stage_totals[:, 2:end], dims=1))
    lines!(ax4, 1:nd4, total4, label=label, linewidth=2, color=col)
end
axislegend(ax4, position=:rt)
save(joinpath(outdir, "refuge_effect.png"), fig4, px_per_unit=2)
println("Saved refuge_effect.png")

# --- Figure 5: Wine region comparison ---
regions = [
    ("Napa Valley", 15.0, 8.0, :blue),
    ("Central Valley", 18.0, 10.0, :orange),
    ("Southern Italy", 18.5, 9.0, :red),
    ("Western Cape SA", 17.0, 6.0, :green),
]

fig5 = Figure(size=(700, 400))
ax5 = Axis(fig5[1, 1], xlabel="Day", ylabel="Total mealybug",
    title="Vine Mealybug Across Wine Regions")

for (label, tmean, tamp, col) in regions
    tvec = [tmean + tamp * sin(2pi * (d - 100) / 365) for d in 1:365]
    e5 = LifeStage(:egg, DistributedDelay(8, 80.0; W0=0.0), egg_dev, 0.01)
    n5 = LifeStage(:nymph, DistributedDelay(15, 350.0; W0=0.0), nymph_dev, 0.015)
    a5 = LifeStage(:adult, DistributedDelay(5, 400.0; W0=0.0), adult_dev, 0.025)
    e5.delay.W[1] = 100.0
    pop5 = Population(:p_ficus, [e5, n5, a5])
    w5 = WeatherSeries([DailyWeather(T) for T in tvec])
    prob5 = PBDMProblem(pop5, w5, (1, 365))
    sol5 = solve(prob5, DirectIteration())
    nd5 = size(sol5.stage_totals, 2) - 1
    total5 = vec(sum(sol5.stage_totals[:, 2:end], dims=1))
    lines!(ax5, 1:nd5, total5, label=label, linewidth=2, color=col)
end
axislegend(ax5, position=:rt)
save(joinpath(outdir, "wine_regions.png"), fig5, px_per_unit=2)
println("Saved wine_regions.png")

println("\nAll 5 figures generated successfully!")
