# Validation script: Processing tomato crop-pest IPM system
# Generates development rate curves and population dynamics for comparison
# with Gutierrez et al. (2006), Jones et al. (1999)

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using PhysiologicallyBasedDemographicModels
using CairoMakie

outdir = joinpath(@__DIR__, "figures", "tomato_ipm")
mkpath(outdir)

# --- Tomato plant development ---
# Literature: base temp 10C, optimal 25-30C, upper ~35C
# DD transplant to harvest: ~1100 above 10C (~80-100 days)
tomato_dev = BriereDevelopmentRate(1.5e-4, 10.0, 38.0)

# Key pests: tomato fruitworm (Helicoverpa zea)
# Literature: base 12C, upper 35C, DD egg-adult ~450
hw_egg_dev = BriereDevelopmentRate(3.0e-4, 12.0, 36.0)
hw_larva_dev = BriereDevelopmentRate(6.0e-5, 12.0, 36.0)
hw_pupa_dev = BriereDevelopmentRate(1.0e-4, 12.0, 36.0)
hw_adult_dev = LinearDevelopmentRate(12.0, 36.0)

# --- Figure 1: Development rate curves ---
temps = range(5.0, 42.0, length=200)
tomato_r = [development_rate(tomato_dev, T) for T in temps]
hw_egg_r = [development_rate(hw_egg_dev, T) for T in temps]
hw_larva_r = [development_rate(hw_larva_dev, T) for T in temps]
hw_pupa_r = [development_rate(hw_pupa_dev, T) for T in temps]

fig1 = Figure(size=(700, 400))
ax1 = Axis(fig1[1, 1],
    xlabel="Temperature (C)", ylabel="Development rate (1/day)",
    title="Tomato and Fruitworm Development Rates")
lines!(ax1, collect(temps), tomato_r, label="Tomato plant", linewidth=2, color=:green, linestyle=:dash)
lines!(ax1, collect(temps), hw_egg_r, label="H. zea egg", linewidth=2, color=:red)
lines!(ax1, collect(temps), hw_larva_r, label="H. zea larva", linewidth=2, color=:orange)
lines!(ax1, collect(temps), hw_pupa_r, label="H. zea pupa", linewidth=2, color=:blue)
# Literature data (Jones et al. 1999)
scatter!(ax1, [20.0, 25.0, 30.0], [0.10, 0.18, 0.22],
    label="Egg (Jones)", marker=:circle, color=:red, markersize=10)
scatter!(ax1, [20.0, 25.0, 30.0], [0.02, 0.04, 0.055],
    label="Larva (Jones)", marker=:diamond, color=:orange, markersize=10)
axislegend(ax1, position=:lt)
save(joinpath(outdir, "devrate_curves.png"), fig1, px_per_unit=2)
println("Saved devrate_curves.png")

# --- Figure 2: Mortality curves ---
function pest_mortality(T, T_opt, mu_min, a_mort)
    return mu_min + a_mort * (T - T_opt)^2
end

egg_mort = [pest_mortality(T, 27.0, 0.02, 0.001) for T in temps]
larva_mort = [pest_mortality(T, 26.0, 0.015, 0.0008) for T in temps]

fig2 = Figure(size=(700, 400))
ax2 = Axis(fig2[1, 1],
    xlabel="Temperature (C)", ylabel="Daily mortality rate",
    title="H. zea Temperature-Dependent Mortality")
lines!(ax2, collect(temps), egg_mort, label="Egg", linewidth=2, color=:red)
lines!(ax2, collect(temps), larva_mort, label="Larva", linewidth=2, color=:orange)
axislegend(ax2, position=:ct)
save(joinpath(outdir, "mortality_curves.png"), fig2, px_per_unit=2)
println("Saved mortality_curves.png")

# --- Figure 3: Fruitworm cohort at 25C ---
e = LifeStage(:egg, DistributedDelay(8, 60.0; W0=0.0), hw_egg_dev, 0.02)
l = LifeStage(:larva, DistributedDelay(15, 280.0; W0=0.0), hw_larva_dev, 0.015)
p = LifeStage(:pupa, DistributedDelay(10, 110.0; W0=0.0), hw_pupa_dev, 0.01)
a = LifeStage(:adult, DistributedDelay(5, 300.0; W0=0.0), hw_adult_dev, 0.04)
e.delay.W[1] = 500.0

pop = Population(:h_zea, [e, l, p, a])
weather = WeatherSeries([DailyWeather(25.0) for _ in 1:180])
prob = PBDMProblem(pop, weather, (1, 180))
sol = solve(prob, DirectIteration())

nd = size(sol.stage_totals, 2) - 1
fig3 = Figure(size=(700, 400))
ax3 = Axis(fig3[1, 1], xlabel="Day", ylabel="Population",
    title="H. zea Cohort at 25C")
lines!(ax3, 1:nd, sol.stage_totals[1, 2:end], label="Egg", linewidth=2, color=:red)
lines!(ax3, 1:nd, sol.stage_totals[2, 2:end], label="Larva", linewidth=2, color=:orange)
lines!(ax3, 1:nd, sol.stage_totals[3, 2:end], label="Pupa", linewidth=2, color=:blue)
lines!(ax3, 1:nd, sol.stage_totals[4, 2:end], label="Adult", linewidth=2, color=:black)
axislegend(ax3, position=:rt)
save(joinpath(outdir, "sim_constant_25C.png"), fig3, px_per_unit=2)
println("Saved sim_constant_25C.png")

# --- Figure 4: IPM comparison (no control vs biocontrol vs chemical) ---
# Central Valley CA climate
cv_temps = [17.0 + 10.0 * sin(2 * pi * (d - 100) / 365) for d in 1:200]

scenarios = [
    ("No control", 0.015, :red),
    ("Biocontrol (Trichogramma)", 0.04, :blue),
    ("Chemical (2 sprays)", 0.06, :green),
    ("IPM (bio + 1 spray)", 0.05, :purple),
]

fig4 = Figure(size=(700, 400))
ax4 = Axis(fig4[1, 1], xlabel="Day", ylabel="Total H. zea",
    title="IPM Strategy Comparison - Central Valley CA")

for (label, mort, col) in scenarios
    e4 = LifeStage(:egg, DistributedDelay(8, 60.0; W0=0.0), hw_egg_dev, 0.02)
    l4 = LifeStage(:larva, DistributedDelay(15, 280.0; W0=0.0), hw_larva_dev, mort)
    p4 = LifeStage(:pupa, DistributedDelay(10, 110.0; W0=0.0), hw_pupa_dev, 0.01)
    a4 = LifeStage(:adult, DistributedDelay(5, 300.0; W0=0.0), hw_adult_dev, 0.04)
    e4.delay.W[1] = 200.0
    pop4 = Population(:h_zea, [e4, l4, p4, a4])
    w4 = WeatherSeries([DailyWeather(T) for T in cv_temps])
    prob4 = PBDMProblem(pop4, w4, (1, 200))
    sol4 = solve(prob4, DirectIteration())
    nd4 = size(sol4.stage_totals, 2) - 1
    total4 = vec(sum(sol4.stage_totals[:, 2:end], dims=1))
    lines!(ax4, 1:nd4, total4, label=label, linewidth=2, color=col)
end
axislegend(ax4, position=:rt)
save(joinpath(outdir, "ipm_comparison.png"), fig4, px_per_unit=2)
println("Saved ipm_comparison.png")

# --- Figure 5: Yield impact bar chart ---
potential_yield = 90.0  # t/ha processing tomato

labels_s = ["No control", "Biocontrol", "Chemical", "IPM"]
morts = [0.015, 0.04, 0.06, 0.05]
colors_s = [:red, :blue, :green, :purple]
yields = Float64[]

for mort in morts
    e5 = LifeStage(:egg, DistributedDelay(8, 60.0; W0=0.0), hw_egg_dev, 0.02)
    l5 = LifeStage(:larva, DistributedDelay(15, 280.0; W0=0.0), hw_larva_dev, mort)
    p5 = LifeStage(:pupa, DistributedDelay(10, 110.0; W0=0.0), hw_pupa_dev, 0.01)
    a5 = LifeStage(:adult, DistributedDelay(5, 300.0; W0=0.0), hw_adult_dev, 0.04)
    e5.delay.W[1] = 200.0
    pop5 = Population(:h_zea, [e5, l5, p5, a5])
    w5 = WeatherSeries([DailyWeather(T) for T in cv_temps])
    prob5 = PBDMProblem(pop5, w5, (1, 200))
    sol5 = solve(prob5, DirectIteration())
    peak = maximum(vec(sum(sol5.stage_totals[:, 2:end], dims=1)))
    loss_frac = min(0.5, peak * 0.0005)
    push!(yields, potential_yield * (1.0 - loss_frac))
end

fig5 = Figure(size=(600, 400))
ax5 = Axis(fig5[1, 1],
    xlabel="Strategy", ylabel="Yield (t/ha)",
    title="Tomato Yield Under Different IPM Strategies",
    xticks=(1:4, labels_s))
barplot!(ax5, 1:4, yields, color=colors_s)
hlines!(ax5, [potential_yield], color=:gray, linestyle=:dash, label="Potential")
axislegend(ax5, position=:rb)
save(joinpath(outdir, "yield_comparison.png"), fig5, px_per_unit=2)
println("Saved yield_comparison.png")

println("\nAll 5 figures generated successfully!")
