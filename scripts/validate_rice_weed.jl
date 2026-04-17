#!/usr/bin/env julia
# Validation script for Rice–Weed Competition (Vignette 15).
#
# System: Oryza sativa (rice) vs Echinochloa crus-galli (barnyard grass).
# Literature: Gutierrez et al. (1984), Graf et al. (1990).
#
# Rice:          T_base 10°C, T_opt ~30-33°C, ~2000 DD to maturity
# Echinochloa:   T_base 8°C, T_opt ~30-35°C (slightly wider thermal niche)

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using PhysiologicallyBasedDemographicModels
using CairoMakie

figdir = joinpath(@__DIR__, "figures", "rice_weed")
mkpath(figdir)

# ============================================================================
# Parameters
# ============================================================================

# Brière development‐rate coefficients (calibrated for reasonable peak rates)
const RICE_a       = 2.0e-5
const RICE_T_LOWER = 10.0
const RICE_T_UPPER = 40.0

const WEED_a       = 1.8e-5
const WEED_T_LOWER = 8.0
const WEED_T_UPPER = 42.0

# Degree‐day requirements above base temperature
const RICE_τ_VEG    = 800.0   # vegetative stage DD
const RICE_τ_REPRO  = 600.0   # reproductive stage DD
const RICE_τ_GRAIN  = 600.0   # grain‐fill stage DD

const WEED_τ_VEG    = 600.0   # vegetative (faster early growth)
const WEED_τ_REPRO  = 500.0   # reproductive
const WEED_τ_SEED   = 400.0   # seed‐set

# Background mortality rates (per degree-day)
const RICE_μ = 0.01
const WEED_μ = 0.01

# Development rate models
rice_dev = BriereDevelopmentRate(RICE_a, RICE_T_LOWER, RICE_T_UPPER)
weed_dev = BriereDevelopmentRate(WEED_a, WEED_T_LOWER, WEED_T_UPPER)

# ============================================================================
# Helper: build populations
# ============================================================================
function build_rice(; n_init=200.0, μ=RICE_μ)
    veg   = LifeStage(:vegetative,
        DistributedDelay(10, RICE_τ_VEG;   W0=n_init / 10), rice_dev, μ)
    repro = LifeStage(:reproductive,
        DistributedDelay(10, RICE_τ_REPRO; W0=0.0), rice_dev, μ)
    grain = LifeStage(:grainfill,
        DistributedDelay(10, RICE_τ_GRAIN; W0=0.0), rice_dev, μ)
    Population(:rice, [veg, repro, grain])
end

function build_weed(; n_init=200.0, μ=WEED_μ)
    veg   = LifeStage(:vegetative,
        DistributedDelay(10, WEED_τ_VEG;   W0=n_init / 10), weed_dev, μ)
    repro = LifeStage(:reproductive,
        DistributedDelay(10, WEED_τ_REPRO; W0=0.0), weed_dev, μ)
    seed  = LifeStage(:seedset,
        DistributedDelay(10, WEED_τ_SEED;  W0=0.0), weed_dev, μ)
    Population(:weed, [veg, repro, seed])
end

# ============================================================================
# Figure 1: Development rate curves
# ============================================================================
Ts = range(0.0, 45.0, length=300)
rice_rates = [development_rate(rice_dev, T) for T in Ts]
weed_rates = [development_rate(weed_dev, T) for T in Ts]

fig1 = Figure(size=(900, 600))
ax1 = Axis(fig1[1, 1],
    xlabel="Temperature (°C)", ylabel="Development rate (1/day)",
    title="Rice vs Echinochloa — Brière Development Rate Curves")
lines!(ax1, collect(Ts), rice_rates, linewidth=2.5, color=:forestgreen,
       label="Rice (T_base=10, T_max=40)")
lines!(ax1, collect(Ts), weed_rates, linewidth=2.5, color=:firebrick,
       linestyle=:dash, label="Echinochloa (T_base=8, T_max=42)")
vlines!(ax1, [RICE_T_LOWER], color=:forestgreen, linestyle=:dot, linewidth=1)
vlines!(ax1, [WEED_T_LOWER], color=:firebrick,   linestyle=:dot, linewidth=1)
axislegend(ax1, position=:lt)
xlims!(ax1, 0, 45)
save(joinpath(figdir, "devrate_curves.png"), fig1, px_per_unit=2)
println("Saved devrate_curves.png — rice peak: ",
        round(maximum(rice_rates), digits=4),
        ", weed peak: ", round(maximum(weed_rates), digits=4))

# ============================================================================
# Figure 2: Biomass accumulation at constant 30°C
# ============================================================================
n_days_growth = 150
weather_30 = WeatherSeries([DailyWeather(30.0) for _ in 1:(n_days_growth + 1)])

pop_rice_alone = build_rice(n_init=200.0)
prob_rice_g = PBDMProblem(pop_rice_alone, weather_30, (1, n_days_growth + 1))
sol_rice_g = solve(prob_rice_g, DirectIteration())
total_rice_g = vec(sum(sol_rice_g.stage_totals[:, 2:end], dims=1))

pop_weed_alone = build_weed(n_init=200.0)
prob_weed_g = PBDMProblem(pop_weed_alone, weather_30, (1, n_days_growth + 1))
sol_weed_g = solve(prob_weed_g, DirectIteration())
total_weed_g = vec(sum(sol_weed_g.stage_totals[:, 2:end], dims=1))

fig2 = Figure(size=(900, 600))
ax2 = Axis(fig2[1, 1],
    xlabel="Day", ylabel="Total biomass (arbitrary units)",
    title="Growth Curves at Constant 30°C — Rice vs Echinochloa")
lines!(ax2, collect(1:n_days_growth), total_rice_g, linewidth=2.5,
       color=:forestgreen, label="Rice alone")
lines!(ax2, collect(1:n_days_growth), total_weed_g, linewidth=2.5,
       color=:firebrick, linestyle=:dash, label="Echinochloa alone")
axislegend(ax2, position=:rt)
xlims!(ax2, 1, n_days_growth)
save(joinpath(figdir, "growth_curves_30C.png"), fig2, px_per_unit=2)
println("Saved growth_curves_30C.png — rice final: ",
        round(total_rice_g[end], digits=1),
        ", weed final: ", round(total_weed_g[end], digits=1))

# ============================================================================
# Figure 3: Competition simulation at 30°C
# ============================================================================
# Competition modeled as increased mortality when both species are present
const COMP_μ = 0.05

# Rice alone
pop_rice_a = build_rice(n_init=200.0)
prob_rice_a = PBDMProblem(pop_rice_a, weather_30, (1, n_days_growth + 1))
sol_rice_a = solve(prob_rice_a, DirectIteration())
total_rice_a = vec(sum(sol_rice_a.stage_totals[:, 2:end], dims=1))

# Weed alone
pop_weed_a = build_weed(n_init=200.0)
prob_weed_a = PBDMProblem(pop_weed_a, weather_30, (1, n_days_growth + 1))
sol_weed_a = solve(prob_weed_a, DirectIteration())
total_weed_a = vec(sum(sol_weed_a.stage_totals[:, 2:end], dims=1))

# Competition: both present with elevated mortality
pop_rice_c = build_rice(n_init=200.0, μ=COMP_μ)
prob_rice_c = PBDMProblem(pop_rice_c, weather_30, (1, n_days_growth + 1))
sol_rice_c = solve(prob_rice_c, DirectIteration())
total_rice_c = vec(sum(sol_rice_c.stage_totals[:, 2:end], dims=1))

pop_weed_c = build_weed(n_init=200.0, μ=COMP_μ)
prob_weed_c = PBDMProblem(pop_weed_c, weather_30, (1, n_days_growth + 1))
sol_weed_c = solve(prob_weed_c, DirectIteration())
total_weed_c = vec(sum(sol_weed_c.stage_totals[:, 2:end], dims=1))

fig3 = Figure(size=(900, 600))
ax3 = Axis(fig3[1, 1],
    xlabel="Day", ylabel="Total biomass (arbitrary units)",
    title="Competition at 30°C — Rice vs Echinochloa")
lines!(ax3, collect(1:n_days_growth), total_rice_a, linewidth=2.5,
       color=:forestgreen, label="Rice alone")
lines!(ax3, collect(1:n_days_growth), total_weed_a, linewidth=2.5,
       color=:firebrick, label="Weed alone")
lines!(ax3, collect(1:n_days_growth), total_rice_c, linewidth=2,
       color=:forestgreen, linestyle=:dash, label="Rice (competing)")
lines!(ax3, collect(1:n_days_growth), total_weed_c, linewidth=2,
       color=:firebrick, linestyle=:dash, label="Weed (competing)")
axislegend(ax3, position=:rt)
xlims!(ax3, 1, n_days_growth)
save(joinpath(figdir, "competition_30C.png"), fig3, px_per_unit=2)
println("Saved competition_30C.png — rice alone peak: ",
        round(maximum(total_rice_a), digits=1),
        ", rice competing peak: ", round(maximum(total_rice_c), digits=1))

# ============================================================================
# Figure 4: Weed density vs rice yield (bar chart)
# ============================================================================
weed_densities = [0, 50, 100, 200, 400]
n_days_yield = 150
weather_28 = WeatherSeries([DailyWeather(28.0) for _ in 1:(n_days_yield + 1)])

rice_yields = Float64[]
for wd in weed_densities
    rice_μ = 0.01 + 0.0001 * wd
    pop_r = build_rice(n_init=200.0, μ=rice_μ)
    prob_r = PBDMProblem(pop_r, weather_28, (1, n_days_yield + 1))
    sol_r = solve(prob_r, DirectIteration())
    final_biomass = sum(sol_r.stage_totals[:, end])
    push!(rice_yields, final_biomass)
end

weed_free_yield = rice_yields[1]
yield_loss_pct = [100.0 * (1.0 - y / weed_free_yield) for y in rice_yields]

fig4 = Figure(size=(900, 600))
ax4a = Axis(fig4[1, 1],
    xlabel="Initial weed density", ylabel="Rice final biomass",
    title="Rice Yield vs Weed Density (150 days, 28°C)",
    xticks=(1:length(weed_densities), string.(weed_densities)))
barplot!(ax4a, 1:length(weed_densities), rice_yields,
         color=:forestgreen, strokewidth=1, strokecolor=:black)

ax4b = Axis(fig4[2, 1],
    xlabel="Initial weed density", ylabel="Yield loss (%)",
    title="Yield Loss Relative to Weed-Free",
    xticks=(1:length(weed_densities), string.(weed_densities)))
barplot!(ax4b, 1:length(weed_densities), yield_loss_pct,
         color=:firebrick, strokewidth=1, strokecolor=:black)

save(joinpath(figdir, "weed_density_yield.png"), fig4, px_per_unit=2)
println("Saved weed_density_yield.png — weed-free yield: ",
        round(weed_free_yield, digits=1),
        ", max yield loss: ", round(yield_loss_pct[end], digits=1), "%")

# ============================================================================
# Figure 5: Climate zone comparison
# ============================================================================
climates = [
    ("Tropical lowland\n(SE Asia, 28±3°C)",    28.0, 3.0),
    ("Subtropical\n(S. China, 24±6°C)",         24.0, 6.0),
    ("Temperate\n(Japan, 20±10°C)",             20.0, 10.0),
]

n_days_clim = 180

function run_climate(T_mean, amp, n_days; species=:rice)
    temps = [T_mean + amp * sin(2π * d / 365) for d in 1:(n_days + 1)]
    weather = WeatherSeries([DailyWeather(T) for T in temps])
    if species == :rice
        pop = build_rice(n_init=200.0)
    else
        pop = build_weed(n_init=200.0)
    end
    prob = PBDMProblem(pop, weather, (1, n_days + 1))
    sol = solve(prob, DirectIteration())
    return sum(sol.stage_totals[:, end])
end

rice_biomass = [run_climate(c[2], c[3], n_days_clim; species=:rice)  for c in climates]
weed_biomass = [run_climate(c[2], c[3], n_days_clim; species=:weed) for c in climates]
labels = [c[1] for c in climates]

fig5 = Figure(size=(900, 600))
ax5 = Axis(fig5[1, 1],
    xlabel="Climate zone", ylabel="Final biomass after 180 days",
    title="Rice vs Echinochloa — Climate Zone Comparison",
    xticks=(1:length(climates), labels))
dodge = 0.35
barplot!(ax5, collect(1:length(climates)) .- dodge/2, rice_biomass,
         width=dodge, color=:forestgreen, strokewidth=1, strokecolor=:black,
         label="Rice")
barplot!(ax5, collect(1:length(climates)) .+ dodge/2, weed_biomass,
         width=dodge, color=:firebrick, strokewidth=1, strokecolor=:black,
         label="Echinochloa")
axislegend(ax5, position=:rt)

save(joinpath(figdir, "climate_comparison.png"), fig5, px_per_unit=2)
println("Saved climate_comparison.png")
for (i, c) in enumerate(climates)
    name = replace(c[1], "\n" => " ")
    println("  $name — rice: ", round(rice_biomass[i], digits=1),
            ", weed: ", round(weed_biomass[i], digits=1))
end

println("\nAll rice-weed competition validation figures saved to: ", figdir)
