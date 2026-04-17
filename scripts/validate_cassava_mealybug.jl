#!/usr/bin/env julia
# Validation script for Cassava Mealybug (Phenacoccus manihoti) PBDM vignette.
#
# Biology and parameters based on:
#   - Gutierrez et al. (1988) Analysis of the cassava system in Africa.
#   - Neuenschwander et al. (2003) Biological control of the cassava mealybug.
#
# The cassava mealybug is the key pest of cassava in sub-Saharan Africa.
# Its parasitoid Epidinocarsis lopezi was introduced for classical biological
# control — one of the most successful biocontrol programs in history.

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using PhysiologicallyBasedDemographicModels
using CairoMakie

figdir = joinpath(@__DIR__, "figures", "cassava_mealybug")
mkpath(figdir)

# ============================================================================
# Parameters — Phenacoccus manihoti (cassava mealybug)
# ============================================================================
# Thermal thresholds
const CM_T_LOWER = 15.0   # lower developmental threshold (°C)
const CM_T_UPPER = 35.0   # upper developmental threshold (°C)

# Degree-day requirements above T_lower
const CM_τ_EGG   = 100.0  # egg stage DD
const CM_τ_NYMPH = 250.0  # crawler + nymph DD
const CM_τ_ADULT = 300.0  # adult longevity DD

# Background mortality rates (per degree-day)
const CM_μ_EGG   = 0.01
const CM_μ_NYMPH = 0.015
const CM_μ_ADULT = 0.03

# Brière development rate parameter 'a' — calibrated so that peak rate
# at ~27-30°C gives biologically reasonable values
const CM_a_EGG   = 2.0e-5
const CM_a_NYMPH = 8.0e-6
const CM_a_ADULT = 6.5e-6

# ============================================================================
# Parameters — Epidinocarsis lopezi (parasitoid)
# ============================================================================
const EL_T_LOWER = 12.0
const EL_T_UPPER = 35.0
const EL_τ       = 200.0  # egg-to-adult DD above 12°C
const EL_a       = 2.5e-5

# ============================================================================
# Temperature-dependent mortality (U-shaped)
# ============================================================================
# μ(T) = α(T - T_opt)² + μ_min, clipped to [μ_min, 1.0]
# Sharp increase above 33°C
function temp_mortality(T, μ_min, T_opt, α)
    return clamp(α * (T - T_opt)^2 + μ_min, μ_min, 1.0)
end

cm_mort_egg(T)   = temp_mortality(T, 0.01, 27.0, 0.0005)
cm_mort_nymph(T) = temp_mortality(T, 0.015, 27.0, 0.0006)
cm_mort_adult(T) = temp_mortality(T, 0.03, 28.0, 0.0008)

# ============================================================================
# Development rate models
# ============================================================================
egg_dev   = BriereDevelopmentRate(CM_a_EGG,   CM_T_LOWER, CM_T_UPPER)
nymph_dev = BriereDevelopmentRate(CM_a_NYMPH, CM_T_LOWER, CM_T_UPPER)
adult_dev = BriereDevelopmentRate(CM_a_ADULT, CM_T_LOWER, CM_T_UPPER)
para_dev  = BriereDevelopmentRate(EL_a,       EL_T_LOWER, EL_T_UPPER)

# Temperature range for curves
Ts = range(5.0, 40.0, length=300)

# ============================================================================
# Figure 1: Development rate curves
# ============================================================================
egg_rates   = [development_rate(egg_dev, T) for T in Ts]
nymph_rates = [development_rate(nymph_dev, T) for T in Ts]
adult_rates = [development_rate(adult_dev, T) for T in Ts]
para_rates  = [development_rate(para_dev, T) for T in Ts]

fig1 = Figure(size=(900, 600))
ax1 = Axis(fig1[1, 1],
    xlabel="Temperature (°C)", ylabel="Development rate (1/day)",
    title="Cassava Mealybug & Parasitoid — Development Rates")
lines!(ax1, collect(Ts), egg_rates,   linewidth=2, color=:firebrick,
       linestyle=:solid, label="CM egg")
lines!(ax1, collect(Ts), nymph_rates, linewidth=2, color=:darkorange,
       linestyle=:dash, label="CM nymph")
lines!(ax1, collect(Ts), adult_rates, linewidth=2, color=:goldenrod,
       linestyle=:dot, label="CM adult")
lines!(ax1, collect(Ts), para_rates,  linewidth=2, color=:steelblue,
       linestyle=:dashdot, label="E. lopezi")
axislegend(ax1, position=:lt)
xlims!(ax1, 5, 40)
save(joinpath(figdir, "devrate_curves.png"), fig1, px_per_unit=2)
println("Saved devrate_curves.png — CM egg peak: ",
        round(maximum(egg_rates), digits=4),
        ", parasitoid peak: ", round(maximum(para_rates), digits=4))

# ============================================================================
# Figure 2: Mortality curves (U-shaped)
# ============================================================================
mort_egg_vals   = [cm_mort_egg(T) for T in Ts]
mort_nymph_vals = [cm_mort_nymph(T) for T in Ts]
mort_adult_vals = [cm_mort_adult(T) for T in Ts]

fig2 = Figure(size=(900, 600))
ax2 = Axis(fig2[1, 1],
    xlabel="Temperature (°C)", ylabel="Mortality rate (per DD)",
    title="Cassava Mealybug — Temperature-Dependent Mortality")
lines!(ax2, collect(Ts), mort_egg_vals,   linewidth=2, color=:firebrick,
       linestyle=:solid, label="Egg")
lines!(ax2, collect(Ts), mort_nymph_vals, linewidth=2, color=:darkorange,
       linestyle=:dash, label="Nymph")
lines!(ax2, collect(Ts), mort_adult_vals, linewidth=2, color=:goldenrod,
       linestyle=:dot, label="Adult")
vlines!(ax2, [33.0], color=:gray, linestyle=:dash, linewidth=1)
text!(ax2, 33.5, 0.15, text="33°C", fontsize=12, color=:gray)
axislegend(ax2, position=:lt)
xlims!(ax2, 5, 40)
ylims!(ax2, 0, 0.25)
save(joinpath(figdir, "mortality_curves.png"), fig2, px_per_unit=2)
println("Saved mortality_curves.png — egg min mort: ",
        round(minimum(mort_egg_vals), digits=4),
        ", adult mort at 35°C: ", round(cm_mort_adult(35.0), digits=4))

# ============================================================================
# Figure 3: Simulation at constant 27°C (optimal temperature)
# ============================================================================
function build_mealybug(; n_eggs=0.0, μ_egg=CM_μ_EGG, μ_nymph=CM_μ_NYMPH,
                          μ_adult=CM_μ_ADULT)
    egg_stage = LifeStage(:egg,
        DistributedDelay(10, CM_τ_EGG; W0=n_eggs / 10),
        egg_dev, μ_egg)
    nymph_stage = LifeStage(:nymph,
        DistributedDelay(15, CM_τ_NYMPH; W0=0.0),
        nymph_dev, μ_nymph)
    adult_stage = LifeStage(:adult,
        DistributedDelay(10, CM_τ_ADULT; W0=0.0),
        adult_dev, μ_adult)
    return Population(:mealybug, [egg_stage, nymph_stage, adult_stage])
end

# Constant 27°C for 200 days — tspan (1, N+1) yields N simulation steps
n_days_const = 200
weather_27 = WeatherSeries([DailyWeather(27.0) for _ in 1:(n_days_const + 1)])
pop_27 = build_mealybug(n_eggs=300.0)
prob_27 = PBDMProblem(pop_27, weather_27, (1, n_days_const + 1))
sol_27 = solve(prob_27, DirectIteration())

days_27 = 1:n_days_const
eggs_27   = sol_27.stage_totals[1, 2:end]
nymphs_27 = sol_27.stage_totals[2, 2:end]
adults_27 = sol_27.stage_totals[3, 2:end]

fig3 = Figure(size=(900, 600))
ax3 = Axis(fig3[1, 1],
    xlabel="Day", ylabel="Population",
    title="Mealybug Cohort at Constant 27°C (300 initial eggs)")
lines!(ax3, collect(days_27), eggs_27,   linewidth=2, color=:firebrick,  label="Eggs")
lines!(ax3, collect(days_27), nymphs_27, linewidth=2, color=:darkorange, label="Nymphs")
lines!(ax3, collect(days_27), adults_27, linewidth=2, color=:goldenrod,  label="Adults")
axislegend(ax3, position=:rt)
xlims!(ax3, 1, n_days_const)
save(joinpath(figdir, "sim_constant_27C.png"), fig3, px_per_unit=2)
println("Saved sim_constant_27C.png — peak nymphs: ",
        round(maximum(nymphs_27), digits=1),
        ", peak adults: ", round(maximum(adults_27), digits=1))

# ============================================================================
# Figure 4: Tritrophic comparison (mealybug alone vs with parasitoid)
# ============================================================================
# Tropical climate: mean 26°C, amplitude 3°C
n_days_trit = 365
tropical_temps = [26.0 + 3.0 * sin(2π * d / 365) for d in 1:(n_days_trit + 1)]
weather_trop = WeatherSeries([DailyWeather(T) for T in tropical_temps])

# Mealybug alone
pop_alone = build_mealybug(n_eggs=300.0)
prob_alone = PBDMProblem(pop_alone, weather_trop, (1, n_days_trit + 1))
sol_alone = solve(prob_alone, DirectIteration())
total_alone = vec(sum(sol_alone.stage_totals[:, 2:end], dims=1))

# Mealybug with parasitoid effect (increased adult mortality)
pop_para = build_mealybug(n_eggs=300.0, μ_adult=0.08)
prob_para = PBDMProblem(pop_para, weather_trop, (1, n_days_trit + 1))
sol_para = solve(prob_para, DirectIteration())
total_para = vec(sum(sol_para.stage_totals[:, 2:end], dims=1))

fig4 = Figure(size=(900, 600))
ax4 = Axis(fig4[1, 1],
    xlabel="Day", ylabel="Total mealybug population",
    title="Tritrophic Comparison — Mealybug ± Parasitoid (E. lopezi)")
lines!(ax4, collect(1:n_days_trit), total_alone, linewidth=2, color=:firebrick,
       label="Mealybug alone (μ_adult=0.03)")
lines!(ax4, collect(1:n_days_trit), total_para,  linewidth=2, color=:steelblue,
       linestyle=:dash, label="With parasitoid (μ_adult=0.08)")
axislegend(ax4, position=:rt)
xlims!(ax4, 1, n_days_trit)
save(joinpath(figdir, "tritrophic_comparison.png"), fig4, px_per_unit=2)
println("Saved tritrophic_comparison.png — peak alone: ",
        round(maximum(total_alone), digits=1),
        ", peak w/ parasitoid: ", round(maximum(total_para), digits=1))

# ============================================================================
# Figure 5: Climate zone comparison
# ============================================================================
function run_climate_zone(T_mean, amplitude, n_days)
    temps = [T_mean + amplitude * sin(2π * d / 365) for d in 1:(n_days + 1)]
    weather = WeatherSeries([DailyWeather(T) for T in temps])
    pop = build_mealybug(n_eggs=300.0)
    prob = PBDMProblem(pop, weather, (1, n_days + 1))
    sol = solve(prob, DirectIteration())
    return vec(sum(sol.stage_totals[:, 2:end], dims=1))
end

n_days_clim = 365
total_humid    = run_climate_zone(26.0, 3.0, n_days_clim)  # humid tropics
total_subhumid = run_climate_zone(25.0, 5.0, n_days_clim)  # sub-humid
total_semiarid = run_climate_zone(28.0, 8.0, n_days_clim)  # semi-arid

fig5 = Figure(size=(900, 600))
ax5 = Axis(fig5[1, 1],
    xlabel="Day", ylabel="Total mealybug population",
    title="Mealybug Dynamics Across African Climate Zones")
lines!(ax5, collect(1:n_days_clim), total_humid, linewidth=2, color=:forestgreen,
       label="Humid tropics (26±3°C)")
lines!(ax5, collect(1:n_days_clim), total_subhumid, linewidth=2, color=:darkorange,
       linestyle=:dash, label="Sub-humid (25±5°C)")
lines!(ax5, collect(1:n_days_clim), total_semiarid, linewidth=2, color=:firebrick,
       linestyle=:dot, label="Semi-arid (28±8°C)")
axislegend(ax5, position=:rt)
xlims!(ax5, 1, n_days_clim)
save(joinpath(figdir, "climate_zones.png"), fig5, px_per_unit=2)
println("Saved climate_zones.png — humid peak: ",
        round(maximum(total_humid), digits=1),
        ", sub-humid peak: ", round(maximum(total_subhumid), digits=1),
        ", semi-arid peak: ", round(maximum(total_semiarid), digits=1))

println("\nAll cassava mealybug validation figures saved to: ", figdir)
