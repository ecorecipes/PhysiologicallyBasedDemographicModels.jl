#!/usr/bin/env julia
# Validate Tuta absoluta PBDM against Ponti et al. (2021).
# Run: julia --project=. scripts/validate_tuta.jl

using PhysiologicallyBasedDemographicModels
using CairoMakie

figdir = joinpath(@__DIR__, "figures", "tuta")
mkpath(figdir)

# ============================================================
# Parameters from Ponti et al. (2021)
# ============================================================

# Eq. 3 coefficients
const A_EP  = 0.0024
const θ_L   = 7.9     # lower thermal threshold (°C)
const θ_V   = 34.95   # upper inflection temperature (°C)
const B_EP  = 3.95
const T_UPPER = 35.0   # upper lethal temperature for Brière model

# Thermal constants (degree-days above θ_L)
const DD_EGG      = 82.3
const DD_LARVA    = 218.8
const DD_PUPA     = 122.9
const DD_ADULT    = 189.2
const DD_PREOVIP  = 30.4
const DD_EP       = DD_EGG + DD_LARVA + DD_PUPA  # 424.0

println("="^65)
println("Tuta absoluta PBDM Validation — Ponti et al. (2021)")
println("="^65)
println()
println("Parameters:")
println("  Eq. 3: a=$A_EP, θ_L=$(θ_L)°C, θ_V=$(θ_V)°C, b=$B_EP")
println("  DD: egg=$DD_EGG, larva=$DD_LARVA, pupa=$DD_PUPA, adult=$DD_ADULT")
println("  DD E-P total: $DD_EP")
println("  Preoviposition: $DD_PREOVIP DD")
println()

# ============================================================
# 1. Development rate functions
# ============================================================

# Ponti et al. (2021) Eq. 3: r(T) = a*(T - θ_L) / (1 + b*(T - θ_V))
# The high-temperature suppression only activates above θ_V.
function tuta_ep_devrate_eq3(T::Real)
    T <= θ_L && return 0.0
    T >= T_UPPER && return 0.0
    denom = 1.0 + B_EP * max(0.0, T - θ_V)
    r = A_EP * (T - θ_L) / denom
    return max(0.0, r)
end

# Brière development rate calibrated to match DD_EP at 25°C.
# At 25°C the E-P rate should be (25 - 7.9)/424 ≈ 0.040 day⁻¹.
# Brière: r = a_B * T * (T - θ_L) * √(T_upper - T)
# Solving: a_B = 0.040 / (25 * 17.1 * √10) ≈ 2.96e-5
const A_BRIERE_EP = 0.040 / (25.0 * (25.0 - θ_L) * sqrt(T_UPPER - 25.0))
briere_ep = BriereDevelopmentRate(A_BRIERE_EP, θ_L, T_UPPER)

println("Development rate models:")
println("-"^70)
println("  T (°C)  |  Eq.3 (day⁻¹)  | Brière (day⁻¹) | E-P time (days)")
println("-"^70)
for T in [10.0, 15.0, 20.0, 25.0, 28.0, 30.0, 33.0, 35.0]
    r_eq3 = tuta_ep_devrate_eq3(T)
    r_br  = development_rate(briere_ep, T)
    r_use = r_br > 0 ? r_br : r_eq3
    d_str = r_use > 0 ? string(lpad(round(1.0/r_use, digits=1), 6)) : "    ∞ "
    println("  $(lpad(T, 5))   |   $(lpad(round(r_eq3, digits=5), 8))    |" *
            "  $(lpad(round(r_br, digits=5), 8))     |  $d_str")
end

# Find peak of Brière development rate
Ts_fine = range(θ_L + 0.1, T_UPPER - 0.1, length=500)
briere_rates_fine = [development_rate(briere_ep, T) for T in Ts_fine]
peak_idx = argmax(briere_rates_fine)
T_peak_br = Ts_fine[peak_idx]
r_peak_br = briere_rates_fine[peak_idx]
println("\nBrière peak: r=$(round(r_peak_br, digits=5)) day⁻¹ at $(round(T_peak_br, digits=1))°C")
println("  (E-P time at peak: $(round(1.0/r_peak_br, digits=1)) days)")
println("  a_Brière = $(round(A_BRIERE_EP, digits=7))")
println()

# ============================================================
# 2. Temperature-dependent mortality
# ============================================================

function tuta_mortality_T(T::Real)
    T_opt = 22.0
    T_cold = 5.0
    T_hot = 35.0
    base_μ = 0.005

    if T < T_cold
        return min(1.0, base_μ + 0.05 * (T_cold - T)^1.5)
    elseif T > T_hot
        return min(1.0, base_μ + 0.08 * (T - T_hot)^1.5)
    else
        cold_effect = T < T_opt ? 0.001 * (T_opt - T)^2 / (T_opt - T_cold)^2 : 0.0
        heat_effect = T > 28.0 ? 0.002 * (T - 28.0)^2 / (T_hot - 28.0)^2 : 0.0
        return base_μ + cold_effect + heat_effect
    end
end

println("Temperature-dependent daily mortality:")
println("-"^50)
println("T (°C)  |  μ(T) per day  |  Daily survival")
println("-"^50)
for T in [-5.0, 0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 28.0, 30.0, 33.0, 35.0, 38.0]
    μ = tuta_mortality_T(T)
    surv = 1.0 - μ
    println("  $(lpad(T, 5))  |    $(lpad(round(μ, digits=4), 7))    |   $(round(surv, digits=4))")
end
println()

# ============================================================
# 3. Fecundity scalar (Ponti et al. 2021, Eq. 7)
# ============================================================

function tuta_fecundity_T(T::Real)
    T <= θ_L && return 0.0
    ϕ = 0.0665 * (T - θ_L) / (1.0 + 2.2 * max(0.0, T - 27.5))
    return clamp(ϕ, 0.0, 1.0)
end

# ============================================================
# Figure 1: Egg-to-pupa development rate vs temperature
# ============================================================

Ts = range(0.0, 42.0, length=300)
rates_eq3_plot = [tuta_ep_devrate_eq3(T) for T in Ts]
rates_br_plot  = [development_rate(briere_ep, T) for T in Ts]

fig1 = Figure(size=(800, 600))
ax1 = Axis(fig1[1, 1],
           xlabel="Temperature (°C)", ylabel="Development rate (day⁻¹)",
           title="Tuta absoluta — Egg-to-Pupa Development Rate",
           xlabelsize=14, ylabelsize=14)
lines!(ax1, collect(Ts), rates_br_plot, color=:black, linewidth=2.5,
       label="Brière model")
lines!(ax1, collect(Ts), rates_eq3_plot, color=:gray, linewidth=1.5,
       linestyle=:dash, label="Eq. 3 (linear)")
vlines!(ax1, [θ_L], color=:blue, linestyle=:dot, linewidth=1,
        label="θ_L = $(θ_L)°C")
scatter!(ax1, [T_peak_br], [r_peak_br], color=:red, markersize=10,
         label="Peak $(round(T_peak_br, digits=1))°C")
xlims!(ax1, 0, 42)
ylims!(ax1, 0, max(maximum(rates_br_plot), maximum(rates_eq3_plot)) * 1.15)
axislegend(ax1, position=:lt)
save(joinpath(figdir, "fig1_ep_devrate.png"), fig1, px_per_unit=2)
println("Saved fig1_ep_devrate.png — Brière peak at $(round(T_peak_br, digits=1))°C")

# ============================================================
# Figure 2: Temperature-dependent mortality
# ============================================================

Ts_mort = range(-10.0, 42.0, length=300)
mort_vals = [tuta_mortality_T(T) for T in Ts_mort]

fig2 = Figure(size=(800, 600))
ax2 = Axis(fig2[1, 1],
           xlabel="Temperature (°C)", ylabel="Daily mortality rate",
           title="Tuta absoluta — Temperature-Dependent Mortality",
           xlabelsize=14, ylabelsize=14)
lines!(ax2, collect(Ts_mort), mort_vals, color=:black, linewidth=2.5)
hlines!(ax2, [0.005], color=:gray, linestyle=:dash, linewidth=1, label="Base μ = 0.005")
vlines!(ax2, [5.0], color=:blue, linestyle=:dot, linewidth=1, label="Cold threshold 5°C")
vlines!(ax2, [35.0], color=:red, linestyle=:dot, linewidth=1, label="Heat threshold 35°C")
xlims!(ax2, -10, 42)
ylims!(ax2, 0, 0.5)
axislegend(ax2, position=:rt)
save(joinpath(figdir, "fig2_mortality.png"), fig2, px_per_unit=2)
println("Saved fig2_mortality.png")

# ============================================================
# Figure 3: Fecundity scalar
# ============================================================

fec_vals = [tuta_fecundity_T(T) for T in Ts]

fig3 = Figure(size=(800, 600))
ax3 = Axis(fig3[1, 1],
           xlabel="Temperature (°C)", ylabel="Fecundity scalar ϕ_fec",
           title="Tuta absoluta — Temperature-Dependent Fecundity (Eq. 7)",
           xlabelsize=14, ylabelsize=14)
lines!(ax3, collect(Ts), fec_vals, color=:black, linewidth=2.5)
xlims!(ax3, 0, 42)
ylims!(ax3, 0, 1.1)
save(joinpath(figdir, "fig3_fecundity.png"), fig3, px_per_unit=2)
println("Saved fig3_fecundity.png")

# ============================================================
# PBDM simulation setup
# ============================================================

# Use LinearDevelopmentRate so DD accumulation (T - θ_L) matches the
# paper's thermal constants directly.
const lin_dev = LinearDevelopmentRate(θ_L, T_UPPER)

function make_tuta_population(; adult_W0=100.0)
    stages = [
        LifeStage(:egg,   DistributedDelay(25, DD_EGG;   W0=0.0),       lin_dev, 0.005),
        LifeStage(:larva, DistributedDelay(25, DD_LARVA;  W0=0.0),       lin_dev, 0.005),
        LifeStage(:pupa,  DistributedDelay(20, DD_PUPA;   W0=0.0),       lin_dev, 0.005),
        LifeStage(:adult, DistributedDelay(15, DD_ADULT;  W0=adult_W0),  lin_dev, 0.008),
    ]
    return Population(:tuta_absoluta, stages)
end

# ============================================================
# Figure 4: Single-species simulation at 25°C (365 days)
# ============================================================

println("\n--- Simulation at constant 25°C for 365 days ---")

pop_25 = make_tuta_population()
ws_25  = WeatherSeries([DailyWeather(25.0) for _ in 1:365]; day_offset=1)
prob_25 = PBDMProblem(pop_25, ws_25, (1, 365))
sol_25  = solve(prob_25, DirectIteration())

tp_25 = total_population(sol_25)
cdd_25 = cumulative_degree_days(sol_25)

println("  DD per day at 25°C:  $(25.0 - θ_L)")
println("  Initial population:  $(round(tp_25[1], digits=1))")
println("  Peak population:     $(round(maximum(tp_25), digits=1))")
println("  Final population:    $(round(tp_25[end], digits=1))")
println("  Annual DD (>7.9°C):  $(round(cdd_25[end], digits=0))")
println("  Est. generations:    $(round(cdd_25[end] / (DD_EP + DD_PREOVIP), digits=1))")

egg_25   = stage_trajectory(sol_25, 1)
larva_25 = stage_trajectory(sol_25, 2)
pupa_25  = stage_trajectory(sol_25, 3)
adult_25 = stage_trajectory(sol_25, 4)
days_25  = sol_25.t

fig4 = Figure(size=(900, 600))
ax4 = Axis(fig4[1, 1],
           xlabel="Day of year", ylabel="Population (individuals)",
           title="Tuta absoluta at Constant 25°C — Stage Dynamics",
           xlabelsize=14, ylabelsize=14)
lines!(ax4, days_25, egg_25,   label="Egg",   color=:gold,  linewidth=1.5)
lines!(ax4, days_25, larva_25, label="Larva", color=:green, linewidth=1.5)
lines!(ax4, days_25, pupa_25,  label="Pupa",  color=:brown, linewidth=1.5)
lines!(ax4, days_25, adult_25, label="Adult", color=:red,   linewidth=2)
axislegend(ax4, position=:rt, framevisible=false)
save(joinpath(figdir, "fig4_sim_25C.png"), fig4, px_per_unit=2)
println("Saved fig4_sim_25C.png")

# ============================================================
# Figure 5: Multi-temperature comparison (15, 20, 25, 30°C)
# ============================================================

println("\n--- Multi-temperature comparison ---")

sim_temps = [15.0, 20.0, 25.0, 30.0]
sim_colors = [:blue, :teal, :orange, :red]
sim_results = Dict{Float64, Any}()

for (T_const, col) in zip(sim_temps, sim_colors)
    pop = make_tuta_population()
    ws  = WeatherSeries([DailyWeather(T_const) for _ in 1:365]; day_offset=1)
    prob = PBDMProblem(pop, ws, (1, 365))
    sol  = solve(prob, DirectIteration())

    tp = total_population(sol)
    cdd = cumulative_degree_days(sol)
    gen = cdd[end] / (DD_EP + DD_PREOVIP)

    sim_results[T_const] = (sol=sol, tp=tp, cdd=cdd, gen=gen)

    println("  $(T_const)°C: peak=$(round(maximum(tp), digits=1)), " *
            "final=$(round(tp[end], digits=1)), " *
            "DD=$(round(cdd[end], digits=0)), " *
            "gen=$(round(gen, digits=1))")
end

fig5 = Figure(size=(1000, 700))

ax5a = Axis(fig5[1, 1],
            xlabel="Day of year", ylabel="Total population",
            title="Total Population at Different Temperatures",
            xlabelsize=13, ylabelsize=13)
for (T_const, col) in zip(sim_temps, sim_colors)
    res = sim_results[T_const]
    lines!(ax5a, res.sol.t, res.tp, color=col, linewidth=2,
           label="$(Int(T_const))°C")
end
axislegend(ax5a, position=:lt, framevisible=false)

ax5b = Axis(fig5[1, 2],
            xlabel="Day of year", ylabel="Cumulative degree-days (>$(θ_L)°C)",
            title="Degree-Day Accumulation",
            xlabelsize=13, ylabelsize=13)
for (T_const, col) in zip(sim_temps, sim_colors)
    res = sim_results[T_const]
    lines!(ax5b, 1:length(res.cdd), res.cdd, color=col, linewidth=2,
           label="$(Int(T_const))°C")
end
hlines!(ax5b, [DD_EP + DD_PREOVIP], color=:gray, linestyle=:dash, linewidth=1,
        label="1 gen ($(DD_EP + DD_PREOVIP) DD)")
axislegend(ax5b, position=:lt, framevisible=false)

# Stage breakdown at 25°C
ax5c = Axis(fig5[2, 1:2],
            xlabel="Day of year", ylabel="Population by stage",
            title="Stage Dynamics at 25°C",
            xlabelsize=13, ylabelsize=13)
sol_25_ref = sim_results[25.0].sol
lines!(ax5c, sol_25_ref.t, stage_trajectory(sol_25_ref, 1),
       label="Egg",   color=:gold,  linewidth=1.5)
lines!(ax5c, sol_25_ref.t, stage_trajectory(sol_25_ref, 2),
       label="Larva", color=:green, linewidth=1.5)
lines!(ax5c, sol_25_ref.t, stage_trajectory(sol_25_ref, 3),
       label="Pupa",  color=:brown, linewidth=1.5)
lines!(ax5c, sol_25_ref.t, stage_trajectory(sol_25_ref, 4),
       label="Adult", color=:red,   linewidth=2)
axislegend(ax5c, position=:rt, framevisible=false)

save(joinpath(figdir, "fig5_multi_temp.png"), fig5, px_per_unit=2)
println("Saved fig5_multi_temp.png")

# ============================================================
# Figure 6 (Combined): All key results in one panel
# ============================================================

println("\n--- Generating combined panel figure ---")

fig_all = Figure(size=(1400, 1000))

# Panel (a): Development rate (Brière + Eq. 3)
ax_a = Axis(fig_all[1, 1],
            xlabel="Temperature (°C)", ylabel="Rate (day⁻¹)",
            title="(a) E-P Development Rate")
lines!(ax_a, collect(Ts), rates_br_plot, color=:black, linewidth=2.5,
       label="Brière")
lines!(ax_a, collect(Ts), rates_eq3_plot, color=:gray, linewidth=1.5,
       linestyle=:dash, label="Eq. 3")
scatter!(ax_a, [T_peak_br], [r_peak_br], color=:red, markersize=8)
xlims!(ax_a, 0, 42)
ylims!(ax_a, 0, max(maximum(rates_br_plot), maximum(rates_eq3_plot)) * 1.15)
axislegend(ax_a, position=:lt, framevisible=false)

# Panel (b): Mortality
ax_b = Axis(fig_all[1, 2],
            xlabel="Temperature (°C)", ylabel="Daily mortality rate",
            title="(b) Temperature-Dependent Mortality")
lines!(ax_b, collect(Ts_mort), mort_vals, color=:black, linewidth=2.5)
hlines!(ax_b, [0.005], color=:gray, linestyle=:dash, linewidth=1)
xlims!(ax_b, -10, 42)
ylims!(ax_b, 0, 0.5)

# Panel (c): Total population across temperatures
ax_c = Axis(fig_all[2, 1],
            xlabel="Day of year", ylabel="Total population",
            title="(c) Population Dynamics")
for (T_const, col) in zip(sim_temps, sim_colors)
    res = sim_results[T_const]
    lines!(ax_c, res.sol.t, res.tp, color=col, linewidth=2,
           label="$(Int(T_const))°C")
end
axislegend(ax_c, position=:lt, framevisible=false)

# Panel (d): Degree-day accumulation
ax_d = Axis(fig_all[2, 2],
            xlabel="Day of year", ylabel="Cumulative DD (>$(θ_L)°C)",
            title="(d) Degree-Day Accumulation")
for (T_const, col) in zip(sim_temps, sim_colors)
    res = sim_results[T_const]
    lines!(ax_d, 1:length(res.cdd), res.cdd, color=col, linewidth=2,
           label="$(Int(T_const))°C")
end
hlines!(ax_d, [DD_EP + DD_PREOVIP], color=:gray, linestyle=:dash, linewidth=1)
axislegend(ax_d, position=:lt, framevisible=false)

Label(fig_all[0, :],
      "Tuta absoluta PBDM Validation — Ponti et al. (2021)",
      fontsize=18, font=:bold)

save(joinpath(figdir, "fig_combined.png"), fig_all, px_per_unit=2)
println("Saved fig_combined.png")

println("\n", "="^65)
println("All Tuta validation figures saved to: $figdir")
println("="^65)
