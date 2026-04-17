#!/usr/bin/env julia
# Validation script for BMSB (Halyomorpha halys) biocontrol PBDM
# matching key results from Gutierrez et al. (2023), BioControl 68:161–181.
#
# Generates figures in scripts/figures/bmsb/ for:
#   1. BMSB nonlinear development rate vs temperature
#   2. Comparative development rates (BMSB, Tj, Tm, As)
#   3. Single-species BMSB simulation at constant 25°C
#   4. Combined panel figure

using PhysiologicallyBasedDemographicModels
using CairoMakie

figdir = joinpath(@__DIR__, "figures", "bmsb")
mkpath(figdir)

# ============================================================
# Table 1 parameters (Gutierrez et al. 2023)
# ============================================================

# --- BMSB ---
const BMSB_T_LOWER   = 12.109  # lower developmental threshold (°C)
const BMSB_T_UPPER   = 33.5    # upper developmental threshold (°C)
const BMSB_DEV_SLOPE = 0.0018  # development rate coefficient
const BMSB_DEV_EXP   = 4.8     # exponent base in nonlinear rate function

const BMSB_DD_EGG     = 105.2   # egg stage (DD)
const BMSB_DD_NYMPH   = 436.88  # nymphal stages 1–5 (DD)
const BMSB_DD_PREOVIP = 370.88  # preoviposition period (DD)
const BMSB_DD_ADULT   = 1293.0  # reproductive adult lifespan (DD)

# --- T. japonicus (Tj) ---
const TJ_T_LOWER     = 11.164  # lower developmental threshold (°C)
const TJ_T_UPPER     = 35.0    # upper threshold in nonlinear function (°C)
const TJ_DEV_SLOPE   = 0.0061  # development rate coefficient
const TJ_DD_IMMATURE = 166.0   # total immature DD
const TJ_DD_ADULT    = 162.0   # adult lifespan DD

# --- T. mitsukurii (Tm) ---
const TM_T_LOWER     = 12.25   # lower developmental threshold (°C)
const TM_T_UPPER     = 35.0    # upper threshold in nonlinear function (°C)
const TM_DEV_SLOPE   = 0.0071  # development rate coefficient
const TM_DD_IMMATURE = 140.2   # total immature DD
const TM_DD_ADULT    = 140.5   # adult lifespan DD

# --- A. sinicus (As) ---
const AS_T_LOWER     = 12.65   # lower developmental threshold (°C)
const AS_T_UPPER     = 35.0    # upper threshold in nonlinear function (°C)
const AS_DEV_SLOPE   = 0.0095  # development rate coefficient
const AS_DD_IMMATURE = 101.9   # total immature DD
const AS_DD_ADULT    = 162.0   # adult lifespan DD

# ============================================================
# Nonlinear development rate (Eq. 3 of Gutierrez et al. 2023)
#   R(T) = slope * (T - T_lower) / (1 + base^(T - T_upper))
# ============================================================

function gutierrez_dev_rate(T, slope, T_lower, T_upper, base)
    T <= T_lower && return 0.0
    r = slope * (T - T_lower) / (1.0 + base^(T - T_upper))
    return max(0.0, r)
end

bmsb_rate(T) = gutierrez_dev_rate(T, BMSB_DEV_SLOPE, BMSB_T_LOWER, BMSB_T_UPPER, BMSB_DEV_EXP)
tj_rate(T)   = gutierrez_dev_rate(T, TJ_DEV_SLOPE, TJ_T_LOWER, TJ_T_UPPER, 4.5)
tm_rate(T)   = gutierrez_dev_rate(T, TM_DEV_SLOPE, TM_T_LOWER, TM_T_UPPER, 4.5)
as_rate(T)   = gutierrez_dev_rate(T, AS_DEV_SLOPE, AS_T_LOWER, AS_T_UPPER, 4.5)

# Print diagnostic output
Ts_diag = [15.0, 20.0, 25.0, 30.0, 33.0]
println("BMSB development rates R(T) = 0.0018*(T-12.109)/(1+4.8^(T-33.5)):")
for T in Ts_diag
    println("  T=$(T)°C: r=$(round(bmsb_rate(T), digits=5)) 1/day")
end

println("\nSpecies lower thresholds:")
println("  BMSB: $(BMSB_T_LOWER)°C")
println("  Tj:   $(TJ_T_LOWER)°C")
println("  Tm:   $(TM_T_LOWER)°C")
println("  As:   $(AS_T_LOWER)°C")

# ============================================================
# Figure 1: BMSB nonlinear development rate vs temperature
# ============================================================

Ts = range(5.0, 42.0, length=300)
bmsb_rates = [bmsb_rate(T) for T in Ts]

# Also compute Brière approximation used by LifeStage
bmsb_briere = BriereDevelopmentRate(0.0000328, BMSB_T_LOWER, BMSB_T_UPPER)
bmsb_briere_rates = [development_rate(bmsb_briere, T) for T in Ts]

# And linear approximation
bmsb_linear = LinearDevelopmentRate(BMSB_T_LOWER, BMSB_T_UPPER)
bmsb_linear_rates = [development_rate(bmsb_linear, T) / BMSB_DD_EGG for T in Ts]

fig1 = Figure(size=(800, 500))
ax1 = Axis(fig1[1, 1],
    title="BMSB Development Rate (Gutierrez et al. 2023, Eq. 3)",
    xlabel="Temperature (°C)",
    ylabel="Development rate (1/day)",
    xlabelsize=14, ylabelsize=14)

lines!(ax1, collect(Ts), bmsb_rates, color=:black, linewidth=2.5,
       label="R(T) = 0.0018(T−12.109) / (1+4.8^(T−33.5))")
lines!(ax1, collect(Ts), bmsb_briere_rates, color=:blue, linewidth=1.5,
       linestyle=:dash, label="Brière approximation")

peak_T = Ts[argmax(bmsb_rates)]
peak_r = maximum(bmsb_rates)
scatter!(ax1, [peak_T], [peak_r], color=:red, markersize=10)
text!(ax1, peak_T + 0.5, peak_r,
    text="peak $(round(peak_r, digits=4)) at $(round(peak_T, digits=1))°C",
    align=(:left, :center), fontsize=10)

xlims!(ax1, 5, 42)
ylims!(ax1, 0, nothing)
axislegend(ax1, position=:lt)

save(joinpath(figdir, "fig1_bmsb_devrate.png"), fig1, px_per_unit=2)
println("\nSaved fig1 — BMSB peak rate: $(round(peak_r, digits=5)) at $(round(peak_T, digits=1))°C")

# ============================================================
# Figure 2: All 4 species development rates on one panel
# ============================================================

tj_rates = [tj_rate(T) for T in Ts]
tm_rates = [tm_rate(T) for T in Ts]
as_rates = [as_rate(T) for T in Ts]

fig2 = Figure(size=(900, 500))
ax2 = Axis(fig2[1, 1],
    title="Temperature-Dependent Development Rates — All Species\n(Gutierrez et al. 2023, Table 1)",
    xlabel="Temperature (°C)",
    ylabel="Development rate (1/day)",
    xlabelsize=14, ylabelsize=14)

lines!(ax2, collect(Ts), bmsb_rates, linewidth=2.5, color=:red,
       label="BMSB (T_L=12.109°C)")
lines!(ax2, collect(Ts), tj_rates, linewidth=2, color=:blue,
       label="T. japonicus (T_L=11.164°C)")
lines!(ax2, collect(Ts), tm_rates, linewidth=2, color=:cyan,
       label="T. mitsukurii (T_L=12.25°C)")
lines!(ax2, collect(Ts), as_rates, linewidth=2, color=:orange,
       linestyle=:dash, label="A. sinicus (T_L=12.65°C)")

vspan!(ax2, 20.0, 30.0, color=(:green, 0.08))
text!(ax2, 25.0, maximum(as_rates) * 1.05,
    text="Optimal range", align=(:center, :bottom), fontsize=10, color=:green)
xlims!(ax2, 5, 42)
ylims!(ax2, 0, nothing)
axislegend(ax2, position=:lt)

save(joinpath(figdir, "fig2_all_devrates.png"), fig2, px_per_unit=2)
println("Saved fig2 — comparative development rates")
println("  Peak rates: BMSB=$(round(maximum(bmsb_rates), digits=5)), " *
        "Tj=$(round(maximum(tj_rates), digits=5)), " *
        "Tm=$(round(maximum(tm_rates), digits=5)), " *
        "As=$(round(maximum(as_rates), digits=5))")

# ============================================================
# Figure 3: Single-species BMSB simulation at constant 25°C
# ============================================================

# Build BMSB population with LinearDevelopmentRate (as in the vignette)
bmsb_dev_lin = LinearDevelopmentRate(BMSB_T_LOWER, BMSB_T_UPPER)

bmsb_stages = [
    LifeStage(:egg,     DistributedDelay(15, BMSB_DD_EGG;     W0=100.0), bmsb_dev_lin, 0.001),
    LifeStage(:nymph,   DistributedDelay(30, BMSB_DD_NYMPH;   W0=0.0),   bmsb_dev_lin, 0.0008),
    LifeStage(:preovip, DistributedDelay(10, BMSB_DD_PREOVIP; W0=0.0),   bmsb_dev_lin, 0.0005),
    LifeStage(:adult,   DistributedDelay(15, BMSB_DD_ADULT;   W0=0.0),   bmsb_dev_lin, 0.0005),
]

bmsb = Population(:bmsb, bmsb_stages)

println("\nBMSB population:")
println("  Life stages: $(n_stages(bmsb))")
println("  Total substages: $(n_substages(bmsb))")
println("  Initial population: $(total_population(bmsb))")

# Constant 25°C weather for 365 days
const_temp = 25.0
dd_per_day = const_temp - BMSB_T_LOWER
println("  DD/day at $(const_temp)°C: $(round(dd_per_day, digits=2))")
println("  Expected egg duration: $(round(BMSB_DD_EGG / dd_per_day, digits=1)) days")
println("  Expected nymph duration: $(round(BMSB_DD_NYMPH / dd_per_day, digits=1)) days")

n_sim_days = 365
temps_const = fill(const_temp, n_sim_days)
weather_const = WeatherSeries(temps_const; day_offset=1)

prob_bmsb = PBDMProblem(bmsb, weather_const, (1, n_sim_days))
sol_bmsb = solve(prob_bmsb, DirectIteration())

println("\nBaseline BMSB simulation (constant 25°C, 365 days):")
stage_names = [:egg, :nymph, :preovip, :adult]
for (i, sname) in enumerate(stage_names)
    traj = stage_trajectory(sol_bmsb, i)
    peak = maximum(traj)
    peak_day = argmax(traj)
    println("  $(sname): peak=$(round(peak, digits=1)) at day $(peak_day)")
end
println("  Return code: $(sol_bmsb.retcode)")
cum_dd = cumulative_degree_days(sol_bmsb)
println("  Cumulative DD at day 365: $(round(cum_dd[end], digits=1))")

fig3 = Figure(size=(900, 600))
ax3 = Axis(fig3[1, 1],
    title="BMSB Population Dynamics at Constant 25°C\n(Single-species, no biocontrol)",
    xlabel="Day",
    ylabel="Population",
    xlabelsize=14, ylabelsize=14)

colors = [:goldenrod, :forestgreen, :steelblue, :firebrick]
for (i, sname) in enumerate(stage_names)
    traj = stage_trajectory(sol_bmsb, i)
    lines!(ax3, sol_bmsb.t, traj, linewidth=2, color=colors[i],
           label=String(sname))
end

axislegend(ax3, position=:rt)
xlims!(ax3, 1, n_sim_days)

save(joinpath(figdir, "fig3_bmsb_constant25.png"), fig3, px_per_unit=2)
println("\nSaved fig3 — BMSB simulation at constant 25°C")

# ============================================================
# Figure 4: Combined panel figure
# ============================================================

fig_all = Figure(size=(1400, 1000))

# Panel (a): BMSB nonlinear development rate
ax_a = Axis(fig_all[1, 1],
    title="(a) BMSB Development Rate",
    xlabel="Temperature (°C)", ylabel="Rate (1/day)")
lines!(ax_a, collect(Ts), bmsb_rates, color=:black, linewidth=2.5)
lines!(ax_a, collect(Ts), bmsb_briere_rates, color=:blue, linewidth=1.5, linestyle=:dash)
xlims!(ax_a, 5, 42); ylims!(ax_a, 0, nothing)

# Panel (b): All species development rates
ax_b = Axis(fig_all[1, 2],
    title="(b) Comparative Development Rates",
    xlabel="Temperature (°C)", ylabel="Rate (1/day)")
lines!(ax_b, collect(Ts), bmsb_rates, linewidth=2.5, color=:red, label="BMSB")
lines!(ax_b, collect(Ts), tj_rates, linewidth=2, color=:blue, label="T. japonicus")
lines!(ax_b, collect(Ts), tm_rates, linewidth=2, color=:cyan, label="T. mitsukurii")
lines!(ax_b, collect(Ts), as_rates, linewidth=2, color=:orange, linestyle=:dash,
       label="A. sinicus")
xlims!(ax_b, 5, 42); ylims!(ax_b, 0, nothing)
axislegend(ax_b, position=:lt, labelsize=10)

# Panel (c): BMSB stage dynamics at 25°C
ax_c = Axis(fig_all[2, 1:2],
    title="(c) BMSB Stage Dynamics at Constant 25°C (365 days)",
    xlabel="Day", ylabel="Population")
for (i, sname) in enumerate(stage_names)
    traj = stage_trajectory(sol_bmsb, i)
    lines!(ax_c, sol_bmsb.t, traj, linewidth=2, color=colors[i],
           label=String(sname))
end
xlims!(ax_c, 1, n_sim_days)
axislegend(ax_c, position=:rt, labelsize=10)

save(joinpath(figdir, "fig4_combined.png"), fig_all, px_per_unit=2)
println("Saved fig4 — combined panel figure")
println("\nAll BMSB validation figures saved to: $(figdir)")
