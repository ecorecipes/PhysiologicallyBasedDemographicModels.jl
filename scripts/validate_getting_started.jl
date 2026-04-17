#!/usr/bin/env julia
# Validation script for the Getting Started tutorial vignette.
# Validates core PBDM mechanics: degree-day accumulation, distributed delay
# (Erlang) dynamics, development rate curves, supply/demand functional
# response, and single-cohort stage progression with conservation check.

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using PhysiologicallyBasedDemographicModels
using CairoMakie
using Statistics

figdir = joinpath(@__DIR__, "figures", "getting_started")
mkpath(figdir)

# ══════════════════════════════════════════════════════════════════════
# 1. Degree-Day Accumulation
# ══════════════════════════════════════════════════════════════════════
println("="^60)
println("1. Degree-Day Accumulation")
println("="^60)

T_base = 10.0
T_upper = 35.0
dev = LinearDevelopmentRate(T_base, T_upper)

temps = [15.0, 20.0, 25.0, 30.0]
n_days = 30
days = 1:n_days

fig1 = Figure(size=(700, 450))
ax1 = Axis(fig1[1, 1];
    xlabel="Day", ylabel="Cumulative Degree-Days",
    title="Degree-Day Accumulation at Constant Temperatures")

for T in temps
    dd_per_day = degree_days(dev, T)
    expected_dd = dd_per_day
    cum_dd = [d * dd_per_day for d in days]
    # Analytical: cumulative DD = d × max(0, T − T_base)
    cum_analytical = [d * max(0.0, T - T_base) for d in days]
    lines!(ax1, collect(days), cum_dd; linewidth=2, label="T=$(Int(T))°C ($(dd_per_day) DD/d)")
    scatter!(ax1, collect(days)[5:5:end], cum_analytical[5:5:end];
        markersize=8, marker=:xcross)

    # Verify the linear formula matches
    max_err = maximum(abs.(cum_dd .- cum_analytical))
    println("  T=$(Int(T))°C: DD/day=$(dd_per_day), 30-day total=$(cum_dd[end]), formula error=$max_err")
end

# Test edge cases
@assert degree_days(dev, 5.0) == 0.0 "Below-threshold should give 0 DD"
dd_40 = degree_days(dev, 40.0)
println("  Edge: T=5°C → $(degree_days(dev, 5.0)) DD (below threshold)")
println("  Edge: T=40°C → $(dd_40) DD (max(0, T-T_base)=$(40.0-T_base))")

axislegend(ax1; position=:lt)
save(joinpath(figdir, "01_degree_day_accumulation.png"), fig1; px_per_unit=2)
println("  → Saved 01_degree_day_accumulation.png\n")

# ══════════════════════════════════════════════════════════════════════
# 2. Erlang Distribution — Distributed Delay Transit Time PDF
# ══════════════════════════════════════════════════════════════════════
println("="^60)
println("2. Erlang Distribution — Distributed Delay")
println("="^60)

τ_mean = 200.0  # mean developmental time in DD
ks = [2, 5, 10, 25]

# Analytical Erlang PDF: f(x; k, θ) where θ = τ/k (scale parameter)
function erlang_pdf(x, k, τ)
    θ = τ / k
    # f(x) = x^(k-1) * exp(-x/θ) / (θ^k * Γ(k))
    # Use log to avoid overflow for large k
    log_f = (k - 1) * log(x) - x / θ - k * log(θ) - sum(log(i) for i in 1:(k-1))
    return exp(log_f)
end

fig2 = Figure(size=(700, 450))
ax2 = Axis(fig2[1, 1];
    xlabel="Physiological Time (DD)", ylabel="Probability Density",
    title="Erlang-Distributed Transit Time (τ=$τ_mean DD)")

x_range = range(0.1, 600.0; length=500)

for k in ks
    σ² = τ_mean^2 / k
    σ = sqrt(σ²)

    pdf_vals = [erlang_pdf(x, k, τ_mean) for x in x_range]
    lines!(ax2, collect(x_range), pdf_vals; linewidth=2, label="k=$k (σ=$(round(σ, digits=1)))")

    println("  k=$k: variance=$(round(σ², digits=1)), σ=$(round(σ, digits=1)), CV=$(round(σ/τ_mean, digits=3))")
end

vlines!(ax2, [τ_mean]; color=:black, linestyle=:dash, linewidth=1, label="mean=$τ_mean DD")
axislegend(ax2; position=:rt)
save(joinpath(figdir, "02_erlang_distribution.png"), fig2; px_per_unit=2)
println("  → Saved 02_erlang_distribution.png\n")

# ══════════════════════════════════════════════════════════════════════
# 3. Development Rate Curves
# ══════════════════════════════════════════════════════════════════════
println("="^60)
println("3. Development Rate Curves")
println("="^60)

linear_dev = LinearDevelopmentRate(10.0, 35.0)
briere_dev = BriereDevelopmentRate(2.0e-4, 10.0, 35.0)
logan_dev = LoganDevelopmentRate(0.02, 0.15, 35.0, 5.0)

T_range = range(0.0, 42.0; length=300)

r_linear = [development_rate(linear_dev, T) for T in T_range]
r_briere = [development_rate(briere_dev, T) for T in T_range]
r_logan  = [development_rate(logan_dev, T) for T in T_range]

# Normalize for visual comparison
r_linear_norm = r_linear ./ maximum(r_linear)
r_briere_norm = r_briere ./ maximum(r_briere)
r_logan_norm  = r_logan ./ maximum(r_logan)

fig3 = Figure(size=(700, 450))
ax3 = Axis(fig3[1, 1];
    xlabel="Temperature (°C)", ylabel="Relative Development Rate",
    title="Development Rate Curves (normalized to peak)")

lines!(ax3, collect(T_range), r_linear_norm; linewidth=2, label="Linear (Tₗ=10, Tᵤ=35)")
lines!(ax3, collect(T_range), r_briere_norm; linewidth=2, label="Brière (a=2e-4, Tₗ=10, Tᵤ=35)")
lines!(ax3, collect(T_range), r_logan_norm; linewidth=2, label="Logan (ψ=0.02, ρ=0.15, Tᵤ=35)")

# Mark optimum temperatures
T_opt_briere = T_range[argmax(r_briere)]
T_opt_logan = T_range[argmax(r_logan)]
vlines!(ax3, [T_opt_briere]; color=Makie.wong_colors()[2], linestyle=:dash, linewidth=1)
vlines!(ax3, [T_opt_logan]; color=Makie.wong_colors()[3], linestyle=:dash, linewidth=1)

println("  Linear: peak at upper bound, rate=$(maximum(r_linear)) DD/day")
println("  Brière: Topt=$(round(T_opt_briere, digits=1))°C, peak rate=$(round(maximum(r_briere), digits=4))")
println("  Logan:  Topt=$(round(T_opt_logan, digits=1))°C, peak rate=$(round(maximum(r_logan), digits=4))")

# Check linear model properties
@assert development_rate(linear_dev, 10.0) == 0.0 "Linear rate should be 0 at T_lower"
@assert development_rate(linear_dev, 5.0) == 0.0 "Linear rate should be 0 below T_lower"
println("  Linear checks: r(10)=0 ✓, r(5)=0 ✓, r(35)=$(development_rate(linear_dev, 35.0)) ✓")

axislegend(ax3; position=:lt)
save(joinpath(figdir, "03_development_rate_curves.png"), fig3; px_per_unit=2)
println("  → Saved 03_development_rate_curves.png\n")

# ══════════════════════════════════════════════════════════════════════
# 4. Supply-Demand Functional Response
# ══════════════════════════════════════════════════════════════════════
println("="^60)
println("4. Supply-Demand Functional Response")
println("="^60)

demand_fixed = 100.0
supply_range = range(0.0, 500.0; length=300)

fig4 = Figure(size=(700, 450))
ax4a = Axis(fig4[1, 1];
    xlabel="Supply", ylabel="Acquisition",
    title="Frazer-Gilbert Functional Response (demand=$demand_fixed)")
ax4b = Axis(fig4[1, 2];
    xlabel="Supply/Demand Ratio", ylabel="φ (Supply-Demand Index)",
    title="Supply-Demand Index φ")

a_values = [0.3, 0.7, 1.5, 3.0]
sd_ratio_range = range(0.0, 5.0; length=300)

for a in a_values
    fr = FraserGilbertResponse(a)

    acq = [acquire(fr, S, demand_fixed) for S in supply_range]
    lines!(ax4a, collect(supply_range), acq; linewidth=2, label="a=$a")

    phi = [supply_demand_ratio(fr, r * demand_fixed, demand_fixed) for r in sd_ratio_range]
    lines!(ax4b, collect(sd_ratio_range), phi; linewidth=2, label="a=$a")
end

# Overlay the demand ceiling on acquisition plot
hlines!(ax4a, [demand_fixed]; color=:black, linestyle=:dash, linewidth=1, label="demand")
hlines!(ax4b, [1.0]; color=:black, linestyle=:dash, linewidth=1)

# Verify key properties
fr_test = FraserGilbertResponse(0.7)
@assert acquire(fr_test, 0.0, 100.0) == 0.0 "Zero supply → zero acquisition"
@assert acquire(fr_test, 100.0, 0.0) == 0.0 "Zero demand → zero acquisition"
acq_large = acquire(fr_test, 1e6, 100.0)
@assert acq_large > 99.9 "Large supply → acquisition ≈ demand"
println("  FG(a=0.7): acquire(0, 100)=$(acquire(fr_test, 0.0, 100.0))")
println("  FG(a=0.7): acquire(100, 100)=$(round(acquire(fr_test, 100.0, 100.0), digits=2))")
println("  FG(a=0.7): acquire(1e6, 100)=$(round(acq_large, digits=4)) → saturates at demand")
println("  FG(a=0.7): φ(50,100)=$(round(supply_demand_ratio(fr_test, 50.0, 100.0), digits=4))")

axislegend(ax4a; position=:rb)
axislegend(ax4b; position=:rb)
save(joinpath(figdir, "04_supply_demand_response.png"), fig4; px_per_unit=2)
println("  → Saved 04_supply_demand_response.png\n")

# ══════════════════════════════════════════════════════════════════════
# 5. Single Cohort Dynamics — 4-Stage Lifecycle at 25°C
# ══════════════════════════════════════════════════════════════════════
println("="^60)
println("5. Single Cohort Dynamics (egg→larva→pupa→adult at 25°C)")
println("="^60)

dev_rate_cohort = LinearDevelopmentRate(10.0, 35.0)
k_sub = 10

# CFL stability requires rate = k*dd/τ < 1 per substage per day.
# At 15°C, dd=5 DD/day. With k=10, need τ > 50 for each stage.
T_const = 15.0
dd_per_day = degree_days(dev_rate_cohort, T_const)  # 5.0

stage_params = [
    (:egg,   k_sub, 100.0,  500.0, 0.0),  # rate = 10*5/100 = 0.5
    (:larva, k_sub, 200.0,  0.0,   0.0),  # rate = 10*5/200 = 0.25
    (:pupa,  k_sub, 80.0,   0.0,   0.0),  # rate = 10*5/80  = 0.625
    (:adult, k_sub, 150.0,  0.0,   0.0),  # rate = 10*5/150 = 0.333
]

println("  CFL check (rate = k·dd/τ < 1 required):")
for (name, k, τ, _, _) in stage_params
    r = k * dd_per_day / τ
    println("    $name: rate = $(round(r, digits=3)) $(r < 1 ? "✓" : "✗ UNSTABLE")")
end

stages = [
    LifeStage(name, DistributedDelay(k, τ; W0=W0), dev_rate_cohort, μ)
    for (name, k, τ, W0, μ) in stage_params
]
pop_cohort = Population(:insect, stages)

n_sim_days = 120  # 120 × 5 = 600 DD (total τ = 530 DD)
weather_cohort = WeatherSeries(fill(T_const, n_sim_days); day_offset=1)

prob = PBDMProblem(pop_cohort, weather_cohort, (1, n_sim_days))
sol = solve(prob, DirectIteration())

# Extract stage trajectories
egg_traj   = stage_trajectory(sol, 1)
larva_traj = stage_trajectory(sol, 2)
pupa_traj  = stage_trajectory(sol, 3)
adult_traj = stage_trajectory(sol, 4)
total_traj = total_population(sol)
cdd = [0.0; cumsum(sol.degree_days)]

# Plot: top = stage dynamics vs DD, bottom = conservation check
fig5 = Figure(size=(700, 600))
ax5a = Axis(fig5[1, 1];
    xlabel="Cumulative Degree-Days", ylabel="Population",
    title="Cohort Progression Through 4 Stages (T=$(Int(T_const))°C, k=$k_sub)")
ax5b = Axis(fig5[2, 1];
    xlabel="Cumulative Degree-Days", ylabel="Total Population",
    title="Conservation Check (no mortality: total should be constant)")

lines!(ax5a, cdd, egg_traj; linewidth=2, label="Egg (τ=100 DD)")
lines!(ax5a, cdd, larva_traj; linewidth=2, label="Larva (τ=200 DD)")
lines!(ax5a, cdd, pupa_traj; linewidth=2, label="Pupa (τ=80 DD)")
lines!(ax5a, cdd, adult_traj; linewidth=2, label="Adult (τ=150 DD)")

# Mark expected stage transition DD
τ_cumulative = cumsum([100.0, 200.0, 80.0, 150.0])
for (i, τc) in enumerate(τ_cumulative)
    if τc <= cdd[end]
        vlines!(ax5a, [τc]; color=:gray, linestyle=:dot, linewidth=1)
    end
end

lines!(ax5b, cdd, total_traj; linewidth=2, color=:blue, label="Population remaining")
# Overlay conservation line: remaining + matured out
initial_pop = total_traj[1]
final_pop = total_traj[end]
cum_maturation = cumsum(sol.maturation)
total_exited = cum_maturation[end]
conserved = [total_traj[1]; [total_traj[i+1] + cum_maturation[i] for i in 1:length(cum_maturation)]]
lines!(ax5b, cdd, conserved; linewidth=2, color=:black, label="Pop + matured out")
hlines!(ax5b, [total_traj[1]]; color=:red, linestyle=:dash, linewidth=1,
    label="Initial ($(round(total_traj[1], digits=0)))")

# True conservation: pop_remaining + matured_out = initial (when μ=0)
conservation_err = abs(final_pop + total_exited - initial_pop) / initial_pop * 100
println("  DD per day at $(Int(T_const))°C: $dd_per_day")
println("  Total DD over $(n_sim_days) days: $(round(cdd[end], digits=1))")
println("  Initial total population: $(round(initial_pop, digits=1))")
println("  Final total population:   $(round(final_pop, digits=1))")
println("  Cumulative maturation out: $(round(total_exited, digits=1))")
println("  Pop + matured out:         $(round(final_pop + total_exited, digits=1))")
println("  Conservation error:        $(round(conservation_err, digits=4))%")
println("  Conservation $(conservation_err < 1.0 ? "PASSED ✓" : "check: $(round(conservation_err, digits=2))% error")")

# Stage peaks
for (i, (name, traj)) in enumerate(zip([:Egg, :Larva, :Pupa, :Adult],
                                        [egg_traj, larva_traj, pupa_traj, adult_traj]))
    peak_idx = argmax(traj)
    peak_dd = cdd[peak_idx]
    println("  $name: peak at DD=$(round(peak_dd, digits=1)), value=$(round(traj[peak_idx], digits=1))")
end

axislegend(ax5a; position=:rt)
axislegend(ax5b; position=:rb)
save(joinpath(figdir, "05_single_cohort_dynamics.png"), fig5; px_per_unit=2)
println("  → Saved 05_single_cohort_dynamics.png\n")

# ══════════════════════════════════════════════════════════════════════
# Summary
# ══════════════════════════════════════════════════════════════════════
println("="^60)
println("All 5 figures generated in: $figdir")
println("="^60)
for f in sort(readdir(figdir))
    println("  $f")
end
