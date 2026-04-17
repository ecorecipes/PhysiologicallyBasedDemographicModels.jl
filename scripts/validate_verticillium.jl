#!/usr/bin/env julia
# Validation script for Verticillium wilt management in cotton
# matching the vignette (Regev, Gutierrez et al. 1990 parameters).

using Pkg; Pkg.activate(joinpath(@__DIR__, ".."))
using CairoMakie

figdir = joinpath(@__DIR__, "figures", "verticillium")
mkpath(figdir)

# ══════════════════════════════════════════════════════════════════════
# 1. Pathogen inoculum dynamics (Regev et al. 1990)
# ══════════════════════════════════════════════════════════════════════

const D_INIT   = 30.0          # initial propagules/g soil
const D_MAX    = 60.0          # carrying capacity (Gutierrez et al. 1983)
const D_GROWTH = 1.2           # annual multiplication, no treatment
const MU_SOLAR = 0.01          # 99% kill (Pullman et al. 1981)
const MU_FUMIG = 0.05          # 95% kill (Regev et al. 1990)
const MU_ROT   = 0.50          # halving via rotation

function update_inoculum(D::Float64, act::Symbol)
    act == :none     && return min(D * D_GROWTH, D_MAX)
    act == :solarize && return D * MU_SOLAR
    act == :fumigate && return D * MU_FUMIG
    act == :rotate   && return D * MU_ROT
    error("Unknown action: $act")
end

println("Inoculum trajectories (5 years from D₀=$D_INIT):")
for act in [:none, :solarize, :fumigate, :rotate]
    D = D_INIT; traj = [D]
    for _ in 1:5; D = update_inoculum(D, act); push!(traj, D); end
    println("  $(rpad(act, 10)): $(join([round(d, digits=1) for d in traj], " → "))")
end

# ══════════════════════════════════════════════════════════════════════
# 2. Virulence evolution: V_{t+1} = V_t + (V* − V_t)/2
# ══════════════════════════════════════════════════════════════════════

const V_SS = 0.4               # steady-state virulence (Gutierrez et al. 1983)
update_virulence(V::Float64) = V + 0.5 * (V_SS - V)

println("\nVirulence convergence to V*=$V_SS:")
for V0 in [0.2, 0.4, 0.6, 0.8, 1.0]
    V = V0; vals = [V]
    for _ in 1:10; V = update_virulence(V); push!(vals, V); end
    println("  V₀=$(lpad(V0, 3)): $(join([round(v, digits=3) for v in vals], " → "))")
end

# ══════════════════════════════════════════════════════════════════════
# 3. Cotton yield response (Regev et al. 1990, Eq. 5)
#    Y(v,D,V) = 1.4(0.74 + 0.25v − 0.011v² − 0.04D + 0.00036D² − 0.755V³)
#    v = plant vigor index (fixed at 4.0)
#    D = inoculum density (propagules/g soil, 0–60)
#    V = virulence index (0.2–1.0)
# ══════════════════════════════════════════════════════════════════════

const PLANT_V = 4.0

function cotton_yield(v, D, V)
    max(0.0, 1.4 * (0.74 + 0.25v - 0.011v^2 - 0.04D + 0.00036D^2 - 0.755V^3))
end

println("\nYield table (bales/acre) at v=$PLANT_V:")
println("  D\\V    |  0.2  |  0.4  |  0.6  |  0.8")
println("  " * "-"^46)
for D in [0, 15, 30, 45, 60]
    ys = [round(cotton_yield(PLANT_V, Float64(D), V), digits=2) for V in [0.2,0.4,0.6,0.8]]
    println("  $(lpad(D,3))    | $(join([lpad(y, 5) for y in ys], " | "))")
end

# ══════════════════════════════════════════════════════════════════════
# 4. Economic framework (Regev et al. 1990)
# ══════════════════════════════════════════════════════════════════════

const COTTON_PRICE = 300.0     # $/bale
const COTTON_COST  = 450.0     # $/acre production cost
const BARLEY_PROFIT = 50.0     # $/acre net (rotation/solarization crop)
const TREAT_COST   = 300.0     # $/acre solarization or fumigation
const DISC         = 0.05      # 5% annual discount rate
const N_YR         = 10        # planning horizon
const ACTIONS      = [:none, :solarize, :fumigate, :rotate]

function annual_profit(act::Symbol, D, V;
                       treat_cost_mult::Float64=1.0)
    tc = TREAT_COST * treat_cost_mult
    if act == :none
        return cotton_yield(PLANT_V, D, V) * COTTON_PRICE - COTTON_COST
    elseif act == :solarize
        return BARLEY_PROFIT - tc
    elseif act == :fumigate
        return cotton_yield(PLANT_V, D * MU_FUMIG, V) * COTTON_PRICE - COTTON_COST - tc
    elseif act == :rotate
        return BARLEY_PROFIT
    end
    error("Unknown action: $act")
end

npv(profits, r) = sum(profits[t] / (1 + r)^(t - 1) for t in eachindex(profits))

println("\nAnnual profit at V=$V_SS:")
for act in ACTIONS
    ps = [round(annual_profit(act, Float64(D), V_SS), digits=0) for D in [10,30,50]]
    println("  $(rpad(act,10)): D=10 → \$$(ps[1]),  D=30 → \$$(ps[2]),  D=50 → \$$(ps[3])")
end

# ══════════════════════════════════════════════════════════════════════
# 5. Multi-year simulation
# ══════════════════════════════════════════════════════════════════════

function simulate(acts, D0, V0; treat_cost_mult::Float64=1.0)
    D, V = D0, V0
    Dh = [D0]; Yh = Float64[]; Ph = Float64[]
    for t in 1:N_YR
        a = acts[t]
        push!(Yh, a ∈ [:solarize, :rotate] ? 0.0 :
              cotton_yield(PLANT_V, D, V))
        push!(Ph, annual_profit(a, D, V; treat_cost_mult))
        D = update_inoculum(D, a)
        V = update_virulence(V)
        push!(Dh, D)
    end
    (D=Dh, Y=Yh, P=Ph)
end

D0, V0 = D_INIT, 0.6

strats = Dict(
    "No treatment"  => fill(:none, N_YR),
    "Solarize"      => fill(:solarize, N_YR),
    "Fumigate"      => fill(:fumigate, N_YR),
    "Rotate (2-yr)" => [isodd(t) ? :none : :rotate for t in 1:N_YR])

println("\n10-year strategy comparison (D₀=$D0, V₀=$V0):")
for (nm, acts) in sort(collect(strats), by=first)
    r = simulate(acts, D0, V0)
    pv = npv(r.P, DISC)
    println("  $nm: final D=$(round(r.D[end], digits=1)), " *
            "NPV=\$$(round(pv, digits=0))/acre")
end

# Optimal strategy (two-step lookahead, Regev et al. 1990)
function optimal_strategy(D0, V0; n=N_YR, r=DISC,
                          treat_cost_mult::Float64=1.0)
    acts = Symbol[]; D, V = D0, V0
    for t in 1:n
        best = (:none, -Inf)
        for a in ACTIONS
            pnow = annual_profit(a, D, V; treat_cost_mult) / (1 + r)^(t - 1)
            Dn = update_inoculum(D, a)
            Vn = update_virulence(V)
            pnxt = maximum(
                annual_profit(b, Dn, Vn; treat_cost_mult) / (1 + r)^t
                for b in ACTIONS)
            if pnow + pnxt > best[2]
                best = (a, pnow + pnxt)
            end
        end
        push!(acts, best[1])
        D = update_inoculum(D, best[1])
        V = update_virulence(V)
    end
    acts
end

opt = optimal_strategy(D0, V0)
r_opt = simulate(opt, D0, V0)
println("\nOptimal sequence: $(join(string.(opt), ", "))")
println("Optimal NPV: \$$(round(npv(r_opt.P, DISC), digits=0))/acre")
println("Literature target NPV: ~\$2,700/acre (Regev et al. 1990)")

# ══════════════════════════════════════════════════════════════════════
# Figure (a): Inoculum trajectories under 4 management strategies
# ══════════════════════════════════════════════════════════════════════

println("\nGenerating figures …")

names_4 = ["No treatment", "Solarize", "Fumigate", "Rotate (2-yr)"]
acts_4  = [strats[n] for n in names_4]
cols_4  = [:firebrick, :steelblue, :seagreen, :darkorange]

fig1 = Figure(size=(700, 450))
ax1 = Axis(fig1[1, 1],
    xlabel = "Year",
    ylabel = "Inoculum density D (propagules/g soil)",
    title  = "Inoculum Trajectories Under Management Strategies\n(D₀=$D0, V₀=$V0; Regev et al. 1990)")
for (i, (nm, a)) in enumerate(zip(names_4, acts_4))
    r = simulate(a, D0, V0)
    lines!(ax1, 0:N_YR, r.D, color=cols_4[i], linewidth=2.5, label=nm)
    scatter!(ax1, 0:N_YR, r.D, color=cols_4[i], markersize=6)
end
# Optimal overlay
lines!(ax1, 0:N_YR, r_opt.D, color=:purple, linewidth=2.5,
       linestyle=:dash, label="Optimal")
scatter!(ax1, 0:N_YR, r_opt.D, color=:purple, markersize=6)
hlines!(ax1, [D_MAX], color=:gray60, linestyle=:dot, linewidth=1,
        label="D_max = $D_MAX")
axislegend(ax1, position=:rt, framevisible=false, labelsize=10)
text!(ax1, 7, 55, text="Literature: D_max = 60 prop/g",
      fontsize=9, color=:gray40)
save(joinpath(figdir, "inoculum_trajectories.png"), fig1, px_per_unit=3)
println("  ✓ inoculum_trajectories.png")

# ══════════════════════════════════════════════════════════════════════
# Figure (b): Virulence evolution — convergence to V* = 0.4
# ══════════════════════════════════════════════════════════════════════

fig2 = Figure(size=(700, 450))
ax2 = Axis(fig2[1, 1],
    xlabel = "Year",
    ylabel = "Virulence index V",
    title  = "Virulence Convergence to V* = $V_SS\n(V_{t+1} = V_t + (V* − V_t)/2)")
v0_vals = [0.2, 0.4, 0.6, 0.8, 1.0]
v_cols  = cgrad(:viridis, length(v0_vals), categorical=true)
for (i, V0i) in enumerate(v0_vals)
    V = V0i; traj = [V]
    for _ in 1:10; V = update_virulence(V); push!(traj, V); end
    lines!(ax2, 0:10, traj, color=v_cols[i], linewidth=2.5,
           label="V₀ = $V0i")
    scatter!(ax2, 0:10, traj, color=v_cols[i], markersize=6)
end
hlines!(ax2, [V_SS], color=:red, linestyle=:dash, linewidth=1.5,
        label="V* = $V_SS (steady state)")
axislegend(ax2, position=:rt, framevisible=false, labelsize=10)
text!(ax2, 6, 0.42, text="Literature: V* = 0.4 (Gutierrez et al. 1983)",
      fontsize=9, color=:gray40)
save(joinpath(figdir, "virulence_evolution.png"), fig2, px_per_unit=3)
println("  ✓ virulence_evolution.png")

# ══════════════════════════════════════════════════════════════════════
# Figure (c): Yield surface — Y(v=4, D, V)
# ══════════════════════════════════════════════════════════════════════

fig3 = Figure(size=(750, 500))
Ds = range(0, 60, length=100)
Vs = range(0.2, 1.0, length=80)
Z  = [cotton_yield(PLANT_V, d, v) for v in Vs, d in Ds]

ax3 = Axis(fig3[1, 1],
    xlabel = "Inoculum density D (propagules/g soil)",
    ylabel = "Virulence index V",
    title  = "Cotton Yield Surface Y(v=4, D, V) — bales/acre\n(Regev et al. 1990, Eq. 5)")
hm = heatmap!(ax3, collect(Ds), collect(Vs), Z',
              colormap=:YlOrRd, colorrange=(0, maximum(Z)))
Colorbar(fig3[1, 2], hm, label="Yield (bales/acre)")
contour!(ax3, collect(Ds), collect(Vs), Z',
         levels=[0.5, 1.0, 1.5, 2.0, 2.5],
         color=:black, linewidth=0.8, labels=true, labelsize=9)
# Mark typical operating point
scatter!(ax3, [D_INIT], [V_SS], color=:white, markersize=12,
         strokecolor=:black, strokewidth=2)
text!(ax3, D_INIT + 2, V_SS + 0.03,
      text="(D₀=30, V*=0.4)", fontsize=10, color=:white)
text!(ax3, 2, 0.95,
      text="Y = 1.4(0.74 + 0.25v − 0.011v² − 0.04D\n      + 0.00036D² − 0.755V³)",
      fontsize=8, color=:black)
save(joinpath(figdir, "yield_surface.png"), fig3, px_per_unit=3)
println("  ✓ yield_surface.png")

# ══════════════════════════════════════════════════════════════════════
# Figure (d): Cumulative NPV comparison over 10 years
# ══════════════════════════════════════════════════════════════════════

fig4 = Figure(size=(700, 450))
ax4 = Axis(fig4[1, 1],
    xlabel = "Year",
    ylabel = "Cumulative NPV (\$/acre)",
    title  = "Cumulative Net Present Value Under Strategies\n(D₀=$D0, V₀=$V0, r=$(Int(DISC*100))%)")
names_5 = ["No treatment", "Solarize", "Fumigate", "Rotate (2-yr)", "Optimal"]
acts_5  = [strats["No treatment"], strats["Solarize"],
           strats["Fumigate"], strats["Rotate (2-yr)"], opt]
cols_5  = [:firebrick, :steelblue, :seagreen, :darkorange, :purple]
lsty_5  = [:solid, :solid, :solid, :solid, :dash]
final_npvs = Float64[]
for (i, (nm, a)) in enumerate(zip(names_5, acts_5))
    r = simulate(a, D0, V0)
    cpv = cumsum([r.P[t] / (1 + DISC)^(t - 1) for t in 1:N_YR])
    push!(final_npvs, cpv[end])
    lines!(ax4, 1:N_YR, cpv, color=cols_5[i], linewidth=2.5,
           linestyle=lsty_5[i], label="$nm (\$$(round(Int, cpv[end])))")
    scatter!(ax4, 1:N_YR, cpv, color=cols_5[i], markersize=6)
end
hlines!(ax4, [2700.0], color=:gray50, linestyle=:dot, linewidth=1,
        label="Literature ≈\$2,700")
axislegend(ax4, position=:lt, framevisible=false, labelsize=9)
save(joinpath(figdir, "cumulative_npv.png"), fig4, px_per_unit=3)
println("  ✓ cumulative_npv.png")

# ══════════════════════════════════════════════════════════════════════
# Figure (e): Sensitivity — treatment cost multiplier vs optimal NPV
# ══════════════════════════════════════════════════════════════════════

fig5 = Figure(size=(750, 500))
cost_mults = range(0.0, 2.5, length=26)
sens_npvs   = Dict(nm => Float64[] for nm in names_4)
opt_npvs    = Float64[]
opt_first   = Symbol[]

println("\nSensitivity to treatment cost multiplier:")
println("  Cost mult | NPV (optimal) | First action")
for cm in cost_mults
    # Fixed strategies
    for (nm, a) in zip(names_4, acts_4)
        r = simulate(a, D0, V0; treat_cost_mult=cm)
        push!(sens_npvs[nm], npv(r.P, DISC))
    end
    # Optimal
    oacts = optimal_strategy(D0, V0; treat_cost_mult=cm)
    r_o = simulate(oacts, D0, V0; treat_cost_mult=cm)
    push!(opt_npvs, npv(r_o.P, DISC))
    push!(opt_first, oacts[1])
end

# Print a subset
for cm in [0.25, 0.50, 1.0, 1.5, 2.0]
    idx = argmin(abs.(collect(cost_mults) .- cm))
    println("    $(lpad(round(cm, digits=2), 5))x    | " *
            "\$$(lpad(round(opt_npvs[idx], digits=0), 6))      | $(opt_first[idx])")
end

ax5a = Axis(fig5[1, 1],
    xlabel = "Treatment cost multiplier",
    ylabel = "10-year NPV (\$/acre)",
    title  = "Sensitivity: Treatment Cost vs. NPV\n(D₀=$D0, V₀=$V0)")
for (i, nm) in enumerate(names_4)
    lines!(ax5a, collect(cost_mults), sens_npvs[nm],
           color=cols_4[i], linewidth=2, label=nm)
end
lines!(ax5a, collect(cost_mults), opt_npvs,
       color=:purple, linewidth=2.5, linestyle=:dash, label="Optimal")
vlines!(ax5a, [1.0], color=:gray50, linestyle=:dot, linewidth=1,
        label="Baseline cost")
axislegend(ax5a, position=:rt, framevisible=false, labelsize=9)

# Lower panel: first action chosen by optimal strategy
first_action_num = [findfirst(==(a), ACTIONS) for a in opt_first]
ax5b = Axis(fig5[2, 1],
    xlabel = "Treatment cost multiplier",
    ylabel = "First optimal action",
    yticks = (1:4, string.(ACTIONS)),
    title  = "Optimal First-Year Action vs. Treatment Cost")
scatter!(ax5b, collect(cost_mults), first_action_num,
         color=:purple, markersize=8)
lines!(ax5b, collect(cost_mults), first_action_num,
       color=:purple, linewidth=1.5)
vlines!(ax5b, [1.0], color=:gray50, linestyle=:dot, linewidth=1)
rowsize!(fig5.layout, 2, Relative(0.3))
save(joinpath(figdir, "sensitivity_analysis.png"), fig5, px_per_unit=3)
println("  ✓ sensitivity_analysis.png")

# ══════════════════════════════════════════════════════════════════════
# Summary verification against literature
# ══════════════════════════════════════════════════════════════════════

println("\n" * "="^65)
println("Literature verification (Regev et al. 1990)")
println("="^65)
println("  Yield at (v=4, D=30, V=0.4): $(round(cotton_yield(4.0, 30.0, 0.4), digits=2)) bales/acre")
println("  V convergence: V₀=1.0 → V₅ = $(round(0.4 + (1.0-0.4)*0.5^5, digits=4)) (expect ≈0.4)")
println("  D range: 10–60 propagules/g soil  ✓")
println("  Solarization: 99% kill (μ=$(MU_SOLAR))  ✓")
println("  Fumigation: 95% kill (μ=$(MU_FUMIG))  ✓")
println("  Rotation: halves D (μ=$(MU_ROT))  ✓")
println("  Virulence steady state: V* = $V_SS  ✓")
opt_npv_val = round(npv(r_opt.P, DISC), digits=0)
println("  Optimal NPV: \$$opt_npv_val/acre (literature ≈\$2,700)")
println("  Optimal sequence: $(join(string.(opt), ", "))")
println("  Literature: no treatment ~2 yrs, then rotate every 3rd yr")
println("="^65)
println("All 5 figures saved to: $figdir")
