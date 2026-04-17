#!/usr/bin/env julia
# Validation script for Rice–Fish Agroecosystem (Vignette 33).
#
# System: Oryza sativa (rice), phytoplankton, zooplankton, Nile tilapia
#         (Oreochromis niloticus) in flooded paddies at Bang Sai, Thailand.
# Literature: d'Oultremont & Gutierrez (2002), Parts I & II, Ecological Modelling.
#
# Key targets from the literature:
#   - Tilapia final mass: 270 g in 150 days (standard), ~100 g (floodwater),
#     270 g (pellet-fed fishpond)
#   - Rice yield: ~14 g/hill (standard), ~21 g/hill (fertilized),
#     ~15.4 g/hill (pellet-fed fish scenario)
#   - Tillers: start 3, peak 19–23 at day 32–35, decline to 8–10
#   - Area ratio: 90 % rice, 10 % pond; fish spend 80 % of time in rice
#   - Rice canopy shading is master variable limiting aquatic productivity
#
# The rice crop is modelled with the Graf et al. (1990) carbon-balance approach:
# organ masses accumulate via photosynthate allocation with phenological
# priority shifts driven by cumulative degree-days.  The vignette's
# DistributedDelay types (designed for insect stage-structured flow) are used
# only for the Q10Respiration and development-rate helpers from the PBDM
# library; organ masses are tracked as simple state variables.

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using PhysiologicallyBasedDemographicModels
using CairoMakie

figdir = joinpath(@__DIR__, "figures", "rice_fish")
mkpath(figdir)

# ============================================================================
# Rice parameters (Yoshida 1981; Graf et al. 1990)
# ============================================================================
const RICE_T_BASE = 10.0
const RICE_T_UPPER = 42.0
rice_dev = LinearDevelopmentRate(RICE_T_BASE, RICE_T_UPPER)

const DD_MAX_TILLER  = 500.0   # peak tillering
const DD_PANICLE_INIT = 1000.0 # panicle initiation
const DD_GRAIN_FILL  = 1200.0  # grain filling onset
const DD_HEADING     = 1500.0  # heading / anthesis
const DD_HARVEST     = 2000.0  # physiological maturity
const SLA_RICE = 0.037         # specific leaf area (m²/g)
const EXT_K = 0.6              # canopy extinction coefficient
const RAD_CONV = 0.65          # radiation use efficiency (g DM / MJ PAR)
const HILLS = 25               # planting density (hills/m²)
const GROWTH_RESP = 0.28       # growth respiration fraction
const TAU_MIN = 0.05           # minimum canopy transmittance (gap fraction)

# Maintenance respiration per organ (Q₁₀ model)
resp_culm  = Q10Respiration(0.010, 2.0, 25.0)
resp_leaf  = Q10Respiration(0.015, 2.0, 25.0)
resp_root  = Q10Respiration(0.008, 2.0, 25.0)
resp_grain = Q10Respiration(0.005, 2.0, 25.0)

# ============================================================================
# Rice sub-models
# ============================================================================

"""Tiller count per hill as function of cumulative degree-days."""
function tiller_count(dd)
    dd < 150 && return 3.0
    dd < DD_MAX_TILLER && return 3.0 + 18.0 * (dd - 150) / (DD_MAX_TILLER - 150)
    return 21.0 - 12.0 * min(1.0, (dd - DD_MAX_TILLER) / (DD_HARVEST - DD_MAX_TILLER))
end

"""Canopy light interception → gross photosynthesis (g DM/hill/day)."""
function rice_photo(lm, rad, til)
    lai = min(10.0, lm * SLA_RICE * HILLS * til / 3.0)
    (1.0 - exp(-EXT_K * lai)) * rad * RAD_CONV / HILLS
end

"""Fraction of above-canopy light reaching the floodwater."""
function canopy_transmittance(lm, til)
    max(TAU_MIN, exp(-EXT_K * min(10.0, lm * SLA_RICE * HILLS * til / 3.0)))
end

"""
Allocation fractions (culm, leaf, root, grain) by phenological phase.
Based on Graf et al. (1990) priority-shift scheme.
"""
function allocation_fractions(cumd)
    if cumd < DD_PANICLE_INIT
        # Vegetative: invest in leaves and roots for canopy establishment
        return (culm=0.25, leaf=0.40, root=0.35, grain=0.00)
    elseif cumd < DD_GRAIN_FILL
        # Reproductive: begin shifting to grain
        return (culm=0.20, leaf=0.20, root=0.10, grain=0.50)
    elseif cumd < DD_HEADING
        # Grain filling: strong grain sink
        return (culm=0.10, leaf=0.05, root=0.05, grain=0.80)
    else
        # Maturation / ripening: nearly all to grain
        return (culm=0.05, leaf=0.02, root=0.03, grain=0.90)
    end
end

"""
Leaf senescence rate (g/day).  Accelerates after grain fill as N
translocates to developing grains.
"""
function leaf_senescence(lm, cumd)
    cumd < DD_GRAIN_FILL && return lm * 0.005
    cumd < DD_HEADING    && return lm * 0.015
    return lm * 0.025
end

# ============================================================================
# Aquatic food web parameters (d'Oultremont & Gutierrez 2002)
# ============================================================================
const PHYTO_RMAX = 1.5    # max growth rate (/day) at 30 °C
const PHYTO_TOPT = 30.0
const PHYTO_LOSS = 0.15   # sinking + non-grazing mortality
const PHYTO_K_LIGHT = 0.3 # Monod half-saturation for light
const PHYTO_CAP = 50.0    # carrying capacity (mg C/L)

const ZOO_RMAX = 0.9      # max growth rate (/day) at 28 °C
const ZOO_TOPT = 28.0
const ZOO_LOSS = 0.10     # non-predation mortality
const ZOO_HALF = 5.0      # half-saturation for phyto grazing
const ZOO_GMAX = 2.0      # max grazing rate
const ZOO_EFF  = 0.30     # assimilation efficiency

phyto_r(T, lf) = T < 5 ? 0.0 :
    PHYTO_RMAX * max(0, 1 - ((T - PHYTO_TOPT) / 15)^2) * lf / (PHYTO_K_LIGHT + lf)
zoo_r(T, P) = T < 8 ? 0.0 :
    ZOO_RMAX * max(0, 1 - ((T - ZOO_TOPT) / 12)^2) * P / (ZOO_HALF + P)

# ============================================================================
# Tilapia parameters (Caulton 1982; d'Oultremont & Gutierrez 2002)
# ============================================================================
const TIL_TBASE = 16.0
const TIL_TOPT  = 28.0
const TIL_TMAX  = 36.0
const TIL_W0    = 12.0    # stocking weight (g)
const TIL_DENS  = 2.0     # stocking density (fish/m²)
const TIL_GMAX  = 0.05    # max specific growth rate (/day) [calibrated]
tilapia_resp = Q10Respiration(0.012, 2.2, 25.0)

function til_tscalar(T)
    (T < TIL_TBASE || T > TIL_TMAX) && return 0.0
    T <= TIL_TOPT ? min(1.0, (T - TIL_TBASE) / (TIL_TOPT - TIL_TBASE)) :
                    max(0.0, 1 - (T - TIL_TOPT) / (TIL_TMAX - TIL_TOPT))
end

# ============================================================================
# System coupling constants
# ============================================================================
const A_RICE       = 0.90  # fraction of area in rice
const FISH_T_RICE  = 0.80  # fraction of time fish spend in rice area
const NUTRIENT_BOOST = 0.10 # N recycling from fish feces → rice yield
const PELLET_RATE  = 0.025 # pellet feed rate (g feed / g fish / day)
const BENTHIC_BASE = 2.0   # persistent benthic/detrital food (simplified from
                            # 14-species model: periphyton, detritus, benthos)

# ============================================================================
# Weather: Bang Sai, Thailand (150-day wet-season crop)
# ============================================================================
const NDAYS = 150
wdays = [let T_m = 28.0 - 1.5d / NDAYS,
             r = max(12.0, 20 + 3sin(2π * d / 45) - 2d / NDAYS)
    DailyWeather(T_m, T_m - 5 + sin(2π * d / 30), T_m + 5 + 1.5sin(2π * d / 30);
                 radiation=r, photoperiod=12.5)
end for d in 1:NDAYS]
weather = WeatherSeries(wdays; day_offset=1)

# ============================================================================
# Simulation engine — carbon-balance approach (Graf et al. 1990)
# ============================================================================
function sim_rf(weather, nd; fish=true, pond=false, pellet=false)
    # Rice organ masses (g DM / hill)
    cm  = zeros(nd + 1)  # culm
    lm  = zeros(nd + 1)  # leaf
    rm_ = zeros(nd + 1)  # root
    gm  = zeros(nd + 1)  # grain
    # Initial transplant seedling mass per hill (3–4 week old seedlings)
    cm[1] = 1.5; lm[1] = 1.5; rm_[1] = 0.8; gm[1] = 0.0

    # Aquatic state
    P  = zeros(nd + 1)   # phytoplankton (mg C/L)
    Z  = zeros(nd + 1)   # zooplankton  (mg C/L)
    fw = zeros(nd + 1)   # fish mass (g/fish)
    P[1] = 8.0; Z[1] = 2.0
    fw[1] = fish ? TIL_W0 : 0.0

    # Tracking arrays
    tr  = zeros(nd)       # canopy transmittance
    cdd = zeros(nd)       # cumulative degree-days
    til = zeros(nd + 1)   # tiller count
    til[1] = 3.0

    cumd = 0.0

    for d in 1:nd
        w = get_weather(weather, min(d, nd))
        T = w.T_mean
        dd = degree_days(rice_dev, T)
        cumd += dd
        cdd[d] = cumd

        tl = tiller_count(cumd)
        til[d+1] = tl

        # ---- Rice carbon balance ----
        gross = rice_photo(lm[d], w.radiation, tl)
        maint = respiration_rate(resp_culm, T)  * cm[d] +
                respiration_rate(resp_leaf, T)  * lm[d] +
                respiration_rate(resp_root, T)  * rm_[d] +
                respiration_rate(resp_grain, T) * gm[d]
        ns = max(0.0, gross - maint) * (1 - GROWTH_RESP)
        # Fish nutrient recycling boost
        pond && pellet && (ns *= 1 + NUTRIENT_BOOST)

        af = allocation_fractions(cumd)
        senes = leaf_senescence(lm[d], cumd)

        cm[d+1]  = max(0.0, cm[d]  + ns * af.culm)
        lm[d+1]  = max(0.0, lm[d]  + ns * af.leaf - senes)
        rm_[d+1] = max(0.0, rm_[d] + ns * af.root)
        gm[d+1]  = max(0.0, gm[d]  + ns * af.grain)

        # Canopy transmittance for aquatic web
        τ = canopy_transmittance(lm[d+1], tl)
        tr[d] = τ

        # ---- Aquatic food web ----
        Pv = P[d]; Zv = Z[d]
        gr = ZOO_GMAX * Zv * Pv / (ZOO_HALF + Pv)
        dP = Pv * (phyto_r(T, τ) * (1 - Pv / PHYTO_CAP) - PHYTO_LOSS) - gr

        fp = 0.0
        if fish
            ft = pond ? FISH_T_RICE : 1.0
            fp = 0.5 * til_tscalar(T) * ft * TIL_DENS * (fw[d] / 100) * Zv / (2 + Zv)
        end
        dZ = Zv * (zoo_r(T, Pv) * ZOO_EFF - ZOO_LOSS) - fp
        P[d+1] = max(0.01, Pv + dP)
        Z[d+1] = max(0.01, Zv + dZ)

        # ---- Tilapia growth ----
        if fish
            W = fw[d]; ts = til_tscalar(T)
            food = Z[d] * 0.5 + P[d] * 0.05 + BENTHIC_BASE
            pond && pellet && (food += PELLET_RATE * W)
            g = min(TIL_GMAX * W * ts, food * TIL_DENS * 0.35)
            fw[d+1] = max(W * 0.95, W + g - respiration_rate(tilapia_resp, T) * W)
        end
    end

    (; culm=cm, leaf=lm, root=rm_, grain=gm, tiller=til,
       phyto=P, zoo=Z, fish=fw, transmit=tr, cum_dd=cdd)
end

# ============================================================================
# Run all three scenarios
# ============================================================================
println("Running rice-fish agroecosystem scenarios...")
r1 = sim_rf(weather, NDAYS; fish=false)
r2 = sim_rf(weather, NDAYS; fish=true, pond=false)
r3 = sim_rf(weather, NDAYS; fish=true, pond=true, pellet=true)

println("── Scenario 1: Rice Only ──")
println("  Grain: $(round(r1.grain[end], digits=2)) g/hill  " *
        "($(round(r1.grain[end] * HILLS / 1000, digits=2)) t/ha equiv.)")
println("  Peak tiller: $(round(maximum(r1.tiller), digits=1))")
println("  Peak phyto:  $(round(maximum(r1.phyto), digits=1)) mg C/L")

println("── Scenario 2: Rice + Tilapia (floodwater) ──")
println("  Grain: $(round(r2.grain[end], digits=2)) g/hill")
println("  Fish:  $(TIL_W0) → $(round(r2.fish[end], digits=1)) g")

println("── Scenario 3: Rice + Fishpond (pellet-fed) ──")
println("  Grain: $(round(r3.grain[end], digits=2)) g/hill")
println("  Fish:  $(TIL_W0) → $(round(r3.fish[end], digits=1)) g")

days = 0:NDAYS

# ============================================================================
# Figure 1: Rice growth curve — tillering and dry-matter allocation
# ============================================================================
fig1 = Figure(size=(1000, 700))

ax1a = Axis(fig1[1, 1],
    xlabel="Days after transplanting",
    ylabel="Tillers per hill",
    title="(a) Tiller Dynamics")
lines!(ax1a, collect(days), r1.tiller, linewidth=2.5, color=:forestgreen, label="Simulated")
hlines!(ax1a, [3.0], color=:grey60, linestyle=:dot, linewidth=1)
hlines!(ax1a, [21.0], color=:grey60, linestyle=:dot, linewidth=1)
hlines!(ax1a, [9.0], color=:grey60, linestyle=:dot, linewidth=1)
text!(ax1a, 5, 3.8; text="Initial: 3", fontsize=11, color=:grey40)
peak_til = maximum(r1.tiller)
peak_day = argmax(r1.tiller) - 1
text!(ax1a, peak_day + 3, peak_til + 0.5;
      text="Peak: $(round(peak_til, digits=1)) @ day $peak_day\n(Lit: 19–23 @ day 32–35)",
      fontsize=11, color=:blue)
text!(ax1a, 110, r1.tiller[end] + 0.8;
      text="Final: $(round(r1.tiller[end], digits=1))\n(Lit: 8–10)",
      fontsize=11, color=:blue)
axislegend(ax1a; position=:rt)

ax1b = Axis(fig1[1, 2],
    xlabel="Days after transplanting",
    ylabel="Dry matter (g/hill)",
    title="(b) Dry Matter Allocation by Organ")
lines!(ax1b, collect(days), r1.leaf, linewidth=2, color=:green, label="Leaf")
lines!(ax1b, collect(days), r1.culm, linewidth=2, color=:saddlebrown, label="Culm")
lines!(ax1b, collect(days), r1.root, linewidth=2, color=:sienna, linestyle=:dash, label="Root")
lines!(ax1b, collect(days), r1.grain, linewidth=2.5, color=:gold, label="Grain")
total_dm = r1.leaf .+ r1.culm .+ r1.root .+ r1.grain
lines!(ax1b, collect(days), total_dm, linewidth=1.5, color=:black, linestyle=:dot, label="Total")
axislegend(ax1b; position=:lt)

ax1c = Axis(fig1[2, 1:2],
    xlabel="Days after transplanting",
    ylabel="Cumulative degree-days (°C·d)",
    title="(c) Phenological Timeline")
lines!(ax1c, collect(1:NDAYS), r1.cum_dd, linewidth=2, color=:firebrick)
hlines!(ax1c, [DD_MAX_TILLER], color=:grey60, linestyle=:dash, linewidth=1)
hlines!(ax1c, [DD_PANICLE_INIT], color=:grey60, linestyle=:dash, linewidth=1)
hlines!(ax1c, [DD_GRAIN_FILL], color=:grey60, linestyle=:dash, linewidth=1)
hlines!(ax1c, [DD_HARVEST], color=:grey60, linestyle=:dash, linewidth=1)
text!(ax1c, 5, DD_MAX_TILLER + 40; text="Max tiller (500 DD)", fontsize=10, color=:grey40)
text!(ax1c, 5, DD_PANICLE_INIT + 40; text="Panicle init (1000 DD)", fontsize=10, color=:grey40)
text!(ax1c, 5, DD_GRAIN_FILL + 40; text="Grain fill (1200 DD)", fontsize=10, color=:grey40)
text!(ax1c, 5, DD_HARVEST + 40; text="Harvest (2000 DD)", fontsize=10, color=:grey40)

save(joinpath(figdir, "rice_growth_curve.png"), fig1, px_per_unit=2)
println("Saved rice_growth_curve.png")

# ============================================================================
# Figure 2: Phytoplankton dynamics — light-dependent growth & canopy shading
# ============================================================================
fig2 = Figure(size=(1000, 600))

ax2a = Axis(fig2[1, 1],
    xlabel="Days after transplanting",
    ylabel="Canopy transmittance (fraction)",
    title="(a) Rice Canopy Light Transmittance")
lines!(ax2a, collect(1:NDAYS), r1.transmit, linewidth=2.5, color=:darkorange, label="Rice only")
lines!(ax2a, collect(1:NDAYS), r2.transmit, linewidth=2, color=:darkorange, linestyle=:dash, label="Rice+fish")
hlines!(ax2a, [0.5], color=:grey60, linestyle=:dot, linewidth=1)
text!(ax2a, 100, 0.53; text="50 % threshold", fontsize=10, color=:grey40)
axislegend(ax2a; position=:rt)

ax2b = Axis(fig2[1, 2],
    xlabel="Days after transplanting",
    ylabel="Biomass (mg C/L)",
    title="(b) Phytoplankton: Canopy Shading Effect")
lines!(ax2b, collect(days), r1.phyto, linewidth=2.5, color=:seagreen, label="No fish")
lines!(ax2b, collect(days), r2.phyto, linewidth=2, color=:seagreen, linestyle=:dash, label="Floodwater fish")
lines!(ax2b, collect(days), r3.phyto, linewidth=2, color=:teal, linestyle=:dashdot, label="Fishpond+pellet")
text!(ax2b, 5, maximum(r1.phyto) * 0.85;
      text="Shading limits phyto\nas canopy closes",
      fontsize=10, color=:grey30)
axislegend(ax2b; position=:rt)

ax2c = Axis(fig2[2, 1:2],
    xlabel="Days after transplanting",
    ylabel="Biomass (mg C/L)",
    title="(c) Zooplankton Response — Trophic Cascade")
lines!(ax2c, collect(days), r1.zoo, linewidth=2, color=:steelblue, label="Zoo (no fish)")
lines!(ax2c, collect(days), r2.zoo, linewidth=2, color=:steelblue, linestyle=:dash, label="Zoo (floodwater)")
lines!(ax2c, collect(days), r3.zoo, linewidth=2, color=:royalblue, linestyle=:dashdot, label="Zoo (fishpond)")
text!(ax2c, 80, maximum(r1.zoo) * 0.7;
      text="Fish grazing depresses zooplankton\n→ cascading from canopy shading",
      fontsize=10, color=:grey30)
axislegend(ax2c; position=:rt)

save(joinpath(figdir, "phytoplankton_dynamics.png"), fig2, px_per_unit=2)
println("Saved phytoplankton_dynamics.png")

# ============================================================================
# Figure 3: Tilapia growth trajectories
# ============================================================================
fig3 = Figure(size=(900, 600))

ax3a = Axis(fig3[1, 1],
    xlabel="Days after stocking",
    ylabel="Body mass (g/fish)",
    title="(a) Tilapia Growth: Floodwater vs Fishpond")
lines!(ax3a, collect(days), r2.fish, linewidth=2.5, color=:royalblue, label="Floodwater only")
lines!(ax3a, collect(days), r3.fish, linewidth=2.5, color=:firebrick, label="Fishpond + pellets")
# Literature target annotations
hlines!(ax3a, [270.0], color=:firebrick, linestyle=:dot, linewidth=1)
hlines!(ax3a, [100.0], color=:royalblue, linestyle=:dot, linewidth=1)
scatter!(ax3a, [150], [270.0], marker=:star5, markersize=14, color=:firebrick)
scatter!(ax3a, [150], [100.0], marker=:star5, markersize=14, color=:royalblue)
text!(ax3a, 95, 280;
      text="Lit. target: 270 g (pellet-fed)",
      fontsize=11, color=:firebrick)
text!(ax3a, 95, 108;
      text="Lit. target: 100 g (light-limited)",
      fontsize=11, color=:royalblue)
text!(ax3a, 5, 18;
      text="Stocking: $(Int(TIL_W0)) g, $(Int(TIL_DENS)) fish/m²",
      fontsize=10, color=:grey40)
axislegend(ax3a; position=:lt)

ax3b = Axis(fig3[1, 2],
    xlabel="Temperature (°C)",
    ylabel="Temperature scalar (0–1)",
    title="(b) Tilapia Thermal Performance")
Ts = range(10.0, 40.0, length=200)
tscalars = [til_tscalar(T) for T in Ts]
lines!(ax3b, collect(Ts), tscalars, linewidth=2.5, color=:darkorange)
vlines!(ax3b, [TIL_TOPT], color=:grey60, linestyle=:dash, linewidth=1)
text!(ax3b, TIL_TOPT + 0.5, 0.95;
      text="T_opt = $(Int(TIL_TOPT)) °C", fontsize=10, color=:grey40)
text!(ax3b, TIL_TBASE + 0.5, 0.05;
      text="T_base = $(Int(TIL_TBASE)) °C", fontsize=10, color=:grey40)

fw_final = round(r2.fish[end], digits=1)
fp_final = round(r3.fish[end], digits=1)
Label(fig3[2, 1:2],
    "Simulated final mass — Floodwater: $(fw_final) g (lit: ~100 g)  |  " *
    "Fishpond+pellet: $(fp_final) g (lit: ~270 g)",
    fontsize=13, halign=:center, padding=(5, 5, 5, 5))

save(joinpath(figdir, "tilapia_growth.png"), fig3, px_per_unit=2)
println("Saved tilapia_growth.png")

# ============================================================================
# Figure 4: Three scenarios — rice yield comparison
# ============================================================================
fig4 = Figure(size=(1000, 700))

ax4a = Axis(fig4[1, 1],
    xlabel="Days after transplanting",
    ylabel="Grain dry matter (g/hill)",
    title="(a) Grain Accumulation by Scenario")
lines!(ax4a, collect(days), r1.grain, linewidth=2.5, color=:forestgreen, label="Rice only")
lines!(ax4a, collect(days), r2.grain, linewidth=2, color=:royalblue, label="Rice + tilapia (flood)")
lines!(ax4a, collect(days), r3.grain, linewidth=2.5, color=:firebrick, label="Rice + fishpond (pellet)")
hlines!(ax4a, [14.0], color=:forestgreen, linestyle=:dot, linewidth=1)
hlines!(ax4a, [15.4], color=:firebrick, linestyle=:dot, linewidth=1)
text!(ax4a, 5, 14.5; text="Lit: 14 g/hill (rice only)", fontsize=10, color=:forestgreen)
text!(ax4a, 5, 15.9; text="Lit: 15.4 g/hill (pellet-fed)", fontsize=10, color=:firebrick)
axislegend(ax4a; position=:lt)

ax4b = Axis(fig4[1, 2],
    xlabel="Scenario",
    ylabel="Final grain (g/hill)",
    title="(b) Final Yield Comparison",
    xticks=(1:3, ["Rice\nonly", "Rice +\ntilapia\n(flood)", "Rice +\nfishpond\n(pellet)"]))
sim_yields = [r1.grain[end], r2.grain[end], r3.grain[end]]
barplot!(ax4b, [1, 2, 3], sim_yields,
         color=[:forestgreen, :royalblue, :firebrick],
         strokewidth=1, strokecolor=:black, label="Simulated")
scatter!(ax4b, [1, 3], [14.0, 15.4],
         marker=:star5, markersize=16, color=:gold,
         strokewidth=1, strokecolor=:black, label="Literature")
axislegend(ax4b; position=:lt)

ax4c = Axis(fig4[2, 1:2],
    xlabel="Scenario",
    ylabel="Fish final mass (g/fish)",
    title="(c) Tilapia Harvest Weight",
    xticks=(1:2, ["Floodwater\n(light-limited)", "Fishpond\n(pellet-fed)"]))
fish_sim = [r2.fish[end], r3.fish[end]]
fish_lit = [100.0, 270.0]
barplot!(ax4c, [1, 2], fish_sim,
         color=[:royalblue, :firebrick],
         strokewidth=1, strokecolor=:black, label="Simulated")
scatter!(ax4c, [1, 2], fish_lit,
         marker=:star5, markersize=16, color=:gold,
         strokewidth=1, strokecolor=:black, label="Lit. target")
axislegend(ax4c; position=:lt)

r1_yield = round(r1.grain[end], digits=2)
r3_yield = round(r3.grain[end], digits=2)
pct_gain = round((r3.grain[end] / max(r1.grain[end], 1e-10) - 1) * 100, digits=1)
Label(fig4[3, 1:2],
    "Rice yield: rice-only $(r1_yield) g/hill (lit: 14)  |  " *
    "Fishpond+pellet $(r3_yield) g/hill (lit: 15.4, +$(pct_gain)% sim.)",
    fontsize=12, halign=:center, padding=(5, 5, 5, 5))

save(joinpath(figdir, "scenario_comparison.png"), fig4, px_per_unit=2)
println("Saved scenario_comparison.png")

# ============================================================================
# Figure 5: Nutrient cycling — nitrogen pool dynamics across scenarios
# ============================================================================
# Tissue N concentrations (Yoshida 1981; Greenland 1997)
const N_LEAF  = 0.030  # leaf N fraction
const N_CULM  = 0.010
const N_ROOT  = 0.008
const N_GRAIN = 0.012
const N_PHYTO = 0.07   # Redfield ratio N:C
const N_ZOO   = 0.10
const N_FISH  = 0.025  # tilapia whole-body N

function compute_n_pools(r; has_fish=false)
    nd = length(r.grain) - 1
    n_rice  = r.leaf * N_LEAF .+ r.culm * N_CULM .+ r.root * N_ROOT .+ r.grain * N_GRAIN
    n_aqua  = r.phyto * N_PHYTO .+ r.zoo * N_ZOO
    n_fish  = has_fish ? r.fish * N_FISH * TIL_DENS : zeros(nd + 1)
    n_total = n_rice .+ n_aqua .+ n_fish
    (; rice=n_rice, aquatic=n_aqua, fish=n_fish, total=n_total)
end

n1 = compute_n_pools(r1; has_fish=false)
n2 = compute_n_pools(r2; has_fish=true)
n3 = compute_n_pools(r3; has_fish=true)

fig5 = Figure(size=(1000, 700))

ax5a = Axis(fig5[1, 1],
    xlabel="Days after transplanting",
    ylabel="N pool (g N / hill equiv.)",
    title="(a) Nitrogen in Rice Biomass")
lines!(ax5a, collect(days), n1.rice, linewidth=2.5, color=:forestgreen, label="Rice only")
lines!(ax5a, collect(days), n2.rice, linewidth=2, color=:royalblue, linestyle=:dash, label="+ Tilapia (flood)")
lines!(ax5a, collect(days), n3.rice, linewidth=2.5, color=:firebrick, linestyle=:dashdot, label="+ Fishpond (pellet)")
text!(ax5a, 100, maximum(n3.rice) * 0.85;
      text="Nutrient boost\nfrom fish feces",
      fontsize=10, color=:firebrick)
axislegend(ax5a; position=:lt)

ax5b = Axis(fig5[1, 2],
    xlabel="Days after transplanting",
    ylabel="N pool (g N / L equiv.)",
    title="(b) Nitrogen in Aquatic Compartments")
lines!(ax5b, collect(days), n1.aquatic, linewidth=2, color=:seagreen, label="Aquatic N (no fish)")
lines!(ax5b, collect(days), n2.aquatic, linewidth=2, color=:seagreen, linestyle=:dash, label="Aquatic N (flood fish)")
lines!(ax5b, collect(days), n2.fish, linewidth=2, color=:steelblue, label="Fish N (flood)")
lines!(ax5b, collect(days), n3.fish, linewidth=2, color=:firebrick, label="Fish N (pellet)")
axislegend(ax5b; position=:lt)

ax5c = Axis(fig5[2, 1:2],
    xlabel="Scenario",
    ylabel="Final N (g N)",
    title="(c) Nitrogen Budget at Harvest (day 150)",
    xticks=(1:3, ["Rice only", "Rice + tilapia\n(floodwater)", "Rice + fishpond\n(pellet-fed)"]))
rice_n_final = [n1.rice[end], n2.rice[end], n3.rice[end]]
aqua_n_final = [n1.aquatic[end], n2.aquatic[end], n3.aquatic[end]]
fish_n_final = [n1.fish[end], n2.fish[end], n3.fish[end]]

barplot!(ax5c, [1, 2, 3], rice_n_final,
         color=:forestgreen, strokewidth=1, strokecolor=:black,
         label="Rice N")
barplot!(ax5c, [1, 2, 3], rice_n_final .+ aqua_n_final,
         color=:seagreen, strokewidth=1, strokecolor=:black,
         label="+ Aquatic N")
barplot!(ax5c, [1, 2, 3], rice_n_final .+ aqua_n_final .+ fish_n_final,
         color=:steelblue, strokewidth=1, strokecolor=:black,
         label="+ Fish N")
axislegend(ax5c; position=:lt)

Label(fig5[3, 1:2],
    "N recycling: fish feces return ~10% of consumed N to paddy " *
    "(d'Oultremont & Gutierrez 2002).  " *
    "Area: 90% rice / 10% pond; fish 80% in paddy.",
    fontsize=11, halign=:center, padding=(5, 5, 5, 5))

save(joinpath(figdir, "nutrient_cycling.png"), fig5, px_per_unit=2)
println("Saved nutrient_cycling.png")

# ============================================================================
# Literature verification summary
# ============================================================================
println("\n" * "="^70)
println("LITERATURE VERIFICATION SUMMARY")
println("="^70)
println("Reference: d'Oultremont & Gutierrez (2002), Ecological Modelling")
println()

peak_til = maximum(r1.tiller)
peak_day = argmax(r1.tiller) - 1
println("Tiller dynamics:")
println("  Peak tiller count: $(round(peak_til, digits=1)) at day $peak_day " *
        "(Lit: 19–23 at day 32–35)")
println("  Final tiller count: $(round(r1.tiller[end], digits=1)) (Lit: 8–10)")
println()

println("Rice yield (g/hill):")
println("  Rice only:           $(round(r1.grain[end], digits=2))  (Lit: ~14)")
println("  Rice + tilapia:      $(round(r2.grain[end], digits=2))")
println("  Rice + fishpond:     $(round(r3.grain[end], digits=2))  (Lit: ~15.4)")
println()

println("Tilapia harvest mass (g/fish):")
println("  Floodwater:          $(round(r2.fish[end], digits=1))  (Lit: ~100)")
println("  Fishpond + pellet:   $(round(r3.fish[end], digits=1))  (Lit: ~270)")
println()

println("Key ecological insight: rice canopy shading → phytoplankton collapse")
println("  → zooplankton decline → fish growth constraint in floodwater.")
println("  Fishpond + pellet feeding decouples fish from this limitation.")
println()
println("All rice-fish validation figures saved to: $figdir")
