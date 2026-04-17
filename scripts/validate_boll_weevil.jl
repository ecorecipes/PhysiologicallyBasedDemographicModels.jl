#!/usr/bin/env julia
# Validation script for the cotton–boll weevil PBDM vignette.
# Generates 5 PNG figures comparing model output to published literature data.
#
# Subcycles the explicit Euler delay stepping to satisfy CFL stability:
#   Δdd per substep < τ/k  for every stage.

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using PhysiologicallyBasedDemographicModels
using CairoMakie

figdir = joinpath(@__DIR__, "figures", "boll_weevil")
mkpath(figdir)

# ══════════════════════════════════════════════════════════════════════
# Parameters from vignette (Gutierrez et al. 1991a,b)
# ══════════════════════════════════════════════════════════════════════

const COTTON_T_BASE  = 12.0
const COTTON_T_UPPER = 40.0
cotton_dev = LinearDevelopmentRate(COTTON_T_BASE, COTTON_T_UPPER)

const DD_FIRST_SQUARE = 400.0
const DD_CUTOUT       = 1100.0
const LEAF_GROWTH_RATE = 0.012

const BW_T_BASE  = 10.8
const BW_T_UPPER = 35.0
bw_dev = LinearDevelopmentRate(BW_T_BASE, BW_T_UPPER)

const BW_FECUNDITY         = 4.0
const BW_FEMALE_FRAC       = 0.50
const BW_OVIPOSITION_SCALE = 0.8
weevil_attack = FraserGilbertResponse(0.5)

const DAMAGE_COEFF   = 0.05
const POTENTIAL_YIELD = 1800.0
bw_damage = LinearDamageFunction(DAMAGE_COEFF)

# Degree-day durations per stage
const EGG_DD   = 57.0
const LARVA_DD = 110.0
const PUPA_DD  = 80.0
const ADULT_DD = 250.0
const TOTAL_IMMATURE_DD = EGG_DD + LARVA_DD + PUPA_DD  # 247

# ── Literature reference data ──
const LIT_BASE_TEMPS = Dict(:egg => 10.9, :larva => 6.6,
                             :pupa => 7.0, :total => 9.0)
const LIT_DD_LO = 248.0;  const LIT_DD_HI = 282.0
const LIT_DEV25_LO = 23.0;  const LIT_DEV25_HI = 24.0
const LIT_T_OPT_LO = 20.0;  const LIT_T_OPT_HI = 30.0
const LIT_UPPER_THERMAL = 39.0

# ══════════════════════════════════════════════════════════════════════
# Subcycled delay stepping (CFL-safe)
# ══════════════════════════════════════════════════════════════════════

"""Compute the number of subcycles needed so Δdd < min(τ/k) for stability."""
function n_subcycles(dd::Real, delays::Vector{<:DistributedDelay})
    min_ratio = minimum(d.τ / d.k for d in delays)
    return max(1, ceil(Int, dd / (0.9 * min_ratio)))  # 0.9 safety margin
end

"""Step a single delay through `dd` degree-days with automatic subcycling."""
function safe_step!(delay::DistributedDelay, dd::Real, inflow::Real;
                    μ::Real=0.0)
    nsub = max(1, ceil(Int, dd / (0.9 * delay.τ / delay.k)))
    dd_sub = dd / nsub
    inflow_sub = inflow / nsub
    total_out = 0.0;  total_atr = 0.0
    for _ in 1:nsub
        r = step_delay!(delay, dd_sub, inflow_sub; μ=μ)
        total_out += r.outflow
        total_atr += r.attrition
    end
    return (outflow=total_out, attrition=total_atr)
end

"""Step a sequential lifecycle (egg→larva→pupa→adult) with subcycling."""
function safe_step_lifecycle!(pop::Population, T_mean::Real)
    ns = n_stages(pop)
    inflow = 0.0
    for j in 1:ns
        s = pop.stages[j]
        dd = degree_days(s.dev_rate, T_mean)
        r = safe_step!(s.delay, dd, inflow; μ=s.μ)
        inflow = r.outflow  # maturation → next stage
    end
    return inflow  # final outflow (adult senescence)
end

# ══════════════════════════════════════════════════════════════════════
# 1. Figure 1 — Development rate vs temperature
# ══════════════════════════════════════════════════════════════════════

println("── Figure 1: Development rate vs temperature ──")

tv = collect(5.0:0.25:42.0)

# Development rate = DD/day / total_DD_required
dev_rate_total  = [degree_days(bw_dev, T) / TOTAL_IMMATURE_DD for T in tv]
dev_rate_egg    = [degree_days(bw_dev, T) / EGG_DD for T in tv]
dev_rate_larva  = [degree_days(bw_dev, T) / LARVA_DD for T in tv]
dev_rate_pupa   = [degree_days(bw_dev, T) / PUPA_DD for T in tv]

fig1 = Figure(size=(900, 600))
ax1 = Axis(fig1[1, 1],
    xlabel="Temperature (°C)",
    ylabel="Development rate (1/day)",
    title="Boll Weevil Development Rate vs Temperature\n(model: linear DD >$(BW_T_BASE)°C; literature: Greenberg 2005)")

lines!(ax1, tv, dev_rate_total, color=:black, linewidth=2.5,
       label="Total immature ($(Int(TOTAL_IMMATURE_DD)) DD)")
lines!(ax1, tv, dev_rate_egg,   color=:gold,   linewidth=1.5, linestyle=:dash,
       label="Egg ($(Int(EGG_DD)) DD)")
lines!(ax1, tv, dev_rate_larva, color=:orange,  linewidth=1.5, linestyle=:dash,
       label="Larva ($(Int(LARVA_DD)) DD)")
lines!(ax1, tv, dev_rate_pupa,  color=:coral,   linewidth=1.5, linestyle=:dash,
       label="Pupa ($(Int(PUPA_DD)) DD)")

vlines!(ax1, [BW_T_BASE], color=:blue, linestyle=:dot, linewidth=1.5,
        label="Model T_base = $(BW_T_BASE)°C")
vlines!(ax1, [LIT_BASE_TEMPS[:total]], color=:red, linestyle=:dot, linewidth=1.5,
        label="Lit. T_base = $(LIT_BASE_TEMPS[:total])°C (Greenberg)")
vspan!(ax1, LIT_T_OPT_LO, LIT_T_OPT_HI, color=(:green, 0.08))
text!(ax1, 25.0, maximum(dev_rate_total) * 0.95,
      text="Optimal 20–30°C", align=(:center, :top), fontsize=11, color=:green4)

# Literature overlay: dev rate at 25°C → 1/23 to 1/24 per day
scatter!(ax1, [25.0, 25.0], [1.0 / LIT_DEV25_HI, 1.0 / LIT_DEV25_LO],
         color=:red, markersize=12, marker=:diamond,
         label="Lit. 25°C: 23–24 d (MDPI 2023)")

model_rate_25 = degree_days(bw_dev, 25.0) / TOTAL_IMMATURE_DD
println("  Model dev rate at 25°C: $(round(model_rate_25, digits=4)) /day → $(round(1/model_rate_25, digits=1)) days")
println("  Literature: $(LIT_DEV25_LO)–$(LIT_DEV25_HI) days")
println("  Model total immature DD: $(Int(TOTAL_IMMATURE_DD))  Lit: $(Int(LIT_DD_LO))–$(Int(LIT_DD_HI))")

axislegend(ax1, position=:lt, framevisible=true, labelsize=10)
save(joinpath(figdir, "devrate_vs_temperature.png"), fig1, px_per_unit=2)
println("  → Saved devrate_vs_temperature.png")

# ══════════════════════════════════════════════════════════════════════
# 2. Figure 2 — Mortality vs temperature
# ══════════════════════════════════════════════════════════════════════

println("\n── Figure 2: Mortality vs temperature ──")

function bw_mortality(T; μ_base=0.003, T_lo=20.0, T_hi=30.0,
                      T_cold=5.0, T_heat=39.0)
    if T < T_lo
        return μ_base + 0.5 * max(0, (T_lo - T) / (T_lo - T_cold))^2
    elseif T > T_hi
        return μ_base + 0.8 * max(0, (T - T_hi) / (T_heat - T_hi))^2
    else
        return μ_base
    end
end

egg_μ(T)   = bw_mortality(T; μ_base=0.005, T_lo=18.0, T_hi=32.0, T_cold=5.0, T_heat=39.0)
larva_μ(T) = bw_mortality(T; μ_base=0.003, T_lo=18.0, T_hi=30.0, T_cold=5.0, T_heat=38.0)
pupa_μ(T)  = bw_mortality(T; μ_base=0.002, T_lo=18.0, T_hi=30.0, T_cold=5.0, T_heat=38.0)
adult_μ(T) = bw_mortality(T; μ_base=0.004, T_lo=18.0, T_hi=32.0, T_cold=3.0, T_heat=39.0)

mtv = collect(5.0:0.25:42.0)

fig2 = Figure(size=(900, 600))
ax2 = Axis(fig2[1, 1],
    xlabel="Temperature (°C)", ylabel="Mortality rate (per DD)",
    title="Boll Weevil Temperature-Dependent Mortality\n(upper thermal limit ≈39°C)")

lines!(ax2, mtv, egg_μ.(mtv),   color=:gold,    linewidth=2, label="Egg (μ₀=0.005)")
lines!(ax2, mtv, larva_μ.(mtv), color=:orange,  linewidth=2, label="Larva (μ₀=0.003)")
lines!(ax2, mtv, pupa_μ.(mtv),  color=:coral,   linewidth=2, label="Pupa (μ₀=0.002)")
lines!(ax2, mtv, adult_μ.(mtv), color=:darkred,  linewidth=2, label="Adult (μ₀=0.004)")

vlines!(ax2, [LIT_UPPER_THERMAL], color=:red, linestyle=:dash, linewidth=2,
        label="Upper thermal limit ≈39°C")
vspan!(ax2, LIT_T_OPT_LO, LIT_T_OPT_HI, color=(:green, 0.08))
text!(ax2, 25.0, 0.005, text="Optimal 20–30°C",
      align=(:center, :bottom), fontsize=11, color=:green4)
hlines!(ax2, [0.005, 0.003, 0.002, 0.004], color=:gray70, linestyle=:dot, linewidth=0.8)

axislegend(ax2, position=:lt, framevisible=true, labelsize=10)
save(joinpath(figdir, "mortality_vs_temperature.png"), fig2, px_per_unit=2)
println("  → Saved mortality_vs_temperature.png")

# ══════════════════════════════════════════════════════════════════════
# 3. Figure 3 — Constant 25°C weevil simulation (180 days)
# ══════════════════════════════════════════════════════════════════════

println("\n── Figure 3: Constant 25°C simulation (180 days) ──")

n_sim = 180
const_T = 25.0
const_weather = WeatherSeries([DailyWeather(const_T) for _ in 1:n_sim])

function build_weevil(; N0=5.0)
    Population(:boll_weevil, [
        LifeStage(:egg,   DistributedDelay(12, EGG_DD;   W0=0.0), bw_dev, 0.005),
        LifeStage(:larva, DistributedDelay(20, LARVA_DD; W0=0.0), bw_dev, 0.003),
        LifeStage(:pupa,  DistributedDelay(15, PUPA_DD;  W0=0.0), bw_dev, 0.002),
        LifeStage(:adult, DistributedDelay(15, ADULT_DD; W0=N0),  bw_dev, 0.004),
    ])
end

function simulate_weevil_const(T_mean, n_days; N0=5.0, K=2000.0)
    pop = build_weevil(N0=N0)
    ns = 4
    totals = zeros(n_days + 1, ns)
    cdd = zeros(n_days + 1)
    for j in 1:ns; totals[1, j] = delay_total(pop.stages[j].delay); end
    cum = 0.0

    for d in 1:n_days
        dd = degree_days(bw_dev, T_mean)
        cum += dd

        # Reproduction (density-limited)
        adults = delay_total(pop.stages[4].delay)
        total_pop = sum(delay_total(pop.stages[j].delay) for j in 1:ns)
        dd_frac = dd / ADULT_DD
        eggs = BW_FECUNDITY * BW_FEMALE_FRAC * adults * dd_frac *
               max(0.0, 1.0 - total_pop / K)
        pop.stages[1].delay.W[1] += max(0.0, eggs)

        # Subcycled lifecycle step
        safe_step_lifecycle!(pop, T_mean)

        for j in 1:ns; totals[d+1, j] = delay_total(pop.stages[j].delay); end
        cdd[d+1] = cum
    end
    return (; totals, cdd)
end

res25 = simulate_weevil_const(const_T, n_sim)

dd_day25 = degree_days(bw_dev, const_T)
dev_time_25 = TOTAL_IMMATURE_DD / dd_day25
println("  DD/day at 25°C: $(dd_day25)")
println("  Model immature dev time: $(round(dev_time_25, digits=1)) days")
println("  Lit: $(LIT_DEV25_LO)–$(LIT_DEV25_HI) days")
println("  Peak adults: $(round(maximum(res25.totals[:, 4]), digits=1))")
println("  Final total: $(round(sum(res25.totals[end, :]), digits=1))")

snames = ["Eggs", "Larvae", "Pupae", "Adults"]
scols  = [:gold, :orange, :coral, :darkred]
days_v = 0:n_sim

fig3 = Figure(size=(900, 700))
ax3a = Axis(fig3[1, 1], xlabel="Days", ylabel="Population",
    title="Boll Weevil at Constant 25°C\n(initial 5 adults; generation ≈ $(round(dev_time_25, digits=0)) days)")

for (j, (nm, cl)) in enumerate(zip(snames, scols))
    lines!(ax3a, days_v, res25.totals[:, j], label=nm, color=cl, linewidth=2)
end
vlines!(ax3a, [dev_time_25], color=:blue, linestyle=:dash, linewidth=1.5,
        label="Model gen. $(round(dev_time_25, digits=0))d")
vlines!(ax3a, [LIT_DEV25_LO], color=:red, linestyle=:dot, linewidth=1.5,
        label="Lit. 23–24 d")
axislegend(ax3a, position=:rt, framevisible=true, labelsize=10)

ax3b = Axis(fig3[2, 1], xlabel="Days",
    ylabel="Cumulative DD (>$(BW_T_BASE)°C)", title="Cumulative Degree-Days")
lines!(ax3b, days_v, res25.cdd, color=:black, linewidth=2, label="CDD at 25°C")
hlines!(ax3b, [TOTAL_IMMATURE_DD], color=:blue, linestyle=:dash, linewidth=1.5,
        label="Total immature = $(Int(TOTAL_IMMATURE_DD)) DD")
hlines!(ax3b, [LIT_DD_LO, LIT_DD_HI], color=:red, linestyle=:dot, linewidth=1.5,
        label="Lit. $(Int(LIT_DD_LO))–$(Int(LIT_DD_HI)) DD")
axislegend(ax3b, position=:lt, framevisible=true, labelsize=10)

save(joinpath(figdir, "sim_constant_25C.png"), fig3, px_per_unit=2)
println("  → Saved sim_constant_25C.png")

# ══════════════════════════════════════════════════════════════════════
# 4. Figure 4 — Cotton with vs without boll weevil
# ══════════════════════════════════════════════════════════════════════

println("\n── Figure 4: Cotton with vs without boll weevil ──")

# Synthetic Londrina weather (from vignette)
n_season = 210
weather_days = DailyWeather{Float64}[]
for d in 1:n_season
    doy = mod(287 + d, 365) + 1
    T_m = 23.0 + 4.5 * sin(2π * (doy - 355) / 365)
    push!(weather_days, DailyWeather(T_m, max(T_m - 4.5, 12.0), min(T_m + 5.0, 36.0);
                                      radiation=18.0 + 5.5*sin(2π*(doy-355)/365),
                                      photoperiod=12.0 + 1.2*sin(2π*(doy-355)/365)))
end
weather = WeatherSeries(weather_days; day_offset=1)

# Cotton organs are PARALLEL pools, not sequential — step each independently.
const LEAF_TAU  = 700.0;  const LEAF_K  = 25
const STEM_TAU  = 2000.0; const STEM_K  = 25
const ROOT_TAU  = 150.0;  const ROOT_K  = 25
const FRUIT_TAU = 800.0;  const FRUIT_K = 25

# Max biomass per organ (g per plant, realistic for IAC-17 cotton)
const MAX_LEAF  = 200.0
const MAX_STEM  = 300.0
const MAX_ROOT  = 100.0
const MAX_FRUIT = 500.0

function simulate_system(weather, n_days; with_weevil=true, colonizer_rate=2.0)
    # Cotton organs (independent delays)
    leaf  = DistributedDelay(LEAF_K,  LEAF_TAU;  W0=0.15)
    stem  = DistributedDelay(STEM_K,  STEM_TAU;  W0=0.10)
    root  = DistributedDelay(ROOT_K,  ROOT_TAU;  W0=0.08)
    fruit = DistributedDelay(FRUIT_K, FRUIT_TAU; W0=0.0)

    # Weevil lifecycle (sequential: egg → larva → pupa → adult)
    bw = build_weevil(N0=0.0)

    cot_out = zeros(n_days + 1, 4)  # leaf, stem, root, fruit (current biomass)
    bw_out  = zeros(n_days + 1, 4)
    cdd_cot = zeros(n_days + 1)
    cum_yield = zeros(n_days + 1)   # cumulative matured boll biomass
    sq_dmg  = zeros(n_days)

    cot_out[1, :] = [delay_total(leaf), delay_total(stem),
                     delay_total(root), delay_total(fruit)]

    cum_cot = 0.0;  cum_bw = 0.0

    for d in 1:n_days
        w = get_weather(weather, d)
        dd_c = degree_days(cotton_dev, w.T_mean)
        dd_b = degree_days(bw_dev, w.T_mean)
        cum_cot += dd_c;  cum_bw += dd_b

        leaf_mass  = delay_total(leaf)
        fruit_mass = delay_total(fruit)

        # Logistic new-leaf growth (photosynthesis → leaves, capped)
        leaf_growth = 0.012 * leaf_mass * dd_c * max(0.0, 1.0 - leaf_mass / MAX_LEAF)

        # Stem & root allocation proportional to leaf growth
        stem_growth = 0.3 * leaf_growth * max(0.0, 1.0 - delay_total(stem) / MAX_STEM)
        root_growth = 0.2 * leaf_growth * max(0.0, 1.0 - delay_total(root) / MAX_ROOT)

        # Fruit initiation once squaring threshold reached
        fruit_growth = 0.0
        if cum_cot >= DD_FIRST_SQUARE && cum_cot < DD_CUTOUT
            fruit_growth = 0.020 * leaf_mass * dd_c * 0.3 *
                           max(0.0, 1.0 - fruit_mass / MAX_FRUIT)
        end

        # ── Weevil damage ──
        if with_weevil
            if cum_cot >= DD_FIRST_SQUARE && cum_cot < DD_FIRST_SQUARE + 100.0
                bw.stages[4].delay.W[1] += colonizer_rate
            end

            adults = delay_total(bw.stages[4].delay)
            if adults > 0.01 && fruit_mass > 0.01
                dd_frac = dd_b / ADULT_DD
                ovi_demand = BW_FECUNDITY * BW_FEMALE_FRAC *
                             BW_OVIPOSITION_SCALE * adults * dd_frac
                eggs_laid = acquire(weevil_attack, fruit_mass, ovi_demand)
                bw.stages[1].delay.W[1] += eggs_laid

                dmg_frac = min(eggs_laid / max(fruit_mass, 1e-10), 0.8)
                for i in 1:fruit.k
                    fruit.W[i] *= (1.0 - dmg_frac)
                end
                sq_dmg[d] = dmg_frac
            end

            total_bw = sum(delay_total(bw.stages[j].delay) for j in 1:4)
            repro = BW_FECUNDITY * BW_FEMALE_FRAC * max(0.0, adults) *
                    dd_b / ADULT_DD * max(0.0, 1.0 - total_bw / 5000.0)
            bw.stages[1].delay.W[1] += max(0.0, repro)

            safe_step_lifecycle!(bw, w.T_mean)
        end

        # Step cotton organs; fruit outflow = matured bolls (= yield)
        safe_step!(leaf,  dd_c, leaf_growth;  μ=0.0008)
        safe_step!(stem,  dd_c, stem_growth;  μ=0.0004)
        safe_step!(root,  dd_c, root_growth;  μ=0.0015)
        fruit_result = safe_step!(fruit, dd_c, fruit_growth; μ=0.0010)

        # Accumulate matured boll biomass as harvestable yield
        cum_yield[d+1] = cum_yield[d] + fruit_result.outflow

        cot_out[d+1, :] = [delay_total(leaf), delay_total(stem),
                           delay_total(root), delay_total(fruit)]
        for j in 1:4; bw_out[d+1, j] = delay_total(bw.stages[j].delay); end
        cdd_cot[d+1] = cum_cot
    end
    return (; cot_out, bw_out, cdd_cot, cum_yield, sq_dmg)
end

res_clean = simulate_system(weather, n_season; with_weevil=false)
res_inf   = simulate_system(weather, n_season; with_weevil=true)

clean_fruit = res_clean.cum_yield[end]
inf_fruit   = res_inf.cum_yield[end]
fruit_loss  = 100.0 * (1.0 - inf_fruit / max(clean_fruit, 1e-10))
println("  Clean cumulative yield: $(round(clean_fruit, digits=2)) g")
println("  Infested cumulative yield: $(round(inf_fruit, digits=2)) g")
println("  Yield reduction: $(round(fruit_loss, digits=1))%")

sd = 0:n_season
fig4 = Figure(size=(950, 850))

ax4a = Axis(fig4[1, 1], xlabel="Days after planting",
    ylabel="Cumulative boll yield (g)",
    title="(a) Cotton Boll Yield: Clean vs Infested")
lines!(ax4a, sd, res_clean.cum_yield, color=:forestgreen, linewidth=2.5,
       label="No weevil")
lines!(ax4a, sd, res_inf.cum_yield, color=:red, linewidth=2.5,
       label="With boll weevil")
sq_day = findfirst(>=(DD_FIRST_SQUARE), res_clean.cdd_cot)
if sq_day !== nothing
    vlines!(ax4a, [sq_day - 1], color=:gray, linestyle=:dash, linewidth=1,
            label="Squaring onset (~400 DD)")
end
axislegend(ax4a, position=:lt)

ax4b = Axis(fig4[2, 1], xlabel="Days after planting",
    ylabel="Biomass (g)", title="(b) Cotton Organ Dynamics (no pest)")
lines!(ax4b, sd, res_clean.cot_out[:, 1], color=:green,  linewidth=2, label="Leaf")
lines!(ax4b, sd, res_clean.cot_out[:, 2], color=:brown,  linewidth=2, label="Stem")
lines!(ax4b, sd, res_clean.cot_out[:, 3], color=:orange, linewidth=2, label="Root")
lines!(ax4b, sd, res_clean.cot_out[:, 4], color=:purple, linewidth=2, label="Fruit")
axislegend(ax4b, position=:lt)

ax4c = Axis(fig4[3, 1], xlabel="Days after planting",
    ylabel="Population", title="(c) Boll Weevil Stage Dynamics (infested)")
lines!(ax4c, sd, res_inf.bw_out[:, 1], color=:gold,    linewidth=2, label="Eggs")
lines!(ax4c, sd, res_inf.bw_out[:, 2], color=:orange,  linewidth=2, label="Larvae")
lines!(ax4c, sd, res_inf.bw_out[:, 3], color=:coral,   linewidth=2, label="Pupae")
lines!(ax4c, sd, res_inf.bw_out[:, 4], color=:darkred,  linewidth=2, label="Adults")
axislegend(ax4c, position=:rt)

save(joinpath(figdir, "cotton_with_vs_without_weevil.png"), fig4, px_per_unit=2)
println("  → Saved cotton_with_vs_without_weevil.png")

# ══════════════════════════════════════════════════════════════════════
# 5. Figure 5 — Economic threshold
# ══════════════════════════════════════════════════════════════════════

println("\n── Figure 5: Economic threshold ──")

price_per_kg    = 1.80
production_cost = 1200.0
breakeven_yield = production_cost / price_per_kg

init_levels = [0.0, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0]
yields  = Float64[]
profits = Float64[]

clean_baseline = simulate_system(weather, n_season; with_weevil=false).cum_yield[end]

for n_init in init_levels
    r = n_init == 0.0 ? 0.0 : n_init / 5.0
    res = simulate_system(weather, n_season; with_weevil=(n_init > 0),
                          colonizer_rate=r)
    # Scale model yield to POTENTIAL_YIELD: model clean → potential
    y = POTENTIAL_YIELD * res.cum_yield[end] / max(clean_baseline, 1e-10)
    push!(yields, y)
    push!(profits, net_profit(y * price_per_kg, production_cost))
end

fig5 = Figure(size=(900, 700))

ax5a = Axis(fig5[1, 1], xlabel="Initial weevil adults",
    ylabel="Lint yield (kg/ha)",
    title="(a) Yield Response to Boll Weevil Pressure")
scatterlines!(ax5a, init_levels, yields, linewidth=2.5, markersize=10,
              color=:firebrick, label="Projected yield")
hlines!(ax5a, [breakeven_yield], color=:black, linestyle=:dash, linewidth=1.5,
        label="Break-even ($(round(Int, breakeven_yield)) kg/ha)")
hlines!(ax5a, [POTENTIAL_YIELD], color=:forestgreen, linestyle=:dot, linewidth=1,
        label="Potential ($(round(Int, POTENTIAL_YIELD)) kg/ha)")
axislegend(ax5a, position=:rt, labelsize=10)

ax5b = Axis(fig5[2, 1], xlabel="Initial weevil adults",
    ylabel="Net profit (USD/ha)",
    title="(b) Profit vs Infestation Pressure\n(cotton \$$(price_per_kg)/kg; cost \$$(round(Int, production_cost))/ha)")
scatterlines!(ax5b, init_levels, profits, linewidth=2.5, markersize=10,
              color=:navy, label="Net profit")
hlines!(ax5b, [0.0], color=:red, linestyle=:dash, linewidth=1.5,
        label="Break-even (profit = 0)")
axislegend(ax5b, position=:rt, labelsize=10)
text!(ax5b, init_levels[end] * 0.55, minimum(profits) * 0.3,
      text="Literature: 30–70% yield loss\nunder heavy infestation\n(Gutierrez et al. 1991)",
      align=(:center, :center), fontsize=10, color=:gray40)

save(joinpath(figdir, "economic_threshold.png"), fig5, px_per_unit=2)
println("  → Saved economic_threshold.png")

# ══════════════════════════════════════════════════════════════════════
# Summary
# ══════════════════════════════════════════════════════════════════════

println("\n══ Validation Summary ══")
println("Model T_base: $(BW_T_BASE)°C  (Lit: egg=$(LIT_BASE_TEMPS[:egg]), larva=$(LIT_BASE_TEMPS[:larva]), pupa=$(LIT_BASE_TEMPS[:pupa]), total=$(LIT_BASE_TEMPS[:total])°C)")
println("Model immature DD: $(Int(TOTAL_IMMATURE_DD))  (Lit: $(Int(LIT_DD_LO))–$(Int(LIT_DD_HI)))")
println("Model dev time 25°C: $(round(TOTAL_IMMATURE_DD / degree_days(bw_dev, 25.0), digits=1)) d  (Lit: $(LIT_DEV25_LO)–$(LIT_DEV25_HI) d)")
println("Fruit loss (standard infestation): $(round(fruit_loss, digits=1))%")
println()

println("══ All figures saved to: $(figdir) ══")
for f in sort(readdir(figdir))
    sz = round(filesize(joinpath(figdir, f)) / 1024, digits=1)
    println("  $f  ($(sz) KB)")
end
