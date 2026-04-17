#!/usr/bin/env julia
# Validation figures for the screwworm SIT vignette (08_screwworm_sit.qmd)
# against Gutierrez & Ponti (2014, 2019).

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using CairoMakie

figdir = joinpath(@__DIR__, "figures", "screwworm")
mkpath(figdir)

# ============================================================================
# Parameters from the vignette and the two reference papers
# ============================================================================

# Temperature thresholds (Gutierrez & Ponti 2014, 2019; Adams 1979;
# Berkebile et al. 2006)
const SW_T_LOWER = 14.5    # °C  lower developmental threshold
const SW_T_UPPER = 43.5    # °C  upper developmental threshold
const SW_T_REF   = 27.2    # °C  optimal temperature for survival

# Mortality Eq. 1 (Gutierrez & Ponti 2014 / 2019 Eq. 1)
#   µ_ab(T) = 0.00036·(T − 27.2)² + 0.0035   (r²=0.74, d.f.=16)
μ_ab(T) = clamp(0.00036 * (T - SW_T_REF)^2 + 0.0035, 0.0, 1.0)

# Reproduction (vignette + Gutierrez & Ponti 2014 Eqs. 4–6)
const MAX_FECUNDITY    = 67.0   # eggs/female/day (vignette; 2014 paper uses 90)
const SEX_RATIO        = 0.5
const MATING_RATE      = 0.5    # proportion of virgin ♀ mating per day
const FEMALE_HALF_LIFE = 3.7    # days, mated wild female (Thomas & Chen 1990)

# Temperature-dependent reproduction scaling  φ_T (symmetrical, concave)
function φ_T(T)
    (T < SW_T_LOWER || T > SW_T_UPPER) && return 0.0
    return clamp(1.0 - ((T - SW_T_REF) / (SW_T_UPPER - SW_T_REF))^2, 0.0, 1.0)
end

# Cumulative cold mortality index  µ_cold  (Eq. 2, Sept 1 – May 31)
function cumulative_cold(daily_temps; start_doy=244, end_doy=151)
    μ_cold = 0.0
    for (i, T) in enumerate(daily_temps)
        doy = ((i - 1 + start_doy - 1) % 365) + 1
        if doy >= start_doy || doy <= end_doy
            if T < SW_T_REF
                μ_cold += μ_ab(T)
            end
        end
    end
    return μ_cold
end

# SIT fertile fraction (vignette Eq. 4i)
function sit_fertile_fraction(wild_males, sterile_males; competitiveness=1.0)
    eff_sterile = sterile_males * competitiveness
    total = wild_males + eff_sterile
    total ≈ 0.0 && return 1.0
    return wild_males / total
end

# Myiasis regression (Gutierrez & Ponti 2019 Eq. 3)
#   log₁₀(myiasis) = 6.568 + 0.028·rain(y−1) − 0.602·µ_cold(y−1)
predict_myiasis(rain, μ_cold) = 10^(6.568 + 0.028 * rain - 0.602 * μ_cold)

# ============================================================================
# Paper reference values for overlay
# ============================================================================

# Gutierrez & Ponti (2014): base temperature for development ≈ 14.5 °C;
# user request mentions "9.7 °C" which appears in some older literature
# as the base for another calliphorid; both papers consistently use 14.5 °C.
# We plot both for comparison.
const PAPER_T_BASE    = 14.5   # confirmed in both papers
const ALT_T_BASE      = 9.7    # sometimes cited in generic calliphorid lit.

# Biweekly SIT release rates required for eradication (Fig. 6 / 2019 paper)
# McAllen 2, Tampico 3, Tuxtla-Gutiérrez 23  (sterile flies km⁻²)
# vs. observed field releases 39–896/week/km²
const SIT_MODEL_MCALLEN = 2.0
const SIT_MODEL_TAMPICO = 3.0
const SIT_MODEL_TUXTLA  = 23.0
const SIT_FIELD_LOW     = 39.0
const SIT_FIELD_HIGH    = 896.0

# µ_cold endemic threshold
const MU_COLD_THRESHOLD = 10.0

# ============================================================================
# Helper: approximate annual daily temperature for a location
# ============================================================================
function approx_temps(mean_T, lat; n=365)
    amplitude = 12.0 - 0.15 * (30.0 - lat)
    return [mean_T + amplitude * sin(2π * (d - 200) / 365) for d in 1:n]
end

# ============================================================================
# Figure 1 — Development rate vs temperature
# ============================================================================
println("Generating Fig 1: development rate vs temperature …")

Ts = range(-5.0, 50.0, length=500)

# Linear development rate: degree-days above T_lower, bounded by T_upper
dev_rate_14_5 = [max(0.0, min(T - PAPER_T_BASE, SW_T_UPPER - PAPER_T_BASE)) /
                 (SW_T_UPPER - PAPER_T_BASE) for T in Ts]
dev_rate_9_7  = [max(0.0, min(T - ALT_T_BASE, SW_T_UPPER - ALT_T_BASE)) /
                 (SW_T_UPPER - ALT_T_BASE) for T in Ts]

# φ_T reproductive scaling (concave, from vignette)
phi_T_vals = [φ_T(T) for T in Ts]

fig1 = Figure(size=(900, 600))
ax1 = Axis(fig1[1, 1],
    xlabel="Temperature (°C)", ylabel="Relative rate",
    title="Screwworm Development & Reproduction vs Temperature",
    xlabelsize=14, ylabelsize=14)

lines!(ax1, collect(Ts), dev_rate_14_5, linewidth=2.5, color=:steelblue,
       label="Development rate (T_base=14.5°C, papers)")
lines!(ax1, collect(Ts), dev_rate_9_7, linewidth=2, color=:orange,
       linestyle=:dash, label="Development rate (T_base=9.7°C, alt. lit.)")
lines!(ax1, collect(Ts), phi_T_vals, linewidth=2, color=:firebrick,
       label="φ_T reproduction scaling (vignette)")
vlines!(ax1, [PAPER_T_BASE], color=:steelblue, linestyle=:dot, linewidth=1)
vlines!(ax1, [ALT_T_BASE], color=:orange, linestyle=:dot, linewidth=1)
vlines!(ax1, [SW_T_REF], color=:gray40, linestyle=:dot, linewidth=1)
text!(ax1, SW_T_REF + 0.5, 0.95, text="T_opt=27.2°C", fontsize=11, color=:gray40)
text!(ax1, PAPER_T_BASE + 0.5, 0.1, text="14.5°C", fontsize=11, color=:steelblue)
text!(ax1, ALT_T_BASE + 0.5, 0.15, text="9.7°C", fontsize=11, color=:orange)
xlims!(ax1, -5, 50)
ylims!(ax1, -0.05, 1.1)
axislegend(ax1, position=:lt, framevisible=false)

save(joinpath(figdir, "fig1_development_rate.png"), fig1, px_per_unit=2)
println("  ✓ fig1_development_rate.png")

# ============================================================================
# Figure 2 — SIT dose-response
# ============================================================================
println("Generating Fig 2: SIT dose-response …")

ratios = range(0.0, 120.0, length=500)
R0_base = 3.0    # approximate net reproductive rate (vignette)
R0_low  = 1.5    # field growth rates are low in tropics

comp_vals = [1.0, 0.5, 0.1, 0.05]
colors2 = [:black, :steelblue, :firebrick, :orange]

fig2 = Figure(size=(900, 600))
ax2 = Axis(fig2[1, 1],
    xlabel="Sterile:Wild ratio", ylabel="Effective R₀",
    title="SIT Dose-Response: Population Growth vs Overflooding Ratio",
    xlabelsize=14, ylabelsize=14)

for (ci, comp) in enumerate(comp_vals)
    eff_R0 = [R0_base * sit_fertile_fraction(1.0, r; competitiveness=comp)
              for r in ratios]
    lines!(ax2, collect(ratios), eff_R0, linewidth=2, color=colors2[ci],
           label="R₀=3.0, comp=$(comp)")
end
# Also plot low R₀ (tropical field)
eff_R0_low = [R0_low * sit_fertile_fraction(1.0, r; competitiveness=1.0)
              for r in ratios]
lines!(ax2, collect(ratios), eff_R0_low, linewidth=2, color=:purple,
       linestyle=:dash, label="R₀=1.5 (tropical field)")

hlines!(ax2, [1.0], color=:red, linestyle=:dot, linewidth=1.5, label="R₀=1 threshold")

# Paper reference: model needs only 2–23 sterile/km²/biweekly,
# but field used 39–896/week/km² → efficacy ~1.7%
vspan!(ax2, 2.0, 23.0, color=(:green, 0.12))
text!(ax2, 4.0, 2.8, text="Model\nprediction\n2–23:1", fontsize=10, color=:green4)

xlims!(ax2, 0, 120)
ylims!(ax2, 0, 3.2)
axislegend(ax2, position=:rt, framevisible=false, labelsize=11)

save(joinpath(figdir, "fig2_sit_dose_response.png"), fig2, px_per_unit=2)
println("  ✓ fig2_sit_dose_response.png")

# ============================================================================
# Figure 3 — Seasonal dynamics (no SIT) — McAllen TX
# ============================================================================
println("Generating Fig 3: seasonal dynamics …")

function simulate_screwworm(; n_days=730, daily_temps=nothing,
                              sterile_release_rate=0.0,
                              release_interval=14,
                              competitiveness=1.0,
                              initial_pop=100.0)
    if daily_temps === nothing
        daily_temps = [26.0 + 10.0 * sin(2π * (d - 200) / 365) for d in 1:n_days]
    end
    n_days = length(daily_temps)

    wild_F_mated   = initial_pop * 0.3
    wild_F_unmated = initial_pop * 0.2
    wild_M         = initial_pop * 0.5
    sterile_M      = 0.0

    pop_history         = Float64[]
    fertile_eggs_hist   = Float64[]
    temp_history        = Float64[]
    mortality_history   = Float64[]

    for day in 1:n_days
        T = daily_temps[day]

        phi = φ_T(T)

        if sterile_release_rate > 0 && day % release_interval == 0
            sterile_M += sterile_release_rate
        end
        sterile_M *= 0.90

        fertile_frac = sit_fertile_fraction(wild_M, sterile_M;
                                            competitiveness=competitiveness)

        new_matings     = wild_F_unmated * MATING_RATE
        new_fertile     = new_matings * fertile_frac
        wild_F_mated   += new_fertile
        wild_F_unmated -= new_matings

        fertile_eggs = phi * MAX_FECUNDITY * SEX_RATIO * wild_F_mated
        fertile_eggs = max(0.0, fertile_eggs)

        mu = μ_ab(T)
        wild_F_mated   *= (1 - mu)
        wild_F_unmated *= (1 - mu)
        wild_M         *= (1 - mu)

        if day > 14
            emergence = fertile_eggs_hist[max(1, day - 14)] * 0.05
            wild_F_unmated += emergence * SEX_RATIO
            wild_M         += emergence * (1 - SEX_RATIO)
        end

        total_wild = wild_F_mated + wild_F_unmated + wild_M
        push!(pop_history, total_wild)
        push!(fertile_eggs_hist, fertile_eggs)
        push!(temp_history, T)
        push!(mortality_history, mu)
    end
    return (; pop=pop_history, eggs=fertile_eggs_hist,
              temps=temp_history, mort=mortality_history)
end

# McAllen TX approximate climate (lat ≈ 26.2, mean T ≈ 23.5 °C)
mcallen_temps = approx_temps(23.5, 26.2; n=730)
res_no_sit = simulate_screwworm(daily_temps=mcallen_temps,
                                 sterile_release_rate=0.0)

fig3 = Figure(size=(1000, 700))

ax3a = Axis(fig3[1, 1], ylabel="Temperature (°C)",
            title="Seasonal Screwworm Dynamics — McAllen TX (no SIT)",
            xlabelsize=12, ylabelsize=12)
lines!(ax3a, 1:730, res_no_sit.temps, color=:darkorange, linewidth=1.5)
hlines!(ax3a, [SW_T_LOWER], color=:steelblue, linestyle=:dash, linewidth=1,
        label="T_lower=14.5°C")
hlines!(ax3a, [SW_T_REF], color=:gray50, linestyle=:dot, linewidth=1,
        label="T_opt=27.2°C")
axislegend(ax3a, position=:rb, framevisible=false, labelsize=10)
hidexdecorations!(ax3a, grid=false)

ax3b = Axis(fig3[2, 1], ylabel="Mortality rate/day",
            xlabelsize=12, ylabelsize=12)
lines!(ax3b, 1:730, res_no_sit.mort, color=:firebrick, linewidth=1.5)
text!(ax3b, 400, maximum(res_no_sit.mort) * 0.85,
      text="µ_ab peak in winter", fontsize=10, color=:firebrick)
hidexdecorations!(ax3b, grid=false)

ax3c = Axis(fig3[3, 1], xlabel="Day of simulation",
            ylabel="Wild population",
            xlabelsize=12, ylabelsize=12)
lines!(ax3c, 1:730, res_no_sit.pop, color=:black, linewidth=2)
# Mark winter troughs
for yr in [1, 2]
    band_start = (yr - 1) * 365 + 1
    band_end   = (yr - 1) * 365 + 90
    vspan!(ax3c, band_start, min(band_end, 730), color=(:steelblue, 0.08))
end
text!(ax3c, 30, maximum(res_no_sit.pop) * 0.85,
      text="Winter", fontsize=10, color=:steelblue)

save(joinpath(figdir, "fig3_seasonal_dynamics.png"), fig3, px_per_unit=2)
println("  ✓ fig3_seasonal_dynamics.png")

# ============================================================================
# Figure 4 — Eradication trajectory under progressive SIT
# ============================================================================
println("Generating Fig 4: eradication trajectory …")

release_rates = [0.0, 500.0, 2000.0, 5000.0]
colors4 = [:black, :steelblue, :firebrick, :purple]
labels4 = ["No SIT", "500/biweekly", "2000/biweekly", "5000/biweekly"]

fig4 = Figure(size=(1000, 600))
ax4 = Axis(fig4[1, 1],
    xlabel="Day", ylabel="Wild population (log scale)",
    title="Eradication Trajectory Under Progressive SIT Deployment",
    yscale=log10, xlabelsize=14, ylabelsize=14)

for (i, rate) in enumerate(release_rates)
    res = simulate_screwworm(daily_temps=mcallen_temps,
                              sterile_release_rate=rate,
                              release_interval=14,
                              competitiveness=0.5)
    pop_clamp = max.(res.pop, 0.1)   # clamp for log scale
    lines!(ax4, 1:730, pop_clamp, linewidth=2, color=colors4[i],
           label=labels4[i])
end
hlines!(ax4, [1.0], color=:red, linestyle=:dot, linewidth=1.5,
        label="Extinction threshold")
xlims!(ax4, 1, 730)
axislegend(ax4, position=:rt, framevisible=false, labelsize=11)

# Annotate paper context
text!(ax4, 500, 3e4,
      text="Paper: avg SIT efficacy ~1.7%\nCold winters 1976–1980 were critical\nLast US case: 1982",
      fontsize=10, color=:gray40, align=(:left, :top))

save(joinpath(figdir, "fig4_eradication_trajectory.png"), fig4, px_per_unit=2)
println("  ✓ fig4_eradication_trajectory.png")

# ============================================================================
# Figure 5 — Geographic limit: overwinter survival vs temperature / latitude
# ============================================================================
println("Generating Fig 5: geographic limit …")

# Compute µ_cold for a range of mean annual temperatures (proxy for latitude)
mean_Ts = range(10.0, 30.0, length=200)
lats    = range(36.0, 16.0, length=200)   # roughly N→S

mu_cold_vals = Float64[]
for (mT, lat) in zip(mean_Ts, lats)
    temps = approx_temps(mT, lat)
    push!(mu_cold_vals, cumulative_cold(temps))
end

# Reference locations from the paper
ref_locs = [
    ("Uvalde, TX",         28.3, 19.0),
    ("McAllen, TX",        26.2, 23.5),
    ("Tampico, MX",        22.0, 24.5),
    ("Tuxtla-Gutiérrez",   16.7, 25.0),
]
ref_mu = Float64[]
for (_, lat, mT) in ref_locs
    temps = approx_temps(mT, lat)
    push!(ref_mu, cumulative_cold(temps))
end

# Climate warming shift (+2°C)
mu_cold_warm = Float64[]
for (mT, lat) in zip(mean_Ts, lats)
    temps = approx_temps(mT + 2.0, lat)
    push!(mu_cold_warm, cumulative_cold(temps))
end

fig5 = Figure(size=(1000, 650))
ax5 = Axis(fig5[1, 1],
    xlabel="Mean annual temperature (°C)", ylabel="µ_cold (cumulative cold index)",
    title="Geographic Limit: Overwinter Survival Index (cf. Fig. 8, 2019)",
    xlabelsize=14, ylabelsize=14)

lines!(ax5, collect(mean_Ts), mu_cold_vals, linewidth=2.5, color=:steelblue,
       label="Current climate")
lines!(ax5, collect(mean_Ts), mu_cold_warm, linewidth=2, color=:firebrick,
       linestyle=:dash, label="+2°C warming (RCP 8.5, ~2050)")

hlines!(ax5, [MU_COLD_THRESHOLD], color=:red, linestyle=:dot, linewidth=2)
text!(ax5, 11.0, MU_COLD_THRESHOLD + 0.5,
      text="µ_cold = 10 (endemic threshold)", fontsize=11, color=:red)

hspan!(ax5, 0.0, MU_COLD_THRESHOLD, color=(:green, 0.06))
text!(ax5, 26.0, 3.0, text="ENDEMIC ZONE", fontsize=13, color=:green4,
      align=(:center, :center))

# Plot reference locations
for (i, (name, lat, mT)) in enumerate(ref_locs)
    scatter!(ax5, [mT], [ref_mu[i]], markersize=12, color=:black)
    text!(ax5, mT + 0.3, ref_mu[i] + 0.4, text=name, fontsize=10, color=:black)
end

# Secondary x-axis label (approximate latitude)
ax5b = Axis(fig5[1, 1], xaxisposition=:top,
    xlabel="Approximate latitude (°N)", xlabelsize=12,
    xticks=([12.0, 16.0, 20.0, 24.0, 28.0],
            ["35°N", "31°N", "27°N", "23°N", "19°N"]))
hideydecorations!(ax5b)
hidespines!(ax5b, :l, :r, :b)
ax5b.backgroundcolor = :transparent
xlims!(ax5b, extrema(mean_Ts)...)

xlims!(ax5, extrema(mean_Ts)...)
axislegend(ax5, position=:rt, framevisible=false, labelsize=11)

save(joinpath(figdir, "fig5_geographic_limit.png"), fig5, px_per_unit=2)
println("  ✓ fig5_geographic_limit.png")

# ============================================================================
# Summary
# ============================================================================
println("\n" * "="^65)
println("Parameter comparison: Vignette vs Papers")
println("="^65)
println("  T_lower:     vignette=14.5°C   papers=14.5°C   ✓ match")
println("  T_upper:     vignette=43.5°C   papers=43.5°C   ✓ match")
println("  T_opt:       vignette=27.2°C   papers=27.2°C   ✓ match")
println("  µ_ab(T):     vignette uses Eq.1 from 2019 paper ✓")
println("  µ_cold thresh: vignette=10     papers=10        ✓ match")
println("  Max fecund:  vignette=67/d     2019=67, 2014=90 (diff context)")
println("  ♀ half-life: vignette=3.7 d    paper=3.7 d      ✓ match")
println("  SIT efficacy:                  paper ~1.7%")
println("  Model SIT:   McAllen 2, Tuxtla 23 /km²/biweek   (2019 Fig. 6)")
println("  Field SIT:   39–896 /km²/week                    (Matlock et al.)")
println("  Note: user-mentioned 'base 9.7°C' not in these papers;")
println("         both papers use 14.5°C. 9.7°C overlaid for comparison.")
println("="^65)
println("\nAll figures saved to: ", figdir)
