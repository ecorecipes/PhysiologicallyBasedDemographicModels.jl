#!/usr/bin/env julia
# Validation script for the Lobesia botrana overwintering vignette (05).
# Compares vignette parameters against values from:
#   - Baumgärtner et al. (2012) J. Ent. Acar. Res. 44(1)
#   - Gutierrez et al. (2017) Agric. For. Meteorol. 232:421–434

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using PhysiologicallyBasedDemographicModels
using CairoMakie

figdir = joinpath(@__DIR__, "figures", "lobesia")
mkpath(figdir)

# ============================================================================
# Paper parameters — Baumgärtner et al. (2012), Table 3
# ============================================================================
# Phase 1 (pre-diapause, mature larvae): α=0.00225, β=5, T_l=8.9, T_u=33.0, ξ=6.0
# Phase 2 (diapause, pupae):             α=3.026e-4, β=1.5, T_l=7.1, T_u=28.5, ξ=1.0
# Phase 3 (post-diapause, pupae):        α=0.00785, β=4.5, T_l=11.5, T_u=33.0, ξ=1.0

# The 2012 paper uses a modified Brière form:
#   r_j(T) = ξ_j · α_j · (T − T_l) / (1 + β_j · (T − T_l))
# This is NOT the standard Brière but a Beddington/rational form.
# We reconstruct it faithfully for comparison.

# Paper's modified development rate (Eq. 7 of Baumgärtner 2012)
function paper_devrate_2012(T, α, β, T_l, T_u, ξ)
    T <= T_l && return 0.0
    T >= T_u && return 0.0
    return ξ * α * (T - T_l) / (1.0 + β * (T - T_l))
end

# Paper phase parameters (Table 3)
const PAPER_PREDIAP  = (α=0.00225,   β=5.0,  T_l=8.9,  T_u=33.0, ξ=6.0)
const PAPER_DIAP     = (α=3.02579e-4, β=1.5, T_l=7.1,  T_u=28.5, ξ=1.0)
const PAPER_POSTDIAP = (α=0.00785,   β=4.5,  T_l=11.5, T_u=33.0, ξ=1.0)

# Vignette parameters (standard Brière form):
# r(T) = a · T · (T − T_lower) · √(T_upper − T)
const VIGNETTE_PREDIAP  = BriereDevelopmentRate(2.0e-4, 6.3, 28.0)
const VIGNETTE_DIAP     = BriereDevelopmentRate(1.0e-4, 6.2, 22.0)
const VIGNETTE_POSTDIAP = BriereDevelopmentRate(2.0e-4, 6.3, 28.0)

# Gutierrez et al. (2017) parameters
# Egg-larval: a=0.00225, T_l=8.9, T_u=33, b=5 (same as 2012 pre-diapause)
# Pupal: a=0.00785, T_l=11.5, T_u=33, b=4.5 (same as 2012 post-diapause)
# Diapause pupae: 115 dd>11.5°C developmental time (85% of 128.5 dd)
# Erlang k=25, σ=23 dd>11.5°C, emergence window: 69–161 dd>11.5°C

# Gutierrez (2017) mortality: μ(T) = c · ((T − T_m) / T_m)  with T_m=21.5
# c=2.2 (eggs/larvae), c=2.0 (pupae/adults)
const T_m_mort = 21.5
const c_EL = 2.2
const c_PA = 2.0

function gutierrez_mortality(T, c, T_m=T_m_mort)
    return clamp(c * ((T - T_m) / T_m), 0.0, 1.0)
end

# ============================================================================
# Photoperiod parameters
# Baumgärtner (2012) Eqs. 1–3:
#   DL_s = 9.83 + 0.1226 · L   (diapause induction starts)
#   DL_e = 7.66 + 0.1226 · L   (diapause induction complete)
# Gutierrez (2017) Eq. 6:
#   diap(dl) = 0.4565 · (14.15 − dl),  clamped [0,1] for dl ≤ 14.15 h
#   Critical day length = 14.15 h
# ============================================================================
const DL_CRIT_2017 = 14.15
const DIAP_COEFF   = 0.4565

DL_start(lat) = 9.83 + 0.1226 * lat
DL_end(lat)   = 7.66 + 0.1226 * lat

function diapause_probability_2017(dl)
    dl >= DL_CRIT_2017 && return 0.0
    return clamp(DIAP_COEFF * (DL_CRIT_2017 - dl), 0.0, 1.0)
end

# Cold hardiness — Gutierrez (2017)
# Diapause pupae supercooling point: −24.5°C
# Non-diapause pupae supercooling point: −22.5°C
const SCP_DIAPAUSE     = -24.5
const SCP_NON_DIAPAUSE = -22.5

# ============================================================================
# FIGURE 1: Development Rate vs Temperature
# ============================================================================
println("Generating Figure 1: Development rates...")

Ts = range(-5.0, 40.0, length=300)

# Paper rates (2012 modified Brière)
paper_pre  = [paper_devrate_2012(T, PAPER_PREDIAP...) for T in Ts]
paper_dia  = [paper_devrate_2012(T, PAPER_DIAP...) for T in Ts]
paper_post = [paper_devrate_2012(T, PAPER_POSTDIAP...) for T in Ts]

# Vignette rates (standard Brière)
vig_pre  = [development_rate(VIGNETTE_PREDIAP, T) for T in Ts]
vig_dia  = [development_rate(VIGNETTE_DIAP, T) for T in Ts]
vig_post = [development_rate(VIGNETTE_POSTDIAP, T) for T in Ts]

fig1 = Figure(size=(1000, 700))
ax1 = Axis(fig1[1,1],
    xlabel="Temperature (°C)", ylabel="Development rate (day⁻¹)",
    title="Lobesia Development Rates — Paper vs Vignette")

# Paper curves (dashed)
lines!(ax1, collect(Ts), paper_pre, color=:firebrick, linewidth=2, linestyle=:dash,
       label="Pre-diapause (2012)")
lines!(ax1, collect(Ts), paper_dia, color=:steelblue, linewidth=2, linestyle=:dash,
       label="Diapause (2012)")
lines!(ax1, collect(Ts), paper_post, color=:forestgreen, linewidth=2, linestyle=:dash,
       label="Post-diapause (2012)")

# Vignette curves (solid)
lines!(ax1, collect(Ts), vig_pre, color=:firebrick, linewidth=2,
       label="Pre-diapause (vignette)")
lines!(ax1, collect(Ts), vig_dia, color=:steelblue, linewidth=2,
       label="Diapause (vignette)")
lines!(ax1, collect(Ts), vig_post, color=:forestgreen, linewidth=2,
       label="Post-diapause (vignette)")

# Annotate paper threshold temperatures
vlines!(ax1, [8.9, 7.1, 11.5], color=:gray60, linestyle=:dot, linewidth=1)
text!(ax1, 9.0, maximum(paper_pre)*0.9, text="T_l=8.9", fontsize=10, color=:firebrick)
text!(ax1, 7.2, maximum(paper_dia)*0.9, text="T_l=7.1", fontsize=10, color=:steelblue)
text!(ax1, 11.6, maximum(paper_post)*0.8, text="T_l=11.5", fontsize=10, color=:forestgreen)

xlims!(ax1, -5, 40)
axislegend(ax1, position=:lt, framevisible=false)

save(joinpath(figdir, "fig1_development_rates.png"), fig1, px_per_unit=2)
println("  Saved fig1_development_rates.png")
println("  Paper pre-diapause peak: ", round(maximum(paper_pre), digits=4),
        " at T≈", round(Ts[argmax(paper_pre)], digits=1), "°C")
println("  Vignette pre-diapause peak: ", round(maximum(vig_pre), digits=4),
        " at T≈", round(Ts[argmax(vig_pre)], digits=1), "°C")

# ============================================================================
# FIGURE 2: Photoperiod Response — Diapause Induction
# ============================================================================
println("\nGenerating Figure 2: Photoperiod response...")

daylengths = range(8.0, 16.0, length=200)
diap_probs = [diapause_probability_2017(dl) for dl in daylengths]

# Latitude-specific critical day lengths (Baumgärtner 2012)
lats_table = [
    ("Pissouri, Cyprus",       34.7),
    ("Sicily",                 37.5),
    ("Bordeaux",               44.8),
    ("Wädenswil, Switzerland", 47.2),
    ("Dresden, Germany",       51.1),
]

fig2 = Figure(size=(1000, 700))
ax2 = Axis(fig2[1,1],
    xlabel="Day length (hours)", ylabel="Diapause induction probability",
    title="Photoperiod-Driven Diapause Induction")

lines!(ax2, collect(daylengths), diap_probs, color=:black, linewidth=2.5,
       label="Gutierrez (2017): 0.4565·(14.15−DL)")

# Mark critical day length
vlines!(ax2, [DL_CRIT_2017], color=:red, linestyle=:dash, linewidth=1.5,
        label="DL_crit = 14.15 h")

# Mark latitude-dependent DL_s and DL_e from 2012
colors_lat = [:purple, :orange, :teal, :royalblue, :brown]
for (i, (name, lat)) in enumerate(lats_table)
    dl_s = DL_start(lat)
    dl_e = DL_end(lat)
    # Shade the induction window
    band!(ax2, [dl_e, dl_s], [0.0, 0.0], [1.0, 1.0],
          color=(colors_lat[i], 0.08))
    vlines!(ax2, [dl_s], color=colors_lat[i], linestyle=:dot, linewidth=1.2)
    vlines!(ax2, [dl_e], color=colors_lat[i], linestyle=:dot, linewidth=1.2)
    text!(ax2, dl_s + 0.05, 0.95 - 0.07*i,
          text="$(name) DLs=$(round(dl_s, digits=1))",
          fontsize=9, color=colors_lat[i])
end

xlims!(ax2, 8, 16)
ylims!(ax2, -0.02, 1.05)
axislegend(ax2, position=:rb, framevisible=false)

save(joinpath(figdir, "fig2_photoperiod_response.png"), fig2, px_per_unit=2)
println("  Saved fig2_photoperiod_response.png")
for (name, lat) in lats_table
    println("  $(name) ($(lat)°N): DL_s=$(round(DL_start(lat), digits=2))h, DL_e=$(round(DL_end(lat), digits=2))h")
end

# ============================================================================
# FIGURE 3: Overwintering Mortality / Cold Survival
# ============================================================================
println("\nGenerating Figure 3: Overwintering mortality...")

Ts_cold = range(-30.0, 35.0, length=300)

# Gutierrez (2017) temperature-dependent mortality
mort_el = [gutierrez_mortality(T, c_EL) for T in Ts_cold]
mort_pa = [gutierrez_mortality(T, c_PA) for T in Ts_cold]

# Survival = 1 − mortality
surv_el = 1.0 .- mort_el
surv_pa = 1.0 .- mort_pa

# Daily survival over 90 winter days at constant temperature
surv_90d_el = [max(0.0, 1.0 - gutierrez_mortality(T, c_EL))^90 for T in Ts_cold]
surv_90d_pa = [max(0.0, 1.0 - gutierrez_mortality(T, c_PA))^90 for T in Ts_cold]

fig3 = Figure(size=(1000, 800))

# Panel a: Daily mortality rate
ax3a = Axis(fig3[1,1],
    xlabel="Temperature (°C)", ylabel="Daily mortality rate",
    title="(a) Temperature-Dependent Mortality (Gutierrez 2017)")
lines!(ax3a, collect(Ts_cold), mort_el, color=:firebrick, linewidth=2,
       label="Eggs/larvae (c=2.2)")
lines!(ax3a, collect(Ts_cold), mort_pa, color=:steelblue, linewidth=2,
       label="Pupae/adults (c=2.0)")
hlines!(ax3a, [1.0], color=:gray, linestyle=:dot)
vlines!(ax3a, [T_m_mort], color=:gray50, linestyle=:dash, linewidth=1,
        label="T_m = 21.5°C")
# Mark supercooling points
vlines!(ax3a, [SCP_DIAPAUSE], color=:navy, linestyle=:dashdot, linewidth=1.5,
        label="SCP diapause = $(SCP_DIAPAUSE)°C")
vlines!(ax3a, [SCP_NON_DIAPAUSE], color=:skyblue3, linestyle=:dashdot, linewidth=1.5,
        label="SCP non-diap = $(SCP_NON_DIAPAUSE)°C")
xlims!(ax3a, -30, 35)
ylims!(ax3a, 0, 1.05)
axislegend(ax3a, position=:rt, framevisible=false, labelsize=10)

# Panel b: Cumulative winter survival (90 days)
ax3b = Axis(fig3[2,1],
    xlabel="Mean winter temperature (°C)",
    ylabel="Proportion surviving 90 days",
    title="(b) 90-Day Winter Survival at Constant Temperature")
lines!(ax3b, collect(Ts_cold), surv_90d_el, color=:firebrick, linewidth=2,
       label="Eggs/larvae")
lines!(ax3b, collect(Ts_cold), surv_90d_pa, color=:steelblue, linewidth=2,
       label="Pupae (diapause)")

# Annotate expected overwinter survival range
band!(ax3b, [5.0, 15.0], [0.0, 0.0], [1.0, 1.0],
      color=(:green, 0.08), label="Typical winter T range")

xlims!(ax3b, -10, 30)
ylims!(ax3b, 0, 1.05)
axislegend(ax3b, position=:lt, framevisible=false, labelsize=10)

save(joinpath(figdir, "fig3_overwintering_mortality.png"), fig3, px_per_unit=2)
println("  Saved fig3_overwintering_mortality.png")
println("  Zero-mortality temperature (pupae): ",
        round(T_m_mort, digits=1), "°C")
println("  Mortality at 5°C (pupae): ",
        round(gutierrez_mortality(5.0, c_PA), digits=3))
println("  90-day survival at 5°C (pupae): ",
        round((1.0 - gutierrez_mortality(5.0, c_PA))^90, digits=3))

# ============================================================================
# FIGURE 4: Seasonal Dynamics — Mediterranean Climate
# ============================================================================
println("\nGenerating Figure 4: Seasonal dynamics...")

fig4 = Figure(size=(1200, 900))

# Simulate at three Mediterranean sites
sites = [
    ("Sicily (37.5°N)",      37.5, :orange),
    ("Bordeaux (44.8°N)",    44.8, :teal),
    ("Switzerland (47.2°N)", 47.2, :royalblue),
]

# Expected emergence days from Baumgärtner (2012)
expected_emergence = Dict(
    37.5 => (day=355, label="Dec 21 ± 10d"),   # Sicily cohort 1
    44.8 => (day=30,  label="Jan 30 ± 15d"),    # Bordeaux cohort 1
    47.2 => (day=54,  label="Feb 23 ± 12d"),    # Switzerland cohort 1
)

for (idx, (name, lat, col)) in enumerate(sites)
    mean_T = 16.0 - 0.15 * (lat - 35.0)
    amplitude = 8.0 + 0.15 * (lat - 35.0)

    sw = SinusoidalWeather(mean_T, amplitude; phase=200.0)

    # Build fresh cohort with vignette parameters
    phases = [
        LifeStage(:prediapause,  DistributedDelay(10, 50.0;  W0=100.0),
                  VIGNETTE_PREDIAP, 0.005),
        LifeStage(:diapause,     DistributedDelay(20, 200.0; W0=0.0),
                  VIGNETTE_DIAP, 0.002),
        LifeStage(:postdiapause, DistributedDelay(10, 50.0;  W0=0.0),
                  VIGNETTE_POSTDIAP, 0.003),
    ]
    cohort = Population(:lobesia, phases)

    n_sim_days = 300
    weather_days = [get_weather(sw, d) for d in 200:(200 + n_sim_days - 1)]
    ws = WeatherSeries(weather_days; day_offset=200)

    prob = PBDMProblem(cohort, ws, (200, 200 + n_sim_days - 1))
    sol = solve(prob, DirectIteration())

    # Extract trajectories
    days = sol.t
    pre_traj  = stage_trajectory(sol, 1)
    dia_traj  = stage_trajectory(sol, 2)
    post_traj = stage_trajectory(sol, 3)
    total     = total_population(sol)

    ax = Axis(fig4[idx, 1],
        xlabel= idx == 3 ? "Day of year" : "",
        ylabel="Population",
        title="$(name)")

    lines!(ax, days, pre_traj,  color=(col, 0.6), linewidth=1.5, linestyle=:dash,
           label="Pre-diapause")
    lines!(ax, days, dia_traj,  color=col, linewidth=2.5,
           label="Diapause")
    lines!(ax, days, post_traj, color=(col, 0.8), linewidth=2, linestyle=:dot,
           label="Post-diapause")
    lines!(ax, days, total,     color=:black, linewidth=1.5, linestyle=:dashdot,
           label="Total")

    # Mark expected emergence from paper
    if haskey(expected_emergence, lat)
        exp = expected_emergence[lat]
        exp_day = exp.day
        # Wrap into simulation range (days 200–499 → handle 355, or 30+365=395, etc.)
        plot_day = exp_day < 200 ? exp_day + 365 : exp_day
        if plot_day >= days[1] && plot_day <= days[end]
            vlines!(ax, [plot_day], color=:red, linestyle=:dash, linewidth=1.5,
                    label="Paper: $(exp.label)")
        end
    end

    # Mark model emergence
    emergence = phenology(sol; threshold=0.5)
    if emergence !== nothing
        vlines!(ax, [emergence], color=:darkgreen, linestyle=:dash, linewidth=1.5,
                label="Model 50% emergence: day $(emergence)")
    end

    xlims!(ax, 200, 500)
    axislegend(ax, position=:rt, framevisible=false, labelsize=9)

    # Temperature panel on right
    ax_t = Axis(fig4[idx, 2],
        xlabel= idx == 3 ? "Day of year" : "",
        ylabel="Temperature (°C)",
        title="Temperature — $(name)")
    temps = [get_weather(sw, d).T_mean for d in days[1]:days[end]]
    lines!(ax_t, collect(days[1]:days[end]), temps, color=col, linewidth=1.5)
    hlines!(ax_t, [0.0], color=:gray, linestyle=:dot)
    xlims!(ax_t, 200, 500)

    # Print summary
    cdd = cumulative_degree_days(sol)
    surv = sum(sol.u[end]) / max(sum(sol.u[1]), 1e-10) * 100
    println("  $(name):")
    println("    Total DD: $(round(cdd[end], digits=0))")
    println("    Model emergence: $(emergence !== nothing ? "day $(emergence)" : "not reached")")
    if haskey(expected_emergence, lat)
        exp = expected_emergence[lat]
        println("    Paper expected:  day $(exp.day) ($(exp.label))")
    end
    println("    Survival: $(round(surv, digits=1))%")
end

save(joinpath(figdir, "fig4_seasonal_dynamics.png"), fig4, px_per_unit=2)
println("  Saved fig4_seasonal_dynamics.png")

# ============================================================================
# FIGURE 5: Climate Comparison — Swiss vs Italian Overwinter Survival
# ============================================================================
println("\nGenerating Figure 5: Climate comparison...")

fig5 = Figure(size=(1200, 800))

# Site definitions
north = (name="Wädenswil, Switzerland (47.2°N)", lat=47.2,
         mean_T=14.6, amplitude=9.8, col=:royalblue)
south = (name="Western Sicily (37.5°N)", lat=37.5,
         mean_T=15.6, amplitude=8.4, col=:orange)

ax5a = Axis(fig5[1,1],
    xlabel="Day of year", ylabel="Population",
    title="Overwintering Cohort — Northern (Switzerland) vs Southern (Sicily)")

ax5b = Axis(fig5[2,1],
    xlabel="Day of year", ylabel="Temperature (°C)",
    title="Temperature Profiles")

# Also track day length
ax5c = Axis(fig5[1,2],
    xlabel="Day of year", ylabel="Day length (hours)",
    title="Photoperiod & Diapause Thresholds")

ax5d = Axis(fig5[2,2],
    xlabel="Day of year", ylabel="Cumulative degree-days",
    title="Cumulative Degree-Day Accumulation")

for site in [north, south]
    sw = SinusoidalWeather(site.mean_T, site.amplitude; phase=200.0)

    phases = [
        LifeStage(:prediapause,  DistributedDelay(10, 50.0;  W0=100.0),
                  VIGNETTE_PREDIAP, 0.005),
        LifeStage(:diapause,     DistributedDelay(20, 200.0; W0=0.0),
                  VIGNETTE_DIAP, 0.002),
        LifeStage(:postdiapause, DistributedDelay(10, 50.0;  W0=0.0),
                  VIGNETTE_POSTDIAP, 0.003),
    ]
    cohort = Population(:lobesia, phases)

    n_sim_days = 300
    weather_days = [get_weather(sw, d) for d in 200:(200 + n_sim_days - 1)]
    ws = WeatherSeries(weather_days; day_offset=200)

    prob = PBDMProblem(cohort, ws, (200, 200 + n_sim_days - 1))
    sol = solve(prob, DirectIteration())

    days = sol.t
    total = total_population(sol)
    cdd = cumulative_degree_days(sol)
    temps = [get_weather(sw, d).T_mean for d in days[1]:days[end]]

    # Compute day lengths
    # Days wrap: 200–365 = year 1, 366–499 = days 1–134 year 2
    local daylengths = Float64[]
    for d in days[1]:days[end]
        doy = d <= 365 ? d : d - 365
        push!(daylengths, photoperiod(site.lat, doy))
    end

    surv = sum(sol.u[end]) / max(sum(sol.u[1]), 1e-10) * 100
    emergence = phenology(sol; threshold=0.5)

    # Population
    lines!(ax5a, days, total, color=site.col, linewidth=2.5,
           label="$(site.name) — surv=$(round(surv, digits=1))%")
    if emergence !== nothing
        scatter!(ax5a, [emergence], [total[findfirst(==(emergence), days)]],
                 color=site.col, markersize=12, marker=:diamond)
    end

    # Temperature
    lines!(ax5b, collect(days[1]:days[end]), temps, color=site.col, linewidth=2,
           label=site.name)

    # Photoperiod
    lines!(ax5c, collect(days[1]:days[end]), daylengths, color=site.col, linewidth=2,
           label="$(site.name)")

    # Cumulative DD
    lines!(ax5d, days[1:length(cdd)], cdd, color=site.col, linewidth=2,
           label=site.name)

    println("  $(site.name):")
    println("    Survival: $(round(surv, digits=1))%")
    println("    Emergence: $(emergence !== nothing ? "day $(emergence)" : "not reached")")
    println("    Final CDD: $(round(cdd[end], digits=0))")
end

# Add reference lines to photoperiod panel
# Mark DL_s thresholds for both latitudes
hlines!(ax5c, [DL_start(47.2)], color=:royalblue, linestyle=:dash, linewidth=1,
        label="DLs Swiss=$(round(DL_start(47.2), digits=1))h")
hlines!(ax5c, [DL_start(37.5)], color=:orange, linestyle=:dash, linewidth=1,
        label="DLs Sicily=$(round(DL_start(37.5), digits=1))h")
hlines!(ax5c, [DL_CRIT_2017], color=:red, linestyle=:dot, linewidth=1,
        label="Gutierrez DL_crit=14.15h")

# Frost line on temperature
hlines!(ax5b, [0.0], color=:gray, linestyle=:dot, linewidth=1)

for ax in [ax5a, ax5b, ax5c, ax5d]
    xlims!(ax, 200, 500)
    axislegend(ax, position=:rt, framevisible=false, labelsize=9)
end

save(joinpath(figdir, "fig5_climate_comparison.png"), fig5, px_per_unit=2)
println("  Saved fig5_climate_comparison.png")

# ============================================================================
# Summary: Parameter comparison table
# ============================================================================
println("\n" * "="^70)
println("PARAMETER COMPARISON: PAPER vs VIGNETTE")
println("="^70)
println()
println("Development rate form:")
println("  Paper (2012):    r(T) = ξ·α·(T−T_l) / (1+β·(T−T_l))  [rational/Beddington]")
println("  Vignette:        r(T) = a·T·(T−T_l)·√(T_u−T)          [standard Brière]")
println()
println("                      Paper (2012)          Vignette")
println("                      T_l    T_u             T_lower  T_upper")
println("  Pre-diapause:       8.9    33.0            6.3      28.0")
println("  Diapause:           7.1    28.5            6.2      22.0")
println("  Post-diapause:     11.5    33.0            6.3      28.0")
println()
println("Key differences:")
println("  ‣ Vignette uses lower T_lower thresholds (6.2–6.3 vs 7.1–11.5°C)")
println("  ‣ Vignette uses lower T_upper for diapause (22 vs 28.5°C)")
println("  ‣ Different functional forms (standard Brière vs rational)")
println("  ‣ 2012 paper has distinct post-diapause T_l=11.5°C; vignette merges with pre-diapause")
println()
println("Photoperiod:")
println("  2012: DL_s = 9.83+0.1226·L, DL_e = 7.66+0.1226·L  [latitude-dependent]")
println("  2017: DL_crit = 14.15 h, slope = 0.4565             [fixed threshold]")
println()
println("Mortality (2017): μ(T) = c·(T−21.5)/21.5  [c=2.2 eggs, 2.0 pupae]")
println("Supercooling points: diapause pupae −24.5°C, non-diapause −22.5°C")
println()
println("All figures saved to: ", figdir)
