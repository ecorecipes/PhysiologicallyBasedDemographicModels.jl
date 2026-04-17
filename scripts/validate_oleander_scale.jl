#!/usr/bin/env julia
# Validation script for the Oleander Scale (Aspidiotus nerii) tritrophic PBDM
# with competing parasitoids Aphytis chilensis and Coccophagoides utilis.
#
# Literature references:
#   Rochat, J. & Gutierrez, A.P. (2001) J. Anim. Ecol. 70:476–490.
#   Gutierrez, A.P. & Pizzamiglio, M.A. (2007) Neotrop. Entomol. 36:70–83.
#   DeBach, P. & Sundby, R.A. (1963) Hilgardia 34:105–166.
#   Murdoch, W.W. et al. (2005) Biological Control (Aphytis melinus proxy).
#
# Generates five PNG figures in scripts/figures/oleander_scale/.

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using CairoMakie
using Statistics

figdir = joinpath(@__DIR__, "figures", "oleander_scale")
mkpath(figdir)

# ============================================================
# Brière development rate:  r(T) = a * T * (T - T_L) * √(T_U - T)
# Returns 0 outside [T_L, T_U].
# ============================================================

function briere_rate(a, T_L, T_U, T)
    (T <= T_L || T >= T_U) && return 0.0
    return a * T * (T - T_L) * sqrt(T_U - T)
end

# ============================================================
# Parameters from Rochat & Gutierrez (2001), Table 1 / text
# ============================================================

# --- Oleander scale (Aspidiotus nerii) ---
const SCALE_TL = 10.0     # lower developmental threshold (°C)
const SCALE_TU = 35.0     # upper developmental threshold (°C)
const SCALE_A  = 0.000035 # Brière coefficient

const SCALE_DD_CRAWLER = 95.0   # DD >10 °C
const SCALE_DD_NYMPH1  = 135.0
const SCALE_DD_NYMPH2  = 175.0
const SCALE_DD_ADULT   = 200.0
const SCALE_PREOVIPOSITION = SCALE_DD_CRAWLER + SCALE_DD_NYMPH1 + SCALE_DD_NYMPH2  # 405 DD

const SCALE_MU_CRAWLER = 0.004  # per-DD background mortality
const SCALE_MU_NYMPH   = 0.001
const SCALE_MU_ADULT   = 0.002

const SCALE_FECUNDITY_MAX   = 4.0   # crawlers/♀/day at optimum
const SCALE_FECUND_TL       = 12.0
const SCALE_FECUND_OPT      = 25.0
const SCALE_FECUND_TU       = 34.0
const SCALE_LIFETIME_FECUND = 100.0

# --- Aphytis chilensis (warm-adapted ectoparasitoid) ---
const AC_TL = 11.0
const AC_TU = 34.0
const AC_A  = 0.000032
const AC_DD_IMMATURE = 200.0   # DD >11 °C
const AC_DD_ADULT    = 125.0
const AC_SEARCH      = 0.15
const AC_SEX_RATIO   = 0.60
const AC_MU_IMM      = 0.002
const AC_MU_ADULT    = 0.003

# --- Coccophagoides utilis (cool-adapted endoparasitoid) ---
const CU_TL = 8.0
const CU_TU = 31.0
const CU_A  = 0.000028
const CU_DD_IMMATURE = 225.0   # DD >8 °C
const CU_DD_ADULT    = 140.0
const CU_SEARCH      = 0.10
const CU_SEX_RATIO   = 0.55
const CU_MU_IMM      = 0.002
const CU_MU_ADULT    = 0.003

# ============================================================
# Helper functions
# ============================================================

function scale_fecundity(T)
    (T < SCALE_FECUND_TL || T > SCALE_FECUND_TU) && return 0.0
    T_mid = (SCALE_FECUND_TU + SCALE_FECUND_TL) / 2.0
    φ = max(0.0, 1.0 - ((T - T_mid) / (T_mid - SCALE_FECUND_TL))^2)
    return SCALE_FECUNDITY_MAX * φ
end

# Temperature-dependent mortality multiplier (increases near thermal limits)
function thermal_mortality(T, T_L, T_U; base=0.01, k_heat=4.0, k_cold=2.0)
    if T <= T_L
        return 1.0  # complete mortality below threshold
    elseif T >= T_U
        return 1.0
    end
    range = T_U - T_L
    mid = (T_L + T_U) / 2.0
    # Asymmetric U-shaped mortality: higher at heat extreme
    frac_cold = max(0.0, (T_L + 0.15 * range - T) / (0.15 * range))
    frac_heat = max(0.0, (T - (T_U - 0.20 * range)) / (0.20 * range))
    return base + (1.0 - base) * (frac_cold^k_cold + frac_heat^k_heat)
end

# Weather generators (return NamedTuple with T_mean and T_max)
function coastal_weather(day)
    doy = mod(day - 1, 365) + 1
    T_mean = 14.5 + 3.5 * sin(2π * (doy - 100) / 365)
    T_max = T_mean + 4.0
    return (T_mean=T_mean, T_max=T_max)
end

function inland_weather(day)
    doy = mod(day - 1, 365) + 1
    T_mean = 18.0 + 10.0 * sin(2π * (doy - 100) / 365)
    T_max = T_mean + 10.0  # Fresno summer maxima reach 38°C
    return (T_mean=T_mean, T_max=T_max)
end

# ============================================================
# Parasitoid emergence uses a discrete delay ring buffer so that
# hosts parasitized on day d yield adult parasitoids after the
# immature development period elapses (tracked in DD).
# ============================================================

mutable struct DelayBuffer
    buf::Vector{Float64}   # parasitized-host cohorts by day
    dd_acc::Vector{Float64} # cumulative DD for each cohort
    dd_required::Float64
    head::Int
    len::Int
end

function DelayBuffer(dd_required; maxlen=120)
    DelayBuffer(zeros(maxlen), zeros(maxlen), dd_required, 1, maxlen)
end

function push_cohort!(db::DelayBuffer, n_parasitized, dd_today)
    db.buf[db.head] = n_parasitized
    db.dd_acc[db.head] = 0.0
    db.head = mod1(db.head + 1, db.len)
    # Accumulate DD on all existing cohorts and collect emerged
    emerged = 0.0
    for i in 1:db.len
        db.dd_acc[i] += dd_today
        if db.buf[i] > 0.0 && db.dd_acc[i] >= db.dd_required
            emerged += db.buf[i]
            db.buf[i] = 0.0
            db.dd_acc[i] = 0.0
        end
    end
    return emerged
end

# ============================================================
# Two-parasitoid simulation with delayed emergence
# ============================================================

function run_two_parasitoids(; weather_fn, n_years=3, K=500.0, S0=50.0, Pa0=2.0, Pc0=2.0)
    n_days = 365 * n_years
    scale_pop = zeros(n_days)
    ac_pop    = zeros(n_days)
    cu_pop    = zeros(n_days)

    S  = S0
    Pa = Pa0
    Pc = Pc0

    db_ac = DelayBuffer(AC_DD_IMMATURE)
    db_cu = DelayBuffer(CU_DD_IMMATURE)

    # Temperature-dependent search efficiency (Gaussian within thermal window)
    function ac_search_eff(T)
        (T <= AC_TL || T >= AC_TU) && return 0.0
        opt = 27.0  # warm-adapted
        σ = 10.0
        return AC_SEARCH * exp(-0.5 * ((T - opt) / σ)^2)
    end
    function cu_search_eff(T)
        (T <= CU_TL || T >= CU_TU) && return 0.0
        opt = 20.0  # cool-adapted
        σ = 8.0
        return CU_SEARCH * exp(-0.5 * ((T - opt) / σ)^2)
    end

    for day in 1:n_days
        wx = weather_fn(day)
        T = wx isa NamedTuple ? wx.T_mean : wx
        T_max = wx isa NamedTuple ? wx.T_max : T + 5.0

        dd_scale = max(0.0, T - SCALE_TL)
        dd_ac    = max(0.0, T - AC_TL)
        dd_cu    = max(0.0, T - CU_TL)

        # Scale growth
        growth = scale_fecundity(T) * 0.1
        S += S * growth * (1.0 - S / K) - SCALE_MU_NYMPH * dd_scale * S

        available = max(0.0, S * 0.6)

        a_ac = ac_search_eff(T)
        a_cu = cu_search_eff(T)

        if a_ac >= a_cu
            attacked_ac = 0.0
            if a_ac > 0 && Pa > 0.01 && available > 0.1
                fr = 1.0 - exp(-a_ac * Pa)
                attacked_ac = min(fr * available, 0.5 * available)
            end
            remaining = max(0.0, available - attacked_ac)
            attacked_cu = 0.0
            if a_cu > 0 && Pc > 0.01 && remaining > 0.1
                fr = 1.0 - exp(-a_cu * Pc)
                attacked_cu = min(fr * remaining, 0.5 * remaining)
            end
        else
            attacked_cu = 0.0
            if a_cu > 0 && Pc > 0.01 && available > 0.1
                fr = 1.0 - exp(-a_cu * Pc)
                attacked_cu = min(fr * available, 0.5 * available)
            end
            remaining = max(0.0, available - attacked_cu)
            attacked_ac = 0.0
            if a_ac > 0 && Pa > 0.01 && remaining > 0.1
                fr = 1.0 - exp(-a_ac * Pa)
                attacked_ac = min(fr * remaining, 0.5 * remaining)
            end
        end

        S = max(0.0, S - attacked_ac - attacked_cu)

        emerged_ac = push_cohort!(db_ac, attacked_ac * AC_SEX_RATIO, dd_ac)
        emerged_cu = push_cohort!(db_cu, attacked_cu * CU_SEX_RATIO, dd_cu)

        # Mortality uses T_max for heat stress (critical for CU upper limit)
        # C. utilis: TU=31°C, severe mortality when T_max > 31°C
        cu_heat = T_max > CU_TU ? 0.25 * ((T_max - CU_TU) / 3.0)^1.5 : 0.0
        # Dormancy below threshold — very low mortality (overwintering in host)
        mort_cu = T < CU_TL ? 0.005 : 0.03 + cu_heat

        # A. chilensis: TU=34°C, tolerates heat well
        ac_heat = T_max > AC_TU ? 0.20 * ((T_max - AC_TU) / 3.0)^1.5 : 0.0
        # Dormancy below threshold — low mortality
        mort_ac = T < AC_TL ? 0.005 : 0.03 + ac_heat

        Pa = max(0.0, Pa + emerged_ac - mort_ac * Pa)
        Pc = max(0.0, Pc + emerged_cu - mort_cu * Pc)
        S  = max(0.0, S)

        scale_pop[day] = S
        ac_pop[day]    = Pa
        cu_pop[day]    = Pc
    end
    return scale_pop, ac_pop, cu_pop
end

# Constant-temperature wrapper
function run_constant_temp(; T_const, n_days=365, kwargs...)
    run_two_parasitoids(; weather_fn=(_)->(T_mean=T_const, T_max=T_const+5.0), n_years=1, kwargs...)
end

# ============================================================
# FIGURE 1: Development Rate Comparison
# ============================================================

println("── Figure 1: Development rate comparison ──")

Ts = range(5.0, 40.0, length=500)
r_scale = [briere_rate(SCALE_A, SCALE_TL, SCALE_TU, T) for T in Ts]
r_ac    = [briere_rate(AC_A, AC_TL, AC_TU, T) for T in Ts]
r_cu    = [briere_rate(CU_A, CU_TL, CU_TU, T) for T in Ts]

# Literature verification points
# A. nerii at 25°C: 48-51 days crawler→adult ⟹ total pre-adult = 405 DD / (25-10) = 27 days
# per-day rate for total organism ≈ 1/48 to 1/51
lit_scale_T = [25.0, 30.0]
lit_scale_r = [1.0 / 49.5, 0.0]  # ~0.020 at 25°C; negligible at 30°C (upper thermal stress)
lit_scale_label = ["A. nerii 25°C:\n48-51 d (lit.)", "A. nerii 30°C:\n~0 (lit.)"]

# Aphytis melinus proxy for A. chilensis: 13-18 days development, base ~13°C
lit_ac_T = [25.0]
lit_ac_r = [1.0 / 15.5]  # midpoint of 13-18 d ≈ 0.065

fig1 = Figure(size=(900, 500))
ax1 = Axis(fig1[1, 1],
    title="Temperature-Dependent Development Rates\n(Oleander Scale and Competing Parasitoids)",
    xlabel="Temperature (°C)",
    ylabel="Development rate (1/day)",
    xlabelsize=14, ylabelsize=14)

lines!(ax1, collect(Ts), r_scale, linewidth=3.0, color=:red,
       label="A. nerii (Tₗ=10°C, Tᵤ=35°C)")
lines!(ax1, collect(Ts), r_ac, linewidth=2.5, color=:blue,
       label="A. chilensis (Tₗ=11°C, Tᵤ=34°C)")
lines!(ax1, collect(Ts), r_cu, linewidth=2.5, color=:forestgreen,
       label="C. utilis (Tₗ=8°C, Tᵤ=31°C)")

# Literature overlays
scatter!(ax1, lit_scale_T, lit_scale_r, color=:red, markersize=14,
         marker=:star5, label="A. nerii lit. (DeBach 1963)")
scatter!(ax1, lit_ac_T, lit_ac_r, color=:blue, markersize=14,
         marker=:diamond, label="Aphytis lit. (Murdoch 2005)")

# Shade thermal niche zones
vspan!(ax1, 15.0, 20.0, color=(:green, 0.08))
text!(ax1, 17.5, maximum(r_scale) * 0.92,
    text="Coastal\nzone", align=(:center, :top), fontsize=9, color=:green)
vspan!(ax1, 25.0, 32.0, color=(:orange, 0.08))
text!(ax1, 28.5, maximum(r_scale) * 0.92,
    text="Inland\nzone", align=(:center, :top), fontsize=9, color=:darkorange)

# Mark C. utilis upper limit
vlines!(ax1, [CU_TU], color=(:forestgreen, 0.6), linestyle=:dash, linewidth=1)
text!(ax1, CU_TU + 0.3, maximum(r_cu) * 0.5,
    text="C. utilis\nTᵤ=31°C", fontsize=8, color=:forestgreen)

xlims!(ax1, 5, 40)
ylims!(ax1, 0, nothing)
axislegend(ax1, position=:lt, framevisible=true, labelsize=9)

save(joinpath(figdir, "devrate_comparison.png"), fig1, px_per_unit=2)
println("  Saved devrate_comparison.png")
println("  Scale peak: $(round(maximum(r_scale), digits=4)) at $(round(collect(Ts)[argmax(r_scale)], digits=1))°C")
println("  A. chilensis peak: $(round(maximum(r_ac), digits=4)) at $(round(collect(Ts)[argmax(r_ac)], digits=1))°C")
println("  C. utilis peak: $(round(maximum(r_cu), digits=4)) at $(round(collect(Ts)[argmax(r_cu)], digits=1))°C")

# Literature check: A. nerii at 25°C
r25_scale = briere_rate(SCALE_A, SCALE_TL, SCALE_TU, 25.0)
dd25 = 25.0 - SCALE_TL  # 15 DD/day
days_to_adult = SCALE_PREOVIPOSITION / dd25
println("  Lit. check — A. nerii at 25°C: $(round(days_to_adult, digits=1)) days crawler→adult (expect 48–51)")
println("  Lit. check — Brière rate at 25°C: $(round(r25_scale, digits=4))")

# ============================================================
# FIGURE 2: Mortality Curves
# ============================================================

println("\n── Figure 2: Temperature-dependent mortality curves ──")

Ts_mort = range(5.0, 40.0, length=400)
mu_scale = [thermal_mortality(T, SCALE_TL, SCALE_TU; base=0.01, k_heat=4.0, k_cold=2.0) for T in Ts_mort]
mu_ac    = [thermal_mortality(T, AC_TL, AC_TU; base=0.01, k_heat=4.0, k_cold=2.0) for T in Ts_mort]
mu_cu    = [thermal_mortality(T, CU_TL, CU_TU; base=0.01, k_heat=3.0, k_cold=2.0) for T in Ts_mort]

# Literature: no significant A. nerii development at 30°C → high stress
lit_mort_T = [30.0]
lit_mort_val = [thermal_mortality(30.0, SCALE_TL, SCALE_TU; base=0.01, k_heat=4.0, k_cold=2.0)]

fig2 = Figure(size=(900, 500))
ax2 = Axis(fig2[1, 1],
    title="Temperature-Dependent Mortality Index\n(Oleander Scale and Parasitoids)",
    xlabel="Temperature (°C)",
    ylabel="Relative mortality (0 = low, 1 = lethal)",
    xlabelsize=14, ylabelsize=14)

lines!(ax2, collect(Ts_mort), mu_scale, linewidth=3.0, color=:red,
       label="A. nerii (10–35°C)")
lines!(ax2, collect(Ts_mort), mu_ac, linewidth=2.5, color=:blue,
       label="A. chilensis (11–34°C)")
lines!(ax2, collect(Ts_mort), mu_cu, linewidth=2.5, color=:forestgreen,
       label="C. utilis (8–31°C)")

# Literature annotation: C. utilis suffers high mortality above 30°C
vspan!(ax2, 30.0, 35.0, color=(:red, 0.06))
text!(ax2, 32.5, 0.55,
    text="C. utilis heat\nstress zone\n(Rochat 2001)", fontsize=8, color=:red,
    align=(:center, :center))

# Annotate: A. chilensis advantage in heat
bracket!(ax2, 25.0, 0.05, 32.0, 0.05, text="A. chilensis\nadvantage",
         orientation=:down, fontsize=8, color=:blue, textcolor=:blue,
         linewidth=1, style=:square)

# Shade climate zones
vspan!(ax2, 8.0, 15.0, color=(:green, 0.06))
text!(ax2, 11.5, 0.85, text="Cool\ncoastal", fontsize=9, color=:green,
      align=(:center, :center))
vspan!(ax2, 25.0, 38.0, color=(:orange, 0.06))
text!(ax2, 31.0, 0.85, text="Hot\ninland", fontsize=9, color=:darkorange,
      align=(:center, :center))

xlims!(ax2, 5, 40)
ylims!(ax2, 0, 1.05)
axislegend(ax2, position=:ct, framevisible=true, labelsize=9)

save(joinpath(figdir, "mortality_curves.png"), fig2, px_per_unit=2)
println("  Saved mortality_curves.png")
println("  C. utilis mortality at 30°C: $(round(thermal_mortality(30.0, CU_TL, CU_TU; base=0.01, k_heat=3.0, k_cold=2.0), digits=3))")
println("  A. chilensis mortality at 30°C: $(round(thermal_mortality(30.0, AC_TL, AC_TU; base=0.01, k_heat=4.0, k_cold=2.0), digits=3))")

# ============================================================
# FIGURE 3: Constant 25°C Simulation — 365 days
# ============================================================

println("\n── Figure 3: Constant 25°C tritrophic simulation ──")

n_days_const = 365
s25, ac25, cu25 = run_constant_temp(T_const=25.0, n_days=n_days_const, S0=50.0, Pa0=2.0, Pc0=2.0)
days_25 = 1:n_days_const

# Literature overlay: fecundity check
# At 25°C, 15 DD/day → adult lifespan = 200/15 ≈ 13.3 days
# Scale fecundity at 25°C:
f25 = scale_fecundity(25.0)
adult_lifespan_days = SCALE_DD_ADULT / (25.0 - SCALE_TL)
lifetime_fecundity_25 = f25 * adult_lifespan_days
println("  Scale fecundity at 25°C: $(round(f25, digits=2)) crawlers/♀/day")
println("  Adult lifespan at 25°C: $(round(adult_lifespan_days, digits=1)) days")
println("  Estimated lifetime fecundity: $(round(lifetime_fecundity_25, digits=1)) crawlers/♀ (lit: ~28–100)")

fig3 = Figure(size=(900, 500))
ax3 = Axis(fig3[1, 1],
    title="Tritrophic Dynamics at Constant 25°C — Scale + Both Parasitoids",
    xlabel="Day",
    ylabel="Population density",
    xlabelsize=14, ylabelsize=14)

lines!(ax3, collect(days_25), s25, linewidth=2.5, color=:red, label="A. nerii (scale)")
lines!(ax3, collect(days_25), ac25, linewidth=2.0, color=:blue, label="A. chilensis")
lines!(ax3, collect(days_25), cu25, linewidth=2.0, color=:forestgreen, label="C. utilis")

# Literature annotations
pk_scale = maximum(s25)
pk_day_scale = argmax(s25)
text!(ax3, min(pk_day_scale + 15, 340), pk_scale * 0.85,
    text="Lit: A. nerii 48–51 d\ncrawler→adult at 25°C\nFecundity ~28 crawlers/♀",
    fontsize=8, color=:gray50, align=(:left, :top))

# Aphytis development time annotation
ac_dev_days = AC_DD_IMMATURE / (25.0 - AC_TL)
text!(ax3, 250, maximum(ac25) * 1.2,
    text="Aphytis dev: $(round(ac_dev_days, digits=0))d\n(lit: 13–18 d, Murdoch)",
    fontsize=8, color=:blue, align=(:center, :bottom))

xlims!(ax3, 1, n_days_const)
axislegend(ax3, position=:rt, framevisible=true, labelsize=10)

save(joinpath(figdir, "sim_constant_25C.png"), fig3, px_per_unit=2)
println("  Saved sim_constant_25C.png")
println("  Scale peak: $(round(pk_scale, digits=1)) at day $(pk_day_scale)")
println("  A. chilensis year-end: $(round(ac25[end], digits=2))")
println("  C. utilis year-end: $(round(cu25[end], digits=2))")

# ============================================================
# FIGURE 4: Coastal vs Inland Comparison
# ============================================================

println("\n── Figure 4: Coastal vs Inland climate comparison ──")

s_coast, ac_coast, cu_coast    = run_two_parasitoids(weather_fn=coastal_weather, n_years=3)
s_inland, ac_inland, cu_inland = run_two_parasitoids(weather_fn=inland_weather, n_years=3)
days_3yr = 1:(365*3)

yr3 = 731:1095
println("  Year 3 means (Coastal):")
println("    Scale: $(round(mean(s_coast[yr3]), digits=1))")
println("    A. chilensis: $(round(mean(ac_coast[yr3]), digits=2))")
println("    C. utilis: $(round(mean(cu_coast[yr3]), digits=2))")
println("  Year 3 means (Inland):")
println("    Scale: $(round(mean(s_inland[yr3]), digits=1))")
println("    A. chilensis: $(round(mean(ac_inland[yr3]), digits=2))")
println("    C. utilis: $(round(mean(cu_inland[yr3]), digits=2))")

fig4 = Figure(size=(1000, 650))

# Top left: Coastal populations
ax4a = Axis(fig4[1, 1], title="(a) Coastal (Berkeley) — Both Parasitoids",
    ylabel="Density", xlabelvisible=false)
lines!(ax4a, collect(days_3yr), s_coast, linewidth=2.5, color=:red, label="Scale")
lines!(ax4a, collect(days_3yr), ac_coast, linewidth=2.0, color=:blue, label="A. chilensis")
lines!(ax4a, collect(days_3yr), cu_coast, linewidth=2.0, color=:forestgreen, label="C. utilis")
axislegend(ax4a, position=:rt, framevisible=true, labelsize=9)

# Top right: Inland populations
ax4b = Axis(fig4[1, 2], title="(b) Inland (Fresno) — Both Parasitoids",
    xlabelvisible=false)
lines!(ax4b, collect(days_3yr), s_inland, linewidth=2.5, color=:red, label="Scale")
lines!(ax4b, collect(days_3yr), ac_inland, linewidth=2.0, color=:blue, label="A. chilensis")
lines!(ax4b, collect(days_3yr), cu_inland, linewidth=2.0, color=:forestgreen, label="C. utilis")
axislegend(ax4b, position=:rt, framevisible=true, labelsize=9)

# Bottom left: Coastal temperature profile
ax4c = Axis(fig4[2, 1], title="(c) Coastal Temperature Profile",
    xlabel="Day", ylabel="Temperature (°C)")
c_temps = [coastal_weather(d).T_mean for d in 1:365]
c_tmax  = [coastal_weather(d).T_max for d in 1:365]
lines!(ax4c, collect(1:365), c_temps, linewidth=2, color=:steelblue, label="T_mean")
lines!(ax4c, collect(1:365), c_tmax, linewidth=1, color=:steelblue, linestyle=:dash, label="T_max")
hlines!(ax4c, [CU_TL], color=:forestgreen, linestyle=:dot, linewidth=1, label="C. utilis Tₗ=8°C")
hlines!(ax4c, [AC_TL], color=:blue, linestyle=:dot, linewidth=1, label="A. chilensis Tₗ=11°C")
axislegend(ax4c, position=:rb, framevisible=true, labelsize=8)

# Bottom right: Inland temperature profile
ax4d = Axis(fig4[2, 2], title="(d) Inland Temperature Profile",
    xlabel="Day", ylabel="Temperature (°C)")
i_temps = [inland_weather(d).T_mean for d in 1:365]
i_tmax  = [inland_weather(d).T_max for d in 1:365]
lines!(ax4d, collect(1:365), i_temps, linewidth=2, color=:coral, label="T_mean")
lines!(ax4d, collect(1:365), i_tmax, linewidth=1, color=:coral, linestyle=:dash, label="T_max")
hlines!(ax4d, [CU_TU], color=:forestgreen, linestyle=:dash, linewidth=1, label="C. utilis Tᵤ=31°C")
hlines!(ax4d, [AC_TU], color=:blue, linestyle=:dash, linewidth=1, label="A. chilensis Tᵤ=34°C")
axislegend(ax4d, position=:rb, framevisible=true, labelsize=8)

# Literature annotation: competitive outcome reversal
Label(fig4[0, :],
    "Rochat & Gutierrez (2001): C. utilis dominates coast; A. chilensis dominates inland",
    fontsize=11, color=:gray40, halign=:center)

save(joinpath(figdir, "coastal_vs_inland.png"), fig4, px_per_unit=2)
println("  Saved coastal_vs_inland.png")

# Verify competitive outcome reversal
coast_ratio = mean(cu_coast[yr3]) / max(mean(ac_coast[yr3]), 1e-6)
inland_ratio = mean(ac_inland[yr3]) / max(mean(cu_inland[yr3]), 1e-6)
println("  Coastal C. utilis/A. chilensis ratio: $(round(coast_ratio, digits=2)) (expect >1)")
println("  Inland A. chilensis/C. utilis ratio: $(round(inland_ratio, digits=2)) (expect >1)")

# ============================================================
# FIGURE 5: Temperature Gradient — Parasitoid Dominance
# ============================================================

println("\n── Figure 5: Temperature gradient analysis ──")

T_mean_range = 12.0:0.25:28.0
scale_eq = Float64[]
ac_eq    = Float64[]
cu_eq    = Float64[]

for T_base in T_mean_range
    function regional_weather(day)
        doy = mod(day - 1, 365) + 1
        amplitude = 2.0 + 0.6 * (T_base - 12.0)
        T_mean = T_base + amplitude * sin(2π * (doy - 100) / 365)
        T_max = T_mean + 3.0 + 0.5 * (T_base - 12.0)
        return (T_mean=T_mean, T_max=T_max)
    end

    s, pa, pc = run_two_parasitoids(weather_fn=regional_weather, n_years=8,
                                     S0=50.0, Pa0=5.0, Pc0=5.0)
    # Average over last 3 years to damp limit-cycle oscillations
    yr_last = (365*5+1):(365*8)
    push!(scale_eq, mean(s[yr_last]))
    push!(ac_eq, mean(pa[yr_last]))
    push!(cu_eq, mean(pc[yr_last]))
end

# 3-point moving average to smooth residual oscillations
function smooth3(v)
    n = length(v)
    s = similar(v)
    s[1] = (2v[1] + v[2]) / 3
    s[n] = (v[n-1] + 2v[n]) / 3
    for i in 2:n-1
        s[i] = (v[i-1] + v[i] + v[i+1]) / 3
    end
    return s
end
ac_eq    = smooth3(smooth3(ac_eq))
cu_eq    = smooth3(smooth3(cu_eq))
scale_eq = smooth3(smooth3(scale_eq))

# Find crossover point
T_vec = collect(T_mean_range)
dominance_diff = ac_eq .- cu_eq
crossover_idx = findfirst(i -> dominance_diff[i] > 0 && dominance_diff[max(1, i-1)] <= 0, 2:length(dominance_diff))
crossover_T = crossover_idx !== nothing ? T_vec[crossover_idx + 1] : NaN
println("  Parasitoid crossover temperature: ~$(round(crossover_T, digits=1))°C")

fig5 = Figure(size=(900, 600))

# Panel a: Parasitoid densities
ax5a = Axis(fig5[1, 1],
    title="(a) Parasitoid Dominance Along Temperature Gradient\n(Rochat & Gutierrez 2001; Gutierrez & Pizzamiglio 2007)",
    xlabel="Mean annual temperature (°C)",
    ylabel="Mean parasitoid density (year 6)",
    xlabelsize=13, ylabelsize=13)

lines!(ax5a, T_vec, ac_eq, linewidth=2.5, color=:blue, label="A. chilensis")
lines!(ax5a, T_vec, cu_eq, linewidth=2.5, color=:forestgreen, label="C. utilis")

# Shade dominance zones
vspan!(ax5a, 12.0, 16.0, color=(:green, 0.08))
vspan!(ax5a, 22.0, 28.0, color=(:orange, 0.08))
text!(ax5a, 14.0, maximum(cu_eq) * 0.9,
    text="Coastal\nC. utilis\ndominates", align=(:center, :top), fontsize=9, color=:green)
text!(ax5a, 25.0, maximum(ac_eq) * 0.9,
    text="Inland\nA. chilensis\ndominates", align=(:center, :top), fontsize=9, color=:darkorange)

if !isnan(crossover_T)
    vlines!(ax5a, [crossover_T], color=:gray40, linestyle=:dash, linewidth=1.5)
    text!(ax5a, crossover_T + 0.3, maximum(vcat(ac_eq, cu_eq)) * 0.6,
        text="Crossover\n≈$(round(crossover_T, digits=1))°C", fontsize=9, color=:gray40)
end

axislegend(ax5a, position=:ct, framevisible=true, labelsize=10)

# Panel b: Scale suppression
ax5b = Axis(fig5[2, 1],
    title="(b) Scale Equilibrium Density Along Temperature Gradient",
    xlabel="Mean annual temperature (°C)",
    ylabel="Mean scale density (year 6)",
    xlabelsize=13, ylabelsize=13)

lines!(ax5b, T_vec, scale_eq, linewidth=2.5, color=:red, label="Oleander scale")

# Which parasitoid dominates? Shade and annotate
for i in eachindex(T_vec)
    dominant = ac_eq[i] > cu_eq[i] ? "Ac" : "Cu"
end

vspan!(ax5b, 12.0, 16.0, color=(:green, 0.08))
vspan!(ax5b, 22.0, 28.0, color=(:orange, 0.08))

# Find minimum scale density (best biocontrol)
min_scale_idx = argmin(scale_eq)
min_scale_T = T_vec[min_scale_idx]
scatter!(ax5b, [min_scale_T], [scale_eq[min_scale_idx]], color=:black,
         markersize=12, marker=:star5, label="Best control ($(round(min_scale_T, digits=1))°C)")

axislegend(ax5b, position=:rt, framevisible=true, labelsize=10)

save(joinpath(figdir, "temperature_gradient.png"), fig5, px_per_unit=2)
println("  Saved temperature_gradient.png")
println("  Min scale density: $(round(minimum(scale_eq), digits=1)) at T=$(min_scale_T)°C")
println("  Scale density range: $(round(minimum(scale_eq), digits=1))–$(round(maximum(scale_eq), digits=1))")

# ============================================================
# Summary
# ============================================================

println("\n" * "="^60)
println("✓ All oleander scale validation figures saved to:")
println("  $(figdir)")
println("="^60)
println("\nLiterature verification summary:")
println("  A. nerii crawler→adult at 25°C: $(round(days_to_adult, digits=1)) d (lit: 48–51 d)")
println("  Aphytis dev time at 25°C: $(round(ac_dev_days, digits=1)) d (lit: 13–18 d)")
println("  Scale Brière rate at 30°C: $(round(briere_rate(SCALE_A, SCALE_TL, SCALE_TU, 30.0), digits=4))")
println("  C. utilis Brière rate at 30°C: $(round(briere_rate(CU_A, CU_TL, CU_TU, 30.0), digits=4)) (near upper limit)")
println("  Coastal: C. utilis dominance ratio = $(round(coast_ratio, digits=2))")
println("  Inland: A. chilensis dominance ratio = $(round(inland_ratio, digits=2))")
