#!/usr/bin/env julia
# Generate medfly validation figures matching Fig. 2 from
# Gutierrez & Ponti (2011) — development and mortality rates.

using PhysiologicallyBasedDemographicModels
using CairoMakie

figdir = joinpath(@__DIR__, "figures", "medfly")
mkpath(figdir)

# === Parameters from vignette / Table 1 ===
const T_LOWER_EGG   = 9.5
const T_LOWER_LARVA = 9.5
const T_LOWER_PUPA  = 9.5
const T_LOWER_ADULT = 9.5

const T_UPPER_EGG   = 35.5
const T_UPPER_LARVA = 35.5
const T_UPPER_PUPA  = 34.5
const T_UPPER_ADULT = 34.5

# DD durations
const τ_EGG   = 31.0
const τ_LARVA = 97.0
const τ_PUPA  = 165.0
const τ_ADULT = 673.0

# Mortality coefficients — Eq. 2
const μ_EL_A = 0.0004;  const μ_EL_B = 0.0145;  const μ_EL_C = 0.1314
const μ_PU_A = 0.0005;  const μ_PU_B = 0.0207;  const μ_PU_C = 0.2142
const μ_AD_A = 0.00049; const μ_AD_B = 0.0187;   const μ_AD_C = 0.1846

# Brière development rate models — calibrated to match Fig 2 peaks
# r(T) = a · T · (T − T_lower) · √(T_upper − T) for T_lower < T < T_upper
egg_dev  = BriereDevelopmentRate(0.0004498, T_LOWER_EGG,  T_UPPER_EGG)
pupa_dev = BriereDevelopmentRate(0.0001055, T_LOWER_PUPA, T_UPPER_PUPA)

# Mortality functions
μ_egg_larval(T) = max(0.0, μ_EL_A * T^2 - μ_EL_B * T + μ_EL_C)
μ_pupal(T)      = max(0.0, μ_PU_A * T^2 - μ_PU_B * T + μ_PU_C)
μ_adult(T)      = max(0.0, μ_AD_A * T^2 - μ_AD_B * T + μ_AD_C)

# Temperature range
Ts = range(0.0, 40.0, length=200)

# === Figure 2a: Egg development rate (1/days) ===
# Paper shows nonlinear Brière curve peaking ~0.65 at ~30°C
egg_rate = [development_rate(egg_dev, T) for T in Ts]

fig1 = Figure(size=(800, 600))
ax1 = Axis(fig1[1,1], xlabel="Temperature (°C)", ylabel="1/days",
           title="(a) Eggs — Development Rate", xlabelsize=14, ylabelsize=14)
lines!(ax1, collect(Ts), egg_rate, color=:black, linewidth=2)
xlims!(ax1, 5, 40)
ylims!(ax1, 0, 0.8)
save(joinpath(figdir, "fig2a_egg_devrate.png"), fig1, px_per_unit=2)
println("Saved fig2a — egg peak rate: ", round(maximum(egg_rate), digits=3))

# === Figure 2b: Pupa development rate ===
pupa_rate = [development_rate(pupa_dev, T) for T in Ts]

fig2 = Figure(size=(800, 600))
ax2 = Axis(fig2[1,1], xlabel="Temperature (°C)", ylabel="1/days",
           title="(b) Pupae — Development Rate", xlabelsize=14, ylabelsize=14)
lines!(ax2, collect(Ts), pupa_rate, color=:black, linewidth=2)
xlims!(ax2, 5, 40)
ylims!(ax2, 0, 0.16)
save(joinpath(figdir, "fig2b_pupa_devrate.png"), fig2, px_per_unit=2)
println("Saved fig2b — pupa peak rate: ", round(maximum(pupa_rate), digits=3))

# === Figure 2c: Egg mortality rate ===
egg_mort = [μ_egg_larval(T) for T in Ts]

fig3 = Figure(size=(800, 600))
ax3 = Axis(fig3[1,1], xlabel="Temperature (°C)", ylabel="Mortality rate/day",
           title="(c) Eggs — Mortality Rate", xlabelsize=14, ylabelsize=14)
lines!(ax3, collect(Ts), egg_mort, color=:black, linewidth=2)
xlims!(ax3, 5, 40)
ylims!(ax3, 0, 0.25)
save(joinpath(figdir, "fig2c_egg_mortality.png"), fig3, px_per_unit=2)
T_opt_egg = μ_EL_B / (2 * μ_EL_A)
println("Saved fig2c — egg min mortality at T=", round(T_opt_egg, digits=1), "°C: ",
        round(μ_egg_larval(T_opt_egg), digits=4))

# === Figure 2d: Pupa mortality rate ===
pupa_mort = [μ_pupal(T) for T in Ts]

fig4 = Figure(size=(800, 600))
ax4 = Axis(fig4[1,1], xlabel="Temperature (°C)", ylabel="Mortality rate/day",
           title="(d) Pupae — Mortality Rate", xlabelsize=14, ylabelsize=14)
lines!(ax4, collect(Ts), pupa_mort, color=:black, linewidth=2)
xlims!(ax4, 5, 40)
ylims!(ax4, 0, 0.14)
save(joinpath(figdir, "fig2d_pupa_mortality.png"), fig4, px_per_unit=2)
T_opt_pupa = μ_PU_B / (2 * μ_PU_A)
println("Saved fig2d — pupa min mortality at T=", round(T_opt_pupa, digits=1), "°C: ",
        round(μ_pupal(T_opt_pupa), digits=4))

# === Figure 2e: All mortality rates combined ===
adult_mort = [μ_adult(T) for T in Ts]
Ts_ext = range(-10.0, 40.0, length=250)
egg_mort_ext   = [μ_egg_larval(T) for T in Ts_ext]
pupa_mort_ext  = [μ_pupal(T) for T in Ts_ext]
adult_mort_ext = [μ_adult(T) for T in Ts_ext]

fig5 = Figure(size=(800, 600))
ax5 = Axis(fig5[1,1], xlabel="Temperature (°C)", ylabel="Mortality rate/day",
           title="(e) Adults — Mortality Rates (all stages)", xlabelsize=14, ylabelsize=14)
lines!(ax5, collect(Ts_ext), egg_mort_ext, color=:gray, linewidth=1.5, linestyle=:dash,
       label="eggs-larvae")
lines!(ax5, collect(Ts_ext), pupa_mort_ext, color=:gray, linewidth=1.5, linestyle=:dot,
       label="pupae")
lines!(ax5, collect(Ts_ext), adult_mort_ext, color=:black, linewidth=2, label="adults")
xlims!(ax5, -10, 40)
ylims!(ax5, 0, 0.4)
axislegend(ax5, position=:rt)
save(joinpath(figdir, "fig2e_all_mortality.png"), fig5, px_per_unit=2)
println("Saved fig2e — adult min mortality at T=",
        round(μ_AD_B / (2 * μ_AD_A), digits=1), "°C")

# === Combined panel figure (matching paper layout) ===
fig_all = Figure(size=(1200, 900))

ax_a = Axis(fig_all[1,1], xlabel="temperature (°C)", ylabel="1/ days", title="a  eggs")
lines!(ax_a, collect(Ts), egg_rate, color=:black, linewidth=2)
xlims!(ax_a, 5, 40); ylims!(ax_a, 0, 0.8)

ax_b = Axis(fig_all[1,2], xlabel="temperature (°C)", ylabel="1/ days", title="b  pupae")
lines!(ax_b, collect(Ts), pupa_rate, color=:black, linewidth=2)
xlims!(ax_b, 5, 40); ylims!(ax_b, 0, 0.16)

ax_c = Axis(fig_all[2,1], xlabel="temperature (°C)", ylabel="mortality rate/day", title="c  eggs")
lines!(ax_c, collect(Ts), egg_mort, color=:black, linewidth=2)
xlims!(ax_c, 5, 40); ylims!(ax_c, 0, 0.25)

ax_d = Axis(fig_all[2,2], xlabel="temperature (°C)", ylabel="mortality rate/day", title="d  pupae")
lines!(ax_d, collect(Ts), pupa_mort, color=:black, linewidth=2)
xlims!(ax_d, 5, 40); ylims!(ax_d, 0, 0.14)

ax_e = Axis(fig_all[3,1:2], xlabel="temperature (°C)", ylabel="mortality rate/day",
            title="e  mortality rate/day")
lines!(ax_e, collect(Ts_ext), egg_mort_ext, color=:gray, linewidth=1.5, linestyle=:dash,
       label="eggs-larvae")
lines!(ax_e, collect(Ts_ext), pupa_mort_ext, color=:gray, linewidth=1.5, linestyle=:dot,
       label="pupae")
lines!(ax_e, collect(Ts_ext), adult_mort_ext, color=:black, linewidth=2, label="adults")
xlims!(ax_e, -10, 40); ylims!(ax_e, 0, 0.4)
axislegend(ax_e, position=:rt)

save(joinpath(figdir, "fig2_combined.png"), fig_all, px_per_unit=2)
println("\nSaved combined Fig 2 panel")
println("All medfly validation figures saved to: ", figdir)
