"""
    temperature_responses.jl

Re-usable temperature- and photoperiod-dependent response functions
extracted from vignettes after a duplication audit. These are
deliberately *plain functions* (not types) — they capture the simple
algebraic forms that recur across many published PBDM applications and
are convenient when scripting illustrative examples or fitting models.

For the type-based variants used inside the solver pipeline, see
`BriereDevelopmentRate`, `LoganDevelopmentRate`, `Q10Respiration`,
`HollingTypeII`, `HollingTypeIII` and `FraserGilbertResponse`.
"""

# -- Triangular thermal scalar --------------------------------------------

"""
    triangular_thermal_scalar(T; θL, θU, Topt = (θL + θU)/2)

Triangular thermal performance scalar in `[0, 1]`, peaking at `Topt`
and zero outside `[θL, θU]`. ASCII alias: `phiT`.

```math
\\phi(T) = \\begin{cases}
0 & T < \\theta_L \\text{ or } T > \\theta_U \\\\
(T-\\theta_L)/(T_{opt}-\\theta_L) & \\theta_L \\le T \\le T_{opt} \\\\
(\\theta_U-T)/(\\theta_U-T_{opt}) & T_{opt} < T \\le \\theta_U
\\end{cases}
```
"""
function triangular_thermal_scalar(T::Real; θL::Real, θU::Real,
                                   Topt::Real = (θL + θU) / 2)
    (T < θL || T > θU) && return zero(float(T))
    if T ≤ Topt
        return (T - θL) / (Topt - θL)
    else
        return (θU - T) / (θU - Topt)
    end
end

const phiT = triangular_thermal_scalar

# -- Briere fecundity / rate ----------------------------------------------

"""
    briere_rate(T; a, θL, θU, m = 2)

Briere-1 form ``a \\cdot T (T-\\theta_L)(\\theta_U-T)^{1/m}``,
zero outside `[θL, θU]`. The exponent `m` defaults to 2 (square root,
the original Brière 1999 form). Used for both development rates and
fecundity functions throughout the corpus.
"""
function briere_rate(T::Real; a::Real, θL::Real, θU::Real, m::Real = 2)
    (T ≤ θL || T ≥ θU) && return zero(float(T))
    return a * T * (T - θL) * (θU - T)^(1 / m)
end

"""
    fecundity_briere(T; a, θL, θU, m = 2)

Alias for `briere_rate` with semantics of "eggs per female per day"
when applied to fecundity data; provided for readability in
fecundity-focussed code.
"""
fecundity_briere(T::Real; kwargs...) = briere_rate(T; kwargs...)

"""
    fecundity_gaussian(T; F_max, T_opt, T_range)

Gaussian fecundity ``F_{max} \\cdot \\exp\\!\\big(-((T-T_{opt})/T_{range})^2\\big)``.
"""
function fecundity_gaussian(T::Real; F_max::Real, T_opt::Real, T_range::Real)
    return F_max * exp(-((T - T_opt) / T_range)^2)
end

# -- Daily mortality (quadratic) ------------------------------------------

"""
    daily_mortality_quadratic(T; a, T_opt, μ_min = 0.0, μ_max = 1.0)

Quadratic-in-temperature daily mortality
``\\mu(T) = a (T - T_{opt})^2 + \\mu_\\text{min}``,
clamped to `[0, μ_max]`. Provides a convenient "U-shape"
mortality envelope around an optimum temperature.
"""
function daily_mortality_quadratic(T::Real; a::Real, T_opt::Real,
                                   μ_min::Real = 0.0, μ_max::Real = 1.0)
    return clamp(a * (T - T_opt)^2 + μ_min, 0.0, μ_max)
end

# -- Diapause induction ---------------------------------------------------

"""
    diapause_fraction_logistic(D; D50, slope)

Logistic photoperiod-driven diapause induction
``1/(1+\\exp(s\\,(D - D_{50})))`` with `D` the day length (hours),
`D50` the inflection day length, and `slope` the steepness.
"""
function diapause_fraction_logistic(D::Real; D50::Real, slope::Real)
    return 1 / (1 + exp(slope * (D - D50)))
end

"""
    diapause_fraction_linear(D, T; D_crit, D_comp, T_mod = -Inf, k = 0.0)

Linear day-length ramp from `D_comp` (full diapause) up to `D_crit`
(no diapause), with optional **temperature modulation** below `T_mod`
that scales the photoperiod fraction by `(1 + k*(T_mod - T))`. The
returned value is clamped to `[0, 1]`.
"""
function diapause_fraction_linear(D::Real, T::Real;
                                  D_crit::Real, D_comp::Real,
                                  T_mod::Real = -Inf, k::Real = 0.0)
    photo = if D ≥ D_crit
        zero(float(D))
    elseif D ≤ D_comp
        one(float(D))
    else
        (D_crit - D) / (D_crit - D_comp)
    end
    if T < T_mod
        photo *= (1 + k * (T_mod - T))
    end
    return clamp(photo, 0.0, 1.0)
end

# -- Gilbert-Fraser ratio-dependent attack --------------------------------

"""
    gilbert_fraser_attack(α, D, H)

Per-period proportion of hosts attacked under the Gilbert–Fraser
ratio-dependent functional response:
``N_a / H = 1 - \\exp(-\\alpha\\, D / H)``,
where `D` is parasitoid demand (e.g. egg load × females × window) and
`H` is host density. Returns 0 when `H ≤ 0`. The parameter `α`
("search rate") is dimensionless: `α = 1` random search, `α > 1`
above-random (chemo-cued).
"""
function gilbert_fraser_attack(α::Real, D::Real, H::Real)
    H ≤ 0 && return zero(float(H))
    return 1 - exp(-α * D / max(H, eps(float(H))))
end
