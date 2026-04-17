"""
Analytical tools for PBDM theory.

Implements the analytical framework from the foundational PBDM papers:
- Gutierrez (1994) "A Physiologically Based Tritrophic Perspective"
- Schreiber & Gutierrez (1998) "A Supply/Demand Perspective"
- Wang & Gutierrez (1980) "Stability Analyses in Population Ecology"
- Gutierrez & Regev (2005) "Bioeconomics of Tritrophic Systems"

Provides: functional response traits, metabolic compensation point,
zero-growth isocline analysis, equilibrium classification, food-web
assembly diagrams, life-history strategy classification, and continuous
ODE generation for PBDM trophic systems.
"""

# ================================================================
# Functional response classification
# ================================================================

"""
    is_ratio_dependent(fr::AbstractFunctionalResponse) -> Bool

Return `true` if the functional response is ratio-dependent (f(S/D)),
`false` if prey-dependent (f(N)).

PBDM theory (Gutierrez 1994) emphasizes ratio-dependent responses
as the biologically correct general form for resource acquisition.
The Fraser-Gilbert response is ratio-dependent; Holling Types II/III
are prey-dependent.
"""
is_ratio_dependent(::FraserGilbertResponse) = true
is_ratio_dependent(::HollingTypeII) = false
is_ratio_dependent(::HollingTypeIII) = false

"""
    apparency(fr::FraserGilbertResponse) -> Real

Return the apparency parameter `a` from a Fraser-Gilbert functional response.

In PBDM theory (Gutierrez 1994, Eq. 5), apparency represents the
effective proportion of a trophic level's biomass available to consumers.
Higher apparency increases acquisition for a given supply/demand ratio.
The parameter controls isocline shape: when `a` < intrinsic growth rate,
the resource isocline is hump-shaped, producing stable coexistence.

See also: [`FraserGilbertResponse`](@ref), [`consumer_isocline`](@ref)
"""
apparency(fr::FraserGilbertResponse) = fr.a

# ================================================================
# Metabolic Compensation Point (MCP)
# ================================================================

"""
    compensation_point(resp::Q10Respiration, T::Real, demand_rate::Real;
                       conversion_efficiency::Real=1.0) -> Real

Compute the metabolic compensation point (MCP) φ* at temperature `T`.

The MCP is the minimum supply/demand ratio required for net positive
growth — the point where resource acquisition exactly balances respiratory
losses (Schreiber & Gutierrez 1998):

    φ* = R(T) / (ε × d)

where `R(T)` is the respiration rate, `d` is the per-capita demand rate,
and `ε` is the conversion efficiency of consumed resources to consumer
biomass.

# Returns
- `φ*` — the compensation point (dimensionless)
- Values ≥ 1.0 indicate the population cannot persist at temperature `T`
- Values < 0 should not occur with valid parameters

# Arguments
- `resp`: Respiration model
- `T`: Temperature (°C)
- `demand_rate`: Per-capita demand rate (`d`)
- `conversion_efficiency`: Fraction of consumed resource converted to
  consumer biomass (default 1.0, i.e., demand already in assimilation units)

See also: [`consumer_isocline`](@ref), [`food_web_assembly`](@ref)
"""
function compensation_point(resp::Q10Respiration, T::Real, demand_rate::Real;
                            conversion_efficiency::Real=1.0)
    demand_rate > 0 || throw(ArgumentError("demand_rate must be positive"))
    conversion_efficiency > 0 || throw(ArgumentError("conversion_efficiency must be positive"))
    return respiration_rate(resp, T) / (conversion_efficiency * demand_rate)
end

"""
    compensation_point(bdf::BiodemographicFunctions, T::Real;
                       demand_rate::Real, conversion_efficiency::Real=1.0)

Compute the MCP from a BDF composite's respiration component.
"""
function compensation_point(bdf::BiodemographicFunctions, T::Real;
                            demand_rate::Real, conversion_efficiency::Real=1.0)
    return compensation_point(respiration_component(bdf), T, demand_rate;
                              conversion_efficiency=conversion_efficiency)
end

# ================================================================
# Life-history strategy classification
# ================================================================

"""
    life_history_strategy(pool::MetabolicPool) -> Symbol

Classify life-history strategy from metabolic pool allocation priorities.

From Gutierrez & Regev (2005), the r–K continuum maps to allocation
priority ordering:
- `:r_selected` — growth prioritized before reproduction (colonizers,
  weeds, high intrinsic rate)
- `:K_selected` — reproduction prioritized before growth (competitors,
  stable environments, investment in offspring quality)
- `:unknown` — allocation labels don't include both `:growth` and
  `:reproduction`

The classification examines which of `:growth` or `:reproduction` appears
first in the priority cascade of the metabolic pool.
"""
function life_history_strategy(pool::MetabolicPool)
    labels = pool.labels
    growth_idx = findfirst(==(:growth), labels)
    repro_idx = findfirst(==(:reproduction), labels)
    (growth_idx === nothing || repro_idx === nothing) && return :unknown
    return growth_idx < repro_idx ? :r_selected : :K_selected
end

function life_history_strategy(model::CoupledPBDMModel)
    return life_history_strategy(allocation_model(model))
end

# ================================================================
# Zero-growth isocline analysis
# ================================================================

"""
    IsoclineResult{T<:Real}

Result of a zero-growth isocline computation in (M₁, M₂) phase space.

# Fields
- `resource_biomass::Vector{T}` — M₁ values (resource axis)
- `consumer_biomass::Vector{T}` — M₂ values (consumer axis)
- `isocline_type::Symbol` — `:consumer` or `:resource`
- `slope::Union{T, Nothing}` — slope if the isocline is linear
  (consumer isocline with Fraser-Gilbert is always a line through the origin)
"""
struct IsoclineResult{T<:Real}
    resource_biomass::Vector{T}
    consumer_biomass::Vector{T}
    isocline_type::Symbol
    slope::Union{T, Nothing}
end

"""
    consumer_isocline(fr::FraserGilbertResponse, resp::Q10Respiration, T::Real;
                      demand_rate::Real, conversion_efficiency::Real=1.0,
                      M1_range::AbstractVector) -> IsoclineResult

Compute the consumer zero-growth isocline (dM₂/dt = 0).

For the Fraser-Gilbert (ratio-dependent) response, the consumer isocline
is a straight line through the origin (Gutierrez 1994, Fig. 3):

    M₂ = [a / (d × ln(1/(1 − φ*)))] × M₁

where φ* is the metabolic compensation point and `a` is the apparency
parameter. The isocline slope depends only on the consumer's physiology,
not on its biomass.

Returns an `IsoclineResult` with `NaN` values if φ* ≥ 1 (population
cannot persist at this temperature).
"""
function consumer_isocline(fr::FraserGilbertResponse, resp::Q10Respiration, T::Real;
                           demand_rate::Real, conversion_efficiency::Real=1.0,
                           M1_range::AbstractVector)
    φ_star = compensation_point(resp, T, demand_rate;
                                conversion_efficiency=conversion_efficiency)
    M1 = collect(Float64, M1_range)
    if φ_star >= 1.0
        return IsoclineResult(M1, fill(NaN, length(M1)), :consumer, nothing)
    end
    s = fr.a / (demand_rate * log(1 / (1 - φ_star)))
    M2 = s .* M1
    return IsoclineResult(M1, M2, :consumer, s)
end

"""
    resource_isocline(fr::FraserGilbertResponse;
                      intrinsic_rate::Real, carrying_capacity::Real,
                      consumer_demand_rate::Real,
                      M1_range::AbstractVector) -> IsoclineResult

Compute the resource zero-growth isocline (dM₁/dt = 0).

For logistic resource growth with Fraser-Gilbert consumer predation,
the isocline is the classic hump-shaped curve:

    r × M₁ × (1 − M₁/K) = d₂ × M₂ × (1 − exp(−a × M₁/(d₂ × M₂)))

The isocline reaches its peak at intermediate M₁ and returns to zero at
M₁ = 0 and M₁ = K. Points above the hump represent net resource decline;
points below represent net resource growth.

Solves implicitly via bisection for each M₁ in `M1_range`.
"""
function resource_isocline(fr::FraserGilbertResponse;
                           intrinsic_rate::Real, carrying_capacity::Real,
                           consumer_demand_rate::Real,
                           M1_range::AbstractVector)
    M1_out = Float64[]
    M2_out = Float64[]

    for M1 in M1_range
        M1 <= 0.0 && continue
        growth = intrinsic_rate * M1 * (1 - M1 / carrying_capacity)
        growth <= 0.0 && continue

        # Target: d₂*M₂*(1 - exp(-a*M₁/(d₂*M₂))) = growth
        # RHS is monotonically increasing in M₂ from 0 to a*M₁
        a_M1 = fr.a * M1
        if growth > a_M1
            continue  # Consumers can never eat enough
        end

        # Bisection
        lo = 1e-12
        hi = max(1e8, 100.0 * carrying_capacity)
        for _ in 1:200
            mid = (lo + hi) / 2
            D = consumer_demand_rate * mid
            consumption = D * (1 - exp(-a_M1 / D))
            if consumption < growth
                lo = mid
            else
                hi = mid
            end
            (hi - lo) < 1e-10 * hi && break
        end
        push!(M1_out, Float64(M1))
        push!(M2_out, (lo + hi) / 2)
    end
    return IsoclineResult(M1_out, M2_out, :resource, nothing)
end

# ================================================================
# Equilibrium analysis
# ================================================================

"""
    EquilibriumResult{T<:Real}

Result of a trophic equilibrium analysis.

# Fields
- `M1_star::T` — resource equilibrium biomass
- `M2_star::T` — consumer equilibrium biomass
- `eigenvalues::Vector{Complex{T}}` — eigenvalues of the Jacobian
- `classification::Symbol` — stability classification:
  `:stable_node`, `:stable_focus`, `:unstable_node`, `:unstable_focus`,
  `:saddle`, `:center`, or `:degenerate`
- `jacobian::Matrix{T}` — the 2×2 Jacobian at equilibrium
"""
struct EquilibriumResult{T<:Real}
    M1_star::T
    M2_star::T
    eigenvalues::Vector{Complex{T}}
    classification::Symbol
    jacobian::Matrix{T}
end

"""
    classify_equilibrium(eigenvalues) -> Symbol

Classify an equilibrium point from its Jacobian eigenvalues (Wang &
Gutierrez 1980).

Returns one of:
- `:stable_node` — both eigenvalues real and negative
- `:stable_focus` — complex conjugate pair with negative real parts
- `:unstable_node` — both eigenvalues real and positive
- `:unstable_focus` — complex conjugate pair with positive real parts
- `:saddle` — real eigenvalues of opposite sign
- `:center` — purely imaginary eigenvalues
- `:degenerate` — at least one zero eigenvalue
"""
function classify_equilibrium(eigenvalues::AbstractVector{<:Complex})
    tol = 1e-10
    re = real.(eigenvalues)
    im_parts = imag.(eigenvalues)

    any(abs.(re) .< tol .&& abs.(im_parts) .< tol) && return :degenerate

    has_imaginary = any(abs.(im_parts) .> tol)

    if has_imaginary
        all(abs.(re) .< tol) && return :center
        all(re .< -tol) && return :stable_focus
        all(re .> tol) && return :unstable_focus
        return :degenerate
    else
        all(re .< -tol) && return :stable_node
        all(re .> tol) && return :unstable_node
        re[1] * re[2] < 0 && return :saddle
        return :degenerate
    end
end

"""
    find_equilibrium(fr::FraserGilbertResponse, resp::Q10Respiration, T::Real;
                     intrinsic_rate::Real, carrying_capacity::Real,
                     consumer_demand_rate::Real,
                     conversion_efficiency::Real=1.0) -> EquilibriumResult

Find and classify the interior equilibrium for a resource–consumer
system with Fraser-Gilbert functional response.

Uses the analytical solution for the intersection of the consumer
isocline (line through origin) and the resource isocline (logistic hump):

    M₁* = K × (1 − a × φ* / (r × ln(1/(1−φ*))))
    M₂* = slope × M₁*

Then computes the Jacobian analytically and classifies via eigenvalues
(Wang & Gutierrez 1980).

Returns `EquilibriumResult` with `NaN` values if no interior equilibrium
exists (MCP ≥ 1 or M₁* ≤ 0).
"""
function find_equilibrium(fr::FraserGilbertResponse, resp::Q10Respiration, T::Real;
                          intrinsic_rate::Real, carrying_capacity::Real,
                          consumer_demand_rate::Real,
                          conversion_efficiency::Real=1.0)
    r = Float64(intrinsic_rate)
    K = Float64(carrying_capacity)
    d₂ = Float64(consumer_demand_rate)
    a = Float64(fr.a)
    ε = Float64(conversion_efficiency)

    φ_star = compensation_point(resp, T, d₂; conversion_efficiency=ε)

    if φ_star >= 1.0
        nan_J = fill(NaN, 2, 2)
        return EquilibriumResult(NaN, NaN, Complex{Float64}[], :degenerate, nan_J)
    end

    log_term = log(1 / (1 - φ_star))  # ln(1/(1-φ*))
    s = a / (d₂ * log_term)           # consumer isocline slope

    # Equilibrium: consumer isocline M₂ = s*M₁ substituted into resource eq
    # r*(1 - M₁/K) = d₂*s*φ*
    M1_star = K * (1 - d₂ * s * φ_star / r)
    M2_star = s * M1_star

    if M1_star <= 0.0 || M2_star <= 0.0
        nan_J = fill(NaN, 2, 2)
        return EquilibriumResult(NaN, NaN, Complex{Float64}[], :degenerate, nan_J)
    end

    # Jacobian at equilibrium (analytical derivatives)
    # u = a*M₁/(d₂*M₂) = a/(d₂*s) = log_term at equilibrium
    exp_neg_u = 1 - φ_star

    # F₁ = r*M₁*(1 - M₁/K) - d₂*M₂*(1 - exp(-u))
    J11 = r * (1 - 2 * M1_star / K) - a * exp_neg_u
    J12 = -(d₂ * φ_star - d₂ * log_term * exp_neg_u)

    # F₂ = (ε*d₂*(1-exp(-u)) - R₂)*M₂
    J21 = ε * a * exp_neg_u
    J22 = -ε * d₂ * log_term * exp_neg_u  # simplified since ε*d₂*φ*=R₂ at equilibrium

    J = [J11 J12; J21 J22]
    eigs = eigvals(J)
    class = classify_equilibrium(complex.(eigs))

    return EquilibriumResult(M1_star, M2_star, complex.(eigs), class, J)
end

# ================================================================
# Food-web assembly
# ================================================================

"""
    SpeciesProfile{T<:Real}

Parameters for a species in a food-web assembly analysis.

# Fields
- `name::Symbol` — species identifier
- `demand_rate::T` — per-capita demand rate
- `resp::Q10Respiration{T}` — respiration model
- `fr::AbstractFunctionalResponse` — functional response (for resource acquisition)
- `conversion_efficiency::T` — assimilation efficiency
"""
struct SpeciesProfile{T<:Real}
    name::Symbol
    demand_rate::T
    resp::Q10Respiration{T}
    fr::AbstractFunctionalResponse
    conversion_efficiency::T
end

function SpeciesProfile(name::Symbol; demand_rate::Real, resp::Q10Respiration,
                        fr::AbstractFunctionalResponse, conversion_efficiency::Real=1.0)
    T = Float64
    SpeciesProfile{T}(name, T(demand_rate), resp, fr, T(conversion_efficiency))
end

"""
    AssemblyResult{T<:Real}

Result of a food-web assembly analysis (Schreiber & Gutierrez 1998).

# Fields
- `species::Vector{Symbol}` — species names in invasion order (lowest MCP first)
- `mcps::Vector{T}` — metabolic compensation points
- `invasion_order::Vector{Int}` — indices into original species list
- `can_persist::BitVector` — whether each species can persist (MCP < 1)
"""
struct AssemblyResult{T<:Real}
    species::Vector{Symbol}
    mcps::Vector{T}
    invasion_order::Vector{Int}
    can_persist::BitVector
end

"""
    food_web_assembly(profiles::Vector{<:SpeciesProfile}, T_celsius::Real) -> AssemblyResult

Compute food-web assembly order from metabolic compensation points
(Schreiber & Gutierrez 1998).

Species with lower MCPs can invade systems where species with higher
MCPs are established. The assembly order (lowest MCP first) predicts
the sequence of successful invasions in a community assembly process.

Species with MCP ≥ 1 cannot persist at the given temperature.

# Arguments
- `profiles`: Vector of `SpeciesProfile` for each candidate species
- `T_celsius`: Temperature (°C)

# Returns
- `AssemblyResult` with species ranked by MCP
"""
function food_web_assembly(profiles::Vector{<:SpeciesProfile}, T_celsius::Real)
    n = length(profiles)
    mcps = Float64[
        compensation_point(sp.resp, T_celsius, sp.demand_rate;
                           conversion_efficiency=sp.conversion_efficiency)
        for sp in profiles
    ]
    order = sortperm(mcps)
    can_persist = mcps .< 1.0

    return AssemblyResult{Float64}(
        [profiles[i].name for i in order],
        mcps[order],
        order,
        can_persist[order]
    )
end


