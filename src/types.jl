"""
Core types for PhysiologicallyBasedDemographicModels.jl

Trait types for dispatch, abstract supertypes, and concrete types
for PBDM model components.
"""

# --- Approach families (MP vs BDF) ---

"""
    AbstractPBDMApproach

High-level supertype for PBDM modeling approaches.
"""
abstract type AbstractPBDMApproach end

"""Approach family for biodemographic functions (BDF)."""
abstract type AbstractBiodemographicApproach <: AbstractPBDMApproach end

"""Approach family for metabolic-pool allocation (MP)."""
abstract type AbstractAllocationApproach <: AbstractPBDMApproach end

"""Approach family for explicit hybrid MP+BDF compositions."""
abstract type AbstractHybridPBDMApproach <: AbstractPBDMApproach end

"""Primitive biodemographic function component."""
abstract type AbstractBiodemographicFunction <: AbstractBiodemographicApproach end

"""Composite biodemographic model."""
abstract type AbstractBiodemographicModel <: AbstractBiodemographicApproach end

"""Allocation model used to partition acquired resources."""
abstract type AbstractAllocationModel <: AbstractAllocationApproach end

"""Temperature- or state-dependent respiration submodel."""
abstract type AbstractRespirationModel <: AbstractBiodemographicFunction end

"""
    approach_family(x)

Return `:bdf`, `:mp`, or `:hybrid` for a package component or composite model.
"""
approach_family(::AbstractBiodemographicApproach) = :bdf
approach_family(::AbstractAllocationApproach) = :mp
approach_family(::AbstractHybridPBDMApproach) = :hybrid

# --- Structure traits (PBDM-specific) ---

abstract type AbstractPBDMStructure <: AbstractProjectionStructure end

"""
    SingleSpeciesPBDM

A single-species PBDM with no trophic interactions.
Population dynamics driven by physiological time and resource allocation.
"""
struct SingleSpeciesPBDM <: AbstractPBDMStructure end

"""
    MultiSpeciesPBDM

Multi-species PBDM structure tag.

This identifies problems that contain multiple populations, but direct
`solve(::PBDMProblem, ::DirectIteration)` support is not yet implemented.
"""
struct MultiSpeciesPBDM <: AbstractPBDMStructure end

# --- Development rate models ---

"""
    AbstractDevelopmentRate

Supertype for temperature-dependent development rate functions.
Development rate maps temperature to daily physiological progress
or daily degree-day accumulation.
"""
abstract type AbstractDevelopmentRate <: AbstractBiodemographicFunction end

"""
    LinearDevelopmentRate(T_lower, T_upper)

Linear degree-day accumulation:
`r(T) = max(0, min(T_upper - T_lower, T - T_lower))`.

This returns daily degree-days above `T_lower`, capped at the interval
`T_upper - T_lower`.
"""
struct LinearDevelopmentRate{T<:Real} <: AbstractDevelopmentRate
    T_lower::T
    T_upper::T

    function LinearDevelopmentRate(T_lower::T, T_upper::T) where {T<:Real}
        T_lower < T_upper || throw(ArgumentError("T_lower must be < T_upper"))
        new{T}(T_lower, T_upper)
    end
end

function LinearDevelopmentRate(T_lower::Real, T_upper::Real)
    T = promote_type(typeof(T_lower), typeof(T_upper))
    LinearDevelopmentRate(T(T_lower), T(T_upper))
end

"""
    BriereDevelopmentRate(a, T_lower, T_upper)

Brière nonlinear development rate:
`r(T) = a * T * (T - T_lower) * sqrt(max(0, T_upper - T))` for T_lower ≤ T ≤ T_upper.
"""
struct BriereDevelopmentRate{T<:Real} <: AbstractDevelopmentRate
    a::T
    T_lower::T
    T_upper::T

    function BriereDevelopmentRate(a::T, T_lower::T, T_upper::T) where {T<:Real}
        T_lower < T_upper || throw(ArgumentError("T_lower must be < T_upper"))
        a > 0 || throw(ArgumentError("a must be positive"))
        new{T}(a, T_lower, T_upper)
    end
end

function BriereDevelopmentRate(a::Real, T_lower::Real, T_upper::Real)
    T = promote_type(typeof(a), typeof(T_lower), typeof(T_upper))
    BriereDevelopmentRate(T(a), T(T_lower), T(T_upper))
end

"""
    LoganDevelopmentRate(ψ, ρ, T_upper, ΔT)

Logan Type III development rate:
`r(T) = ψ * (exp(ρ * T) - exp(ρ * T_upper - (T_upper - T) / ΔT))`.
"""
struct LoganDevelopmentRate{T<:Real} <: AbstractDevelopmentRate
    ψ::T
    ρ::T
    T_upper::T
    ΔT::T

    function LoganDevelopmentRate(ψ::T, ρ::T, T_upper::T, ΔT::T) where {T<:Real}
        ψ > 0 || throw(ArgumentError("ψ must be positive"))
        ρ > 0 || throw(ArgumentError("ρ must be positive"))
        ΔT > 0 || throw(ArgumentError("ΔT must be positive"))
        new{T}(ψ, ρ, T_upper, ΔT)
    end
end

function LoganDevelopmentRate(ψ::Real, ρ::Real, T_upper::Real, ΔT::Real)
    T = promote_type(typeof(ψ), typeof(ρ), typeof(T_upper), typeof(ΔT))
    LoganDevelopmentRate(T(ψ), T(ρ), T(T_upper), T(ΔT))
end

"""
    development_rate(model, T)

Compute instantaneous development rate at temperature `T`.
"""
function development_rate end

function development_rate(m::LinearDevelopmentRate, T::Real)
    T <= m.T_lower && return zero(T)
    T >= m.T_upper && return (m.T_upper - m.T_lower)
    return T - m.T_lower
end

function development_rate(m::BriereDevelopmentRate, T::Real)
    (T <= m.T_lower || T >= m.T_upper) && return zero(T)
    return m.a * T * (T - m.T_lower) * sqrt(m.T_upper - T)
end

function development_rate(m::LoganDevelopmentRate, T::Real)
    return m.ψ * (exp(m.ρ * T) - exp(m.ρ * m.T_upper - (m.T_upper - T) / m.ΔT))
end

"""
    degree_days(model, T)

Compute degree-day accumulation for one day at temperature `T`.
For `LinearDevelopmentRate`, this is `max(0, T - T_lower)`.
"""
function degree_days(m::LinearDevelopmentRate, T::Real)
    return max(zero(T), T - m.T_lower)
end

function degree_days(m::AbstractDevelopmentRate, T::Real)
    return development_rate(m, T)
end

# --- Distributed Delay (Erlang k-substage) ---

"""
    DistributedDelay{T<:Real}

Manetsch/Vansickle k-substage distributed delay for a single life stage.
Approximates the distribution of developmental times with an Erlang distribution.

The delay has `k` substages, mean developmental time `τ` (in physiological time units),
and variance `σ² = τ²/k`. The flow rate from substage `i` is `r_i = (k/τ) * W_i`.

# Fields
- `k::Int`: Number of substages (higher k → tighter distribution)
- `τ::T`: Mean developmental time in degree-days
- `W::Vector{T}`: Substage state vector (mass, numbers, or other quantity)
"""
struct DistributedDelay{T<:Real}
    k::Int
    τ::T
    W::Vector{T}

    function DistributedDelay(k::Int, τ::T, W::Vector{T}) where {T<:Real}
        k > 0 || throw(ArgumentError("k must be positive"))
        τ > 0 || throw(ArgumentError("τ must be positive"))
        length(W) == k || throw(DimensionMismatch("W must have length k=$k"))
        new{T}(k, τ, W)
    end
end

function DistributedDelay(k::Int, τ::Real; W0::Real=zero(τ))
    T = typeof(float(τ))
    W = fill(T(W0), k)
    DistributedDelay(k, T(τ), W)
end

"""Variance of developmental times."""
delay_variance(d::DistributedDelay) = d.τ^2 / d.k

"""Per-substage flow rate coefficient."""
delay_rate(d::DistributedDelay) = d.k / d.τ

"""Total content across all substages."""
delay_total(d::DistributedDelay) = sum(d.W)

# --- Functional Response (Supply/Demand) ---

"""
    AbstractFunctionalResponse

Supertype for supply/demand functional response models.
"""
abstract type AbstractFunctionalResponse <: AbstractBiodemographicFunction end

"""
    FraserGilbertResponse(a)

Frazer-Gilbert demand-driven functional response (ratio-dependent):
`acquisition = demand * (1 - exp(-a * supply / demand))`

Used for photosynthesis, nutrient uptake, and herbivore consumption.
Parameter `a` is the apparency (search/interception rate): the effective
proportion of the resource trophic level's biomass available to consumers
(Gutierrez 1994, Eq. 5). Higher `a` increases acquisition for a given
supply/demand ratio and controls isocline shape in phase-plane analysis.

See also: [`apparency`](@ref), [`acquire`](@ref), [`is_ratio_dependent`](@ref)
"""
struct FraserGilbertResponse{T<:Real} <: AbstractFunctionalResponse
    a::T

    function FraserGilbertResponse(a::T) where {T<:Real}
        a > 0 || throw(ArgumentError("a must be positive"))
        new{T}(a)
    end
end

"""
    acquire(response, supply, demand)

Compute actual acquisition given available supply and demand.
Returns a value in `[0, demand]`.
"""
function acquire(r::FraserGilbertResponse, supply::Real, demand::Real)
    demand <= 0 && return zero(supply)
    supply <= 0 && return zero(demand)
    return demand * (1 - exp(-r.a * supply / demand))
end

"""
    supply_demand_ratio(response, supply, demand)

Compute the supply/demand index φ ∈ [0, 1].
"""
function supply_demand_ratio(r::FraserGilbertResponse, supply::Real, demand::Real)
    demand <= 0 && return one(supply)
    return acquire(r, supply, demand) / demand
end

# --- Metabolic Pool ---

"""
    MetabolicPool{T<:Real}

Metabolic pool representing resource allocation with priority-based demand.

# Fields
- `supply::T`: Total resource supply (e.g., photosynthate production)
- `demands::Vector{T}`: Demands in priority order (highest priority first)
- `labels::Vector{Symbol}`: Labels for each demand component
"""
struct MetabolicPool{T<:Real} <: AbstractAllocationModel
    supply::T
    demands::Vector{T}
    labels::Vector{Symbol}

    function MetabolicPool(supply::T, demands::Vector{T}, labels::Vector{Symbol}) where {T<:Real}
        length(demands) == length(labels) ||
            throw(DimensionMismatch("demands and labels must have same length"))
        new{T}(supply, demands, labels)
    end
end

function MetabolicPool(supply::Real, demands::Vector{<:Real}, labels::Vector{Symbol})
    T = promote_type(typeof(supply), eltype(demands))
    MetabolicPool(T(supply), T.(demands), labels)
end

"""
    allocate(pool::MetabolicPool)

Allocate supply to demands in priority order. Returns a vector of actual
allocations (same length as demands). Higher-priority demands are filled first;
remaining supply cascades to lower-priority demands.
"""
function allocate(pool::MetabolicPool{T}) where {T}
    n = length(pool.demands)
    alloc = zeros(T, n)
    remaining = pool.supply
    for i in 1:n
        alloc[i] = min(remaining, pool.demands[i])
        remaining -= alloc[i]
        remaining <= 0 && break
    end
    return alloc
end

"""
    supply_demand_index(pool::MetabolicPool)

Compute overall supply/demand ratio φ ∈ [0, 1].
"""
function supply_demand_index(pool::MetabolicPool)
    total_demand = sum(pool.demands)
    total_demand <= 0 && return one(pool.supply)
    return min(one(pool.supply), pool.supply / total_demand)
end

# --- Explicit MP/BDF compositions ---

"""
    BiodemographicFunctions(development, acquisition, respiration; label=:bdf)

Bundle a development-rate model, acquisition/functional-response model, and
respiration model into an explicit BDF layer.
"""
struct BiodemographicFunctions{
    D<:AbstractDevelopmentRate,
    F<:AbstractFunctionalResponse,
    R<:AbstractRespirationModel
} <: AbstractBiodemographicModel
    development::D
    acquisition::F
    respiration::R
    label::Symbol
end

function BiodemographicFunctions(development::D, acquisition::F, respiration::R;
                                 label::Symbol=:bdf) where {
                                     D<:AbstractDevelopmentRate,
                                     F<:AbstractFunctionalResponse,
                                     R<:AbstractRespirationModel
                                 }
    BiodemographicFunctions{D,F,R}(development, acquisition, respiration, label)
end

development_component(model::BiodemographicFunctions) = model.development
acquisition_component(model::BiodemographicFunctions) = model.acquisition
respiration_component(model::BiodemographicFunctions) = model.respiration

"""
    CoupledPBDMModel(bdf, allocation; label=:hybrid)

Explicit composition of a BDF layer and a metabolic-pool allocation layer.
"""
struct CoupledPBDMModel{
    B<:AbstractBiodemographicModel,
    A<:AbstractAllocationModel
} <: AbstractHybridPBDMApproach
    bdf::B
    allocation::A
    label::Symbol
end

function CoupledPBDMModel(bdf::B, allocation::A;
                          label::Symbol=:hybrid) where {
                              B<:AbstractBiodemographicModel,
                              A<:AbstractAllocationModel
                          }
    CoupledPBDMModel{B,A}(bdf, allocation, label)
end

biodemography(model::CoupledPBDMModel) = model.bdf
allocation_model(model::CoupledPBDMModel) = model.allocation

# --- Respiration ---

"""
    Q10Respiration{T<:Real}

Temperature-dependent respiration rate using Q₁₀ scaling.
`R(T) = R_ref * Q₁₀^((T - T_ref) / 10)`

# Fields
- `R_ref::T`: Reference respiration rate (fraction of dry mass per day)
- `Q10::T`: Q₁₀ factor (typically 2.0–2.5)
- `T_ref::T`: Reference temperature (°C)
"""
struct Q10Respiration{T<:Real} <: AbstractRespirationModel
    R_ref::T
    Q10::T
    T_ref::T

    function Q10Respiration(R_ref::T, Q10::T, T_ref::T) where {T<:Real}
        R_ref >= 0 || throw(ArgumentError("R_ref must be non-negative"))
        Q10 > 0 || throw(ArgumentError("Q10 must be positive"))
        new{T}(R_ref, Q10, T_ref)
    end
end

function Q10Respiration(R_ref::Real, Q10::Real, T_ref::Real)
    T = promote_type(typeof(R_ref), typeof(Q10), typeof(T_ref))
    Q10Respiration(T(R_ref), T(Q10), T(T_ref))
end

"""
    respiration_rate(model, T)

Compute respiration rate at temperature `T`.
"""
function respiration_rate(m::Q10Respiration, T::Real)
    return m.R_ref * m.Q10^((T - m.T_ref) / 10)
end

# --- Life Stage (composite type per biological stage) ---

"""
    LifeStage{T<:Real}

A single biological life stage (e.g., egg, larva, pupa, adult, leaf, fruit)
modeled with a distributed delay and optional attrition.

# Fields
- `name::Symbol`: Stage identifier
- `delay::DistributedDelay{T}`: Erlang k-substage delay for maturation
- `dev_rate::AbstractDevelopmentRate`: Temperature-dependent development rate
- `μ::T`: Background mortality rate (per degree-day)
"""
struct LifeStage{T<:Real, D<:AbstractDevelopmentRate}
    name::Symbol
    delay::DistributedDelay{T}
    dev_rate::D
    μ::T

    function LifeStage(name::Symbol, delay::DistributedDelay{T},
                       dev_rate::D, μ::T) where {T<:Real, D<:AbstractDevelopmentRate}
        μ >= 0 || throw(ArgumentError("μ must be non-negative"))
        new{T,D}(name, delay, dev_rate, μ)
    end
end

function LifeStage(name::Symbol, delay::DistributedDelay{T},
                   dev_rate::AbstractDevelopmentRate, μ::Real) where {T}
    LifeStage(name, delay, dev_rate, T(μ))
end

# --- Population (collection of life stages) ---

"""
    Population{T<:Real}

A population of a single species, consisting of ordered life stages
connected by maturation flows.

# Fields
- `name::Symbol`: Species/population name
- `stages::Vector{LifeStage{T}}`: Life stages in developmental order
"""
struct Population{T<:Real}
    name::Symbol
    stages::Vector{<:LifeStage{T}}

    function Population(name::Symbol, stages::Vector{<:LifeStage{T}}) where {T<:Real}
        length(stages) > 0 || throw(ArgumentError("must have at least one life stage"))
        new{T}(name, stages)
    end
end

n_stages(pop::Population) = length(pop.stages)

"""Total number of substages across all life stages."""
function n_substages(pop::Population)
    return sum(s.delay.k for s in pop.stages)
end

"""Total population across all substages and stages."""
function total_population(pop::Population)
    return sum(delay_total(s.delay) for s in pop.stages)
end
