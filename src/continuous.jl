"""
Continuous-time PBDM formulations.

Provides problem types for expressing PBDM dynamics as continuous differential
equations rather than discrete daily steps. Three formulations are supported:

1. **Aggregate ODE** (`ContinuousPBDMProblem`): Aggregate-biomass ODEs using
   the linear chain trick (Erlang substages → ODE state variables).
   Follows Gutierrez (1994) and Schreiber & Gutierrez (1998).
   Solved via OrdinaryDiffEq.jl (package extension).

2. **Delay DE** (`DelayPBDMProblem`): Delay differential equations where
   maturation is modeled as an explicit distributed delay with history.
   Solved via DelayDiffEq.jl (package extension).

3. **Size-structured PDE** (`PSPMProblem`): McKendrick–von Foerster PDE
   discretized via method of lines (FMU, EBT, or CM methods; see
   Joshi et al. 2023). Solved as a large ODE system via OrdinaryDiffEq.jl.

All problem types reuse the existing package types (FraserGilbertResponse,
Q10Respiration, development rates, etc.) so parameterization is consistent
with the discrete-time models.
"""

# ============================================================================
# Abstract supertypes
# ============================================================================

"""
    AbstractContinuousPBDM

Supertype for continuous-time PBDM formulations.
"""
abstract type AbstractContinuousPBDM end

"""
    AbstractPSPMMethod

Supertype for numerical methods that discretize the McKendrick–von Foerster
PDE into an ODE system (method of lines).
"""
abstract type AbstractPSPMMethod end

"""
    AbstractPSPMSpecies{T<:Real}

Supertype for species in PSPM formulations.

Subtypes:
- `PSPMSpecies`: single PDE for one continuous state variable.
- `StagedPSPMSpecies`: multiple sequentially-coupled PDEs representing
  distinct life stages (e.g., egg → larva → pupa → adult).
"""
abstract type AbstractPSPMSpecies{T<:Real} end

# ============================================================================
# Trophic species specification (shared across all continuous formulations)
# ============================================================================

"""
    ContinuousSpecies{T<:Real}

Parameterization of one biological species for continuous-time PBDM dynamics.

Bundles the physiological parameters needed to write the continuous ODE/DDE/PDE:
development rate, functional response, respiration, and life-stage structure.

# Fields
- `name::Symbol` — species identifier
- `n_stages::Int` — number of life stages
- `k::Vector{Int}` — substages per life stage (Erlang shape parameters)
- `τ::Vector{T}` — mean developmental time per stage (degree-days)
- `μ::Vector{T}` — background mortality per stage (per degree-day)
- `dev_rate::Vector{<:AbstractDevelopmentRate}` — development rate model per stage
- `fr::AbstractFunctionalResponse` — functional response for resource acquisition
- `resp::Q10Respiration{T}` — respiration model
- `demand_rate::T` — per-capita demand rate (d)
- `conversion_efficiency::T` — assimilation efficiency (ε)
- `intrinsic_rate::T` — intrinsic growth rate (basal resources only; 0 otherwise)
- `carrying_capacity::T` — carrying capacity (basal resources only; Inf otherwise)
"""
struct ContinuousSpecies{T<:Real}
    name::Symbol
    n_stages::Int
    k::Vector{Int}
    τ::Vector{T}
    μ::Vector{T}
    dev_rate::Vector{<:AbstractDevelopmentRate}
    fr::AbstractFunctionalResponse
    resp::Q10Respiration{T}
    demand_rate::T
    conversion_efficiency::T
    intrinsic_rate::T
    carrying_capacity::T
end

"""
    ContinuousSpecies(pop::Population; fr, resp, demand_rate,
                      conversion_efficiency=1.0, intrinsic_rate=0.0,
                      carrying_capacity=Inf)

Construct a `ContinuousSpecies` from an existing discrete `Population`,
extracting stage structure automatically.
"""
function ContinuousSpecies(pop::Population{T};
                           fr::AbstractFunctionalResponse,
                           resp::Q10Respiration,
                           demand_rate::Real,
                           conversion_efficiency::Real=1.0,
                           intrinsic_rate::Real=0.0,
                           carrying_capacity::Real=Inf) where {T}
    ns = n_stages(pop)
    ks = [s.delay.k for s in pop.stages]
    τs = T[s.delay.τ for s in pop.stages]
    μs = T[s.μ for s in pop.stages]
    devs = [s.dev_rate for s in pop.stages]
    ContinuousSpecies{T}(pop.name, ns, ks, τs, μs, devs,
                         fr, resp, T(demand_rate), T(conversion_efficiency),
                         T(intrinsic_rate), T(carrying_capacity))
end

"""
    ContinuousSpecies(name::Symbol; k, τ, μ, dev_rate, fr, resp, demand_rate,
                      conversion_efficiency=1.0, intrinsic_rate=0.0,
                      carrying_capacity=Inf)

Construct a `ContinuousSpecies` directly from parameter vectors.
"""
function ContinuousSpecies(name::Symbol;
                           k::Vector{Int},
                           τ::Vector{<:Real},
                           μ::Vector{<:Real},
                           dev_rate::Vector{<:AbstractDevelopmentRate},
                           fr::AbstractFunctionalResponse,
                           resp::Q10Respiration,
                           demand_rate::Real,
                           conversion_efficiency::Real=1.0,
                           intrinsic_rate::Real=0.0,
                           carrying_capacity::Real=Inf)
    T = Float64
    ns = length(k)
    length(τ) == ns || throw(DimensionMismatch("τ must have length $(ns)"))
    length(μ) == ns || throw(DimensionMismatch("μ must have length $(ns)"))
    length(dev_rate) == ns || throw(DimensionMismatch("dev_rate must have length $(ns)"))
    ContinuousSpecies{T}(name, ns, k, T.(τ), T.(μ), dev_rate,
                         fr, resp, T(demand_rate), T(conversion_efficiency),
                         T(intrinsic_rate), T(carrying_capacity))
end

"""
Total number of ODE state variables for one species (sum of all substages).
"""
_total_substages(sp::ContinuousSpecies) = sum(sp.k)

# ============================================================================
# Trophic interaction specification
# ============================================================================

"""
    ContinuousTrophicLink

Directed trophic interaction between two continuous species.

# Fields
- `consumer::Symbol` — name of the consuming species
- `resource::Symbol` — name of the resource species
"""
struct ContinuousTrophicLink
    consumer::Symbol
    resource::Symbol
end

# ============================================================================
# Aggregate ODE problem (linear chain trick)
# ============================================================================

"""
    ContinuousPBDMProblem{T<:Real}

Aggregate-biomass ODE formulation of a PBDM trophic system.

Each life stage's Erlang k-substage distributed delay is represented as k
ODE state variables (linear chain trick). The resulting system is:

For substage ℓ of stage j of species i:
    dW_{i,j,ℓ}/dt = r_{i,j} · g_j(T) · W_{i,j,ℓ-1}
                   - r_{i,j} · g_j(T) · W_{i,j,ℓ}
                   - μ_{i,j} · g_j(T) · W_{i,j,ℓ}
                   - σ_i(t) · W_{i,j,ℓ}  (stress-driven attrition)

where r_{i,j} = k_j / τ_j, g_j(T) = development_rate(dev_j, T(t)),
and σ_i depends on the supply/demand ratio via trophic interactions.

For basal resources, logistic growth replaces trophic acquisition.

# Fields
- `species::Vector{ContinuousSpecies{T}}` — species in trophic order
- `links::Vector{ContinuousTrophicLink}` — trophic interactions
- `u0::Vector{T}` — initial state vector (flattened substages)
- `tspan::Tuple{T,T}` — time span
- `T_forcing` — temperature forcing: constant `T` or callable `T(t)`
- `p` — additional user parameters

# Requires
OrdinaryDiffEq.jl package extension for `solve()`.
"""
struct ContinuousPBDMProblem{T<:Real, TF} <: AbstractContinuousPBDM
    species::Vector{ContinuousSpecies{T}}
    links::Vector{ContinuousTrophicLink}
    u0::Vector{T}
    tspan::Tuple{T,T}
    T_forcing::TF
    p::Any
end

function ContinuousPBDMProblem(;
        species::Vector{<:ContinuousSpecies},
        links::Vector{ContinuousTrophicLink}=ContinuousTrophicLink[],
        u0::Union{Nothing,Vector{<:Real}}=nothing,
        tspan::Tuple{<:Real,<:Real},
        T_forcing=25.0,
        p=nothing)
    T = Float64

    # Build initial state vector from species' substage masses if not provided
    if u0 === nothing
        u0_vec = T[]
        for sp in species
            for j in 1:sp.n_stages
                append!(u0_vec, fill(T(0), sp.k[j]))
            end
        end
        # Default: distribute W0 evenly across substages of stage 1
        # User should provide u0 for non-trivial initialization
    else
        u0_vec = T.(u0)
    end

    n_expected = sum(_total_substages(sp) for sp in species)
    length(u0_vec) == n_expected ||
        throw(DimensionMismatch("u0 has length $(length(u0_vec)), expected $(n_expected) (sum of all substages)"))

    ContinuousPBDMProblem{T, typeof(T_forcing)}(
        species, links, u0_vec, (T(tspan[1]), T(tspan[2])), T_forcing, p)
end

"""
    species_state_ranges(prob::ContinuousPBDMProblem)

Return a vector of `(species_name, stage_index, UnitRange)` tuples mapping
each life stage to its indices in the flattened state vector.
"""
function species_state_ranges(prob::ContinuousPBDMProblem)
    ranges = Tuple{Symbol, Int, UnitRange{Int}}[]
    idx = 1
    for sp in prob.species
        for j in 1:sp.n_stages
            r = idx:(idx + sp.k[j] - 1)
            push!(ranges, (sp.name, j, r))
            idx += sp.k[j]
        end
    end
    return ranges
end

"""
    species_total_ranges(prob::ContinuousPBDMProblem)

Return a Dict mapping each species name to the full range of indices
in the state vector for that species.
"""
function species_total_ranges(prob::ContinuousPBDMProblem)
    result = Dict{Symbol, UnitRange{Int}}()
    idx = 1
    for sp in prob.species
        n = _total_substages(sp)
        result[sp.name] = idx:(idx + n - 1)
        idx += n
    end
    return result
end

"""
    _get_temperature(T_forcing, t)

Get temperature at time `t` from forcing (constant or callable).
"""
_get_temperature(T_f::Real, t) = T_f
_get_temperature(T_f, t) = T_f(t)

"""
    flatten_population(pop::Population)

Extract the current substage masses from a `Population` into a flat vector,
suitable for initializing a `ContinuousPBDMProblem`.
"""
function flatten_population(pop::Population{T}) where {T}
    u = T[]
    for stage in pop.stages
        append!(u, stage.delay.W)
    end
    return u
end

"""
    flatten_populations(pops::Vector{<:Population})

Flatten multiple populations into a single state vector.
"""
function flatten_populations(pops::Vector{<:Population{T}}) where {T}
    u = T[]
    for pop in pops
        append!(u, flatten_population(pop))
    end
    return u
end

"""
Placeholder for solve dispatch — implemented by OrdinaryDiffEq extension.
"""
function solve_continuous end

# ============================================================================
# Delay DE problem
# ============================================================================

"""
    DelayPBDMProblem{T<:Real}

Delay differential equation formulation of a PBDM system.

Instead of the linear chain trick, each life stage is represented as a
single state variable with an explicit distributed delay modeled via the
history function. The Erlang-distributed maturation time τ with shape k
is handled via the gamma-distributed delay kernel.

State variables are aggregate biomass per species per stage (no substages).
Delay terms appear as `h(t - τ_j)` weighted by the gamma kernel.

For practical computation, DelayDiffEq.jl's `MethodOfSteps` solvers are used,
and the distributed delay is either:
- Approximated by k discrete lags (taps from the Erlang distribution), or
- Computed via quadrature over the history.

# Fields
- `species::Vector{ContinuousSpecies{T}}` — species parameterizations
- `links::Vector{ContinuousTrophicLink}` — trophic interactions
- `u0::Vector{T}` — initial state (one value per species per stage)
- `h0` — history function `h(p, t) -> Vector` for t < t0
- `tspan::Tuple{T,T}` — time span
- `T_forcing` — temperature forcing
- `p` — additional parameters

# Requires
DelayDiffEq.jl package extension for `solve()`.
"""
struct DelayPBDMProblem{T<:Real, TF, H} <: AbstractContinuousPBDM
    species::Vector{ContinuousSpecies{T}}
    links::Vector{ContinuousTrophicLink}
    u0::Vector{T}
    h0::H
    tspan::Tuple{T,T}
    T_forcing::TF
    p::Any
end

function DelayPBDMProblem(;
        species::Vector{<:ContinuousSpecies},
        links::Vector{ContinuousTrophicLink}=ContinuousTrophicLink[],
        u0::Vector{<:Real},
        h0=nothing,
        tspan::Tuple{<:Real,<:Real},
        T_forcing=25.0,
        p=nothing)
    T = Float64

    # One state variable per species per stage
    n_expected = sum(sp.n_stages for sp in species)
    length(u0) == n_expected ||
        throw(DimensionMismatch("u0 has length $(length(u0)), expected $(n_expected) (one per species per stage)"))

    # Default history: constant at u0
    if h0 === nothing
        u0_copy = T.(u0)
        h0_fn = (p, t) -> u0_copy
    else
        h0_fn = h0
    end

    DelayPBDMProblem{T, typeof(T_forcing), typeof(h0_fn)}(
        species, links, T.(u0), h0_fn,
        (T(tspan[1]), T(tspan[2])), T_forcing, p)
end

"""
Placeholder — implemented by DelayDiffEq extension.
"""
function solve_delay end

# ============================================================================
# PSPM / McKendrick–von Foerster PDE problem
# ============================================================================

"""
    FixedMeshUpwind <: AbstractPSPMMethod

Fixed mesh upwind (FMU) discretization of the McKendrick–von Foerster PDE.
First-order upwind finite differences on a fixed size grid.
Fast but diffusive. From Zhang et al. (2017), Joshi et al. (2023).
"""
struct FixedMeshUpwind <: AbstractPSPMMethod
    n_mesh::Int
end
FixedMeshUpwind(; n_mesh::Int=100) = FixedMeshUpwind(n_mesh)

"""
    EscalatorBoxcarTrain <: AbstractPSPMMethod

Escalator Boxcar Train (EBT) method: tracks cohorts as moving "boxcars"
with position (mean size) and number. Highly accurate for sharp
distributions. From de Roos (1988), Joshi et al. (2023).
"""
struct EscalatorBoxcarTrain <: AbstractPSPMMethod
    max_cohorts::Int
    cohort_interval::Float64
end
EscalatorBoxcarTrain(; max_cohorts::Int=200, cohort_interval::Float64=1.0) =
    EscalatorBoxcarTrain(max_cohorts, cohort_interval)

"""
    CharacteristicMethod <: AbstractPSPMMethod

Method of characteristics: tracks cohort boundaries along
characteristic curves rather than mean sizes (as in EBT). Each
cohort tracks (x_lower, x_upper, N_count). The boundaries
propagate at the local growth rate, giving exact characteristic
tracking for smooth growth functions. From Joshi et al. (2023).
"""
struct CharacteristicMethod <: AbstractPSPMMethod
    max_cohorts::Int
end
CharacteristicMethod(; max_cohorts::Int=200) = CharacteristicMethod(max_cohorts)

"""
    ImplicitFixedMeshUpwind <: AbstractPSPMMethod

Implicit Fixed Mesh Upwind (IFMU): operator-split semi-implicit
scheme treating mortality implicitly for stiff problems. Advection
uses FMU (explicit upwind), mortality is solved implicitly via a
stiff ODE solver. From Zhang et al. (2017), Joshi et al. (2023).
"""
struct ImplicitFixedMeshUpwind <: AbstractPSPMMethod
    n_mesh::Int
end
ImplicitFixedMeshUpwind(; n_mesh::Int=100) = ImplicitFixedMeshUpwind(n_mesh)

"""
    ImplicitEscalatorBoxcarTrain <: AbstractPSPMMethod

Implicit EBT (IEBT): EBT cohort tracking with implicit mortality
treatment for stiff problems. Uses a stiff ODE solver so that
high-mortality regimes remain stable without tiny time steps.
From Joshi et al. (2023).
"""
struct ImplicitEscalatorBoxcarTrain <: AbstractPSPMMethod
    max_cohorts::Int
    cohort_interval::Float64
end
ImplicitEscalatorBoxcarTrain(; max_cohorts::Int=200, cohort_interval::Float64=1.0) =
    ImplicitEscalatorBoxcarTrain(max_cohorts, cohort_interval)

"""
    ImplicitCharacteristicMethod <: AbstractPSPMMethod

Implicit CM (ICM): characteristic method with implicit mortality
for stiff problems. Uses a stiff ODE solver.
From Joshi et al. (2023).
"""
struct ImplicitCharacteristicMethod <: AbstractPSPMMethod
    max_cohorts::Int
end
ImplicitCharacteristicMethod(; max_cohorts::Int=200) = ImplicitCharacteristicMethod(max_cohorts)

"""
    LaxFriedrichsUpwind <: AbstractPSPMMethod

Lax-Friedrichs upwind (FMU-LF): fixed mesh with Lax-Friedrichs
numerical flux. Adds controlled numerical diffusion for stability
on coarser meshes, at the cost of smearing sharp features.
The flux is: F_{j+1/2} = (g·n_j + g·n_{j+1})/2 − (α/2)(n_{j+1} − n_j)
where α = max|g(x)| over the domain.
From Joshi et al. (2023).
"""
struct LaxFriedrichsUpwind <: AbstractPSPMMethod
    n_mesh::Int
end
LaxFriedrichsUpwind(; n_mesh::Int=100) = LaxFriedrichsUpwind(n_mesh)

"""
    ImplicitLaxFriedrichsUpwind <: AbstractPSPMMethod

Implicit Lax-Friedrichs upwind (IFMU-LF): Lax-Friedrichs flux
with implicit mortality for stiff problems. Combines the diffusive
stability of the LF flux with implicit treatment of mortality.
From Joshi et al. (2023).
"""
struct ImplicitLaxFriedrichsUpwind <: AbstractPSPMMethod
    n_mesh::Int
end
ImplicitLaxFriedrichsUpwind(; n_mesh::Int=100) = ImplicitLaxFriedrichsUpwind(n_mesh)

"""
    PSPMSpecies{T<:Real}

Size-structured species for PSPM formulation.

The McKendrick–von Foerster PDE for this species is:
    ∂n/∂t + ∂(g·n)/∂x = -μ·n

where x is the physiological state variable (e.g., body mass),
g(x,E,t) is the growth rate, μ(x,E,t) is the mortality rate,
and n(x,t) is the size density.

Boundary condition at x = x_b:
    g(x_b) · n(x_b, t) = ∫ β(x,E,t) · n(x,t) dx

where β is the fecundity rate.

# Fields
- `name::Symbol` — species identifier
- `x_birth::T` — birth size (lower boundary)
- `x_max::T` — maximum size (upper boundary)
- `growth_rate` — `g(x, E, t)` → growth rate at size x
- `mortality_rate` — `μ(x, E, t)` → mortality rate at size x
- `fecundity_rate` — `β(x, E, t)` → fecundity at size x
- `init_density` — `n₀(x)` → initial size distribution
"""
struct PSPMSpecies{T<:Real, G, M, F, I} <: AbstractPSPMSpecies{T}
    name::Symbol
    x_birth::T
    x_max::T
    growth_rate::G
    mortality_rate::M
    fecundity_rate::F
    init_density::I
end

function PSPMSpecies(name::Symbol;
                     x_birth::Real, x_max::Real,
                     growth_rate, mortality_rate, fecundity_rate,
                     init_density)
    T = Float64
    PSPMSpecies{T, typeof(growth_rate), typeof(mortality_rate),
                typeof(fecundity_rate), typeof(init_density)}(
        name, T(x_birth), T(x_max),
        growth_rate, mortality_rate, fecundity_rate, init_density)
end

# ============================================================================
# Multi-stage PSPM species
# ============================================================================

"""
    PSPMStage{G, M}

One life stage in a staged PSPM species.  Each stage has its own
McKendrick–von Foerster PDE over physiological age ``x \\in [x_{\\min}, x_{\\max}]``.

Individuals enter the stage at `x_birth` and leave at `x_max`.
Inter-stage boundary conditions are handled automatically by
`StagedPSPMSpecies`.

# Fields
- `name::Symbol` — stage identifier (e.g., `:egg`, `:larva`)
- `x_birth::Float64` — lower boundary of physiological age
- `x_max::Float64` — upper boundary of physiological age
- `growth_rate` — `g(x, E, t)` → development (advection) rate
- `mortality_rate` — `μ(x, E, t)` → mortality rate
"""
struct PSPMStage{G, M}
    name::Symbol
    x_birth::Float64
    x_max::Float64
    growth_rate::G
    mortality_rate::M
end

function PSPMStage(name::Symbol; x_birth::Real=0.0, x_max::Real=1.0,
                   growth_rate, mortality_rate)
    PSPMStage{typeof(growth_rate), typeof(mortality_rate)}(
        name, Float64(x_birth), Float64(x_max), growth_rate, mortality_rate)
end

"""
    StagedPSPMSpecies{T<:Real, R, I} <: AbstractPSPMSpecies{T}

Multi-stage physiologically structured species where individuals pass
sequentially through life stages (e.g., egg → larva → pupa → adult).

Each stage ``i`` has its own McKendrick–von Foerster PDE.  Boundary
conditions are coupled automatically:

- **Stage 1**: influx = `reproduction_flux(E, t, stage_totals)`
  where `stage_totals` is a vector of total abundance per stage.
- **Stage ``i > 1``**: influx = outflow from stage ``i-1`` at
  ``x = x_{\\max}``, i.e., ``g_{i-1}(x_{\\max}) \\cdot n_{i-1}(x_{\\max})``.

# Fields
- `name::Symbol` — species identifier
- `stages::Vector{<:PSPMStage}` — ordered life stages
- `reproduction_flux` — `(E, t, stage_totals) → flux` — birth flux into stage 1
- `init_density` — `(stage_idx, x) → density` — initial size distribution per stage

# Example
```julia
stages = [
    PSPMStage(:egg;   growth_rate=(x,E,t)->0.5, mortality_rate=(x,E,t)->0.01),
    PSPMStage(:larva; growth_rate=(x,E,t)->0.3, mortality_rate=(x,E,t)->0.02),
    PSPMStage(:adult; growth_rate=(x,E,t)->0.1, mortality_rate=(x,E,t)->0.05),
]
sp = StagedPSPMSpecies(:insect;
    stages = stages,
    reproduction_flux = (E, t, totals) -> 10.0 * totals[3],
    init_density = (stage, x) -> stage == 3 ? 1.0 : 0.0)
```
"""
struct StagedPSPMSpecies{T<:Real, R, I} <: AbstractPSPMSpecies{T}
    name::Symbol
    stages::Vector{<:PSPMStage}
    reproduction_flux::R
    init_density::I
end

function StagedPSPMSpecies(name::Symbol;
                           stages::Vector{<:PSPMStage},
                           reproduction_flux,
                           init_density)
    T = Float64
    StagedPSPMSpecies{T, typeof(reproduction_flux), typeof(init_density)}(
        name, stages, reproduction_flux, init_density)
end

"""
    n_pspm_stages(sp::StagedPSPMSpecies)

Number of life stages in a staged PSPM species.
"""
n_pspm_stages(sp::StagedPSPMSpecies) = length(sp.stages)

"""
    PSPMProblem{T<:Real}

Physiologically Structured Population Model problem.

Discretizes the McKendrick–von Foerster PDE using the specified method
(FMU, EBT, or CM) and solves the resulting ODE system.

Accepts both single-PDE species (`PSPMSpecies`) and multi-stage species
(`StagedPSPMSpecies`).

# Fields
- `species::Vector{<:AbstractPSPMSpecies{T}}` — structured species
- `environment` — environment computation: `E(u, t)` or constant
- `method::AbstractPSPMMethod` — discretization method
- `tspan::Tuple{T,T}` — time span
- `p` — additional parameters

# Requires
OrdinaryDiffEq.jl package extension for `solve()`.
"""
struct PSPMProblem{T<:Real, E} <: AbstractContinuousPBDM
    species::Vector{<:AbstractPSPMSpecies{T}}
    environment::E
    method::AbstractPSPMMethod
    tspan::Tuple{T,T}
    p::Any
end

function PSPMProblem(;
        species::Vector{<:AbstractPSPMSpecies},
        environment=nothing,
        method::AbstractPSPMMethod=FixedMeshUpwind(),
        tspan::Tuple{<:Real,<:Real},
        p=nothing)
    T = Float64
    PSPMProblem{T, typeof(environment)}(
        species, environment, method, (T(tspan[1]), T(tspan[2])), p)
end

"""
Placeholder — implemented by OrdinaryDiffEq extension.
"""
function solve_pspm end

"""
Placeholder — implemented by OrdinaryDiffEq extension.
"""
function staged_species_stage_totals end

# ============================================================================
# Continuous solution type
# ============================================================================

"""
    ContinuousPBDMSolution{T<:Real}

Result of solving a continuous-time PBDM problem.

Wraps the raw SciML solution and adds species-aware accessors.

# Fields
- `t::Vector{T}` — time points
- `u::Matrix{T}` — state matrix (n_states × n_timepoints)
- `species_names::Vector{Symbol}` — species names
- `species_ranges::Dict{Symbol, UnitRange{Int}}` — index ranges per species
- `retcode::Symbol` — solver return code
- `raw_sol` — underlying SciML solution object (for interpolation etc.)
"""
struct ContinuousPBDMSolution{T<:Real, S}
    t::Vector{T}
    u::Matrix{T}
    species_names::Vector{Symbol}
    species_ranges::Dict{Symbol, UnitRange{Int}}
    retcode::Symbol
    raw_sol::S
end

"""
    species_trajectory(sol::ContinuousPBDMSolution, name::Symbol)

Extract the total biomass trajectory for a named species.
"""
function species_trajectory(sol::ContinuousPBDMSolution, name::Symbol)
    haskey(sol.species_ranges, name) ||
        throw(KeyError("species $name not found"))
    r = sol.species_ranges[name]
    return vec(sum(sol.u[r, :]; dims=1))
end

"""
    stage_trajectory(sol::ContinuousPBDMSolution, name::Symbol, stage_idx::Int,
                     prob::ContinuousPBDMProblem)

Extract the biomass trajectory for a specific stage of a named species.
"""
function stage_trajectory(sol::ContinuousPBDMSolution, name::Symbol,
                          stage_idx::Int, prob::ContinuousPBDMProblem)
    ranges = species_state_ranges(prob)
    for (sn, si, r) in ranges
        if sn == name && si == stage_idx
            return vec(sum(sol.u[r, :]; dims=1))
        end
    end
    throw(KeyError("stage $stage_idx of species $name not found"))
end

function Base.show(io::IO, sol::ContinuousPBDMSolution)
    nt = length(sol.t)
    ns = length(sol.species_names)
    nv = size(sol.u, 1)
    print(io, "ContinuousPBDMSolution($nt timepoints, $ns species, $nv state vars, retcode=$(sol.retcode))")
end
