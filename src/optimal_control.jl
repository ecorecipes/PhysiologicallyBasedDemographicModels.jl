"""
Optimal management control types for PBDM systems.

Provides problem types for bioeconomic optimization following the approach
of Regev, Gutierrez & Schreiber (1998) and Gutierrez & Regev (2005).

The optimization uses Euler-discretized PBDM dynamics with JuMP/Ipopt,
following the pattern established in EpiPolicies (Frost & Montes 2025).

Core types are defined here (no JuMP dependency). The `optimize_management`
solver is provided by the JuMP package extension (`ext/JuMPExt.jl`).
"""

# ================================================================
# Trophic level parameterization (simplified, for optimization)
# ================================================================

"""
    TrophicLevel{T<:Real}

Simplified aggregate-biomass parameterization of one trophic level for
management optimization. This is NOT a full distributed-delay population;
it captures only the key mass-balance parameters needed for the
discrete-time optimal control formulation.

# Fields
- `name::Symbol` — level identifier (e.g., `:plant`, `:herbivore`, `:predator`)
- `demand_rate::T` — per-capita demand rate (d)
- `fr::AbstractFunctionalResponse` — functional response for resource acquisition
- `resp::Q10Respiration{T}` — respiration model
- `conversion_efficiency::T` — assimilation efficiency (ε)
- `intrinsic_rate::T` — intrinsic growth rate (only for basal resource; 0 otherwise)
- `carrying_capacity::T` — carrying capacity (only for basal resource; Inf otherwise)
"""
struct TrophicLevel{T<:Real}
    name::Symbol
    demand_rate::T
    fr::AbstractFunctionalResponse
    resp::Q10Respiration{T}
    conversion_efficiency::T
    intrinsic_rate::T
    carrying_capacity::T
end

function TrophicLevel(name::Symbol; demand_rate::Real,
                      fr::AbstractFunctionalResponse, resp::Q10Respiration,
                      conversion_efficiency::Real=1.0,
                      intrinsic_rate::Real=0.0,
                      carrying_capacity::Real=Inf)
    T = Float64
    TrophicLevel{T}(name, T(demand_rate), fr, resp,
                    T(conversion_efficiency), T(intrinsic_rate),
                    T(carrying_capacity))
end

# ================================================================
# Management control actions
# ================================================================

"""
    AbstractManagementAction

Supertype for management interventions in a PBDM system.
Each action modifies the dynamics of one or more trophic levels.
"""
abstract type AbstractManagementAction end

"""
    PesticideControl(; name, target, max_rate, efficacy, cost_weight)

Chemical control that adds mortality to a target trophic level.

During optimization, the control variable `u(t) ∈ [0, max_rate]`
applies additional mortality `efficacy × u(t)` to the target population.

# Fields
- `name::Symbol` — control identifier
- `target::Int` — index of the target trophic level
- `max_rate::Float64` — maximum application rate per time step
- `efficacy::Float64` — mortality per unit application rate
- `cost_weight::Float64` — quadratic cost penalty weight in objective
"""
struct PesticideControl <: AbstractManagementAction
    name::Symbol
    target::Int
    max_rate::Float64
    efficacy::Float64
    cost_weight::Float64
end

function PesticideControl(; name::Symbol, target::Int, max_rate::Real,
                          efficacy::Real=1.0, cost_weight::Real=1.0)
    PesticideControl(name, target, Float64(max_rate), Float64(efficacy),
                     Float64(cost_weight))
end

"""
    BiologicalReleaseControl(; name, target, max_rate, cost_weight)

Augmentative biological control that adds biomass to a target
(natural enemy) trophic level.

During optimization, the control variable `u(t) ∈ [0, max_rate]`
injects `u(t)` units of biomass into the target predator population.

# Fields
- `name::Symbol` — control identifier
- `target::Int` — index of the augmented trophic level
- `max_rate::Float64` — maximum release rate per time step
- `cost_weight::Float64` — quadratic cost penalty weight
"""
struct BiologicalReleaseControl <: AbstractManagementAction
    name::Symbol
    target::Int
    max_rate::Float64
    cost_weight::Float64
end

function BiologicalReleaseControl(; name::Symbol, target::Int, max_rate::Real,
                                  cost_weight::Real=1.0)
    BiologicalReleaseControl(name, target, Float64(max_rate), Float64(cost_weight))
end

"""
    HarvestControl(; name, target, max_rate, revenue_per_unit, cost_weight)

Harvesting action that removes biomass from a target trophic level,
generating revenue.

# Fields
- `name::Symbol` — control identifier
- `target::Int` — index of the harvested trophic level
- `max_rate::Float64` — maximum harvest rate per time step
- `revenue_per_unit::Float64` — revenue generated per unit harvested
- `cost_weight::Float64` — quadratic cost penalty weight
"""
struct HarvestControl <: AbstractManagementAction
    name::Symbol
    target::Int
    max_rate::Float64
    revenue_per_unit::Float64
    cost_weight::Float64
end

function HarvestControl(; name::Symbol, target::Int, max_rate::Real,
                        revenue_per_unit::Real=0.0, cost_weight::Real=1.0)
    HarvestControl(name, target, Float64(max_rate), Float64(revenue_per_unit),
                   Float64(cost_weight))
end

# ================================================================
# Management objectives
# ================================================================

"""
    AbstractManagementObjective

Supertype for optimization objectives.
"""
abstract type AbstractManagementObjective end

"""
    MinimizeDamage(; damage_weights, control_cost_weight)

Minimize weighted pest damage plus quadratic control costs.

Objective: `min Σₜ [Σᵢ wᵢ × Mᵢ(t) + Σⱼ cⱼ × uⱼ(t)²]`

where `wᵢ` are damage weights per trophic level and `cⱼ` are
control cost weights per action.

# Fields
- `damage_weights::Vector{Float64}` — weight for each trophic level
  (positive = penalize high biomass, e.g., for pests; zero = ignore)
- `control_cost_weight::Float64` — global scaling of control costs
"""
struct MinimizeDamage <: AbstractManagementObjective
    damage_weights::Vector{Float64}
    control_cost_weight::Float64
end

function MinimizeDamage(; damage_weights::Vector{<:Real},
                        control_cost_weight::Real=1.0)
    MinimizeDamage(Float64.(damage_weights), Float64(control_cost_weight))
end

"""
    MaximizeProfit(; resource_level, price_per_unit, control_cost_weight, discount_rate)

Maximize discounted profit from a resource trophic level.

Objective: `max Σₜ δᵗ × [p × M_resource(t) − Σⱼ cⱼ × uⱼ(t)²]`

# Fields
- `resource_level::Int` — index of the revenue-generating trophic level
- `price_per_unit::Float64` — revenue per unit biomass
- `control_cost_weight::Float64` — global scaling of control costs
- `discount_rate::Float64` — per-period discount factor δ = 1/(1+r)
"""
struct MaximizeProfit <: AbstractManagementObjective
    resource_level::Int
    price_per_unit::Float64
    control_cost_weight::Float64
    discount_rate::Float64
end

function MaximizeProfit(; resource_level::Int, price_per_unit::Real,
                        control_cost_weight::Real=1.0, discount_rate::Real=0.0)
    MaximizeProfit(resource_level, Float64(price_per_unit),
                   Float64(control_cost_weight), Float64(discount_rate))
end

# ================================================================
# Control problem
# ================================================================

"""
    PBDMControlProblem{T<:Real}

Optimal management problem for a PBDM trophic system.

Bundles trophic dynamics, management controls, objective function, and
simulation parameters. Solve with `optimize_management(prob)` (requires
the JuMP package extension).

# Fields
- `levels::Vector{TrophicLevel{T}}` — trophic levels (resource → predator)
- `controls::Vector{<:AbstractManagementAction}` — available management actions
- `objective::AbstractManagementObjective` — optimization objective
- `u0::Vector{T}` — initial biomass for each level
- `tspan::Tuple{T,T}` — time span (start, end)
- `dt::T` — Euler discretization time step
- `T_celsius::T` — temperature (°C, constant or reference)

# Example
```julia
using PhysiologicallyBasedDemographicModels

plant = TrophicLevel(:plant;
    demand_rate=0.5, intrinsic_rate=0.1, carrying_capacity=1000.0,
    fr=FraserGilbertResponse(1.0),
    resp=Q10Respiration(0.01, 2.0, 25.0))

pest = TrophicLevel(:pest;
    demand_rate=0.3,
    fr=FraserGilbertResponse(0.8),
    resp=Q10Respiration(0.02, 2.0, 25.0),
    conversion_efficiency=0.4)

spray = PesticideControl(name=:spray, target=2, max_rate=0.5,
                         efficacy=0.8, cost_weight=10.0)

obj = MinimizeDamage(damage_weights=[0.0, 1.0], control_cost_weight=1.0)

prob = PBDMControlProblem(
    levels=[plant, pest],
    controls=[spray],
    objective=obj,
    u0=[500.0, 50.0],
    tspan=(0.0, 365.0),
    dt=1.0,
    T_celsius=25.0)

# Solve (requires `using JuMP, Ipopt`):
# sol = optimize_management(prob)
```
"""
struct PBDMControlProblem{T<:Real}
    levels::Vector{TrophicLevel{T}}
    controls::Vector{<:AbstractManagementAction}
    objective::AbstractManagementObjective
    u0::Vector{T}
    tspan::Tuple{T,T}
    dt::T
    T_celsius::T
end

function PBDMControlProblem(; levels::Vector{<:TrophicLevel},
                            controls::Vector{<:AbstractManagementAction},
                            objective::AbstractManagementObjective,
                            u0::Vector{<:Real},
                            tspan::Tuple{<:Real,<:Real},
                            dt::Real, T_celsius::Real)
    T = Float64
    PBDMControlProblem{T}(
        levels, controls, objective, T.(u0),
        (T(tspan[1]), T(tspan[2])), T(dt), T(T_celsius))
end

# ================================================================
# Management solution
# ================================================================

"""
    ManagementSolution{T<:Real}

Result of an optimal management computation.

# Fields
- `t::Vector{T}` — time points
- `state::Matrix{T}` — optimal state trajectories (n_levels × n_timesteps)
- `controls::Matrix{T}` — optimal control trajectories (n_controls × n_timesteps)
- `objective_value::T` — optimal objective function value
- `state_labels::Vector{Symbol}` — labels for state variables
- `control_labels::Vector{Symbol}` — labels for control variables
- `termination_status::Symbol` — solver status
"""
struct ManagementSolution{T<:Real}
    t::Vector{T}
    state::Matrix{T}
    controls::Matrix{T}
    objective_value::T
    state_labels::Vector{Symbol}
    control_labels::Vector{Symbol}
    termination_status::Symbol
end

"""
    optimize_management(prob::PBDMControlProblem; solver=nothing, silent=true)

Solve an optimal management problem using nonlinear programming.

Requires the JuMP package extension: `using JuMP, Ipopt` must be called
before this function is available.

See [`PBDMControlProblem`](@ref) for problem setup.
"""
function optimize_management end
