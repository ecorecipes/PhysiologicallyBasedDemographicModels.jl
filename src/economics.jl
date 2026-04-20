"""
Economic analysis types and functions for PBDM bioeconomic models.

Provides cost, revenue, and profit computation for agricultural systems,
including crop damage functions, input costs, and net present value.
"""

# --- Abstract types ---

"""
    AbstractCostFunction

Supertype for per-hectare cost components in a PBDM bioeconomic model
(e.g. [`FixedCost`](@ref), [`VariableCost`](@ref), [`InputCostBundle`](@ref)).
Concrete subtypes must support `total_cost(::AbstractCostFunction, ctx)`.
"""
abstract type AbstractCostFunction end

"""
    AbstractRevenueFunction

Supertype for revenue components in a bioeconomic model (e.g.
[`CropRevenue`](@ref)). Concrete subtypes implement
`revenue(::AbstractRevenueFunction, yield)`.
"""
abstract type AbstractRevenueFunction end

"""
    AbstractDamageFunction

Supertype for crop-damage models that map a pest pressure to a fractional
yield loss (e.g. [`LinearDamageFunction`](@ref),
[`ExponentialDamageFunction`](@ref)). Concrete subtypes implement
`yield_loss(::AbstractDamageFunction, pressure)`.
"""
abstract type AbstractDamageFunction end

# --- Cost functions ---

"""
    FixedCost(amount, label)

A fixed per-hectare cost (e.g., land rent, equipment).
"""
struct FixedCost{T<:Real} <: AbstractCostFunction
    amount::T
    label::Symbol
end

"""
    VariableCost(rate, label)

A cost that scales with an input quantity (e.g., USD/kg pesticide applied).
"""
struct VariableCost{T<:Real} <: AbstractCostFunction
    rate::T
    label::Symbol
end

"""
    InputCostBundle(costs)

A bundle of input costs for a production system (seed, pesticide, fertilizer, labor).
"""
struct InputCostBundle{T<:Real}
    costs::Vector{Pair{Symbol, T}}
end

function InputCostBundle(; kwargs...)
    T = Float64
    pairs = Pair{Symbol, T}[]
    for (k, v) in kwargs
        push!(pairs, k => T(v))
    end
    InputCostBundle{T}(pairs)
end

"""
    total_cost(bundle::InputCostBundle)

Sum all costs in the bundle.
"""
function total_cost(bundle::InputCostBundle)
    return sum(last(p) for p in bundle.costs)
end

"""
    total_cost(costs::Vector{<:AbstractCostFunction}, quantities::Dict{Symbol, <:Real})

Compute total cost given cost functions and input quantities.
"""
function total_cost(costs::Vector{<:AbstractCostFunction},
                    quantities::Dict{Symbol, <:Real}=Dict{Symbol,Float64}())
    total = 0.0
    for c in costs
        if c isa FixedCost
            total += c.amount
        elseif c isa VariableCost
            q = get(quantities, c.label, 0.0)
            total += c.rate * q
        end
    end
    return total
end

# --- Revenue / Yield functions ---

"""
    CropRevenue(price_per_unit, unit_label)

Revenue from crop sales at a given market price.
"""
struct CropRevenue{T<:Real} <: AbstractRevenueFunction
    price_per_unit::T
    unit_label::Symbol
end

"""
    revenue(cr::CropRevenue, yield_quantity)

Compute gross revenue from a crop yield.
"""
function revenue(cr::CropRevenue, yield_quantity::Real)
    return cr.price_per_unit * yield_quantity
end

# --- Damage functions ---

"""
    LinearDamageFunction(loss_per_pest)

Yield loss proportional to pest density: `loss = loss_per_pest × pest_density`.
"""
struct LinearDamageFunction{T<:Real} <: AbstractDamageFunction
    loss_per_pest::T
end

"""
    ExponentialDamageFunction(a)

Yield loss as fraction: `damage_fraction = 1 - exp(-a × pest_density)`.
"""
struct ExponentialDamageFunction{T<:Real} <: AbstractDamageFunction
    a::T
end

"""
    yield_loss(df::AbstractDamageFunction, pest_density, potential_yield)

Compute yield loss from pest damage.
"""
function yield_loss(df::LinearDamageFunction, pest_density::Real, potential_yield::Real)
    loss = df.loss_per_pest * pest_density
    return min(loss, potential_yield)
end

function yield_loss(df::ExponentialDamageFunction, pest_density::Real, potential_yield::Real)
    damage_frac = 1.0 - exp(-df.a * pest_density)
    return damage_frac * potential_yield
end

"""
    actual_yield(df::AbstractDamageFunction, pest_density, potential_yield)

Compute actual yield after pest damage.
"""
function actual_yield(df::AbstractDamageFunction, pest_density::Real, potential_yield::Real)
    return potential_yield - yield_loss(df, pest_density, potential_yield)
end

# --- Profit calculation ---

"""
    net_profit(gross_revenue, total_costs)

Simple net profit.
"""
function net_profit(gross_revenue::Real, total_costs::Real)
    return gross_revenue - total_costs
end

"""
    net_profit(revenue_fn, yield_qty, cost_bundle)

Compute net profit from revenue function, yield, and costs.
"""
function net_profit(rev::CropRevenue, yield_qty::Real, costs::InputCostBundle)
    return revenue(rev, yield_qty) - total_cost(costs)
end

"""
    daily_income(annual_profit; days=365)

Convert annual profit to daily income.
"""
function daily_income(annual_profit::Real; days::Int=365)
    return annual_profit / days
end

# --- Net Present Value ---

"""
    npv(cash_flows, discount_rate)

Compute net present value of a vector of periodic cash flows.
`discount_rate` is per-period (e.g., 0.05 for 5% per period).
"""
function npv(cash_flows::AbstractVector{<:Real}, discount_rate::Real)
    total = 0.0
    for (t, cf) in enumerate(cash_flows)
        total += cf / (1 + discount_rate)^t
    end
    return total
end

"""
    benefit_cost_ratio(benefits, costs)

Compute benefit-cost ratio.
"""
function benefit_cost_ratio(benefits::Real, costs::Real)
    costs <= 0 && return Inf
    return benefits / costs
end

# --- Climate-yield regression helpers ---

"""
    RainfallYieldModel(a, b, c)

Quadratic rainfall-yield model: `yield = a × mm² + b × mm + c`.
"""
struct RainfallYieldModel{T<:Real}
    a::T
    b::T
    c::T
end

"""
    predict_yield(m, args...)

Apply a yield-prediction model to weather/management drivers.
Methods exist for [`RainfallYieldModel`](@ref) (rainfall only) and
[`WeatherYieldModel`](@ref) (degree-days × rainfall).
"""
function predict_yield(m::RainfallYieldModel, rainfall_mm::Real)
    return m.a * rainfall_mm^2 + m.b * rainfall_mm + m.c
end

"""
    WeatherYieldModel(β_dd, β_rain, β_interaction, intercept)

Weather-driven yield model: `yield = intercept + β_dd × dd + β_rain × rain + β_int × dd × rain`.
"""
struct WeatherYieldModel{T<:Real}
    β_dd::T
    β_rain::T
    β_interaction::T
    intercept::T
end

function predict_yield(m::WeatherYieldModel, degree_days::Real, rainfall_mm::Real)
    return m.intercept + m.β_dd * degree_days + m.β_rain * rainfall_mm +
           m.β_interaction * degree_days * rainfall_mm
end
