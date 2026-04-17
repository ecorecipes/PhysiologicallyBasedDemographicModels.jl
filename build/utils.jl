"""
Utility functions for PhysiologicallyBasedDemographicModels.jl
"""

"""
    photoperiod(latitude, day_of_year)

Estimate day length (hours) using the CBM model.
`latitude` in degrees, `day_of_year` ∈ [1, 365].
"""
function photoperiod(latitude::Real, day_of_year::Int)
    θ = 0.2163108 + 2 * atan(0.9671396 * tan(0.00860 * (day_of_year - 186)))
    ϕ = asin(0.39795 * cos(θ))
    lat_rad = deg2rad(latitude)
    arg = (sin(deg2rad(-0.8333)) + sin(lat_rad) * sin(ϕ)) / (cos(lat_rad) * cos(ϕ))
    arg = clamp(arg, -1.0, 1.0)
    return 24 - (24 / π) * acos(arg)
end

"""
    degree_days_sine(T_min, T_max, T_lower; T_upper=Inf)

Compute daily degree-day accumulation using the single sine method.
More accurate than simple mean-based degree-days.
"""
function degree_days_sine(T_min::Real, T_max::Real, T_lower::Real;
                          T_upper::Real=Inf)
    T_mean = (T_min + T_max) / 2
    amplitude = (T_max - T_min) / 2

    if T_max <= T_lower
        return 0.0
    elseif T_min >= T_lower
        if T_upper == Inf || T_max <= T_upper
            return T_mean - T_lower
        else
            # Both thresholds matter
            return _dd_sine_both(T_mean, amplitude, T_lower, T_upper)
        end
    else
        # Partial day above threshold
        return _dd_sine_partial(T_mean, amplitude, T_lower, T_upper)
    end
end

function _dd_sine_partial(T_mean, amplitude, T_lower, T_upper)
    amplitude == 0 && return max(0.0, T_mean - T_lower)
    θ1 = asin(clamp((T_lower - T_mean) / amplitude, -1.0, 1.0))
    dd = ((T_mean - T_lower) * (π / 2 - θ1) + amplitude * cos(θ1)) / π
    if T_upper < Inf && T_mean + amplitude > T_upper
        θ2 = asin(clamp((T_upper - T_mean) / amplitude, -1.0, 1.0))
        dd -= ((T_mean - T_upper) * (π / 2 - θ2) + amplitude * cos(θ2)) / π
    end
    return max(0.0, dd)
end

function _dd_sine_both(T_mean, amplitude, T_lower, T_upper)
    θ1 = asin(clamp((T_lower - T_mean) / amplitude, -1.0, 1.0))
    θ2 = asin(clamp((T_upper - T_mean) / amplitude, -1.0, 1.0))
    dd = ((T_mean - T_lower) * (π / 2 - θ1) + amplitude * cos(θ1)) / π
    dd -= ((T_mean - T_upper) * (π / 2 - θ2) + amplitude * cos(θ2)) / π
    return max(0.0, dd)
end

"""
    make_population(name, stage_specs, dev_rate; mortality=0.0)

Convenience constructor for a Population from a vector of (name, k, τ) tuples.
All stages share the same development rate model.

# Arguments
- `name::Symbol`: Population name
- `stage_specs`: Vector of `(name::Symbol, k::Int, τ::Real)` tuples
- `dev_rate`: Development rate model applied to all stages
- `mortality`: Constant mortality rate for all stages (default 0)
"""
function make_population(name::Symbol,
                         stage_specs::Vector{<:Tuple{Symbol, Int, <:Real}},
                         dev_rate::AbstractDevelopmentRate;
                         mortality::Real=0.0)
    T = typeof(float(stage_specs[1][3]))
    stages = LifeStage{T, typeof(dev_rate)}[]
    for (sname, k, τ) in stage_specs
        delay = DistributedDelay(k, T(τ))
        push!(stages, LifeStage(sname, delay, dev_rate, T(mortality)))
    end
    Population(name, stages)
end
