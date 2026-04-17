"""
Distributed delay dynamics — the core PBDM simulation engine.

Implements the Manetsch/Vansickle k-substage distributed delay
difference equations for stage-structured population dynamics
driven by physiological time.
"""

"""
    step_delay!(delay, dd, inflow; μ=0.0, stress=0.0)

Advance a distributed delay by one physiological time step.

# Arguments
- `delay::DistributedDelay`: The delay to update (modified in-place)
- `dd::Real`: Degree-days accumulated this time step
- `inflow::Real`: Input to the first substage
- `μ::Real`: Background mortality rate (per degree-day)
- `stress::Real`: Stress-induced attrition rate (from supply/demand deficit)

# Returns
Named tuple `(outflow, attrition)`:
- `outflow`: quantity flowing out of the last substage (maturation)
- `attrition`: total quantity lost to mortality and stress
"""
function step_delay!(delay::DistributedDelay{T}, dd::Real, inflow::Real;
                     μ::Real=zero(T), stress::Real=zero(T)) where {T}
    dd >= 0 || throw(ArgumentError("dd must be non-negative"))
    inflow >= 0 || throw(ArgumentError("inflow must be non-negative"))
    μ >= 0 || throw(ArgumentError("μ must be non-negative"))
    (zero(T) <= stress <= one(T)) ||
        throw(ArgumentError("stress must be in [0, 1], got $(stress)"))

    k = delay.k
    τ = delay.τ
    W = delay.W

    total_attrition = zero(T)
    total_outflow = zero(T)

    # Split large daily steps into stable substeps so no substage can send more
    # mass onward than it contains.
    total_fraction = ((k / τ) + μ + stress) * dd
    n_substeps = max(1, ceil(Int, total_fraction))
    dd_sub = T(dd / n_substeps)
    flow_fraction = (k / τ) * dd_sub
    attrition_fraction = (μ + stress) * dd_sub

    substep_inflow = T(inflow)
    for _ in 1:n_substeps
        # Process substages from first to last so inflow cascades forward.
        prev_flow = substep_inflow
        substep_inflow = zero(T)

        for i in 1:k
            w_i = W[i]
            loss = attrition_fraction * w_i
            flow_out = flow_fraction * w_i
            remaining = w_i - flow_out - loss

            if remaining < zero(T)
                # Preserve mass balance when floating-point error or an
                # unexpectedly aggressive parameter set would drive a substage
                # slightly negative.
                total_attrition += loss - remaining
                remaining = zero(T)
            else
                total_attrition += loss
            end

            W[i] = remaining + prev_flow
            prev_flow = flow_out
        end

        total_outflow += prev_flow
    end

    return (outflow=total_outflow, attrition=total_attrition)
end

"""
    step_population!(pop, weather_day)

Advance a single-species Population by one calendar day.
Computes degree-days from the weather, then steps each life stage's
distributed delay. Maturation outflow from stage j becomes inflow to stage j+1.
The last stage's outflow is returned as the population's reproductive output.

# Returns
Named tuple with fields:
- `degree_days`: DD accumulated this day for the first stage (legacy summary)
- `stage_degree_days`: DD accumulated this day for each stage
- `maturation`: outflow from the final life stage
- `total_attrition`: total mortality across all stages
- `stage_totals`: vector of total population per stage
"""
function step_population!(pop::Population{T}, weather_day::DailyWeather;
                          inflow::Real=zero(T),
                          stage_stress=nothing) where {T}
    ns = n_stages(pop)

    total_attrition = zero(T)
    stage_totals = Vector{T}(undef, ns)
    stage_dds = Vector{T}(undef, ns)

    stress_vec = stage_stress === nothing ? fill(zero(T), ns) : T.(stage_stress)
    length(stress_vec) == ns ||
        throw(DimensionMismatch("stage_stress must have length $(ns)"))
    all(zero(T) .<= stress_vec .<= one(T)) ||
        throw(ArgumentError("stage_stress values must lie in [0, 1]"))

    stage_inflow = T(inflow)  # External inflow to first stage (set by reproduction)
    maturation = zero(T)

    for j in 1:ns
        stage = pop.stages[j]
        dd = degree_days(stage.dev_rate, weather_day.T_mean)
        stage_dds[j] = dd

        result = step_delay!(stage.delay, dd, stage_inflow;
                             μ=stage.μ, stress=stress_vec[j])

        total_attrition += result.attrition
        stage_inflow = result.outflow  # Maturation feeds into next stage
        stage_totals[j] = delay_total(stage.delay)
    end

    maturation = stage_inflow  # Outflow from last stage

    return (degree_days=stage_dds[1], stage_degree_days=stage_dds,
            maturation=maturation,
            total_attrition=total_attrition, stage_totals=stage_totals,
            stage_stress=stress_vec)
end

"""
    step_system!(pop, weather_day, approach=nothing; inflow=0.0, kwargs...)

Advance a population using an explicit approach layer. Legacy stepping uses the
plain distributed-delay update; BDF/MP/hybrid approaches can inject stress or
resource-allocation effects before the stage update.
"""
function step_system!(pop::Population{T}, weather_day::DailyWeather,
                      ::Nothing=nothing; inflow::Real=zero(T), kwargs...) where {T}
    result = step_population!(pop, weather_day; inflow=inflow)
    return (; result..., approach_family=:legacy)
end

function _population_demands(pop::Population{T}) where {T}
    T[delay_total(stage.delay) for stage in pop.stages]
end

function _population_labels(pop::Population)
    [stage.name for stage in pop.stages]
end

function _allocation_stress(pool::MetabolicPool{T}) where {T}
    alloc = allocate(pool)
    stress = similar(alloc)
    for i in eachindex(alloc)
        demand = pool.demands[i]
        stress[i] = demand > 0 ? max(zero(T), one(T) - alloc[i] / demand) : zero(T)
    end
    return alloc, stress
end

function _bdf_budget(pop::Population{T}, weather_day::DailyWeather,
                     model::BiodemographicFunctions) where {T}
    stage_demands = _population_demands(pop)
    total_demand = sum(stage_demands)
    gross_supply = acquire(model.acquisition, weather_day.radiation, total_demand)
    resp = respiration_rate(model.respiration, weather_day.T_mean) * total_demand
    net_supply = max(zero(T), T(gross_supply - resp))
    return (gross_supply=T(gross_supply), respiration=T(resp),
            net_supply=net_supply, demand=total_demand)
end

function step_system!(pop::Population{T}, weather_day::DailyWeather,
                      model::BiodemographicFunctions;
                      inflow::Real=zero(T), kwargs...) where {T}
    budget = _bdf_budget(pop, weather_day, model)
    deficit = budget.demand > 0 ? max(zero(T), one(T) - budget.net_supply / budget.demand) : zero(T)
    ratio = budget.demand > 0 ? min(one(T), budget.net_supply / budget.demand) : one(T)
    stage_stress = fill(deficit, n_stages(pop))
    result = step_population!(pop, weather_day; inflow=inflow, stage_stress=stage_stress)
    return (; result..., gross_supply=budget.gross_supply, respiration=budget.respiration,
            net_supply=budget.net_supply, demand=budget.demand,
            supply_demand=ratio, supply_demand_ratio=ratio,
            approach_family=:bdf, approach_label=model.label)
end

function _instantiate_pool(template::MetabolicPool{T}, pop::Population{T};
                           supply::Real=template.supply) where {T}
    stage_totals = _population_demands(pop)
    if length(template.demands) == length(stage_totals)
        demands = T.(template.demands) .* stage_totals
        labels = template.labels
    else
        demands = stage_totals
        labels = _population_labels(pop)
    end
    return MetabolicPool(T(supply), demands, labels)
end

function step_system!(pop::Population{T}, weather_day::DailyWeather,
                      pool::MetabolicPool{T};
                      inflow::Real=zero(T), kwargs...) where {T}
    realized_pool = _instantiate_pool(pool, pop)
    alloc, stress = _allocation_stress(realized_pool)
    ratio = supply_demand_index(realized_pool)
    result = step_population!(pop, weather_day; inflow=inflow, stage_stress=stress)
    return (; result..., allocations=alloc, demand=sum(realized_pool.demands),
            supply=realized_pool.supply,
            supply_demand=ratio, supply_demand_ratio=ratio,
            approach_family=:mp)
end

function step_system!(pop::Population{T}, weather_day::DailyWeather,
                      model::CoupledPBDMModel;
                      inflow::Real=zero(T), kwargs...) where {T}
    budget = _bdf_budget(pop, weather_day, biodemography(model))
    realized_pool = _instantiate_pool(allocation_model(model), pop; supply=budget.net_supply)
    alloc, stress = _allocation_stress(realized_pool)
    ratio = supply_demand_index(realized_pool)
    result = step_population!(pop, weather_day; inflow=inflow, stage_stress=stress)
    return (; result..., gross_supply=budget.gross_supply, respiration=budget.respiration,
            net_supply=budget.net_supply, demand=sum(realized_pool.demands),
            allocations=alloc,
            supply_demand=ratio, supply_demand_ratio=ratio,
            approach_family=:hybrid, approach_label=model.label)
end

function step_system!(pop::Population{T}, weather_day::DailyWeather,
                      ::AbstractBiodemographicApproach;
                      inflow::Real=zero(T), kwargs...) where {T}
    throw(ArgumentError("No step_system! implementation exists for this biodemographic approach type. Define a specialized method for $(typeof(pop)) and the concrete approach subtype."))
end

function step_system!(pop::Population{T}, weather_day::DailyWeather,
                      ::AbstractAllocationApproach;
                      inflow::Real=zero(T), kwargs...) where {T}
    throw(ArgumentError("No step_system! implementation exists for this allocation approach type. Define a specialized method for $(typeof(pop)) and the concrete approach subtype."))
end
