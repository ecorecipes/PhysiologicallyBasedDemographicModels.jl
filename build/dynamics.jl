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
    k = delay.k
    τ = delay.τ
    W = delay.W

    # Fractional development this step: proportion of stage completed
    # Rate per substage per degree-day = k / τ
    rate = (k / τ) * dd

    # Combined attrition per degree-day
    attrition_rate = (μ + stress) * dd

    total_attrition = zero(T)
    outflow = zero(T)

    # Process substages from last to first (reverse to avoid overwriting)
    prev_flow = T(inflow)

    for i in 1:k
        # Attrition from this substage
        loss = attrition_rate * W[i]
        total_attrition += loss

        # Flow out of this substage
        flow_out = rate * W[i]

        # Update: gain from previous substage, lose to next and to attrition
        W[i] = W[i] + prev_flow - flow_out - loss

        # Clamp to non-negative
        if W[i] < 0
            W[i] = zero(T)
        end

        prev_flow = flow_out
    end

    outflow = prev_flow  # Flow out of the last substage

    return (outflow=outflow, attrition=total_attrition)
end

"""
    step_population!(pop, weather_day)

Advance a single-species Population by one calendar day.
Computes degree-days from the weather, then steps each life stage's
distributed delay. Maturation outflow from stage j becomes inflow to stage j+1.
The last stage's outflow is returned as the population's reproductive output.

# Returns
Named tuple with fields:
- `degree_days`: DD accumulated this day
- `maturation`: outflow from the final life stage
- `total_attrition`: total mortality across all stages
- `stage_totals`: vector of total population per stage
"""
function step_population!(pop::Population{T}, weather_day::DailyWeather) where {T}
    ns = n_stages(pop)

    # Compute degree-days from first stage's development rate
    # (in practice each stage may have its own dev_rate)
    total_dd = zero(T)
    total_attrition = zero(T)
    stage_totals = Vector{T}(undef, ns)

    inflow = zero(T)  # External inflow to first stage (set by reproduction)
    maturation = zero(T)

    for j in 1:ns
        stage = pop.stages[j]
        dd = degree_days(stage.dev_rate, weather_day.T_mean)

        if j == 1
            total_dd = dd
        end

        result = step_delay!(stage.delay, dd, inflow;
                             μ=stage.μ)

        total_attrition += result.attrition
        inflow = result.outflow  # Maturation feeds into next stage
        stage_totals[j] = delay_total(stage.delay)
    end

    maturation = inflow  # Outflow from last stage

    return (degree_days=total_dd, maturation=maturation,
            total_attrition=total_attrition, stage_totals=stage_totals)
end
