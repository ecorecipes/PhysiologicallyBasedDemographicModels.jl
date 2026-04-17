"""
Analysis functions for PBDM solutions.

Extends StructuredPopulationCore analysis where applicable and provides
PBDM-specific analyses.
"""

"""
    cumulative_degree_days(sol::PBDMSolution; stage_idx=1)

Return cumulative degree-day accumulation over the simulation for a selected
life stage.
"""
function cumulative_degree_days(sol::PBDMSolution; stage_idx::Int=1)
    return cumsum(stage_degree_days(sol, stage_idx))
end

"""
    stage_degree_days(sol::PBDMSolution, stage_idx::Int)

Extract the daily degree-day trajectory for a single stage.
"""
function stage_degree_days(sol::PBDMSolution, stage_idx::Int)
    return sol.stage_degree_days[stage_idx, :]
end

"""
    stage_trajectory(sol::PBDMSolution, stage_idx::Int)

Extract the population trajectory for a single stage.
Returns a vector of length `length(sol.t)`.
"""
function stage_trajectory(sol::PBDMSolution, stage_idx::Int)
    return sol.stage_totals[stage_idx, :]
end

"""
    total_population(sol::PBDMSolution)

Compute total population at each time step (sum across all stages).
"""
function total_population(sol::PBDMSolution)
    return [sum(u) for u in sol.u]
end

"""
    net_growth_rate(sol::PBDMSolution)

Compute mean daily growth rate λ over the simulation.
Only considers days with positive lambda.
"""
function net_growth_rate(sol::PBDMSolution)
    pos = filter(l -> l > 0, sol.lambdas)
    isempty(pos) && return 0.0
    return exp(mean(log.(pos)))
end

"""
    phenology(sol::PBDMSolution; threshold=0.5)

Estimate phenological event timing from maturation output.
Returns the calendar day when cumulative maturation first exceeds
`threshold` fraction of total maturation.
"""
function phenology(sol::PBDMSolution; threshold::Real=0.5)
    cum = cumsum(sol.maturation)
    total = cum[end]
    total <= 0 && return nothing
    idx = findfirst(c -> c >= threshold * total, cum)
    idx === nothing && return nothing
    return sol.t[idx]
end
