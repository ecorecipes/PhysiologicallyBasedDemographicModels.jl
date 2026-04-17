"""
JuMP extension for PhysiologicallyBasedDemographicModels.jl

Provides `optimize_management(prob::PBDMControlProblem)` using JuMP
with Euler-discretized PBDM dynamics, following the pattern established
in EpiPolicies (Frost & Montes 2025).

Activated automatically when `using JuMP` is called alongside the
main package.
"""
module JuMPExt

using PhysiologicallyBasedDemographicModels
using JuMP

import PhysiologicallyBasedDemographicModels: optimize_management,
    PBDMControlProblem, ManagementSolution,
    TrophicLevel, PesticideControl, BiologicalReleaseControl, HarvestControl,
    MinimizeDamage, MaximizeProfit,
    acquire, respiration_rate

"""
    optimize_management(prob::PBDMControlProblem;
                        solver=nothing, silent=true, max_iter=5000)

Solve an optimal management problem using nonlinear programming.

Euler-discretizes the PBDM trophic dynamics as JuMP nonlinear constraints
and optimizes control trajectories using the specified solver (defaults
to Ipopt if available).

# Arguments
- `prob`: The control problem specification
- `solver`: JuMP-compatible optimizer (e.g., `Ipopt.Optimizer`)
- `silent`: Suppress solver output (default: true)
- `max_iter`: Maximum solver iterations (default: 5000)

# Returns
- `ManagementSolution` with optimal state and control trajectories
"""
function PhysiologicallyBasedDemographicModels.optimize_management(
        prob::PBDMControlProblem;
        solver=nothing, silent::Bool=true, max_iter::Int=5000)

    if solver === nothing
        error("No solver specified. Use e.g. `optimize_management(prob; solver=Ipopt.Optimizer)`")
    end

    levels = prob.levels
    n_levels = length(levels)
    n_controls = length(prob.controls)
    dt = prob.dt
    T_c = prob.T_celsius
    n_steps = Int(ceil((prob.tspan[2] - prob.tspan[1]) / dt))
    ts = [prob.tspan[1] + i * dt for i in 0:n_steps]

    # Precompute respiration rates at constant temperature
    R = [respiration_rate(lev.resp, T_c) for lev in levels]

    # Build JuMP model
    model = Model(solver)
    set_optimizer_attribute(model, "max_iter", max_iter)
    silent && set_silent(model)

    # State variables: biomass for each level at each timestep
    @variable(model, M[1:n_levels, 1:(n_steps+1)] >= 0)

    # Control variables: one per action at each timestep
    ctrl_vars = []
    for (j, action) in enumerate(prob.controls)
        max_val = if action isa PesticideControl
            action.max_rate
        elseif action isa BiologicalReleaseControl
            action.max_rate
        elseif action isa HarvestControl
            action.max_rate
        else
            1.0
        end
        cv = @variable(model, [1:(n_steps+1)], lower_bound=0, upper_bound=max_val)
        push!(ctrl_vars, cv)
    end

    # Initial conditions
    for i in 1:n_levels
        @constraint(model, M[i, 1] == prob.u0[i])
    end

    # Euler-discretized dynamics as constraints
    for t_idx in 1:n_steps
        for i in 1:n_levels
            lev = levels[i]
            Mi = M[i, t_idx]

            if i == 1
                # Basal resource: logistic growth
                growth_expr = lev.intrinsic_rate * Mi * (1 - Mi / lev.carrying_capacity)
            else
                # Consumer: demand-driven acquisition from level below
                # For JuMP compatibility, inline the Fraser-Gilbert equation
                M_resource = M[i-1, t_idx]
                demand = lev.demand_rate * Mi
                # acquisition = demand * (1 - exp(-a * supply / demand))
                # = lev.demand_rate * Mi * (1 - exp(-a * M[i-1] / (lev.demand_rate * Mi)))
                a_val = _get_apparency(lev.fr)
                # Guard against division by zero: when Mi ≈ 0, acquisition ≈ 0
                growth_expr = lev.conversion_efficiency * lev.demand_rate * Mi *
                    (1 - exp(-a_val * M_resource / (lev.demand_rate * Mi + 1e-10)))
            end

            # Respiration
            resp_expr = R[i] * Mi

            # Consumption by level above
            if i < n_levels
                lev_above = levels[i+1]
                a_above = _get_apparency(lev_above.fr)
                consumption_expr = lev_above.demand_rate * M[i+1, t_idx] *
                    (1 - exp(-a_above * Mi / (lev_above.demand_rate * M[i+1, t_idx] + 1e-10)))
            else
                consumption_expr = 0.0
            end

            # Apply control actions affecting this level
            control_effect = AffExpr(0.0)
            for (j, action) in enumerate(prob.controls)
                u_j = ctrl_vars[j][t_idx]
                if action isa PesticideControl && action.target == i
                    control_effect += action.efficacy * u_j * Mi
                elseif action isa BiologicalReleaseControl && action.target == i
                    control_effect -= u_j  # Release adds biomass (subtracted from loss)
                elseif action isa HarvestControl && action.target == i
                    control_effect += u_j * Mi
                end
            end

            # Euler step: M[i, t+1] = M[i, t] + dt * (growth - resp - consumption - control)
            @constraint(model,
                M[i, t_idx+1] == Mi + dt * (growth_expr - resp_expr - consumption_expr) - dt * control_effect)
        end
    end

    # Objective function
    obj = prob.objective
    if obj isa MinimizeDamage
        obj_expr = AffExpr(0.0)
        for t_idx in 1:(n_steps+1)
            for i in 1:n_levels
                if i <= length(obj.damage_weights) && obj.damage_weights[i] > 0
                    obj_expr += obj.damage_weights[i] * M[i, t_idx]
                end
            end
            for (j, action) in enumerate(prob.controls)
                cw = _get_cost_weight(action) * obj.control_cost_weight
                obj_expr += cw * ctrl_vars[j][t_idx]^2
            end
        end
        @objective(model, Min, obj_expr)

    elseif obj isa MaximizeProfit
        obj_expr = AffExpr(0.0)
        for t_idx in 1:(n_steps+1)
            discount = 1.0 / (1.0 + obj.discount_rate)^((t_idx-1) * dt)
            obj_expr += discount * obj.price_per_unit * M[obj.resource_level, t_idx]
            for (j, action) in enumerate(prob.controls)
                cw = _get_cost_weight(action) * obj.control_cost_weight
                obj_expr -= discount * cw * ctrl_vars[j][t_idx]^2
            end
        end
        @objective(model, Max, obj_expr)
    end

    # Solve
    optimize!(model)

    # Extract results
    status = termination_status(model)
    status_sym = Symbol(string(status))

    state_matrix = zeros(n_levels, n_steps + 1)
    for i in 1:n_levels, t_idx in 1:(n_steps+1)
        state_matrix[i, t_idx] = value(M[i, t_idx])
    end

    control_matrix = zeros(n_controls, n_steps + 1)
    for j in 1:n_controls, t_idx in 1:(n_steps+1)
        control_matrix[j, t_idx] = value(ctrl_vars[j][t_idx])
    end

    obj_val = objective_value(model)

    return ManagementSolution{Float64}(
        ts, state_matrix, control_matrix, obj_val,
        [lev.name for lev in levels],
        [_get_name(a) for a in prob.controls],
        status_sym)
end

# Helper functions to extract fields from different action types
_get_apparency(fr::PhysiologicallyBasedDemographicModels.FraserGilbertResponse) = fr.a
_get_apparency(fr) = error("JuMP extension currently supports only FraserGilbertResponse")

_get_cost_weight(a::PesticideControl) = a.cost_weight
_get_cost_weight(a::BiologicalReleaseControl) = a.cost_weight
_get_cost_weight(a::HarvestControl) = a.cost_weight

_get_name(a::PesticideControl) = a.name
_get_name(a::BiologicalReleaseControl) = a.name
_get_name(a::HarvestControl) = a.name

end # module JuMPExt
