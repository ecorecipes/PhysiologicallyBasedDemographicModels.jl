"""
DelayDiffEq extension for PhysiologicallyBasedDemographicModels.jl

Provides `solve_delay` for `DelayPBDMProblem` using delay differential
equations. Each life stage is a single state variable with an explicit
gamma-distributed delay for maturation, solved via `MethodOfSteps`.

Activated automatically when `using DelayDiffEq` is called alongside
the main package.
"""
module DelayDiffEqExt

using PhysiologicallyBasedDemographicModels
using DelayDiffEq
using OrdinaryDiffEq: Tsit5

import PhysiologicallyBasedDemographicModels:
    solve_delay,
    DelayPBDMProblem, ContinuousPBDMSolution,
    ContinuousSpecies, ContinuousTrophicLink,
    _get_temperature,
    FraserGilbertResponse, HollingTypeII,
    development_rate, respiration_rate, acquire

# ============================================================================
# DelayPBDMProblem → DDEProblem
# ============================================================================

"""
    _build_dde_rhs(prob::DelayPBDMProblem)

Build the DDE right-hand-side function.

State variables: one biomass per species per stage.
State layout: species 1 stage 1, species 1 stage 2, ..., species N stage M.

Maturation from stage j to stage j+1 uses a discrete-lag approximation of
the Erlang-k distributed delay. With shape parameter k and mean delay τ,
the maturation outflow is approximated as:

    maturation_j(t) ≈ (k/τ)^k · h(t - τ) · τ / Γ(k) · W_j(t)

For integer k (Erlang), this simplifies to the k-fold convolution of
exponential waiting times. We approximate by using k discrete taps
equally spaced across the delay distribution.
"""
function _build_dde_rhs(prob::DelayPBDMProblem)
    species = prob.species
    links = prob.links
    T_forcing = prob.T_forcing
    n_sp = length(species)

    # State layout: species[i].n_stages variables per species
    species_offsets = Int[]
    stage_indices = Vector{Vector{Int}}()
    idx = 1
    for sp in species
        push!(species_offsets, idx)
        si = Int[]
        for j in 1:sp.n_stages
            push!(si, idx)
            idx += 1
        end
        push!(stage_indices, si)
    end

    # Name → index
    name_to_idx = Dict{Symbol, Int}(sp.name => i for (i, sp) in enumerate(species))

    # Trophic links
    consumers_of = Dict{Int, Vector{Int}}()
    resources_of = Dict{Int, Vector{Int}}()
    for link in links
        ci = get(name_to_idx, link.consumer, 0)
        ri = get(name_to_idx, link.resource, 0)
        ci > 0 && ri > 0 || continue
        push!(get!(consumers_of, ri, Int[]), ci)
        push!(get!(resources_of, ci, Int[]), ri)
    end

    # Pre-compute delay lags (in calendar days) for each stage
    # For each stage, the mean lag τ is in degree-days; convert at reference T
    # For DDE, we use a single representative lag per stage.
    # At constant temperature T, lag_days = τ / degree_days_per_day(T).

    function dde_rhs!(du, u, h, p, t)
        T_c = _get_temperature(T_forcing, t)

        # Compute totals per species
        totals = zeros(n_sp)
        for i in 1:n_sp
            for j in 1:species[i].n_stages
                totals[i] += max(0.0, u[stage_indices[i][j]])
            end
        end

        # Trophic acquisition/consumption
        acquisition = zeros(n_sp)
        consumed_from = zeros(n_sp)
        for i in 1:n_sp
            sp = species[i]
            if haskey(resources_of, i)
                demand = sp.demand_rate * totals[i]
                for ri in resources_of[i]
                    supply = totals[ri]
                    acq = _dde_acquire(sp.fr, supply, demand)
                    acquisition[i] += sp.conversion_efficiency * acq
                    consumed_from[ri] += acq
                end
            end
        end

        for i in 1:n_sp
            sp = species[i]
            total_i = totals[i]
            pred_rate = total_i > 0 ? consumed_from[i] / total_i : 0.0
            R_i = respiration_rate(sp.resp, T_c)

            for j in 1:sp.n_stages
                si = stage_indices[i][j]
                W = u[si]
                g_j = development_rate(sp.dev_rate[j], T_c)

                # Maturation uses aggregate rate (1/τ)*g, not the per-substage
                # rate (k/τ)*g from the linear chain trick.  k controls delay
                # shape (Erlang order) but the DDE already handles the delay
                # explicitly via the history function.
                aggregate_rate = g_j > 0 ? g_j / sp.τ[j] : 0.0

                # Lag in calendar days for this stage
                lag_days = g_j > 0 ? sp.τ[j] / g_j : 1e6
                lag_days = min(lag_days, 1e4)

                # Maturation outflow — only for non-terminal stages
                if j < sp.n_stages
                    # Use delayed state for outflow (consistent with DDE paradigm)
                    h_self = h(p, t - lag_days)
                    W_delayed = max(0.0, h_self[si])
                    mat_out = aggregate_rate * W_delayed
                else
                    # Terminal stage: no maturation exit (biomass stays until
                    # consumed by mortality, predation, or respiration)
                    mat_out = 0.0
                end

                # Maturation inflow to this stage
                if j == 1
                    if sp.intrinsic_rate > 0
                        # Basal resource: logistic growth
                        inflow = sp.intrinsic_rate * total_i * (1.0 - total_i / sp.carrying_capacity)
                        inflow = max(0.0, inflow)
                    else
                        # Consumer: acquisition
                        acq_rate = total_i > 0 ? acquisition[i] / total_i : 0.0
                        inflow = acq_rate * total_i
                    end
                else
                    # Delayed maturation from previous stage
                    prev_si = stage_indices[i][j-1]
                    prev_g = development_rate(sp.dev_rate[j-1], T_c)
                    prev_lag = prev_g > 0 ? sp.τ[j-1] / prev_g : 1e6
                    prev_lag = min(prev_lag, 1e4)

                    h_val = h(p, t - prev_lag)
                    W_prev_delayed = max(0.0, h_val[prev_si])
                    prev_agg_rate = prev_g > 0 ? prev_g / sp.τ[j-1] : 0.0
                    inflow = prev_agg_rate * W_prev_delayed
                end

                # dW/dt = inflow - outflow - mortality - predation - respiration
                du[si] = inflow - mat_out - sp.μ[j] * g_j * W - pred_rate * W - R_i * W
            end
        end
    end

    # Collect all distinct lags needed
    function compute_lags(T_c)
        lags = Float64[]
        for sp in species
            for j in 1:sp.n_stages
                g_j = development_rate(sp.dev_rate[j], T_c)
                lag = g_j > 0 ? sp.τ[j] / g_j : 1e6
                push!(lags, min(lag, 1e4))
            end
        end
        return unique(sort(lags))
    end

    return dde_rhs!, compute_lags
end

function _dde_acquire(fr::FraserGilbertResponse, supply, demand)
    (demand <= 0 || supply <= 0) && return 0.0
    return demand * (1.0 - exp(-fr.a * supply / demand))
end

function _dde_acquire(fr, supply, demand)
    return acquire(fr, supply, demand)
end

"""
    PhysiologicallyBasedDemographicModels.solve_delay(
        prob::DelayPBDMProblem; alg=nothing, kwargs...)

Solve a delay-differential-equation PBDM problem.

# Arguments
- `prob`: The delay PBDM problem
- `alg`: DDE solver algorithm (default: `MethodOfSteps(Tsit5())`)
- `kwargs...`: Passed to DelayDiffEq.solve

# Returns
`ContinuousPBDMSolution` with time series and species trajectories.
"""
function PhysiologicallyBasedDemographicModels.solve_delay(
        prob::DelayPBDMProblem; alg=nothing, kwargs...)

    dde_rhs!, compute_lags = _build_dde_rhs(prob)

    # Compute representative lags at a reference temperature
    T_ref = _get_temperature(prob.T_forcing, prob.tspan[1])
    lags = compute_lags(T_ref)

    # Filter out extremely large lags (non-developing stages)
    lags = filter(l -> l < 1e4, lags)

    # Build DDE problem
    h_fn = prob.h0
    dde_prob = DDEProblem(dde_rhs!, prob.u0, h_fn, prob.tspan, prob.p;
                          constant_lags=lags)

    if alg === nothing
        alg = MethodOfSteps(Tsit5())
    end

    sol = DelayDiffEq.solve(dde_prob, alg; kwargs...)

    # Convert to ContinuousPBDMSolution
    t_vec = Float64.(sol.t)
    n_states = length(prob.u0)
    n_t = length(t_vec)

    u_mat = zeros(n_states, n_t)
    for (j, ui) in enumerate(sol.u)
        u_mat[:, j] = ui
    end

    # Build species ranges (one variable per stage)
    sp_ranges = Dict{Symbol, UnitRange{Int}}()
    idx = 1
    for sp in prob.species
        sp_ranges[sp.name] = idx:(idx + sp.n_stages - 1)
        idx += sp.n_stages
    end

    sp_names = [sp.name for sp in prob.species]
    retcode = Symbol(string(sol.retcode))

    return ContinuousPBDMSolution{Float64, typeof(sol)}(
        t_vec, u_mat, sp_names, sp_ranges, retcode, sol)
end

end # module DelayDiffEqExt
