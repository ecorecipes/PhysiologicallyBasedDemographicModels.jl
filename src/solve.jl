"""
    PBDMSolution

Result of solving a PBDMProblem.

Follows the AbstractProjectionSolution interface from StructuredPopulationCore.jl,
extended with PBDM-specific outputs (degree-day accumulation, supply/demand).

# Fields
- `t`: calendar day labels for each completed simulation day
- `u`: population state after each simulated day
- `kernel_matrices`: nothing (PBDMs don't use matrix kernels)
- `eigenanalysis`: nothing (PBDMs are non-linear; no eigendecomposition)
- `retcode`: :Success or :Failure
- `lambdas`: per-step growth rates
- `degree_days`: daily degree-days for the first stage (legacy summary)
- `stage_degree_days`: matrix (n_stages × n_days) of per-stage daily degree-days
- `stage_totals`: matrix (n_stages × n_days) of per-stage populations
- `maturation`: per-day maturation output from the final stage
"""
struct PBDMSolution{T<:Real} <: AbstractProjectionSolution
    t::Vector{Int}
    u::Vector{Vector{T}}
    kernel_matrices::Nothing
    eigenanalysis::Nothing
    retcode::Symbol
    lambdas::Vector{T}
    degree_days::Vector{T}
    stage_degree_days::Matrix{T}
    stage_totals::Matrix{T}
    maturation::Vector{T}
end

function Base.show(io::IO, sol::PBDMSolution)
    nt = length(sol.t)
    ns = size(sol.stage_totals, 1)
    print(io, "PBDMSolution($nt days, $ns stages, retcode=$(sol.retcode))")
end

# --- Main solve dispatch ---

function CommonSolve.solve(prob::PBDMProblem, alg::DirectIteration=DirectIteration();
                           kwargs...)
    _solve(prob.structure, prob.density, prob.stochasticity, prob.approach, prob, alg; kwargs...)
end

function CommonSolve.solve(prob::PBDMProblem, ::EigenAnalysis; kwargs...)
    throw(ArgumentError("EigenAnalysis is not supported for nonlinear PBDM problems. Use DirectIteration() instead."))
end

# --- Single-species, density-independent, deterministic ---

function _solve(::SingleSpeciesPBDM, ::DensityIndependent, ::Deterministic,
                approach,
                prob::PBDMProblem, ::DirectIteration; kwargs...)
    pop = prob.populations
    weather = prob.weather
    t0, tf = prob.tspan
    n_days = tf - t0 + 1

    T = eltype(pop.stages[1].delay.W)
    ns = n_stages(pop)

    # Pre-allocate outputs
    u = Vector{Vector{T}}(undef, n_days)
    lambdas = Vector{T}(undef, n_days)
    dd_accum = Vector{T}(undef, n_days)
    stage_dds = Matrix{T}(undef, ns, n_days)
    stage_tots = Matrix{T}(undef, ns, n_days)
    matur = Vector{T}(undef, n_days)

    prev_totals = T[delay_total(s.delay) for s in pop.stages]

    for d in 1:n_days
        day = t0 + d - 1
        w = get_weather(weather, day)
        result = step_system!(pop, w, approach)

        dd_accum[d] = result.degree_days
        stage_dds[:, d] = result.stage_degree_days
        matur[d] = result.maturation

        stage_tots[:, d] = result.stage_totals

        total_now = sum(result.stage_totals)
        total_prev = sum(prev_totals)
        lambdas[d] = total_prev > 0 ? total_now / total_prev : zero(T)
        u[d] = copy(result.stage_totals)
        prev_totals = u[d]
    end

    ts = collect(t0:tf)
    return PBDMSolution(ts, u, nothing, nothing, :Success,
                        lambdas, dd_accum, stage_dds, stage_tots, matur)
end

# --- Single-species, density-dependent, deterministic ---

function _solve(::SingleSpeciesPBDM, ::DensityDependent, ::Deterministic,
                approach,
                prob::PBDMProblem, ::DirectIteration;
                reproduction_fn=nothing, kwargs...)
    pop = prob.populations
    weather = prob.weather
    t0, tf = prob.tspan
    n_days = tf - t0 + 1

    T = eltype(pop.stages[1].delay.W)
    ns = n_stages(pop)

    u = Vector{Vector{T}}(undef, n_days)
    lambdas = Vector{T}(undef, n_days)
    dd_accum = Vector{T}(undef, n_days)
    stage_dds = Matrix{T}(undef, ns, n_days)
    stage_tots = Matrix{T}(undef, ns, n_days)
    matur = Vector{T}(undef, n_days)

    prev_totals = T[delay_total(s.delay) for s in pop.stages]

    for d in 1:n_days
        day = t0 + d - 1
        w = get_weather(weather, day)

        offspring = reproduction_fn !== nothing ? reproduction_fn(pop, w, prob.p, day) : zero(T)

        result = step_system!(pop, w, approach; inflow=offspring)

        dd_accum[d] = result.degree_days
        stage_dds[:, d] = result.stage_degree_days
        matur[d] = result.maturation

        stage_tots[:, d] = result.stage_totals

        total_now = sum(result.stage_totals)
        total_prev = sum(prev_totals)
        lambdas[d] = total_prev > 0 ? total_now / total_prev : zero(T)
        u[d] = copy(result.stage_totals)
        prev_totals = u[d]
    end

    ts = collect(t0:tf)
    return PBDMSolution(ts, u, nothing, nothing, :Success,
                        lambdas, dd_accum, stage_dds, stage_tots, matur)
end

function _solve(::MultiSpeciesPBDM, ::AbstractDensityDependence, ::AbstractStochasticity,
                approach,
                prob::PBDMProblem, ::DirectIteration; kwargs...)
    throw(ArgumentError("Direct solve support for MultiSpeciesPBDM is not implemented yet. Solve each population separately or use an explicit custom interaction simulation."))
end
