"""
    PBDMSolution

Result of solving a PBDMProblem.

Follows the AbstractProjectionSolution interface from StructuredPopulationCore.jl,
extended with PBDM-specific outputs (degree-day accumulation, supply/demand).

# Fields
- `t`: calendar day time steps
- `u`: population state at each time step (vector of stage totals per day)
- `kernel_matrices`: nothing (PBDMs don't use matrix kernels)
- `eigenanalysis`: nothing (PBDMs are non-linear; no eigendecomposition)
- `retcode`: :Success or :Failure
- `lambdas`: per-step growth rates
- `degree_days`: cumulative degree-day accumulation
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
    _solve(prob.structure, prob.density, prob.stochasticity, prob, alg; kwargs...)
end

# --- Single-species, density-independent, deterministic ---

function _solve(::SingleSpeciesPBDM, ::DensityIndependent, ::Deterministic,
                prob::PBDMProblem, ::DirectIteration; kwargs...)
    pop = prob.populations
    weather = prob.weather
    t0, tf = prob.tspan
    n_days = tf - t0

    T = eltype(pop.stages[1].delay.W)
    ns = n_stages(pop)

    # Pre-allocate outputs
    u = Vector{Vector{T}}(undef, n_days + 1)
    lambdas = Vector{T}(undef, n_days)
    dd_accum = Vector{T}(undef, n_days)
    stage_tots = Matrix{T}(undef, ns, n_days + 1)
    matur = Vector{T}(undef, n_days)

    # Initial state
    init_totals = T[delay_total(s.delay) for s in pop.stages]
    u[1] = copy(init_totals)
    for j in 1:ns
        stage_tots[j, 1] = init_totals[j]
    end

    for d in 1:n_days
        day = t0 + d - 1
        w = get_weather(weather, day)
        result = step_population!(pop, w)

        dd_accum[d] = result.degree_days
        matur[d] = result.maturation

        # Compute per-stage totals
        for j in 1:ns
            stage_tots[j, d + 1] = result.stage_totals[j]
        end

        total_now = sum(result.stage_totals)
        total_prev = sum(u[d])
        lambdas[d] = total_prev > 0 ? total_now / total_prev : zero(T)
        u[d + 1] = copy(result.stage_totals)
    end

    ts = collect(t0:tf)
    return PBDMSolution(ts, u, nothing, nothing, :Success,
                        lambdas, dd_accum, stage_tots, matur)
end

# --- Single-species, density-dependent, deterministic ---

function _solve(::SingleSpeciesPBDM, ::DensityDependent, ::Deterministic,
                prob::PBDMProblem, ::DirectIteration;
                reproduction_fn=nothing, kwargs...)
    pop = prob.populations
    weather = prob.weather
    t0, tf = prob.tspan
    n_days = tf - t0

    T = eltype(pop.stages[1].delay.W)
    ns = n_stages(pop)

    u = Vector{Vector{T}}(undef, n_days + 1)
    lambdas = Vector{T}(undef, n_days)
    dd_accum = Vector{T}(undef, n_days)
    stage_tots = Matrix{T}(undef, ns, n_days + 1)
    matur = Vector{T}(undef, n_days)

    init_totals = T[delay_total(s.delay) for s in pop.stages]
    u[1] = copy(init_totals)
    for j in 1:ns
        stage_tots[j, 1] = init_totals[j]
    end

    for d in 1:n_days
        day = t0 + d - 1
        w = get_weather(weather, day)

        # If reproduction function is provided, inject offspring into stage 1
        if reproduction_fn !== nothing
            offspring = reproduction_fn(pop, w, prob.p, day)
            pop.stages[1].delay.W[1] += offspring
        end

        result = step_population!(pop, w)

        dd_accum[d] = result.degree_days
        matur[d] = result.maturation

        for j in 1:ns
            stage_tots[j, d + 1] = result.stage_totals[j]
        end

        total_now = sum(result.stage_totals)
        total_prev = sum(u[d])
        lambdas[d] = total_prev > 0 ? total_now / total_prev : zero(T)
        u[d + 1] = copy(result.stage_totals)
    end

    ts = collect(t0:tf)
    return PBDMSolution(ts, u, nothing, nothing, :Success,
                        lambdas, dd_accum, stage_tots, matur)
end
