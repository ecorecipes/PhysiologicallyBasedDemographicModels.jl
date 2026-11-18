# Resilience of PBDMs through DynamicalResilience.jl and ResilienceCore.jl.
#
# - State vectors: the complete numerical state of a population, a system or a problem as one
#   vector, and back.
# - An integrator for `PBDMProblem` (`init`, `step!`, `reinit!`, `set_u!`, …), so that
#   DynamicalResilience's simulation-based metrics (`recovery`, …) run on PBDMs unchanged.
# - The annual map of a model under a repeated weather year, its fixed point (the annual
#   cycle) and its Jacobian, whose eigenvalues are the cycle's Floquet multipliers.
# - `linearise` for the consumer–resource equilibria of `theory.jl`, and the same model as an
#   `ODEProblem` for DynamicalResilience's continuation along temperature.

# ============================================================================
# State vectors
# ============================================================================

_state_variable_names(sys::PopulationSystem) =
    sort!([name for (name, sv) in sys.state if sv isa Union{ScalarState, ArrayState}])

"""
    full_state(x) -> Vector{Float64}

The numerical state of `x` as one vector:

- a `Population`: the contents of every substage of every stage's distributed delay, stage
  by stage;
- a `BulkPopulation`: its value;
- a `PopulationSystem`: its components in `keys(sys)` order, then its `ScalarState` and
  `ArrayState` variables in order of name;
- a `PBDMProblem`: its populations.

Other state variables (`DictState`, `GenomeState`, `DiapauseState`, …) belong to the system
and are copied and restored with it, but are not part of the vector: they are neither perturbed
nor linearised. [`set_full_state!`](@ref) is the inverse and [`full_state_labels`](@ref) names the
entries.
"""
full_state(pop::Population) = Float64[x for s in pop.stages for x in s.delay.W]
full_state(bp::BulkPopulation) = Float64[bp.value[]]
function full_state(sys::PopulationSystem)
    u = Float64[]
    for name in keys(sys)
        append!(u, full_state(sys[name].population))
    end
    for name in _state_variable_names(sys)
        sv = sys.state[name]
        sv isa ScalarState ? push!(u, sv.value[]) : append!(u, sv.value)
    end
    return u
end
full_state(prob::PBDMProblem) = full_state(prob.populations)

_state_length(pop::Population) = sum(s.delay.k for s in pop.stages)
_state_length(::BulkPopulation) = 1
function _state_length(sys::PopulationSystem)
    n = sum(_state_length(sys[name].population) for name in keys(sys))
    for name in _state_variable_names(sys)
        sv = sys.state[name]
        n += sv isa ScalarState ? 1 : length(sv.value)
    end
    return n
end

function _check_state_length(x, u)
    n = _state_length(x)
    length(u) == n || throw(DimensionMismatch("the state of this $(nameof(typeof(x))) has $n entries, got $(length(u))"))
    return nothing
end

"""
    set_full_state!(x, u) -> x

Write the state vector `u`, laid out as by [`full_state`](@ref), into the population, system
or problem `x`. Values are stored as given: a negative entry is not clipped (perturbations used
to linearise a map can be negative), but a `BulkPopulation` above its carrying capacity is.
"""
function set_full_state!(pop::Population, u::AbstractVector{<:Real})
    _check_state_length(pop, u)
    i = 0
    for s in pop.stages, j in eachindex(s.delay.W)
        s.delay.W[j] = u[i += 1]
    end
    return pop
end
function set_full_state!(bp::BulkPopulation, u::AbstractVector{<:Real})
    _check_state_length(bp, u)
    bp.value[] = min(u[1], bp.K)
    return bp
end
function set_full_state!(sys::PopulationSystem, u::AbstractVector{<:Real})
    _check_state_length(sys, u)
    i = 0
    for name in keys(sys)
        pop = sys[name].population
        n = _state_length(pop)
        set_full_state!(pop, view(u, (i + 1):(i + n)))
        i += n
    end
    for name in _state_variable_names(sys)
        sv = sys.state[name]
        if sv isa ScalarState
            sv.value[] = u[i += 1]
        else
            n = length(sv.value)
            sv.value .= view(u, (i + 1):(i + n))
            i += n
        end
    end
    return sys
end
set_full_state!(prob::PBDMProblem, u::AbstractVector{<:Real}) = (set_full_state!(prob.populations, u); prob)

"""
    full_state_labels(x) -> Vector{String}

Names for the entries of [`full_state`](@ref)`(x)`: `"pest.larva[3]"` for the third
substage of stage `larva` of component (or population) `pest`, `"scale"` for a
`BulkPopulation`, `"state.temperature"` or `"state.pools[2]"` for state variables.
"""
full_state_labels(pop::Population) = [string(pop.name, ".", s.name, "[", j, "]") for s in pop.stages for j in 1:s.delay.k]
full_state_labels(bp::BulkPopulation) = [string(bp.name)]
function full_state_labels(sys::PopulationSystem)
    labels = String[]
    for name in keys(sys)
        pop = sys[name].population
        sub = full_state_labels(pop)
        pop isa BulkPopulation ? push!(labels, string(name)) :
            append!(labels, [string(name, chopprefix(l, string(pop.name)))
                             for l in sub])
    end
    for name in _state_variable_names(sys)
        sv = sys.state[name]
        sv isa ScalarState ? push!(labels, string("state.", name)) :
            append!(labels, [string("state.", name, "[", j, "]") for j in eachindex(sv.value)])
    end
    return labels
end
full_state_labels(prob::PBDMProblem) = full_state_labels(prob.populations)

# ============================================================================
# The integrator: SciML's integrator interface over daily steps
# ============================================================================

"""
    PBDMIntegrator

The integrator returned by `init(prob::PBDMProblem, DirectIteration())`. It steps a copy of the
problem's populations one day at a time, exactly as `solve` does, and implements the part of
SciML's integrator interface that DynamicalResilience.jl uses:

- `integ.t`: the last completed calendar day (`first(prob.tspan) - 1` before the first step);
- `integ.u`: [`full_state`](@ref) of the populations after day `integ.t`;
- `step!(integ)` advances one day; `step!(integ, dt, true)` advances `dt` days (`dt` must be a
  whole number of days);
- `reinit!(integ[, u0])` restarts from the problem's populations (or from the state `u0`) at
  `first(prob.tspan) - 1`. PBDMs are deterministic and weather is a function of the day, so a
  run replays exactly;
- `set_u!(integ, u)` writes a new state; or change `integ.u` in place and then call
  `derivative_discontinuity!(integ, true)`;
- `terminate!(integ)` stops further steps; `check_error(integ)` returns `ReturnCode.Failure`
  once a state entry is not finite, `ReturnCode.Terminated` after `terminate!`, and
  `ReturnCode.Success` otherwise;
- `integ.p` (the problem's parameters) may be reassigned between steps.

Covered: single-species problems (`SingleSpeciesPBDM`, density-independent or density-dependent,
deterministic; give a density-dependent problem's `reproduction_fn` to `init`, as to `solve`)
and coupled `PopulationSystem` problems (`MultiTypePBDM`, `MultiSpeciesPBDMNew`, …). Each day
runs PRE_EVENT callbacks, events, state variables, stress, stepping, rules and END_OF_DAY
callbacks; observables are not recorded. Events and rules must not keep their own mutable
state between days, or a replay will not be exact.

`solve(prob)` steps `prob.populations` in place, so `init` a problem before solving it, or
`remake` it, if the integrator should start from the original state.
"""
mutable struct PBDMIntegrator{P<:PBDMProblem, S, A}
    prob::P
    sys::S
    alg::A
    t::Int
    u::Vector{Float64}
    p::Any
    reproduction_fn::Any
    parts::Any
    comp_names::Vector{Symbol}
    retcode::SciMLBase.ReturnCode.T
end

function Base.setproperty!(integ::PBDMIntegrator, name::Symbol, value)
    if name === :p
        setfield!(integ, :p, value)
        integ.prob.populations isa PopulationSystem && setfield!(integ, :parts, _coupled_parts(value))
        return value
    end
    return setfield!(integ, name, convert(fieldtype(typeof(integ), name), value))
end

function Base.show(io::IO, integ::PBDMIntegrator)
    print(io, "PBDMIntegrator(day $(integ.t), $(length(integ.u)) state entries, $(integ.retcode))")
end

_check_integrable(prob::PBDMProblem) = _check_integrable(prob.structure, prob.density, prob.stochasticity, prob.populations)
_check_integrable(::SingleSpeciesPBDM, ::AbstractDensityDependence, ::Deterministic, ::Population) = nothing
_check_integrable(::MultiPopulationPBDM, ::DensityDependent, ::Deterministic, ::PopulationSystem) = nothing
_check_integrable(structure, density, stochasticity, populations) = throw(ArgumentError(
    "an integrator is available for deterministic single-species problems with a Population and for coupled " *
    "PopulationSystem problems; got $(nameof(typeof(structure))), $(nameof(typeof(density))), " *
    "$(nameof(typeof(stochasticity))) with $(typeof(populations))"))

"""
    init(prob::PBDMProblem, alg = DirectIteration(); reproduction_fn = nothing, seed = nothing) -> PBDMIntegrator

An integrator over the days of `prob`, starting from a copy of `prob.populations` (see
[`PBDMIntegrator`](@ref)). `seed` is accepted for compatibility with stochastic SciML problems
and ignored: the problem types covered are deterministic.
"""
function CommonSolve.init(prob::PBDMProblem, alg::DirectIteration = DirectIteration();
                          reproduction_fn = nothing, seed = nothing, kwargs...)
    _check_integrable(prob)
    sys = deepcopy(prob.populations)
    coupled = sys isa PopulationSystem
    parts = coupled ? _coupled_parts(prob.p) : nothing
    names = coupled ? _check_coupled_components(sys) : Symbol[]
    return PBDMIntegrator(prob, sys, alg, first(prob.tspan) - 1, full_state(sys), prob.p, reproduction_fn,
                          parts, names, SciMLBase.ReturnCode.Default)
end

function _step_day!(integ::PBDMIntegrator)
    day = integ.t + 1
    w = get_weather(integ.prob.weather, day)
    sys = integ.sys
    if sys isa PopulationSystem
        _coupled_day!(sys, w, day, integ.prob.approach, integ.parts, integ.comp_names)
        user_p = integ.parts.user_p
        for cb in integ.parts.cb_end_of_day
            cb.fn(sys, w, day, user_p)
        end
    else
        _single_day!(integ.prob.density, sys, w, integ.prob.approach, integ.p, day, integ.reproduction_fn)
    end
    setfield!(integ, :t, day)
    setfield!(integ, :u, full_state(sys))
    all(isfinite, integ.u) || setfield!(integ, :retcode, SciMLBase.ReturnCode.Failure)
    return integ
end

_running(integ::PBDMIntegrator) = integ.retcode != SciMLBase.ReturnCode.Failure &&
                                  integ.retcode != SciMLBase.ReturnCode.Terminated

function CommonSolve.step!(integ::PBDMIntegrator)
    _running(integ) && _step_day!(integ)
    return integ
end

function CommonSolve.step!(integ::PBDMIntegrator, dt::Real, stop_at_tdt::Bool = false)
    n = round(Int, dt)
    (abs(dt - n) <= 1e-9 && n >= 0) || throw(ArgumentError(
        "a PBDM steps whole days: dt must be a non-negative whole number, got $dt"))
    for _ in 1:n
        _running(integ) || break
        _step_day!(integ)
    end
    return integ
end

function SciMLBase.reinit!(integ::PBDMIntegrator, u0 = nothing; t0::Integer = first(integ.prob.tspan), kwargs...)
    sys = deepcopy(integ.prob.populations)
    u0 === nothing || set_full_state!(sys, u0)
    setfield!(integ, :sys, sys)
    setfield!(integ, :t, Int(t0) - 1)
    setfield!(integ, :u, full_state(sys))
    setfield!(integ, :retcode, SciMLBase.ReturnCode.Default)
    return integ
end

function SciMLBase.set_u!(integ::PBDMIntegrator, u::AbstractVector{<:Real})
    set_full_state!(integ.sys, u)
    setfield!(integ, :u, full_state(integ.sys))
    return integ
end

function SciMLBase.derivative_discontinuity!(integ::PBDMIntegrator, modified::Bool)
    modified && SciMLBase.set_u!(integ, copy(integ.u))
    return integ
end

function SciMLBase.terminate!(integ::PBDMIntegrator, retcode = SciMLBase.ReturnCode.Terminated)
    setfield!(integ, :retcode, retcode)
    return integ
end

SciMLBase.check_error(integ::PBDMIntegrator) =
    integ.retcode == SciMLBase.ReturnCode.Default ? SciMLBase.ReturnCode.Success : integ.retcode

# ============================================================================
# The annual map and its fixed point
# ============================================================================

"""
    annual_map(prob; period = 365, t0 = first(prob.tspan), reproduction_fn = nothing) -> F

The map `F(u)` that takes a state `u` (laid out as by [`full_state`](@ref)) at the start of
day `t0` to the state `period` days later, by stepping `prob`'s model through days
`t0, …, t0 + period - 1` with [`init`](@ref)`/step!`. Every evaluation reuses the same weather
days, so `F` is the year-to-year map under that weather year repeated: with weather that is
periodic in the day (`SinusoidalWeather`, or a `WeatherSeries` built from a periodic function)
its fixed point is the model's annual cycle. Scheduled events fire on their calendar days in
every evaluation, so a one-off event (a single release) becomes an annual one.
"""
function annual_map(prob::PBDMProblem; period::Integer = 365, t0::Integer = first(prob.tspan),
                    reproduction_fn = nothing)
    period > 0 || throw(ArgumentError("period must be positive"))
    base = remake(prob; tspan = (Int(t0), Int(t0) + Int(period) - 1))
    return function F(u)
        integ = CommonSolve.init(remake(base; u0 = u), DirectIteration(); reproduction_fn)
        CommonSolve.step!(integ, period, true)
        return copy(integ.u)
    end
end

"""
    annual_fixed_point(prob; period = 365, u0 = full_state(prob), spinup = 0,
                       method = :finite_differences, h = 1e-6, tol = 1e-9, maxiters = 50,
                       reproduction_fn = nothing)
        -> (u, multipliers, jacobian, converged, iterations)

The annual cycle of `prob` under its repeated weather year: a fixed point of
[`annual_map`](@ref), found by DynamicalResilience's `fixed_point` (Newton's method) after
`spinup` plain iterations of the map from `u0`. `jacobian` is the `DiscreteJacobian` of the
annual map there and `multipliers` its eigenvalues, the cycle's Floquet multipliers: the cycle
is stable when every multiplier has modulus below one, and the leading modulus is the factor by
which a small perturbation shrinks per year.

The Jacobian is taken by finite differences, since PBDM's populations store `Float64` and do
not carry dual numbers. Each column uses a central difference with step `h · max(|u_j|, 1)`,
except where that step would make the entry negative: there it uses a forward difference.
Populations cannot be negative and PBDM code clips them at zero, so at a cycle on which a
species is absent the map has a kink, and a central difference would average the true growth
of a rare invader with the clipped side (halving it, for instance). The forward difference
gives the multiplier of the absent species correctly: its growth factor per year when rare,
so the cycle resists invasion when it is below one.

Newton converges to the cycle near its starting point whether it is stable or not; a model that
runs away (unbounded growth) has no cycle to find.
"""
function annual_fixed_point(prob::PBDMProblem; period::Integer = 365, u0::AbstractVector = full_state(prob),
                            spinup::Integer = 0, method::Symbol = :finite_differences, h::Real = 1e-6,
                            tol::Real = 1e-9, maxiters::Integer = 50, reproduction_fn = nothing)
    F = annual_map(prob; period, reproduction_fn)
    u = Vector{Float64}(u0)
    for _ in 1:spinup
        u = F(u)
    end
    fp = DynamicalResilience.fixed_point(F, u; method, h, tol, maxiters)
    method === :finite_differences || return fp
    J = _annual_jacobian(F, fp.u, h)
    return (u = fp.u, multipliers = eigvals(J), jacobian = ResilienceCore.DiscreteJacobian(J, fp.u),
            converged = fp.converged, iterations = fp.iterations)
end

# Finite-difference Jacobian of `F` at `u`: central, or forward where a central step would make
# the entry negative (an absent species, whose population code clips at zero).
function _annual_jacobian(F, u::AbstractVector{Float64}, h::Real)
    n = length(u)
    J = Matrix{Float64}(undef, n, n)
    Fu = nothing
    for j in 1:n
        s = h * max(abs(u[j]), 1.0)
        e = zeros(n)
        e[j] = s
        if u[j] - s < 0
            Fu === nothing && (Fu = F(u))
            J[:, j] = (F(u .+ e) .- Fu) ./ s
        else
            J[:, j] = (F(u .+ e) .- F(u .- e)) ./ (2s)
        end
    end
    return J
end

"""
    linearise(prob::PBDMProblem; period = 365, check = true, kwargs...) -> DiscreteJacobian

The `DiscreteJacobian` of the annual map of `prob` at its annual cycle, from
[`annual_fixed_point`](@ref) (`kwargs` are passed to it). ResilienceCore's metrics then apply
per year: `resilience` is `−log ρ`, with `ρ` the leading Floquet multiplier's modulus. With
`check = true` an error is thrown if Newton's method did not converge.
"""
function ResilienceCore.linearise(prob::PBDMProblem; period::Integer = 365, check::Bool = true, kwargs...)
    fp = annual_fixed_point(prob; period, kwargs...)
    check && !fp.converged && throw(ArgumentError(
        "no annual cycle found: Newton's method did not converge in $(fp.iterations) steps; " *
        "try `spinup` years of plain iteration first, or check that the model neither runs away nor dies out"))
    return fp.jacobian
end

# ============================================================================
# Consumer–resource theory: equilibria as Jacobians, and the model as an ODEProblem
# ============================================================================

"""
    linearise(r::EquilibriumResult) -> ContinuousJacobian

The Jacobian of the consumer–resource model at the interior equilibrium found by
[`find_equilibrium`](@ref), with the equilibrium `(M₁*, M₂*)`, as ResilienceCore's
`ContinuousJacobian`. Its eigenvalues are `r.eigenvalues`. Throws if there is no interior
equilibrium (the consumer cannot persist).
"""
function ResilienceCore.linearise(r::EquilibriumResult)
    (isfinite(r.M1_star) && isfinite(r.M2_star) && all(isfinite, r.jacobian)) || throw(ArgumentError(
        "no interior equilibrium: the consumer cannot persist (metabolic compensation point ≥ 1, or M₁* ≤ 0)"))
    return ResilienceCore.ContinuousJacobian(r.jacobian, [r.M1_star, r.M2_star])
end

function _consumer_resource_rhs(u, p, t)
    M1, M2 = u[1], u[2]
    acquired = M2 > 0 ? 1 - exp(-p.a * M1 / (p.d2 * M2)) : one(M1)   # fraction of demand met
    R2 = respiration_rate(p.resp, p.T)
    return [p.r * M1 * (1 - M1 / p.K) - p.d2 * M2 * acquired,
            (p.eps * p.d2 * acquired - R2) * M2]
end

"""
    consumer_resource_problem(fr::FraserGilbertResponse, resp::Q10Respiration, T;
                              intrinsic_rate, carrying_capacity, consumer_demand_rate,
                              conversion_efficiency = 1.0, u0 = nothing, tspan = (0.0, 100.0))
        -> ODEProblem

The resource–consumer model behind [`find_equilibrium`](@ref) as a SciML `ODEProblem`:

    dM₁/dt = r M₁ (1 − M₁/K) − d₂ M₂ (1 − exp(−a M₁/(d₂ M₂)))
    dM₂/dt = (ε d₂ (1 − exp(−a M₁/(d₂ M₂))) − R₂(T)) M₂

with `R₂(T)` the consumer's `Q10Respiration` at temperature `T`. The parameters `p` are a
`NamedTuple` with fields `a, r, K, d2, eps, resp, T`, so DynamicalResilience's
`resilience_along(prob, (p, T) -> merge(p, (T = T,)), temperatures)` follows the equilibrium
along temperature. `u0` defaults to the equilibrium from `find_equilibrium` when it exists.
"""
function consumer_resource_problem(fr::FraserGilbertResponse, resp::Q10Respiration, T::Real;
                                   intrinsic_rate::Real, carrying_capacity::Real, consumer_demand_rate::Real,
                                   conversion_efficiency::Real = 1.0, u0 = nothing, tspan = (0.0, 100.0))
    p = (a = Float64(fr.a), r = Float64(intrinsic_rate), K = Float64(carrying_capacity),
         d2 = Float64(consumer_demand_rate), eps = Float64(conversion_efficiency), resp = resp, T = Float64(T))
    if u0 === nothing
        eq = find_equilibrium(fr, resp, T; intrinsic_rate, carrying_capacity, consumer_demand_rate,
                              conversion_efficiency)
        u0 = isfinite(eq.M1_star) ? [eq.M1_star, eq.M2_star] : [Float64(carrying_capacity), 1.0]
    end
    return SciMLBase.ODEProblem(_consumer_resource_rhs, Vector{Float64}(u0), tspan, p)
end
