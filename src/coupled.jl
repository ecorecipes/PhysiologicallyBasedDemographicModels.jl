"""
Coupled population systems — multi-type, multi-species, and metapopulation PBDMs.

Provides `PopulationSystem`, transition/interaction rules, scheduled events,
an extended `PBDMSolution` for coupled systems, and a `DirectIteration` solver
that orchestrates daily stepping, event application, and rule evaluation
without hand-rolled vignette loops.
"""

# ============================================================================
# Structure tags
# ============================================================================

"""
    MultiPopulationPBDM <: AbstractPBDMStructure

Solver-facing umbrella tag for any problem with multiple named `Population`
components.  Semantic subtypes exist for documentation and convenience
dispatch but all share the same solver path.
"""
abstract type MultiPopulationPBDM <: AbstractPBDMStructure end

"""
    MultiTypePBDM <: MultiPopulationPBDM

One biological species split into distinct types (sex, mating status,
genotype, sterile vs wild, etc.).  Example: screwworm SIT with wild males,
sterile males, mated females, unmated females, sterile-mated females.
"""
struct MultiTypePBDM <: MultiPopulationPBDM end

"""
    MultiSpeciesPBDMNew <: MultiPopulationPBDM

Multiple biological species interacting through trophic links,
competition, parasitism, or mutualism.  Replaces the old placeholder
`MultiSpeciesPBDM` with an actual solver implementation.
"""
struct MultiSpeciesPBDMNew <: MultiPopulationPBDM end

"""
    MetapopulationPBDM <: MultiPopulationPBDM

Spatially replicated populations connected by dispersal.
"""
struct MetapopulationPBDM <: MultiPopulationPBDM end

# ============================================================================
# StateVariable — mutable auxiliary state in the solver loop
# ============================================================================

"""
    AbstractStateVariable

Supertype for named mutable state variables that live inside a
`PopulationSystem` and are accessible to rules, events, and observables.

State variables can optionally have an `update_fn` that the solver calls
automatically each day (after events, before stepping).
"""
abstract type AbstractStateVariable end

"""
    ScalarState{T}(name, init; update=nothing)

A single mutable scalar value.

`update` signature: `(current_value, system, weather_day, day, p) -> new_value`.
If `nothing`, the variable is only modified manually by rules/events.
"""
struct ScalarState{T<:Real, F} <: AbstractStateVariable
    name::Symbol
    value::Ref{T}
    update_fn::F
end

function ScalarState(name::Symbol, init::T; update=nothing) where {T<:Real}
    ScalarState{T, typeof(update)}(name, Ref(init), update)
end

"""
    ArrayState{T}(name, init; update=nothing)

A mutable vector of values (e.g., per-stage carbon pools, per-patch water).

`update` signature: `(vec, system, weather_day, day, p) -> nothing` (mutate in-place).
"""
struct ArrayState{T<:Real, F} <: AbstractStateVariable
    name::Symbol
    value::Vector{T}
    update_fn::F
end

function ArrayState(name::Symbol, init::Vector{T}; update=nothing) where {T<:Real}
    ArrayState{T, typeof(update)}(name, copy(init), update)
end

"""
    DictState{K,V}(name, init; update=nothing)

A mutable dictionary (e.g., per-component or per-patch state maps).

`update` signature: `(dict, system, weather_day, day, p) -> nothing` (mutate in-place).
"""
struct DictState{K, V, F} <: AbstractStateVariable
    name::Symbol
    value::Dict{K, V}
    update_fn::F
end

function DictState(name::Symbol, init::Dict{K,V}; update=nothing) where {K, V}
    DictState{K, V, typeof(update)}(name, copy(init), update)
end

"""Get the current value of a state variable."""
get_state(sv::ScalarState) = sv.value[]
get_state(sv::ArrayState) = sv.value
get_state(sv::DictState) = sv.value

"""Set the value of a scalar state variable."""
set_state!(sv::ScalarState, v) = (sv.value[] = v)

"""Check whether a state variable has an automatic update function."""
has_auto_update(sv::AbstractStateVariable) = sv.update_fn !== nothing

"""
    update_state!(sv, system, weather_day, day, p)

Call the state variable's automatic update function (no-op if `update_fn` is nothing).
"""
function update_state!(sv::ScalarState, sys, w, day::Int, p)
    sv.update_fn === nothing && return
    sv.value[] = sv.update_fn(sv.value[], sys, w, day, p)
end

function update_state!(sv::ArrayState, sys, w, day::Int, p)
    sv.update_fn === nothing && return
    sv.update_fn(sv.value, sys, w, day, p)
end

function update_state!(sv::DictState, sys, w, day::Int, p)
    sv.update_fn === nothing && return
    sv.update_fn(sv.value, sys, w, day, p)
end

"""Snapshot the current value for history recording."""
snapshot(sv::ScalarState) = sv.value[]
snapshot(sv::ArrayState) = copy(sv.value)
snapshot(sv::DictState) = copy(sv.value)

function Base.show(io::IO, sv::ScalarState)
    print(io, "ScalarState(:$(sv.name), $(sv.value[]))")
end
function Base.show(io::IO, sv::ArrayState)
    print(io, "ArrayState(:$(sv.name), $(length(sv.value))-element)")
end
function Base.show(io::IO, sv::DictState)
    print(io, "DictState(:$(sv.name), $(length(sv.value)) entries)")
end

# ============================================================================
# BulkPopulation — unstructured scalar population
# ============================================================================

"""
    BulkPopulation{T}(name, init; growth_fn=nothing, K=Inf)

An unstructured (non-stage-structured) population represented by a single
scalar value. Participates in `PopulationSystem` alongside distributed-delay
`Population`s.

Useful for logistic/exponential growth models, soil pathogen inoculum,
phytoplankton density, or any quantity that doesn't need substage tracking.

# Fields
- `name::Symbol`: population identifier
- `value::Ref{T}`: current population size (mutable)
- `growth_fn`: `(N, weather_day, day, p) -> N_new`; `nothing` = rule-managed only
- `K::T`: carrying capacity (`Inf` = unlimited)
"""
struct BulkPopulation{T<:Real, F}
    name::Symbol
    value::Ref{T}
    growth_fn::F
    K::T
end

function BulkPopulation(name::Symbol, init::Real; growth_fn=nothing, K::Real=Inf)
    T = Float64
    BulkPopulation{T, typeof(growth_fn)}(name, Ref(T(init)), growth_fn, T(K))
end

total_population(bp::BulkPopulation) = max(bp.value[], 0.0)
n_stages(::BulkPopulation) = 1

"""Set the absolute value of a BulkPopulation."""
function set_value!(bp::BulkPopulation, val::Real)
    bp.value[] = clamp(val, 0.0, bp.K)
end

function inject!(bp::BulkPopulation, amount::Real)
    bp.value[] = min(bp.value[] + amount, bp.K)
end

function inject!(bp::BulkPopulation, ::Int, amount::Real)
    inject!(bp, amount)
end

function remove_fraction!(bp::BulkPopulation, frac::Real)
    0 <= frac <= 1 || throw(ArgumentError("fraction must be in [0, 1]"))
    bp.value[] *= (1.0 - frac)
end

function remove_fraction!(bp::BulkPopulation, ::Int, frac::Real)
    remove_fraction!(bp, frac)
end

"""
    step_bulk!(bp, weather_day, day, p)

Advance a `BulkPopulation` by one day using its `growth_fn`.
No-op if `growth_fn` is `nothing`.
"""
function step_bulk!(bp::BulkPopulation, w, day::Int, p)
    bp.growth_fn === nothing && return
    bp.value[] = clamp(bp.growth_fn(bp.value[], w, day, p), 0.0, bp.K)
end

function Base.show(io::IO, bp::BulkPopulation)
    print(io, "BulkPopulation(:$(bp.name), $(bp.value[]))")
end

# ============================================================================
# PopulationSystem — named collection of Population components
# ============================================================================

"""
    PopulationComponent{T}

A single named component in a `PopulationSystem`.

# Fields
- `population`: a `Population{T}` (stage-structured) or `BulkPopulation{T}` (scalar)
- `species::Symbol`: biological species label (for grouping)
- `type::Symbol`: sub-type within a species (`:default` if not applicable)
- `patch::Symbol`: spatial patch label (`:default` if not applicable)
"""
struct PopulationComponent{T<:Real}
    population::Union{Population{T}, BulkPopulation{T}}
    species::Symbol
    type::Symbol
    patch::Symbol
end

function PopulationComponent(pop::Population{T};
                             species::Symbol=pop.name,
                             type::Symbol=:default,
                             patch::Symbol=:default) where {T}
    PopulationComponent{T}(pop, species, type, patch)
end

function PopulationComponent(pop::BulkPopulation{T};
                             species::Symbol=pop.name,
                             type::Symbol=:default,
                             patch::Symbol=:default) where {T}
    PopulationComponent{T}(pop, species, type, patch)
end

"""
    PopulationSystem{T}

An ordered, named collection of `PopulationComponent`s with optional
`StateVariable`s for auxiliary mutable state.

Access components by name: `system[:wild_male]`.
Access state variables: `get_state(system, :allele_freq)`.
Iterate with `for (name, comp) in pairs(system)`.
"""
struct PopulationSystem{T<:Real}
    names::Vector{Symbol}
    components::Dict{Symbol, PopulationComponent{T}}
    order::Vector{Symbol}  # deterministic iteration order
    state::Dict{Symbol, AbstractStateVariable}
end

function PopulationSystem(pairs::Pair{Symbol}...; state::Vector=AbstractStateVariable[])
    T = Float64
    names = Symbol[]
    components = Dict{Symbol, PopulationComponent{T}}()
    for (name, val) in pairs
        name in names && throw(ArgumentError("duplicate component name: $name"))
        push!(names, name)
        if val isa Population
            components[name] = PopulationComponent(val)
        elseif val isa BulkPopulation
            components[name] = PopulationComponent(val)
        elseif val isa PopulationComponent
            components[name] = val
        else
            throw(ArgumentError("expected Population, BulkPopulation, or PopulationComponent, got $(typeof(val))"))
        end
    end
    length(names) > 0 || throw(ArgumentError("PopulationSystem must have ≥ 1 component"))
    state_dict = Dict{Symbol, AbstractStateVariable}(sv.name => sv for sv in state)
    PopulationSystem{T}(names, components, copy(names), state_dict)
end

Base.getindex(sys::PopulationSystem, name::Symbol) = sys.components[name]
Base.haskey(sys::PopulationSystem, name::Symbol) = haskey(sys.components, name)
Base.keys(sys::PopulationSystem) = sys.order
Base.length(sys::PopulationSystem) = length(sys.order)
Base.pairs(sys::PopulationSystem) = ((n, sys.components[n]) for n in sys.order)

function Base.show(io::IO, sys::PopulationSystem)
    ns = length(sys.state)
    state_str = ns > 0 ? ", $ns state vars" : ""
    print(io, "PopulationSystem($(length(sys)) components: $(join(sys.order, ", "))$state_str)")
end

"""Get a state variable's current value by name."""
get_state(sys::PopulationSystem, name::Symbol) = get_state(sys.state[name])

"""Set a scalar state variable's value by name."""
set_state!(sys::PopulationSystem, name::Symbol, v) = set_state!(sys.state[name], v)

"""Check if a state variable exists."""
has_state(sys::PopulationSystem, name::Symbol) = haskey(sys.state, name)

"""Total population summed across all components."""
function total_population(sys::PopulationSystem)
    return sum(total_population(c.population) for (_, c) in pairs(sys))
end

"""Total population for a single named component."""
function component_total(sys::PopulationSystem, name::Symbol)
    return total_population(sys[name].population)
end

"""Named totals as a Dict."""
function component_totals(sys::PopulationSystem{T}) where {T}
    return Dict{Symbol, T}(n => total_population(c.population) for (n, c) in pairs(sys))
end

"""Filter components by species."""
function by_species(sys::PopulationSystem, species::Symbol)
    return [(n, c) for (n, c) in pairs(sys) if c.species == species]
end

"""Filter components by type."""
function by_type(sys::PopulationSystem, type::Symbol)
    return [(n, c) for (n, c) in pairs(sys) if c.type == type]
end

"""Filter components by patch."""
function by_patch(sys::PopulationSystem, patch::Symbol)
    return [(n, c) for (n, c) in pairs(sys) if c.patch == patch]
end

"""Inject inflow into a specific stage of a named component."""
function inject!(sys::PopulationSystem, name::Symbol, stage_idx::Int, amount::Real)
    pop = sys[name].population
    if pop isa BulkPopulation
        inject!(pop, amount)
    else
        1 <= stage_idx <= n_stages(pop) ||
            throw(BoundsError("stage $stage_idx out of range for $name ($(n_stages(pop)) stages)"))
        pop.stages[stage_idx].delay.W[1] += amount
    end
    return nothing
end

"""Inject inflow into the first stage of a named component."""
function inject!(sys::PopulationSystem, name::Symbol, amount::Real)
    inject!(sys, name, 1, amount)
end

"""Remove a fraction of all substages of a named component."""
function remove_fraction!(sys::PopulationSystem, name::Symbol, fraction::Real)
    0 <= fraction <= 1 || throw(ArgumentError("fraction must be in [0, 1]"))
    pop = sys[name].population
    if pop isa BulkPopulation
        remove_fraction!(pop, fraction)
    else
        for stage in pop.stages
            stage.delay.W .*= (1.0 - fraction)
        end
    end
    return nothing
end

"""Remove a fraction from a specific stage of a named component."""
function remove_fraction!(sys::PopulationSystem, name::Symbol, stage_idx::Int, fraction::Real)
    0 <= fraction <= 1 || throw(ArgumentError("fraction must be in [0, 1]"))
    pop = sys[name].population
    if pop isa BulkPopulation
        remove_fraction!(pop, fraction)
    else
        pop.stages[stage_idx].delay.W .*= (1.0 - fraction)
    end
    return nothing
end

# ============================================================================
# Abstract rule and event interfaces
# ============================================================================

"""
    AbstractInteractionRule

Supertype for rules evaluated each day after all components are stepped.
Rules can transfer mass between components, apply mortality, compute
reproduction, etc.

Implement `apply_rule!(rule, system, weather_day, day, p) -> NamedTuple`
where the returned tuple contains any per-day observables the rule wants
to record.
"""
abstract type AbstractInteractionRule end

"""
    AbstractScheduledEvent

Supertype for discrete events applied on specific days (releases,
introductions, sprays, immigration pulses).

Implement `apply_event!(event, system, weather_day, day, p) -> nothing`.
The event should return `true` if it fired on this day, `false` otherwise.
"""
abstract type AbstractScheduledEvent end

# ============================================================================
# Concrete events
# ============================================================================

"""
    PulseRelease(target, stage_idx, amount, interval; start_day=1, end_day=typemax(Int))

Release `amount` individuals into `target` component (stage `stage_idx`)
every `interval` days, starting at `start_day`.
"""
struct PulseRelease <: AbstractScheduledEvent
    target::Symbol
    stage_idx::Int
    amount::Float64
    interval::Int
    start_day::Int
    end_day::Int
end

function PulseRelease(target::Symbol, stage_idx::Int, amount::Real, interval::Int;
                      start_day::Int=1, end_day::Int=typemax(Int))
    PulseRelease(target, stage_idx, Float64(amount), interval, start_day, end_day)
end

function PulseRelease(target::Symbol, amount::Real, interval::Int; kwargs...)
    PulseRelease(target, 1, amount, interval; kwargs...)
end

function apply_event!(ev::PulseRelease, sys::PopulationSystem, ::Any, day::Int, ::Any)
    day < ev.start_day && return false
    day > ev.end_day && return false
    (day - ev.start_day) % ev.interval == 0 || return false
    inject!(sys, ev.target, ev.stage_idx, ev.amount)
    return true
end

"""
    SingleDayRelease(target, stage_idx, amount, day)

Release `amount` individuals into `target` on exactly one day.
"""
struct SingleDayRelease <: AbstractScheduledEvent
    target::Symbol
    stage_idx::Int
    amount::Float64
    day::Int
end

function SingleDayRelease(target::Symbol, amount::Real, day::Int; stage_idx::Int=1)
    SingleDayRelease(target, stage_idx, Float64(amount), day)
end

function apply_event!(ev::SingleDayRelease, sys::PopulationSystem, ::Any, day::Int, ::Any)
    day == ev.day || return false
    inject!(sys, ev.target, ev.stage_idx, ev.amount)
    return true
end

"""
    SprayEvent(targets, kill_fractions, days)

Apply kill fractions to named components on specified days.
`targets` and `kill_fractions` should be parallel vectors.
`days` is the set of calendar days on which the spray fires.
"""
struct SprayEvent <: AbstractScheduledEvent
    targets::Vector{Symbol}
    kill_fractions::Vector{Float64}
    days::Vector{Int}
end

function apply_event!(ev::SprayEvent, sys::PopulationSystem, ::Any, day::Int, ::Any)
    day in ev.days || return false
    for (t, kf) in zip(ev.targets, ev.kill_fractions)
        haskey(sys, t) && remove_fraction!(sys, t, kf)
    end
    return true
end

"""
    ConditionalRelease(target, stage_idx, amount_fn, interval; start_day, end_day)

Like `PulseRelease` but the release amount is computed dynamically by
`amount_fn(system, weather_day, day, p) -> Float64`.
"""
struct ConditionalRelease{F} <: AbstractScheduledEvent
    target::Symbol
    stage_idx::Int
    amount_fn::F
    interval::Int
    start_day::Int
    end_day::Int
end

function ConditionalRelease(target::Symbol, amount_fn, interval::Int;
                            stage_idx::Int=1, start_day::Int=1,
                            end_day::Int=typemax(Int))
    ConditionalRelease(target, stage_idx, amount_fn, interval, start_day, end_day)
end

function apply_event!(ev::ConditionalRelease, sys::PopulationSystem, w, day::Int, p)
    day < ev.start_day && return false
    day > ev.end_day && return false
    (day - ev.start_day) % ev.interval == 0 || return false
    amount = ev.amount_fn(sys, w, day, p)
    amount > 0 && inject!(sys, ev.target, ev.stage_idx, amount)
    return true
end

"""
    WeatherConditionalEvent(name, predicate, action)

An event that fires when a weather/state predicate is true.
Unlike schedule-based events, this checks conditions every day.

- `predicate(weather_day, day, system, p) -> Bool`
- `action(system, weather_day, day, p) -> Bool` (true if it acted)

# Examples
```julia
# Diapause entry when photoperiod drops below threshold
WeatherConditionalEvent(:diapause_entry,
    (w, d, sys, p) -> w.photoperiod < 12.5,
    (sys, w, d, p) -> begin
        # transfer adults to diapause pool
        ...
        true
    end
)
```
"""
struct WeatherConditionalEvent{F, G} <: AbstractScheduledEvent
    name::Symbol
    predicate::F
    action::G
end

function apply_event!(ev::WeatherConditionalEvent, sys::PopulationSystem, w, day::Int, p)
    ev.predicate(w, day, sys, p) || return false
    return ev.action(sys, w, day, p)
end

# ============================================================================
# Concrete interaction rules
# ============================================================================

"""
    TransferRule(source, source_stage, target, target_stage, fraction_fn)

Each day, transfer a fraction of source → target.
`fraction_fn(system, weather_day, day, p) -> Float64` in [0, 1].
"""
struct TransferRule{F} <: AbstractInteractionRule
    source::Symbol
    source_stage::Int
    target::Symbol
    target_stage::Int
    fraction_fn::F
end

function TransferRule(source::Symbol, target::Symbol, fraction_fn;
                      source_stage::Int=1, target_stage::Int=1)
    TransferRule(source, source_stage, target, target_stage, fraction_fn)
end

function apply_rule!(rule::TransferRule, sys::PopulationSystem, w, day::Int, p)
    frac = rule.fraction_fn(sys, w, day, p)
    frac = clamp(frac, 0.0, 1.0)
    source_pop = sys[rule.source].population
    amount = delay_total(source_pop.stages[rule.source_stage].delay) * frac
    remove_fraction!(sys, rule.source, rule.source_stage, frac)
    inject!(sys, rule.target, rule.target_stage, amount)
    return (transferred=amount,)
end

"""
    ReproductionRule(parents, offspring_target, offspring_stage, reproduction_fn)

Each day, compute offspring from the state of parent components and inject
them into the offspring target.

`reproduction_fn(system, weather_day, day, p) -> Float64` returns the
number of new individuals to inject.
"""
struct ReproductionRule{F} <: AbstractInteractionRule
    offspring_target::Symbol
    offspring_stage::Int
    reproduction_fn::F
end

function ReproductionRule(target::Symbol, reproduction_fn; stage::Int=1)
    ReproductionRule(target, stage, reproduction_fn)
end

function apply_rule!(rule::ReproductionRule, sys::PopulationSystem, w, day::Int, p)
    offspring = max(0.0, rule.reproduction_fn(sys, w, day, p))
    inject!(sys, rule.offspring_target, rule.offspring_stage, offspring)
    return (offspring=offspring,)
end

"""
    MortalityRule(target, mortality_fn)

Apply additional daily mortality to a named component.
`mortality_fn(system, weather_day, day, p) -> Float64` in [0, 1].
"""
struct MortalityRule{F} <: AbstractInteractionRule
    target::Symbol
    mortality_fn::F
end

function apply_rule!(rule::MortalityRule, sys::PopulationSystem, w, day::Int, p)
    frac = clamp(rule.mortality_fn(sys, w, day, p), 0.0, 1.0)
    remove_fraction!(sys, rule.target, frac)
    return (mortality=frac,)
end

"""
    PredationRule(predator, prey, response, conversion; predator_stage, prey_stage)

Holling-style predation/parasitism: predator adults attack prey,
removing prey and optionally injecting converted biomass back into
the predator.
"""
struct PredationRule{FR<:AbstractFunctionalResponse} <: AbstractInteractionRule
    predator::Symbol
    predator_stage::Int
    prey::Symbol
    prey_stage::Int
    response::FR
    conversion::Float64
end

function PredationRule(predator::Symbol, prey::Symbol, response::AbstractFunctionalResponse;
                       conversion::Float64=0.0,
                       predator_stage::Int=-1,  # -1 means last (adult)
                       prey_stage::Int=0)       # 0 means all stages
    PredationRule(predator, predator_stage, prey, prey_stage, response, conversion)
end

function apply_rule!(rule::PredationRule, sys::PopulationSystem, w, day::Int, p)
    pred_pop = sys[rule.predator].population
    prey_pop = sys[rule.prey].population

    pred_stage_idx = rule.predator_stage == -1 ? n_stages(pred_pop) : rule.predator_stage
    n_predators = delay_total(pred_pop.stages[pred_stage_idx].delay)

    if rule.prey_stage == 0
        n_prey = total_population(prey_pop)
    else
        n_prey = delay_total(prey_pop.stages[rule.prey_stage].delay)
    end

    n_predators <= 0 && return (consumed=0.0,)
    n_prey <= 0 && return (consumed=0.0,)

    per_capita = functional_response(rule.response, n_prey)
    consumed = min(per_capita * n_predators, n_prey * 0.9)

    # Remove consumed prey
    if rule.prey_stage == 0
        frac = consumed / max(n_prey, 1e-10)
        remove_fraction!(sys, rule.prey, clamp(frac, 0.0, 1.0))
    else
        frac = consumed / max(n_prey, 1e-10)
        remove_fraction!(sys, rule.prey, rule.prey_stage, clamp(frac, 0.0, 1.0))
    end

    # Convert consumed to predator biomass
    if rule.conversion > 0 && consumed > 0
        inject!(sys, rule.predator, 1, consumed * rule.conversion)
    end

    return (consumed=consumed,)
end

"""
    CompetitionRule(competitors, resource_fn, allocation_fn)

Generic competition: `resource_fn` computes shared resource,
`allocation_fn` partitions it among competitors and applies stress.
"""
struct CompetitionRule{RF, AF} <: AbstractInteractionRule
    competitors::Vector{Symbol}
    resource_fn::RF
    allocation_fn::AF
end

function apply_rule!(rule::CompetitionRule, sys::PopulationSystem, w, day::Int, p)
    resource = rule.resource_fn(sys, w, day, p)
    stresses = rule.allocation_fn(sys, rule.competitors, resource, w, day, p)
    return (resource=resource, stresses=stresses)
end

"""
    CustomRule(name, apply_fn)

Escape hatch: `apply_fn(system, weather_day, day, p) -> NamedTuple`
can do anything.
"""
struct CustomRule{F} <: AbstractInteractionRule
    name::Symbol
    apply_fn::F
end

function apply_rule!(rule::CustomRule, sys::PopulationSystem, w, day::Int, p)
    return rule.apply_fn(sys, w, day, p)
end

# ============================================================================
# Pre-step stress rules
# ============================================================================

"""
    AbstractStressRule

Rules evaluated *before* stepping to compute per-component `stage_stress`
vectors that are passed to `step_population!`.
"""
abstract type AbstractStressRule end

"""
    StressRule(name, stress_fn)

Compute stress vectors for components before stepping.
`stress_fn(system, weather_day, day, p) -> Dict{Symbol, Vector{Float64}}`
mapping component names to per-stage stress vectors in [0, 1].

Components not present in the returned Dict receive zero stress.
"""
struct StressRule{F} <: AbstractStressRule
    name::Symbol
    stress_fn::F
end

function compute_stress(rule::StressRule, sys::PopulationSystem, w, day::Int, p)
    return rule.stress_fn(sys, w, day, p)
end

# ============================================================================
# Observable specification
# ============================================================================

"""
    Observable(name, fn)

Records a scalar value each day.
`fn(system, weather_day, day, p) -> Real`.
"""
struct Observable{F}
    name::Symbol
    fn::F
end

# ============================================================================
# PhaseCallback — hook into specific solver phases
# ============================================================================

"""
    CallbackPhase

Enum specifying when a `PhaseCallback` fires during the solver loop.

- `PRE_EVENT`:  before scheduled events
- `POST_EVENT`: after events, before state variable auto-update
- `PRE_STEP`:   after state update + stress computation, before stepping
- `POST_STEP`:  after stepping, before interaction rules
- `END_OF_DAY`: after observables are recorded
"""
@enum CallbackPhase begin
    PRE_EVENT
    POST_EVENT
    PRE_STEP
    POST_STEP
    END_OF_DAY
end

"""
    PhaseCallback(name, phase, fn)

A lightweight hook that runs at a specific solver phase.
`fn(system, weather_day, day, p) -> nothing` (can mutate system/state).

More general than `CustomRule` (POST_STEP only) or `StressRule` (PRE_STEP only)
for patterns that need to run at unconventional times.
"""
struct PhaseCallback{F}
    name::Symbol
    phase::CallbackPhase
    fn::F
end

# ============================================================================
# CoupledPBDMSolution
# ============================================================================

"""
    CoupledPBDMSolution{T}

Result of solving a coupled (multi-population) PBDMProblem.

# Fields
- `t`: calendar day labels
- `component_names`: ordered names of population components
- `component_totals`: Dict of name => Vector{T} daily totals
- `component_stage_totals`: Dict of name => Matrix{T} (stages × days)
- `state_history`: Dict of state_name => Vector (daily snapshots)
- `observables`: Dict of name => Vector{T} custom observable trajectories
- `event_log`: Vector of (day, event_type) tuples
- `rule_log`: Dict of rule_name => Vector{NamedTuple} per-day outputs
- `retcode`: :Success or :Failure
"""
struct CoupledPBDMSolution{T<:Real}
    t::Vector{Int}
    component_names::Vector{Symbol}
    component_totals::Dict{Symbol, Vector{T}}
    component_stage_totals::Dict{Symbol, Matrix{T}}
    state_history::Dict{Symbol, Vector}
    observables::Dict{Symbol, Vector{T}}
    event_log::Vector{Tuple{Int, Symbol}}
    rule_log::Dict{Symbol, Vector}
    retcode::Symbol
end

function Base.show(io::IO, sol::CoupledPBDMSolution)
    nt = length(sol.t)
    nc = length(sol.component_names)
    ns = length(sol.state_history)
    no = length(sol.observables)
    state_str = ns > 0 ? ", $ns state vars" : ""
    print(io, "CoupledPBDMSolution($nt days, $nc components$state_str, $no observables, retcode=$(sol.retcode))")
end

"""Get daily total trajectory for a component."""
Base.getindex(sol::CoupledPBDMSolution, name::Symbol) = sol.component_totals[name]

# ============================================================================
# GenomeState — allele frequency tracking for resistance genetics
# ============================================================================

"""
    GenomeState(name, locus; update=nothing)

Track allele frequencies at a diallelic locus within the coupled solver.
Wraps a `DialleleicLocus` and provides automatic Hardy-Weinberg genotype
computation.  The optional `update` function has signature
`(locus, system, weather_day, day, p) -> nothing` and is called in Phase 2.

# Examples
```julia
gs = GenomeState(:resistance, DialleleicLocus(0.01, 0.5))
fRR, fRS, fSS = get_genotypes(gs)
```
"""
struct GenomeState{F} <: AbstractStateVariable
    name::Symbol
    locus::DialleleicLocus
    update_fn::F
end

function GenomeState(name::Symbol, locus::DialleleicLocus; update=nothing)
    GenomeState{typeof(update)}(name, locus, update)
end

get_state(gs::GenomeState) = gs.locus.R
set_state!(gs::GenomeState, v) = (gs.locus.R = clamp(v, 0.0, 1.0))
has_auto_update(gs::GenomeState) = gs.update_fn !== nothing
snapshot(gs::GenomeState) = gs.locus.R

function update_state!(gs::GenomeState, sys, w, day, p)
    if gs.update_fn !== nothing
        gs.update_fn(gs.locus, sys, w, day, p)
        gs.locus.R = clamp(gs.locus.R, 0.0, 1.0)
    end
    nothing
end

"""Return (fRR, fRS, fSS) genotype frequencies from Hardy-Weinberg."""
get_genotypes(gs::GenomeState) = genotype_frequencies(gs.locus)

"""Return the `DialleleicLocus` object."""
get_locus(gs::GenomeState) = gs.locus

"""
    SelectionRule(genome_state, fitness_fn; name)

Apply one generation of single-locus Hardy-Weinberg selection to the
`GenomeState` named `genome_state`. `fitness_fn(sys, w, day, p) -> GenotypeFitness`
is called each day to compute the current genotype-specific relative fitness.

The rule updates the `GenomeState`'s allele frequency via `selection_step!`
and returns `(R_before=..., R_after=..., w_bar=...)` metrics. Use this rule
to couple daily spray/stress pressure to resistance evolution without writing
bespoke selection loops in the vignette.
"""
struct SelectionRule{F} <: AbstractInteractionRule
    name::Symbol
    genome_state::Symbol
    fitness_fn::F
end

function SelectionRule(genome_state::Symbol, fitness_fn;
        name::Symbol = Symbol("selection_", genome_state))
    return SelectionRule(name, genome_state, fitness_fn)
end

function apply_rule!(rule::SelectionRule, sys::PopulationSystem, w, day::Int, p)
    has_state(sys, rule.genome_state) || throw(ArgumentError(
        "SelectionRule: no GenomeState named :$(rule.genome_state) in system"))
    gs = sys.state[rule.genome_state]
    gs isa GenomeState || throw(ArgumentError(
        "SelectionRule: state :$(rule.genome_state) is not a GenomeState"))
    fitness = rule.fitness_fn(sys, w, day, p)
    fitness isa GenotypeFitness || throw(ArgumentError(
        "SelectionRule fitness_fn must return a GenotypeFitness, got $(typeof(fitness))"))

    R_before = gs.locus.R
    freq = genotype_frequencies(gs.locus)
    w_bar = freq.SS * fitness.w_SS + freq.SR * fitness.w_SR + freq.RR * fitness.w_RR
    selection_step!(gs.locus, fitness)
    return (R_before = R_before, R_after = gs.locus.R, w_bar = w_bar)
end

# ============================================================================
# DiapauseState + DiapauseRule — photoperiod-driven dormancy
# ============================================================================

"""
    DiapauseState(name, pool; cold_dd=0.0, induction_fn, emergence_fn, cold_survival_fn)

Track a dormant/diapausing population pool.  `induction_fn(daylength, T) -> Float64`
gives the fraction of active individuals entering diapause on a given day.
`emergence_fn(cum_dd) -> Float64` gives the fraction of the pool emerging.
`cold_survival_fn(cold_dd) -> Float64` gives survivorship as a function of
accumulated cold degree-days.
"""
struct DiapauseState{F1, F2, F3} <: AbstractStateVariable
    name::Symbol
    pool::Base.RefValue{Float64}
    cold_dd::Base.RefValue{Float64}
    induction_fn::F1
    emergence_fn::F2
    cold_survival_fn::F3
end

function DiapauseState(name::Symbol, pool::Real=0.0;
                       cold_dd::Real=0.0,
                       induction_fn=(dl, T) -> 0.0,
                       emergence_fn=(dd) -> 0.0,
                       cold_survival_fn=(cdd) -> 1.0)
    DiapauseState(name, Ref(Float64(pool)), Ref(Float64(cold_dd)),
                  induction_fn, emergence_fn, cold_survival_fn)
end

get_state(ds::DiapauseState) = ds.pool[]
set_state!(ds::DiapauseState, v) = (ds.pool[] = v)
has_auto_update(::DiapauseState) = false
snapshot(ds::DiapauseState) = (pool=ds.pool[], cold_dd=ds.cold_dd[])

function update_state!(::DiapauseState, sys, w, day, p)
    nothing  # Managed by DiapauseRule
end

"""
    DiapauseRule(source, diapause_state; T_cold_base=10.0)

Manage diapause dynamics: induction (active→dormant), cold degree-day
accumulation, cold mortality, and spring emergence (dormant→active).

Fires in the interaction-rule phase (Phase 4).  The `source` is the name
of the active `BulkPopulation`, and `diapause_state` names a `DiapauseState`
in the system.
"""
struct DiapauseRule <: AbstractInteractionRule
    source::Symbol
    diapause_state::Symbol
    T_cold_base::Float64
end

function DiapauseRule(source::Symbol, diapause_state::Symbol; T_cold_base::Float64=10.0)
    DiapauseRule(source, diapause_state, T_cold_base)
end

function apply_rule!(rule::DiapauseRule, sys::PopulationSystem, w, day::Int, p)
    ds = sys.state[rule.diapause_state]::DiapauseState
    T = w.T_mean
    dl = hasfield(typeof(w), :photoperiod) ? w.photoperiod : 12.0

    # Accumulate cold degree-days (below base)
    cold_dd_today = max(0.0, rule.T_cold_base - T)
    ds.cold_dd[] += cold_dd_today

    # Induction: fraction of active population entering diapause
    induction_frac = clamp(ds.induction_fn(dl, T), 0.0, 1.0)
    entering = 0.0
    if induction_frac > 0
        active_N = total_population(sys[rule.source].population)
        entering = active_N * induction_frac
        remove_fraction!(sys, rule.source, induction_frac)
        ds.pool[] += entering
    end

    # Cold mortality on diapausing pool
    lx = clamp(ds.cold_survival_fn(ds.cold_dd[]), 0.0, 1.0)
    ds.pool[] *= lx

    # Emergence: fraction of pool becoming active
    dd_today = max(0.0, T - rule.T_cold_base)
    emergence_frac = clamp(ds.emergence_fn(dd_today > 0 ? ds.cold_dd[] : 0.0), 0.0, 1.0)
    emerged = ds.pool[] * emergence_frac
    ds.pool[] -= emerged
    inject!(sys, rule.source, emerged)

    return (pool=ds.pool[], cold_dd=ds.cold_dd[], induction=induction_frac,
            entering=entering, emerged=emerged, survival=lx)
end

# ============================================================================
# PhenologyState — crop/host developmental milestones
# ============================================================================

"""
    PhenologyState(name, milestones; base_temp=10.0, init_dd=0.0)

Track cumulative degree-days and phenological phase based on milestone
thresholds.  `milestones` is a vector of `(name, dd_threshold)` pairs
in ascending order.

# Examples
```julia
ps = PhenologyState(:cotton, [(:squaring, 450.0), (:flowering, 700.0),
                               (:boll_fill, 1000.0), (:open_boll, 1400.0)])
```
"""
struct PhenologyState{F} <: AbstractStateVariable
    name::Symbol
    cum_dd::Base.RefValue{Float64}
    milestones::Vector{Tuple{Symbol, Float64}}
    current_phase::Base.RefValue{Symbol}
    base_temp::Float64
    dd_fn::F
end

function PhenologyState(name::Symbol, milestones::Vector;
                        base_temp::Float64=10.0, init_dd::Float64=0.0,
                        dd_fn=nothing)
    ms = [(Symbol(n), Float64(t)) for (n, t) in milestones]
    sort!(ms, by=last)
    initial_phase = :pre
    for (mname, mdd) in ms
        if init_dd >= mdd
            initial_phase = mname
        end
    end
    PhenologyState{typeof(dd_fn)}(name, Ref(init_dd), ms, Ref(initial_phase),
                                  base_temp, dd_fn)
end

get_state(ps::PhenologyState) = ps.cum_dd[]
set_state!(ps::PhenologyState, v) = (ps.cum_dd[] = v)
has_auto_update(::PhenologyState) = true
snapshot(ps::PhenologyState) = (cum_dd=ps.cum_dd[], phase=ps.current_phase[])

function update_state!(ps::PhenologyState, sys, w, day, p)
    # `dd_fn(w, day, p) -> Float64` overrides the default linear rule if
    # provided — useful for nonlinear development-rate functions (Brière,
    # Logan, Sharpe–DeMichele, etc.).
    dd = if ps.dd_fn === nothing
        max(0.0, w.T_mean - ps.base_temp)
    else
        max(0.0, Float64(ps.dd_fn(w, day, p)))
    end
    ps.cum_dd[] += dd
    # Advance phase
    for (mname, mdd) in ps.milestones
        if ps.cum_dd[] >= mdd
            ps.current_phase[] = mname
        end
    end
    nothing
end

"""Get the current phenological phase symbol."""
get_phase(ps::PhenologyState) = ps.current_phase[]
get_phase(sys::PopulationSystem, name::Symbol) = get_phase(sys.state[name]::PhenologyState)

"""Check whether a milestone has been reached."""
past_milestone(ps::PhenologyState, milestone::Symbol) =
    any(m -> m[1] == milestone && ps.cum_dd[] >= m[2], ps.milestones)

"""Get availability index (0–1) for a named resource (linear ramp from milestone DD)."""
function resource_availability(ps::PhenologyState, milestone::Symbol; ramp_dd::Float64=100.0)
    for (mname, mdd) in ps.milestones
        if mname == milestone
            return clamp((ps.cum_dd[] - mdd) / ramp_dd, 0.0, 1.0)
        end
    end
    return 0.0
end

# ============================================================================
# SoilState — soil-borne pathogen inoculum
# ============================================================================

"""
    SoilState(name; inoculum=0.0, virulence=0.0, update=nothing)

Track soil-borne pathogen inoculum density and virulence for disease models.
The `update` function signature is `(inoculum, virulence, sys, w, day, p) -> (new_inoc, new_vir)`.
Fires in Phase 2 if provided.
"""
struct SoilState{F} <: AbstractStateVariable
    name::Symbol
    inoculum::Base.RefValue{Float64}
    virulence::Base.RefValue{Float64}
    update_fn::F
end

function SoilState(name::Symbol; inoculum::Real=0.0, virulence::Real=0.0, update=nothing)
    SoilState{typeof(update)}(name, Ref(Float64(inoculum)), Ref(Float64(virulence)), update)
end

get_state(ss::SoilState) = (inoculum=ss.inoculum[], virulence=ss.virulence[])
set_state!(ss::SoilState, v::Tuple) = (ss.inoculum[] = v[1]; ss.virulence[] = v[2])
set_state!(ss::SoilState, v::NamedTuple) = (ss.inoculum[] = v.inoculum; ss.virulence[] = v.virulence)
has_auto_update(ss::SoilState) = ss.update_fn !== nothing
snapshot(ss::SoilState) = (inoculum=ss.inoculum[], virulence=ss.virulence[])

function update_state!(ss::SoilState, sys, w, day, p)
    if ss.update_fn !== nothing
        new_inoc, new_vir = ss.update_fn(ss.inoculum[], ss.virulence[], sys, w, day, p)
        ss.inoculum[] = new_inoc
        ss.virulence[] = new_vir
    end
    nothing
end

"""Get inoculum density from a SoilState."""
get_inoculum(ss::SoilState) = ss.inoculum[]
get_inoculum(sys::PopulationSystem, name::Symbol) = (sys.state[name]::SoilState).inoculum[]

"""Get virulence from a SoilState."""
get_virulence(ss::SoilState) = ss.virulence[]
get_virulence(sys::PopulationSystem, name::Symbol) = (sys.state[name]::SoilState).virulence[]

# ============================================================================
# SpatialGrid + DispersalRule — metapopulation dynamics
# ============================================================================

"""
    SpatialGrid(patches; connectivity=nothing)

A collection of `PopulationSystem` patches connected by dispersal.
`connectivity` is an optional NxN matrix giving dispersal weights between
patches (defaults to uniform).

Each patch is a full `PopulationSystem` that can be stepped independently.
`DispersalRule` coordinates inter-patch movement.
"""
struct SpatialGrid{T}
    patches::Vector{PopulationSystem{T}}
    patch_names::Vector{Symbol}
    connectivity::Matrix{Float64}
end

function SpatialGrid(patches::Vector{<:Pair}; connectivity=nothing)
    names = [first(p) for p in patches]
    systems = [last(p) for p in patches]
    n = length(patches)
    conn = connectivity === nothing ? ones(n, n) - I : connectivity
    T = Float64
    SpatialGrid{T}(systems, names, conn)
end

function SpatialGrid(patches::Vector{PopulationSystem{T}}; connectivity=nothing) where T
    names = [Symbol("patch_", i) for i in eachindex(patches)]
    n = length(patches)
    conn = connectivity === nothing ? ones(n, n) - I : connectivity
    SpatialGrid{T}(patches, names, conn)
end

Base.length(sg::SpatialGrid) = length(sg.patches)
Base.getindex(sg::SpatialGrid, i::Int) = sg.patches[i]
Base.getindex(sg::SpatialGrid, name::Symbol) = sg.patches[findfirst(==(name), sg.patch_names)]

"""
    DispersalRule(species; emigration_fn, patch_finding=0.5)

Rule that moves individuals of `species` between patches in a `SpatialGrid`.
`emigration_fn(N, supply_demand, w, day, p) -> fraction_emigrating` determines
how many leave each patch.  `patch_finding` is the fraction of emigrants that
successfully find a new patch.  Distribution among patches is proportional to
connectivity weights.
"""
struct DispersalRule{F} <: AbstractInteractionRule
    species::Symbol
    emigration_fn::F
    patch_finding::Float64
end

function DispersalRule(species::Symbol; emigration_fn=(N, phi, w, day, p) -> 0.0,
                       patch_finding::Float64=0.5)
    DispersalRule{typeof(emigration_fn)}(species, emigration_fn, patch_finding)
end

"""
    apply_dispersal!(rule, grid, w, day, p)

Apply dispersal across all patches in a SpatialGrid.
"""
function apply_dispersal!(rule::DispersalRule, grid::SpatialGrid, w, day, p)
    n = length(grid)
    species = rule.species
    emigrants = zeros(n)
    remaining = zeros(n)

    # Compute emigration from each patch
    for i in 1:n
        patch = grid.patches[i]
        haskey(patch, species) || continue
        N = total_population(patch[species].population)
        phi = 1.0  # default supply/demand ratio
        frac = clamp(rule.emigration_fn(N, phi, w, day, p), 0.0, 1.0)
        emigrants[i] = N * frac
        remove_fraction!(patch, species, frac)
        remaining[i] = total_population(patch[species].population)
    end

    # Distribute survivors among patches
    total_em = sum(emigrants)
    survivors = total_em * rule.patch_finding
    if survivors > 0 && n > 1
        for j in 1:n
            weight_sum = sum(grid.connectivity[i, j] for i in 1:n if i != j)
            weight_sum <= 0 && continue
            incoming = sum(emigrants[i] * grid.connectivity[i, j] / weight_sum
                          for i in 1:n if i != j) * rule.patch_finding
            inject!(grid.patches[j], species, incoming)
        end
    end

    return (total_emigrants=total_em, survivors=survivors)
end

# ============================================================================
# EnsemblePBDMProblem — SciML-inspired parameter sweeps
# ============================================================================

"""
    EnsemblePBDMProblem(prob; prob_func, output_func, reduction, u_init)

SciML-style ensemble problem for running multiple PBDM trajectories with
different parameters, initial conditions, or weather.

*   `prob_func(prob, i, repeat) -> new_prob`: Modify the template problem for trajectory `i`.
*   `output_func(sol, i) -> (output, rerun)`: Extract output from each solution.
*   `reduction(u, data, I) -> (u, converged)`: Reduce batch results.
*   `u_init`: Initial accumulator for reduction.

# Example
```julia
prob = PBDMProblem(system, weather, (1, 365); p=(r=0.1,))
ens = EnsemblePBDMProblem(prob;
    prob_func = (prob, i, _) -> remake(prob; p=(r=0.05*i,)),
    output_func = (sol, i) -> (sol[:pest][end], false))
results = solve(ens, DirectIteration(); trajectories=10)
```
"""
struct EnsemblePBDMProblem{P, PF, OF, RF, U}
    prob::P
    prob_func::PF
    output_func::OF
    reduction::RF
    u_init::U
end

function EnsemblePBDMProblem(prob;
    prob_func = (prob, i, repeat) -> prob,
    output_func = (sol, i) -> (sol, false),
    reduction = (u, data, I) -> (append!(u, data), false),
    u_init = [])
    EnsemblePBDMProblem(prob, prob_func, output_func, reduction, u_init)
end

"""
    EnsemblePBDMSolution(solutions, converged)

Result of solving an `EnsemblePBDMProblem`.
"""
struct EnsemblePBDMSolution{T}
    u::T
    converged::Bool
end

Base.length(es::EnsemblePBDMSolution) = length(es.u)
Base.getindex(es::EnsemblePBDMSolution, i) = es.u[i]
Base.iterate(es::EnsemblePBDMSolution, args...) = iterate(es.u, args...)

function Base.show(io::IO, es::EnsemblePBDMSolution)
    print(io, "EnsemblePBDMSolution($(length(es.u)) trajectories, converged=$(es.converged))")
end

function CommonSolve.solve(ep::EnsemblePBDMProblem, alg=DirectIteration(); trajectories::Int, kwargs...)
    u = deepcopy(ep.u_init)
    converged = false
    batch_data = []

    for i in 1:trajectories
        converged && break
        repeat = 1
        rerun = true
        while rerun
            modified_prob = ep.prob_func(deepcopy(ep.prob), i, repeat)
            sol = solve(modified_prob, alg; kwargs...)
            output, rerun = ep.output_func(sol, i)
            push!(batch_data, output)
            repeat += 1
            repeat > 10 && break
        end
    end

    u, converged = ep.reduction(u, batch_data, 1:trajectories)
    return EnsemblePBDMSolution(u, converged)
end

# ============================================================================
# PBDMProblem convenience constructors for coupled systems
# ============================================================================

function PBDMProblem(structure::MultiPopulationPBDM,
                     system::PopulationSystem, weather, tspan::Tuple{Int,Int};
                     approach=nothing, p=nothing,
                     rules::Vector=AbstractInteractionRule[],
                     events::Vector=AbstractScheduledEvent[],
                     observables::Vector=Observable[],
                     stress_rules::Vector=AbstractStressRule[],
                     callbacks::Vector=PhaseCallback[])
    t0, tf = tspan
    t0 <= tf || throw(ArgumentError("tspan must satisfy start <= end, got $(tspan)"))
    combined_p = (;
        user_p = p,
        rules = rules,
        events = events,
        observables = observables,
        stress_rules = stress_rules,
        callbacks = callbacks
    )
    PBDMProblem(structure, DensityDependent(), Deterministic(),
                approach, system, weather, tspan, combined_p)
end

# Convenience: infer MultiTypePBDM when not specified
function PBDMProblem(system::PopulationSystem, weather, tspan::Tuple{Int,Int}; kwargs...)
    PBDMProblem(MultiTypePBDM(), system, weather, tspan; kwargs...)
end

# ============================================================================
# Coupled solver
# ============================================================================

function _solve(::S, ::DensityDependent, ::Deterministic,
                approach,
                prob::PBDMProblem, ::DirectIteration;
                kwargs...) where {S<:MultiPopulationPBDM}
    sys = prob.populations
    sys isa PopulationSystem || throw(ArgumentError(
        "MultiPopulationPBDM solve requires a PopulationSystem, got $(typeof(sys))"))

    weather = prob.weather
    t0, tf = prob.tspan
    n_days = tf - t0 + 1
    T = Float64

    # Unpack coupled parameters
    cp = prob.p
    rules = hasfield(typeof(cp), :rules) ? cp.rules : AbstractInteractionRule[]
    events = hasfield(typeof(cp), :events) ? cp.events : AbstractScheduledEvent[]
    observables = hasfield(typeof(cp), :observables) ? cp.observables : Observable[]
    stress_rules = hasfield(typeof(cp), :stress_rules) ? cp.stress_rules : AbstractStressRule[]
    callbacks = hasfield(typeof(cp), :callbacks) ? cp.callbacks : PhaseCallback[]
    user_p = hasfield(typeof(cp), :user_p) ? cp.user_p : nothing

    # Partition callbacks by phase
    cb_pre_event  = [cb for cb in callbacks if cb.phase == PRE_EVENT]
    cb_post_event = [cb for cb in callbacks if cb.phase == POST_EVENT]
    cb_pre_step   = [cb for cb in callbacks if cb.phase == PRE_STEP]
    cb_post_step  = [cb for cb in callbacks if cb.phase == POST_STEP]
    cb_end_of_day = [cb for cb in callbacks if cb.phase == END_OF_DAY]

    comp_names = collect(keys(sys))

    # Pre-allocate per-component outputs
    comp_totals = Dict{Symbol, Vector{T}}(
        n => zeros(T, n_days) for n in comp_names)
    comp_stage_tots = Dict{Symbol, Matrix{T}}(
        n => zeros(T, n_stages(sys[n].population), n_days) for n in comp_names)

    # Observable storage
    obs_storage = Dict{Symbol, Vector{T}}(
        o.name => zeros(T, n_days) for o in observables)

    # State variable history storage
    state_names = collect(Base.keys(sys.state))
    state_history = Dict{Symbol, Vector}(
        sn => Vector{Any}(undef, n_days) for sn in state_names)

    # Event log
    event_log = Tuple{Int, Symbol}[]

    # Rule log
    rule_names = Symbol[]
    seen_names = Set{Symbol}()
    for (i, r) in enumerate(rules)
        nm = if hasproperty(r, :name)
            r.name
        else
            # Derive a readable default from the rule type (e.g. :DiapauseRule)
            Symbol(nameof(typeof(r)))
        end
        # Guarantee uniqueness: suffix with index if collision.
        if nm in seen_names
            nm = Symbol(nm, "_", i)
        end
        push!(seen_names, nm)
        push!(rule_names, nm)
    end
    rule_log = Dict{Symbol, Vector}(nm => [] for nm in rule_names)

    for d in 1:n_days
        day = t0 + d - 1
        w = get_weather(weather, day)

        # Phase 0: PRE_EVENT callbacks
        for cb in cb_pre_event
            cb.fn(sys, w, day, user_p)
        end

        # Phase 1: Apply scheduled events
        for ev in events
            fired = apply_event!(ev, sys, w, day, user_p)
            if fired
                if ev isa WeatherConditionalEvent
                    push!(event_log, (day, ev.name))
                else
                    ev_name = Symbol(typeof(ev).name.name)
                    push!(event_log, (day, ev_name))
                end
            end
        end

        # Phase 1.5: POST_EVENT callbacks
        for cb in cb_post_event
            cb.fn(sys, w, day, user_p)
        end

        # Phase 2: Update auto state variables
        for (_, sv) in sys.state
            update_state!(sv, sys, w, day, user_p)
        end

        # Phase 2.5: PRE_STEP callbacks + compute stress
        for cb in cb_pre_step
            cb.fn(sys, w, day, user_p)
        end

        component_stress = Dict{Symbol, Vector{T}}()
        for sr in stress_rules
            stress_dict = compute_stress(sr, sys, w, day, user_p)
            for (cname, svec) in stress_dict
                if haskey(component_stress, cname)
                    component_stress[cname] .= max.(component_stress[cname], svec)
                else
                    component_stress[cname] = T.(svec)
                end
            end
        end

        # Phase 3: Step each component
        for name in comp_names
            comp = sys[name]
            pop = comp.population
            if pop isa BulkPopulation
                step_bulk!(pop, w, day, user_p)
            elseif approach === nothing
                ss = get(component_stress, name, nothing)
                step_population!(pop, w; stage_stress=ss)
            else
                step_system!(pop, w, approach)
            end
        end

        # Phase 3.5: POST_STEP callbacks
        for cb in cb_post_step
            cb.fn(sys, w, day, user_p)
        end

        # Phase 4: Apply interaction rules
        for (i, rule) in enumerate(rules)
            result = apply_rule!(rule, sys, w, day, user_p)
            push!(rule_log[rule_names[i]], result)
        end

        # Phase 5: Record component state + state variable history
        for name in comp_names
            pop = sys[name].population
            comp_totals[name][d] = total_population(pop)
            if pop isa BulkPopulation
                comp_stage_tots[name][1, d] = total_population(pop)
            else
                for j in 1:n_stages(pop)
                    comp_stage_tots[name][j, d] = delay_total(pop.stages[j].delay)
                end
            end
        end

        for sn in state_names
            state_history[sn][d] = snapshot(sys.state[sn])
        end

        # Phase 6: Record observables
        for obs in observables
            obs_storage[obs.name][d] = obs.fn(sys, w, day, user_p)
        end

        # Phase 6.5: END_OF_DAY callbacks
        for cb in cb_end_of_day
            cb.fn(sys, w, day, user_p)
        end
    end

    ts = collect(t0:tf)
    return CoupledPBDMSolution{T}(
        ts, comp_names, comp_totals, comp_stage_tots,
        state_history, obs_storage, event_log, rule_log, :Success)
end
