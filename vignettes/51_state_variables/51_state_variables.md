# State Variables, Bulk Populations, and Phase Callbacks
Simon Frost

- [Overview](#overview)
- [Setup](#setup)
- [State variables](#state-variables)
- [Bulk populations](#bulk-populations)
- [Lower-level dynamics](#lower-level-dynamics)
- [Phase callbacks](#phase-callbacks)
- [Phenology and soil state](#phenology-and-soil-state)
- [Genetic state introspection](#genetic-state-introspection)
- [`remake`](#remake)
- [Summary](#summary)

## Overview

The PBDM coupled API gives every `PopulationSystem` access to typed
state containers (`ScalarState`, `ArrayState`, `DictState`) and exposes
a five-phase solver loop hookable via `PhaseCallback`. This vignette
walks through:

- the three concrete `AbstractStateVariable` implementations,
- the `update_state!` / `has_auto_update` interface,
- `BulkPopulation` and its in-place advance `step_bulk!`,
- the `CallbackPhase` enum (`PRE_EVENT`, `POST_EVENT`, `PRE_STEP`,
  `POST_STEP`, `END_OF_DAY`) and `PhaseCallback` hooks,
- the lower-level `step_population!` / `step_system!` building blocks,
- `remake` for problem variants.

## Setup

``` julia
using PhysiologicallyBasedDemographicModels
```

## State variables

``` julia
ss = ScalarState(:resource, 100.0)
println("typeof(ss).name.name = ", typeof(ss).name.name)
println("ss isa AbstractStateVariable = ", ss isa AbstractStateVariable)
println("get_state(ss) = ", get_state(ss))
set_state!(ss, 90.0)
println("after set_state!  = ", get_state(ss))
println("has_auto_update(ss) = ", has_auto_update(ss))
```

    typeof(ss).name.name = ScalarState
    ss isa AbstractStateVariable = true
    get_state(ss) = 100.0
    after set_state!  = 90.0
    has_auto_update(ss) = false

``` julia
as = ArrayState(:demand, [1.0, 2.0, 3.0])
println("typeof(as).name.name = ", typeof(as).name.name)
println("as isa AbstractStateVariable = ", as isa AbstractStateVariable)
println("get_state(as) = ", get_state(as))
```

    typeof(as).name.name = ArrayState
    as isa AbstractStateVariable = true
    get_state(as) = [1.0, 2.0, 3.0]

``` julia
ds = DictState(:patches, Dict(:north => 0.0, :south => 0.0))
println("typeof(ds).name.name = ", typeof(ds).name.name)
println("ds isa AbstractStateVariable = ", ds isa AbstractStateVariable)
println("get_state(ds) = ", get_state(ds))
```

    typeof(ds).name.name = DictState
    ds isa AbstractStateVariable = true
    get_state(ds) = Dict(:south => 0.0, :north => 0.0)

A state variable with an automatic update function:

``` julia
ss_auto = ScalarState(:t, 0.0;
    update = (cur, sys, w, day, p) -> cur + 1.0)
println("has_auto_update(ss_auto) = ", has_auto_update(ss_auto))
```

    has_auto_update(ss_auto) = true

`update_state!(sv, system, weather_day, day, p)` is what the solver
invokes — these methods are dispatched on the concrete state type.

``` julia
println("update_state! method count = ", length(methods(update_state!)))
```

    update_state! method count = 7

## Bulk populations

A `BulkPopulation` is a single scalar pool (e.g., resource biomass) with
an optional growth function and carrying capacity.

``` julia
bp = BulkPopulation(:plant_biomass, 50.0;
    growth_fn = (v, w, day, p) -> v + 0.5,
    K = 200.0)
println(bp)
step_bulk!(bp, nothing, 1, NamedTuple())
println("after step_bulk! = ", bp)
```

    BulkPopulation(:plant_biomass, 50.0)
    after step_bulk! = BulkPopulation(:plant_biomass, 50.5)

A no-growth bulk population (no `growth_fn`) is a pure carrier.

``` julia
bp2 = BulkPopulation(:soil_water, 25.0)
step_bulk!(bp2, nothing, 1, NamedTuple())
println(bp2, "  unchanged")
```

    BulkPopulation(:soil_water, 25.0)  unchanged

## Lower-level dynamics

`step_population!` and `step_system!` are the per-day building blocks
called by `solve` on a `PBDMProblem`.

``` julia
println("step_population! method count = ", length(methods(step_population!)))
println("step_system!     method count = ", length(methods(step_system!)))
```

    step_population! method count = 1
    step_system!     method count = 7

## Phase callbacks

The solver advances each day through five phases (in order):

1.  `PRE_EVENT` — before scheduled events fire
2.  `POST_EVENT` — after events, before state auto-update
3.  `PRE_STEP` — after state update + stress, before stepping
4.  `POST_STEP` — after stepping, before interaction rules
5.  `END_OF_DAY` — after observables are recorded

``` julia
println("Enum members: ", instances(CallbackPhase))
println("PRE_EVENT  isa CallbackPhase = ", PRE_EVENT  isa CallbackPhase)
println("POST_EVENT isa CallbackPhase = ", POST_EVENT isa CallbackPhase)
println("PRE_STEP   isa CallbackPhase = ", PRE_STEP   isa CallbackPhase)
println("POST_STEP  isa CallbackPhase = ", POST_STEP  isa CallbackPhase)
println("END_OF_DAY isa CallbackPhase = ", END_OF_DAY isa CallbackPhase)
```

    Enum members: (PRE_EVENT, POST_EVENT, PRE_STEP, POST_STEP, END_OF_DAY)
    PRE_EVENT  isa CallbackPhase = true
    POST_EVENT isa CallbackPhase = true
    PRE_STEP   isa CallbackPhase = true
    POST_STEP  isa CallbackPhase = true
    END_OF_DAY isa CallbackPhase = true

A `PhaseCallback` couples a name + phase + a `(sys, w, day, p)`
function:

``` julia
cb_log = PhaseCallback(:counter, END_OF_DAY,
    (sys, w, day, p) -> println("  end of day $day"))
println("typeof(cb_log).name.name = ", typeof(cb_log).name.name)
println("cb_log.phase             = ", cb_log.phase)
```

    typeof(cb_log).name.name = PhaseCallback
    cb_log.phase             = END_OF_DAY

## Phenology and soil state

The package also ships specialised state types whose update / accessor
helpers are part of the public API.

``` julia
println("get_phase           method count = ", length(methods(get_phase)))
println("past_milestone      method count = ", length(methods(past_milestone)))
println("resource_availability method count = ", length(methods(resource_availability)))
println("get_virulence       method count = ", length(methods(get_virulence)))
```

    get_phase           method count = 2
    past_milestone      method count = 1
    resource_availability method count = 1
    get_virulence       method count = 2

## Genetic state introspection

``` julia
println("get_genotypes method count = ", length(methods(get_genotypes)))
println("get_locus     method count = ", length(methods(get_locus)))
```

    get_genotypes method count = 1
    get_locus     method count = 1

## `remake`

`remake(prob; ...)` produces a copy of a `PBDMProblem` with overrides —
the canonical pattern for parameter sweeps.

``` julia
println("remake method count = ", length(methods(remake)))
```

    remake method count = 1

## Summary

- State containers: `ScalarState`, `ArrayState`, `DictState`, all
  `<: AbstractStateVariable`. `update_state!` and `has_auto_update` form
  the auto-update protocol.
- `BulkPopulation` + `step_bulk!` are scalar (non-stage-structured)
  populations.
- `step_population!` and `step_system!` are the legacy-style per-day
  drivers.
- `CallbackPhase` enum (`PRE_EVENT`, `POST_EVENT`, `PRE_STEP`,
  `POST_STEP`, `END_OF_DAY`) + `PhaseCallback` are the user-extensible
  solver hooks.
- `get_phase`, `past_milestone`, `resource_availability`,
  `get_virulence`, `get_genotypes`, `get_locus` are the typed accessors
  for `PhenologyState`, `SoilState`, and `GenomeState`.
- `remake` builds problem variants.
