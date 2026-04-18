# Extended Rules and Scheduled Events
Simon Frost

- [Overview](#overview)
- [Setup](#setup)
- [TransferRule](#transferrule)
- [PredationRule](#predationrule)
- [CompetitionRule](#competitionrule)
- [ConditionalRelease](#conditionalrelease)
- [WeatherConditionalEvent](#weatherconditionalevent)
- [SprayEvent](#sprayevent)
- [DispersalRule and SpatialGrid](#dispersalrule-and-spatialgrid)
- [`apply_event!` dispatch surface](#apply_event-dispatch-surface)
- [Summary](#summary)

## Overview

The earlier corpus vignettes use the basic `TransferRule` /
`ReproductionRule` / `MortalityRule` / `PulseRelease` /
`SingleDayRelease` patterns. This vignette covers the **extended** rule
and event types introduced for ecological interactions, conditional
triggers, and spatial dispersal:

- `TransferRule` (also covered here as the rule prototype),
- `PredationRule` and `CompetitionRule`,
- `ConditionalRelease` and `WeatherConditionalEvent`,
- `SprayEvent`,
- `DispersalRule` + `apply_dispersal!`,
- the `apply_event!` dispatch surface.

## Setup

``` julia
using PhysiologicallyBasedDemographicModels
```

## TransferRule

``` julia
tr = TransferRule(:eggs, :larvae, (sys, w, day, p) -> 0.05;
    source_stage = 1, target_stage = 1)
println("typeof(tr).name.name = ", typeof(tr).name.name)
println("tr isa AbstractInteractionRule = ", tr isa AbstractInteractionRule)
```

    typeof(tr).name.name = TransferRule
    tr isa AbstractInteractionRule = true

## PredationRule

A `PredationRule` couples a predator and prey component via a functional
response (`HollingTypeII`, `HollingTypeIII`, or
`FraserGilbertResponse`).

``` julia
fr_pred = HollingTypeIII(0.8, 0.5, 0.1)
pr = PredationRule(:wasp, :aphid, fr_pred;
    conversion = 0.1, predator_stage = -1, prey_stage = 0)
println("typeof(pr).name.name = ", typeof(pr).name.name)
println("pr isa AbstractInteractionRule = ", pr isa AbstractInteractionRule)
```

    typeof(pr).name.name = PredationRule
    pr isa AbstractInteractionRule = true

## CompetitionRule

``` julia
cr = CompetitionRule(
    [:species_a, :species_b],
    (sys, w, day, p) -> 100.0,
    (demands, supply) -> demands .* (supply / max(sum(demands), eps())))
println("typeof(cr).name.name = ", typeof(cr).name.name)
println("cr isa AbstractInteractionRule = ", cr isa AbstractInteractionRule)
```

    typeof(cr).name.name = CompetitionRule
    cr isa AbstractInteractionRule = true

## ConditionalRelease

A scheduled release of natural enemies that fires every `interval` days
provided the amount returned by `amount_fn` is positive.

``` julia
cr2 = ConditionalRelease(:parasitoid,
    (sys, w, day, p) -> day > 30 ? 50.0 : 0.0,
    7; stage_idx = 1, start_day = 1)
println("typeof(cr2).name.name = ", typeof(cr2).name.name)
println("cr2 isa AbstractScheduledEvent = ", cr2 isa AbstractScheduledEvent)
```

    typeof(cr2).name.name = ConditionalRelease
    cr2 isa AbstractScheduledEvent = true

## WeatherConditionalEvent

Fires whenever a user-supplied `(sys, w, day, p) -> Bool` predicate
returns true. A common idiom: trigger diapause / emergence on a
photoperiod or temperature threshold.

``` julia
wce = WeatherConditionalEvent(
    :cold_snap_kill,
    (sys, w, day, p) -> (w !== nothing && hasproperty(w, :T_min) && w.T_min < -2.0),
    (sys, w, day, p) -> remove_fraction!(sys, :pest, 0.5))
println("typeof(wce).name.name = ", typeof(wce).name.name)
println("wce isa AbstractScheduledEvent = ", wce isa AbstractScheduledEvent)
```

    typeof(wce).name.name = WeatherConditionalEvent
    wce isa AbstractScheduledEvent = true

## SprayEvent

Multiplicative mortality on listed days for one or more components.

``` julia
spray = SprayEvent([:pest], [0.85], [60, 90])
println("typeof(spray).name.name = ", typeof(spray).name.name)
println("spray isa AbstractScheduledEvent = ", spray isa AbstractScheduledEvent)
```

    typeof(spray).name.name = SprayEvent
    spray isa AbstractScheduledEvent = true

## DispersalRule and SpatialGrid

Dispersal connects components by name across patches in a `SpatialGrid`.

``` julia
disp = DispersalRule(:pest;
    emigration_fn = (N, phi, w, day, p) -> 0.05,
    patch_finding = 0.5)
println("typeof(disp).name.name = ", typeof(disp).name.name)
println("disp isa AbstractInteractionRule = ", disp isa AbstractInteractionRule)
```

    typeof(disp).name.name = DispersalRule
    disp isa AbstractInteractionRule = true

`apply_dispersal!(rule, grid, w, day, p)` is the per-day dispersal
driver called by the spatial solver loop.

``` julia
println("apply_dispersal! method count = ", length(methods(apply_dispersal!)))
```

    apply_dispersal! method count = 1

A trivial `SpatialGrid` shows the constructor surface (full multi-patch
demos live in the corpus vignettes).

``` julia
sg_methods = length(methods(SpatialGrid))
println("SpatialGrid constructor count = ", sg_methods)
```

    SpatialGrid constructor count = 3

## `apply_event!` dispatch surface

``` julia
println("apply_event! method count = ", length(methods(apply_event!)))
```

    apply_event! method count = 5

Each of `PulseRelease`, `SingleDayRelease`, `SprayEvent`,
`ConditionalRelease`, `WeatherConditionalEvent` (and any user-defined
`AbstractScheduledEvent`) is dispatched by this generic.

## Summary

- Extended interaction rules: `TransferRule`, `PredationRule` (driven by
  any `AbstractFunctionalResponse` such as `HollingTypeIII`),
  `CompetitionRule`, `DispersalRule`. All `<: AbstractInteractionRule`.
- Extended scheduled events: `ConditionalRelease`,
  `WeatherConditionalEvent`, `SprayEvent`. All
  `<: AbstractScheduledEvent`.
- Spatial dispersal: `SpatialGrid` + `DispersalRule` +
  `apply_dispersal!`.
- The `apply_event!` generic is the per-day dispatch hook for every
  scheduled event subtype.
