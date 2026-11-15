# Dynamics

```@docs
step_delay!
step_population!
step_system!
```

## Interaction & Lifecycle Rules

The rule types `TransferRule`, `ReproductionRule`, `MortalityRule` and `CustomRule` are
defined in StructuredPopulationCore and shared with MatrixProjectionModels.jl. Here a rule
is applied in place with `apply_rule!(rule, system, weather_day, day, p)`, and the
functions the rules hold are called as `fn(system, weather_day, day, p)`. Fractions
returned by those functions are clamped to `[0, 1]`, and negative offspring count as zero.
In a `CoupledPBDMSolution` these three rule types are logged under their type name
(`:TransferRule`, …), whatever their `name` field; other rules are logged under their
`name`, or their type name if they have none.

```@docs
apply_rule!
AbstractInteractionRule
CompetitionRule
PredationRule
MortalityRule
ReproductionRule
DiapauseRule
TransferRule
SelectionRule
StressRule
AbstractStressRule
DispersalRule
CustomRule
apply_dispersal!
```

## Scheduled Events & Phase Callbacks

`PulseRelease`, `SingleDayRelease` and `SprayEvent` are defined in StructuredPopulationCore
and shared with MatrixProjectionModels.jl. `PulseRelease` starts on day 0 unless
`start_day` is given, release amounts must be nonnegative, kill fractions must lie in
`[0, 1]`, and a `SprayEvent` naming a component that is not in the system throws when it
fires. The solver calls `apply_event!(event, system, weather_day, day, p)`; for these three
events that falls back to the shared `apply_event!(event, system, day, p)`.

```@docs
apply_event!
AbstractScheduledEvent
SprayEvent
SingleDayRelease
PulseRelease
ConditionalRelease
WeatherConditionalEvent
PhaseCallback
CallbackPhase
step_bulk!
by_patch
by_species
by_type
```
