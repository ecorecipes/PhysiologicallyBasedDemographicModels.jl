# Resilience

Resilience metrics for PBDMs come from
[ResilienceCore.jl](https://github.com/ecorecipes/ResilienceCore.jl) (operators and their
metrics) and [DynamicalResilience.jl](https://github.com/ecorecipes/DynamicalResilience.jl)
(continuation, fixed points of maps, and simulation-based recovery). This package supplies:

- an integrator for `PBDMProblem` (`init`, `step!`, `reinit!`, `set_u!`, …), so that
  DynamicalResilience's `recovery` runs on a PBDM unchanged;
- the state of a population, system or problem as one vector;
- the annual map of a model under a repeated weather year, its fixed point (the annual cycle)
  and the cycle's Floquet multipliers;
- `linearise` for the consumer–resource equilibria of the theory helpers, and the same model
  as an `ODEProblem` for continuation along temperature.

See the tutorial *Resilience of PBDMs* (vignette 65).

## Integrator

```@docs
PBDMIntegrator
init(::PBDMProblem, ::DirectIteration)
```

## State vectors

```@docs
full_state
set_full_state!
full_state_labels
```

## Annual cycles

```@docs
annual_map
annual_fixed_point
linearise(::PBDMProblem)
```

## Consumer–resource theory

```@docs
linearise(::EquilibriumResult)
consumer_resource_problem
```
