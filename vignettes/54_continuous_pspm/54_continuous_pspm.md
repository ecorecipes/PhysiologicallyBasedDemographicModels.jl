# Continuous-Time and PSPM Formulations
Simon Frost

- [Overview](#overview)
- [Setup](#setup)
- [Supertypes](#supertypes)
- [Building a `ContinuousSpecies`](#building-a-continuousspecies)
- [Trophic links](#trophic-links)
- [Aggregate-biomass ODE problem](#aggregate-biomass-ode-problem)
- [Delay DDE problem](#delay-dde-problem)
- [PSPM discretisation methods](#pspm-discretisation-methods)
- [Flattening helpers (discrete → continuous
  bridge)](#flattening-helpers-discrete--continuous-bridge)
- [`ContinuousPBDMSolution` and
  `species_trajectory`](#continuouspbdmsolution-and-species_trajectory)
- [Summary](#summary)

## Overview

`PhysiologicallyBasedDemographicModels.jl` ships three continuous-time
formulations of the same physiologically-based dynamics:

1.  **`ContinuousPBDMProblem`** — the linear-chain-trick ODE expansion
    of the discrete-substage model.
2.  **`DelayPBDMProblem`** — a delay differential equation formulation
    with explicit distributed-delay history.
3.  **`PSPMProblem`** — Physiologically Structured Population Model
    (McKendrick–von Foerster PDE) discretised by a chosen
    `AbstractPSPMMethod`.

All three share the `AbstractContinuousPBDM` supertype, and the
`solve_continuous` / `solve_delay` generics are wired up by the
`OrdinaryDiffEq` and `DelayDiffEq` package extensions.

This vignette enumerates the API surface and demonstrates introspection
helpers that work without loading the SciML extensions.

## Setup

``` julia
using PhysiologicallyBasedDemographicModels
```

## Supertypes

``` julia
println("AbstractContinuousPBDM isabstract = ", isabstracttype(AbstractContinuousPBDM))
println("AbstractPSPMMethod     isabstract = ", isabstracttype(AbstractPSPMMethod))
println("AbstractPSPMSpecies    isabstract = ", isabstracttype(AbstractPSPMSpecies))
```

    AbstractContinuousPBDM isabstract = true
    AbstractPSPMMethod     isabstract = true
    AbstractPSPMSpecies    isabstract = true

## Building a `ContinuousSpecies`

``` julia
fr   = FraserGilbertResponse(1.0)
resp = Q10Respiration(0.05, 2.0, 25.0)
dev  = LinearDevelopmentRate(10.0, 30.0)

sp_consumer = ContinuousSpecies(:consumer;
    k = [3, 5, 4],
    τ = [50.0, 100.0, 80.0],
    μ = [0.001, 0.0005, 0.001],
    dev_rate = [dev, dev, dev],
    fr = fr,
    resp = resp,
    demand_rate = 1.0,
    conversion_efficiency = 0.6,
    intrinsic_rate = 0.0,
    carrying_capacity = Inf)

sp_resource = ContinuousSpecies(:resource;
    k = [1],
    τ = [10.0],
    μ = [0.0],
    dev_rate = [dev],
    fr = fr,
    resp = resp,
    demand_rate = 1.0,
    conversion_efficiency = 1.0,
    intrinsic_rate = 0.5,
    carrying_capacity = 100.0)

println("typeof(sp_consumer).name.name = ", typeof(sp_consumer).name.name)
```

    typeof(sp_consumer).name.name = ContinuousSpecies

## Trophic links

``` julia
link = ContinuousTrophicLink(:consumer, :resource)
println("typeof(link).name.name = ", typeof(link).name.name)
println("link.consumer = ", link.consumer)
println("link.resource = ", link.resource)
```

    typeof(link).name.name = ContinuousTrophicLink
    link.consumer = consumer
    link.resource = resource

## Aggregate-biomass ODE problem

``` julia
total_substages = sum(sum(sp.k) for sp in (sp_consumer, sp_resource))
u0 = zeros(total_substages)
u0[1] = 50.0   # seed the first consumer substage
u0[end] = 50.0 # seed the resource substage

cprob = ContinuousPBDMProblem(
    species = [sp_consumer, sp_resource],
    links   = [link],
    u0      = u0,
    tspan   = (0.0, 50.0),
    T_forcing = 25.0)
println("typeof(cprob).name.name           = ", typeof(cprob).name.name)
println("cprob isa AbstractContinuousPBDM  = ", cprob isa AbstractContinuousPBDM)
println("cprob.tspan                       = ", cprob.tspan)
```

    typeof(cprob).name.name           = ContinuousPBDMProblem
    cprob isa AbstractContinuousPBDM  = true
    cprob.tspan                       = (0.0, 50.0)

`species_state_ranges` and `species_total_ranges` map the flat state
vector back to per-species / per-stage views.

``` julia
ranges = species_state_ranges(cprob)
println("species_state_ranges (first 3) = ", ranges[1:min(3, end)])
totals = species_total_ranges(cprob)
println("species_total_ranges = ", totals)
```

    species_state_ranges (first 3) = Tuple{Symbol, Int64, UnitRange{Int64}}[(:consumer, 1, 1:3), (:consumer, 2, 4:8), (:consumer, 3, 9:12)]
    species_total_ranges = Dict{Symbol, UnitRange{Int64}}(:consumer => 1:12, :resource => 13:13)

`solve_continuous` is the dispatch hook implemented by the
`OrdinaryDiffEq` extension.

``` julia
println("solve_continuous method count = ", length(methods(solve_continuous)))
```

    solve_continuous method count = 0

## Delay DDE problem

``` julia
n_stage_total = sum(sp.n_stages for sp in (sp_consumer, sp_resource))
u0_delay = fill(10.0, n_stage_total)

dprob = DelayPBDMProblem(
    species = [sp_consumer, sp_resource],
    links   = [link],
    u0      = u0_delay,
    tspan   = (0.0, 50.0),
    T_forcing = 25.0)
println("typeof(dprob).name.name           = ", typeof(dprob).name.name)
println("dprob isa AbstractContinuousPBDM  = ", dprob isa AbstractContinuousPBDM)
println("solve_delay method count          = ", length(methods(solve_delay)))
```

    typeof(dprob).name.name           = DelayPBDMProblem
    dprob isa AbstractContinuousPBDM  = true
    solve_delay method count          = 0

## PSPM discretisation methods

``` julia
methods_pspm = [
    FixedMeshUpwind, EscalatorBoxcarTrain, CharacteristicMethod,
    ImplicitFixedMeshUpwind, ImplicitEscalatorBoxcarTrain,
    ImplicitCharacteristicMethod, LaxFriedrichsUpwind,
    ImplicitLaxFriedrichsUpwind,
]
for M in methods_pspm
    println(rpad(string(M), 32), " <: AbstractPSPMMethod = ", M <: AbstractPSPMMethod)
end
```

    PhysiologicallyBasedDemographicModels.FixedMeshUpwind <: AbstractPSPMMethod = true
    PhysiologicallyBasedDemographicModels.EscalatorBoxcarTrain <: AbstractPSPMMethod = true
    PhysiologicallyBasedDemographicModels.CharacteristicMethod <: AbstractPSPMMethod = true
    PhysiologicallyBasedDemographicModels.ImplicitFixedMeshUpwind <: AbstractPSPMMethod = true
    PhysiologicallyBasedDemographicModels.ImplicitEscalatorBoxcarTrain <: AbstractPSPMMethod = true
    PhysiologicallyBasedDemographicModels.ImplicitCharacteristicMethod <: AbstractPSPMMethod = true
    PhysiologicallyBasedDemographicModels.LaxFriedrichsUpwind <: AbstractPSPMMethod = true
    PhysiologicallyBasedDemographicModels.ImplicitLaxFriedrichsUpwind <: AbstractPSPMMethod = true

## Flattening helpers (discrete → continuous bridge)

`flatten_population` / `flatten_populations` extract the substage masses
from a discrete `Population` to seed `u0` for a continuous problem.

``` julia
println("flatten_population  method count = ", length(methods(flatten_population)))
println("flatten_populations method count = ", length(methods(flatten_populations)))
```

    flatten_population  method count = 1
    flatten_populations method count = 1

## `ContinuousPBDMSolution` and `species_trajectory`

The wrapper type emitted by `solve_continuous`;
`species_trajectory(sol, :name)` returns the per-species total-biomass
trajectory.

``` julia
println("ContinuousPBDMSolution exists  = ", isdefined(@__MODULE__, :ContinuousPBDMSolution))
println("species_trajectory  method count = ", length(methods(species_trajectory)))
```

    ContinuousPBDMSolution exists  = true
    species_trajectory  method count = 1

## Summary

- Supertypes: `AbstractContinuousPBDM`, `AbstractPSPMMethod`,
  `AbstractPSPMSpecies`.
- Species: `ContinuousSpecies`, `ContinuousTrophicLink`.
- Problems: `ContinuousPBDMProblem`, `DelayPBDMProblem` (PSPMProblem
  also exists for the PDE formulation).
- PSPM discretisations include `ImplicitCharacteristicMethod` and
  `ImplicitEscalatorBoxcarTrain` alongside the explicit variants.
- Discrete → continuous bridge: `flatten_population`,
  `flatten_populations`.
- Inspection helpers: `species_state_ranges`, `species_total_ranges`,
  `species_trajectory`.
- Solver hooks: `solve_continuous`, `solve_delay` (provided by SciML
  extensions).
- Solution wrapper: `ContinuousPBDMSolution`.
