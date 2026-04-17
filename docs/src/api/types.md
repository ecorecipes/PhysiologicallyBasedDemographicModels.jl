# Types & Traits

## Structure Traits

`MultiSpeciesPBDM` currently acts as a structure tag and container shape.
Direct `solve(..., DirectIteration())` support remains single-species only.

```@docs
AbstractPBDMStructure
SingleSpeciesPBDM
MultiSpeciesPBDM
```

## Approach Types

```@docs
AbstractPBDMApproach
AbstractBiodemographicApproach
AbstractAllocationApproach
AbstractHybridPBDMApproach
AbstractBiodemographicFunction
AbstractBiodemographicModel
AbstractAllocationModel
AbstractRespirationModel
BiodemographicFunctions
CoupledPBDMModel
approach_family
development_component
acquisition_component
respiration_component
biodemography
allocation_model
```

## Development Rates

`LinearDevelopmentRate` returns daily degree-day accumulation rather than a
normalized 0–1 fraction.

```@docs
AbstractDevelopmentRate
LinearDevelopmentRate
BriereDevelopmentRate
LoganDevelopmentRate
development_rate
degree_days
```

## Distributed Delay

```@docs
DistributedDelay
delay_variance
delay_rate
delay_total
```

## Functional Response

```@docs
AbstractFunctionalResponse
FraserGilbertResponse
acquire
supply_demand_ratio
```

## Metabolic Pool

```@docs
MetabolicPool
allocate
supply_demand_index
```

## Respiration

```@docs
Q10Respiration
respiration_rate
```

## Life Stages & Populations

```@docs
LifeStage
Population
n_stages
n_substages
total_population
```
