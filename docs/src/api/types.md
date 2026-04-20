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
apparency
is_ratio_dependent
consumer_isocline
resource_isocline
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
BulkPopulation
n_stages
n_substages
```

## Structure Tags & Density Dependence

```@docs
AbstractProjectionStructure
AbstractDensityDependence
DensityIndependent
DensityDependent
AbstractStochasticity
Deterministic
StochasticKernelResampled
StochasticParameterResampled
DirectIteration
```

## State Variables & Containers

```@docs
AbstractStateVariable
ScalarState
ArrayState
DictState
DiapauseState
SoilState
PhenologyState
GenomeState
PopulationComponent
PopulationSystem
component_total
component_totals
flatten_population
flatten_populations
get_state
set_state!
has_state
get_genotypes
get_locus
get_phase
get_inoculum
get_virulence
resource_availability
past_milestone
has_auto_update
inject!
remove_fraction!
set_value!
snapshot
update_state!
Observable
```

## Coupled Multi-Species Structures

```@docs
MultiSpeciesPBDMNew
MultiPopulationPBDM
MultiTypePBDM
MetapopulationPBDM
SpatialGrid
CoupledPBDMSolution
Base.getindex(::CoupledPBDMSolution, ::Symbol)
```
