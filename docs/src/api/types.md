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

`PopulationComponent`, `PopulationSystem`, `AbstractStateVariable`, `Observable` and the
functions that read and change a system (`get_state`, `set_state!`, `has_state`,
`component_total`, `component_totals`, `inject!`, `remove_fraction!`) are defined in
StructuredPopulationCore and shared with MatrixProjectionModels.jl, so they mean the same
thing when both packages are loaded. This package exports them and adds the methods for
its own population carriers (`Population`, `BulkPopulation`) and state variable types.

A component made from a `Population` or `BulkPopulation` holds it, by reference, in
`component.population`; `component.model` is `nothing` and `component.stage_names` are the
names of its stages. A system's state is a dictionary from names to state variables (held
by reference) or plain values. `get_state(system, name)` returns the value of a state
variable, or a plain value as it is; plain values are not updated automatically by the
solver.

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
