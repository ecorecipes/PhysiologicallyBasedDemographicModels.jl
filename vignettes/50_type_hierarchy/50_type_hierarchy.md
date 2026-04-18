# Type Hierarchy Reference
Simon Frost

- [Overview](#overview)
- [Setup](#setup)
- [Approach hierarchy](#approach-hierarchy)
- [Structure hierarchy](#structure-hierarchy)
- [Continuous and PSPM hierarchies](#continuous-and-pspm-hierarchies)
- [Density dependence and stochasticity
  (re-exports)](#density-dependence-and-stochasticity-re-exports)
- [Solution supertype](#solution-supertype)
- [Disease and genotype roots](#disease-and-genotype-roots)
- [Weather supertype](#weather-supertype)
- [Summary](#summary)

## Overview

`PhysiologicallyBasedDemographicModels.jl` ships a deep abstract-type
hierarchy organised along three orthogonal axes:

1.  **Approach** (`AbstractPBDMApproach`) — biodemographic functions vs.
    metabolic-pool allocation vs. hybrid combinations.
2.  **Structure** (`AbstractPBDMStructure`) — single-species,
    multi-type, multi-species, metapopulation.
3.  **Density / stochasticity** — re-exported from the shared
    `StructuredPopulationCore` package.

This short reference vignette enumerates each layer so downstream code
knows which supertype to dispatch on.

## Setup

``` julia
using PhysiologicallyBasedDemographicModels
```

## Approach hierarchy

``` julia
println("AbstractPBDMApproach            isabstract = ", isabstracttype(AbstractPBDMApproach))
println("AbstractBiodemographicApproach <: AbstractPBDMApproach = ",
        AbstractBiodemographicApproach <: AbstractPBDMApproach)
println("AbstractAllocationApproach     <: AbstractPBDMApproach = ",
        AbstractAllocationApproach     <: AbstractPBDMApproach)
println("AbstractHybridPBDMApproach     <: AbstractPBDMApproach = ",
        AbstractHybridPBDMApproach     <: AbstractPBDMApproach)
```

    AbstractPBDMApproach            isabstract = true
    AbstractBiodemographicApproach <: AbstractPBDMApproach = true
    AbstractAllocationApproach     <: AbstractPBDMApproach = true
    AbstractHybridPBDMApproach     <: AbstractPBDMApproach = true

Within the biodemographic and allocation branches:

``` julia
println("AbstractBiodemographicFunction <: AbstractBiodemographicApproach = ",
        AbstractBiodemographicFunction <: AbstractBiodemographicApproach)
println("AbstractBiodemographicModel    <: AbstractBiodemographicApproach = ",
        AbstractBiodemographicModel    <: AbstractBiodemographicApproach)
println("AbstractAllocationModel        <: AbstractAllocationApproach    = ",
        AbstractAllocationModel        <: AbstractAllocationApproach)
println("AbstractRespirationModel       <: AbstractBiodemographicFunction = ",
        AbstractRespirationModel       <: AbstractBiodemographicFunction)
println("AbstractFunctionalResponse     <: AbstractBiodemographicFunction = ",
        AbstractFunctionalResponse     <: AbstractBiodemographicFunction)
```

    AbstractBiodemographicFunction <: AbstractBiodemographicApproach = true
    AbstractBiodemographicModel    <: AbstractBiodemographicApproach = true
    AbstractAllocationModel        <: AbstractAllocationApproach    = true
    AbstractRespirationModel       <: AbstractBiodemographicFunction = true
    AbstractFunctionalResponse     <: AbstractBiodemographicFunction = true

## Structure hierarchy

``` julia
println("AbstractPBDMStructure  <: AbstractProjectionStructure = ",
        AbstractPBDMStructure <: AbstractProjectionStructure)
println("SingleSpeciesPBDM      <: AbstractPBDMStructure       = ",
        SingleSpeciesPBDM     <: AbstractPBDMStructure)
println("MultiPopulationPBDM    <: AbstractPBDMStructure       = ",
        MultiPopulationPBDM   <: AbstractPBDMStructure)
println("MetapopulationPBDM     <: MultiPopulationPBDM         = ",
        MetapopulationPBDM    <: MultiPopulationPBDM)
```

    AbstractPBDMStructure  <: AbstractProjectionStructure = true
    SingleSpeciesPBDM      <: AbstractPBDMStructure       = true
    MultiPopulationPBDM    <: AbstractPBDMStructure       = true
    MetapopulationPBDM     <: MultiPopulationPBDM         = true

Concrete instantiations:

``` julia
println("SingleSpeciesPBDM()    = ", SingleSpeciesPBDM())
println("MetapopulationPBDM()   = ", MetapopulationPBDM())
```

    SingleSpeciesPBDM()    = SingleSpeciesPBDM()
    MetapopulationPBDM()   = MetapopulationPBDM()

## Continuous and PSPM hierarchies

The continuous-time formulations have their own root.

``` julia
println("AbstractContinuousPBDM isabstract = ", isabstracttype(AbstractContinuousPBDM))
println("AbstractPSPMMethod     isabstract = ", isabstracttype(AbstractPSPMMethod))
println("AbstractPSPMSpecies    isabstract = ", isabstracttype(AbstractPSPMSpecies))
```

    AbstractContinuousPBDM isabstract = true
    AbstractPSPMMethod     isabstract = true
    AbstractPSPMSpecies    isabstract = true

## Density dependence and stochasticity (re-exports)

These are re-exported from `StructuredPopulationCore` so PBDM users do
not need to import the core package directly.

``` julia
println("AbstractDensityDependence isabstract = ", isabstracttype(AbstractDensityDependence))
println("DensityIndependent()      isa AbstractDensityDependence = ",
        DensityIndependent() isa AbstractDensityDependence)

println("AbstractStochasticity     isabstract = ", isabstracttype(AbstractStochasticity))
println("Deterministic()           isa AbstractStochasticity = ",
        Deterministic() isa AbstractStochasticity)
println("StochasticKernelResampled    isa Type = ", StochasticKernelResampled    isa Type)
println("StochasticParameterResampled isa Type = ", StochasticParameterResampled isa Type)
```

    AbstractDensityDependence isabstract = true
    DensityIndependent()      isa AbstractDensityDependence = true
    AbstractStochasticity     isabstract = true
    Deterministic()           isa AbstractStochasticity = true
    StochasticKernelResampled    isa Type = true
    StochasticParameterResampled isa Type = true

## Solution supertype

``` julia
println("AbstractProjectionSolution isabstract = ", isabstracttype(AbstractProjectionSolution))
```

    AbstractProjectionSolution isabstract = true

## Disease and genotype roots

``` julia
println("AbstractDiseaseModel  isabstract = ", isabstracttype(AbstractDiseaseModel))
println("AbstractGenotypeModel isabstract = ", isabstracttype(AbstractGenotypeModel))
```

    AbstractDiseaseModel  isabstract = true
    AbstractGenotypeModel isabstract = true

## Weather supertype

``` julia
println("AbstractWeather isabstract = ", isabstracttype(AbstractWeather))
```

    AbstractWeather isabstract = true

## Summary

- **Approach:** `AbstractPBDMApproach` →
  `AbstractBiodemographicApproach` / `AbstractAllocationApproach` /
  `AbstractHybridPBDMApproach`, with `AbstractBiodemographicFunction`,
  `AbstractBiodemographicModel`, `AbstractAllocationModel`,
  `AbstractRespirationModel`, `AbstractFunctionalResponse` underneath.
- **Structure:** `AbstractPBDMStructure` (a subtype of
  `AbstractProjectionStructure`) → `SingleSpeciesPBDM`,
  `MultiPopulationPBDM` → `MetapopulationPBDM` (and `MultiTypePBDM`,
  `MultiSpeciesPBDMNew`).
- **Continuous-time:** `AbstractContinuousPBDM`, `AbstractPSPMMethod`,
  `AbstractPSPMSpecies`.
- **Cross-cutting:** `AbstractDensityDependence` (`DensityIndependent`),
  `AbstractStochasticity` (`Deterministic`, `StochasticKernelResampled`,
  `StochasticParameterResampled`), `AbstractProjectionSolution`,
  `AbstractDiseaseModel`, `AbstractGenotypeModel`, `AbstractWeather`.
