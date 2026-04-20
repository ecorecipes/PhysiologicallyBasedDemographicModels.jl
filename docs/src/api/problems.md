# Problem & Solution

`DirectIteration()` is the supported solver algorithm for `PBDMProblem`.
`EigenAnalysis()` is intentionally unsupported for nonlinear PBDM dynamics,
and direct multi-species `solve` support is not implemented yet.

```@docs
PBDMProblem
PBDMSolution
```

## Continuous-Time & PSPM Problems

```@docs
AbstractContinuousPBDM
ContinuousPBDMProblem
ContinuousPBDMSolution
ContinuousSpecies
ContinuousTrophicLink
DelayPBDMProblem
PSPMProblem
PSPMSpecies
PSPMStage
StagedPSPMSpecies
AbstractPSPMMethod
AbstractPSPMSpecies
n_pspm_stages
staged_species_stage_totals
species_state_ranges
species_total_ranges
species_trajectory
variables
solve_continuous
solve_delay
solve_pspm
```

## Internal Helpers

```@docs
PhysiologicallyBasedDemographicModels._total_substages
PhysiologicallyBasedDemographicModels._get_temperature
```

## Numerical Methods

```@docs
CharacteristicMethod
EscalatorBoxcarTrain
FixedMeshUpwind
LaxFriedrichsUpwind
ImplicitCharacteristicMethod
ImplicitEscalatorBoxcarTrain
ImplicitFixedMeshUpwind
ImplicitLaxFriedrichsUpwind
```

## Ensemble & Scenario

```@docs
EnsemblePBDMProblem
EnsemblePBDMSolution
run_scenarios
remake
```
