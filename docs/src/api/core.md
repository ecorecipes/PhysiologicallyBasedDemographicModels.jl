# Core Abstractions (StructuredPopulationCore)

## Type Hierarchy

```@docs
AbstractProjectionSolution
AbstractStateDomain
ContinuousDomain
DiscreteDomain
AbstractStateSemantics
ContinuousState
FiniteState
GeneralContinuousState
SimpleContinuousState
AbstractTimeSemantics
ContinuousTime
DiscreteTime
AbstractContinuousStateStructure
AbstractIPMStructure
SimpleIPM
GeneralIPM
StateBlockLayout
blocknames
blockrange
blockranges
bounds
meshpoints
n_states
step_size
combine_state
split_state
extract_population
```

## Eigen & Demographic Analysis

```@docs
EigenAnalysis
eigenanalysis_full
eigenanalysis_power
lambda
stable_distribution
reproductive_value
damping_ratio
sensitivity
elasticity
net_repro_rate_lagged
generation_time_lagged
stochastic_growth_rate
is_ergodic
is_irreducible
is_primitive
mean_kernel
```

## Delay & Lag Utilities

```@docs
DelayGeneratorTerm
TimeLagStructure
augment_population
expand_lag_matrix
extract_lag_components
area_under_curve
```
