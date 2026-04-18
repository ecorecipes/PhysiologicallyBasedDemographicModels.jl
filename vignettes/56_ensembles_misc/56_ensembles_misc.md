# Ensembles, Filters, and Weather Helpers
Simon Frost

- [Overview](#overview)
- [Setup](#setup)
- [Population system filters](#population-system-filters)
- [Weather: single-sine degree days](#weather-single-sine-degree-days)
- [Per-stage degree days from a solved
  `PBDMProblem`](#per-stage-degree-days-from-a-solved-pbdmproblem)
- [Ensemble parameter sweep](#ensemble-parameter-sweep)
- [Summary](#summary)

## Overview

Final corner of the PBDM API:

- `EnsemblePBDMProblem` / `EnsemblePBDMSolution` — SciML-style parameter
  sweeps on top of any `PBDMProblem`.
- `by_species`, `by_type`, `by_patch` — filter helpers for
  `PopulationSystem`.
- `degree_days_sine` — single-sine daily degree-day estimator.
- `stage_degree_days` — per-stage daily degree-days, accessible from a
  `PBDMSolution`.

## Setup

``` julia
using PhysiologicallyBasedDemographicModels
using CommonSolve: solve
```

## Population system filters

Build a small system spanning two species, two types, and two patches:

``` julia
make_pop(name) = begin
    dev = LinearDevelopmentRate(10.0, 35.0)
    Population(name, [
        LifeStage(:juv,   DistributedDelay(5, 30.0;  W0=10.0), dev, 0.01),
        LifeStage(:adult, DistributedDelay(1, 100.0; W0=0.0),  dev, 0.02),
    ])
end

pestN  = PopulationComponent(make_pop(:pest);  species=:pest,    type=:wild,    patch=:north)
pestS  = PopulationComponent(make_pop(:pest);  species=:pest,    type=:wild,    patch=:south)
predN  = PopulationComponent(make_pop(:pred);  species=:predator, type=:natural, patch=:north)

sys = PopulationSystem(
    :pest_n => pestN, :pest_s => pestS, :pred_n => predN)

println("by_species(sys, :pest)  count = ", length(by_species(sys, :pest)))
println("by_species(sys, :predator) count = ", length(by_species(sys, :predator)))
println("by_type(sys, :wild)     count = ", length(by_type(sys, :wild)))
println("by_patch(sys, :north)   count = ", length(by_patch(sys, :north)))
println("by_patch(sys, :south)   count = ", length(by_patch(sys, :south)))
```

    by_species(sys, :pest)  count = 2
    by_species(sys, :predator) count = 1
    by_type(sys, :wild)     count = 2
    by_patch(sys, :north)   count = 2
    by_patch(sys, :south)   count = 1

## Weather: single-sine degree days

`degree_days_sine(T_min, T_max, T_lower; T_upper=Inf)` returns the daily
heat units accumulated above the lower threshold using the single-sine
approximation.

``` julia
println("DD(T_min=8, T_max=24, T_lower=10) = ",
        round(degree_days_sine(8.0, 24.0, 10.0); digits = 3))
println("DD(T_min=8, T_max=24, T_lower=10, T_upper=22) = ",
        round(degree_days_sine(8.0, 24.0, 10.0; T_upper = 22.0); digits = 3))
println("DD(T_min=5, T_max=9,  T_lower=10) = ",
        round(degree_days_sine(5.0, 9.0, 10.0); digits = 3))
```

    DD(T_min=8, T_max=24, T_lower=10) = 6.304
    DD(T_min=8, T_max=24, T_lower=10, T_upper=22) = 6.0
    DD(T_min=5, T_max=9,  T_lower=10) = 0.0

## Per-stage degree days from a solved `PBDMProblem`

`PBDMSolution` exposes `stage_degree_days` as a `(n_stages × n_days)`
matrix of per-stage degree-days accumulated each day, useful for
post-hoc phenology accounting.

``` julia
pop_pest = make_pop(:pest)
weather  = WeatherSeries(fill(25.0, 60); day_offset=1)
prob     = PBDMProblem(pop_pest, weather, (1, 60))

sol_solo = solve(prob, DirectIteration())
println("typeof(sol_solo).name.name = ", typeof(sol_solo).name.name)
println("size(sol_solo.stage_degree_days) = ", size(sol_solo.stage_degree_days))
println("total stage-degree-days summed   = ",
        round(sum(sol_solo.stage_degree_days); digits = 2))
```

    typeof(sol_solo).name.name = PBDMSolution
    size(sol_solo.stage_degree_days) = (2, 60)
    total stage-degree-days summed   = 1800.0

## Ensemble parameter sweep

`EnsemblePBDMProblem` wraps a template problem with three SciML-style
hooks: `prob_func` (modify the template per trajectory), `output_func`
(extract per-solution output), and `reduction` (accumulate batch
results).

``` julia
ens = EnsemblePBDMProblem(prob;
    prob_func = (p, i, _) -> p,
    output_func = (sol, i) -> (i, false))

ens_sol = solve(ens, DirectIteration(); trajectories = 4)
println("typeof(ens).name.name      = ", typeof(ens).name.name)
println("typeof(ens_sol).name.name  = ", typeof(ens_sol).name.name)
println("length(ens_sol)            = ", length(ens_sol))
println("ens_sol.converged          = ", ens_sol.converged)
println("ens_sol[1]                 = ", ens_sol[1])
```

    typeof(ens).name.name      = EnsemblePBDMProblem
    typeof(ens_sol).name.name  = EnsemblePBDMSolution
    length(ens_sol)            = 4
    ens_sol.converged          = false
    ens_sol[1]                 = 1

## Summary

- Filters on `PopulationSystem`: `by_species`, `by_type`, `by_patch`.
- Weather helper: `degree_days_sine` (with optional upper threshold).
- Phenology bookkeeping: `stage_degree_days` field on `PBDMSolution`.
- Ensembles: `EnsemblePBDMProblem`, `EnsemblePBDMSolution`, both
  iterable / indexable for downstream analysis.
