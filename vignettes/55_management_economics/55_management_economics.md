# Management Optimisation and Economics
Simon Frost

- [Overview](#overview)
- [Setup](#setup)
- [Trophic level](#trophic-level)
- [Management actions](#management-actions)
- [Objectives](#objectives)
- [Control problem](#control-problem)
- [Solver hook](#solver-hook)
- [Economics — costs](#economics--costs)
- [Sterile insect technique (SIT)](#sterile-insect-technique-sit)
- [Economics — revenue and damage](#economics--revenue-and-damage)
- [Summary](#summary)

## Overview

The optimal-control surface in
`PhysiologicallyBasedDemographicModels.jl` ships:

- A simplified `TrophicLevel` parameterisation suitable for the
  Euler-discretised optimisation (no distributed-delay substages).
- A handful of management actions (`PesticideControl`,
  `BiologicalReleaseControl`, `HarvestControl`) all
  `<: AbstractManagementAction`.
- Two objectives (`MaximizeProfit`, `MinimizeDamage`) both
  `<: AbstractManagementObjective`.
- The `PBDMControlProblem` and `ManagementSolution` containers and the
  `optimize_management` solver hook (provided by the JuMP package
  extension).

The economics layer adds `AbstractCostFunction`,
`AbstractRevenueFunction`, `AbstractDamageFunction`, with `FixedCost` /
`VariableCost` as the primary cost-function concretes.

## Setup

``` julia
using PhysiologicallyBasedDemographicModels
```

## Trophic level

``` julia
fr = FraserGilbertResponse(1.0)
resp = Q10Respiration(0.05, 2.0, 25.0)

plant = TrophicLevel(:plant;
    demand_rate = 0.5, fr = fr, resp = resp,
    intrinsic_rate = 0.1, carrying_capacity = 1000.0)
pest  = TrophicLevel(:pest;
    demand_rate = 0.3, fr = fr, resp = resp,
    conversion_efficiency = 0.4)

println("typeof(plant).name.name = ", typeof(plant).name.name)
```

    typeof(plant).name.name = TrophicLevel

## Management actions

``` julia
spray = PesticideControl(name = :spray, target = 2,
    max_rate = 0.5, efficacy = 1.0, cost_weight = 1.0)
release = BiologicalReleaseControl(name = :predator_release, target = 3,
    max_rate = 100.0, cost_weight = 0.5)
harvest = HarvestControl(name = :crop_harvest, target = 1,
    max_rate = 100.0, revenue_per_unit = 2.5, cost_weight = 0.1)

for a in (spray, release, harvest)
    println(rpad(string(typeof(a).name.name), 28),
            " <: AbstractManagementAction = ", a isa AbstractManagementAction)
end
```

    PesticideControl             <: AbstractManagementAction = true
    BiologicalReleaseControl     <: AbstractManagementAction = true
    HarvestControl               <: AbstractManagementAction = true

## Objectives

``` julia
obj_min = MinimizeDamage(damage_weights = [0.0, 1.0],
    control_cost_weight = 1.0)
obj_max = MaximizeProfit(resource_level = 1, price_per_unit = 5.0,
    control_cost_weight = 1.0, discount_rate = 0.05)

println("typeof(obj_min).name.name = ", typeof(obj_min).name.name)
println("typeof(obj_max).name.name = ", typeof(obj_max).name.name)
println("obj_min isa AbstractManagementObjective = ",
        obj_min isa AbstractManagementObjective)
println("obj_max isa AbstractManagementObjective = ",
        obj_max isa AbstractManagementObjective)
```

    typeof(obj_min).name.name = MinimizeDamage
    typeof(obj_max).name.name = MaximizeProfit
    obj_min isa AbstractManagementObjective = true
    obj_max isa AbstractManagementObjective = true

## Control problem

``` julia
prob = PBDMControlProblem(
    levels = [plant, pest],
    controls = AbstractManagementAction[spray],
    objective = obj_min,
    u0 = [500.0, 50.0],
    tspan = (0.0, 120.0),
    dt = 1.0,
    T_celsius = 25.0)
println("typeof(prob).name.name = ", typeof(prob).name.name)
```

    typeof(prob).name.name = PBDMControlProblem

## Solver hook

`optimize_management(prob; ...)` is provided by the JuMP package
extension; it returns a `ManagementSolution`.

``` julia
println("ManagementSolution    isdefined = ", isdefined(@__MODULE__, :ManagementSolution))
println("optimize_management   method count = ", length(methods(optimize_management)))
```

    ManagementSolution    isdefined = true
    optimize_management   method count = 0

## Economics — costs

``` julia
fixed = FixedCost(150.0, :rent)
variable = VariableCost(12.0, :pesticide)
println("fixed    isa AbstractCostFunction = ", fixed    isa AbstractCostFunction)
println("variable isa AbstractCostFunction = ", variable isa AbstractCostFunction)

bundle = InputCostBundle(seed = 200.0, pesticide = 60.0, labour = 350.0)
println("total_cost(bundle) = ", total_cost(bundle))
```

    fixed    isa AbstractCostFunction = true
    variable isa AbstractCostFunction = true
    total_cost(bundle) = 610.0

``` julia
println("total_cost on (FixedCost, VariableCost) = ",
        total_cost([fixed, variable], Dict(:pesticide => 5.0)))
```

    total_cost on (FixedCost, VariableCost) = 210.0

## Sterile insect technique (SIT)

`SITRelease` packages the parameters for a sterile-insect release
programme and `fertile_mating_fraction` returns the fraction of matings
that remain fertile under continued release.

``` julia
sit = SITRelease(10_000.0, 0.6, 7)  # 10k sterile males, 60% competitive, weekly
println("typeof(sit).name.name = ", typeof(sit).name.name)
println("fertile_mating_fraction(1000, sit, 0)  = ",
        round(fertile_mating_fraction(1000.0, sit, 0); digits = 4))
println("fertile_mating_fraction(1000, sit, 3)  = ",
        round(fertile_mating_fraction(1000.0, sit, 3); digits = 4))
```

    typeof(sit).name.name = SITRelease
    fertile_mating_fraction(1000, sit, 0)  = 0.1429
    fertile_mating_fraction(1000, sit, 3)  = 0.1861

## Economics — revenue and damage

``` julia
println("AbstractRevenueFunction isabstract = ",
        isabstracttype(AbstractRevenueFunction))
println("AbstractDamageFunction  isabstract = ",
        isabstracttype(AbstractDamageFunction))
```

    AbstractRevenueFunction isabstract = true
    AbstractDamageFunction  isabstract = true

## Summary

- Bioeconomic optimisation surface: `TrophicLevel`,
  `AbstractManagementAction` (`PesticideControl`,
  `BiologicalReleaseControl`, `HarvestControl`),
  `AbstractManagementObjective` (`MinimizeDamage`, `MaximizeProfit`),
  `PBDMControlProblem`, `ManagementSolution`, `optimize_management`.
- Economics surface: `AbstractCostFunction` (`FixedCost`,
  `VariableCost`), `AbstractRevenueFunction`, `AbstractDamageFunction`.
