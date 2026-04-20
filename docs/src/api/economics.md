# Economics

## Cost Functions

```@docs
AbstractCostFunction
FixedCost
VariableCost
InputCostBundle
total_cost
```

## Revenue

```@docs
AbstractRevenueFunction
CropRevenue
revenue
```

## Damage Functions

```@docs
AbstractDamageFunction
LinearDamageFunction
ExponentialDamageFunction
yield_loss
actual_yield
```

## Profit & Valuation

```@docs
net_profit
daily_income
npv
benefit_cost_ratio
```

## Yield Models

```@docs
RainfallYieldModel
WeatherYieldModel
predict_yield
```

## Surrogate Models

[`LogLinearSurrogate`](@ref) packages an offline-fit log-linear regression
surrogate for PBDM-derived yield, pest pressure, or treatment-utility
maps. It is used by vignettes 58 (CBB bioeconomics) and 60
(Tuta absoluta invasion) to summarise multi-parameter sweeps.

```@docs
LogLinearSurrogate
predict_log
predict
marginal_effects
```

## Management & Optimal Control

```@docs
AbstractManagementAction
AbstractManagementObjective
BiologicalReleaseControl
HarvestControl
PesticideControl
MaximizeProfit
MinimizeDamage
PBDMControlProblem
ManagementSolution
optimize_management
```
