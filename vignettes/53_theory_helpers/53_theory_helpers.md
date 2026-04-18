# Theoretical Helpers — Compensation, Isoclines, and Assembly
Simon Frost

- [Overview](#overview)
- [Setup](#setup)
- [Building blocks](#building-blocks)
- [Functional-response
  introspection](#functional-response-introspection)
- [Distributed-delay rate](#distributed-delay-rate)
- [Composite biodemography
  accessors](#composite-biodemography-accessors)
- [Compensation point and life-history
  strategy](#compensation-point-and-life-history-strategy)
- [Isoclines](#isoclines)
- [Equilibrium and stability](#equilibrium-and-stability)
- [Food-web assembly](#food-web-assembly)
- [Cross-population helpers](#cross-population-helpers)
- [Summary](#summary)

## Overview

`PhysiologicallyBasedDemographicModels.jl` ships an analytical toolbox
for investigating the steady-state and topological properties of trophic
PBDMs without running full simulations:

- **Compensation point** of a respiration model — the supply/demand
  ratio at which net energetic balance is zero.
- **Life-history strategy** classification from allocation-pool
  ordering.
- **Zero-growth isoclines** of resource–consumer systems.
- **Equilibrium discovery** and **stability classification**.
- **Food-web assembly** under the temperature-mediated functional
  response.

Plus the small set of accessor / introspection helpers
(`is_ratio_dependent`, `apparency`, `delay_rate`, `biodemography`,
`allocation_model`, `acquisition_component`, `respiration_component`,
`development_component`, `compute_stress`, `fertile_mating_fraction`).

## Setup

``` julia
using PhysiologicallyBasedDemographicModels
```

## Building blocks

``` julia
fr   = FraserGilbertResponse(1.0)
resp = Q10Respiration(0.05, 2.0, 25.0)
dev  = LinearDevelopmentRate(10.0, 30.0)
acq  = fr
delay = DistributedDelay(5, 10.0)
println("fr   = ", fr)
println("resp = ", resp)
```

    fr   = FraserGilbertResponse{Float64}(1.0)
    resp = Q10Respiration{Float64}(0.05, 2.0, 25.0)

## Functional-response introspection

``` julia
println("is_ratio_dependent(FraserGilbertResponse) = ", is_ratio_dependent(fr))
println("is_ratio_dependent(HollingTypeII)         = ", is_ratio_dependent(HollingTypeII(1.0, 1.0)))
println("is_ratio_dependent(HollingTypeIII)        = ", is_ratio_dependent(HollingTypeIII(1.0, 1.0, 0.1)))
println("apparency(fr)                              = ", apparency(fr))
```

    is_ratio_dependent(FraserGilbertResponse) = true
    is_ratio_dependent(HollingTypeII)         = false
    is_ratio_dependent(HollingTypeIII)        = false
    apparency(fr)                              = 1.0

## Distributed-delay rate

``` julia
println("delay_rate(DistributedDelay(τ=10, k=5)) = ", delay_rate(delay))
```

    delay_rate(DistributedDelay(τ=10, k=5)) = 0.5

## Composite biodemography accessors

``` julia
bdf   = BiodemographicFunctions(dev, acq, resp)
pool  = MetabolicPool(2.0, [1.0, 0.5, 0.5], [:respiration, :growth, :reproduction])
model = CoupledPBDMModel(bdf, pool)

println("development_component(bdf)  isa typeof(dev)  = ",
        development_component(bdf)  isa typeof(dev))
println("acquisition_component(bdf)  isa typeof(acq)  = ",
        acquisition_component(bdf)  isa typeof(acq))
println("respiration_component(bdf)  isa typeof(resp) = ",
        respiration_component(bdf) isa typeof(resp))
println("biodemography(model)        isa typeof(bdf)  = ",
        biodemography(model)        isa typeof(bdf))
println("allocation_model(model)     isa typeof(pool) = ",
        allocation_model(model)     isa typeof(pool))
```

    development_component(bdf)  isa typeof(dev)  = true
    acquisition_component(bdf)  isa typeof(acq)  = true
    respiration_component(bdf)  isa typeof(resp) = true
    biodemography(model)        isa typeof(bdf)  = true
    allocation_model(model)     isa typeof(pool) = true

## Compensation point and life-history strategy

``` julia
phi_star = compensation_point(resp, 25.0, 1.0; conversion_efficiency = 0.8)
println("metabolic compensation point φ* = ", round(phi_star, digits=4))

println("life_history_strategy(pool)  = ", life_history_strategy(pool))
println("life_history_strategy(model) = ", life_history_strategy(model))
```

    metabolic compensation point φ* = 0.0625
    life_history_strategy(pool)  = r_selected
    life_history_strategy(model) = r_selected

## Isoclines

``` julia
M1_grid = 0.0:1.0:50.0
ci = consumer_isocline(fr, resp, 25.0;
    demand_rate = 1.0, conversion_efficiency = 0.8,
    M1_range = M1_grid)
println("typeof(ci).name.name      = ", typeof(ci).name.name)
println("ci.isocline_type          = ", ci.isocline_type)
println("ci.slope                  = ", round(ci.slope, digits=4))
println("first 4 (M₁, M₂) pairs    = ", collect(zip(ci.resource_biomass[1:4], ci.consumer_biomass[1:4])))
```

    typeof(ci).name.name      = IsoclineResult
    ci.isocline_type          = consumer
    ci.slope                  = 15.4946
    first 4 (M₁, M₂) pairs    = [(0.0, 0.0), (1.0, 15.494622163225383), (2.0, 30.989244326450766), (3.0, 46.483866489676146)]

``` julia
ri = resource_isocline(fr;
    intrinsic_rate = 0.5, carrying_capacity = 100.0,
    consumer_demand_rate = 1.0,
    M1_range = 1.0:1.0:99.0)
println("typeof(ri).name.name = ", typeof(ri).name.name)
println("ri.isocline_type     = ", ri.isocline_type)
println("ri sample length     = ", length(ri.resource_biomass))
```

    typeof(ri).name.name = IsoclineResult
    ri.isocline_type     = resource
    ri sample length     = 99

## Equilibrium and stability

``` julia
eq = find_equilibrium(fr, resp, 25.0;
    intrinsic_rate = 0.5, carrying_capacity = 100.0,
    consumer_demand_rate = 1.0, conversion_efficiency = 0.8)
println("typeof(eq).name.name = ", typeof(eq).name.name)
println("M₁* = ", round(eq.M1_star, digits=4))
println("M₂* = ", round(eq.M2_star, digits=4))
println("classification = ", eq.classification)
println("eigenvalues    = ", round.(eq.eigenvalues, digits=4))
```

    typeof(eq).name.name = EquilibriumResult
    M₁* = NaN
    M₂* = NaN
    classification = degenerate
    eigenvalues    = ComplexF64[]

`classify_equilibrium` is also usable directly on a vector of complex
eigenvalues.

``` julia
println("classify_equilibrium([-1+0im, -2+0im]) = ",
        classify_equilibrium(Complex{Float64}[-1.0+0im, -2.0+0im]))
println("classify_equilibrium([1+1im, 1-1im])   = ",
        classify_equilibrium(Complex{Float64}[1.0+1.0im, 1.0-1.0im]))
```

    classify_equilibrium([-1+0im, -2+0im]) = stable_node
    classify_equilibrium([1+1im, 1-1im])   = unstable_focus

## Food-web assembly

``` julia
sp1 = SpeciesProfile(:invader1; demand_rate = 1.0, resp = resp,
    fr = fr, conversion_efficiency = 0.8)
sp2 = SpeciesProfile(:invader2; demand_rate = 1.5, resp = resp,
    fr = FraserGilbertResponse(1.5), conversion_efficiency = 0.7)

ar = food_web_assembly([sp1, sp2], 25.0)
println("typeof(ar).name.name = ", typeof(ar).name.name)
println("assembly result      = ", ar)
```

    typeof(ar).name.name = AssemblyResult
    assembly result      = AssemblyResult{Float64}([:invader2, :invader1], [0.04761904761904763, 0.0625], [2, 1], Bool[1, 1])

## Cross-population helpers

``` julia
println("fertile_mating_fraction method count = ",
        length(methods(fertile_mating_fraction)))
println("compute_stress           method count = ",
        length(methods(compute_stress)))
```

    fertile_mating_fraction method count = 1
    compute_stress           method count = 1

## Summary

- Accessors: `is_ratio_dependent`, `apparency`, `delay_rate`,
  `biodemography`, `allocation_model`, `acquisition_component`,
  `respiration_component`, `development_component`.
- `compensation_point` + `life_history_strategy` for steady-state
  diagnostics.
- `consumer_isocline` / `resource_isocline` → `IsoclineResult`.
- `find_equilibrium` → `EquilibriumResult`; `classify_equilibrium` is
  the pure-eigenvalue helper.
- `SpeciesProfile` + `food_web_assembly` → `AssemblyResult`.
- `fertile_mating_fraction` (SIT) and `compute_stress` (`StressRule`).
