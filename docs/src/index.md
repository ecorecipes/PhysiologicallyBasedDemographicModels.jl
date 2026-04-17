# PhysiologicallyBasedDemographicModels.jl

A Julia package for building and analyzing Physiologically Based Demographic Models (PBDMs).

## Overview

Physiologically Based Demographic Models are mechanistic, process-based simulation models
that link individual physiology to population dynamics through temperature-driven development,
resource acquisition (supply/demand), and distributed maturation delays. This package provides:

- **Distributed delay dynamics** (Manetsch/Vansickle k-substage Erlang model)
- **Temperature-driven development** (linear, Brière, Logan rate functions)
- **Supply/demand functional responses** (Frazer-Gilbert, Holling Type II/III)
- **Metabolic pool allocation** with priority-based resource partitioning
- **SciML-compatible** Problem/solve interface via StructuredPopulationCore.jl
- **Economics module** for bioeconomic analysis (costs, revenue, NPV, damage functions)
- **Genetics module** for resistance evolution (Hardy-Weinberg, selection, refuges)
- **Epidemiology module** for vector-borne disease dynamics (SIR, vector-host)
- **Multi-species interactions** (trophic webs, SIT mating models)

## Quick Start

```julia
using PhysiologicallyBasedDemographicModels

# Define a 3-stage insect lifecycle
dev = LinearDevelopmentRate(10.0, 35.0)  # base 10°C, upper 35°C

egg   = LifeStage(:egg,   DistributedDelay(10, 100.0; W0=50.0), dev, 0.05)
larva = LifeStage(:larva, DistributedDelay(15, 200.0; W0=0.0),  dev, 0.03)
adult = LifeStage(:adult, DistributedDelay(8,  150.0; W0=0.0),  dev, 0.02)

pop = Population(:insect, [egg, larva, adult])

# Simulate over 180 days with sinusoidal temperature
weather = SinusoidalWeather(22.0, 8.0)
prob = PBDMProblem(pop, weather, (1, 180))
sol = solve(prob)

println("Peak population: ", maximum(total_population(sol)))
println("Phenology (50% maturation): day ", phenology(sol))
```

## Package Architecture

PhysiologicallyBasedDemographicModels.jl follows the same trait-based dispatch pattern
as its sibling packages in the projection model ecosystem:

| Package | Domain |
|---------|--------|
| [StructuredPopulationCore.jl](https://ecorecipes.github.io/StructuredPopulationCore.jl) | Shared abstractions and analysis |
| [MatrixProjectionModels.jl](https://ecorecipes.github.io/MatrixProjectionModels.jl) | Discrete matrix population models |
| [IntegralProjectionModels.jl](https://ecorecipes.github.io/IntegralProjectionModels.jl) | Continuous-state integral projection models |
| [CategoricalPopulationDynamics.jl](https://ecorecipes.github.io/CategoricalPopulationDynamics.jl) | Categorical composition and Kan extensions |
| **PhysiologicallyBasedDemographicModels.jl** | **Process-based physiological models** |

## Tutorials

The tutorials cover a range of PBDM applications from basic concepts to
full bioeconomic systems:

These tutorial pages are generated from the canonical `vignettes/*/*.qmd`
sources during the docs build.

| # | Tutorial | System | Key Concepts |
|---|----------|--------|-------------|
| 1 | [Getting Started](tutorials/01_getting_started.md) | Generic insect | Degree-days, distributed delays, supply/demand |
| 2 | [Cotton Plant](tutorials/02_cotton_plant.md) | *Gossypium* | Multi-organ plant model, explicit BDF/MP hybrid |
| 3 | [Coffee Berry Borer](tutorials/03_coffee_berry_borer.md) | *Hypothenemus hampei* | Lactin development, 7-stage lifecycle |
| 4 | [Grapevine C/N](tutorials/04_grapevine.md) | *Vitis vinifera* | Carbon/nitrogen allocation, Q₁₀ respiration |
| 5 | [Lobesia Overwintering](tutorials/05_lobesia_overwintering.md) | *Lobesia botrana* | 3-phase diapause, photoperiod |
| 6 | [Bt Cotton Resistance](tutorials/06_bt_cotton_resistance.md) | Bt cotton + bollworm | Hardy-Weinberg genetics, refuges |
| 7 | [Pesticide Resistance](tutorials/07_pesticide_resistance.md) | Generic pest | Economic optimization + resistance |
| 8 | [Screwworm SIT](tutorials/08_screwworm_sit.md) | *Cochliomyia* | Sterile insect technique, mating competition |
| 9 | [Indian Bt Cotton](tutorials/09_bt_cotton_india.md) | Indian smallholders | Bioeconomics, rainfall-yield, NPV |
| 10 | [Tsetse Ecosocial](tutorials/10_tsetse_ecosocial.md) | Tsetse-cattle-human | Vector-borne disease, trophic coupling, welfare |
| 11 | [Olive Climate](tutorials/11_olive_climate.md) | Mediterranean olive | Climate scenarios, regional economics |
