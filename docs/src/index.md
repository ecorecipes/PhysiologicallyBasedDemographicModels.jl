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
| 12 | [Cassava Mealybug Biocontrol](tutorials/12_cassava_mealybug.md) | Cassava–mealybug–parasitoid | Biocontrol, trophic interactions, functional response |
| 13 | [Medfly Invasion](tutorials/13_medfly_invasion.md) | Mediterranean fruit fly | Climate suitability, dispersal, establishment risk |
| 14 | [East Coast Fever](tutorials/14_east_coast_fever.md) | Cattle–tick–pathogen | Disease transmission, livestock demography, control |
| 15 | [Rice–Weed Competition](tutorials/15_rice_weed_competition.md) | Rice–weed system | Competition, crop yield, management |
| 16 | [Aedes albopictus Invasion](tutorials/16_aedes_albopictus.md) | Asian tiger mosquito | Climate limits, invasion risk, phenology |
| 17 | [Pink Bollworm Climate Limits](tutorials/17_pink_bollworm.md) | *Pectinophora gossypiella* | Temperature effects, voltinism, range limits |
| 18 | [Vine Mealybug Biocontrol](tutorials/18_vine_mealybug.md) | Grapevine–mealybug–parasitoid | Biocontrol, host–parasitoid dynamics, phenology |
| 19 | [BMSB Biocontrol](tutorials/19_bmsb_biocontrol.md) | Stink bug–parasitoid | Biological control, host–parasitoid interactions |
| 20 | [Cabbage Root Fly Diapause](tutorials/20_cabbage_maggot.md) | *Delia radicum* | Diapause, phenology, temperature-driven development |
| 21 | [Tomato IPM](tutorials/21_tomato_ipm.md) | Processing tomato | IPM, pest management, crop yield |
| 22 | [Tuta absoluta Invasion](tutorials/22_tuta_absoluta.md) | Tomato leafminer | Invasion risk, climate suitability, spread |
| 23 | [Yellow Starthistle Biocontrol](tutorials/23_yellow_starthistle.md) | Weed–agent system | Biocontrol, plant demography, herbivory |
| 24 | [Asian Citrus Psyllid](tutorials/24_asian_citrus_psyllid.md) | Citrus–psyllid–disease | Vector-borne disease, transmission, host dynamics |
| 25 | [Bemisia tabaci Invasion](tutorials/25_bemisia_tabaci.md) | Whitefly in Europe | Climate suitability, invasion risk, development |
| 26 | [Light Brown Apple Moth](tutorials/26_light_brown_apple_moth.md) | *Epiphyas postvittana* | Invasion potential, climate limits, phenology |
| 27 | [Spotted Alfalfa Aphid](tutorials/27_spotted_alfalfa_aphid.md) | Alfalfa–aphid–natural enemy | Biological control, predator–prey, crop damage |
| 28 | [Olive–Olive Fly Climate](tutorials/28_olive_bactrocera.md) | Olive–*Bactrocera oleae* | Climate warming, pest phenology, crop losses |
| 29 | [Cowpea–Thrips](tutorials/29_cowpea_thrips.md) | Cowpea–thrips in West Africa | Agroecosystem, pest damage, yield loss |
| 30 | [Bean Growth Types](tutorials/30_bean_growth.md) | Common bean types I–III | Growth types, yield prediction, phenology |
| 31 | [Cotton–Boll Weevil](tutorials/31_cotton_boll_weevil.md) | Cotton–weevil in Brazil | Pest pressure, crop dynamics, management |
| 32 | [Oleander Scale Parasitoids](tutorials/32_oleander_scale.md) | Scale–parasitoid system | Competing parasitoids, regulation, biocontrol |
| 33 | [Rice–Fish Agroecosystem](tutorials/33_rice_fish_agroecosystem.md) | Rice–fish integration | Agroecology, coupled production, resource partitioning |
| 34 | [Plant–Aphid–Parasitoid](tutorials/34_plant_aphid_parasite.md) | Tritrophic system | Trophic interactions, host quality, parasitism |
| 35 | [Apple Tree Growth](tutorials/35_apple_tree.md) | Golden Delicious apple | Carbon allocation, phenology, growth |
| 36 | [Fusarium–Nematode in Cotton](tutorials/36_fusarium_nematode.md) | Cotton disease–nematode | Soil-borne disease, nematodes, crop loss |
| 37 | [Whitefly–Autoparasitoid](tutorials/37_whitefly_autoparasitoid.md) | Cotton–whitefly–autoparasitoid | Autoparasitism, population dynamics, biocontrol |
| 38 | [Tropical Fruit Flies](tutorials/38_tropical_fruit_flies.md) | Tropical fruit fly invasion | Climate change, invasion potential, establishment |
| 39 | [China Bt Cotton Economics](tutorials/39_china_bt_cotton.md) | Bt cotton in China | Bioeconomics, pest control, adoption |
| 40 | [Fall Armyworm in Europe](tutorials/40_spodoptera_frugiperda.md) | *Spodoptera frugiperda* | Establishment risk, climate suitability, spread |
| 41 | [Olive Fly ODE](tutorials/41_olive_fly_ode.md) | Olive fruit fly | ODEs, phenology, population dynamics |
| 42 | [Spodoptera exigua Biocontrol](tutorials/42_spodoptera_biocontrol.md) | Beet armyworm–parasitoid | Biological control, host–parasitoid, mortality |
| 43 | [Philaenus Phenology](tutorials/43_philaenus_phenology.md) | *Philaenus spumarius* | Phenology, development rates, climate |
| 44 | [Physiological Risk Index](tutorials/44_risk_index.md) | Pest establishment | Thermal physiology, risk index, establishment |
| 45 | [D. suzukii PDE](tutorials/45_dsuzukii_pde.md) | *Drosophila suzukii* | PDEs, stage structure, survival |
| 46 | [Medfly Kolmogorov](tutorials/46_medfly_kolmogorov.md) | Medfly climate change | Distribution dynamics, dispersal, climate |
| 47 | [Consumer–Resource Dynamics](tutorials/47_consumer_resource.md) | Stage-structured consumer–resource | Stage structure, consumer–resource interactions |
| 48 | [Xylella Eco-Epidemiology](tutorials/48_xylella_ecoepi.md) | Olive grove–*Xylella* | Disease transmission, vectors, olive production |
| 49 | [Bombus HSP Tolerance](tutorials/49_bombus_hsp.md) | *Bombus terrestris* | Genetic modification, heat shock proteins, thermal tolerance |
| 50 | [Type Hierarchy Reference](tutorials/50_type_hierarchy.md) | Core types | Abstract types, inheritance, model structure |
| 51 | [State Variables & Callbacks](tutorials/51_state_variables.md) | State machinery | State variables, bulk populations, callbacks |
| 52 | [Extended Rules & Events](tutorials/52_extended_rules_events.md) | Rules framework | Events, scheduling, model updates |
| 53 | [Theory Helpers](tutorials/53_theory_helpers.md) | Analysis utilities | Compensation, isoclines, assembly rules |
| 54 | [Continuous-Time & PSPM](tutorials/54_continuous_pspm.md) | Continuous models | Continuous time, PSPM, transport equations |
| 55 | [Management & Economics](tutorials/55_management_economics.md) | Optimization | Economics, optimization, control tactics |
| 56 | [Ensembles & Weather](tutorials/56_ensembles_misc.md) | Simulation utilities | Ensembles, filters, weather data |
| 57 | [Verticillium DP Management](tutorials/57_verticillium_dp.md) | Verticillium wilt | Dynamic programming, control, seasonality |
| 58 | [CBB Bioeconomics](tutorials/58_cbb_bioeconomics.md) | Coffee berry borer | Bioeconomics, control tactics, optimization |
| 59 | [Olive–Fly Bioeconomics](tutorials/59_olive_climate.md) | Olive–olive fly | Climate warming, bioeconomics, pest management |
| 60 | [Tuta absoluta Mechanistic PBDM](tutorials/60_tuta_absoluta_invasion.md) | *Tuta absoluta* invasion | Invasion risk, mechanistic PBDM, climate |
| 61 | [Lobesia Voltinism Shifts](tutorials/61_lobesia_voltinism.md) | *Lobesia botrana* | Voltinism, climate warming, phenology |
| 62 | [Screwworm SIT Overflooding](tutorials/62_screwworm_sit.md) | *Cochliomyia* SIT | SIT, overflooding ratio, eradication |
| 63 | [Bayesian Stage Mortality](tutorials/63_bayesian_mortality.md) | Stage-structured inference | Bayesian inference, mortality, stage structure |
| 64 | [BMSB Tritrophic Biocontrol](tutorials/64_bmsb_tritrophic.md) | Tritrophic stink bug | Tritrophic interactions, biocontrol, parasitoids |
