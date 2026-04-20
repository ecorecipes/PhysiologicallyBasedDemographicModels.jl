# PhysiologicallyBasedDemographicModels.jl

A Julia package for **Physiologically Based Demographic Models** (PBDMs):
physiologically mechanistic, weather-driven, stage-structured population models
of plants, arthropods, and their interactions, following the framework
popularised by Gutierrez, Ponti, and collaborators.

This package is the **application and reference suite** for the projection-model
workspace. It is the exactness oracle that validates the shared abstractions in
[StructuredPopulationCore.jl](https://github.com/ecorecipes/StructuredPopulationCore.jl)
and the continuous-time backends in
[FiniteStatePopulationDynamics.jl](https://github.com/ecorecipes/FiniteStatePopulationDynamics.jl)
and
[ContinuousStatePopulationDynamics.jl](https://github.com/ecorecipes/ContinuousStatePopulationDynamics.jl)
by reproducing published PBDM studies end-to-end.

## Features

### Core abstractions

- `BulkPopulation`, `Population`, `PopulationComponent`, `PopulationSystem` —
  scalar-pool and stage-structured population carriers with optional scalar
  auxiliary state variables (`ScalarState`, `StateVariable`).
- Distributed-delay stage structure with temperature-driven development, per-
  stage mortality, and mass/biomass tracking.
- `MetabolicPool` allocation for supply/demand energy budgeting.

### Coupled multi-species / multi-rule API

- `PBDMProblem` + `DirectIteration` solver with explicit phases:
  weather → stress → step → reproduction → events → callbacks.
- `ReproductionRule`, `MortalityRule`, `StressRule`, `CustomRule`,
  `PhaseCallback`, and weather-conditional events.
- `MultiSpeciesPBDMNew` orchestrates trophic webs with interaction rules between
  populations.
- `run_scenarios` + `compare_metrics` for factorial scenario sweeps.

### Extensions

- **Genome**: single- and multi-locus genotype tracking across populations.
- **Diapause**: latitude- and photoperiod-driven diapause entry/exit.
- **Phenology**: degree-day stage progression with user-defined thresholds.
- **Soil coupling**: soil moisture/temperature state feeding crop stress.
- **Multispecies**: predator–prey, plant–herbivore, tritrophic interactions.

### Vignettes

Sixty-four vignettes cover the full PBDM corpus, ranging from introductory
tutorials through classical reference models (cotton, coffee berry borer,
grapevine, olive fly, screwworm SIT) to recent applications (Bt resistance,
Spodoptera frugiperda, Drosophila suzukii PDE, Xylella eco-epidemiology). All
handrolled simulations have been migrated onto the coupled API and are rendered
via Quarto to HTML, PDF, and GitHub-flavoured Markdown.

## Installation

This package is not yet registered in the Julia General registry. Install
directly from GitHub:

```julia
using Pkg
Pkg.add(url="https://github.com/ecorecipes/PhysiologicallyBasedDemographicModels.jl")
```

## Quick start

```julia
using PhysiologicallyBasedDemographicModels

# Build a single-species stage-structured population.
pop = make_cotton_plant()             # example builder from the vignettes
sys = PopulationSystem(:cotton => PopulationComponent(pop))

prob = PBDMProblem(
    SingleSpeciesPBDM(),
    sys,
    (1, 200);
    weather = make_daily_weather(),
)
sol = solve(prob, DirectIteration())
```

See the [`vignettes/`](vignettes/) directory for fully worked examples.
The annotated index — title, source paper, and post-audit fidelity grade
for all 64 vignettes — is in [`vignettes/README.md`](vignettes/README.md).
A browseable [Documenter.jl site](docs/) (HTML + KaTeX, with all 64
vignettes as tutorial pages and an API reference) can be built locally
with `bash scripts/build_docs.sh`; output lands in `docs/build/`.
Recent fixes and refactors are listed in [`CHANGELOG.md`](CHANGELOG.md).

## Related

- [StructuredPopulationCore.jl](https://github.com/ecorecipes/StructuredPopulationCore.jl)
  — shared abstractions used by every package in this stack
- [MatrixProjectionModels.jl](https://github.com/ecorecipes/MatrixProjectionModels.jl)
  — discrete-stage, discrete-time matrix projection models
- [IntegralProjectionModels.jl](https://github.com/ecorecipes/IntegralProjectionModels.jl)
  — continuous-state, discrete-time integral projection models
- [FiniteStatePopulationDynamics.jl](https://github.com/ecorecipes/FiniteStatePopulationDynamics.jl)
  — finite-state continuous-time dynamics
- [ContinuousStatePopulationDynamics.jl](https://github.com/ecorecipes/ContinuousStatePopulationDynamics.jl)
  — continuous-state continuous-time dynamics
- [CategoricalPopulationDynamics.jl](https://github.com/ecorecipes/CategoricalPopulationDynamics.jl)
  — compositional categorical front-end
