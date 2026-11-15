# Changelog

All notable user-visible changes to `PhysiologicallyBasedDemographicModels.jl`
are recorded here. The format is loosely based on
[Keep a Changelog](https://keepachangelog.com/) and the project follows
semantic versioning once a `1.0` is tagged.

## [0.2.0]

### Changed (breaking): coupled systems are shared with MatrixProjectionModels.jl

The data layer of coupled systems now lives in StructuredPopulationCore
(`src/coupled_systems.jl`), and this package and MatrixProjectionModels.jl both build on
it. The two packages can be loaded together: the names they both export are the same
types and functions. This package still exports all of them, and the solver and its
results are unchanged, but the following differs from 0.1.

- **`PopulationComponent`** has the fields `model`, `population`, `stage_names`,
  `species`, `type`, `patch`. For a `Population` or `BulkPopulation`, `model` is
  `nothing` and `population` is the carrier, held by reference as before. Build
  components with `PopulationComponent(pop; species, type, patch)` or by giving the
  carrier to `PopulationSystem`; the positional
  `PopulationComponent{T}(pop, species, type, patch)` is gone.
- **`PopulationSystem`** has the fields `order`, `components`, `state`. The `names` field
  is gone (use `keys(system)` or `system.order`). `state` is a `Dict{Symbol, Any}`: a
  vector of state variables is stored by reference as before, and a `NamedTuple`,
  dictionary or vector of pairs of plain values is now accepted too (plain values are
  copied in, are not updated automatically, and are recorded in `state_history`).
  `set_state!(system, name, value)` on a name that does not exist now creates a plain
  entry instead of throwing a `KeyError`.
- **`PulseRelease`** starts on day 0 by default (`start_day = 1` before). Pass
  `start_day = 1` to keep the old schedule.
- **`PulseRelease`, `SingleDayRelease`, `SprayEvent`** have a `name` field (first), and
  the amount is the element type of a parametric struct. Constructors check their
  arguments: `interval > 0`, `stage_idx >= 1`, a finite nonnegative amount, kill
  fractions in `[0, 1]`, and as many kill fractions as targets. The positional forms
  `PulseRelease(target, stage_idx, amount, interval; ...)` and
  `SingleDayRelease(target, stage_idx, amount, day)` still work; the six-argument
  `PulseRelease(target, stage_idx, amount, interval, start_day, end_day)` does not (use
  the keywords).
- **`SprayEvent`** with a target that is not a component of the system throws an
  `ArgumentError` when it fires. It was skipped silently before.
- **`TransferRule`, `ReproductionRule`, `MortalityRule`** have a `name` field (first).
  `ReproductionRule`'s fields are `target` and `stage` (were `offspring_target`,
  `offspring_stage`), and `stage` must be at least 1. `MortalityRule` has a `stage` field
  (`0`, the default, means all stages, as before; a positive value limits the mortality to
  that stage). The positional forms `TransferRule(source, source_stage, target,
  target_stage, fn)` and `ReproductionRule(target, stage, fn)` still work. The solver
  logs these three rule types under their type name, as before, whatever their `name`.
- **`remake`** is `SciMLBase.remake` with a method for `PBDMProblem`, not a function of
  this package. SciMLBase is now a direct dependency.
- **`FixedMeshUpwind`** is the type defined in StructuredPopulationCore and shared with
  ContinuousStatePopulationDynamics.jl. `FixedMeshUpwind(n_mesh = 50)` works as before;
  `FixedMeshUpwind()` has `n_mesh == nothing`, which PSPM solves read as 100, the old
  default. `AbstractPSPMMethod` is an alias of the core's
  `AbstractTransportDiscretization`.

### Added

- `apply_event!(event, system, day, p)`, without the weather argument, for the three
  shared events.
- Carrier methods `inject!(pop, stage, amount)` and `remove_fraction!(pop, [stage,]
  fraction)` for `Population`.

## [Unreleased]

### Added

- **`docs/`** — Documenter.jl site that surfaces all 64 vignettes as
  tutorial pages alongside the API reference. The driver
  (`docs/make.jl`) syncs the pre-rendered GFM markdown and figure
  directories from `vignettes/<n>/` into `docs/src/tutorials/`,
  rewrites Quarto's multi-line `<img …>` blocks to standard markdown
  image syntax, and produces an HTML site under `docs/build/`. Set
  `PBDM_BUILD_PDF=true` to additionally produce a LaTeX/PDF build,
  or `PBDM_RERENDER=true` to re-run `quarto render` per vignette.
  Convenience wrapper: `bash scripts/build_docs.sh`.
- **`vignettes/README.md`** — index of all 64 vignettes with title, source
  paper bib-key, and post-audit fidelity grade.
- **`src/temperature_responses.jl`** — eight reusable temperature-response
  helpers extracted from the vignettes during the comprehensive review:
  `triangular_thermal_scalar` / `phiT`, `briere_rate`, `fecundity_briere`,
  `fecundity_gaussian`, `daily_mortality_quadratic`,
  `diapause_fraction_logistic`, `diapause_fraction_linear`,
  `gilbert_fraser_attack`. Covered by
  `test/test_temperature_responses.jl` (7 testsets).
- Multi-format render driver: `scripts/render_all_vignettes.sh` now accepts
  an output-format argument (`gfm`, `html`, or `pdf`; default `gfm`).
  All 64 vignettes render cleanly in all three formats.

### Changed (vignette refactors)

- **#62 screwworm SIT** and **#64 BMSB tritrophic** ported onto the new
  shared temperature-response helpers; behaviour preserved.
- **`scripts/render_all_vignettes.sh`** parameterised on output format
  (writes per-format logs to `/tmp/vign_run_log_${FORMAT}.csv`).

### Fixed (vignette fidelity audit)

A workspace-wide audit (recorded in `vign_audit_01_30.md`,
`vign_audit_31_56.md`, `vign_audit_57_64.md`) compared every vignette
against its source paper. The following fixes were applied to bring
quantitative results back in line with the published figures:

- **#14 East Coast Fever** — corrected two parameter transcription errors:
  `ϑ` 0.004087 → 0.008087, `ρ` 1.8869 → 1.18869. Lifts grade C+ → A.
- **#29 Cowpea thrips** — corrected stage degree-day totals
  (DD_pupa 63 → 91.4, DD_adult 250 → 571) and added explicit
  pedagogical-simplification disclaimer for fecundity / mortality
  forms. C → B.
- **#30 Bean growth** — removed spurious respiration divisor; restored
  type-specific organ lifetimes per Loomis & Ng. C → B.
- **#34 Plant–aphid–parasitoid** — corrected aphid and parasitoid
  thermal thresholds (`T_BASE`).
- **#38 Tropical fruit flies** — added linear density-dependence note
  for the egg→larva transition.
- **#46 Medfly Kolmogorov** — added pedagogical simplification
  disclaimer for the moment-closure derivation.
- **#48 Xylella eco-epidemiology** — added DDE/ODE simplification
  disclaimer.
- **#49 Bombus HSP** — explicitly flagged scenarios as speculative
  rather than calibrated.
- **#57 Verticillium DP** — fixed 10× decimal error in the alternative
  inoculation density (`π_alt` 31.4 → 314.0). C → A.
- **#58 CBB bioeconomics** — corrected Colombia coffee-utility
  coefficient (−0.6228 → −0.7607). B → A.
- **#59 Olive climate** — fixed `\max` / `\min` LaTeX subscripting
  (`Y_\max` → `Y_{\max}`) so the vignette renders to PDF under
  `lualatex`.
- **#60 Tuta absoluta invasion** — added `cap` keyword argument to
  microclimate adjustment routine.
- **#61 Lobesia voltinism** — restored published oviposition `Tm` /
  `c`, diapause-pupa thermal sum (`ΔDP` = 115), and stage-specific
  Briere constants. B → A.
- **#62 Screwworm SIT** — corrected stray `0.5` factor in the
  extinction-equation prose.
- **#63 Bayesian mortality** — re-anchored citations and explicitly
  flagged the log-Gaussian likelihood as a pedagogical approximation
  of the paper's Gamma. C → B.
- **#64 BMSB tritrophic** — raised standing densities to outbreak
  level so the `α` calibration matches the paper's response surface.
  C → B.

### Documentation

- HTML and PDF render outputs are produced under `vignettes/<n>/<n>.html`
  and `vignettes/<n>/<n>.pdf`; both are gitignored. Only the `.qmd`
  source and the GFM `.md` outputs are tracked in git.

## [0.1.0]

Initial worked-corpus release. See `git log` for the full development
history; key milestones are captured in the workspace
`~/.copilot/session-state/.../checkpoints/` index.

[Unreleased]: ./
[0.1.0]: ./
