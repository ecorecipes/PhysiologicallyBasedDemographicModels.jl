# Changelog

All notable user-visible changes to `PhysiologicallyBasedDemographicModels.jl`
are recorded here. The format is loosely based on
[Keep a Changelog](https://keepachangelog.com/) and the project follows
semantic versioning once a `1.0` is tagged.

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
