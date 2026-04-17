# VLM-based figure review — design note

## Context

Each PBDM vignette renders figures that should resemble figures in the source
paper. Numerical regressions catch bugs in state trajectories and summary
metrics, but they do not catch visual regressions — e.g. a phenology plot whose
shape has drifted away from the published figure, or a coupled dynamics plot
whose relative peaks no longer align.

## Goal

Build an advisory pipeline that flags likely figure regressions by comparing
rendered vignette figures against the corresponding published-paper figures.

## Why "advisory"

- Published figures are under publisher copyright; we cannot vendor them into
  this repository.
- Even with access, figures differ in rendering (fonts, axis extents, color
  palette, aspect), so an exact pixel match is not a useful criterion.
- The only robust automated signal is *qualitative correspondence*: same
  trajectory shape, same ordering of peaks, same sign of response to a
  perturbation.

So the VLM step should produce a **human-readable judgement**, not a hard
pass/fail gate.

## Recommended pipeline

1. **Figure registry (per vignette).** Add a lightweight YAML/JSON file at
   `vignettes/<NN>_<slug>/figures.yaml` listing:
   - rendered figure filename (relative to the vignette dir)
   - a short natural-language description of what the figure should show
   - reference citation key pointing into `papers/refs.txt`
   - reference figure label (e.g. "Fig. 3a")
   - qualitative acceptance criteria (e.g. "larval peak precedes adult peak",
     "crop-load × yield curve is monotone")

2. **Reference-figure handling.** Reference images are **not** committed.
   Instead, store a local-only path override (e.g. an env var
   `PBDM_REFERENCE_FIGURE_DIR`) that points at a user-provided directory of
   scanned/extracted paper figures. CI runs skip the VLM step when that
   directory is absent.

3. **VLM prompt template.** For each (rendered, reference) pair, send the two
   images plus the description and acceptance criteria to a vision-capable
   model with a prompt such as:

   > You are reviewing a scientific figure regression. Image A is the rendered
   > figure from our simulation. Image B is the reference figure from the
   > original paper. Ignore styling differences (colors, fonts, axis labels).
   > For each acceptance criterion below, answer PASS / FAIL / UNSURE and give
   > a one-sentence justification.

4. **Report artifact.** Emit a Markdown table at
   `vignettes/_reports/figure-review.md` with one row per criterion. A CI job
   can fail only if the VLM reports FAIL on a criterion that was previously
   PASS.

5. **Rendered-only fallback.** When no reference image is available, run a
   **self-consistency** prompt: ask the VLM whether the rendered figure is
   consistent with its own caption and with the vignette's narrative text
   immediately above it. This catches the common regression where a plot
   renders correctly but shows the wrong variable.

## Current limitations

- No VLM is wired into the Julia test process, and Julia's test runner is not a
  good place to call an external VLM API (non-determinism, flakiness, API key
  management).
- A shell-level harness outside `Pkg.test()` is the right layer: a
  `scripts/review_figures.jl` driver that can be invoked manually or from CI,
  reading the figure registry and writing the report.
- Reference figures remain a manual, user-provided input for copyright reasons.

## Minimum viable next step

Pick 3 vignettes with well-known reference figures, write their
`figures.yaml` registry entries by hand, and stand up `scripts/review_figures.jl`
as a prompt-only stub (no API call yet). Once the schema is proven, bolt on an
actual VLM client as a separate change.

## Non-goals

- Pixel-diffing. Style and rendering differences dominate, so pixel diffs are
  not informative.
- Committing reference figures. Copyright on the source papers prevents this.
- Making VLM review a required gate in `Pkg.test()`. It should stay advisory.
