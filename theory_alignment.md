# Theory–Implementation Correspondence

This document maps the foundational constructs from the PBDM theory literature to their implementations in `PhysiologicallyBasedDemographicModels.jl`.

## Foundational Papers

| # | Citation | Key Contributions |
|---|---------|-------------------|
| 1 | Gutierrez, A.P. (1994). "A Physiologically Based Tritrophic Perspective on Bottom-Up–Top-Down Regulation of Populations" | Mass dynamics dM/dt, metabolic pool model, demand-driven functional response ("hunting equation"), zero-growth isoclines, apparency parameter |
| 2 | Gutierrez, A.P. et al. (1984). "A General Distributed Delay Time Varying Life Table Plant Population Model" | Erlang k-stage distributed delay, time-varying life table, supply/demand ratio, temperature-dependent development rates |
| 3 | Schreiber, S.J. & Gutierrez, A.P. (1998). "A Supply/Demand Perspective of Species Invasions and Coexistence" | Metabolic compensation point (MCP), ecological compensation point (ECP), invasion criteria, food-chain assembly diagrams |
| 4 | Wang, Y.H. & Gutierrez, A.P. (1980). "An Assessment of the Use of Stability Analyses in Population Ecology" | Jacobian eigenvalue analysis, equilibrium-space classification, age-structured stability conditions |
| 5 | Regev, U., Gutierrez, A.P. & Schreiber, S.J. (1998). "Biological and Economic Foundations of Renewable Resource Exploitation" | Bioeconomic optimization, optimal control, discount rate effects, technology parameter |
| 6 | Gutierrez, A.P. & Regev, U. (2005). "The Bioeconomics of Tritrophic Systems: Applications to Invasive Species" | Energy allocation priority scheme (wastes → respiration → reproduction → growth), r–K strategy analysis, constrained optimization |

---

## ✅ Strong Correspondence

Core constructs that are faithfully implemented.

### Distributed Delay (Erlang k-Stage)

- **Paper**: Gutierrez et al. (1984), Eq. 3–5 — the Manetsch/Vansickle k-substage distributed delay model with physiological-time stepping and attrition.
- **Implementation**: `DistributedDelay(k, τ; W0)` + `step_delay!(delay, rate, dt; attrition)`
- **Assessment**: Direct and faithful. The delay struct stores k substages, mean duration τ, and mass W0. Stepping advances cohorts through substages at a rate driven by physiological time (degree-days). Attrition is applied per-substage.

### Demand-Driven Functional Response

- **Paper**: Gutierrez (1994), Eq. 5 — the "hunting equation": acquisition = D(1 − exp(−aS/D)), where D = demand, S = supply, a = apparency.
- **Implementation**: `FraserGilbertResponse(a)` + `acquire(fr, supply, demand)`
- **Assessment**: Exact equation. The name "Fraser-Gilbert" is the original attribution that Gutierrez uses. Returns the fraction of demand satisfied given available supply.

### Metabolic Pool Allocation

- **Paper**: Gutierrez (1994), Gutierrez & Wang (1977) — production is allocated by priority: wastes/respiration first, then reproduction, then growth (or the reverse for r-selected species).
- **Implementation**: `MetabolicPool` + `allocate(pool)` — priority-ordered cascade allocation.
- **Assessment**: Faithful. The priority ordering is user-specified, matching the theory's flexibility to model both r-selected (growth-first) and K-selected (reproduction-first) strategies. Each allocation step draws from a shrinking residual budget.

### Supply/Demand Ratio (φ)

- **Paper**: Central integrating variable across all PBDM papers. φ = acquisition/demand; drives growth, reproduction, stress, and movement.
- **Implementation**: `supply_demand_ratio(fr, S, D)` returns φ ∈ [0, 1]; also `supply_demand_index()` for mapping φ to demographic multipliers.
- **Assessment**: Pervasive in the API, as in the theory. φ appears in every vignette and drives the same cascade of demographic consequences described in the papers.

### Temperature-Driven Development

- **Paper**: Gutierrez et al. (1984), Section 2 — physiological time via degree-day accumulation, with nonlinear corrections at temperature extremes.
- **Implementation**: Three development-rate models:
  - `LinearDevelopmentRate(T_base, T_opt)` — simplest degree-day model
  - `BriereDevelopmentRate(a, T_L, T_M)` — nonlinear insect development curve
  - `LoganDevelopmentRate(ψ, ρ, T_max, ΔT)` — Type III sigmoid with high-temperature decline
- **Assessment**: Covers the spectrum of forms used across the applied papers, from linear (plants) to nonlinear (insects).

### Q₁₀ Respiration

- **Paper**: Gutierrez (1994), Section 3.1 — exponential temperature scaling of metabolic cost: R(T) = R_ref · Q₁₀^((T − T_ref)/10).
- **Implementation**: `Q10Respiration(R_ref, Q10, T_ref)` + `respiration_rate(resp, T)`
- **Assessment**: Exact Q₁₀ formulation from the papers.

### Biodemographic Functions (BDF)

- **Paper**: Gutierrez (1992, 1994) — the composite of development, acquisition, and respiration that defines a species' physiology.
- **Implementation**: `BiodemographicFunctions(dev, acq, resp)` — bundles a development rate model, a functional response, and a respiration model.
- **Assessment**: Clean mapping. BDF objects serve as the physiological "identity card" for each population, exactly as described in the theory.

### Tritrophic Structure

- **Paper**: Gutierrez (1994), Fig. 1 — the canonical plant → herbivore → natural enemy food chain, each level with its own BDF.
- **Implementation**: `TrophicLink` + `TrophicWeb` + `PredationRule` — coupled solver handles N-trophic interactions.
- **Assessment**: Generalizes beyond tritrophic to arbitrary food webs. Each link specifies consumer, resource, and functional response.

### Sterile Insect Technique

- **Paper**: Gutierrez et al. (various applied papers on screwworm and medfly).
- **Implementation**: `SITRelease` + `fertile_mating_fraction(sir, N_wild)` — competitiveness-weighted fertile mating fraction.
- **Assessment**: Faithful to the Knipling-style SIT model used in the PBDM applications.

### Genetics and Resistance Evolution

- **Paper**: Gutierrez et al. (Bt cotton papers) — Hardy-Weinberg selection dynamics for insecticide resistance.
- **Implementation**: `DialleleicLocus`, `GenotypeFitness`, `selection_step!`, `TwoLocusResistance`.
- **Assessment**: Covers single- and two-locus resistance genetics as used in the Bt cotton analyses.

### Bioeconomic Analysis

- **Paper**: Regev, Gutierrez & Schreiber (1998); Gutierrez & Regev (2005) — yield loss, profit analysis, discounting.
- **Implementation**: Economics layer: `CropRevenue`, damage functions, `npv()`, `benefit_cost_ratio()`.
- **Assessment**: Covers the static bioeconomic analysis (costs, revenue, NPV) used in the applied papers.

---

## ⚠️ Partial Correspondence

Present in the package but could be strengthened or made more explicit.

### Metapopulation Dispersal

- **Theory**: Supply/demand-driven emigration, where individuals leave patches where φ < threshold (Gutierrez et al. 1999 cassava; various other papers).
- **Current state**: `SpatialGrid` + `DispersalRule` types exist but the coupled solver doesn't integrate them into the main `solve()` path.
- **Gap**: Spatial models currently require hand-rolled loops (see `pbdm/inprogress/40_cassava_metapopulation/`).
- **Recommendation**: Wire `SpatialGrid` dispersal into the coupled solver as an optional post-interaction phase.

### Energy Currency

- **Theory**: All PBDM papers consistently use grams dry weight as the mass currency, with explicit conversions (e.g., 1 g dry ≈ 7000 calories).
- **Current state**: The package is unit-agnostic — delays track mass but units are not enforced.
- **Gap**: Could add documentation conventions or optional unit annotations.

---

## ❌ Not Yet Implemented

Theory constructs that don't have a corresponding implementation.

*(None — all foundational constructs now have implementations.)*

---

## Recently Implemented (moved from ⚠️/❌ → ✅)

### Continuous ODE/DDE/PSPM Formulations — was ❌, now ✅

- **Paper**: Gutierrez (1994), Eq. 2a — the canonical tritrophic ODE: dM₂/dt = h[f(M₁/M₂)M₂] − g(M₂/M₃)M₃.
- **Implementation**: Three continuous-time formulations via package extensions:
  1. **ODE (Linear Chain Trick)** — `ContinuousPBDMProblem` + `solve_continuous()` via `OrdinaryDiffEqExt`. Each Erlang k-substage becomes an ODE state variable; adaptive Tsit5 time-stepping.
  2. **DDE (Delay Differential Equations)** — `DelayPBDMProblem` + `solve_delay()` via `DelayDiffEqExt`. One state per species per stage with history-based maturation delays solved by `MethodOfSteps(Tsit5())`.
  3. **PSPM (McKendrick–von Foerster PDE)** — `PSPMProblem` + `solve_pspm()` via `OrdinaryDiffEqExt`. Size-structured PDE discretized using FMU (Fixed Mesh Upwind) or EBT (Escalator Boxcar Train) methods, following Joshi et al. (2023).
- **Assessment**: Full continuous-time capability. All three formulations support multi-trophic coupling, temperature-dependent development, Q10 respiration, and Fraser-Gilbert/Holling functional responses. Compatible with the SciML ecosystem for sensitivity analysis, adjoint methods, and ensemble simulations.

The following constructs were added to `src/theory.jl` and `src/optimal_control.jl`:

### Apparency Parameter (*a*) — was ⚠️, now ✅

- **Implementation**: `apparency(fr)` returns the apparency parameter; `is_ratio_dependent(fr)` trait distinguishes ratio-dependent (Fraser-Gilbert) from prey-dependent (Holling) functional responses.
- **Assessment**: Apparency is now an explicit, named API concept. Docstring on `FraserGilbertResponse` cross-references Gutierrez (1994) Eq. 5.

### Metabolic Compensation Point (MCP) — was ⚠️, now ✅

- **Implementation**: `compensation_point(resp, T, demand_rate; conversion_efficiency)` → φ_MCP.
- **Assessment**: Computes φ* = R(T)/(ε·d) exactly as in Schreiber & Gutierrez (1998). Central to coexistence predictions and assembly analysis.

### Zero-Growth Isoclines — was ⚠️, now ✅

- **Implementation**: `consumer_isocline(fr, φ_star, d)` and `resource_isocline(fr, d, φ_star, r, K; M1_range, n_points)` → `IsoclineResult`.
- **Assessment**: Consumer isocline is the analytical straight line M₂ = s·M₁; resource isocline is computed via bisection, yielding the classic hump shape from Gutierrez (1994) Fig. 3–5.

### Ratio-Dependent vs. Prey-Dependent — was ⚠️, now ✅

- **Implementation**: `is_ratio_dependent(fr)` trait returns `true` for `FraserGilbertResponse`, `false` for `HollingTypeII`/`HollingTypeIII`.
- **Assessment**: The API now explicitly flags which paradigm each functional response implements.

### Stability/Eigenvalue Analysis — was ❌, now ✅

- **Implementation**: `find_equilibrium(fr, resp, T; ...)` → `EquilibriumResult` with M₁*, M₂*, Jacobian, eigenvalues, and classification (stable node/focus, unstable node/focus, saddle, center, degenerate). `classify_equilibrium(eigenvalues)` for standalone eigenvalue classification.
- **Assessment**: Full analytical equilibrium computation following Wang & Gutierrez (1980). Jacobian entries are closed-form at the interior equilibrium.

### Assembly Diagrams — was ❌, now ✅

- **Implementation**: `SpeciesProfile(name, fr, resp, T, demand_rate, ε)` + `food_web_assembly(species)` → `AssemblyResult` with invasion order, coexistence species, and MCP values.
- **Assessment**: MCP-ordered species invasion following Schreiber & Gutierrez (1998). Species persist only if φ* < 1; assembly order is MCP-ascending.

### r–K Strategy Parameterization — was ❌, now ✅

- **Implementation**: `life_history_strategy(pool; threshold)` → `:r_selected` | `:K_selected` | `:intermediate`.
- **Assessment**: Classifies based on allocation priority ordering exactly as in Gutierrez & Regev (2005): growth-first = r-selected, reproduction-first = K-selected.

### Optimal Control / Bioeconomic Optimization — was ❌, now ✅

- **Implementation**: `PBDMControlProblem` bundles `TrophicLevel` parameterizations, `AbstractManagementAction` subtypes (`PesticideControl`, `BiologicalReleaseControl`, `HarvestControl`), and `AbstractManagementObjective` subtypes (`MinimizeDamage`, `MaximizeProfit`). The `optimize_management(prob; solver)` function is provided by the JuMP package extension (`ext/JuMPExt.jl`), which Euler-discretizes the PBDM dynamics following the EpiPolicies pattern (Frost & Montes 2025).
- **Assessment**: Full discrete-time optimal control with Ipopt as default NLP solver. Supports quadratic control costs, discounted profit maximization, and arbitrary trophic chain length.

---

## Summary

| Category | Count | Assessment |
|----------|-------|------------|
| ✅ Strong correspondence | 20 constructs | Core simulation + analytical theory tools + continuous formulations |
| ⚠️ Partial correspondence | 0 constructs | — |
| ❌ Not yet implemented | 0 constructs | — |

The package provides both a **simulation framework** and **analytical theory tools** that together faithfully implement the PBDM literature. The core simulation machinery (distributed delays, demand-driven functional responses, metabolic pool allocation, supply/demand ratios, temperature-driven physiology, multi-trophic coupling) is complemented by analytical utilities (compensation points, isocline computation, equilibrium/stability analysis, food-web assembly, life-history classification, optimal control) and continuous-time formulations (ODE via linear chain trick, DDE via explicit delays, PSPM via McKendrick–von Foerster PDE discretization). The 39 vignettes comprehensively cover the applied PBDM literature.

Future extensions: **metapopulation solver integration** (spatial models currently require hand-rolled loops) and **proper Characteristic Method** for PSPM (currently delegates to EBT). These represent natural next steps beyond the foundational theory.
