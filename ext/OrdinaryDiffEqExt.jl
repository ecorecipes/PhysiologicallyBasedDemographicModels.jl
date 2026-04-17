"""
OrdinaryDiffEq extension for PhysiologicallyBasedDemographicModels.jl

Provides `solve` methods for:
1. `ContinuousPBDMProblem` — aggregate-biomass ODEs (linear chain trick)
2. `PSPMProblem` — McKendrick–von Foerster PDE via method of lines

Activated automatically when `using OrdinaryDiffEq` is called alongside
the main package.
"""
module OrdinaryDiffEqExt

using PhysiologicallyBasedDemographicModels
using OrdinaryDiffEq

import PhysiologicallyBasedDemographicModels:
    solve_continuous, solve_pspm, staged_species_stage_totals,
    ContinuousPBDMProblem, ContinuousPBDMSolution,
    PSPMProblem, PSPMSpecies, AbstractPSPMSpecies,
    PSPMStage, StagedPSPMSpecies, n_pspm_stages,
    FixedMeshUpwind, EscalatorBoxcarTrain, CharacteristicMethod,
    ImplicitFixedMeshUpwind, ImplicitEscalatorBoxcarTrain, ImplicitCharacteristicMethod,
    LaxFriedrichsUpwind, ImplicitLaxFriedrichsUpwind,
    ContinuousSpecies, ContinuousTrophicLink,
    _total_substages, _get_temperature, species_total_ranges,
    FraserGilbertResponse, HollingTypeII, HollingTypeIII,
    development_rate, respiration_rate, acquire

# ============================================================================
# 1. ContinuousPBDMProblem → ODEProblem (Linear Chain Trick)
# ============================================================================

"""
    _build_ode_rhs(prob::ContinuousPBDMProblem)

Build the ODE right-hand-side function for the linear chain trick formulation.

State vector layout: for each species, for each stage j, k_j substage masses
W_{j,1}, ..., W_{j,k_j} are concatenated.

The ODE for substage ℓ of stage j of species i is:

  dW/dt = r_j · g_j(T) · W_{prev} - (r_j · g_j(T) + μ_j · g_j(T) + σ_i) · W

where:
- r_j = k_j / τ_j (Erlang flow rate)
- g_j(T) = development_rate(dev_j, T(t)) (degree-day rate)
- μ_j = background mortality per degree-day
- σ_i = trophic stress (consumption by predators)
- W_{prev} is the previous substage (or inflow for ℓ=1)

For basal resources (species with intrinsic_rate > 0), the first substage
additionally receives logistic growth input.
"""
function _build_ode_rhs(prob::ContinuousPBDMProblem)
    species = prob.species
    links = prob.links
    T_forcing = prob.T_forcing
    n_sp = length(species)

    # Pre-compute index ranges for each species and each stage
    # species_offsets[i] = starting index for species i
    # stage_ranges[i][j] = UnitRange for substages of stage j of species i
    species_offsets = Int[]
    stage_ranges = Vector{Vector{UnitRange{Int}}}()
    idx = 1
    for sp in species
        push!(species_offsets, idx)
        sr = UnitRange{Int}[]
        for j in 1:sp.n_stages
            push!(sr, idx:(idx + sp.k[j] - 1))
            idx += sp.k[j]
        end
        push!(stage_ranges, sr)
    end

    # Build name → index lookup
    name_to_idx = Dict{Symbol, Int}(sp.name => i for (i, sp) in enumerate(species))

    # Build consumer/resource link lookup
    # consumers_of[i] = [(consumer_idx, consumer_species), ...]
    consumers_of = Dict{Int, Vector{Int}}()
    resources_of = Dict{Int, Vector{Int}}()
    for link in links
        ci = get(name_to_idx, link.consumer, 0)
        ri = get(name_to_idx, link.resource, 0)
        ci > 0 && ri > 0 || continue
        push!(get!(consumers_of, ri, Int[]), ci)
        push!(get!(resources_of, ci, Int[]), ri)
    end

    function rhs!(du, u, p, t)
        T_c = _get_temperature(T_forcing, t)

        # Phase 1: Compute total biomass per species
        totals = zeros(n_sp)
        for i in 1:n_sp
            offset = species_offsets[i]
            n_sub = _total_substages(species[i])
            for idx in offset:(offset + n_sub - 1)
                totals[i] += u[idx]
            end
        end

        # Phase 2: Compute trophic acquisition and consumption
        # acquisition[i] = total biomass acquired by species i from resources
        # consumed_from[i] = total biomass consumed from species i by predators
        acquisition = zeros(n_sp)
        consumed_from = zeros(n_sp)

        for i in 1:n_sp
            sp = species[i]
            if haskey(resources_of, i) && !isempty(resources_of[i])
                # Consumer: compute demand-driven acquisition from each resource
                demand = sp.demand_rate * totals[i]
                for ri in resources_of[i]
                    supply = totals[ri]
                    acq = _ode_acquire(sp.fr, supply, demand)
                    acquisition[i] += sp.conversion_efficiency * acq
                    consumed_from[ri] += acq
                end
            end
        end

        # Phase 3: Compute du/dt for each substage
        for i in 1:n_sp
            sp = species[i]
            total_i = totals[i]

            # Per-capita consumption rate from predators
            predation_rate = total_i > 0 ? consumed_from[i] / total_i : 0.0

            # Per-capita acquisition rate (for consumers)
            acq_rate = total_i > 0 ? acquisition[i] / total_i : 0.0

            for j in 1:sp.n_stages
                r = stage_ranges[i][j]
                k_j = sp.k[j]
                r_j = k_j / sp.τ[j]  # Erlang flow rate
                g_j = development_rate(sp.dev_rate[j], T_c)  # DD/day at current T
                μ_j = sp.μ[j]

                # Respiration rate (per-capita per day)
                R_i = respiration_rate(sp.resp, T_c)

                for ℓ in 1:k_j
                    idx = r[ℓ]
                    W = u[idx]

                    # Developmental flow
                    flow_rate = r_j * g_j
                    mortality = μ_j * g_j

                    # Inflow from previous substage
                    if ℓ == 1 && j == 1
                        # First substage of first stage
                        if sp.intrinsic_rate > 0
                            # Basal resource: logistic growth as inflow
                            inflow = sp.intrinsic_rate * total_i * (1.0 - total_i / sp.carrying_capacity)
                            inflow = max(0.0, inflow)
                        else
                            # Consumer: acquisition feeds into first substage
                            inflow = acq_rate * total_i
                        end
                    elseif ℓ == 1
                        # First substage of stage j: receives outflow from last substage of stage j-1
                        prev_r = stage_ranges[i][j-1]
                        prev_last = prev_r[end]
                        prev_flow_rate = sp.k[j-1] / sp.τ[j-1]
                        prev_g = development_rate(sp.dev_rate[j-1], T_c)
                        inflow = prev_flow_rate * prev_g * u[prev_last]
                    else
                        # Internal substage: receives from previous substage
                        inflow = flow_rate * u[idx - 1]
                    end

                    # Outflow + mortality + predation + respiration
                    du[idx] = inflow - flow_rate * W - mortality * W - predation_rate * W - R_i * W
                end
            end
        end
    end

    return rhs!
end

"""
Acquire biomass (ODE-safe, no branching on zero for AD compatibility).
"""
function _ode_acquire(fr::FraserGilbertResponse, supply, demand)
    (demand <= 0 || supply <= 0) && return 0.0
    return demand * (1.0 - exp(-fr.a * supply / demand))
end

function _ode_acquire(fr::HollingTypeII, supply, demand)
    supply <= 0 && return 0.0
    return fr.a * supply / (1.0 + fr.a * fr.h * supply) * (demand > 0 ? demand / (demand + 1e-10) : 0.0)
end

function _ode_acquire(fr, supply, demand)
    return acquire(fr, supply, demand)
end

"""
    PhysiologicallyBasedDemographicModels.solve_continuous(
        prob::ContinuousPBDMProblem; alg=nothing, kwargs...)

Solve an aggregate-biomass ODE PBDM problem.

# Arguments
- `prob`: The continuous PBDM problem
- `alg`: ODE solver algorithm (default: auto-selected by OrdinaryDiffEq)
- `kwargs...`: Passed to OrdinaryDiffEq.solve (e.g., `abstol`, `reltol`, `saveat`)

# Returns
`ContinuousPBDMSolution` with time series, species trajectories, and the raw SciML solution.
"""
function PhysiologicallyBasedDemographicModels.solve_continuous(
        prob::ContinuousPBDMProblem;
        alg=nothing, kwargs...)

    rhs! = _build_ode_rhs(prob)

    ode_prob = ODEProblem(rhs!, prob.u0, prob.tspan, prob.p)

    if alg === nothing
        sol = OrdinaryDiffEq.solve(ode_prob; kwargs...)
    else
        sol = OrdinaryDiffEq.solve(ode_prob, alg; kwargs...)
    end

    # Convert to ContinuousPBDMSolution
    t_vec = Float64.(sol.t)
    n_states = length(prob.u0)
    n_t = length(t_vec)

    u_mat = zeros(n_states, n_t)
    for (j, ui) in enumerate(sol.u)
        u_mat[:, j] = ui
    end

    sp_ranges = species_total_ranges(prob)
    sp_names = [sp.name for sp in prob.species]
    retcode = Symbol(string(sol.retcode))

    return ContinuousPBDMSolution{Float64, typeof(sol)}(
        t_vec, u_mat, sp_names, sp_ranges, retcode, sol)
end

# ============================================================================
# 2. PSPMProblem → ODEProblem (McKendrick–von Foerster PDE discretization)
# ============================================================================

"""
    PhysiologicallyBasedDemographicModels.solve_pspm(
        prob::PSPMProblem; alg=nothing, kwargs...)

Solve a size-structured population model by discretizing the
McKendrick–von Foerster PDE and solving the resulting ODE system.

Dispatches on `prob.method` — eight methods following Joshi et al. (2023):

**Explicit methods** (default `Tsit5()` solver):
- `FixedMeshUpwind`: first-order upwind finite differences on a fixed grid
- `EscalatorBoxcarTrain`: cohort-tracking method (de Roos 1988)
- `CharacteristicMethod`: boundary-tracking along characteristics
- `LaxFriedrichsUpwind`: Lax-Friedrichs flux for extra diffusive stability

**Semi-implicit methods** (default `Rosenbrock23()` stiff solver):
- `ImplicitFixedMeshUpwind`: IFMU with implicit mortality
- `ImplicitEscalatorBoxcarTrain`: IEBT with implicit mortality
- `ImplicitCharacteristicMethod`: ICM with implicit mortality
- `ImplicitLaxFriedrichsUpwind`: IFMU-LF with implicit mortality

# Returns
`ContinuousPBDMSolution` with size-density state trajectories.
"""
function PhysiologicallyBasedDemographicModels.solve_pspm(
        prob::PSPMProblem; alg=nothing, kwargs...)
    return _solve_pspm(prob, prob.method; alg=alg, kwargs...)
end

# ============================================================================
# Helper: wrap ODE solution into ContinuousPBDMSolution
# ============================================================================

function _wrap_pspm_solution(sol, u0, prob, sp_ranges)
    t_vec = Float64.(sol.t)
    n_states = length(u0)
    n_t = length(t_vec)
    u_mat = zeros(n_states, n_t)
    for (j, ui) in enumerate(sol.u)
        u_mat[:, j] = ui
    end
    sp_names = [sp.name for sp in prob.species]
    retcode = Symbol(string(sol.retcode))
    return ContinuousPBDMSolution{Float64, typeof(sol)}(
        t_vec, u_mat, sp_names, sp_ranges, retcode, sol)
end

# ============================================================================
# Helper: solve an ODEProblem with default algorithm selection
# ============================================================================

function _solve_ode(ode_prob, alg, stiff::Bool; kwargs...)
    if alg !== nothing
        return OrdinaryDiffEq.solve(ode_prob, alg; kwargs...)
    elseif stiff
        return OrdinaryDiffEq.solve(ode_prob, Rosenbrock23(); kwargs...)
    else
        return OrdinaryDiffEq.solve(ode_prob; kwargs...)
    end
end

# ============================================================================
# 2a. Fixed Mesh Upwind (FMU) — explicit
# ============================================================================

function _solve_pspm(prob::PSPMProblem, method::FixedMeshUpwind;
                     alg=nothing, kwargs...)
    n_mesh = method.n_mesh
    if _has_staged(prob)
        meshes, u0, sp_ranges = _init_staged_fmu_state(prob, n_mesh)
        rhs! = _build_staged_fmu_rhs(prob, meshes, n_mesh, :upwind)
    else
        meshes, u0, sp_ranges = _init_fmu_state(prob, n_mesh)
        rhs! = _build_fmu_rhs(prob, meshes, n_mesh, :upwind)
    end
    ode_prob = ODEProblem(rhs!, u0, prob.tspan, prob.p)
    sol = _solve_ode(ode_prob, alg, false; kwargs...)
    return _wrap_pspm_solution(sol, u0, prob, sp_ranges)
end

# ============================================================================
# 2b. Implicit Fixed Mesh Upwind (IFMU) — semi-implicit / stiff
# ============================================================================

function _solve_pspm(prob::PSPMProblem, method::ImplicitFixedMeshUpwind;
                     alg=nothing, kwargs...)
    n_mesh = method.n_mesh
    if _has_staged(prob)
        meshes, u0, sp_ranges = _init_staged_fmu_state(prob, n_mesh)
        rhs! = _build_staged_fmu_rhs(prob, meshes, n_mesh, :upwind)
    else
        meshes, u0, sp_ranges = _init_fmu_state(prob, n_mesh)
        rhs! = _build_fmu_rhs(prob, meshes, n_mesh, :upwind)
    end
    ode_prob = ODEProblem(rhs!, u0, prob.tspan, prob.p)
    sol = _solve_ode(ode_prob, alg, true; kwargs...)
    return _wrap_pspm_solution(sol, u0, prob, sp_ranges)
end

# ============================================================================
# 2c. Lax-Friedrichs Upwind (FMU-LF) — explicit with numerical diffusion
# ============================================================================

function _solve_pspm(prob::PSPMProblem, method::LaxFriedrichsUpwind;
                     alg=nothing, kwargs...)
    n_mesh = method.n_mesh
    if _has_staged(prob)
        meshes, u0, sp_ranges = _init_staged_fmu_state(prob, n_mesh)
        rhs! = _build_staged_fmu_rhs(prob, meshes, n_mesh, :lax_friedrichs)
    else
        meshes, u0, sp_ranges = _init_fmu_state(prob, n_mesh)
        rhs! = _build_fmu_rhs(prob, meshes, n_mesh, :lax_friedrichs)
    end
    ode_prob = ODEProblem(rhs!, u0, prob.tspan, prob.p)
    sol = _solve_ode(ode_prob, alg, false; kwargs...)
    return _wrap_pspm_solution(sol, u0, prob, sp_ranges)
end

# ============================================================================
# 2d. Implicit Lax-Friedrichs Upwind (IFMU-LF) — semi-implicit LF
# ============================================================================

function _solve_pspm(prob::PSPMProblem, method::ImplicitLaxFriedrichsUpwind;
                     alg=nothing, kwargs...)
    n_mesh = method.n_mesh
    if _has_staged(prob)
        meshes, u0, sp_ranges = _init_staged_fmu_state(prob, n_mesh)
        rhs! = _build_staged_fmu_rhs(prob, meshes, n_mesh, :lax_friedrichs)
    else
        meshes, u0, sp_ranges = _init_fmu_state(prob, n_mesh)
        rhs! = _build_fmu_rhs(prob, meshes, n_mesh, :lax_friedrichs)
    end
    ode_prob = ODEProblem(rhs!, u0, prob.tspan, prob.p)
    sol = _solve_ode(ode_prob, alg, true; kwargs...)
    return _wrap_pspm_solution(sol, u0, prob, sp_ranges)
end

# ============================================================================
# 2e. Escalator Boxcar Train (EBT) — explicit
# ============================================================================

function _solve_pspm(prob::PSPMProblem, method::EscalatorBoxcarTrain;
                     alg=nothing, kwargs...)
    max_c = method.max_cohorts
    if _has_staged(prob)
        u0, sp_ranges = _init_staged_ebt_state(prob, max_c)
        rhs! = _build_staged_ebt_rhs(prob, max_c)
    else
        u0, sp_ranges = _init_ebt_state(prob, max_c)
        rhs! = _build_ebt_rhs(prob, max_c)
    end
    ode_prob = ODEProblem(rhs!, u0, prob.tspan, prob.p)
    sol = _solve_ode(ode_prob, alg, false; kwargs...)
    return _wrap_pspm_solution(sol, u0, prob, sp_ranges)
end

# ============================================================================
# 2f. Implicit EBT (IEBT) — semi-implicit / stiff
# ============================================================================

function _solve_pspm(prob::PSPMProblem, method::ImplicitEscalatorBoxcarTrain;
                     alg=nothing, kwargs...)
    max_c = method.max_cohorts
    if _has_staged(prob)
        u0, sp_ranges = _init_staged_ebt_state(prob, max_c)
        rhs! = _build_staged_ebt_rhs(prob, max_c)
    else
        u0, sp_ranges = _init_ebt_state(prob, max_c)
        rhs! = _build_ebt_rhs(prob, max_c)
    end
    ode_prob = ODEProblem(rhs!, u0, prob.tspan, prob.p)
    sol = _solve_ode(ode_prob, alg, true; kwargs...)
    return _wrap_pspm_solution(sol, u0, prob, sp_ranges)
end

# ============================================================================
# 2g. Characteristic Method (CM) — explicit, boundary-tracking
# ============================================================================

function _solve_pspm(prob::PSPMProblem, method::CharacteristicMethod;
                     alg=nothing, kwargs...)
    max_c = method.max_cohorts
    if _has_staged(prob)
        u0, sp_ranges = _init_staged_cm_state(prob, max_c)
        rhs! = _build_staged_cm_rhs(prob, max_c)
    else
        u0, sp_ranges = _init_cm_state(prob, max_c)
        rhs! = _build_cm_rhs(prob, max_c)
    end
    ode_prob = ODEProblem(rhs!, u0, prob.tspan, prob.p)
    sol = _solve_ode(ode_prob, alg, false; kwargs...)
    return _wrap_pspm_solution(sol, u0, prob, sp_ranges)
end

# ============================================================================
# 2h. Implicit Characteristic Method (ICM) — semi-implicit
# ============================================================================

function _solve_pspm(prob::PSPMProblem, method::ImplicitCharacteristicMethod;
                     alg=nothing, kwargs...)
    max_c = method.max_cohorts
    if _has_staged(prob)
        u0, sp_ranges = _init_staged_cm_state(prob, max_c)
        rhs! = _build_staged_cm_rhs(prob, max_c)
    else
        u0, sp_ranges = _init_cm_state(prob, max_c)
        rhs! = _build_cm_rhs(prob, max_c)
    end
    ode_prob = ODEProblem(rhs!, u0, prob.tspan, prob.p)
    sol = _solve_ode(ode_prob, alg, true; kwargs...)
    return _wrap_pspm_solution(sol, u0, prob, sp_ranges)
end

# ============================================================================
# Shared initialisation: FMU-family (FMU, IFMU, FMU-LF, IFMU-LF)
# ============================================================================

function _init_fmu_state(prob, n_mesh)
    meshes = []
    u0 = Float64[]
    sp_ranges = Dict{Symbol, UnitRange{Int}}()
    idx = 1
    for sp in prob.species
        dx = (sp.x_max - sp.x_birth) / n_mesh
        x_centers = [sp.x_birth + (i - 0.5) * dx for i in 1:n_mesh]
        push!(meshes, (x_centers=x_centers, dx=dx, x_b=sp.x_birth, x_m=sp.x_max))
        n0 = [sp.init_density(x) for x in x_centers]
        append!(u0, n0)
        sp_ranges[sp.name] = idx:(idx + n_mesh - 1)
        idx += n_mesh
    end
    return meshes, u0, sp_ranges
end

"""
    _build_fmu_rhs(prob, meshes, n_mesh, flux_type)

Build the RHS for fixed-mesh methods.

`flux_type` is `:upwind` (pure first-order upwind) or `:lax_friedrichs`
(Lax-Friedrichs flux with numerical diffusion α = max|g|).
"""
function _build_fmu_rhs(prob, meshes, n_mesh, flux_type::Symbol)
    function fmu_rhs!(du, u, p, t)
        E = prob.environment !== nothing ? prob.environment(u, t) : nothing
        offset = 0
        for (si, sp) in enumerate(prob.species)
            mesh = meshes[si]
            dx = mesh.dx
            n_m = n_mesh

            # Compute birth flux (boundary condition)
            birth_flux = 0.0
            for i in 1:n_m
                x = mesh.x_centers[i]
                birth_flux += sp.fecundity_rate(x, E, t) * u[offset + i] * dx
            end

            if flux_type == :lax_friedrichs
                # Compute max wave speed α = max|g(x)| across the mesh
                α = 0.0
                for i in 1:n_m
                    α = max(α, abs(sp.growth_rate(mesh.x_centers[i], E, t)))
                end

                for i in 1:n_m
                    x = mesh.x_centers[i]
                    g_i = sp.growth_rate(x, E, t)
                    μ_i = sp.mortality_rate(x, E, t)
                    n_i = u[offset + i]

                    # Left interface: exact upwind at birth boundary,
                    # Lax-Friedrichs at internal interfaces
                    if i == 1
                        # Exact boundary condition: g(x_b)·n(x_b) = birth_flux
                        g_b = sp.growth_rate(sp.x_birth, E, t)
                        flux_in = g_b > 0 ? birth_flux / (g_b + 1e-30) : 0.0
                        f_left = g_i * flux_in
                    else
                        x_l = mesh.x_centers[i-1]
                        g_l = sp.growth_rate(x_l, E, t)
                        n_l = u[offset + i - 1]
                        f_left = 0.5 * (g_l * n_l + g_i * n_i) - 0.5 * α * (n_i - n_l)
                    end

                    # Right interface: outflow boundary or LF flux
                    if i == n_m
                        f_right = g_i * n_i
                    else
                        x_r = mesh.x_centers[i+1]
                        g_r = sp.growth_rate(x_r, E, t)
                        n_r = u[offset + i + 1]
                        f_right = 0.5 * (g_i * n_i + g_r * n_r) - 0.5 * α * (n_r - n_i)
                    end

                    du[offset + i] = -(f_right - f_left) / dx - μ_i * n_i
                end
            else  # :upwind
                for i in 1:n_m
                    x = mesh.x_centers[i]
                    g = sp.growth_rate(x, E, t)
                    μ = sp.mortality_rate(x, E, t)

                    if i == 1
                        g_b = sp.growth_rate(sp.x_birth, E, t)
                        flux_in = g_b > 0 ? birth_flux / (g_b + 1e-30) : 0.0
                        flux_left = g * flux_in
                    else
                        flux_left = g * u[offset + i - 1]
                    end
                    flux_right = g * u[offset + i]

                    du[offset + i] = -(flux_right - flux_left) / dx - μ * u[offset + i]
                end
            end
            offset += n_m
        end
    end
    return fmu_rhs!
end

# ============================================================================
# Shared initialisation: EBT-family (EBT, IEBT)
# ============================================================================

function _init_ebt_state(prob, max_c)
    u0 = Float64[]
    sp_ranges = Dict{Symbol, UnitRange{Int}}()
    idx = 1
    for sp in prob.species
        dx = (sp.x_max - sp.x_birth) / max_c
        for c in 1:max_c
            x_c = sp.x_birth + (c - 0.5) * dx
            n_c = sp.init_density(x_c) * dx
            push!(u0, x_c)   # position (mean size)
            push!(u0, n_c)   # number
        end
        sp_ranges[sp.name] = idx:(idx + 2 * max_c - 1)
        idx += 2 * max_c
    end
    return u0, sp_ranges
end

function _build_ebt_rhs(prob, max_c)
    function ebt_rhs!(du, u, p, t)
        E = prob.environment !== nothing ? prob.environment(u, t) : nothing
        offset = 0
        for (si, sp) in enumerate(prob.species)
            for c in 1:max_c
                xi = offset + 2 * (c - 1) + 1
                ni = offset + 2 * (c - 1) + 2
                x = u[xi]
                n = u[ni]
                g = sp.growth_rate(x, E, t)
                μ = sp.mortality_rate(x, E, t)
                du[xi] = g          # dx/dt = growth rate
                du[ni] = -μ * n     # dn/dt = -mortality * n
            end
            offset += 2 * max_c
        end
    end
    return ebt_rhs!
end

# ============================================================================
# Characteristic Method (CM) — boundary-tracking
#
# State per cohort: (x_lower, x_upper, N_count) — 3 variables.
# Boundaries propagate at their local growth rates, preserving the
# characteristic structure of the PDE. The density within a cohort
# is n ≈ N / (x_upper − x_lower). Unlike EBT (which tracks mean
# size), CM tracks the full support of each cohort.
# ============================================================================

function _init_cm_state(prob, max_c)
    u0 = Float64[]
    sp_ranges = Dict{Symbol, UnitRange{Int}}()
    idx = 1
    for sp in prob.species
        dx = (sp.x_max - sp.x_birth) / max_c
        for c in 1:max_c
            x_lo = sp.x_birth + (c - 1) * dx
            x_hi = sp.x_birth + c * dx
            x_mid = 0.5 * (x_lo + x_hi)
            n_c = sp.init_density(x_mid) * dx
            push!(u0, x_lo)  # lower boundary
            push!(u0, x_hi)  # upper boundary
            push!(u0, n_c)   # number in cohort
        end
        sp_ranges[sp.name] = idx:(idx + 3 * max_c - 1)
        idx += 3 * max_c
    end
    return u0, sp_ranges
end

function _build_cm_rhs(prob, max_c)
    function cm_rhs!(du, u, p, t)
        E = prob.environment !== nothing ? prob.environment(u, t) : nothing
        offset = 0
        for (si, sp) in enumerate(prob.species)
            for c in 1:max_c
                lo_i = offset + 3 * (c - 1) + 1
                hi_i = offset + 3 * (c - 1) + 2
                n_i  = offset + 3 * (c - 1) + 3
                x_lo = u[lo_i]
                x_hi = u[hi_i]
                N    = u[n_i]

                # Boundaries propagate at their local growth rate
                g_lo = sp.growth_rate(x_lo, E, t)
                g_hi = sp.growth_rate(x_hi, E, t)

                # Mortality evaluated at the cohort midpoint
                x_mid = 0.5 * (x_lo + x_hi)
                μ_mid = sp.mortality_rate(x_mid, E, t)

                du[lo_i] = g_lo       # dx_lower/dt = g(x_lower)
                du[hi_i] = g_hi       # dx_upper/dt = g(x_upper)
                du[n_i]  = -μ_mid * N # dN/dt = -μ·N
            end
            offset += 3 * max_c
        end
    end
    return cm_rhs!
end

# ============================================================================
# 3. Staged PSPM species support
# ============================================================================

"""Check whether any species in the problem is a StagedPSPMSpecies."""
_has_staged(prob::PSPMProblem) = any(sp -> sp isa StagedPSPMSpecies, prob.species)

# ============================================================================
# 3a. Staged FMU initialisation
# ============================================================================

"""
    _init_staged_fmu_state(prob, n_mesh)

Initialise the FMU state vector for problems with `StagedPSPMSpecies`.

For a staged species with `ns` stages, each stage gets `n_mesh` cells,
giving `ns × n_mesh` state variables per species.  The species range
covers all stages.
"""
function _init_staged_fmu_state(prob, n_mesh)
    meshes = []      # list of (per-species) lists of per-stage mesh info
    u0 = Float64[]
    sp_ranges = Dict{Symbol, UnitRange{Int}}()
    idx = 1

    for sp in prob.species
        if sp isa StagedPSPMSpecies
            stage_meshes = []
            n_stages = length(sp.stages)
            sp_start = idx
            for (si, stage) in enumerate(sp.stages)
                dx = (stage.x_max - stage.x_birth) / n_mesh
                x_centers = [stage.x_birth + (i - 0.5) * dx for i in 1:n_mesh]
                push!(stage_meshes, (x_centers=x_centers, dx=dx,
                                     x_b=stage.x_birth, x_m=stage.x_max))
                n0 = [sp.init_density(si, x) for x in x_centers]
                append!(u0, n0)
                idx += n_mesh
            end
            push!(meshes, stage_meshes)
            sp_ranges[sp.name] = sp_start:(idx - 1)
        else
            # Regular PSPMSpecies — single PDE
            dx = (sp.x_max - sp.x_birth) / n_mesh
            x_centers = [sp.x_birth + (i - 0.5) * dx for i in 1:n_mesh]
            push!(meshes, [(x_centers=x_centers, dx=dx,
                            x_b=sp.x_birth, x_m=sp.x_max)])
            n0 = [sp.init_density(x) for x in x_centers]
            append!(u0, n0)
            sp_ranges[sp.name] = idx:(idx + n_mesh - 1)
            idx += n_mesh
        end
    end
    return meshes, u0, sp_ranges
end

# ============================================================================
# 3b. Staged FMU RHS builder
# ============================================================================

"""
    _build_staged_fmu_rhs(prob, meshes, n_mesh, flux_type)

Build the RHS for FMU methods with `StagedPSPMSpecies`.

For each staged species, the inter-stage boundary conditions are:
- Stage 1: influx = `sp.reproduction_flux(E, t, stage_totals)`
- Stage i>1: influx = outflow from stage i-1 at x = x_max

`flux_type` is `:upwind` or `:lax_friedrichs`.
"""
function _build_staged_fmu_rhs(prob, meshes, n_mesh, flux_type::Symbol)
    function staged_fmu_rhs!(du, u, p, t)
        E = prob.environment !== nothing ? prob.environment(u, t) : nothing
        offset = 0

        for (si, sp) in enumerate(prob.species)
            sp_meshes = meshes[si]

            if sp isa StagedPSPMSpecies
                n_stages = length(sp.stages)
                n_m = n_mesh

                # Compute stage totals for reproduction_flux and
                # density-dependent mortality (AD-safe: use eltype(u))
                stage_totals = zeros(eltype(u), n_stages)
                for s in 1:n_stages
                    s_off = offset + (s - 1) * n_m
                    dx_s = sp_meshes[s].dx
                    for j in 1:n_m
                        stage_totals[s] += u[s_off + j] * dx_s
                    end
                end

                for s in 1:n_stages
                    stage = sp.stages[s]
                    mesh = sp_meshes[s]
                    dx = mesh.dx
                    s_off = offset + (s - 1) * n_m

                    # Compute boundary influx
                    if s == 1
                        influx = sp.reproduction_flux(E, t, stage_totals)
                    else
                        # Outflow from previous stage: g_{s-1}(x_max) · n_{s-1}(x_max)
                        prev_stage = sp.stages[s - 1]
                        prev_mesh = sp_meshes[s - 1]
                        prev_off = offset + (s - 2) * n_m
                        g_prev_max = prev_stage.growth_rate(prev_mesh.x_m, E, t)
                        n_prev_max = u[prev_off + n_m]  # last cell of previous stage
                        influx = g_prev_max * n_prev_max
                    end

                    if flux_type == :lax_friedrichs
                        # Compute max wave speed for this stage
                        α = 0.0
                        for j in 1:n_m
                            x = mesh.x_centers[j]
                            α = max(α, abs(stage.growth_rate(x, E, t)))
                        end

                        for j in 1:n_m
                            x = mesh.x_centers[j]
                            g_j = stage.growth_rate(x, E, t)
                            μ_j = stage.mortality_rate(x, E, t)
                            n_j = u[s_off + j]

                            if j == 1
                                # Boundary cell: influx / dx
                                g_b = stage.growth_rate(mesh.x_b, E, t)
                                n_b = g_b > 0 ? influx / (g_b + 1e-30) : 0.0
                                f_left = g_j * n_b
                            else
                                x_l = mesh.x_centers[j - 1]
                                g_l = stage.growth_rate(x_l, E, t)
                                n_l = u[s_off + j - 1]
                                f_left = 0.5 * (g_l * n_l + g_j * n_j) -
                                         0.5 * α * (n_j - n_l)
                            end

                            if j == n_m
                                f_right = g_j * n_j
                            else
                                x_r = mesh.x_centers[j + 1]
                                g_r = stage.growth_rate(x_r, E, t)
                                n_r = u[s_off + j + 1]
                                f_right = 0.5 * (g_j * n_j + g_r * n_r) -
                                          0.5 * α * (n_r - n_j)
                            end

                            du[s_off + j] = -(f_right - f_left) / dx - μ_j * n_j
                        end
                    else  # :upwind
                        for j in 1:n_m
                            x = mesh.x_centers[j]
                            g = stage.growth_rate(x, E, t)
                            μ = stage.mortality_rate(x, E, t)

                            if j == 1
                                flux_left = influx
                            else
                                flux_left = g * u[s_off + j - 1]
                            end
                            flux_right = g * u[s_off + j]

                            du[s_off + j] = (flux_left - flux_right) / dx -
                                            μ * u[s_off + j]
                        end
                    end
                end
                offset += n_stages * n_m

            else
                # Regular PSPMSpecies — delegate to single-species FMU logic
                mesh = sp_meshes[1]
                dx = mesh.dx
                n_m = n_mesh

                birth_flux = 0.0
                for j in 1:n_m
                    x = mesh.x_centers[j]
                    birth_flux += sp.fecundity_rate(x, E, t) * u[offset + j] * dx
                end

                if flux_type == :lax_friedrichs
                    α = 0.0
                    for j in 1:n_m
                        α = max(α, abs(sp.growth_rate(mesh.x_centers[j], E, t)))
                    end
                    for j in 1:n_m
                        x = mesh.x_centers[j]
                        g_j = sp.growth_rate(x, E, t)
                        μ_j = sp.mortality_rate(x, E, t)
                        n_j = u[offset + j]
                        if j == 1
                            g_b = sp.growth_rate(sp.x_birth, E, t)
                            n_b = g_b > 0 ? birth_flux / (g_b + 1e-30) : 0.0
                            f_left = g_j * n_b
                        else
                            x_l = mesh.x_centers[j - 1]
                            g_l = sp.growth_rate(x_l, E, t)
                            n_l = u[offset + j - 1]
                            f_left = 0.5 * (g_l * n_l + g_j * n_j) -
                                     0.5 * α * (n_j - n_l)
                        end
                        if j == n_m
                            f_right = g_j * n_j
                        else
                            x_r = mesh.x_centers[j + 1]
                            g_r = sp.growth_rate(x_r, E, t)
                            n_r = u[offset + j + 1]
                            f_right = 0.5 * (g_j * n_j + g_r * n_r) -
                                      0.5 * α * (n_r - n_j)
                        end
                        du[offset + j] = -(f_right - f_left) / dx - μ_j * n_j
                    end
                else
                    for j in 1:n_m
                        x = mesh.x_centers[j]
                        g = sp.growth_rate(x, E, t)
                        μ = sp.mortality_rate(x, E, t)
                        if j == 1
                            g_b = sp.growth_rate(sp.x_birth, E, t)
                            flux_in = g_b > 0 ? birth_flux / (g_b + 1e-30) : 0.0
                            flux_left = g * flux_in
                        else
                            flux_left = g * u[offset + j - 1]
                        end
                        flux_right = g * u[offset + j]
                        du[offset + j] = -(flux_right - flux_left) / dx -
                                         μ * u[offset + j]
                    end
                end
                offset += n_m
            end
        end
    end
    return staged_fmu_rhs!
end

# ============================================================================
# 3c. Staged EBT initialisation and RHS
# ============================================================================

function _init_staged_ebt_state(prob, max_c)
    u0 = Float64[]
    sp_ranges = Dict{Symbol, UnitRange{Int}}()
    idx = 1

    for sp in prob.species
        if sp isa StagedPSPMSpecies
            sp_start = idx
            for (si, stage) in enumerate(sp.stages)
                dx = (stage.x_max - stage.x_birth) / max_c
                for c in 1:max_c
                    x_c = stage.x_birth + (c - 0.5) * dx
                    n_c = sp.init_density(si, x_c) * dx
                    push!(u0, x_c)   # position
                    push!(u0, n_c)   # number
                end
            end
            n_total = length(sp.stages) * 2 * max_c
            sp_ranges[sp.name] = sp_start:(sp_start + n_total - 1)
            idx = sp_start + n_total
        else
            dx = (sp.x_max - sp.x_birth) / max_c
            for c in 1:max_c
                x_c = sp.x_birth + (c - 0.5) * dx
                n_c = sp.init_density(x_c) * dx
                push!(u0, x_c)
                push!(u0, n_c)
            end
            sp_ranges[sp.name] = idx:(idx + 2 * max_c - 1)
            idx += 2 * max_c
        end
    end
    return u0, sp_ranges
end

function _build_staged_ebt_rhs(prob, max_c)
    function staged_ebt_rhs!(du, u, p, t)
        E = prob.environment !== nothing ? prob.environment(u, t) : nothing
        offset = 0

        for (si, sp) in enumerate(prob.species)
            if sp isa StagedPSPMSpecies
                n_stages = length(sp.stages)
                vars_per_stage = 2 * max_c

                for s in 1:n_stages
                    stage = sp.stages[s]
                    s_off = offset + (s - 1) * vars_per_stage
                    for c in 1:max_c
                        xi = s_off + 2 * (c - 1) + 1
                        ni = s_off + 2 * (c - 1) + 2
                        x = u[xi]
                        n = u[ni]
                        g = stage.growth_rate(x, E, t)
                        μ = stage.mortality_rate(x, E, t)
                        du[xi] = g
                        du[ni] = -μ * n
                    end
                end
                offset += n_stages * vars_per_stage
            else
                for c in 1:max_c
                    xi = offset + 2 * (c - 1) + 1
                    ni = offset + 2 * (c - 1) + 2
                    x = u[xi]
                    n = u[ni]
                    g = sp.growth_rate(x, E, t)
                    μ = sp.mortality_rate(x, E, t)
                    du[xi] = g
                    du[ni] = -μ * n
                end
                offset += 2 * max_c
            end
        end
    end
    return staged_ebt_rhs!
end

# ============================================================================
# 3d. Staged CM initialisation and RHS
# ============================================================================

function _init_staged_cm_state(prob, max_c)
    u0 = Float64[]
    sp_ranges = Dict{Symbol, UnitRange{Int}}()
    idx = 1

    for sp in prob.species
        if sp isa StagedPSPMSpecies
            sp_start = idx
            for (si, stage) in enumerate(sp.stages)
                dx = (stage.x_max - stage.x_birth) / max_c
                for c in 1:max_c
                    x_lo = stage.x_birth + (c - 1) * dx
                    x_hi = stage.x_birth + c * dx
                    x_mid = 0.5 * (x_lo + x_hi)
                    n_c = sp.init_density(si, x_mid) * dx
                    push!(u0, x_lo)
                    push!(u0, x_hi)
                    push!(u0, n_c)
                end
            end
            n_total = length(sp.stages) * 3 * max_c
            sp_ranges[sp.name] = sp_start:(sp_start + n_total - 1)
            idx = sp_start + n_total
        else
            dx = (sp.x_max - sp.x_birth) / max_c
            for c in 1:max_c
                x_lo = sp.x_birth + (c - 1) * dx
                x_hi = sp.x_birth + c * dx
                x_mid = 0.5 * (x_lo + x_hi)
                n_c = sp.init_density(x_mid) * dx
                push!(u0, x_lo)
                push!(u0, x_hi)
                push!(u0, n_c)
            end
            sp_ranges[sp.name] = idx:(idx + 3 * max_c - 1)
            idx += 3 * max_c
        end
    end
    return u0, sp_ranges
end

function _build_staged_cm_rhs(prob, max_c)
    function staged_cm_rhs!(du, u, p, t)
        E = prob.environment !== nothing ? prob.environment(u, t) : nothing
        offset = 0

        for (si, sp) in enumerate(prob.species)
            if sp isa StagedPSPMSpecies
                n_stages = length(sp.stages)
                vars_per_stage = 3 * max_c

                for s in 1:n_stages
                    stage = sp.stages[s]
                    s_off = offset + (s - 1) * vars_per_stage
                    for c in 1:max_c
                        lo_i = s_off + 3 * (c - 1) + 1
                        hi_i = s_off + 3 * (c - 1) + 2
                        n_i  = s_off + 3 * (c - 1) + 3
                        x_lo = u[lo_i]
                        x_hi = u[hi_i]
                        N    = u[n_i]
                        g_lo = stage.growth_rate(x_lo, E, t)
                        g_hi = stage.growth_rate(x_hi, E, t)
                        x_mid = 0.5 * (x_lo + x_hi)
                        μ_mid = stage.mortality_rate(x_mid, E, t)
                        du[lo_i] = g_lo
                        du[hi_i] = g_hi
                        du[n_i]  = -μ_mid * N
                    end
                end
                offset += n_stages * vars_per_stage
            else
                for c in 1:max_c
                    lo_i = offset + 3 * (c - 1) + 1
                    hi_i = offset + 3 * (c - 1) + 2
                    n_i  = offset + 3 * (c - 1) + 3
                    x_lo = u[lo_i]
                    x_hi = u[hi_i]
                    N    = u[n_i]
                    g_lo = sp.growth_rate(x_lo, E, t)
                    g_hi = sp.growth_rate(x_hi, E, t)
                    x_mid = 0.5 * (x_lo + x_hi)
                    μ_mid = sp.mortality_rate(x_mid, E, t)
                    du[lo_i] = g_lo
                    du[hi_i] = g_hi
                    du[n_i]  = -μ_mid * N
                end
                offset += 3 * max_c
            end
        end
    end
    return staged_cm_rhs!
end

# ============================================================================
# 3e. Staged PSPM solution accessors
# ============================================================================

"""
    staged_species_stage_totals(sol::ContinuousPBDMSolution,
                                 sp::StagedPSPMSpecies, n_mesh::Int)

Extract per-stage total abundance trajectories from a staged PSPM solution.

Returns a matrix `(n_timepoints × n_stages)` where each column is the total
abundance in that life stage over time.
"""
function PhysiologicallyBasedDemographicModels.staged_species_stage_totals(
        sol::ContinuousPBDMSolution,
        sp::StagedPSPMSpecies, n_mesh::Int)
    r = sol.species_ranges[sp.name]
    n_stages = length(sp.stages)
    n_t = length(sol.t)
    totals = zeros(n_t, n_stages)
    for s in 1:n_stages
        stage = sp.stages[s]
        dx = (stage.x_max - stage.x_birth) / n_mesh
        s_start = r.start + (s - 1) * n_mesh
        for k in 1:n_t
            for j in 1:n_mesh
                totals[k, s] += sol.u[s_start + j - 1, k] * dx
            end
        end
    end
    return totals
end

end # module OrdinaryDiffEqExt
