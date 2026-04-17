using Test
using PhysiologicallyBasedDemographicModels
import PhysiologicallyBasedDemographicModels: _total_substages

@testset "Continuous-time PBDM formulations" begin

    @testset "ContinuousSpecies construction" begin
        # From scratch
        sp = ContinuousSpecies(:plant;
            k=[3, 3], τ=[100.0, 200.0], μ=[0.01, 0.02],
            dev_rate=[LinearDevelopmentRate(10.0, 35.0), LinearDevelopmentRate(10.0, 35.0)],
            fr=FraserGilbertResponse(0.8),
            resp=Q10Respiration(0.01, 2.0, 25.0),
            demand_rate=0.5,
            intrinsic_rate=0.3, carrying_capacity=5000.0)

        @test sp.name == :plant
        @test sp.n_stages == 2
        @test sp.k == [3, 3]
        @test _total_substages(sp) == 6
        @test sp.intrinsic_rate == 0.3
        @test sp.carrying_capacity == 5000.0

        # From a Population object
        pop = Population(:test_pop, [
            LifeStage(:egg, DistributedDelay(3, 50.0; W0=10.0),
                      LinearDevelopmentRate(10.0, 35.0), 0.01),
            LifeStage(:larva, DistributedDelay(5, 150.0; W0=5.0),
                      LinearDevelopmentRate(10.0, 35.0), 0.02)
        ])

        sp2 = ContinuousSpecies(pop;
            fr=FraserGilbertResponse(0.8),
            resp=Q10Respiration(0.02, 2.0, 25.0),
            demand_rate=0.3)

        @test sp2.name == :test_pop
        @test sp2.n_stages == 2
        @test sp2.k == [3, 5]
        @test sp2.τ == [50.0, 150.0]
        @test sp2.μ == [0.01, 0.02]
        @test _total_substages(sp2) == 8
    end

    @testset "ContinuousTrophicLink" begin
        link = ContinuousTrophicLink(:herbivore, :plant)
        @test link.consumer == :herbivore
        @test link.resource == :plant
    end

    @testset "flatten_population" begin
        pop = Population(:test, [
            LifeStage(:s1, DistributedDelay(3, 50.0; W0=10.0),
                      LinearDevelopmentRate(10.0, 35.0), 0.0),
            LifeStage(:s2, DistributedDelay(2, 100.0; W0=5.0),
                      LinearDevelopmentRate(10.0, 35.0), 0.0)
        ])
        u = flatten_population(pop)
        @test length(u) == 5  # 3 + 2 substages
        @test u[1] ≈ 10.0
        @test u[4] ≈ 5.0
    end

    @testset "ContinuousPBDMProblem construction" begin
        plant = ContinuousSpecies(:plant;
            k=[3], τ=[100.0], μ=[0.0],
            dev_rate=[LinearDevelopmentRate(10.0, 35.0)],
            fr=FraserGilbertResponse(0.8),
            resp=Q10Respiration(0.01, 2.0, 25.0),
            demand_rate=0.5,
            intrinsic_rate=0.5, carrying_capacity=1000.0)

        prob = ContinuousPBDMProblem(
            species=[plant],
            u0=[100.0, 100.0, 100.0],
            tspan=(0.0, 365.0),
            T_forcing=25.0)

        @test length(prob.u0) == 3
        @test prob.tspan == (0.0, 365.0)

        ranges = species_state_ranges(prob)
        @test length(ranges) == 1
        @test ranges[1] == (:plant, 1, 1:3)

        tr = species_total_ranges(prob)
        @test tr[:plant] == 1:3
    end

    @testset "DelayPBDMProblem construction" begin
        sp = ContinuousSpecies(:grass;
            k=[3], τ=[80.0], μ=[0.01],
            dev_rate=[LinearDevelopmentRate(10.0, 35.0)],
            fr=FraserGilbertResponse(0.8),
            resp=Q10Respiration(0.01, 2.0, 25.0),
            demand_rate=0.3,
            intrinsic_rate=0.5, carrying_capacity=2000.0)

        prob = DelayPBDMProblem(
            species=[sp],
            u0=[500.0],
            tspan=(0.0, 100.0),
            T_forcing=25.0)

        @test length(prob.u0) == 1
        @test prob.tspan == (0.0, 100.0)
        @test prob.h0 isa Function  # default history function
    end

    @testset "PSPMSpecies and PSPMProblem construction" begin
        sp = PSPMSpecies(:daphnia;
            x_birth=0.1, x_max=10.0,
            growth_rate=(x, E, t) -> 0.5 * (1 - x/10.0),
            mortality_rate=(x, E, t) -> 0.01 + 0.001 * x,
            fecundity_rate=(x, E, t) -> x > 3.0 ? 0.1 * x : 0.0,
            init_density=(x) -> x < 5.0 ? 1.0 : 0.0)

        @test sp.name == :daphnia
        @test sp.x_birth ≈ 0.1
        @test sp.x_max ≈ 10.0

        # FMU method
        prob_fmu = PSPMProblem(
            species=[sp],
            method=FixedMeshUpwind(n_mesh=50),
            tspan=(0.0, 100.0))

        @test prob_fmu.method isa FixedMeshUpwind
        @test prob_fmu.method.n_mesh == 50

        # EBT method
        prob_ebt = PSPMProblem(
            species=[sp],
            method=EscalatorBoxcarTrain(max_cohorts=100),
            tspan=(0.0, 100.0))

        @test prob_ebt.method isa EscalatorBoxcarTrain

        # CM method
        prob_cm = PSPMProblem(
            species=[sp],
            method=CharacteristicMethod(max_cohorts=100),
            tspan=(0.0, 100.0))

        @test prob_cm.method isa CharacteristicMethod
    end

    @testset "ContinuousPBDMSolution accessors" begin
        # Create a mock solution
        t = collect(0.0:1.0:10.0)
        u = rand(6, 11)
        sp_ranges = Dict(:plant => 1:3, :pest => 4:6)
        sol = ContinuousPBDMSolution{Float64, Nothing}(
            t, u, [:plant, :pest], sp_ranges, :Success, nothing)

        traj_plant = species_trajectory(sol, :plant)
        @test length(traj_plant) == 11
        @test traj_plant[1] ≈ sum(u[1:3, 1])

        traj_pest = species_trajectory(sol, :pest)
        @test length(traj_pest) == 11
        @test traj_pest[end] ≈ sum(u[4:6, end])

        @test_throws KeyError species_trajectory(sol, :nonexistent)

        # show method
        buf = IOBuffer()
        show(buf, sol)
        s = String(take!(buf))
        @test occursin("11 timepoints", s)
        @test occursin("2 species", s)
    end

    # ====================================================================
    # Integration tests with actual solvers (require extensions)
    # ====================================================================

    @testset "ODE solve (linear chain trick)" begin
        # Single basal resource with logistic growth, 3 substages
        plant = ContinuousSpecies(:plant;
            k=[3], τ=[100.0], μ=[0.0],
            dev_rate=[LinearDevelopmentRate(10.0, 35.0)],
            fr=FraserGilbertResponse(0.8),
            resp=Q10Respiration(0.005, 2.0, 25.0),
            demand_rate=0.5,
            intrinsic_rate=0.5, carrying_capacity=1000.0)

        prob = ContinuousPBDMProblem(
            species=[plant],
            u0=[100.0, 100.0, 100.0],
            tspan=(0.0, 200.0),
            T_forcing=25.0)

        try
            using OrdinaryDiffEq
            sol = solve_continuous(prob; reltol=1e-6, abstol=1e-8)

            @test sol.retcode == :Success || sol.retcode == :ReturnCode_Success
            @test length(sol.t) > 1
            @test size(sol.u, 1) == 3
            @test size(sol.u, 2) == length(sol.t)

            # Plant should grow toward carrying capacity
            plant_traj = species_trajectory(sol, :plant)
            @test plant_traj[end] > plant_traj[1]  # grows
            @test plant_traj[end] < 1100.0  # bounded by K

            # Raw solution accessible
            @test sol.raw_sol !== nothing
        catch e
            if isa(e, ArgumentError) && occursin("OrdinaryDiffEq", string(e))
                @warn "OrdinaryDiffEq not loaded, skipping ODE solve test"
            else
                rethrow(e)
            end
        end
    end

    @testset "ODE solve — tritrophic system" begin
        dev = LinearDevelopmentRate(10.0, 35.0)
        fr = FraserGilbertResponse(0.8)

        plant = ContinuousSpecies(:plant;
            k=[3], τ=[100.0], μ=[0.0],
            dev_rate=[dev], fr=fr,
            resp=Q10Respiration(0.005, 2.0, 25.0),
            demand_rate=0.5,
            intrinsic_rate=1.0, carrying_capacity=2000.0)

        herbivore = ContinuousSpecies(:herbivore;
            k=[3], τ=[80.0], μ=[0.01],
            dev_rate=[dev], fr=fr,
            resp=Q10Respiration(0.02, 2.0, 25.0),
            demand_rate=0.3,
            conversion_efficiency=0.4)

        predator = ContinuousSpecies(:predator;
            k=[3], τ=[120.0], μ=[0.02],
            dev_rate=[dev], fr=fr,
            resp=Q10Respiration(0.03, 2.0, 25.0),
            demand_rate=0.2,
            conversion_efficiency=0.3)

        links = [
            ContinuousTrophicLink(:herbivore, :plant),
            ContinuousTrophicLink(:predator, :herbivore)
        ]

        prob = ContinuousPBDMProblem(
            species=[plant, herbivore, predator],
            links=links,
            u0=[200.0, 200.0, 200.0,   # plant substages
                30.0, 30.0, 30.0,       # herbivore substages
                5.0, 5.0, 5.0],         # predator substages
            tspan=(0.0, 365.0),
            T_forcing=25.0)

        try
            using OrdinaryDiffEq
            sol = solve_continuous(prob; reltol=1e-6, abstol=1e-8)

            @test sol.retcode == :Success || sol.retcode == :ReturnCode_Success
            @test :plant in sol.species_names
            @test :herbivore in sol.species_names
            @test :predator in sol.species_names

            # All species should have non-negative biomass
            plant_traj = species_trajectory(sol, :plant)
            herb_traj = species_trajectory(sol, :herbivore)
            pred_traj = species_trajectory(sol, :predator)

            @test all(plant_traj .>= -1e-6)  # allow tiny numerical noise
            @test all(herb_traj .>= -1e-6)
            @test all(pred_traj .>= -1e-6)

            # Tritrophic regulation: plant > herbivore > predator in biomass
            # (typical for biomass pyramids)
            @test plant_traj[end] > herb_traj[end]
        catch e
            if isa(e, ArgumentError) && occursin("OrdinaryDiffEq", string(e))
                @warn "OrdinaryDiffEq not loaded, skipping tritrophic ODE test"
            else
                rethrow(e)
            end
        end
    end

    @testset "DDE solve — single species logistic" begin
        sp = ContinuousSpecies(:grass;
            k=[3], τ=[80.0], μ=[0.0],
            dev_rate=[LinearDevelopmentRate(10.0, 35.0)],
            fr=FraserGilbertResponse(0.8),
            resp=Q10Respiration(0.005, 2.0, 25.0),
            demand_rate=0.3,
            intrinsic_rate=0.5, carrying_capacity=2000.0)

        prob = DelayPBDMProblem(
            species=[sp],
            u0=[500.0],
            tspan=(0.0, 100.0),
            T_forcing=25.0)

        try
            using DelayDiffEq
            sol = solve_delay(prob; reltol=1e-4, abstol=1e-6)

            @test sol.retcode == :Success || sol.retcode == :ReturnCode_Success
            @test length(sol.t) > 1

            grass_traj = species_trajectory(sol, :grass)
            @test grass_traj[end] > grass_traj[1]  # grows toward K
        catch e
            if isa(e, ArgumentError) && occursin("DelayDiffEq", string(e))
                @warn "DelayDiffEq not loaded, skipping DDE solve test"
            else
                rethrow(e)
            end
        end
    end

    @testset "PSPM solve (FMU) — size-structured population" begin
        sp = PSPMSpecies(:daphnia;
            x_birth=0.5, x_max=10.0,
            growth_rate=(x, E, t) -> 0.2 * (1 - x / 10.0),
            mortality_rate=(x, E, t) -> 0.01,
            fecundity_rate=(x, E, t) -> x > 3.0 ? 0.05 * x : 0.0,
            init_density=(x) -> x < 5.0 ? 10.0 : 0.0)

        prob = PSPMProblem(
            species=[sp],
            method=FixedMeshUpwind(n_mesh=50),
            tspan=(0.0, 50.0))

        try
            using OrdinaryDiffEq
            sol = solve_pspm(prob; reltol=1e-4, abstol=1e-6)

            @test sol.retcode == :Success || sol.retcode == :ReturnCode_Success
            @test length(sol.t) > 1
            @test size(sol.u, 1) == 50  # 50 mesh points
            @test :daphnia in sol.species_names
        catch e
            if isa(e, ArgumentError) && occursin("OrdinaryDiffEq", string(e))
                @warn "OrdinaryDiffEq not loaded, skipping PSPM test"
            else
                rethrow(e)
            end
        end
    end

    @testset "PSPM solve (EBT) — size-structured population" begin
        sp = PSPMSpecies(:daphnia;
            x_birth=0.5, x_max=10.0,
            growth_rate=(x, E, t) -> 0.2 * (1 - x / 10.0),
            mortality_rate=(x, E, t) -> 0.01,
            fecundity_rate=(x, E, t) -> x > 3.0 ? 0.05 * x : 0.0,
            init_density=(x) -> x < 5.0 ? 10.0 : 0.0)

        prob = PSPMProblem(
            species=[sp],
            method=EscalatorBoxcarTrain(max_cohorts=30),
            tspan=(0.0, 50.0))

        try
            using OrdinaryDiffEq
            sol = solve_pspm(prob; reltol=1e-4, abstol=1e-6)

            @test sol.retcode == :Success || sol.retcode == :ReturnCode_Success
            @test length(sol.t) > 1
            @test size(sol.u, 1) == 60  # 30 cohorts × 2 (position + number)
        catch e
            if isa(e, ArgumentError) && occursin("OrdinaryDiffEq", string(e))
                @warn "OrdinaryDiffEq not loaded, skipping EBT test"
            else
                rethrow(e)
            end
        end
    end

    # =========================================================================
    # Remaining PSPM methods (Joshi et al. 2023)
    # =========================================================================

    # Shared species definition for all PSPM method tests
    function _make_pspm_test_species()
        PSPMSpecies(:daphnia;
            x_birth=0.5, x_max=10.0,
            growth_rate=(x, E, t) -> 0.2 * (1 - x / 10.0),
            mortality_rate=(x, E, t) -> 0.01,
            fecundity_rate=(x, E, t) -> x > 3.0 ? 0.05 * x : 0.0,
            init_density=(x) -> x < 5.0 ? 10.0 : 0.0)
    end

    # Stiff species with very high size-dependent mortality (needs implicit)
    function _make_stiff_pspm_species()
        PSPMSpecies(:stiff_sp;
            x_birth=0.1, x_max=5.0,
            growth_rate=(x, E, t) -> 0.5 * (1 - x / 5.0),
            mortality_rate=(x, E, t) -> 0.5 + 10.0 * max(0.0, x - 3.0),
            fecundity_rate=(x, E, t) -> x > 1.5 ? 0.2 * x : 0.0,
            init_density=(x) -> exp(-0.5 * (x - 1.0)^2))
    end

    @testset "PSPM solve (IFMU) — implicit fixed mesh upwind" begin
        sp = _make_pspm_test_species()
        prob = PSPMProblem(
            species=[sp],
            method=ImplicitFixedMeshUpwind(n_mesh=50),
            tspan=(0.0, 50.0))

        try
            using OrdinaryDiffEq
            sol = solve_pspm(prob; reltol=1e-4, abstol=1e-6)
            @test sol.retcode == :Success || sol.retcode == :ReturnCode_Success
            @test length(sol.t) > 1
            @test size(sol.u, 1) == 50
        catch e
            if isa(e, ArgumentError) && occursin("OrdinaryDiffEq", string(e))
                @warn "OrdinaryDiffEq not loaded, skipping IFMU test"
            else
                rethrow(e)
            end
        end
    end

    @testset "PSPM solve (IFMU) — stiff mortality" begin
        sp = _make_stiff_pspm_species()
        prob = PSPMProblem(
            species=[sp],
            method=ImplicitFixedMeshUpwind(n_mesh=40),
            tspan=(0.0, 20.0))

        try
            using OrdinaryDiffEq
            sol = solve_pspm(prob; reltol=1e-4, abstol=1e-6)
            @test sol.retcode == :Success || sol.retcode == :ReturnCode_Success
            @test length(sol.t) > 1
        catch e
            if isa(e, ArgumentError) && occursin("OrdinaryDiffEq", string(e))
                @warn "OrdinaryDiffEq not loaded, skipping IFMU stiff test"
            else
                rethrow(e)
            end
        end
    end

    @testset "PSPM solve (FMU-LF) — Lax-Friedrichs upwind" begin
        sp = _make_pspm_test_species()
        prob = PSPMProblem(
            species=[sp],
            method=LaxFriedrichsUpwind(n_mesh=50),
            tspan=(0.0, 50.0))

        try
            using OrdinaryDiffEq
            sol = solve_pspm(prob; reltol=1e-4, abstol=1e-6)
            @test sol.retcode == :Success || sol.retcode == :ReturnCode_Success
            @test length(sol.t) > 1
            @test size(sol.u, 1) == 50
        catch e
            if isa(e, ArgumentError) && occursin("OrdinaryDiffEq", string(e))
                @warn "OrdinaryDiffEq not loaded, skipping FMU-LF test"
            else
                rethrow(e)
            end
        end
    end

    @testset "PSPM solve (IFMU-LF) — implicit Lax-Friedrichs" begin
        sp = _make_stiff_pspm_species()
        prob = PSPMProblem(
            species=[sp],
            method=ImplicitLaxFriedrichsUpwind(n_mesh=40),
            tspan=(0.0, 20.0))

        try
            using OrdinaryDiffEq
            sol = solve_pspm(prob; reltol=1e-4, abstol=1e-6)
            @test sol.retcode == :Success || sol.retcode == :ReturnCode_Success
            @test length(sol.t) > 1
        catch e
            if isa(e, ArgumentError) && occursin("OrdinaryDiffEq", string(e))
                @warn "OrdinaryDiffEq not loaded, skipping IFMU-LF test"
            else
                rethrow(e)
            end
        end
    end

    @testset "PSPM solve (IEBT) — implicit escalator boxcar train" begin
        sp = _make_stiff_pspm_species()
        prob = PSPMProblem(
            species=[sp],
            method=ImplicitEscalatorBoxcarTrain(max_cohorts=30),
            tspan=(0.0, 20.0))

        try
            using OrdinaryDiffEq
            sol = solve_pspm(prob; reltol=1e-4, abstol=1e-6)
            @test sol.retcode == :Success || sol.retcode == :ReturnCode_Success
            @test length(sol.t) > 1
            @test size(sol.u, 1) == 60  # 30 cohorts × 2
        catch e
            if isa(e, ArgumentError) && occursin("OrdinaryDiffEq", string(e))
                @warn "OrdinaryDiffEq not loaded, skipping IEBT test"
            else
                rethrow(e)
            end
        end
    end

    @testset "PSPM solve (CM) — characteristic method" begin
        sp = _make_pspm_test_species()
        prob = PSPMProblem(
            species=[sp],
            method=CharacteristicMethod(max_cohorts=30),
            tspan=(0.0, 50.0))

        try
            using OrdinaryDiffEq
            sol = solve_pspm(prob; reltol=1e-4, abstol=1e-6)
            @test sol.retcode == :Success || sol.retcode == :ReturnCode_Success
            @test length(sol.t) > 1
            @test size(sol.u, 1) == 90  # 30 cohorts × 3 (lo, hi, N)
        catch e
            if isa(e, ArgumentError) && occursin("OrdinaryDiffEq", string(e))
                @warn "OrdinaryDiffEq not loaded, skipping CM test"
            else
                rethrow(e)
            end
        end
    end

    @testset "PSPM solve (ICM) — implicit characteristic method" begin
        sp = _make_stiff_pspm_species()
        prob = PSPMProblem(
            species=[sp],
            method=ImplicitCharacteristicMethod(max_cohorts=30),
            tspan=(0.0, 20.0))

        try
            using OrdinaryDiffEq
            sol = solve_pspm(prob; reltol=1e-4, abstol=1e-6)
            @test sol.retcode == :Success || sol.retcode == :ReturnCode_Success
            @test length(sol.t) > 1
            @test size(sol.u, 1) == 90  # 30 cohorts × 3
        catch e
            if isa(e, ArgumentError) && occursin("OrdinaryDiffEq", string(e))
                @warn "OrdinaryDiffEq not loaded, skipping ICM test"
            else
                rethrow(e)
            end
        end
    end

    @testset "PSPM method consistency — FMU-LF mesh convergence" begin
        # LF adds numerical diffusion that interacts with fecundity,
        # so it differs from FMU for coupled problems. Instead, test
        # that LF converges as mesh refines (total pop stabilizes).
        sp = _make_pspm_test_species()

        try
            using OrdinaryDiffEq
            pops = Float64[]
            for nm in [50, 100, 200]
                prob_lf = PSPMProblem(
                    species=[sp],
                    method=LaxFriedrichsUpwind(n_mesh=nm),
                    tspan=(0.0, 30.0))
                sol = solve_pspm(prob_lf; reltol=1e-6, abstol=1e-8)
                dx = (sp.x_max - sp.x_birth) / nm
                push!(pops, sum(sol.u[:, end]) * dx)
            end
            # Refinement should reduce the difference between successive meshes
            diff1 = abs(pops[2] - pops[1])
            diff2 = abs(pops[3] - pops[2])
            @test diff2 < diff1  # convergence: finer meshes agree more
        catch e
            if isa(e, ArgumentError) && occursin("OrdinaryDiffEq", string(e))
                @warn "OrdinaryDiffEq not loaded, skipping LF convergence test"
            else
                rethrow(e)
            end
        end
    end

    @testset "PSPM method consistency — FMU vs IFMU give same result" begin
        sp = _make_pspm_test_species()

        prob_fmu = PSPMProblem(
            species=[sp],
            method=FixedMeshUpwind(n_mesh=50),
            tspan=(0.0, 30.0))

        prob_ifmu = PSPMProblem(
            species=[sp],
            method=ImplicitFixedMeshUpwind(n_mesh=50),
            tspan=(0.0, 30.0))

        try
            using OrdinaryDiffEq
            sol_fmu  = solve_pspm(prob_fmu;  reltol=1e-6, abstol=1e-8)
            sol_ifmu = solve_pspm(prob_ifmu; reltol=1e-6, abstol=1e-8)

            # Both use the same RHS — only the ODE solver differs
            # Total population should agree closely for non-stiff problems
            dx = (sp.x_max - sp.x_birth) / 50
            pop_fmu  = sum(sol_fmu.u[:, end]) * dx
            pop_ifmu = sum(sol_ifmu.u[:, end]) * dx

            @test abs(pop_fmu - pop_ifmu) / max(abs(pop_fmu), 1e-10) < 0.05
        catch e
            if isa(e, ArgumentError) && occursin("OrdinaryDiffEq", string(e))
                @warn "OrdinaryDiffEq not loaded, skipping FMU/IFMU consistency test"
            else
                rethrow(e)
            end
        end
    end

    @testset "PSPM construction — new method types" begin
        sp = _make_pspm_test_species()

        p1 = PSPMProblem(species=[sp], method=ImplicitFixedMeshUpwind(n_mesh=60), tspan=(0.0, 10.0))
        @test p1.method isa ImplicitFixedMeshUpwind
        @test p1.method.n_mesh == 60

        p2 = PSPMProblem(species=[sp], method=LaxFriedrichsUpwind(n_mesh=80), tspan=(0.0, 10.0))
        @test p2.method isa LaxFriedrichsUpwind

        p3 = PSPMProblem(species=[sp], method=ImplicitLaxFriedrichsUpwind(n_mesh=70), tspan=(0.0, 10.0))
        @test p3.method isa ImplicitLaxFriedrichsUpwind

        p4 = PSPMProblem(species=[sp], method=ImplicitEscalatorBoxcarTrain(max_cohorts=50), tspan=(0.0, 10.0))
        @test p4.method isa ImplicitEscalatorBoxcarTrain

        p5 = PSPMProblem(species=[sp], method=ImplicitCharacteristicMethod(max_cohorts=40), tspan=(0.0, 10.0))
        @test p5.method isa ImplicitCharacteristicMethod
    end

    @testset "StagedPSPMSpecies construction" begin
        stages = [
            PSPMStage(:egg;   growth_rate=(x,E,t)->0.5, mortality_rate=(x,E,t)->0.01),
            PSPMStage(:larva; growth_rate=(x,E,t)->0.3, mortality_rate=(x,E,t)->0.02),
            PSPMStage(:adult; x_birth=0.0, x_max=1.0,
                      growth_rate=(x,E,t)->0.1, mortality_rate=(x,E,t)->0.05),
        ]

        sp = StagedPSPMSpecies(:insect;
            stages = stages,
            reproduction_flux = (E, t, totals) -> 5.0 * totals[3],
            init_density = (stage, x) -> stage == 3 ? 1.0 : 0.0)

        @test sp.name == :insect
        @test n_pspm_stages(sp) == 3
        @test sp.stages[1].name == :egg
        @test sp.stages[2].name == :larva
        @test sp.stages[3].name == :adult
        @test sp.stages[1].x_birth ≈ 0.0
        @test sp.stages[1].x_max ≈ 1.0

        # Type hierarchy
        @test sp isa AbstractPSPMSpecies{Float64}

        # Construction via PSPMProblem
        prob = PSPMProblem(
            species = [sp],
            method = FixedMeshUpwind(n_mesh=20),
            tspan = (0.0, 50.0))
        @test prob.method isa FixedMeshUpwind
        @test length(prob.species) == 1
    end

    @testset "StagedPSPMSpecies solve — FMU family" begin
        # 2-stage species with known analytic behavior:
        # constant growth and mortality, no density dependence
        stages = [
            PSPMStage(:juvenile; growth_rate=(x,E,t)->0.5, mortality_rate=(x,E,t)->0.01),
            PSPMStage(:adult;    growth_rate=(x,E,t)->0.2, mortality_rate=(x,E,t)->0.02),
        ]

        sp = StagedPSPMSpecies(:test_sp;
            stages = stages,
            reproduction_flux = (E, t, totals) -> 2.0 * totals[2],
            init_density = (stage, x) -> stage == 2 ? 5.0 : 0.0)

        # FMU explicit
        prob = PSPMProblem(species=[sp], method=FixedMeshUpwind(n_mesh=30),
                           tspan=(0.0, 20.0))
        sol = solve_pspm(prob; reltol=1e-6, abstol=1e-8, saveat=1.0)
        @test sol.retcode == :Success
        @test length(sol.t) == 21
        @test :test_sp in sol.species_names

        totals = staged_species_stage_totals(sol, sp, 30)
        @test size(totals) == (21, 2)
        @test totals[1, 2] > 0   # adults present initially
        @test totals[end, 1] > 0  # juveniles produced via reproduction

        # IFMU (stiff solver)
        prob2 = PSPMProblem(species=[sp], method=ImplicitFixedMeshUpwind(n_mesh=30),
                            tspan=(0.0, 20.0))
        sol2 = solve_pspm(prob2; reltol=1e-6, abstol=1e-8, saveat=1.0)
        @test sol2.retcode == :Success

        # Lax-Friedrichs
        prob3 = PSPMProblem(species=[sp], method=LaxFriedrichsUpwind(n_mesh=30),
                            tspan=(0.0, 20.0))
        sol3 = solve_pspm(prob3; reltol=1e-6, abstol=1e-8, saveat=1.0)
        @test sol3.retcode == :Success

        # Implicit Lax-Friedrichs
        prob4 = PSPMProblem(species=[sp], method=ImplicitLaxFriedrichsUpwind(n_mesh=30),
                            tspan=(0.0, 20.0))
        sol4 = solve_pspm(prob4; reltol=1e-6, abstol=1e-8, saveat=1.0)
        @test sol4.retcode == :Success
    end

    @testset "StagedPSPMSpecies solve — EBT and CM" begin
        stages = [
            PSPMStage(:juvenile; growth_rate=(x,E,t)->0.5, mortality_rate=(x,E,t)->0.01),
            PSPMStage(:adult;    growth_rate=(x,E,t)->0.2, mortality_rate=(x,E,t)->0.02),
        ]

        sp = StagedPSPMSpecies(:test_sp;
            stages = stages,
            reproduction_flux = (E, t, totals) -> 2.0 * totals[2],
            init_density = (stage, x) -> stage == 2 ? 5.0 : 0.0)

        # EBT
        prob_ebt = PSPMProblem(species=[sp], method=EscalatorBoxcarTrain(max_cohorts=20),
                               tspan=(0.0, 20.0))
        sol_ebt = solve_pspm(prob_ebt; reltol=1e-6, abstol=1e-8, saveat=1.0)
        @test sol_ebt.retcode == :Success

        # Implicit EBT
        prob_iebt = PSPMProblem(species=[sp], method=ImplicitEscalatorBoxcarTrain(max_cohorts=20),
                                tspan=(0.0, 20.0))
        sol_iebt = solve_pspm(prob_iebt; reltol=1e-6, abstol=1e-8, saveat=1.0)
        @test sol_iebt.retcode == :Success

        # CM
        prob_cm = PSPMProblem(species=[sp], method=CharacteristicMethod(max_cohorts=20),
                              tspan=(0.0, 20.0))
        sol_cm = solve_pspm(prob_cm; reltol=1e-6, abstol=1e-8, saveat=1.0)
        @test sol_cm.retcode == :Success

        # Implicit CM
        prob_icm = PSPMProblem(species=[sp], method=ImplicitCharacteristicMethod(max_cohorts=20),
                               tspan=(0.0, 20.0))
        sol_icm = solve_pspm(prob_icm; reltol=1e-6, abstol=1e-8, saveat=1.0)
        @test sol_icm.retcode == :Success
    end

    @testset "StagedPSPMSpecies with environment" begin
        # Use environment for temperature-dependent rates
        stages = [
            PSPMStage(:egg;   growth_rate=(x,E,t) -> E.T > 10.0 ? 0.5 : 0.0,
                              mortality_rate=(x,E,t) -> 0.01 + 0.001 * (E.T - 20.0)^2),
            PSPMStage(:adult; growth_rate=(x,E,t) -> E.T > 10.0 ? 0.2 : 0.0,
                              mortality_rate=(x,E,t) -> 0.02 + 0.001 * (E.T - 20.0)^2),
        ]

        sp = StagedPSPMSpecies(:env_sp;
            stages = stages,
            reproduction_flux = (E, t, totals) -> 3.0 * totals[2],
            init_density = (stage, x) -> stage == 2 ? 1.0 : 0.0)

        env_func(u, t) = (T = 20.0 + 5.0 * sin(2π * t / 365),)

        prob = PSPMProblem(species=[sp], environment=env_func,
                           method=FixedMeshUpwind(n_mesh=20),
                           tspan=(0.0, 365.0))
        sol = solve_pspm(prob; reltol=1e-6, abstol=1e-8, saveat=10.0)
        @test sol.retcode == :Success
        @test length(sol.t) > 1

        totals = staged_species_stage_totals(sol, sp, 20)
        @test size(totals, 2) == 2
    end
end
