@testset "Core Types" begin
    @testset "LinearDevelopmentRate" begin
        ldr = LinearDevelopmentRate(10.0, 35.0)
        @test ldr.T_lower == 10.0
        @test ldr.T_upper == 35.0
        @test development_rate(ldr, 5.0) == 0.0    # Below threshold
        @test development_rate(ldr, 10.0) == 0.0   # At threshold
        @test development_rate(ldr, 20.0) == 10.0   # 20 - 10 = 10 DD
        @test development_rate(ldr, 40.0) == 25.0   # Capped at T_upper - T_lower
        @test degree_days(ldr, 20.0) == 10.0
        @test degree_days(ldr, 5.0) == 0.0
        @test_throws ArgumentError LinearDevelopmentRate(35.0, 10.0)  # lower > upper
    end

    @testset "BriereDevelopmentRate" begin
        bdr = BriereDevelopmentRate(1e-4, 6.0, 35.0)
        @test development_rate(bdr, 5.0) == 0.0
        @test development_rate(bdr, 36.0) == 0.0
        @test development_rate(bdr, 20.0) > 0.0
        @test_throws ArgumentError BriereDevelopmentRate(-1.0, 6.0, 35.0)
    end

    @testset "LoganDevelopmentRate" begin
        ldr = LoganDevelopmentRate(0.02, 0.15, 35.0, 5.0)
        r = development_rate(ldr, 20.0)
        @test r isa Float64
        @test_throws ArgumentError LoganDevelopmentRate(-0.02, 0.15, 35.0, 5.0)
        @test_throws ArgumentError LoganDevelopmentRate(0.02, 0.0, 35.0, 5.0)
        @test_throws ArgumentError LoganDevelopmentRate(0.02, 0.15, 35.0, 0.0)
    end

    @testset "DistributedDelay" begin
        dd = DistributedDelay(30, 750.0)
        @test dd.k == 30
        @test dd.τ == 750.0
        @test length(dd.W) == 30
        @test all(dd.W .== 0.0)
        @test delay_variance(dd) ≈ 750.0^2 / 30
        @test delay_rate(dd) ≈ 30 / 750.0
        @test delay_total(dd) == 0.0

        # With initial value
        dd2 = DistributedDelay(5, 100.0; W0=10.0)
        @test delay_total(dd2) == 50.0

        @test_throws ArgumentError DistributedDelay(0, 100.0)
        @test_throws ArgumentError DistributedDelay(5, -1.0)
    end

    @testset "FraserGilbertResponse" begin
        fr = FraserGilbertResponse(0.5)
        @test fr.a == 0.5

        # When supply >> demand, acquisition ≈ demand
        acq = acquire(fr, 1000.0, 10.0)
        @test acq ≈ 10.0 atol=0.01

        # When supply << demand, acquisition ≈ a * supply (approximately)
        acq2 = acquire(fr, 1.0, 1000.0)
        @test acq2 < 1.0  # Can't acquire more than supply

        # Zero cases
        @test acquire(fr, 0.0, 10.0) == 0.0
        @test acquire(fr, 10.0, 0.0) == 0.0

        # Supply/demand ratio
        φ = supply_demand_ratio(fr, 100.0, 100.0)
        @test 0.0 <= φ <= 1.0

        @test_throws ArgumentError FraserGilbertResponse(-1.0)
    end

    @testset "MetabolicPool" begin
        pool = MetabolicPool(100.0,
                             [30.0, 40.0, 50.0],
                             [:respiration, :growth, :reserves])

        alloc = allocate(pool)
        @test alloc[1] == 30.0   # Full respiration
        @test alloc[2] == 40.0   # Full growth
        @test alloc[3] == 30.0   # Only 30 left for reserves

        @test supply_demand_index(pool) ≈ 100.0 / 120.0

        # Deficit scenario
        pool2 = MetabolicPool(20.0,
                              [30.0, 40.0],
                              [:resp, :growth])
        alloc2 = allocate(pool2)
        @test alloc2[1] == 20.0
        @test alloc2[2] == 0.0
    end

    @testset "Q10Respiration" begin
        resp = Q10Respiration(0.03, 2.3, 25.0)
        @test respiration_rate(resp, 25.0) ≈ 0.03  # At reference
        @test respiration_rate(resp, 35.0) ≈ 0.03 * 2.3  # 10°C above
        @test respiration_rate(resp, 15.0) < 0.03  # Below reference
    end

    @testset "MP vs BDF Abstractions" begin
        dev = LinearDevelopmentRate(10.0, 35.0)
        acq = FraserGilbertResponse(0.5)
        resp = Q10Respiration(0.03, 2.3, 25.0)
        pool = MetabolicPool(100.0, [40.0, 20.0], [:growth, :reserves])

        bdf = BiodemographicFunctions(dev, acq, resp; label=:cotton_bdf)
        hybrid = CoupledPBDMModel(bdf, pool; label=:cotton_hybrid)

        @test dev isa AbstractBiodemographicFunction
        @test acq isa AbstractBiodemographicFunction
        @test resp isa AbstractBiodemographicFunction
        @test pool isa AbstractAllocationModel
        @test bdf isa AbstractBiodemographicModel
        @test hybrid isa AbstractHybridPBDMApproach

        @test approach_family(dev) == :bdf
        @test approach_family(pool) == :mp
        @test approach_family(bdf) == :bdf
        @test approach_family(hybrid) == :hybrid

        @test development_component(bdf) === dev
        @test acquisition_component(bdf) === acq
        @test respiration_component(bdf) === resp
        @test biodemography(hybrid) === bdf
        @test allocation_model(hybrid) === pool
        @test hybrid.label == :cotton_hybrid
    end

    @testset "LifeStage" begin
        delay = DistributedDelay(20, 500.0; W0=1.0)
        dr = LinearDevelopmentRate(10.0, 35.0)
        stage = LifeStage(:larva, delay, dr, 0.01)
        @test stage.name == :larva
        @test stage.μ == 0.01
        @test_throws ArgumentError LifeStage(:bad, delay, dr, -0.1)
    end

    @testset "Population" begin
        dr = LinearDevelopmentRate(10.0, 35.0)
        stages = [
            LifeStage(:egg, DistributedDelay(10, 60.0; W0=100.0), dr, 0.01),
            LifeStage(:larva, DistributedDelay(20, 140.0), dr, 0.02),
            LifeStage(:pupa, DistributedDelay(10, 70.0), dr, 0.005),
            LifeStage(:adult, DistributedDelay(5, 200.0), dr, 0.05)
        ]
        pop = Population(:pest, stages)
        @test n_stages(pop) == 4
        @test n_substages(pop) == 45
        @test total_population(pop) > 0
    end
end
