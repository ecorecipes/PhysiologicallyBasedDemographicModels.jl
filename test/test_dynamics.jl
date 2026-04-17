@testset "Dynamics" begin
    @testset "step_delay! basic" begin
        # Simple delay: k=3, τ=30 DD, initially W=[10,10,10]
        dd = DistributedDelay(3, 30.0; W0=10.0)
        @test delay_total(dd) == 30.0

        # Step with 10 degree-days and 5 units inflow
        result = step_delay!(dd, 10.0, 5.0)

        # Should have outflow from last substage and some attrition=0
        @test result.outflow >= 0.0
        @test result.attrition == 0.0
        @test all(dd.W .>= 0.0)
    end

    @testset "step_delay! with mortality" begin
        dd = DistributedDelay(5, 50.0; W0=100.0)
        total_before = delay_total(dd)

        result = step_delay!(dd, 5.0, 0.0; μ=0.01)

        total_after = delay_total(dd)
        # Total should decrease (mortality + outflow, no inflow)
        @test total_after < total_before
        @test result.attrition > 0.0
    end

    @testset "step_delay! mass conservation (no mortality)" begin
        dd = DistributedDelay(10, 100.0; W0=0.0)
        total_in = 0.0
        total_out = 0.0

        for _ in 1:200
            inflow = 1.0
            total_in += inflow
            result = step_delay!(dd, 5.0, inflow)
            total_out += result.outflow
        end

        # Total in should equal total in delay + total out
        @test total_in ≈ delay_total(dd) + total_out atol=1e-6
    end

    @testset "step_delay! remains conservative for large daily DD" begin
        dd = DistributedDelay(1, 50.0; W0=1.0)
        result = step_delay!(dd, 100.0, 0.0)
        @test all(dd.W .>= 0.0)
        @test delay_total(dd) + result.outflow + result.attrition ≈ 1.0 atol=1e-6
    end

    @testset "step_population!" begin
        dr = LinearDevelopmentRate(10.0, 35.0)
        stages = [
            LifeStage(:egg, DistributedDelay(10, 60.0; W0=100.0), dr, 0.001),
            LifeStage(:larva, DistributedDelay(15, 150.0), dr, 0.002),
            LifeStage(:pupa, DistributedDelay(10, 70.0), dr, 0.001),
        ]
        pop = Population(:insect, stages)

        w = DailyWeather(25.0)
        result = step_population!(pop, w)

        @test result.degree_days == 15.0  # 25 - 10
        @test result.stage_degree_days == [15.0, 15.0, 15.0]
        @test result.maturation >= 0.0
        @test result.total_attrition >= 0.0
        @test length(result.stage_totals) == 3

        # Population should still be mostly in egg stage after one day
        @test result.stage_totals[1] > result.stage_totals[2]
    end

    @testset "step_population! with inflow and stress" begin
        dr = LinearDevelopmentRate(10.0, 35.0)
        stages = [
            LifeStage(:egg, DistributedDelay(5, 50.0; W0=10.0), dr, 0.0),
            LifeStage(:larva, DistributedDelay(5, 100.0), dr, 0.0),
        ]
        pop = Population(:test, stages)

        result = step_population!(pop, DailyWeather(25.0); inflow=5.0, stage_stress=[0.1, 0.2])
        @test length(result.stage_stress) == 2
        @test result.stage_stress[1] == 0.1
        @test result.total_attrition > 0.0
        @test_throws ArgumentError step_population!(pop, DailyWeather(25.0); stage_stress=[1.2, 0.2])
    end

    @testset "step_population! records stage-specific degree-days" begin
        stages = [
            LifeStage(:juvenile, DistributedDelay(5, 50.0; W0=10.0), LinearDevelopmentRate(10.0, 35.0), 0.0),
            LifeStage(:adult, DistributedDelay(5, 50.0; W0=0.0), LinearDevelopmentRate(15.0, 35.0), 0.0),
        ]
        pop = Population(:heterogeneous, stages)
        result = step_population!(pop, DailyWeather(20.0))
        @test result.degree_days == 10.0
        @test result.stage_degree_days == [10.0, 5.0]
    end

    @testset "step_system! approach-aware dispatch" begin
        dr = LinearDevelopmentRate(10.0, 35.0)
        acq = FraserGilbertResponse(0.6)
        resp = Q10Respiration(0.02, 2.0, 25.0)
        pool = MetabolicPool(15.0, [1.0, 0.5], [:egg, :larva])
        bdf = BiodemographicFunctions(dr, acq, resp)
        hybrid = CoupledPBDMModel(bdf, pool)

        stages = [
            LifeStage(:egg, DistributedDelay(10, 60.0; W0=50.0), dr, 0.001),
            LifeStage(:larva, DistributedDelay(10, 120.0), dr, 0.002),
        ]
        pop_legacy = Population(:insect_legacy, deepcopy(stages))
        pop_bdf = Population(:insect_bdf, deepcopy(stages))
        pop_mp = Population(:insect_mp, deepcopy(stages))
        pop_hybrid = Population(:insect_hybrid, deepcopy(stages))
        w = DailyWeather(25.0, 20.0, 30.0; radiation=18.0)

        legacy = step_system!(pop_legacy, w)
        bdf_step = step_system!(pop_bdf, w, bdf)
        mp_step = step_system!(pop_mp, w, pool)
        hybrid_step = step_system!(pop_hybrid, w, hybrid)

        @test legacy.approach_family == :legacy
        @test bdf_step.approach_family == :bdf
        @test mp_step.approach_family == :mp
        @test hybrid_step.approach_family == :hybrid

        @test haskey(pairs(bdf_step), :gross_supply)
        @test haskey(pairs(mp_step), :allocations)
        @test haskey(pairs(hybrid_step), :allocations)
        @test 0.0 <= bdf_step.supply_demand <= 1.0
        @test bdf_step.supply_demand_ratio == bdf_step.supply_demand
        @test 0.0 <= mp_step.supply_demand <= 1.0
        @test mp_step.supply_demand_ratio == mp_step.supply_demand
        @test 0.0 <= hybrid_step.supply_demand <= 1.0
        @test hybrid_step.supply_demand_ratio == hybrid_step.supply_demand
        @test length(hybrid_step.stage_stress) == 2
    end

    @testset "step_system! rejects unspecialized abstract approaches" begin
        struct DummyBDF <: AbstractBiodemographicApproach end
        struct DummyMP <: AbstractAllocationApproach end

        dr = LinearDevelopmentRate(10.0, 35.0)
        pop = Population(:test, [LifeStage(:egg, DistributedDelay(5, 50.0; W0=10.0), dr, 0.0)])
        w = DailyWeather(20.0, 15.0, 25.0; radiation=10.0)

        @test_throws ArgumentError step_system!(pop, w, DummyBDF())
        @test_throws ArgumentError step_system!(pop, w, DummyMP())
    end
end
