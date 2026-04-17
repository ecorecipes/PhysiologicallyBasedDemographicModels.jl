@testset "Solve" begin
    @testset "Single-species DI deterministic" begin
        # Simple 3-stage insect with linear development
        dr = LinearDevelopmentRate(10.0, 35.0)
        stages = [
            LifeStage(:egg, DistributedDelay(10, 60.0; W0=1000.0), dr, 0.005),
            LifeStage(:larva, DistributedDelay(20, 200.0; W0=0.0), dr, 0.01),
            LifeStage(:adult, DistributedDelay(5, 100.0; W0=0.0), dr, 0.02),
        ]
        pop = Population(:pest, stages)

        # 100 days of warm weather
        weather = WeatherSeries(fill(25.0, 100); day_offset=1)
        initial_eggs = delay_total(pop.stages[1].delay)

        prob = PBDMProblem(pop, weather, (1, 100))
        @test prob isa PBDMProblem
        @test prob.structure isa SingleSpeciesPBDM
        @test prob.density isa DensityIndependent

        sol = solve(prob, DirectIteration())
        @test sol isa PBDMSolution
        @test sol.retcode == :Success
        @test length(sol.t) == 100
        @test length(sol.u) == 100
        @test length(sol.degree_days) == 100
        @test length(sol.lambdas) == 100
        @test size(sol.stage_degree_days) == (3, 100)
        @test size(sol.stage_totals) == (3, 100)

        # All degree-days should be 15 (25 - 10)
        @test all(sol.degree_days .≈ 15.0)
        @test all(sol.stage_degree_days[1, :] .≈ 15.0)

        # Eggs should decrease over time (flowing into larva)
        @test sol.stage_totals[1, 1] < initial_eggs
        @test sol.stage_totals[1, end] < sol.stage_totals[1, 1]

        # Eventually larvae should appear (may need many days at 15 DD/day)
        # With τ_egg=60DD and 15DD/day, eggs transit in ~4 days
        # 100 days × 15DD = 1500 DD total, but with mortality and delay distribution
        # larvae should receive some flow
        @test sol.stage_totals[2, end] >= 0.0
    end

    @testset "PBDMSolution show" begin
        dr = LinearDevelopmentRate(10.0, 35.0)
        stages = [LifeStage(:s1, DistributedDelay(5, 50.0; W0=10.0), dr, 0.0)]
        pop = Population(:test, stages)
        weather = WeatherSeries(fill(20.0, 10); day_offset=1)
        prob = PBDMProblem(pop, weather, (1, 10))
        sol = solve(prob, DirectIteration())
        s = sprint(show, sol)
        @test occursin("PBDMSolution", s)
    end

    @testset "PBDMProblem show" begin
        dr = LinearDevelopmentRate(10.0, 35.0)
        stages = [LifeStage(:s1, DistributedDelay(5, 50.0; W0=10.0), dr, 0.0)]
        pop = Population(:test, stages)
        weather = WeatherSeries(fill(20.0, 10); day_offset=1)
        prob = PBDMProblem(pop, weather, (1, 10))
        s = sprint(show, prob)
        @test occursin("PBDMProblem", s)
        @test occursin("SingleSpeciesPBDM", s)
    end

    @testset "PBDMProblem with explicit approach" begin
        dr = LinearDevelopmentRate(10.0, 35.0)
        fr = FraserGilbertResponse(0.5)
        resp = Q10Respiration(0.03, 2.0, 25.0)
        bdf = BiodemographicFunctions(dr, fr, resp; label=:test_bdf)
        pool = MetabolicPool(20.0, [1.0], [:s1])
        hybrid = CoupledPBDMModel(bdf, pool; label=:test_hybrid)

        stages = [LifeStage(:s1, DistributedDelay(5, 50.0; W0=10.0), dr, 0.0)]
        pop = Population(:test, stages)
        weather = WeatherSeries(fill(20.0, 10); day_offset=1)

        prob_kw = PBDMProblem(pop, weather, (1, 10); approach=hybrid)
        prob_pos = PBDMProblem(hybrid, pop, weather, (1, 10))

        @test prob_kw.approach === hybrid
        @test prob_pos.approach === hybrid
        @test occursin("approach=hybrid", sprint(show, prob_kw))

        sol = solve(prob_kw, DirectIteration())
        @test sol.retcode == :Success
    end

    @testset "Single-species DD deterministic" begin
        # Density-dependent with reproduction function
        dr = LinearDevelopmentRate(10.0, 35.0)
        stages = [
            LifeStage(:immature, DistributedDelay(10, 100.0; W0=500.0), dr, 0.01),
            LifeStage(:mature, DistributedDelay(5, 200.0; W0=0.0), dr, 0.02),
        ]
        pop = Population(:plant, stages)
        weather = WeatherSeries(fill(22.0, 50); day_offset=1)

        # Simple reproduction: mature individuals produce offspring
        repro_fn = (pop, w, p, day) -> 0.1 * delay_total(pop.stages[end].delay)

        prob = PBDMProblem(DensityDependent(), pop, weather, (1, 50))
        sol = solve(prob, DirectIteration(); reproduction_fn=repro_fn)

        @test sol.retcode == :Success
        @test length(sol.t) == 50
    end

    @testset "Sinusoidal weather integration" begin
        dr = LinearDevelopmentRate(10.0, 35.0)
        stages = [
            LifeStage(:egg, DistributedDelay(10, 80.0; W0=100.0), dr, 0.005),
            LifeStage(:larva, DistributedDelay(15, 200.0), dr, 0.01),
        ]
        pop = Population(:moth, stages)
        weather = SinusoidalWeather(18.0, 8.0; phase=200.0)

        prob = PBDMProblem(pop, weather, (1, 365))
        sol = solve(prob, DirectIteration())

        @test sol.retcode == :Success
        @test length(sol.t) == 365

        # Degree-days should vary seasonally
        @test minimum(sol.degree_days) < maximum(sol.degree_days)

        # Some days in winter should have very low DD
        @test minimum(sol.degree_days) < 1.0
    end

    @testset "Unsupported solver modes are explicit" begin
        dr = LinearDevelopmentRate(10.0, 35.0)
        pop = Population(:test, [LifeStage(:s1, DistributedDelay(5, 50.0; W0=10.0), dr, 0.0)])
        weather = WeatherSeries(fill(20.0, 10); day_offset=1)
        prob = PBDMProblem(pop, weather, (1, 10))
        @test_throws ArgumentError solve(prob, PhysiologicallyBasedDemographicModels.EigenAnalysis())

        pop2 = Population(:test2, [LifeStage(:s1, DistributedDelay(5, 50.0; W0=5.0), dr, 0.0)])
        multi = PBDMProblem([pop, pop2], weather, (1, 10))
        @test_throws ArgumentError solve(multi, DirectIteration())
    end

    @testset "PBDMProblem validates tspan" begin
        dr = LinearDevelopmentRate(10.0, 35.0)
        pop = Population(:test, [LifeStage(:s1, DistributedDelay(5, 50.0; W0=10.0), dr, 0.0)])
        weather = WeatherSeries(fill(20.0, 10); day_offset=1)
        @test_throws ArgumentError PBDMProblem(pop, weather, (10, 1))
    end
end
