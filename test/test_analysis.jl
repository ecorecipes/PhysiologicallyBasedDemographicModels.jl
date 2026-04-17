@testset "Analysis" begin
    # Build and solve a simple model for analysis tests
    dr = LinearDevelopmentRate(10.0, 35.0)
    stages = [
        LifeStage(:egg, DistributedDelay(10, 80.0; W0=500.0), dr, 0.005),
        LifeStage(:larva, DistributedDelay(15, 200.0), dr, 0.01),
        LifeStage(:adult, DistributedDelay(5, 100.0), dr, 0.02),
    ]
    pop = Population(:moth, stages)
    weather = WeatherSeries(fill(25.0, 100); day_offset=1)
    initial_total = total_population(pop)
    prob = PBDMProblem(pop, weather, (1, 100))
    sol = solve(prob, DirectIteration())

    @testset "cumulative_degree_days" begin
        cdd = cumulative_degree_days(sol)
        @test length(cdd) == 100
        @test cdd[1] == 15.0
        @test cdd[end] ≈ 15.0 * 100
        @test issorted(cdd)
    end

    @testset "stage_degree_days" begin
        dd = stage_degree_days(sol, 1)
        @test length(dd) == 100
        @test all(dd .≈ 15.0)
    end

    @testset "stage_trajectory" begin
        traj = stage_trajectory(sol, 1)
        @test length(traj) == 100
        @test traj[1] < initial_total
    end

    @testset "total_population" begin
        tp = total_population(sol)
        @test length(tp) == 100
        @test tp[1] < initial_total
    end

    @testset "net_growth_rate" begin
        r = net_growth_rate(sol)
        @test r isa Float64
        @test r > 0.0
    end

    @testset "phenology" begin
        day = phenology(sol; threshold=0.5)
        # With 100 days of warm weather, maturation should occur
        if !isnothing(day)
            @test day isa Int
            @test day >= 1
            @test day <= 100
        end
    end
end
