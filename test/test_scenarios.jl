using Test
using PhysiologicallyBasedDemographicModels

@testset "scenarios — run_scenarios and compare_metrics" begin
    # Trivial simulator
    simulator = s -> (; yield = s.dose * 10.0, pest = 100.0 - s.dose * 5.0)

    @testset "Dict of scenarios" begin
        scenarios = Dict("low" => (; dose = 1.0),
                         "high" => (; dose = 5.0))
        results = run_scenarios(simulator, scenarios)
        @test results isa AbstractDict
        @test length(results) == 2
        @test results["low"].yield ≈ 10.0
        @test results["high"].yield ≈ 50.0
    end

    @testset "Vector of Pair scenarios" begin
        scenarios = [:a => (; dose = 2.0), :b => (; dose = 3.0)]
        results = run_scenarios(simulator, scenarios)
        @test results[:a].yield ≈ 20.0
        @test results[:b].yield ≈ 30.0
    end

    @testset "post hook adds derived metrics" begin
        scenarios = Dict("s1" => (; dose = 2.0))
        results = run_scenarios(simulator, scenarios;
            post = r -> merge(r, (; profit = r.yield - 5.0)))
        @test results["s1"].profit ≈ 15.0
        @test results["s1"].yield ≈ 20.0
    end

    @testset "compare_metrics — symbol list" begin
        scenarios = Dict("low" => (; dose = 1.0), "high" => (; dose = 5.0))
        results = run_scenarios(simulator, scenarios)
        rows = compare_metrics(results, [:yield, :pest])
        @test length(rows) == 2
        @test all(haskey(r, :scenario) for r in rows)
        @test all(haskey(r, :yield) for r in rows)
        @test all(haskey(r, :pest) for r in rows)
        low_row = first(filter(r -> r.scenario == "low", rows))
        @test low_row.yield ≈ 10.0
    end

    @testset "compare_metrics — extractor pairs" begin
        scenarios = Dict("s1" => (; dose = 4.0))
        results = run_scenarios(simulator, scenarios)
        rows = compare_metrics(results,
            [:yield => r -> r.yield,
             :net   => r -> r.yield - 2 * r.pest])
        @test length(rows) == 1
        @test rows[1].yield ≈ 40.0
        @test rows[1].net ≈ 40.0 - 2 * 80.0
    end

    @testset "empty results" begin
        results = run_scenarios(simulator, Dict{String, Any}())
        @test isempty(results)
        @test compare_metrics(results, [:yield]) == NamedTuple[]
    end
end
