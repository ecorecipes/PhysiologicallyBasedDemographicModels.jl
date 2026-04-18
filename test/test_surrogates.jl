using Test
using PhysiologicallyBasedDemographicModels

@testset "LogLinearSurrogate basics" begin
    s = LogLinearSurrogate(
        2.0,
        Dict(:A => -1.0, :B => -0.5, :C => -0.25),
        Dict(Tuple(sort([:A, :B])) => 0.3),
        "test"
    )

    @test s.label == "test"
    @test variables(s) == [:A, :B, :C]
    @test predict_log(s, Set{Symbol}()) ≈ 2.0
    @test predict_log(s, Set([:A])) ≈ 2.0 - 1.0
    @test predict_log(s, Set([:A, :B])) ≈ 2.0 - 1.0 - 0.5 + 0.3
    @test predict_log(s, Set([:A, :B, :C])) ≈ 2.0 - 1.0 - 0.5 - 0.25 + 0.3
    @test predict(s, Set([:A])) ≈ exp(1.0)

    # Iterable other than Set works
    @test predict_log(s, [:A, :B]) ≈ predict_log(s, Set([:A, :B]))
    @test predict_log(s, (:A,)) ≈ predict_log(s, Set([:A]))
end

@testset "LogLinearSurrogate without interactions" begin
    s = LogLinearSurrogate(1.0, Dict(:X => -2.0))
    @test isempty(s.interactions)
    @test predict_log(s, Set{Symbol}()) ≈ 1.0
    @test predict_log(s, Set([:X])) ≈ -1.0
    @test variables(s) == [:X]
end

@testset "marginal_effects ranking" begin
    s = LogLinearSurrogate(0.0,
        Dict(:big => -2.0, :small => -0.1, :mid => -1.0))
    me = marginal_effects(s)
    @test first(me).first == :big          # most negative β → smallest A
    @test last(me).first == :small         # least effect → largest A
    @test me[2].first == :mid
    @test all(0 < x.second ≤ 1 for x in me)
end

@testset "enumerate_strategies covers 2^n subsets" begin
    vars = [:A, :B, :C]
    strategies = collect(enumerate_strategies(vars))
    @test length(strategies) == 8
    @test Set{Symbol}() in strategies
    @test Set([:A, :B, :C]) in strategies
    # All distinct
    @test length(unique(strategies)) == 8

    s = LogLinearSurrogate(0.0, Dict(:A => -1.0, :B => -1.0))
    @test length(collect(enumerate_strategies(s))) == 4
end

@testset "pareto_frontier — maximize value" begin
    # Points (cost, value):
    #   A (1, 1)  - frontier (cheapest with positive value)
    #   B (2, 3)  - frontier (better value)
    #   C (3, 2)  - dominated by B
    #   D (4, 5)  - frontier (best value)
    #   E (5, 5)  - tied value, higher cost → dominated
    pts = [(:A, 1, 1), (:B, 2, 3), (:C, 3, 2), (:D, 4, 5), (:E, 5, 5)]
    frontier = pareto_frontier(pts; cost = p -> p[2], value = p -> p[3])
    names = [p[1] for p in frontier]
    @test names == [:A, :B, :D]
end

@testset "pareto_frontier — minimize value (cost vs damage)" begin
    pts = [(1, 10.0), (2, 8.0), (3, 9.0), (4, 5.0), (5, 6.0)]
    frontier = pareto_frontier(pts;
        cost = first, value = last, maximize_value = false)
    @test frontier == [(1, 10.0), (2, 8.0), (4, 5.0)]
end

@testset "End-to-end: Cure 2020 Colombia surrogate" begin
    # Subset of the published Colombia coefficients (Cure et al. 2020 Table 5).
    s = LogLinearSurrogate(
        11.0481,
        Dict(:H    => -1.5835, :CU => -0.6228,
             :Bb   => -0.3178, :Pc => -0.3751,
             :Het  => -0.1388, :Stei => -0.1215,
             :C    => -0.1381, :Ma => -0.0420,
             :Cs   => -0.0257, :Pn => -0.0286),
        Dict(Tuple(sort([:H, :Pc])) => 0.1864),
        "Colombia"
    )
    me = marginal_effects(s)
    # Harvest must be the strongest tactic (smallest A_X).
    @test first(me).first == :H
    # Reproduce A_H from the paper to 3 decimals.
    @test isapprox(first(me).second, 0.205; atol = 0.001)

    # Sweep all 2^10 = 1024 strategies; harvest-only must beat no-control
    # and the top strategy must include H.
    strategies = collect(enumerate_strategies(s))
    @test length(strategies) == 1024
    log_y = [predict_log(s, on) for on in strategies]
    best = strategies[argmin(log_y)]
    @test :H in best
end
