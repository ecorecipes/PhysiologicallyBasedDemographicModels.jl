#!/usr/bin/env julia

include(joinpath(@__DIR__, "runtime_validation_utils.jl"))

println("Running runtime validation for 36_apple_tree.qmd")
info = run_vignette_runtime("36_apple_tree.qmd"; min_figures=0)

function simulate_apple(; temp_offset::Float64=0.0, fruit_factor::Float64=1.0)
    weather_days = if iszero(temp_offset)
        Main.weather.days
    else
        [DailyWeather(Main.temps_mean[d] + temp_offset,
                      Main.temps_min[d] + temp_offset,
                      Main.temps_max[d] + temp_offset;
                      radiation=Main.rads[d])
         for d in 1:Main.n_days]
    end
    weather = WeatherSeries(weather_days; day_offset=1)

    apple = Population(:golden_delicious, [
        LifeStage(:leaf, DistributedDelay(Main.k, 700.0; W0=0.3), Main.dev_rate, Main.DEMAND_LEAF),
        LifeStage(:shoot, DistributedDelay(Main.k, 1200.0; W0=0.15), Main.dev_rate, Main.DEMAND_SHOOT),
        LifeStage(:fruit, DistributedDelay(Main.k, 1100.0; W0=0.0), Main.dev_rate, Main.DEMAND_FRUIT * fruit_factor),
        LifeStage(:root, DistributedDelay(Main.k, 200.0; W0=0.05), Main.dev_rate, Main.DEMAND_ROOT),
    ])

    prob = PBDMProblem(Main.apple_hybrid, apple, weather, (1, Main.n_days))
    sol = solve(prob, DirectIteration())
    return (;
        sol,
        cdd = cumulative_degree_days(sol)[end],
        λ = net_growth_rate(sol),
        fruit = stage_trajectory(sol, 3)[end],
    )
end

baseline = simulate_apple()
cool = simulate_apple(temp_offset=-2.0)
warm = simulate_apple(temp_offset=2.0)
light = simulate_apple(fruit_factor=0.5)
heavy = simulate_apple(fruit_factor=1.5)

@testset "Apple tree runtime validation" begin
    @test info.fig_count >= 0
    @test baseline.cdd > 2000
    @test maximum(stage_trajectory(Main.sol, 1)) > 0
    @test cool.cdd < baseline.cdd < warm.cdd
    @test cool.λ > baseline.λ > warm.λ
    @test light.λ > heavy.λ
    @test isfinite(baseline.fruit)
end

println("Validated apple tree runtime ($(info.fig_count) figures, baseline λ=$(round(baseline.λ, digits=4)))")
