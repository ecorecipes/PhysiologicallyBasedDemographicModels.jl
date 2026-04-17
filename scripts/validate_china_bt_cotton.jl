#!/usr/bin/env julia

include(joinpath(@__DIR__, "runtime_validation_utils.jl"))

println("Running runtime validation for 41_china_bt_cotton.qmd")
info = run_vignette_runtime("41_china_bt_cotton.qmd"; min_figures=0)

@testset "China Bt cotton runtime validation" begin
    @test info.fig_count >= 0
    @test length(Main.strategies) == 6
    @test length(Main.results) == 6
    @test length(Main.profits) == 6
    @test all(isfinite, Main.profits)
    @test Main.results["S5: Bt(high), no spray"].yield_kg_ha >= Main.results["S3: Bt(low), no spray"].yield_kg_ha
    @test Main.results["S6: Bt(high) + 3 mirid"].yield_kg_ha >= 0
    @test Main.results["S1: Natural control"].yield_kg_ha >= 0
end

println("Validated China Bt cotton runtime ($(info.fig_count) figures)")
