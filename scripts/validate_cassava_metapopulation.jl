#!/usr/bin/env julia

include(joinpath(@__DIR__, "runtime_validation_utils.jl"))

println("Running runtime validation for 40_cassava_metapopulation.qmd")
info = run_vignette_runtime("40_cassava_metapopulation.qmd"; min_figures=0)

@testset "Cassava metapopulation runtime validation" begin
    @test info.fig_count >= 0
    @test length(Main.result.mean_tuber) == Main.N_DAYS
    @test length(Main.result.mean_cm) == Main.N_DAYS
    @test length(Main.result.mean_cgm) == Main.N_DAYS
    @test Main.result.mean_tuber[end] > 0
    @test maximum(Main.result.mean_cm) > 0
    @test maximum(Main.result.mean_cgm) > 0
    @test length(Main.phi_values) == length(Main.cm_peaks) == length(Main.tuber_yields)
    @test all(isfinite, Main.cm_peaks)
    @test all(y -> y >= 0 && isfinite(y), Main.tuber_yields)
end

println("Validated cassava metapopulation runtime ($(info.fig_count) figures)")
