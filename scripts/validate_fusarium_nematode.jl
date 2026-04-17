#!/usr/bin/env julia

include(joinpath(@__DIR__, "runtime_validation_utils.jl"))

println("Running runtime validation for 37_fusarium_nematode.qmd")
info = run_vignette_runtime("37_fusarium_nematode.qmd"; min_figures=0)

@testset "Fusarium-nematode runtime validation" begin
    @test info.fig_count >= 0
    @test Main.field_disease[:A].prevalence[end] > Main.field_disease[:B].prevalence[end] > Main.field_disease[:C].prevalence[end]
    @test Main.metric_results[:A].actual_yield < Main.metric_results[:B].actual_yield <= Main.metric_results[:C].actual_yield
    @test Main.metric_results[:B].profit > Main.metric_results[:A].profit
    @test Main.metric_results[:C].profit > Main.metric_results[:A].profit
    @test Main.metric_results[:A].symptom.day < Main.metric_results[:B].symptom.day
end

println("Validated fusarium-nematode runtime ($(info.fig_count) figures)")
