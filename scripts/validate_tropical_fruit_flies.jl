#!/usr/bin/env julia

include(joinpath(@__DIR__, "runtime_validation_utils.jl"))
using Statistics

println("Running runtime validation for 39_tropical_fruit_flies.qmd")
info = run_vignette_runtime("39_tropical_fruit_flies.qmd"; min_figures=0)

mean_current = mean(d.T_mean for d in Main.days_current)
mean_rcp85 = mean(d.T_mean for d in Main.days_rcp85)

@testset "Tropical fruit flies runtime validation" begin
    @test info.fig_count >= 0
    @test length(Main.results) == length(Main.all_species) == 4
    @test mean_rcp85 > mean_current
    for sp in Main.all_species
        current = Main.results[sp.abbrev][:current]
        warm = Main.results[sp.abbrev][:rcp85]
        @test length(current.traj_pupae) == Main.N_DAYS
        @test length(warm.traj_pupae) == Main.N_DAYS
        @test all(isfinite, current.traj_pupae)
        @test all(isfinite, warm.traj_pupae)
        @test sum(warm.traj_pupae) >= sum(current.traj_pupae)
        @test sum(Main.fi_rcp85[sp.abbrev] .> 0.5) >= sum(Main.fi_current[sp.abbrev] .> 0.5)
    end
end

println("Validated tropical fruit fly runtime ($(info.fig_count) figures)")
