#!/usr/bin/env julia

include(joinpath(@__DIR__, "runtime_validation_utils.jl"))

println("Running runtime validation for 38_whitefly_autoparasitoid.qmd")
info = run_vignette_runtime("38_whitefly_autoparasitoid.qmd"; min_figures=0)

@testset "Whitefly autoparasitoid runtime validation" begin
    @test info.fig_count >= 0
    @test length(Main.results) == length(Main.scenarios) == 7
    @test length(Main.res_all.traj_wf) == Main.N_DAYS + 1
    @test length(Main.res_all.traj_sr_p2) == Main.N_DAYS
    @test all(x -> isfinite(x) && x >= 0, Main.res_all.traj_wf)
    @test all(x -> 0.0 <= x <= 1.0, Main.res_all.traj_sr_p2)
    @test all(x -> 0.0 <= x <= 1.0, Main.res_all.traj_sr_p3)
    for result in values(Main.results)
        @test length(result.traj_wf) == Main.N_DAYS + 1
        @test all(isfinite, result.traj_wf)
        @test all(isfinite, result.traj_p1)
        @test all(isfinite, result.traj_p2)
        @test all(isfinite, result.traj_p3)
    end
end

println("Validated whitefly-autoparasitoid runtime ($(info.fig_count) figures)")
