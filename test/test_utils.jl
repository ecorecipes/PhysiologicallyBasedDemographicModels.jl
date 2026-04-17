@testset "Utilities" begin
    @testset "photoperiod" begin
        # Summer solstice at 45°N should be > 15 hours
        dl_summer = photoperiod(45.0, 172)  # ~June 21
        @test dl_summer > 15.0
        @test dl_summer < 24.0

        # Winter solstice at 45°N should be < 10 hours
        dl_winter = photoperiod(45.0, 355)  # ~Dec 21
        @test dl_winter < 10.0
        @test dl_winter > 0.0

        # Equator should be ~12 hours year-round
        dl_equator = photoperiod(0.0, 80)
        @test dl_equator ≈ 12.0 atol=0.5
    end

    @testset "degree_days_sine" begin
        # All above threshold
        dd = degree_days_sine(15.0, 25.0, 10.0)
        @test dd ≈ 10.0  # Mean(15,25)=20 - 10 = 10

        # All below threshold
        dd_cold = degree_days_sine(5.0, 8.0, 10.0)
        @test dd_cold == 0.0

        # Partial day
        dd_partial = degree_days_sine(8.0, 22.0, 10.0)
        @test dd_partial > 0.0
        @test dd_partial < 10.0  # Less than full contribution

        # With upper threshold
        dd_capped = degree_days_sine(20.0, 40.0, 10.0; T_upper=30.0)
        @test dd_capped < degree_days_sine(20.0, 40.0, 10.0)
    end

    @testset "make_population" begin
        dr = LinearDevelopmentRate(10.0, 35.0)
        pop = make_population(:moth,
            [(:egg, 10, 80.0), (:larva, 20, 200.0), (:pupa, 10, 70.0)],
            dr; mortality=0.01)

        @test pop.name == :moth
        @test n_stages(pop) == 3
        @test pop.stages[1].name == :egg
        @test pop.stages[1].delay.k == 10
        @test pop.stages[1].μ == 0.01
        @test pop.stages[2].delay.τ == 200.0
    end
end
