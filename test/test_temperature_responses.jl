using Test
using PhysiologicallyBasedDemographicModels

@testset "Temperature responses" begin

    @testset "triangular_thermal_scalar" begin
        # zero outside [θL, θU]
        @test triangular_thermal_scalar(5.0; θL = 10.0, θU = 30.0) == 0
        @test triangular_thermal_scalar(35.0; θL = 10.0, θU = 30.0) == 0
        # peak at midpoint by default = 1
        @test triangular_thermal_scalar(20.0; θL = 10.0, θU = 30.0) ≈ 1.0
        # ASCII alias
        @test phiT(20.0; θL = 10.0, θU = 30.0) == triangular_thermal_scalar(20.0; θL = 10.0, θU = 30.0)
        # custom Topt: linear interp on each side
        f = triangular_thermal_scalar
        @test f(15.0; θL = 10.0, θU = 30.0, Topt = 25.0) ≈ (15-10)/(25-10)
        @test f(28.0; θL = 10.0, θU = 30.0, Topt = 25.0) ≈ (30-28)/(30-25)
    end

    @testset "briere_rate / fecundity_briere" begin
        # zero outside [θL, θU]
        @test briere_rate(5.0;  a = 1e-4, θL = 10.0, θU = 35.0) == 0
        @test briere_rate(40.0; a = 1e-4, θL = 10.0, θU = 35.0) == 0
        # positive in interior
        @test briere_rate(25.0; a = 1e-4, θL = 10.0, θU = 35.0) > 0
        # alias matches
        @test fecundity_briere(25.0; a = 1e-4, θL = 10.0, θU = 35.0) ==
              briere_rate(25.0; a = 1e-4, θL = 10.0, θU = 35.0)
        # exponent m = 2 → sqrt form: a*T*(T-θL)*sqrt(θU-T)
        T = 25.0
        @test briere_rate(T; a = 1.0, θL = 10.0, θU = 35.0, m = 2) ≈
              T * (T - 10.0) * sqrt(35.0 - T)
    end

    @testset "fecundity_gaussian" begin
        @test fecundity_gaussian(25.0; F_max = 10.0, T_opt = 25.0, T_range = 5.0) ≈ 10.0
        @test fecundity_gaussian(20.0; F_max = 10.0, T_opt = 25.0, T_range = 5.0) ≈
              10 * exp(-1)
    end

    @testset "daily_mortality_quadratic" begin
        @test daily_mortality_quadratic(25.0; a = 0.001, T_opt = 25.0, μ_min = 0.01) ≈ 0.01
        @test daily_mortality_quadratic(40.0; a = 0.001, T_opt = 25.0, μ_min = 0.01) ≈
              0.001*15^2 + 0.01
        # clamps to μ_max
        @test daily_mortality_quadratic(100.0; a = 1.0, T_opt = 25.0, μ_max = 0.5) == 0.5
        # never negative
        @test daily_mortality_quadratic(25.0; a = -1.0, T_opt = 25.0) ≥ 0.0
    end

    @testset "diapause_fraction_logistic" begin
        # at D = D50, sigmoid is 0.5
        @test diapause_fraction_logistic(13.0; D50 = 13.0, slope = 5.0) ≈ 0.5
        # short days (D < D50) with positive slope → near 1
        @test diapause_fraction_logistic(10.0; D50 = 13.0, slope = 5.0) > 0.99
        # long days (D > D50) → near 0
        @test diapause_fraction_logistic(16.0; D50 = 13.0, slope = 5.0) < 0.01
    end

    @testset "diapause_fraction_linear" begin
        # above critical day length: zero diapause
        @test diapause_fraction_linear(15.0, 20.0; D_crit = 13.5, D_comp = 12.0) == 0.0
        # below complete day length: full diapause
        @test diapause_fraction_linear(11.0, 20.0; D_crit = 13.5, D_comp = 12.0) == 1.0
        # midway between: half (clamp confirms)
        @test diapause_fraction_linear(12.75, 20.0;
                                       D_crit = 13.5, D_comp = 12.0) ≈ 0.5
        # cold modulation amplifies diapause fraction
        cold = diapause_fraction_linear(12.75, 5.0;
                                        D_crit = 13.5, D_comp = 12.0,
                                        T_mod = 10.0, k = 0.1)
        @test cold > 0.5
        @test cold ≤ 1.0
    end

    @testset "gilbert_fraser_attack" begin
        # zero hosts → zero attack
        @test gilbert_fraser_attack(1.0, 5.0, 0.0) == 0.0
        # at α*D/H = 1, expect 1 - exp(-1) ≈ 0.632
        @test gilbert_fraser_attack(1.0, 1.0, 1.0) ≈ 1 - exp(-1)
        # increasing α monotonically increases attack
        @test gilbert_fraser_attack(2.0, 1.0, 1.0) >
              gilbert_fraser_attack(1.0, 1.0, 1.0)
        # bounded above by 1
        @test gilbert_fraser_attack(100.0, 100.0, 1.0) ≤ 1.0
        @test gilbert_fraser_attack(100.0, 100.0, 1.0) > 0.999
    end

end
