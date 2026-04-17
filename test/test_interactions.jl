@testset "Interactions" begin
    @testset "HollingTypeII" begin
        fr = HollingTypeII(0.5, 0.1)
        @test functional_response(fr, 0.0) ≈ 0.0
        # At low density, approximately linear: a * N
        @test functional_response(fr, 0.01) ≈ 0.5 * 0.01 atol=1e-4
        # Asymptotes at 1/h = 10
        @test functional_response(fr, 1e6) < 1.0 / fr.h + 0.01

        @test_throws ArgumentError HollingTypeII(0.0, 0.1)
        @test_throws ArgumentError HollingTypeII(0.5, -0.1)
    end

    @testset "HollingTypeIII" begin
        fr = HollingTypeIII(0.5, 10.0, 0.1)
        @test functional_response(fr, 0.0) ≈ 0.0
        # Type III is sigmoid: starts slowly at low density
        low = functional_response(fr, 1.0)
        mid = functional_response(fr, 10.0)
        high = functional_response(fr, 100.0)
        @test low < mid < high
    end

    @testset "FraserGilbert via functional_response" begin
        fg = FraserGilbertResponse(1.0)
        acq = functional_response(fg, 100.0, 50.0)
        @test acq > 0
        @test acq ≤ 50.0
    end

    @testset "TrophicWeb" begin
        web = TrophicWeb()
        link1 = TrophicLink(:predator, :prey, HollingTypeII(0.5, 0.1), 0.3)
        link2 = TrophicLink(:predator, :prey2, HollingTypeII(0.3, 0.2), 0.2)
        add_link!(web, link1)
        add_link!(web, link2)

        pred_links = find_links(web, :predator)
        @test length(pred_links) == 2

        prey_predators = find_predators(web, :prey)
        @test length(prey_predators) == 1
        @test prey_predators[1].predator_name == :predator

        empty_links = find_links(web, :herbivore)
        @test length(empty_links) == 0
    end

    @testset "SIT mating" begin
        sit = SITRelease(10000.0, 0.8, 7)

        # On release day
        frac = fertile_mating_fraction(1000.0, sit, 7)
        @test frac < 1.0
        @test frac > 0.0

        # Without SIT (effectively zero steriles in distant future)
        frac_no = fertile_mating_fraction(1000.0, SITRelease(0.0, 0.8, 7), 0)
        @test frac_no ≈ 1.0

        # Zero wild males
        frac_zero = fertile_mating_fraction(0.0, sit, 7)
        @test frac_zero ≈ 0.0
    end
end
