@testset "Genetics" begin
    @testset "DialleleicLocus" begin
        locus = DialleleicLocus(0.01, 0.0)
        @test locus.R ≈ 0.01
        @test locus.dominance ≈ 0.0

        # Boundary checks
        @test_throws ArgumentError DialleleicLocus(-0.1, 0.0)
        @test_throws ArgumentError DialleleicLocus(0.5, 1.5)
    end

    @testset "Hardy-Weinberg frequencies" begin
        freq = genotype_frequencies(DialleleicLocus(0.1))
        @test freq.SS ≈ 0.81
        @test freq.SR ≈ 0.18
        @test freq.RR ≈ 0.01

        # Sum to 1
        @test freq.SS + freq.SR + freq.RR ≈ 1.0

        # Pure susceptible
        freq0 = genotype_frequencies(0.0)
        @test freq0.SS ≈ 1.0
        @test freq0.RR ≈ 0.0
    end

    @testset "Selection step" begin
        locus = DialleleicLocus(0.01)
        # RR has higher fitness → R should increase
        fit = GenotypeFitness(0.5, 0.5, 1.0)
        R_new = selection_step!(locus, fit)
        @test R_new > 0.01

        # SS has higher fitness → R should decrease
        locus2 = DialleleicLocus(0.5)
        fit2 = GenotypeFitness(1.0, 0.5, 0.1)
        R_new2 = selection_step!(locus2, fit2)
        @test R_new2 < 0.5

        # Equal fitness → no change
        locus3 = DialleleicLocus(0.3)
        fit3 = GenotypeFitness(1.0, 1.0, 1.0)
        R_new3 = selection_step!(locus3, fit3)
        @test R_new3 ≈ 0.3
    end

    @testset "Allele frequency from counts" begin
        # 100 SS, 0 SR, 0 RR → R = 0
        @test allele_frequency_from_adults(100.0, 0.0, 0.0) ≈ 0.0
        # 0 SS, 0 SR, 100 RR → R = 1
        @test allele_frequency_from_adults(0.0, 0.0, 100.0) ≈ 1.0
        # 25 SS, 50 SR, 25 RR → R = 0.5
        @test allele_frequency_from_adults(25.0, 50.0, 25.0) ≈ 0.5
    end

    @testset "Dose-response" begin
        dr = DoseResponse(1.0, 4.0)
        @test mortality_probability(dr, 0.0) ≈ 0.0
        @test mortality_probability(dr, 1.0) ≈ 0.5
        # High dose → high mortality
        @test mortality_probability(dr, 10.0) > 0.99
    end

    @testset "Refuge dilution" begin
        @test refuge_dilution(0.5, 0.0, 0.1) ≈ 0.45
        @test refuge_dilution(0.5, 0.5, 0.5) ≈ 0.5  # No dilution when equal
    end

    @testset "Two-locus resistance" begin
        l1 = DialleleicLocus(0.1)
        l2 = DialleleicLocus(0.2)
        tlr = TwoLocusResistance(l1, l2)
        prob = probability_fully_resistant(tlr)
        @test prob ≈ 0.01 * 0.04  # RR₁ × RR₂
    end
end
