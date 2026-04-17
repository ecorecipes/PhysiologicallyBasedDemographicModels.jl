@testset "Epidemiology" begin
    @testset "SIRDisease construction" begin
        d = SIRDisease(0.3, 0.1, 0.01)
        @test d.β ≈ 0.3
        @test d.γ ≈ 0.1
        @test d.μ_d ≈ 0.01

        @test_throws ArgumentError SIRDisease(-0.1, 0.1, 0.01)
        @test_throws ArgumentError SIRDisease(0.3, -0.1, 0.01)
    end

    @testset "DiseaseState" begin
        ds = DiseaseState(990.0, 10.0)
        @test total_alive(ds) ≈ 1000.0
        @test prevalence(ds) ≈ 0.01
        @test ds.R ≈ 0.0
        @test ds.D ≈ 0.0
    end

    @testset "SIR step_disease!" begin
        state = DiseaseState(990.0, 10.0)
        disease = SIRDisease(0.3, 0.1, 0.01)
        result = step_disease!(state, disease, 1.0)

        @test result.new_infections >= 0
        @test result.recoveries >= 0
        @test result.deaths >= 0
        @test total_alive(state) ≤ 1000.0  # no births
        @test state.S >= 0
        @test state.I >= 0
        @test state.R >= 0
    end

    @testset "SIR conservation" begin
        state = DiseaseState(990.0, 10.0)
        disease = SIRDisease(0.3, 0.1, 0.0)  # no disease death
        initial = total_alive(state)

        for _ in 1:100
            step_disease!(state, disease, 1.0)
        end

        @test total_alive(state) ≈ initial atol=1e-6
    end

    @testset "R0" begin
        d = SIRDisease(0.5, 0.1, 0.01)
        r0 = R0(d, 1.0)
        @test r0 ≈ 0.5 / 0.11
    end

    @testset "VectorBorneDisease" begin
        vbd = VectorBorneDisease(0.4, 0.3, 0.1, 0.01, 10.0)
        host = DiseaseState(990.0, 10.0)
        vec = VectorState(1000.0)

        # Introduce some infected vectors
        vec.I = 50.0
        vec.S = 950.0

        result = step_vector_disease!(host, vec, vbd, 0.5)
        @test result.host_infections >= 0
        @test result.vector_infections >= 0
        @test host.S >= 0
        @test vec.S >= 0
    end

    @testset "Empty populations" begin
        state = DiseaseState(0.0, 0.0)
        disease = SIRDisease(0.3, 0.1, 0.01)
        result = step_disease!(state, disease, 1.0)
        @test result.new_infections ≈ 0.0
    end
end
