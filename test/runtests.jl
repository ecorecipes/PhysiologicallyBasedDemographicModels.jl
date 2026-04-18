using PhysiologicallyBasedDemographicModels
using Test
using LinearAlgebra
using Statistics

@testset "PhysiologicallyBasedDemographicModels.jl" begin
    include("test_types.jl")
    include("test_weather.jl")
    include("test_dynamics.jl")
    include("test_solve.jl")
    include("test_analysis.jl")
    include("test_utils.jl")
    include("test_temperature_responses.jl")
    include("test_interactions.jl")
    include("test_economics.jl")
    include("test_genetics.jl")
    include("test_epidemiology.jl")
    include("test_coupled.jl")
    include("test_scenarios.jl")
    include("test_vignette_outputs.jl")
    include("test_theory.jl")
    include("test_continuous.jl")
    include("test_surrogates.jl")
end
