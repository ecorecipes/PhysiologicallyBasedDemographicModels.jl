module PhysiologicallyBasedDemographicModels

using CommonSolve
using LinearAlgebra
using StructuredPopulationCore
using Statistics

# --- Layer 1: Core types ---
include("types.jl")
export AbstractPBDMStructure, SingleSpeciesPBDM, MultiSpeciesPBDM
# Re-export shared types from StructuredPopulationCore
export AbstractProjectionStructure
export AbstractDensityDependence, DensityIndependent, DensityDependent
export AbstractStochasticity, Deterministic, StochasticKernelResampled, StochasticParameterResampled
export DirectIteration, EigenAnalysis

# Development rates
export AbstractDevelopmentRate
export LinearDevelopmentRate, BriereDevelopmentRate, LoganDevelopmentRate
export development_rate, degree_days

# Distributed delay
export DistributedDelay
export delay_variance, delay_rate, delay_total

# Functional response
export AbstractFunctionalResponse, FraserGilbertResponse
export acquire, supply_demand_ratio

# Metabolic pool
export MetabolicPool, allocate, supply_demand_index

# Respiration
export Q10Respiration, respiration_rate

# Life stages and populations
export LifeStage, Population
export n_stages, n_substages, total_population

# --- Layer 2: Utilities ---
include("utils.jl")
export photoperiod, degree_days_sine
export make_population

# --- Layer 4: Weather/forcing ---
include("weather.jl")
export AbstractWeather, DailyWeather, WeatherSeries, SinusoidalWeather
export get_weather

# --- Layer 5: Dynamics ---
include("dynamics.jl")
export step_delay!, step_population!

# --- Layer 6: Problem/Solution ---
include("problems.jl")
export PBDMProblem, remake

include("solve.jl")
using CommonSolve: solve
export solve, PBDMSolution

# --- Layer 7: Analysis ---
include("analysis.jl")
# Re-export from StructuredPopulationCore
export AbstractProjectionSolution
export cumulative_degree_days, stage_trajectory
export net_growth_rate, phenology

# --- Layer 8: Multi-species interactions ---
include("interactions.jl")
export HollingTypeII, HollingTypeIII
export functional_response
export TrophicLink, TrophicWeb, add_link!, find_links, find_predators
export SITRelease, fertile_mating_fraction

# --- Layer 9: Economics ---
include("economics.jl")
export AbstractCostFunction, AbstractRevenueFunction, AbstractDamageFunction
export FixedCost, VariableCost, InputCostBundle, total_cost
export CropRevenue, revenue
export LinearDamageFunction, ExponentialDamageFunction
export yield_loss, actual_yield
export net_profit, daily_income, npv, benefit_cost_ratio
export RainfallYieldModel, WeatherYieldModel, predict_yield

# --- Layer 10: Genetics ---
include("genetics.jl")
export AbstractGenotypeModel, DialleleicLocus
export genotype_frequencies, GenotypeFitness, selection_step!
export allele_frequency_from_adults
export DoseResponse, mortality_probability
export refuge_dilution
export TwoLocusResistance, probability_fully_resistant

# --- Layer 11: Epidemiology ---
include("epidemiology.jl")
export AbstractDiseaseModel, SIRDisease, DiseaseState
export total_alive, prevalence, step_disease!, R0
export VectorBorneDisease, VectorState, total_vectors, step_vector_disease!

end # module PhysiologicallyBasedDemographicModels
