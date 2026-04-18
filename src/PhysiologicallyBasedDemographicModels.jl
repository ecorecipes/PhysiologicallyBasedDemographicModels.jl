module PhysiologicallyBasedDemographicModels

using CommonSolve
using LinearAlgebra
using StructuredPopulationCore
using Statistics

# --- Layer 1: Core types ---
include("types.jl")
export AbstractPBDMStructure, SingleSpeciesPBDM, MultiSpeciesPBDM
export AbstractPBDMApproach, AbstractBiodemographicApproach, AbstractAllocationApproach
export AbstractHybridPBDMApproach, AbstractBiodemographicFunction
export AbstractBiodemographicModel, AbstractAllocationModel, AbstractRespirationModel
export BiodemographicFunctions, CoupledPBDMModel
export approach_family, development_component, acquisition_component, respiration_component
export biodemography, allocation_model
# Re-export shared types from StructuredPopulationCore
export AbstractProjectionStructure
export AbstractDensityDependence, DensityIndependent, DensityDependent
export AbstractStochasticity, Deterministic, StochasticKernelResampled, StochasticParameterResampled
export DirectIteration

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

include("temperature_responses.jl")
export triangular_thermal_scalar, phiT
export briere_rate, fecundity_briere, fecundity_gaussian
export daily_mortality_quadratic
export diapause_fraction_logistic, diapause_fraction_linear
export gilbert_fraser_attack
export make_population

# --- Layer 4: Weather/forcing ---
include("weather.jl")
export AbstractWeather, DailyWeather, WeatherSeries, SinusoidalWeather
export get_weather

# --- Layer 5: Dynamics ---
include("dynamics.jl")
export step_delay!, step_population!, step_system!

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
export cumulative_degree_days, stage_degree_days, stage_trajectory
export net_growth_rate, phenology

# --- Layer 8: Multi-species interactions ---
include("interactions.jl")
export HollingTypeII, HollingTypeIII
export functional_response
export TrophicLink, TrophicWeb, add_link!, find_links, find_predators
export SITRelease, fertile_mating_fraction

# --- Layer 8a: Genetics (before coupled, which depends on DialleleicLocus) ---
include("genetics.jl")
export AbstractGenotypeModel, DialleleicLocus
export genotype_frequencies, GenotypeFitness, selection_step!
export allele_frequency_from_adults
export DoseResponse, mortality_probability
export refuge_dilution
export TwoLocusResistance, probability_fully_resistant

# --- Layer 8b: Coupled population systems ---
include("coupled.jl")
export MultiPopulationPBDM, MultiTypePBDM, MultiSpeciesPBDMNew, MetapopulationPBDM
export AbstractStateVariable, ScalarState, ArrayState, DictState
export get_state, set_state!, has_auto_update, update_state!, snapshot, has_state
export BulkPopulation, step_bulk!
export PopulationComponent, PopulationSystem
export component_total, component_totals, by_species, by_type, by_patch
export inject!, remove_fraction!, set_value!
export AbstractInteractionRule, AbstractScheduledEvent
export PulseRelease, SingleDayRelease, SprayEvent, ConditionalRelease
export WeatherConditionalEvent
export TransferRule, ReproductionRule, MortalityRule, PredationRule
export CompetitionRule, CustomRule
export AbstractStressRule, StressRule, compute_stress
export Observable
export CallbackPhase, PRE_EVENT, POST_EVENT, PRE_STEP, POST_STEP, END_OF_DAY
export PhaseCallback
export CoupledPBDMSolution
export apply_event!, apply_rule!
# GenomeState
export GenomeState, get_genotypes, get_locus, SelectionRule
# DiapauseState + DiapauseRule
export DiapauseState, DiapauseRule
# PhenologyState
export PhenologyState, get_phase, past_milestone, resource_availability
# SoilState
export SoilState, get_inoculum, get_virulence
# SpatialGrid + DispersalRule
export SpatialGrid, DispersalRule, apply_dispersal!
# EnsemblePBDMProblem
export EnsemblePBDMProblem, EnsemblePBDMSolution

# --- Layer 9: Theory / analytical tools ---
include("theory.jl")
export is_ratio_dependent, apparency
export compensation_point
export life_history_strategy
export IsoclineResult, consumer_isocline, resource_isocline
export EquilibriumResult, classify_equilibrium, find_equilibrium
export SpeciesProfile, AssemblyResult, food_web_assembly

# --- Layer 10: Optimal control types ---
include("optimal_control.jl")
export AbstractManagementAction, PesticideControl, BiologicalReleaseControl, HarvestControl
export AbstractManagementObjective, MinimizeDamage, MaximizeProfit
export TrophicLevel, PBDMControlProblem, ManagementSolution
export optimize_management

# --- Layer 10b: Continuous-time formulations ---
include("continuous.jl")
export AbstractContinuousPBDM, AbstractPSPMMethod, AbstractPSPMSpecies
export ContinuousSpecies, ContinuousTrophicLink
export ContinuousPBDMProblem, DelayPBDMProblem, PSPMProblem
export PSPMSpecies, PSPMStage, StagedPSPMSpecies, n_pspm_stages
export FixedMeshUpwind, EscalatorBoxcarTrain, CharacteristicMethod
export ImplicitFixedMeshUpwind, ImplicitEscalatorBoxcarTrain, ImplicitCharacteristicMethod
export LaxFriedrichsUpwind, ImplicitLaxFriedrichsUpwind
export ContinuousPBDMSolution
export species_state_ranges, species_total_ranges
export flatten_population, flatten_populations
export species_trajectory
export solve_continuous, solve_delay, solve_pspm
export staged_species_stage_totals

# --- Layer 11: Economics ---
include("economics.jl")
export AbstractCostFunction, AbstractRevenueFunction, AbstractDamageFunction
export FixedCost, VariableCost, InputCostBundle, total_cost
export CropRevenue, revenue
export LinearDamageFunction, ExponentialDamageFunction
export yield_loss, actual_yield
export net_profit, daily_income, npv, benefit_cost_ratio
export RainfallYieldModel, WeatherYieldModel, predict_yield

# --- Layer 12: Epidemiology ---
include("epidemiology.jl")
export AbstractDiseaseModel, SIRDisease, DiseaseState
export total_alive, prevalence, step_disease!, R0
export VectorBorneDisease, VectorState, total_vectors, step_vector_disease!

# --- Layer 13: Scenario comparison helpers ---
include("scenarios.jl")
export run_scenarios, compare_metrics

# --- Layer 14: Surrogate predictive models ---
include("surrogates.jl")
export LogLinearSurrogate, predict_log, predict, marginal_effects
export variables, enumerate_strategies, pareto_frontier

end # module PhysiologicallyBasedDemographicModels
