"""
    PBDMProblem

Central problem type for physiologically based demographic models.
Parameterized by trait types following the StructuredPopulationCore convention.

# Type parameters
- `S`: AbstractPBDMStructure (SingleSpeciesPBDM or MultiSpeciesPBDM)
- `D`: AbstractDensityDependence (DensityIndependent or DensityDependent)
- `T`: AbstractStochasticity (Deterministic or StochasticParameterResampled)
- `A`: Approach type (`Nothing` for legacy stepping, or an `AbstractPBDMApproach`)
- `Pop`: Population or Vector{Population} type
- `W`: Weather/forcing type
- `P`: Parameters type
"""
struct PBDMProblem{S<:AbstractPBDMStructure, D<:AbstractDensityDependence,
                   St<:AbstractStochasticity, A, Pop, W, P}
    structure::S
    density::D
    stochasticity::St
    approach::A
    populations::Pop        # Population or Vector{Population}
    weather::W              # WeatherSeries or function
    tspan::Tuple{Int,Int}   # (day_start, day_end) in calendar days
    p::P                    # Additional parameters
end

# Master constructor
function PBDMProblem(structure::AbstractPBDMStructure,
                     density::AbstractDensityDependence,
                     stochasticity::AbstractStochasticity,
                     populations, weather, tspan;
                     approach=nothing,
                     p=nothing)
    t0, tf = tspan
    t0 <= tf || throw(ArgumentError("tspan must satisfy start <= end, got $(tspan)"))
    PBDMProblem(structure, density, stochasticity,
                approach, populations, weather, tspan, p)
end

# Convenience: single population, DI, deterministic
function PBDMProblem(population::Population, weather, tspan::Tuple{Int,Int}; kwargs...)
    PBDMProblem(SingleSpeciesPBDM(), DensityIndependent(), Deterministic(),
                population, weather, tspan; kwargs...)
end

# Convenience: single population, density-dependent
function PBDMProblem(density::DensityDependent,
                     population::Population, weather, tspan::Tuple{Int,Int}; kwargs...)
    PBDMProblem(SingleSpeciesPBDM(), density, Deterministic(),
                population, weather, tspan; kwargs...)
end

# Convenience: vector of populations → multi-species
function PBDMProblem(populations::AbstractVector{<:Population}, weather,
                     tspan::Tuple{Int,Int}; kwargs...)
    PBDMProblem(MultiSpeciesPBDM(), DensityDependent(), Deterministic(),
                populations, weather, tspan; kwargs...)
end

# Convenience with explicit approach
function PBDMProblem(approach::AbstractPBDMApproach,
                     population::Population, weather, tspan::Tuple{Int,Int}; kwargs...)
    PBDMProblem(SingleSpeciesPBDM(), DensityIndependent(), Deterministic(),
                population, weather, tspan; approach=approach, kwargs...)
end

function PBDMProblem(approach::AbstractPBDMApproach, density::DensityDependent,
                     population::Population, weather, tspan::Tuple{Int,Int}; kwargs...)
    PBDMProblem(SingleSpeciesPBDM(), density, Deterministic(),
                population, weather, tspan; approach=approach, kwargs...)
end

function PBDMProblem(approach::AbstractPBDMApproach,
                     populations::AbstractVector{<:Population}, weather,
                     tspan::Tuple{Int,Int}; kwargs...)
    PBDMProblem(MultiSpeciesPBDM(), DensityDependent(), Deterministic(),
                populations, weather, tspan; approach=approach, kwargs...)
end

"""
    remake(prob::PBDMProblem; kwargs...)

Create a new PBDMProblem with selected fields replaced.
"""
function remake(prob::PBDMProblem;
                structure=prob.structure, density=prob.density,
                stochasticity=prob.stochasticity,
                approach=prob.approach,
                populations=prob.populations, weather=prob.weather,
                tspan=prob.tspan, p=prob.p)
    PBDMProblem(structure, density, stochasticity,
                populations, weather, tspan; approach=approach, p=p)
end

function Base.show(io::IO, prob::PBDMProblem)
    sname = typeof(prob.structure).name.name
    dname = typeof(prob.density).name.name
    stname = typeof(prob.stochasticity).name.name
    if prob.approach === nothing
        print(io, "PBDMProblem($sname, $dname, $stname, tspan=$(prob.tspan))")
    else
        print(io, "PBDMProblem($sname, $dname, $stname, approach=$(approach_family(prob.approach)), tspan=$(prob.tspan))")
    end
end
