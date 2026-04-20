"""
Multi-species trophic interaction module.

Provides functional responses and coupling mechanisms for multi-species
PBDM systems (plant-herbivore, predator-prey, tritrophic).
"""

# --- Functional responses ---

"""
    HollingTypeII(a, h)

Holling Type II functional response (disc equation):
`f(N) = a × N / (1 + a × h × N)`

- `a`: attack rate (area searched per predator per time)
- `h`: handling time per prey item
"""
struct HollingTypeII{T<:Real} <: AbstractFunctionalResponse
    a::T
    h::T

    function HollingTypeII(a::T, h::T) where {T<:Real}
        a > 0 || throw(ArgumentError("a must be positive"))
        h >= 0 || throw(ArgumentError("h must be non-negative"))
        new{T}(a, h)
    end
end

function HollingTypeII(a::Real, h::Real)
    T = promote_type(typeof(a), typeof(h))
    HollingTypeII(T(a), T(h))
end

"""
    functional_response(fr, prey_density)

Compute per-predator consumption rate from functional response.
"""
function functional_response(fr::HollingTypeII, prey_density::Real)
    prey_density <= 0 && return 0.0
    return fr.a * prey_density / (1 + fr.a * fr.h * prey_density)
end

"""
    HollingTypeIII(a_max, b, h)

Holling Type III (sigmoid) functional response:
`f(N) = a_max × N² / (b² + N²) / (1 + a_max × h × N² / (b² + N²))`

Simplified form: `f(N) = a_max × N² / (b² + N²) × 1/(1 + h × ...)`
"""
struct HollingTypeIII{T<:Real} <: AbstractFunctionalResponse
    a_max::T
    b::T
    h::T
end

function functional_response(fr::HollingTypeIII, prey_density::Real)
    prey_density <= 0 && return 0.0
    attack = fr.a_max * prey_density^2 / (fr.b^2 + prey_density^2)
    return attack / (1 + attack * fr.h)
end

# Extend FraserGilbertResponse with functional_response dispatch
function functional_response(fr::FraserGilbertResponse, supply::Real, demand::Real)
    return acquire(fr, supply, demand)
end

# --- Trophic interaction ---

"""
    TrophicLink(predator_name, prey_name, response, conversion_efficiency)

A directed trophic link between two populations.
`conversion_efficiency` is the fraction of consumed prey biomass converted to predator biomass.
"""
struct TrophicLink{T<:Real, FR<:AbstractFunctionalResponse}
    predator_name::Symbol
    prey_name::Symbol
    response::FR
    conversion_efficiency::T
end

"""
    TrophicWeb(links)

A collection of trophic links defining a food web.
"""
struct TrophicWeb{T<:Real}
    links::Vector{TrophicLink{T}}
end

TrophicWeb() = TrophicWeb{Float64}(TrophicLink{Float64}[])

"""
    add_link!(web::TrophicWeb, link::TrophicLink) -> TrophicWeb

Append a [`TrophicLink`](@ref) to `web` in place and return the web.
"""
function add_link!(web::TrophicWeb, link::TrophicLink)
    push!(web.links, link)
    return web
end

"""
    find_links(web, predator_name)

Find all trophic links where the given species is the predator.
"""
function find_links(web::TrophicWeb, predator_name::Symbol)
    return filter(l -> l.predator_name == predator_name, web.links)
end

"""
    find_predators(web, prey_name)

Find all trophic links where the given species is the prey.
"""
function find_predators(web::TrophicWeb, prey_name::Symbol)
    return filter(l -> l.prey_name == prey_name, web.links)
end

# --- SIT mating model ---

"""
    SITRelease(sterile_males, competitiveness, release_interval)

Sterile insect technique release parameters.
- `sterile_males`: number released per event
- `competitiveness`: mating competitiveness relative to wild males (0-1)
- `release_interval`: days between releases
"""
struct SITRelease{T<:Real}
    sterile_males::T
    competitiveness::T
    release_interval::Int
end

"""
    fertile_mating_fraction(wild_males, sit::SITRelease, day)

Compute the fraction of matings that are fertile given SIT releases.
"""
function fertile_mating_fraction(wild_males::Real, sit::SITRelease, day::Int)
    # Sterile males present if within release cycle
    sterile = (day % sit.release_interval == 0) ? sit.sterile_males : sit.sterile_males * 0.9^(day % sit.release_interval)
    effective_sterile = sterile * sit.competitiveness
    total = wild_males + effective_sterile
    total ≈ 0.0 && return 1.0
    return wild_males / total
end
