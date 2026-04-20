"""
Epidemiology module for vector-borne and livestock disease dynamics.

Provides SIR-type compartmental disease models that couple with
PBDM population dynamics for host-vector-pathogen systems.
"""

# --- Types ---

"""
    AbstractDiseaseModel

Supertype for compartmental disease models in the epidemiology module
(e.g. [`SIRDisease`](@ref), [`VectorBorneDisease`](@ref)). Concrete
subtypes are advanced via [`step_disease!`](@ref) /
[`step_vector_disease!`](@ref).
"""
abstract type AbstractDiseaseModel end

"""
    SIRDisease(β, γ, μ_d)

Simple SIR disease model with:
- `β`: transmission rate (per contact per day)
- `γ`: recovery rate (per day)
- `μ_d`: disease-induced mortality rate (per day)
"""
struct SIRDisease{T<:Real} <: AbstractDiseaseModel
    β::T
    γ::T
    μ_d::T

    function SIRDisease(β::T, γ::T, μ_d::T) where {T<:Real}
        β >= 0 || throw(ArgumentError("β must be non-negative"))
        γ >= 0 || throw(ArgumentError("γ must be non-negative"))
        μ_d >= 0 || throw(ArgumentError("μ_d must be non-negative"))
        new{T}(β, γ, μ_d)
    end
end

function SIRDisease(β::Real, γ::Real, μ_d::Real)
    T = promote_type(typeof(β), typeof(γ), typeof(μ_d))
    SIRDisease(T(β), T(γ), T(μ_d))
end

"""
    DiseaseState(S, I, R, D)

Compartmental disease state: Susceptible, Infected, Recovered, Dead-from-disease.
"""
mutable struct DiseaseState{T<:Real}
    S::T
    I::T
    R::T
    D::T
end

DiseaseState(S::Real, I::Real) = DiseaseState(float(S), float(I), 0.0, 0.0)
DiseaseState(N::Real) = DiseaseState(float(N), 0.0, 0.0, 0.0)

"""
    total_alive(ds::DiseaseState) -> Real

Sum of `S + I + R` (excludes disease-induced deaths `D`).
"""
total_alive(ds::DiseaseState) = ds.S + ds.I + ds.R

"""
    prevalence(ds::DiseaseState) -> Real

Fraction of the alive sub-population currently infected: `I / (S+I+R)`.
Returns `0.0` if the alive sub-population is empty.
"""
prevalence(ds::DiseaseState) = total_alive(ds) > 0 ? ds.I / total_alive(ds) : 0.0

"""
    step_disease!(state, disease, contact_rate; dt=1.0)

Advance disease dynamics by one time step.
Returns named tuple of new infections, recoveries, and deaths.
"""
function step_disease!(state::DiseaseState, disease::SIRDisease,
                       contact_rate::Real; dt::Real=1.0)
    N = total_alive(state)
    N <= 0 && return (new_infections=0.0, recoveries=0.0, deaths=0.0)

    # Force of infection
    foi = disease.β * contact_rate * state.I / N
    new_infections = foi * state.S * dt
    new_infections = min(new_infections, state.S)

    # Recovery and disease death
    recoveries = disease.γ * state.I * dt
    deaths = disease.μ_d * state.I * dt
    total_loss = recoveries + deaths
    if total_loss > state.I
        scale = state.I / total_loss
        recoveries *= scale
        deaths *= scale
    end

    # Update compartments
    state.S -= new_infections
    state.I += new_infections - recoveries - deaths
    state.R += recoveries
    state.D += deaths

    # Clamp
    state.S = max(0.0, state.S)
    state.I = max(0.0, state.I)
    state.R = max(0.0, state.R)

    return (new_infections=new_infections, recoveries=recoveries, deaths=deaths)
end

# --- Vector-borne transmission ---

"""
    VectorBorneDisease(β_vh, β_hv, γ_h, μ_h, extrinsic_incubation)

Vector-borne disease transmission model.
- `β_vh`: vector-to-host transmission probability per bite
- `β_hv`: host-to-vector transmission probability per bite
- `γ_h`: host recovery rate
- `μ_h`: host disease-induced mortality
- `extrinsic_incubation`: days for pathogen to develop in vector
"""
struct VectorBorneDisease{T<:Real} <: AbstractDiseaseModel
    β_vh::T
    β_hv::T
    γ_h::T
    μ_h::T
    extrinsic_incubation::T
end

"""
    VectorState(S, E, I)

Vector (insect) disease state: Susceptible, Exposed (incubating), Infectious.
"""
mutable struct VectorState{T<:Real}
    S::T
    E::T
    I::T
end

VectorState(N::Real) = VectorState(float(N), 0.0, 0.0)

"""
    total_vectors(vs::VectorState) -> Real

Total vector population `S + E + I` for a [`VectorState`](@ref).
"""
total_vectors(vs::VectorState) = vs.S + vs.E + vs.I

"""
    step_vector_disease!(host_state, vector_state, disease, bite_rate; dt=1.0)

Advance vector-borne disease dynamics for one time step.
"""
function step_vector_disease!(host::DiseaseState, vector::VectorState,
                              disease::VectorBorneDisease, bite_rate::Real;
                              dt::Real=1.0)
    N_h = total_alive(host)
    N_v = total_vectors(vector)
    (N_h <= 0 || N_v <= 0) && return (host_infections=0.0, vector_infections=0.0)

    # Vector → Host transmission
    infective_bites = bite_rate * vector.I * dt
    host_infections = disease.β_vh * infective_bites * (host.S / N_h)
    host_infections = min(host_infections, host.S)

    # Host → Vector transmission
    vector_bites_on_infected = bite_rate * (host.I / N_h) * vector.S * dt
    vector_infections = disease.β_hv * vector_bites_on_infected
    vector_infections = min(vector_infections, vector.S)

    # Host disease dynamics
    recoveries = disease.γ_h * host.I * dt
    deaths = disease.μ_h * host.I * dt

    host.S -= host_infections
    host.I += host_infections - recoveries - deaths
    host.R += recoveries
    host.D += deaths
    host.S = max(0.0, host.S)
    host.I = max(0.0, host.I)

    # Vector dynamics (E → I transition based on incubation)
    ei_rate = dt / max(1.0, disease.extrinsic_incubation)
    becoming_infectious = vector.E * ei_rate

    vector.S -= vector_infections
    vector.E += vector_infections - becoming_infectious
    vector.I += becoming_infectious
    vector.S = max(0.0, vector.S)
    vector.E = max(0.0, vector.E)

    return (host_infections=host_infections, vector_infections=vector_infections)
end

# --- Basic reproduction number ---

"""
    R0(disease::SIRDisease, contact_rate)

Basic reproduction number for SIR model.
"""
function R0(disease::SIRDisease, contact_rate::Real)
    return disease.β * contact_rate / (disease.γ + disease.μ_d)
end
