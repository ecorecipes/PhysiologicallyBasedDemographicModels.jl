"""
Genetics module for tracking resistance allele frequencies.

Implements Hardy-Weinberg single-locus genetics for modeling pesticide
and Bt toxin resistance evolution in pest populations.
"""

# --- Types ---

"""
    AbstractGenotypeModel

Supertype for genetic tracking models.
"""
abstract type AbstractGenotypeModel end

"""
    DialleleicLocus(R, dominance)

A single diallelic locus with allele frequency R (resistant) and 1-R (susceptible).
`dominance` ∈ [0,1]: 0 = fully recessive resistance, 0.5 = codominant, 1 = fully dominant.
"""
mutable struct DialleleicLocus{T<:Real} <: AbstractGenotypeModel
    R::T
    dominance::T

    function DialleleicLocus(R::T, dominance::T) where {T<:Real}
        0 <= R <= 1 || throw(ArgumentError("R must be in [0,1]"))
        0 <= dominance <= 1 || throw(ArgumentError("dominance must be in [0,1]"))
        new{T}(R, dominance)
    end
end

function DialleleicLocus(R::Real, dominance::Real)
    T = promote_type(typeof(float(R)), typeof(float(dominance)))
    DialleleicLocus(T(R), T(dominance))
end

DialleleicLocus(R::Real) = DialleleicLocus(R, zero(R))

"""
    genotype_frequencies(locus::DialleleicLocus)

Return Hardy-Weinberg genotype frequencies (SS, SR, RR).
"""
function genotype_frequencies(locus::DialleleicLocus)
    R = locus.R
    S = 1 - R
    return (SS=S^2, SR=2*S*R, RR=R^2)
end

"""
    genotype_frequencies(R::Real)

Return Hardy-Weinberg genotype frequencies from allele frequency R.
"""
function genotype_frequencies(R::Real)
    S = 1 - R
    return (SS=S^2, SR=2*S*R, RR=R^2)
end

# --- Fitness and selection ---

"""
    GenotypeFitness(w_SS, w_SR, w_RR)

Relative fitness values for each genotype.
"""
struct GenotypeFitness{T<:Real}
    w_SS::T
    w_SR::T
    w_RR::T
end

"""
    selection_step!(locus, fitness)

Update allele frequency after one generation of selection.
Uses the standard population genetics recursion.
"""
function selection_step!(locus::DialleleicLocus, fitness::GenotypeFitness)
    freq = genotype_frequencies(locus)
    w_bar = freq.SS * fitness.w_SS + freq.SR * fitness.w_SR + freq.RR * fitness.w_RR
    w_bar <= 0 && return locus.R
    R_new = (freq.SR * fitness.w_SR * 0.5 + freq.RR * fitness.w_RR) / w_bar
    locus.R = clamp(R_new, 0.0, 1.0)
    return locus.R
end

"""
    allele_frequency_from_adults(n_SS, n_SR, n_RR)

Compute resistance allele frequency from genotype counts.
"""
function allele_frequency_from_adults(n_SS::Real, n_SR::Real, n_RR::Real)
    n_total = n_SS + n_SR + n_RR
    n_total ≈ 0.0 && return 0.0
    return (2 * n_RR + n_SR) / (2 * n_total)
end

# --- Dose-response ---

"""
    DoseResponse(ld50, slope)

Probit/logistic dose-response curve. `ld50` is the dose killing 50%
of susceptible individuals; `slope` controls steepness.
"""
struct DoseResponse{T<:Real}
    ld50::T
    slope::T
end

"""
    mortality_probability(dr::DoseResponse, dose)

Compute mortality probability at a given dose.
Uses a logistic model: `P(death) = 1 / (1 + (ld50/dose)^slope)`.
"""
function mortality_probability(dr::DoseResponse, dose::Real)
    dose <= 0 && return 0.0
    return 1.0 / (1.0 + (dr.ld50 / dose)^dr.slope)
end

# --- Refuge dynamics ---

"""
    refuge_dilution(R_field, R_refuge, refuge_fraction)

Dilute resistance allele frequency by gene flow from a refuge.
"""
function refuge_dilution(R_field::Real, R_refuge::Real, refuge_fraction::Real)
    return (1 - refuge_fraction) * R_field + refuge_fraction * R_refuge
end

# --- Multi-locus (two-toxin pyramid) ---

"""
    TwoLocusResistance(locus1, locus2)

Independent two-locus resistance model (e.g., Cry1Ac + Cry2Ab).
"""
struct TwoLocusResistance{T<:Real} <: AbstractGenotypeModel
    locus1::DialleleicLocus{T}
    locus2::DialleleicLocus{T}
end

"""
    probability_fully_resistant(tlr::TwoLocusResistance)

Probability of being homozygous resistant at both loci (independent).
"""
function probability_fully_resistant(tlr::TwoLocusResistance)
    f1 = genotype_frequencies(tlr.locus1)
    f2 = genotype_frequencies(tlr.locus2)
    return f1.RR * f2.RR
end
