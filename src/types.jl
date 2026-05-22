export OrbitalSpace
export AbstractOrbitalAlgorithm
export AbstractOrbitalSelector, AbstractOrbitalGenerator, AbstractOrbitalTransformer

"""
    OrbitalSpace{B, T, R <: Real}

A generic representation of an orbital manifold. Decouples the basis and orbitals
from stateful host objects like DFTK's `scfres`.

# Fields
- `basis::B`: The underlying basis (e.g., PlaneWaveBasis).
- `ψ::Vector{AbstractMatrix{T}}`: The orbital coefficients per k-point.
- `eigenvalues::Vector{Vector{R}}`: The energies per k-point.
- `occupations::Vector{Vector{R}}`: The fractional occupations per k-point.
- `εF::R`: The Fermi energy of the system.
"""
struct OrbitalSpace{B, T, R <: Real}
    basis::B
    ψ::Vector{AbstractMatrix{T}}
    eigenvalues::Vector{Vector{R}}
    occupations::Vector{Vector{R}}
    εF::R
end

# Extract the space directly from a converged DFTK scfres
function OrbitalSpace(scfres)
    # The scfres.ψ is usually a Vector{Matrix{T}}. We convert the signature to AbstractMatrix
    # to allow mixing dense matrices with BlockArrays later.
    B = typeof(scfres.basis)
    T = eltype(scfres.ψ[1])
    R = eltype(scfres.eigenvalues[1])
    
    ψ_abs = AbstractMatrix{T}[scfres.ψ[ik] for ik in 1:length(scfres.ψ)]
    
    return OrbitalSpace{B, T, R}(
        scfres.basis,
        ψ_abs,
        scfres.eigenvalues,
        scfres.occupation,
        scfres.εF
    )
end

# Algorithm Traits
abstract type AbstractOrbitalAlgorithm end
abstract type AbstractOrbitalSelector <: AbstractOrbitalAlgorithm end
abstract type AbstractOrbitalGenerator <: AbstractOrbitalAlgorithm end
abstract type AbstractOrbitalTransformer <: AbstractOrbitalAlgorithm end
