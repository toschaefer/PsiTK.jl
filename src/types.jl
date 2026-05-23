export OrbitalSpace


"""
    OrbitalSpace{B, T, R <: Real}

A generic representation of an orbital manifold. Decouples the basis and orbitals
from stateful host objects like DFTK's `scfres`.

# Fields
- `basis::B`: The underlying basis (e.g., PlaneWaveBasis).
- `ψ::Vector{Matrix{T}}`: The orbital coefficients per k-point.
- `eigenvalues::Vector{Vector{R}}`: The energies per k-point.
- `occupations::Vector{Vector{R}}`: The fractional occupations per k-point.
- `εF::R`: The Fermi energy of the system.
"""
struct OrbitalSpace{B,T,R<:Real}
    basis::B
    ψ::Vector{Matrix{T}}
    eigenvalues::Vector{Vector{R}}
    occupations::Vector{Vector{R}}
    εF::R
end

# Extract the space directly from a converged DFTK scfres
function OrbitalSpace(scfres)
    B = typeof(scfres.basis)
    T = eltype(scfres.ψ[1])
    R = eltype(scfres.eigenvalues[1])

    return OrbitalSpace{B,T,R}(
        scfres.basis,
        scfres.ψ,
        scfres.eigenvalues,
        scfres.occupation,
        scfres.εF,
    )
end
