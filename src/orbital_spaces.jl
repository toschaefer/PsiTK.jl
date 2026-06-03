export merge_spaces
export OrbitalSpace, diagonalize_orbitals
export split_space_occupied_virtual, extract_occupied_space, extract_virtual_space
export select_orbitals

using LinearAlgebra


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
    is_orthonormal::Bool
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
        true # Converged SCF states are always orthonormal
    )
end

"""
    merge_spaces(spaces::OrbitalSpace...)

Merges multiple `OrbitalSpace`s into a single `OrbitalSpace`.
"""
function merge_spaces(spaces::OrbitalSpace{B,T,R}...) where {B,T,R}
    @assert length(spaces) > 0 "Must provide at least one OrbitalSpace to merge."

    basis = spaces[1].basis
    nkpt = length(basis.kpoints)
    is_orthonormal_merged = all(s.is_orthonormal for s in spaces)

    ψ_merged = Matrix{T}[]
    eigenvalues_merged = Vector{R}[]
    occupations_merged = Vector{R}[]

    for ik = 1:nkpt
        # Standard horizontal concatenation of wavefunctions
        ψ_k_blocks = [s.ψ[ik] for s in spaces]
        push!(ψ_merged, reduce(hcat, ψ_k_blocks))

        # Standard concatenation for 1D scalar arrays (negligible memory)
        push!(eigenvalues_merged, vcat([s.eigenvalues[ik] for s in spaces]...))
        push!(occupations_merged, vcat([s.occupations[ik] for s in spaces]...))
    end

    return OrbitalSpace{B,T,R}(
        basis,
        ψ_merged,
        eigenvalues_merged,
        occupations_merged,
        spaces[1].εF,
        is_orthonormal_merged
    )
end

"""
    diagonalize_orbitals(space::OrbitalSpace, hamiltonian)

Diagonalizes the Fock Hamiltonian in the subspace defined by `space.ψ` to
return canonicalized, orthogonal eigenstates.
"""
function diagonalize_orbitals(space::OrbitalSpace{B,T,R}, hamiltonian) where {B,T,R}
    ψ_canon = Matrix{T}[]
    eigenvalues_canon = Vector{R}[]

    for ik = 1:length(space.basis.kpoints)
        X = space.ψ[ik]
        # H is the Hamiltonian applied to the subspace
        HX = hamiltonian[ik] * X
        # Project into the subspace to form a small dense matrix
        h_sub = Hermitian(Matrix(X' * HX))

        if space.is_orthonormal
            res = eigen(h_sub)
        else
            S_sub = Hermitian(Matrix(X' * X))
            res = eigen(h_sub, S_sub)
        end

        # Rotate original wavefunctions to diagonalize 
        push!(ψ_canon, X * res.vectors)
        push!(eigenvalues_canon, res.values)
    end

    return OrbitalSpace{B,T,R}(
        space.basis,
        ψ_canon,
        eigenvalues_canon,
        space.occupations,
        space.εF,
        true 
    )
end



"""
    split_space_occupied_virtual(space::OrbitalSpace; threshold=1e-6)

Splits the given `space` into an occupied and an empty (virtual) `OrbitalSpace` based on the 
fractional occupation `threshold`.
Returns a tuple `(occupied_space, virtual_space)`. Note that this allocates two new `OrbitalSpace` objects and copies the data. If you only need one, use `extract_occupied_space` or `extract_virtual_space`.
"""
function split_space_occupied_virtual(space::OrbitalSpace; threshold=1e-6)
    masks = DFTK.occupied_empty_masks(space.occupations, threshold)
    occ_space = select_orbitals(space, masks.mask_occ)
    empty_space = select_orbitals(space, masks.mask_empty)
    return occ_space, empty_space
end

"""
    extract_occupied_space(space::OrbitalSpace; threshold=1e-6)

Extracts the occupied subspace from the given `space` based on the fractional occupation `threshold`.
"""
function extract_occupied_space(space::OrbitalSpace; threshold=1e-6)
    masks = DFTK.occupied_empty_masks(space.occupations, threshold)
    return select_orbitals(space, masks.mask_occ)
end

"""
    extract_virtual_space(space::OrbitalSpace; threshold=1e-6)

Extracts the empty (virtual) subspace from the given `space` based on the fractional occupation `threshold`.
"""
function extract_virtual_space(space::OrbitalSpace; threshold=1e-6)
    masks = DFTK.occupied_empty_masks(space.occupations, threshold)
    return select_orbitals(space, masks.mask_empty)
end

"""
    select_orbitals(space::OrbitalSpace, indices)

Extracts the orbitals from the given `space` at the specific indices.
`indices` can either be a single `AbstractVector{Int}` (applied to all k-points),
or an `AbstractVector{<:AbstractVector{Int}}` providing a specific list of indices for each k-point.
"""
function select_orbitals(space::OrbitalSpace{B,T,R}, indices::AbstractVector{<:AbstractVector{Int}}) where {B,T,R}
    @assert length(indices) == length(space.ψ)
    ψ_sel = Matrix{T}[]
    eigenvalues_sel = Vector{R}[]
    occupations_sel = Vector{R}[]

    for ik = 1:length(space.ψ)
        push!(ψ_sel, space.ψ[ik][:, indices[ik]])
        push!(eigenvalues_sel, space.eigenvalues[ik][indices[ik]])
        push!(occupations_sel, space.occupations[ik][indices[ik]])
    end

    return OrbitalSpace{B,T,R}(
        space.basis,
        ψ_sel,
        eigenvalues_sel,
        occupations_sel,
        space.εF,
        space.is_orthonormal
    )
end

function select_orbitals(space::OrbitalSpace, indices::AbstractVector{Int})
    # Apply the same indices to all k-points
    return select_orbitals(space, [indices for _ in 1:length(space.ψ)])
end
