"""
    OrbitalSpace{B,T,R<:Real}

A generic representation of an orbital manifold. Decouples the basis and orbitals
from stateful host objects like DFTK's `scfres`.

# Fields
- `basis::B`: The underlying basis (e.g., PlaneWaveBasis).
- `ψ::Vector{Matrix{T}}`: The orbital coefficients per k-point.
- `eigenvalues::Vector{Vector{R}}`: The energies per k-point.
- `occupation::Vector{Vector{R}}`: The fractional occupations per k-point.
- `εF::R`: The Fermi energy of the system.
"""
struct OrbitalSpace{B,T,R<:Real}
    basis::B
    ψ::Vector{Matrix{T}}
    eigenvalues::Vector{Vector{R}}
    occupation::Vector{Vector{R}}
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
        true, # Converged SCF states are always orthonormal
    )
end

"""
    merge_spaces(spaces::OrbitalSpace...)

Merges multiple `OrbitalSpace`s into a single `OrbitalSpace`.
"""
function merge_spaces(
    space::OrbitalSpace{B,T,R},
    more_spaces::OrbitalSpace{B,T,R}...,
) where {B,T,R}
    spaces = (space, more_spaces...)
    basis = space.basis
    is_orthonormal_merged = all(s.is_orthonormal for s in spaces)

    ψ_merged = Matrix{T}[]
    eigenvalues_merged = Vector{R}[]
    occupation_merged = Vector{R}[]

    for ik in eachindex(basis.kpoints)
        # Standard horizontal concatenation of wavefunctions
        ψ_k_blocks = [s.ψ[ik] for s in spaces]
        push!(ψ_merged, reduce(hcat, ψ_k_blocks))

        # Standard concatenation for 1D scalar arrays (negligible memory)
        push!(eigenvalues_merged, vcat([s.eigenvalues[ik] for s in spaces]...))
        push!(occupation_merged, vcat([s.occupation[ik] for s in spaces]...))
    end

    return OrbitalSpace{B,T,R}(
        basis,
        ψ_merged,
        eigenvalues_merged,
        occupation_merged,
        spaces[1].εF,
        is_orthonormal_merged,
    )
end

"""
    canonicalize_orbitals(space::OrbitalSpace, hamiltonian)

Diagonalizes the Fock Hamiltonian in the subspace defined by `space.ψ` to
return canonicalized, orthogonal eigenstates.
"""
function canonicalize_orbitals(space::OrbitalSpace{B,T,R}, hamiltonian) where {B,T,R}
    ψ_canon = Matrix{T}[]
    eigenvalues_canon = Vector{R}[]

    for ik in eachindex(space.basis.kpoints)
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
        space.occupation,
        space.εF,
        true,
    )
end

"""
    split_occupied_virtual(space::OrbitalSpace; threshold=1e-6)

Split `space` into its occupied and virtual orbitals, i.e. those with fractional
occupation above resp. below `threshold`. Returns a tuple `(occupied_space, virtual_space)`.
"""
function split_occupied_virtual(space::OrbitalSpace; threshold=1e-6)
    masks = DFTK.occupied_empty_masks(space.occupation, threshold)
    return select_orbitals(space, masks.mask_occ), select_orbitals(space, masks.mask_empty)
end

"""
    select_orbitals(space::OrbitalSpace, indices)

Extracts the orbitals from the given `space` at the specific indices.
`indices` can either be a single `AbstractVector{Int}` (applied to all k-points),
or an `AbstractVector{<:AbstractVector{Int}}` providing a specific list of indices for
each k-point.
"""
function select_orbitals(
    space::OrbitalSpace{B,T,R},
    indices::AbstractVector{<:AbstractVector{Int}},
) where {B,T,R}
    @assert length(indices) == length(space.ψ)
    ψ_sel = Matrix{T}[]
    eigenvalues_sel = Vector{R}[]
    occupation_sel = Vector{R}[]

    for ik in eachindex(space.ψ)
        push!(ψ_sel, space.ψ[ik][:, indices[ik]])
        push!(eigenvalues_sel, space.eigenvalues[ik][indices[ik]])
        push!(occupation_sel, space.occupation[ik][indices[ik]])
    end

    return OrbitalSpace{B,T,R}(
        space.basis,
        ψ_sel,
        eigenvalues_sel,
        occupation_sel,
        space.εF,
        space.is_orthonormal,
    )
end

function select_orbitals(space::OrbitalSpace, indices::AbstractVector{Int})
    # Apply the same indices to all k-points
    return select_orbitals(space, [indices for _ in eachindex(space.ψ)])
end
