"""
    OrbitalSpace{B,T,R<:Real}
    OrbitalSpace(scfres)

A generic representation of an orbital manifold. Decouples the basis and orbitals
from stateful host objects like DFTK's `scfres`. `OrbitalSpace(scfres)` takes all orbitals
of a converged DFTK SCF, occupied and unoccupied (see [`split_occupied_virtual`](@ref)).

# Fields
- `basis::B`: The underlying basis (e.g., PlaneWaveBasis).
- `ψ::Vector{Matrix{T}}`: The orbital coefficients per k-point.
- `eigenvalues::Vector{Vector{R}}`: The energies per k-point.
- `occupation::Vector{Vector{R}}`: The fractional occupations per k-point.
- `εF::R`: The Fermi energy of the system.
- `is_orthonormal::Bool`: Whether the orbitals are orthonormal, which e.g.
  [`DensitySpecificVirtuals`](@ref) are not (see [`canonicalize_orbitals`](@ref)).
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
    canonicalize_orbitals(space::OrbitalSpace, hamiltonian; occupation_tol=1e-4)

Diagonalizes the Fock Hamiltonian in the subspace defined by `space.ψ` to return
canonicalized, orthonormal eigenstates sorted by energy. The occupations of `space` are
redistributed by the aufbau principle: the largest occupations are assigned to the lowest
eigenvalues. The Fermi level `εF` is placed in the middle of the gap between occupied and
empty orbitals (kept as is if there is no gap).

The aufbau occupation reproduces the input one only if the diagonalization does not mix
orbitals of different occupation, e.g. for a converged SCF canonicalized with its own
Hamiltonian. Otherwise, e.g. for DFT orbitals canonicalized with a Fock operator, it
defines a different determinant, which is not self-consistent with `hamiltonian`. A
warning is issued if the input density matrix, expressed in the canonical orbitals,
deviates from the aufbau occupation by more than `occupation_tol` in any element.

A deviation `d` corresponds to an occupied-virtual mixing angle `θ ≈ d / f` (with `f` the
occupation, e.g. 2) and changes the energy of the determinant by the order of `θ²` times
the orbital energy gap. The default thus flags changes beyond those left by an SCF
converged to about 1e-9 Hartree in the total energy.
"""
function canonicalize_orbitals(
    space::OrbitalSpace{B,T,R},
    hamiltonian;
    occupation_tol=1e-4,
) where {B,T,R}
    ψ_canon = Matrix{T}[]
    eigenvalues_canon = Vector{R}[]
    occupation_canon = Vector{R}[]

    for ik in eachindex(space.basis.kpoints)
        X = space.ψ[ik]
        f = space.occupation[ik]
        # H is the Hamiltonian applied to the subspace
        HX = hamiltonian[ik] * X
        # Project into the subspace to form a small dense matrix
        h_sub = Hermitian(Matrix(X' * HX))

        # Eigenvectors C, ascending in energy, with S C: the overlap of X with the
        # canonical orbitals X C (S = X'X = I for orthonormal X)
        if space.is_orthonormal
            res = eigen(h_sub)
            SC = res.vectors
        else
            S_sub = Hermitian(Matrix(X' * X))
            res = eigen(h_sub, S_sub)
            SC = S_sub * res.vectors
        end

        # Aufbau: largest occupations to the lowest eigenvalues
        f_aufbau = sort(f; rev=true)

        # Input density matrix Σ_i f_i |x_i⟩⟨x_i| in the canonical orbitals. It equals
        # diag(f_aufbau) iff the rotation only mixes orbitals of equal occupation, otherwise
        # the deviation grows linearly with the occupied-virtual mixing angle
        γ = SC' * Diagonal(f) * SC
        deviation = maximum(abs, γ - Diagonal(f_aufbau))
        if deviation > occupation_tol
            @warn "canonicalize_orbitals mixed orbitals of different occupation " *
                  "(max. deviation $deviation at k-point $ik): the aufbau occupation " *
                  "defines a new determinant, not self-consistent with the Hamiltonian."
        end

        push!(ψ_canon, X * res.vectors)
        push!(eigenvalues_canon, res.values)
        push!(occupation_canon, f_aufbau)
    end

    return OrbitalSpace{B,T,R}(
        space.basis,
        ψ_canon,
        eigenvalues_canon,
        occupation_canon,
        _fermi_level_in_gap(eigenvalues_canon, occupation_canon, space.εF),
        true,
    )
end

"""
    split_occupied_virtual(space::OrbitalSpace; threshold=1e-6)

Split `space` into its occupied and virtual orbitals, i.e. those with fractional
occupation above resp. below `threshold`. Returns a tuple `(occupied_space, virtual_space)`.
"""
function split_occupied_virtual(space::OrbitalSpace; threshold=1e-6)
    (; occupied, empty) = _occupied_empty_indices(space.occupation, threshold)
    return select_orbitals(space, occupied), select_orbitals(space, empty)
end

# Indices of the occupied and empty orbitals per k-point. Unlike DFTK.occupied_empty_masks,
# this does not assume the occupied orbitals to come first (e.g. in merged spaces)
function _occupied_empty_indices(occupation, threshold)
    occupied = [findall(f -> abs(f) > threshold, occ) for occ in occupation]
    empty = [findall(f -> abs(f) <= threshold, occ) for occ in occupation]
    return (; occupied, empty)
end

# Fermi level in the middle of the gap between occupied and empty orbitals. The given εF
# is kept if there is no such gap: only occupied or only empty orbitals, or fractional
# occupations around the Fermi level
function _fermi_level_in_gap(eigenvalues, occupation, εF; threshold=1e-6)
    (; occupied, empty) = _occupied_empty_indices(occupation, threshold)
    ε_occ = mapreduce((ε, i) -> ε[i], vcat, eigenvalues, occupied)
    ε_empty = mapreduce((ε, i) -> ε[i], vcat, eigenvalues, empty)
    (isempty(ε_occ) || isempty(ε_empty)) && return εF
    ε_homo, ε_lumo = maximum(ε_occ), minimum(ε_empty)
    return ε_homo < ε_lumo ? (ε_homo + ε_lumo) / 2 : εF
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
