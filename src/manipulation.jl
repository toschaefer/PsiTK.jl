export merge_spaces
export canonicalize_orbitals

using LinearAlgebra

"""
    merge_spaces(spaces::OrbitalSpace...)

Merges multiple `OrbitalSpace`s into a single `OrbitalSpace`.
"""
function merge_spaces(spaces::OrbitalSpace{B, T, R}...) where {B, T, R}
    @assert length(spaces) > 0 "Must provide at least one OrbitalSpace to merge."
    
    basis = spaces[1].basis
    nkpt = length(basis.kpoints)
    
    ψ_merged = Matrix{T}[]
    eigenvalues_merged = Vector{R}[]
    occupations_merged = Vector{R}[]

    for ik in 1:nkpt
        # Standard horizontal concatenation of wavefunctions
        ψ_k_blocks = [s.ψ[ik] for s in spaces]
        push!(ψ_merged, reduce(hcat, ψ_k_blocks))
        
        # Standard concatenation for 1D scalar arrays (negligible memory)
        push!(eigenvalues_merged, vcat([s.eigenvalues[ik] for s in spaces]...))
        push!(occupations_merged, vcat([s.occupations[ik] for s in spaces]...))
    end
    
    return OrbitalSpace{B, T, R}(
        basis,
        ψ_merged,
        eigenvalues_merged,
        occupations_merged,
        spaces[1].εF
    )
end

"""
    canonicalize_orbitals(space::OrbitalSpace, hamiltonian)

Diagonalizes the Fock Hamiltonian in the subspace defined by `space.ψ` to
return canonicalized, orthogonal eigenstates.
"""
function canonicalize_orbitals(space::OrbitalSpace{B, T, R}, hamiltonian) where {B, T, R}
    ψ_canon = Matrix{T}[]
    eigenvalues_canon = Vector{R}[]
    
    for ik in 1:length(space.basis.kpoints)
        X = space.ψ[ik]
        # H is the Hamiltonian applied to the subspace
        HX = hamiltonian[ik] * X
        # Project into the subspace to form a small dense matrix
        h_sub = Hermitian(Matrix(X' * HX))
        
        res = eigen(h_sub)
        
        # Rotate original wavefunctions to diagonalize
        push!(ψ_canon, X * res.vectors)
        push!(eigenvalues_canon, res.values)
    end
    
    return OrbitalSpace{B, T, R}(
        space.basis,
        ψ_canon,
        eigenvalues_canon,
        space.occupations,
        space.εF
    )
end
