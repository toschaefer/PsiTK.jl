export OccupiedOrbitals
export select_orbitals

"""
    OccupiedOrbitals(; threshold=1e-6)

Algorithm to select the occupied subspace from an `OrbitalSpace` based on a given
fractional occupation threshold.
"""
Base.@kwdef struct OccupiedOrbitals <: AbstractOrbitalSelector
    threshold::Float64 = 1e-6
end

"""
    select_orbitals(space::OrbitalSpace, alg::OccupiedOrbitals)

Extracts the occupied orbitals from the given `space` where occupation > `alg.threshold`.
"""
function select_orbitals(space::OrbitalSpace{B, T, R}, alg::OccupiedOrbitals) where {B, T, R}
    ψ_occ = Matrix{T}[]
    eigenvalues_occ = Vector{R}[]
    occupations_occ = Vector{R}[]

    for ik in 1:length(space.ψ)
        occ = space.occupations[ik]
        
        # Determine how many bands are occupied
        n_occ = count(x -> x > alg.threshold, occ)
        
        push!(ψ_occ, space.ψ[ik][:, 1:n_occ])
        push!(eigenvalues_occ, space.eigenvalues[ik][1:n_occ])
        push!(occupations_occ, space.occupations[ik][1:n_occ])
    end

    return OrbitalSpace{B, T, R}(
        space.basis,
        ψ_occ,
        eigenvalues_occ,
        occupations_occ,
        space.εF
    )
end
