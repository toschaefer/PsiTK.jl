export DownfoldDSV, CanonicalVirtuals
export generate_orbitals

using LinearAlgebra
using DFTK

function construct_stochastic_orbitals(N, kpt, orbitalType)
    NG = length(kpt.G_vectors)
    radius = rand(NG, N)
    phase = cis.(2π .* rand(NG, N))
    ϕk = zeros(orbitalType, NG, N)
    ϕk .= radius .* phase
    for a in 1:N
        ϕk[:, a] ./= norm(ϕk[:, a])
    end
    qr_decomp = qr(ϕk)
    return Matrix(qr_decomp.Q)
end

"""
    DownfoldDSV(; n_orbitals=100, tol=1e-5, maxiter=500)

Algorithm to compute compressed virtual orbitals (Downfolded Singular Vectors)
by diagonalizing the exact exchange operator in the full virtual space.
"""
Base.@kwdef struct DownfoldDSV <: AbstractOrbitalGenerator
    n_orbitals::Int = 100
    tol::Float64 = 1e-5
    maxiter::Int = 500
end

"""
    CanonicalVirtuals(; n_orbitals=100)

Algorithm to compute canonical virtual orbitals by diagonalizing the Fock
Hamiltonian. Set `n_orbitals = :all` to compute the full virtual plane-wave space.
"""
Base.@kwdef struct CanonicalVirtuals <: AbstractOrbitalGenerator
    n_orbitals::Union{Int, Symbol} = 100
end

"""
    generate_orbitals(alg::DownfoldDSV, scfres, occ_space)

Generates DSVs using DFTK's `LOBPCG` on the exact exchange operator.
"""
function generate_orbitals(alg::DownfoldDSV, scfres, occ_space::OrbitalSpace{B, T, R}) where {B, T, R}
    basis = scfres.basis
    Ecut = basis.Ecut

    # Reconstruct the exchange operator K
    ExactExchangeTerm = only([term for term in basis.terms if term isa DFTK.TermExactExchange])
    _, K = DFTK.ene_ops(ExactExchangeTerm, basis, scfres.ψ, scfres.occupation)
    
    ψ_dsv = Matrix{T}[]
    eigenvalues_dsv = Vector{R}[]
    occupations_dsv = Vector{R}[]

    for ik in 1:length(basis.kpoints)
        kpt = basis.kpoints[ik]
        Kk = K[ik]
        ψocck = occ_space.ψ[ik]
        
        # Stochastic guess
        ϕk = construct_stochastic_orbitals(alg.n_orbitals, kpt, T)
        
        # Build LevelShifted operators
        ε_homo = maximum(occ_space.eigenvalues[ik])
        ham_hf_levelshifted = LevelShiftedOperator(scfres.ham[ik], ψocck, ε_homo, 1e-5, 2 * Ecut)
        
        # Minimum HF eigenvalue over all k-points for safe shift
        shift = abs(minimum(minimum.(scfres.eigenvalues))) + 2.0
        Kk_virt = ProjectedShiftedOperator(Kk, ψocck, shift)
        
        kinetic_preconditioner = DFTK.PreconditionerTPA(scfres.ham[ik].basis, kpt)
        
        # LOBPCG
        dsv = DFTK.LOBPCG(
            Kk_virt, ϕk, ham_hf_levelshifted, kinetic_preconditioner,
            alg.tol, alg.maxiter, callback=DFTK.DefaultLobpcgCallback()
        )
        
        # Recanonicalize DSVs
        qr_decomp = qr(dsv.X)
        X_ortho = Matrix(qr_decomp.Q)
        h_dsv = Hermitian(X_ortho' * (scfres.ham[ik] * X_ortho))
        canonical_dsv_res = eigen(h_dsv)
        
        push!(ψ_dsv, X_ortho * canonical_dsv_res.vectors)
        push!(eigenvalues_dsv, canonical_dsv_res.values)
        push!(occupations_dsv, zeros(R, alg.n_orbitals))
    end
    
    return OrbitalSpace{B, T, R}(basis, ψ_dsv, eigenvalues_dsv, occupations_dsv, occ_space.εF)
end

"""
    generate_orbitals(alg::CanonicalVirtuals, scfres, occ_space)

Generates canonical virtuals from the Fock Hamiltonian.
"""
function generate_orbitals(alg::CanonicalVirtuals, scfres, occ_space::OrbitalSpace{B, T, R}) where {B, T, R}
    basis = scfres.basis
    Ecut = basis.Ecut
    
    ψ_virt = Matrix{T}[]
    eigenvalues_virt = Vector{R}[]
    occupations_virt = Vector{R}[]

    for ik in 1:length(basis.kpoints)
        kpt = basis.kpoints[ik]
        ψocck = occ_space.ψ[ik]
        Nfull = length(kpt.G_vectors)
        
        N_virt = alg.n_orbitals === :all ? Nfull - size(ψocck, 2) : alg.n_orbitals
        
        ϕk_canon = construct_stochastic_orbitals(N_virt, kpt, T)
        ϕk_canon .-= ψocck * (ψocck' * ϕk_canon)   # Project out occupied
        ϕk_canon = Matrix(qr(ϕk_canon).Q)
        
        ε_homo = maximum(occ_space.eigenvalues[ik])
        ham_hf_levelshifted = LevelShiftedOperator(scfres.ham[ik], ψocck, ε_homo, 1e-5, 2 * Ecut)
        kinetic_preconditioner = DFTK.PreconditionerTPA(scfres.ham[ik].basis, kpt)
        
        canon_res = DFTK.LOBPCG(
            ham_hf_levelshifted, ϕk_canon, IdentityOperator(Nfull, T), 
            kinetic_preconditioner, 1e-5, 500, callback=DFTK.DefaultLobpcgCallback()
        )
        
        ε_virtk = canon_res.λ .+ ham_hf_levelshifted.ε_homo .- ham_hf_levelshifted.safe_shift
        
        push!(ψ_virt, canon_res.X)
        push!(eigenvalues_virt, ε_virtk)
        push!(occupations_virt, zeros(R, N_virt))
    end
    
    return OrbitalSpace{B, T, R}(basis, ψ_virt, eigenvalues_virt, occupations_virt, occ_space.εF)
end
