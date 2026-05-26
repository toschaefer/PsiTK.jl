export DensitySpecificVirtuals, CanonicalVirtuals
export generate_orbitals

using LinearAlgebra
using DFTK

function construct_stochastic_orbitals(N, kpt, orbitalType)
    NG = length(kpt.G_vectors)
    radius = rand(NG, N)
    phase = cis.(2π .* rand(NG, N))
    ϕk = zeros(orbitalType, NG, N)
    ϕk .= radius .* phase
    for a = 1:N
        ϕk[:, a] ./= norm(ϕk[:, a])
    end
    qr_decomp = qr(ϕk)
    return Matrix(qr_decomp.Q)
end

# --- Targets ---

@doc raw"""
    DensitySpecificVirtuals(; n_orbitals)

Target to compute compressed virtual orbitals (Density Specific Virtuals)
by solving a generalized eigenvalue problem in the full virtual space:
```math
\mathcal K \varphi  =  \lambda h \varphi 
```
where $\mathcal K$ and $h$ are the Fock exchange operator and the Fock Hamiltonian, respectively.

When generated, the orbitals returned ($\psi_i$) will be strictly orthogonalized spanning the DSV subspace.
The associated orbital energies will strictly track the generalized Rayleigh quotient ($\lambda_i = \langle \psi_i | \mathcal K | \psi_i \rangle / \langle \psi_i | \mathcal h | \psi_i \rangle$) evaluated exactly within the orthogonalized states.
"""
Base.@kwdef struct DensitySpecificVirtuals
    n_orbitals::Int
end

"""
    CanonicalVirtuals(; n_orbitals)

Target to compute canonical virtual orbitals by diagonalizing the Fock
Hamiltonian. Set `n_orbitals = :all` to compute the full virtual plane-wave space.
"""
Base.@kwdef struct CanonicalVirtuals
    n_orbitals::Union{Int,Symbol}
end

# --- Generators ---

"""
    generate_orbitals(target::DensitySpecificVirtuals, scfres, occ_space, solver::LOBPCGEigensolver)

Generates DSVs using an iterative LOBPCG solver.
"""
function generate_orbitals(
    target::DensitySpecificVirtuals,
    scfres,
    occ_space::OrbitalSpace{B,T,R},
    solver::LOBPCGEigensolver,
) where {B,T,R}
    basis = scfres.basis
    Ecut = basis.Ecut

    # Reconstruct the exchange operator K
    ExactExchangeTerm =
        only([term for term in basis.terms if term isa DFTK.TermExactExchange])
    _, K = DFTK.ene_ops(ExactExchangeTerm, basis, scfres.ψ, scfres.occupation)

    ψ_dsv = Matrix{T}[]
    eigenvalues_dsv = Vector{R}[]
    occupations_dsv = Vector{R}[]

    for ik = 1:length(basis.kpoints)
        kpt = basis.kpoints[ik]
        Kk = K[ik]
        ψocck = occ_space.ψ[ik]
        Nfull = length(kpt.G_vectors)

        if target.n_orbitals > 0.1 * Nfull
            @warn "DSV n_orbitals ($(target.n_orbitals)) is > 10% of plane waves ($Nfull). Full diagonalization might be faster."
        end

        # Stochastic guess
        ϕk = construct_stochastic_orbitals(target.n_orbitals, kpt, T)

        # Build LevelShifted operators
        ε_homo = maximum(occ_space.eigenvalues[ik])
        ham_hf_levelshifted =
            LevelShiftedOperator(scfres.ham[ik], ψocck, ε_homo, 1e-5, 2 * Ecut)

        # Minimum HF eigenvalue over all k-points for safe shift
        shift = abs(minimum(minimum.(scfres.eigenvalues))) + 2.0
        Kk_virt = ProjectedShiftedOperator(Kk, ψocck, shift)

        kinetic_preconditioner = DFTK.PreconditionerTPA(scfres.ham[ik].basis, kpt)

        # LOBPCG
        dsv = DFTK.LOBPCG(
            Kk_virt,
            ϕk,
            ham_hf_levelshifted,
            kinetic_preconditioner,
            solver.tol,
            solver.maxiter,
            callback = DFTK.DefaultLobpcgCallback(),
        )

        # Orthogonalize DSVs (re-canonicalization is left to the user if needed)
        qr_decomp = qr(dsv.X)
        X_ortho = Matrix(qr_decomp.Q)
        
        # Estimate eigenvalues via the generalized Rayleigh quotient for the orthogonalized orbitals
        # y_i = <φ_i | A | φ_i> / <φ_i | B | φ_i> where A = Kk_virt and B = ham_hf_levelshifted
        # Since A X = B X Λ and X' B X = I, we have X' A X = Λ
        # With X = Q R  =>  Q = X R^{-1}, letting invR = R \ I, we get Q' A Q = invR' Λ invR and Q' B Q = invR' invR
        invR = qr_decomp.R \ I
        num = real(diag(invR' * Diagonal(dsv.λ) * invR))
        den = real(diag(invR' * invR))
        h_dsv_diag = num ./ den

        push!(ψ_dsv, X_ortho)
        push!(eigenvalues_dsv, h_dsv_diag)
        push!(occupations_dsv, zeros(R, target.n_orbitals))
    end

    return OrbitalSpace{B,T,R}(basis, ψ_dsv, eigenvalues_dsv, occupations_dsv, occ_space.εF)
end

"""
    generate_orbitals(target::CanonicalVirtuals, scfres, occ_space, solver::LOBPCGEigensolver)

Generates canonical virtuals using an iterative LOBPCG solver on the Fock operator.
"""
function generate_orbitals(
    target::CanonicalVirtuals,
    scfres,
    occ_space::OrbitalSpace{B,T,R},
    solver::LOBPCGEigensolver,
) where {B,T,R}
    basis = scfres.basis
    Ecut = basis.Ecut

    ψ_virt = Matrix{T}[]
    eigenvalues_virt = Vector{R}[]
    occupations_virt = Vector{R}[]

    for ik = 1:length(basis.kpoints)
        kpt = basis.kpoints[ik]
        ψocck = occ_space.ψ[ik]
        Nfull = length(kpt.G_vectors)

        N_virt = target.n_orbitals === :all ? (Nfull - size(ψocck, 2)) : target.n_orbitals

        if N_virt > 0.1 * Nfull
            @warn "CanonicalVirtuals n_orbitals ($N_virt) is > 10% of plane waves ($Nfull). FullDiagonalizationEigensolver might be faster."
        end

        ϕk_canon = construct_stochastic_orbitals(N_virt, kpt, T)
        ϕk_canon .-= ψocck * (ψocck' * ϕk_canon)   # Project out occupied
        ϕk_canon = Matrix(qr(ϕk_canon).Q)

        ε_homo = maximum(occ_space.eigenvalues[ik])
        ham_hf_levelshifted =
            LevelShiftedOperator(scfres.ham[ik], ψocck, ε_homo, 1e-5, 2 * Ecut)
        kinetic_preconditioner = DFTK.PreconditionerTPA(scfres.ham[ik].basis, kpt)

        canon_res = DFTK.LOBPCG(
            ham_hf_levelshifted,
            ϕk_canon,
            IdentityOperator(Nfull, T),
            kinetic_preconditioner,
            solver.tol,
            solver.maxiter,
            callback = DFTK.DefaultLobpcgCallback(),
        )

        ε_virtk =
            canon_res.λ .+ ham_hf_levelshifted.ε_homo .- ham_hf_levelshifted.safe_shift

        push!(ψ_virt, canon_res.X)
        push!(eigenvalues_virt, ε_virtk)
        push!(occupations_virt, zeros(R, N_virt))
    end

    return OrbitalSpace{B,T,R}(
        basis,
        ψ_virt,
        eigenvalues_virt,
        occupations_virt,
        occ_space.εF,
    )
end

"""
    generate_orbitals(target::CanonicalVirtuals, scfres, occ_space, solver::FullDiagonalizationEigensolver)

Generates canonical virtuals using full dense exact diagonalization on the Fock operator.
"""
function generate_orbitals(
    target::CanonicalVirtuals,
    scfres,
    occ_space::OrbitalSpace{B,T,R},
    solver::FullDiagonalizationEigensolver,
) where {B,T,R}
    basis = scfres.basis

    ψ_virt = Matrix{T}[]
    eigenvalues_virt = Vector{R}[]
    occupations_virt = Vector{R}[]

    for ik = 1:length(basis.kpoints)
        kpt = basis.kpoints[ik]
        ψocck = occ_space.ψ[ik]
        Nfull = length(kpt.G_vectors)
        N_virt = target.n_orbitals === :all ? (Nfull - size(ψocck, 2)) : target.n_orbitals

        # Shift the occupied states using LevelShiftedOperator, identical to LOBPCGEigensolver
        ε_homo = maximum(occ_space.eigenvalues[ik])
        Ecut = basis.Ecut
        ham_hf_levelshifted = LevelShiftedOperator(scfres.ham[ik], ψocck, ε_homo, 1e-5, 2 * Ecut)

        # Build the full identity matrix for this k-point space to compute the dense Hamiltonian
        I_mat = Matrix{T}(I, Nfull, Nfull)
        H_dense = Hermitian(Matrix(ham_hf_levelshifted * I_mat))

        # Diagonalize the full matrix
        eigen_res = eigen(H_dense)

        # Extract the lowest N_virt eigenvalues and eigenvectors (the true virtual states)
        indices = 1:N_virt

        # Shift the eigenvalues back since LevelShiftedOperator shifts the whole spectrum
        ε_virtk = eigen_res.values[indices] .+ ham_hf_levelshifted.ε_homo .- ham_hf_levelshifted.safe_shift

        push!(ψ_virt, eigen_res.vectors[:, indices])
        push!(eigenvalues_virt, ε_virtk)
        push!(occupations_virt, zeros(R, N_virt))
    end

    return OrbitalSpace{B,T,R}(
        basis,
        ψ_virt,
        eigenvalues_virt,
        occupations_virt,
        occ_space.εF,
    )
end

# --- Fallbacks ---

function generate_orbitals(target::DensitySpecificVirtuals, scfres, occ_space)
    generate_orbitals(target, scfres, occ_space, LOBPCGEigensolver())
end

function generate_orbitals(target::CanonicalVirtuals, scfres, occ_space)
    solver = target.n_orbitals === :all ? FullDiagonalizationEigensolver() : LOBPCGEigensolver()
    generate_orbitals(target, scfres, occ_space, solver)
end
