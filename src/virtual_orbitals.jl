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

When generated, the orbitals returned ($\varphi_i$) will directly span the DSV subspace. 
These orbitals will NOT be strictly orthogonalized initially, ensuring their orbital energies strictly match the generalized Rayleigh quotient ($\lambda_i = \langle \varphi_i | \mathcal K | \varphi_i \rangle / \langle \varphi_i | \mathcal h | \varphi_i \rangle$).
They can be orthonormalized and canonicalized later using `diagonalize_orbitals`.
"""
struct DensitySpecificVirtuals{TH, TK}
    n_orbitals::Int
    ham::TH  # The DFTK Hamiltonian
    K::TK    # The Fock exchange operator
end

function DensitySpecificVirtuals(scfres, occ_space::OrbitalSpace; n_orbitals::Int)
    basis = scfres.basis
    ham = scfres.ham
    ExactExchangeTerm = only([term for term in basis.terms if term isa DFTK.TermExactExchange])
    _, K = DFTK.ene_ops(ExactExchangeTerm, basis, occ_space.ψ, occ_space.occupations)
    return DensitySpecificVirtuals(n_orbitals, ham, K)
end

"""
    CanonicalVirtuals(; n_orbitals)

Target to compute canonical virtual orbitals by diagonalizing the Fock
Hamiltonian. Set `n_orbitals = :all` to compute the full virtual plane-wave space.
"""
struct CanonicalVirtuals{TH}
    n_orbitals::Union{Int,Symbol}
    ham::TH  # The DFTK Hamiltonian
end

function CanonicalVirtuals(scfres; n_orbitals)
    return CanonicalVirtuals(n_orbitals, scfres.ham)
end

@doc raw"""
    MaximalExchangeVirtuals(; n_orbitals)

Target to compute virtual orbitals that maximize exchange, by solving the eigenvalue problem:
```math
\mathcal K \varphi  =  \lambda \varphi 
```
where $\mathcal K$ is the Fock exchange operator.
"""
struct MaximalExchangeVirtuals{TK}
    n_orbitals::Int
    K::TK    # The Fock exchange operator
end

function MaximalExchangeVirtuals(scfres, occ_space::OrbitalSpace; n_orbitals::Int)
    basis = scfres.basis
    ExactExchangeTerm = only([term for term in basis.terms if term isa DFTK.TermExactExchange])
    _, K = DFTK.ene_ops(ExactExchangeTerm, basis, occ_space.ψ, occ_space.occupations)
    return MaximalExchangeVirtuals(n_orbitals, K)
end

# --- Generators ---

"""
    generate_orbitals(target::DensitySpecificVirtuals, occ_space, solver::LOBPCGEigensolver)

Generates DSVs using an iterative LOBPCG solver.
"""
function generate_orbitals(
    target::DensitySpecificVirtuals,
    occ_space::OrbitalSpace{B,T,R},
    solver::LOBPCGEigensolver,
) where {B,T,R}
    ham = target.ham
    K = target.K
    basis = ham.basis
    Ecut = basis.Ecut

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
            LevelShiftedOperator(ham[ik], ψocck, ε_homo, 1e-5, 2 * Ecut)

        # Minimum HF eigenvalue over all k-points for safe shift
        shift = abs(minimum(minimum.(occ_space.eigenvalues))) + 2.0
        Kk_virt = ProjectedShiftedOperator(Kk, ψocck, shift)

        kinetic_preconditioner = DFTK.PreconditionerTPA(ham[ik].basis, kpt)

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

        push!(ψ_dsv, dsv.X)
        push!(eigenvalues_dsv, dsv.λ)
        push!(occupations_dsv, zeros(R, target.n_orbitals))
    end

    return OrbitalSpace{B,T,R}(basis, ψ_dsv, eigenvalues_dsv, occupations_dsv, occ_space.εF, false)
end

"""
    generate_orbitals(target::CanonicalVirtuals, occ_space, solver::LOBPCGEigensolver)

Generates canonical virtuals using an iterative LOBPCG solver on the Fock operator.
"""
function generate_orbitals(
    target::CanonicalVirtuals,
    occ_space::OrbitalSpace{B,T,R},
    solver::LOBPCGEigensolver,
) where {B,T,R}
    ham = target.ham
    basis = ham.basis
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
            LevelShiftedOperator(ham[ik], ψocck, ε_homo, 1e-5, 2 * Ecut)
        kinetic_preconditioner = DFTK.PreconditionerTPA(ham[ik].basis, kpt)

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
        true # canonical virtuals are strictly orthonormal
    )
end

"""
    generate_orbitals(target::CanonicalVirtuals, occ_space, solver::FullDiagonalizationEigensolver)

Generates canonical virtuals using full dense exact diagonalization on the Fock operator.
"""
function generate_orbitals(
    target::CanonicalVirtuals,
    occ_space::OrbitalSpace{B,T,R},
    solver::FullDiagonalizationEigensolver,
) where {B,T,R}
    ham = target.ham
    basis = ham.basis

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
        ham_hf_levelshifted = LevelShiftedOperator(ham[ik], ψocck, ε_homo, 1e-5, 2 * Ecut)

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
        true # canonical virtuals are strictly orthonormal
    )
end

"""
    generate_orbitals(target::MaximalExchangeVirtuals, occ_space, solver::BlockDavidsonEigensolver)

Generates Maximal Exchange Virtuals by solving the eigenvalue problem K φ = λ φ using a BlockDavidsonEigensolver solver.
"""
function generate_orbitals(
    target::MaximalExchangeVirtuals,
    occ_space::OrbitalSpace{B,T,R},
    solver::BlockDavidsonEigensolver,
) where {B,T,R}
    # TODO: Mereto
    error("Not implemented yet.")
end

# --- Fallbacks ---

function generate_orbitals(target::DensitySpecificVirtuals, occ_space)
    generate_orbitals(target, occ_space, LOBPCGEigensolver())
end

function generate_orbitals(target::CanonicalVirtuals, occ_space)
    solver = target.n_orbitals === :all ? FullDiagonalizationEigensolver() : LOBPCGEigensolver()
    generate_orbitals(target, occ_space, solver)
end
