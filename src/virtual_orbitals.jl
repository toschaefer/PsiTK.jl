export DensitySpecificVirtuals, CanonicalVirtuals, MaximalExchangeVirtuals
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

# --- Targets: which eigenvalue problem defines the virtual orbitals (physics) ---

@doc raw"""
    DensitySpecificVirtuals(; n_orbitals)

Target for [`generate_orbitals`](@ref): compressed virtual orbitals (Density Specific
Virtuals) from the lowest eigenpairs of the generalized eigenvalue problem in the virtual space
```math
\mathcal K \varphi  =  \lambda h \varphi
```
where $\mathcal K$ and $h$ are the Fock exchange operator and the Fock Hamiltonian, respectively.

The generated orbitals are NOT orthonormal (they are $h$-orthonormal), and their `eigenvalues`
are the generalized Rayleigh quotients $\lambda_i$, not orbital energies. Use
[`canonicalize_orbitals`](@ref) to obtain orthonormal orbitals with Fock energies.
"""
Base.@kwdef struct DensitySpecificVirtuals
    n_orbitals::Int
end

"""
    CanonicalVirtuals(; n_orbitals=:all)

Target for [`generate_orbitals`](@ref): canonical virtual orbitals, i.e. the lowest
eigenpairs of the Fock Hamiltonian in the virtual space. `n_orbitals = :all` yields the
complete virtual plane-wave space.
"""
Base.@kwdef struct CanonicalVirtuals
    n_orbitals::Union{Int,Symbol} = :all
end

@doc raw"""
    MaximalExchangeVirtuals(; n_orbitals)

Target for [`generate_orbitals`](@ref): virtual orbitals maximizing the exchange
interaction with the occupied space, i.e. the lowest (most negative) eigenpairs of
```math
\mathcal K \varphi  =  \lambda \varphi
```
where $\mathcal K$ is the Fock exchange operator.
"""
Base.@kwdef struct MaximalExchangeVirtuals
    n_orbitals::Int
end

const VirtualOrbitalTarget =
    Union{DensitySpecificVirtuals,CanonicalVirtuals,MaximalExchangeVirtuals}

# --- Eigenvalue problems: per k-point (; A, B, ε_offset) such that the lowest eigenpairs
#     of A φ = λ B φ are the wanted orbitals with eigenvalue λ + ε_offset ---

# Fock operator with the occupied space pushed above the virtual spectrum
function _levelshifted_fock(ham, occ_space, ik)
    ε_homo = maximum(occ_space.eigenvalues[ik])
    safe_shift = 1e-5            # keeps the shifted virtual spectrum strictly positive
    penalty = 2 * ham.basis.Ecut # lifts the occupied space above all plane-wave energies
    op = LevelShiftedOperator(ham[ik], occ_space.ψ[ik], ε_homo, safe_shift, penalty)
    return op, ε_homo - safe_shift
end

# Fock exchange operator of the occupied orbitals, one block per k-point
function _exchange_operator(ham, occ_space)
    basis = ham.basis
    term = only(t for t in basis.terms if t isa DFTK.TermExactExchange)
    _, K = DFTK.ene_ops(term, basis, occ_space.ψ, occ_space.occupations)
    return K
end

# Exchange operator restricted to the virtual space: occupied components are shifted to
# positive energies, so they are never among the lowest (negative) exchange eigenvalues
function _projected_exchange(K, occ_space, ik)
    shift = abs(minimum(minimum.(occ_space.eigenvalues))) + 2.0
    return ProjectedShiftedOperator(K[ik], occ_space.ψ[ik], shift)
end

function _eigenproblems(::CanonicalVirtuals, occ_space, ham)
    return map(eachindex(ham.basis.kpoints)) do ik
        A, ε_offset = _levelshifted_fock(ham, occ_space, ik)
        (; A, B = I, ε_offset)
    end
end

function _eigenproblems(::DensitySpecificVirtuals, occ_space, ham)
    K = _exchange_operator(ham, occ_space)
    return map(eachindex(ham.basis.kpoints)) do ik
        B, ε_offset = _levelshifted_fock(ham, occ_space, ik)
        # eigenvalues are Rayleigh quotients λ = <φ|K|φ>/<φ|h|φ>, reported as they are
        (; A = _projected_exchange(K, occ_space, ik), B, ε_offset = zero(ε_offset))
    end
end

function _eigenproblems(::MaximalExchangeVirtuals, occ_space, ham)
    K = _exchange_operator(ham, occ_space)
    return map(eachindex(ham.basis.kpoints)) do ik
        A = _projected_exchange(K, occ_space, ik)
        (; A, B = I, ε_offset = zero(eltype(occ_space.eigenvalues[ik])))
    end
end

function _n_orbitals(target::VirtualOrbitalTarget, occ_space, ik)
    if target.n_orbitals === :all
        return length(occ_space.basis.kpoints[ik].G_vectors) - size(occ_space.ψ[ik], 2)
    end
    return target.n_orbitals
end

# Generalized eigenvectors (B ≠ I) are only B-orthonormal
_is_orthonormal(::CanonicalVirtuals) = true
_is_orthonormal(::DensitySpecificVirtuals) = false
_is_orthonormal(::MaximalExchangeVirtuals) = true

_default_solver(::VirtualOrbitalTarget) = LOBPCG()
_default_solver(target::CanonicalVirtuals) =
    target.n_orbitals === :all ? FullDiagonalization() : LOBPCG()

# --- Generator ---

"""
    generate_orbitals(target, occ_space::OrbitalSpace, ham; solver=LOBPCG())

Generate virtual orbitals orthogonal to the occupied space `occ_space` by solving the
eigenvalue problem defined by `target` with the eigensolver `solver`.

# Arguments
- `target`: [`CanonicalVirtuals`](@ref), [`DensitySpecificVirtuals`](@ref) or
  [`MaximalExchangeVirtuals`](@ref)
- `occ_space`: the occupied orbitals (define the exchange operator and are projected out)
- `ham`: the DFTK Hamiltonian of the converged Hartree-Fock calculation (`scfres.ham`)
- `solver`: [`LOBPCG`](@ref), [`FullDiagonalization`](@ref) or [`BlockDavidson`](@ref).
  Defaults to `LOBPCG()`, or `FullDiagonalization()` for `CanonicalVirtuals(n_orbitals=:all)`

# Returns
An `OrbitalSpace` with `target.n_orbitals` orbitals per k-point.
"""
function generate_orbitals(
    target::VirtualOrbitalTarget,
    occ_space::OrbitalSpace{TB,T,R},
    ham;
    solver = _default_solver(target),
) where {TB,T,R}
    basis = ham.basis
    problems = _eigenproblems(target, occ_space, ham)

    ψ_virt = Matrix{T}[]
    eigenvalues_virt = Vector{R}[]
    occupations_virt = Vector{R}[]

    for (ik, kpt) in enumerate(basis.kpoints)
        n_orbitals = _n_orbitals(target, occ_space, ik)
        ψocck = occ_space.ψ[ik]
        Nfull = length(kpt.G_vectors)

        if solver isa LOBPCG && n_orbitals > 0.1 * Nfull
            @warn "n_orbitals ($n_orbitals) is > 10% of plane waves ($Nfull). " *
                  "FullDiagonalization() might be faster."
        end

        # Stochastic initial guess in the orthogonal complement of the occupied space
        X0 = construct_stochastic_orbitals(n_orbitals, kpt, T)
        X0 .-= ψocck * (ψocck' * X0)
        X0 = Matrix(qr(X0).Q)

        preconditioner = DFTK.PreconditionerTPA(basis, kpt)
        (; A, B, ε_offset) = problems[ik]
        (; λ, X) = _solve(A, B, X0, preconditioner, solver)

        push!(ψ_virt, X)
        push!(eigenvalues_virt, λ .+ ε_offset)
        push!(occupations_virt, zeros(R, n_orbitals))
    end

    return OrbitalSpace{TB,T,R}(
        basis,
        ψ_virt,
        eigenvalues_virt,
        occupations_virt,
        occ_space.εF,
        _is_orthonormal(target),
    )
end
