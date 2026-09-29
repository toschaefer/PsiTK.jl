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

# --- Targets: which eigenvalue problem defines the virtual orbitals (physics) ---

@doc raw"""
    DensitySpecificVirtuals(; n_orbitals)

Target for [`generate_orbitals`](@ref): compressed virtual orbitals (Density Specific
Virtuals) from the lowest eigenpairs of the generalized eigenvalue problem in the virtual
space
```math
\mathcal K \varphi  =  \lambda h \varphi
```
where ``\mathcal K`` and ``h`` are the Fock exchange operator and the Fock Hamiltonian,
respectively.

The generated orbitals are NOT orthonormal (they are ``h``-orthonormal), and their
`eigenvalues` are the generalized Rayleigh quotients ``\lambda_i``, not orbital energies.
Use [`canonicalize_orbitals`](@ref) to obtain orthonormal orbitals with Fock energies.
"""
Base.@kwdef struct DensitySpecificVirtuals
    n_orbitals::Int
end

"""
    CanonicalVirtuals(; n_orbitals=:all)

Target for [`generate_orbitals`](@ref): canonical virtual orbitals, i.e. the lowest
eigenpairs of the Fock Hamiltonian in the virtual space. `n_orbitals=:all` yields the
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
where ``\mathcal K`` is the Fock exchange operator.
"""
Base.@kwdef struct MaximalExchangeVirtuals
    n_orbitals::Int
end

"""
Union of the virtual-orbital targets accepted by [`generate_orbitals`](@ref). A new target
is added here and implements [`_eigenproblems`](@ref) and [`_is_orthonormal`](@ref).
"""
const VirtualOrbitalTarget =
    Union{DensitySpecificVirtuals,CanonicalVirtuals,MaximalExchangeVirtuals}

# --- Eigenvalue problems ---

"""
    _eigenproblems(target, occ_space, ham)

Extension point for virtual-orbital targets: the eigenvalue problem that defines the
virtual orbitals of `target`, as a vector over k-points of `(; A, B, ε_offset)`. The lowest
eigenpairs of `A φ = λ B φ` (`B = I` for a standard eigenvalue problem) are the wanted
orbitals with eigenvalues `λ + ε_offset`. The eigenvalue problem is solved by
[`_solve`](@ref), so targets and eigensolvers combine freely.
"""
function _eigenproblems end

# Fock operator with the occupied space pushed above the virtual spectrum
function _levelshifted_fock(ham_k, ψocc_k, εocc_k)
    ε_homo = maximum(εocc_k)
    safe_shift = 1e-5              # keeps the shifted virtual spectrum strictly positive
    penalty = 2 * ham_k.basis.Ecut # lifts the occupied space above all plane-wave energies
    op = LevelShiftedOperator(ham_k, ψocc_k, ε_homo, safe_shift, penalty)
    return op, ε_homo - safe_shift
end

# Fock exchange operator of the occupied orbitals restricted to the virtual space, one block
# per k-point: occupied components are shifted to positive energies, so they are never
# among the lowest (negative) exchange eigenvalues
function _projected_exchange(basis, ψocc, occupation, εocc)
    term = only(t for t in basis.terms if t isa DFTK.TermExactExchange)
    _, K = DFTK.ene_ops(term, basis, ψocc, occupation)
    shift = abs(minimum(minimum, εocc)) + 2.0
    return [ProjectedShiftedOperator(K[ik], ψocc[ik], shift) for ik in eachindex(K)]
end

function _eigenproblems(::CanonicalVirtuals, occ_space, ham)
    (; ψ, eigenvalues) = occ_space
    return map(eachindex(ham.basis.kpoints)) do ik
        A, ε_offset = _levelshifted_fock(ham[ik], ψ[ik], eigenvalues[ik])
        (; A, B=I, ε_offset)
    end
end

function _eigenproblems(::DensitySpecificVirtuals, occ_space, ham)
    (; ψ, eigenvalues, occupation) = occ_space
    K = _projected_exchange(ham.basis, ψ, occupation, eigenvalues)
    return map(eachindex(ham.basis.kpoints)) do ik
        B, ε_offset = _levelshifted_fock(ham[ik], ψ[ik], eigenvalues[ik])
        # eigenvalues are Rayleigh quotients λ = <φ|K|φ>/<φ|h|φ>, reported as they are
        (; A=K[ik], B, ε_offset=zero(ε_offset))
    end
end

function _eigenproblems(::MaximalExchangeVirtuals, occ_space, ham)
    (; ψ, eigenvalues, occupation) = occ_space
    K = _projected_exchange(ham.basis, ψ, occupation, eigenvalues)
    return map(eachindex(ham.basis.kpoints)) do ik
        (; A=K[ik], B=I, ε_offset=zero(eltype(eigenvalues[ik])))
    end
end

function _n_orbitals(target::VirtualOrbitalTarget, n_G, n_occ)
    target.n_orbitals === :all && return n_G - n_occ
    return target.n_orbitals
end

"""
    _is_orthonormal(target)

Whether the orbitals generated for `target` are orthonormal. Solutions of a generalized
eigenvalue problem (`B ≠ I` in [`_eigenproblems`](@ref)) are only `B`-orthonormal.
"""
function _is_orthonormal end

_is_orthonormal(::CanonicalVirtuals) = true
_is_orthonormal(::DensitySpecificVirtuals) = false
_is_orthonormal(::MaximalExchangeVirtuals) = true

_default_solver(::VirtualOrbitalTarget) = LOBPCG()
function _default_solver(target::CanonicalVirtuals)
    return target.n_orbitals === :all ? FullDiagonalization() : LOBPCG()
end

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
  Defaults to `LOBPCG()`, or `FullDiagonalization()` for
  `CanonicalVirtuals(n_orbitals=:all)`

# Returns
An `OrbitalSpace` with `target.n_orbitals` orbitals per k-point.
"""
function generate_orbitals(
    target::VirtualOrbitalTarget,
    occ_space::OrbitalSpace{TB,T,R},
    ham;
    solver=_default_solver(target),
) where {TB,T,R}
    basis = ham.basis
    problems = _eigenproblems(target, occ_space, ham)

    ψ_virt = Matrix{T}[]
    eigenvalues_virt = Vector{R}[]
    occupation_virt = Vector{R}[]

    for (ik, kpt) in enumerate(basis.kpoints)
        ψocck = occ_space.ψ[ik]
        Nfull = length(kpt.G_vectors)
        n_orbitals = _n_orbitals(target, Nfull, size(ψocck, 2))

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
        push!(occupation_virt, zeros(R, n_orbitals))
    end

    return OrbitalSpace{TB,T,R}(
        basis,
        ψ_virt,
        eigenvalues_virt,
        occupation_virt,
        occ_space.εF,
        _is_orthonormal(target),
    )
end
