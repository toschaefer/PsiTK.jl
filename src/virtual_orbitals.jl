function construct_stochastic_orbitals(N, kpt, orbitalType)
    Npw = length(kpt.G_vectors)
    radius = rand(Npw, N)
    phase = cis.(2π .* rand(Npw, N))
    ϕk = zeros(orbitalType, Npw, N)
    ϕk .= radius .* phase
    for a in 1:N
        ϕk[:, a] ./= norm(ϕk[:, a])
    end
    qr_decomp = qr(ϕk)
    return Matrix(qr_decomp.Q)
end

# --- Targets: which eigenvalue problem defines the virtual orbitals (physics) ---

@doc raw"""
    DensitySpecificVirtuals(; n_orbitals, interaction_kernel)

Target for [`generate_orbitals`](@ref): compressed virtual orbitals (Density Specific
Virtuals) from the lowest eigenpairs of the generalized eigenvalue problem in the virtual
space
```math
\mathcal K \varphi  =  \lambda h \varphi
```
where ``\mathcal K`` is the exchange operator of the active occupied orbitals with the
DFTK `InteractionKernel` `interaction_kernel` (at full strength and uncompressed) and
``h`` the mean-field Hamiltonian of the reference (the Fock operator for Hartree-Fock, the
Kohn-Sham Hamiltonian for DFT). The kernel has no default, since the orbitals depend on it
(see [`generate_orbitals`](@ref) for taking it from the SCF).

The generated orbitals are NOT orthonormal (they are ``h``-orthonormal), and their
`eigenvalues` are the generalized Rayleigh quotients ``\lambda_i``, not orbital energies.
Use [`canonicalize_orbitals`](@ref) to obtain orthonormal orbitals with energies of ``h``.

# Cost
Per operator application to one vector, on top of the solver's cost:
- ``\mathcal K``: ``O(N_k N_\text{act} N_r \log N_r)``, two FFTs per active orbital and
  k-point
- ``h``: ``O(N_k N_\text{occ} N_r \log N_r)`` with exact exchange (two FFTs per occupied
  orbital and k-point), ``O(N_r \log N_r)`` without
"""
Base.@kwdef struct DensitySpecificVirtuals{K}
    n_orbitals::Int
    interaction_kernel::K
end

@doc raw"""
    CanonicalVirtuals(; n_orbitals=:all)

Target for [`generate_orbitals`](@ref): canonical virtual orbitals, i.e. the lowest
eigenpairs of the mean-field Hamiltonian ``h`` of the reference (the Fock operator for
Hartree-Fock, the Kohn-Sham Hamiltonian for DFT) in the virtual space. `n_orbitals=:all`
yields the complete virtual plane-wave space.

# Cost
Per operator application to one vector, on top of the solver's cost:
``O(N_k N_\text{occ} N_r \log N_r)`` for ``h`` with exact exchange (two FFTs per occupied
orbital and k-point), ``O(N_r \log N_r)`` without.
"""
Base.@kwdef struct CanonicalVirtuals
    n_orbitals::Union{Int,Symbol} = :all
end

@doc raw"""
    MaximalExchangeVirtuals(; n_orbitals, interaction_kernel)

Target for [`generate_orbitals`](@ref): virtual orbitals maximizing the exchange
interaction with the active occupied orbitals, i.e. the lowest (most negative) eigenpairs
of
```math
\mathcal K \varphi  =  \lambda \varphi
```
where ``\mathcal K`` is the exchange operator of the active occupied orbitals with the
DFTK `InteractionKernel` `interaction_kernel` (at full strength and uncompressed). The
kernel has no default, since the orbitals depend on it (see [`generate_orbitals`](@ref)
for taking it from the SCF).

# Cost
Per operator application to one vector, on top of the solver's cost:
``O(N_k N_\text{act} N_r \log N_r)`` for ``\mathcal K``, two FFTs per active orbital and
k-point.
"""
Base.@kwdef struct MaximalExchangeVirtuals{K}
    n_orbitals::Int
    interaction_kernel::K
end

"""
Union of the virtual-orbital targets accepted by [`generate_orbitals`](@ref). A new target
is a configuration struct with an `n_orbitals` field, added here, with methods of
[`_eigenproblems`](@ref) and [`_is_orthonormal`](@ref). It documents the cost of one
application of its operators in its docstring.
"""
const VirtualOrbitalTarget =
    Union{DensitySpecificVirtuals,CanonicalVirtuals,MaximalExchangeVirtuals}

# --- Eigenvalue problems ---

"""
    _eigenproblems(target, occ_space, active, ham)

Extension point for virtual-orbital targets: the eigenvalue problem that defines the
virtual orbitals of `target`, as a vector over k-points of `(; A, B, ε_offset)`. The lowest
eigenpairs of `A φ = λ B φ` (`B = I` for a standard eigenvalue problem) are the wanted
orbitals with eigenvalues `λ + ε_offset`. The eigenvalue problem is solved by
[`_solve`](@ref), so targets and eigensolvers combine freely (see
[`VirtualOrbitalTarget`](@ref) for adding a target).

The arguments are checked by [`generate_orbitals`](@ref):
- `occ_space`: all occupied orbitals of the reference, which the operators must push out
  of the low end of the spectrum
- `active`: the occupied orbitals to correlate, a subset of `occ_space`, which define
  exchange-based targets
- `ham`: the mean-field Hamiltonian of `occ_space`, with an uncompressed exchange operator
"""
function _eigenproblems end

# Mean-field Hamiltonian restricted to the virtual space, with the occupied space pushed
# above the virtual spectrum
function _projected_hamiltonian(ham_k, ψocc_k, εocc_k)
    ε_homo = maximum(εocc_k)
    safe_shift = 1e-5              # keeps the shifted virtual spectrum strictly positive
    penalty = 2 * ham_k.basis.Ecut # lifts the occupied space above all plane-wave energies
    op = ProjectedShiftedOperator(ham_k, ψocc_k, safe_shift - ε_homo, penalty)
    return op, ε_homo - safe_shift
end

# Exchange operator of the active orbitals restricted to the virtual space, one block per
# k-point: the occupied space gets a positive penalty, so it is never among the lowest
# (negative) exchange eigenvalues. The term is built at full strength (scaling 1) and
# applied uncompressed, independent of the exchange term of the SCF model.
function _projected_exchange(interaction_kernel, occ_space, active)
    basis = occ_space.basis
    term = DFTK.ExactExchange(; kernel=interaction_kernel)(basis)
    _, K = DFTK.ene_ops(
        term,
        basis,
        active.ψ,
        active.occupation;
        exxalg=DFTK.VanillaExx(),
    )
    penalty = abs(minimum(minimum, occ_space.eigenvalues)) + 2.0
    return map(eachindex(K)) do ik
        ProjectedShiftedOperator(K[ik], occ_space.ψ[ik], 0.0, penalty)
    end
end

function _eigenproblems(::CanonicalVirtuals, occ_space, active, ham)
    (; ψ, eigenvalues) = occ_space
    return map(eachindex(ham.basis.kpoints)) do ik
        A, ε_offset = _projected_hamiltonian(ham[ik], ψ[ik], eigenvalues[ik])
        (; A, B=I, ε_offset)
    end
end

function _eigenproblems(target::DensitySpecificVirtuals, occ_space, active, ham)
    (; ψ, eigenvalues) = occ_space
    K = _projected_exchange(target.interaction_kernel, occ_space, active)
    return map(eachindex(ham.basis.kpoints)) do ik
        B, ε_offset = _projected_hamiltonian(ham[ik], ψ[ik], eigenvalues[ik])
        # eigenvalues are Rayleigh quotients λ = <φ|K|φ>/<φ|h|φ>, reported as they are
        (; A=K[ik], B, ε_offset=zero(ε_offset))
    end
end

function _eigenproblems(target::MaximalExchangeVirtuals, occ_space, active, ham)
    K = _projected_exchange(target.interaction_kernel, occ_space, active)
    return map(eachindex(ham.basis.kpoints)) do ik
        (; A=K[ik], B=I, ε_offset=zero(eltype(occ_space.eigenvalues[ik])))
    end
end

function _n_orbitals(target::VirtualOrbitalTarget, Npw, Nocc)
    target.n_orbitals === :all && return Npw - Nocc
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

# The active orbitals must lie in the span of the orthonormal occupied orbitals
function _check_active_space(active, occ_space; tol=1e-8)
    active.basis === occ_space.basis ||
        error("`active` and `occ_space` must be defined on the same basis.")
    for (Xa, Xo) in zip(active.ψ, occ_space.ψ)
        deviation = norm(Xa - Xo * (Xo' * Xa))
        deviation < tol * max(1, norm(Xa)) ||
            error("The active orbitals do not lie in the occupied space (deviation " *
                  "$deviation).")
    end
    return nothing
end

@doc raw"""
    generate_orbitals(target, occ_space::OrbitalSpace, ham; active=occ_space, solver)

Generate virtual orbitals orthogonal to the occupied space `occ_space` by solving the
eigenvalue problem defined by `target` with the eigensolver `solver`.

# Arguments
- `target`: [`CanonicalVirtuals`](@ref), [`DensitySpecificVirtuals`](@ref) or
  [`MaximalExchangeVirtuals`](@ref)
- `occ_space`: all occupied orbitals of the reference, which are projected out
- `ham`: the DFTK mean-field Hamiltonian of `occ_space` (the Fock operator for
  Hartree-Fock, the Kohn-Sham Hamiltonian for DFT), e.g. `scfres.ham`. Its exchange
  operator, if any, must be uncompressed: an ACE exchange (`exxalg=AceExx()`) is exact
  only on the orbitals it was sketched with and is rejected with an error
- `active`: the occupied orbitals to correlate, a subset of `occ_space` (e.g. via
  [`select_orbitals`](@ref) for embedding). They define the exchange operator
  ``\mathcal K`` of [`DensitySpecificVirtuals`](@ref) and
  [`MaximalExchangeVirtuals`](@ref); [`CanonicalVirtuals`](@ref) ignore them
- `solver`: [`LOBPCG`](@ref), [`FullDiagonalization`](@ref) or [`BlockDavidson`](@ref).
  Defaults to `LOBPCG()`, or `FullDiagonalization()` for
  `CanonicalVirtuals(n_orbitals=:all)`

The interaction kernel of ``\mathcal K`` is set by `target`. To use the one of a
Hartree-Fock or hybrid SCF, take it from the model; the exact-exchange fraction of a
hybrid is stored separately and not part of the kernel:
```julia
kernel = only(t.kernel for t in scfres.basis.model.term_types if t isa ExactExchange)
target = DensitySpecificVirtuals(; n_orbitals=50, interaction_kernel=kernel)
```

# Returns
An `OrbitalSpace` with `target.n_orbitals` orbitals per k-point.

# Cost
Per k-point, with ``N_\text{virt}`` = `target.n_orbitals`: that of `solver` (see its
docstring), with the operator applications of `target` (see its docstring), plus
``O(N_\text{pw} N_\text{virt}^2)`` time for the initial guess. Exchange-based targets
build their exchange operator once for all k-points:
- time: ``O(N_k^2 N_\text{act}^2 N_r \log N_r)`` for the exchange energy DFTK evaluates
  along
- memory: ``N_k N_\text{act} N_r`` complex numbers for the active orbitals in real space
"""
function generate_orbitals(
    target::VirtualOrbitalTarget,
    occ_space::OrbitalSpace{TB,T,R},
    ham;
    active::OrbitalSpace=occ_space,
    solver=_default_solver(target),
) where {TB,T,R}
    _check_uncompressed_exchange(ham)
    ham.basis === occ_space.basis ||
        error("`ham` and `occ_space` must be defined on the same basis.")
    _check_active_space(active, occ_space)

    basis = ham.basis
    problems = _eigenproblems(target, occ_space, active, ham)

    ψ_virt = Matrix{T}[]
    eigenvalues_virt = Vector{R}[]
    occupation_virt = Vector{R}[]

    for (ik, kpt) in enumerate(basis.kpoints)
        ψocck = occ_space.ψ[ik]
        Npw = length(kpt.G_vectors)
        n_orbitals = _n_orbitals(target, Npw, size(ψocck, 2))

        if solver isa LOBPCG && n_orbitals > 0.1 * Npw
            @warn "n_orbitals ($n_orbitals) is > 10% of plane waves ($Npw). " *
                  "FullDiagonalization() might be faster."
        end

        # Stochastic initial guess in the orthogonal complement of the occupied space
        n_guess = min(_n_guess(solver, n_orbitals), Npw - size(ψocck, 2))
        X0 = construct_stochastic_orbitals(n_guess, kpt, T)
        X0 .-= ψocck * (ψocck' * X0)
        X0 = Matrix(qr(X0).Q)

        preconditioner = DFTK.PreconditionerTPA(basis, kpt)
        (; A, B, ε_offset) = problems[ik]
        (; λ, X) = _solve(A, B, X0, n_orbitals, preconditioner, solver)

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
