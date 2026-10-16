@doc raw"""
    LOBPCG(; tol=1e-6, maxiter=200, extra_fraction=0.03, n_extra_min=13,
           callback=DefaultLobpcgCallback())

Iterative eigensolver (`lobpcg` from LOBPCGEigensolver.jl) computing only the requested
number of lowest eigenpairs. Suited when `n_orbitals` is small compared to the number of
plane waves.

To avoid stalling on a (near-)degenerate multiplet cut by the highest requested eigenpair,
the block is enlarged by ``\max(n_\text{extra,min}, \lceil f_\text{extra} n \rceil)``
vectors, but only the lowest ``n`` eigenpairs must converge (`n_conv_check`) and are
returned.

# Fields
- `tol::Float64`: the residual tolerance for convergence
- `maxiter::Int`: maximum number of iterations before the solver aborts
- `extra_fraction::Float64`: fraction ``f_\text{extra}`` of extra block vectors
- `n_extra_min::Int`: minimal number ``n_\text{extra,min}`` of extra block vectors
- `callback`: called with the solver state after each iteration; the default prints the
  convergence progress, `callback=identity` silences it

# Cost
Per iteration, with ``n`` the block size (requested eigenpairs plus extra vectors):
- time: ``n`` applications of the operators (see the target's docstring) plus
  ``O(N_\text{pw} n^2 + n^3)`` for the Rayleigh-Ritz step
- memory: ``O(N_\text{pw} n)`` complex numbers for a few blocks of ``n`` vectors
"""
Base.@kwdef struct LOBPCG
    tol::Float64 = 1e-6
    maxiter::Int = 200
    extra_fraction::Float64 = 0.03
    n_extra_min::Int = 13
    callback = DefaultLobpcgCallback()
end

@doc raw"""
    FullDiagonalization()

Dense exact diagonalization (`LinearAlgebra.eigen`) of the operator in the full
plane-wave basis. Yields all eigenpairs at once and is preferable to
[`LOBPCG`](@ref) when a large fraction of the spectrum is requested,
e.g. `CanonicalVirtuals(n_orbitals=:all)`.

# Cost
- time: ``N_\text{pw}`` applications of the operators to build them as dense matrices
  (see the target's docstring), plus ``O(N_\text{pw}^3)`` for the diagonalization
- memory: ``O(N_\text{pw}^2)`` complex numbers for the dense matrices
"""
struct FullDiagonalization end

include("davidson.jl")

"""
    _n_guess(solver, n)

Number of guess vectors `solver` wants in `X0` for computing the lowest `n` eigenpairs in
[`_solve`](@ref). Defaults to `n`; a solver converging better with a larger block (e.g.
[`LOBPCG`](@ref)) returns more.
"""
_n_guess(solver, n) = n
function _n_guess(solver::LOBPCG, n)
    return n + max(solver.n_extra_min, ceil(Int, solver.extra_fraction * n))
end

"""
    _solve(A, B, X0, n, preconditioner, solver) -> (; λ, X)

Extension point for eigensolvers: the lowest `n` eigenpairs of `A φ = λ B φ` (`B` may be
`I`) with the eigenvalues `λ` and the eigenvectors as columns of `X`, starting from the
guess `X0` with `size(X0, 2) ≥ n` columns (see [`_n_guess`](@ref)). The eigenvalue problems
come from [`_eigenproblems`](@ref). Each solver documents its cost per operator application
in its docstring.
"""
function _solve end

function _solve(A, B, X0, n, preconditioner, solver::LOBPCG)
    res = lobpcg(A, X0, B, preconditioner, solver.tol, solver.maxiter;
                 n_conv_check=n, solver.callback)
    return (; λ=res.λ[1:n], X=res.X[:, 1:n])
end

function _solve(A, B, X0, n, preconditioner, ::FullDiagonalization)
    N = size(X0, 1)
    T = eltype(X0)
    A_dense = Hermitian(Matrix(A * Matrix{T}(I, N, N)))
    if B === I
        res = eigen(A_dense)
    else
        res = eigen(A_dense, Hermitian(Matrix(B * Matrix{T}(I, N, N))))
    end
    return (; λ=res.values[1:n], X=res.vectors[:, 1:n])
end

function _solve(A, B, X0, n, preconditioner, solver::BlockDavidson)
    return davidson(A, X0, solver)
end
