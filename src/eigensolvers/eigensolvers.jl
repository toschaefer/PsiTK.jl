@doc raw"""
    LOBPCG(; tol=1e-6, maxiter=200, callback=DefaultLobpcgCallback())

Iterative eigensolver (`lobpcg` from LOBPCGEigensolver.jl) computing only the requested
number of lowest eigenpairs. Suited when `n_orbitals` is small compared to the number of
plane waves.

# Fields
- `tol::Float64`: the residual tolerance for convergence
- `maxiter::Int`: maximum number of iterations before the solver aborts
- `callback`: called with the solver state after each iteration; the default prints the
  convergence progress, `callback=identity` silences it

# Cost
Per iteration, with ``n`` the number of requested eigenpairs:
- time: ``n`` applications of the operators (see the target's docstring) plus
  ``O(N_\text{pw} n^2 + n^3)`` for the Rayleigh-Ritz step
- memory: ``O(N_\text{pw} n)`` complex numbers for a few blocks of ``n`` vectors
"""
Base.@kwdef struct LOBPCG
    tol::Float64 = 1e-6
    maxiter::Int = 200
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
    _solve(A, B, X0, preconditioner, solver) -> (; λ, X)

Extension point for eigensolvers: the lowest `size(X0, 2)` eigenpairs of `A φ = λ B φ`
(`B` may be `I`) with the eigenvalues `λ` and the eigenvectors as columns of `X`, starting
from the guess `X0`. The eigenvalue problems come from [`_eigenproblems`](@ref). Each
solver documents its cost per operator application in its docstring.
"""
function _solve end

function _solve(A, B, X0, preconditioner, solver::LOBPCG)
    res = lobpcg(A, X0, B, preconditioner, solver.tol, solver.maxiter; solver.callback)
    return (; λ=res.λ, X=res.X)
end

function _solve(A, B, X0, preconditioner, ::FullDiagonalization)
    N, n = size(X0)
    T = eltype(X0)
    A_dense = Hermitian(Matrix(A * Matrix{T}(I, N, N)))
    if B === I
        res = eigen(A_dense)
    else
        res = eigen(A_dense, Hermitian(Matrix(B * Matrix{T}(I, N, N))))
    end
    return (; λ=res.values[1:n], X=res.vectors[:, 1:n])
end

function _solve(A, B, X0, preconditioner, solver::BlockDavidson)
    return davidson(A, X0, solver)
end
