export LOBPCG, FullDiagonalization, BlockDavidson

"""
    LOBPCG(; tol=1e-6, maxiter=200)

Iterative eigensolver (`lobpcg` from LOBPCGEigensolver.jl) computing only the requested
number of lowest eigenpairs. Suited when `n_orbitals` is small compared to the number of plane waves.

# Fields
- `tol::Float64`: the residual tolerance for convergence
- `maxiter::Int`: maximum number of iterations before the solver aborts
"""
Base.@kwdef struct LOBPCG
    tol::Float64 = 1e-6
    maxiter::Int = 200
end

"""
    FullDiagonalization()

Dense exact diagonalization (`LinearAlgebra.eigen`) of the operator in the full
plane-wave basis. Yields all eigenpairs at once and is preferable to
[`LOBPCG`](@ref) when a large fraction of the spectrum is requested,
e.g. `CanonicalVirtuals(n_orbitals=:all)`.
"""
struct FullDiagonalization end

include("davidson.jl")
