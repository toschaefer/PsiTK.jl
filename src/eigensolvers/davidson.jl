"""
    BlockDavidsonEigensolver(; tol=1e-6, use_jabobi_davidson=false, maxiter=200)

Configuration struct for the Block Davidson iterative eigensolver.

# Fields
- `tol::Float64`: The residual tolerance for convergence.
- `use_jabobi_davidson::Bool`: If true, applies a Jacobi-Davidson style correction equation... 
- `maxiter::Int`: Maximum number of iterations before the solver aborts.
"""
Base.@kwdef struct BlockDavidsonEigensolver
    tol::Float64 = 1e-6
    use_jabobi_davidson::Bool = false
    maxiter::Int = 200
end

"""
    davidson(A, X0, solver::BlockDavidsonEigensolver; kwargs...)

Iteratively computes the lowest eigenvalues and eigenvectors of the linear operator `A` 
using a Block Davidson method.

TODO: some more details would be nice

References:
Francesco Mereto, Master's Thesis, 2026 (link?)

# Arguments
- `A`: The linear operator or matrix to diagonalize.
- `X0`: The initial guess for the eigenvectors.
- `solver`: The `BlockDavidsonEigensolver` configuration (tolerance, max iterations, etc.).
"""
function davidson(A, X0, solver::BlockDavidsonEigensolver; kwargs...)
    # TODO: Mereto
    error("Block Davidson solver not implemented yet.")
    
    # Expected to return a named tuple or struct, for example:
    # return (λ = eigenvalues, X = eigenvectors, residual = residuals)
end
