"""
    davidson(A, X0, solver::BlockDavidson; kwargs...)

Iteratively computes the lowest eigenvalues and eigenvectors of the linear operator `A` 
using a Block Davidson method.

TODO: some more details would be nice

References:
Francesco Mereto, Master's Thesis, 2026 (link?)

# Arguments
- `A`: The linear operator or matrix to diagonalize.
- `X0`: The initial guess for the eigenvectors.
- `solver`: The `BlockDavidson` configuration (tolerance, max iterations, etc.).
"""
function davidson(A, X0, solver::BlockDavidson; kwargs...)
    # TODO: Mereto
    error("Block Davidson solver not implemented yet.")
    
    # Expected to return a named tuple or struct, for example:
    # return (λ = eigenvalues, X = eigenvectors, residual = residuals)
end
