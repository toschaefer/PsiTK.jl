export LOBPCGEigensolver, FullDiagonalizationEigensolver

# Wrapper for DFTK's LOBPCG
Base.@kwdef struct LOBPCGEigensolver
    tol::Float64 = 1e-6
    maxiter::Int = 200
end

# Wrapper for LinearAlgebra.eigen
struct FullDiagonalizationEigensolver end

# ---------------------------------------------------------
# FUTURE SOLVERS (e.g. Davidson)
# ---------------------------------------------------------
# export DavidsonEigensolver
# include("davidson.jl")
