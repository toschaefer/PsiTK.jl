export LOBPCGEigensolver, FullDiagonalizationEigensolver, BlockDavidson

# Wrapper for DFTK's LOBPCG
Base.@kwdef struct LOBPCGEigensolver
    tol::Float64 = 1e-6
    maxiter::Int = 200
end

# Wrapper for LinearAlgebra.eigen
struct FullDiagonalizationEigensolver end

Base.@kwdef struct BlockDavidson
    tol::Float64 = 1e-6
    use_jabobi_davidson::Bool = false
    maxiter::Int = 200
    # TODO: add further arguments here
end

include("davidson.jl")
