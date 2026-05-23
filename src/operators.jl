using LinearAlgebra

"""
    LevelShiftedOperator(base_op, V, ε_ref, safe_shift, penalty)

A wrapper that applies a constant energy shift and a strong penalty to a projected subspace.
When applied to a vector `X`, it computes:
    (base_op - ε_ref + safe_shift) * X + penalty * V * (V' * X)

This effectively shifts the entire spectrum by `safe_shift - ε_ref`, while artificially 
pushing the eigenvalues of the subspace defined by `V` up by `penalty`. 
Useful for isolating virtual manifolds or forcing iterative eigensolvers away from an occupied subspace.
"""
struct LevelShiftedOperator{TH,TV}
    base_op::TH                      # the base operator
    V::TV                            # basis for penalty projector
    ε_homo::Float64                  # reference energy level
    safe_shift::Float64              # shift to avoid zero/negative eigenvalues
    penalty::Float64                 # penalty energy for the subspace
end

function LevelShiftedOperator(base_op, V, ε_homo, safe_shift, penalty)
    return LevelShiftedOperator{typeof(base_op),typeof(V)}(
        base_op,
        V,
        ε_homo,
        safe_shift,
        penalty,
    )
end

function LinearAlgebra.mul!(Y, op::LevelShiftedOperator, X)
    mul!(Y, op.base_op, X)
    Y .+= (-op.ε_homo + op.safe_shift) .* X

    # Projector penalty
    coeffs = op.V' * X
    mul!(Y, op.V, coeffs, op.penalty, 1.0)
    return Y
end

function Base.:*(op::LevelShiftedOperator, X::AbstractMatrix)
    Y = similar(X)
    mul!(Y, op, X)
    return Y
end

Base.size(op::LevelShiftedOperator, args...) = size(op.base_op, args...)
Base.eltype(op::LevelShiftedOperator) = eltype(op.base_op)
LinearAlgebra.ishermitian(op::LevelShiftedOperator) = ishermitian(op.base_op)


"""
    ProjectedShiftedOperator(base_op, V, shift)

A wrapper that projects the operator into the orthogonal complement of the subspace `V`, 
while applying a constant energy `shift` directly to the `V` subspace components.

Mathematically, it acts as:
    (1 - V * V') * base_op * (1 - V * V') * X + shift * V * (V' * X)

This ensures that the eigensolver remains within the orthogonal complement of `V`, 
by penalizing any components inside `V` with the positive energy `shift`.
"""
struct ProjectedShiftedOperator{TOp,TV}
    base_op::TOp       # the base operator
    V::TV              # basis for projector
    shift::Float64     # energy penalty shift for the projected subspace
end

function ProjectedShiftedOperator(base_op, V, shift)
    return ProjectedShiftedOperator{typeof(base_op),typeof(V)}(base_op, V, shift)
end

function LinearAlgebra.mul!(Y, op::ProjectedShiftedOperator, X)
    workX = similar(X)
    copy!(workX, X)
    coeffs_in = op.V' * X
    mul!(workX, op.V, coeffs_in, -1.0, 1.0) # workX -= op.V * coeffs_in
    mul!(Y, op.base_op, workX)              #     Y  = op.base_op * workX
    coeffs_out = op.V' * Y
    mul!(Y, op.V, coeffs_out, -1.0, 1.0)    #     Y -= op.V * coeffs_out
    mul!(Y, op.V, coeffs_in, op.shift, 1.0) #     Y += shift*coeffs_in
    return Y
end

function Base.:*(op::ProjectedShiftedOperator, X::AbstractMatrix)
    Y = similar(X)
    mul!(Y, op, X)
    return Y
end

Base.size(op::ProjectedShiftedOperator, args...) = size(op.base_op, args...)
Base.eltype(op::ProjectedShiftedOperator) = eltype(op.base_op)
LinearAlgebra.ishermitian(op::ProjectedShiftedOperator) = ishermitian(op.base_op)


# Identity operator — used as B = I in the standard eigenproblem for canonical virtuals
struct IdentityOperator{T}
    n::Int
end

IdentityOperator(n::Int, ::Type{T}) where {T} = IdentityOperator{T}(n)

function LinearAlgebra.mul!(Y, ::IdentityOperator, X)
    copy!(Y, X)
    return Y
end

Base.:*(::IdentityOperator, X::AbstractMatrix) = copy(X)
Base.size(op::IdentityOperator, args...) = (op.n, op.n)[args...]
Base.eltype(::IdentityOperator{T}) where {T} = T
LinearAlgebra.ishermitian(::IdentityOperator) = true
