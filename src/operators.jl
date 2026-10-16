@doc raw"""
    ProjectedShiftedOperator(base_op, V, shift, penalty)

A wrapper that projects the operator into the orthogonal complement of the subspace
spanned by the orthonormal columns of ``V``, shifts it there by a constant, and assigns a
penalty energy to the subspace itself. Applied to a vector ``X``, it computes
```math
(1 - V V^\dagger) (A + σ) (1 - V V^\dagger) X + p \, V V^\dagger X
```
with ``A`` = `base_op`, ``σ`` = `shift` and ``p`` = `penalty`. The subspace of ``V`` is an
exact eigenspace with eigenvalue ``p``, even if it is not an invariant subspace of ``A``
(e.g. orbitals of a loosely converged SCF), so the eigenvectors in the complement are
exactly orthogonal to ``V``. A penalty above the wanted part of the spectrum keeps
eigensolvers away from the subspace.
"""
struct ProjectedShiftedOperator{TOp,TV}
    base_op::TOp       # the base operator
    V::TV              # basis of the projected-out subspace
    shift::Float64     # energy shift in the orthogonal complement
    penalty::Float64   # eigenvalue of the projected-out subspace
end

function ProjectedShiftedOperator(base_op, V, shift, penalty)
    return ProjectedShiftedOperator{typeof(base_op),typeof(V)}(base_op, V, shift, penalty)
end

function LinearAlgebra.mul!(Y, op::ProjectedShiftedOperator, X)
    workX = similar(X)
    copy!(workX, X)
    coeffs_in = op.V' * X
    mul!(workX, op.V, coeffs_in, -1.0, 1.0)   # workX -= op.V * coeffs_in
    mul!(Y, op.base_op, workX)                #     Y  = op.base_op * workX
    Y .+= op.shift .* workX
    coeffs_out = op.V' * Y
    mul!(Y, op.V, coeffs_out, -1.0, 1.0)      #     Y -= op.V * coeffs_out
    mul!(Y, op.V, coeffs_in, op.penalty, 1.0) #     Y += penalty*coeffs_in
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
