using LinearAlgebra

# Level shifted Operator (Shift + Penalty)
struct LevelShiftedOperator{TH, TV}
    base_op::TH                      # the operator
    V::TV                            # basis for penalty projector
    ε_homo::Float64                  # energy of HOMO           
    safe_shift::Float64              # avoid zero gap (~ 1e-5 Ha)
    penalty::Float64                 # should be 2*Ecut
end

function LevelShiftedOperator(base_op, V, ε_homo, safe_shift, penalty)
    return LevelShiftedOperator{typeof(base_op), typeof(V)}(
        base_op, V, ε_homo, safe_shift, penalty
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
 

# Projected and shifted operator (using (1-VV') as projector)
struct ProjectedShiftedOperator{TOp, TV}
    base_op::TOp       # the operator
    V::TV              # basis for projector
    shift::Float64     # should be a lilttle larger than the lower epsilon_hf
end

function ProjectedShiftedOperator(base_op, V, shift)
    return ProjectedShiftedOperator{typeof(base_op), typeof(V)}(
        base_op, V, shift
    )
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

IdentityOperator(n::Int, ::Type{T}) where T = IdentityOperator{T}(n)

function LinearAlgebra.mul!(Y, ::IdentityOperator, X)
    copy!(Y, X)
    return Y
end

Base.:*(::IdentityOperator, X::AbstractMatrix) = copy(X)
Base.size(op::IdentityOperator, args...) = (op.n, op.n)[args...]
Base.eltype(::IdentityOperator{T}) where T = T
LinearAlgebra.ishermitian(::IdentityOperator) = true
