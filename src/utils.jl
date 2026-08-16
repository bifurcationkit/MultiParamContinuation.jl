struct TrivialWeight; end

struct Weight{T}
    weights::T
end
@inline get_weights(w::Weight) = w.weights
@inline apply_T(w::Weight, T, du) = _apply_T(get_weights(w), T, du)
@inline _apply_T(::TrivialWeight, T, du) = T' * du
@inline _apply_T(w::AbstractVector, T, du) = T' * (w .* du)

@inline weighted_norm(w::Weight, u) = _weighted_norm(get_weights(w), u)
@inline _weighted_norm(::TrivialWeight, u) = norm(u)
@inline _weighted_norm(w::AbstractVector, u) = norm(sqrt.(w) .* u)

"""
Squared euclidean distance. This version is non allocating compared to `norm(u1 - u2, 2)^2` albeit perhaps less performant for large dimensions.
"""
@inline dist2(w::Weight, u1, u2) = dist2(get_weights(w), u1, u2)
@inline dist2(::TrivialWeight, u1, u2) = mapreduce(x -> abs2(x[1] - x[2]), +, zip(u1, u2))
@inline dist2(w::AbstractVector, u1, u2) = mapreduce(x -> x[3] * abs2(x[1] - x[2]), +, zip(u1, u2, w))


abstract type AbstractJacobianType end

"""
Struct to specify the use of jacobian-free method based on ForwardDiff.jl
"""
struct JacobianFreeFD <: AbstractJacobianType end

"""
Struct to specify the use of jacobian-free method based on user passed code. The constraint corresponding to the tangent space is appended.
"""
struct MyJacobianFree <: AbstractJacobianType end

"""
Struct to specify the use of jacobian-free method based on ForwardDiff.jl
"""
struct MyJacobian <: AbstractJacobianType end

