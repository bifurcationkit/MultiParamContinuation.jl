struct TrivialWeight; end

struct Weight{T}
    weights::T
end
@inline get_weights(w::Weight) = w.weights
@inline apply_T(w::Weight, T, du) = _apply_T(get_weights(w), T, du)
@inline _apply_T(::TrivialWeight, T, du) = T' * du
@inline _apply_T(w::AbstractVector, T, du) = T' * (w .* du)
#━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
@inline weighted_norm(w::Weight, u) = _weighted_norm(get_weights(w), u)
@inline _weighted_norm(::TrivialWeight, u) = norm(u)
@inline _weighted_norm(w::AbstractVector, u) = norm(sqrt.(w) .* u)
#━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
"""
Squared euclidean distance. This version is non allocating compared to `norm(u1 - u2, 2)^2` albeit perhaps less performant for large dimensions.
"""
@inline dist2(w::Weight, u1, u2) = dist2(get_weights(w), u1, u2)
@inline dist2(::TrivialWeight, u1, u2) = mapreduce(x -> abs2(x[1] - x[2]), +, zip(u1, u2))
@inline dist2(w::AbstractVector, u1, u2) = mapreduce(x -> x[3] * abs2(x[1] - x[2]), +, zip(u1, u2, w))
#━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
"""
Measure of the difference between the tangent spaces spanned by the (n × k) matrices
`A` and `B`. It returns the minimum over the columns `bᵢ` of `B` of the squared weighted
norm of the projection of `bᵢ` onto `span(A)`,

    dot = minᵢ ‖Aᵀ D bᵢ‖²

with `D = Diagonal(w)`. For `k = 1` it is `cos²θ`, `θ` being the angle between the two
tangent lines. Charts whose tangent spaces have `dot < dotmin` are considered not to
intersect (see `CoveringPar.dotmin`).
"""
@inline matrix_dot(w::Weight, A, B) = _matrix_dot(get_weights(w), A, B)
@inline _matrix_dot(::TrivialWeight, A, B) = minimum(sum(abs2, A' * B; dims = 1))
@inline _matrix_dot(w::AbstractVector, A, B) = minimum(sum(abs2, A' * Diagonal(w) * B; dims = 1))
#━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
"""
Orthonormalize the columns of `T` (n × k) in the weighted metric D = Diagonal(w),
i.e. return Φ = T * R⁻¹ with Φ' * D * Φ = I.
"""
function _weighted_orthonormalize(T, w::AbstractVector)
    G = T' * (Diagonal(w) * T)
    R = cholesky(Hermitian(G)).U
    return T / R
end

weighted_orthonormalize(T, w::Weight) = _weighted_orthonormalize(T, get_weights(w))
_weighted_orthonormalize(T, ::TrivialWeight) = Matrix(qr(T).Q)
_weighted_orthonormalize(T::StA.StaticArray, ::TrivialWeight) = SMatrix(qr(T).Q)
#━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
_myrand(::AbstractArray{𝒯}, n, m) where {𝒯} = rand(𝒯, n, m)
_myrand(::StA.StaticVector{N, 𝒯}, n, m) where {N, 𝒯} = rand(StA.SMatrix{n, m, 𝒯, n*m})

@inline __myrhs(𝒯, m, dim) = vcat(zeros(𝒯, m, dim), I(dim))
_myrhs(::AbstractArray{𝒯}, m, dim) where {𝒯} = __myrhs(𝒯, m, dim)
_myrhs(::StA.StaticVector{N, 𝒯}, m, dim) where {N, 𝒯} = StA.SMatrix{m+dim, dim}(__myrhs(𝒯, m, dim))
#━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
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

