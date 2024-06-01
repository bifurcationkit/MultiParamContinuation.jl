"""
Squared euclidean distance. This version is non allocating compared to `norm(u1 - u2, 2)^2` albeit perhaps less performant for large dimensions.
"""
@inline dist2(u1, u2) = mapreduce(x -> abs2(x[1] - x[2]), +, zip(u1, u2))

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

