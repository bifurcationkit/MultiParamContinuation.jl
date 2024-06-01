using MultiParamContinuation
using Test

@testset "aabb" begin
    include("aabb.jl")
end

@testset "plane" begin
    include("plane.jl")
end

@testset "sphere" begin
    include("sphere.jl")
end

@testset "4d" begin
    include("4d.jl")
end

@testset "sphere BifurcationKit" begin
    include("sphere_bk.jl")
end
