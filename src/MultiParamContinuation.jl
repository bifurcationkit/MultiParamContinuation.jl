module MultiParamContinuation
    using CircularArrays: CircularVector
    using ForwardDiff
    using LinearAlgebra, Parameters
    using Base.Iterators: reverse
    using ProgressLogging: @progress, @withprogress, @logprogress
    using  NonlinearSolve
    using StaticArrays

    using Reexport
    @reexport using NonlinearSolve: solve, NewtonRaphson

    using DocStringExtensions
    ###########
    include("utils.jl")
    include("space.jl")
    include("parameters.jl")
    include("atlas.jl")
    include("bvhtree.jl")
    include("problem.jl")
    include("plot.jl")
    include("henderson.jl")

    export NonLinearSolveSpec
    export CoveringPar
    export Cube, ProductSpace
    export continuation, step!, ManifoldProblem, ManifoldProblem_BK
    export Henderson, HendersonEllipse
    export Atlas, Chart
    export AABB, BVHNode, add!
end
