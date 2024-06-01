module GLMakieExt
    using GLMakie, MultiParamContinuation

    using DocStringExtensions, StaticArrays

    import MultiParamContinuation: plot2d, plotd, plot2d_improved, plotd_improved, plotcenters, is_on_boundary, get_tangent

    import MultiParamContinuation: plotkd, plotkd!

    import MultiParamContinuation: BVHNode, npoints, is_leaf, ρ

    include("plot.jl")
end
