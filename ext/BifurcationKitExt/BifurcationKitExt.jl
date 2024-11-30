module BifurcationKitExt
    using BifurcationKit, MultiParamContinuation
    using DocStringExtensions, LinearAlgebra
    using ForwardDiff
    using Accessors: @set # to modify fields of NewtonPar
    const BK = BifurcationKit
    
    import MultiParamContinuation: correct_guess,
                                    project_on_M,
                                    ManifoldProblemBK,
                                    ManifoldProblem_BK,
                                    get_radius_default,
                                    event_default,
                                    finalize_default,
                                    _has_projection,
                                    get_tangent,
                                    _has_tangent_computation,
                                    jacobian,
                                    d2F,
                                    continuation,
                                    CoveringPar

    include("problems.jl")
end
