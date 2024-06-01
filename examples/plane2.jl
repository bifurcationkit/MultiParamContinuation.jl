using Revise
cd(@__DIR__)
using Pkg
pkg"activate ."
using MultiParamContinuation

using GLMakie
Makie.inline!(false)
Makie.inline!(true)


using Test, LinearAlgebra
const MPC = MultiParamContinuation

F(u,p) = [u[1] + u[2] - u[3]]

prob = ManifoldProblem(F, zeros(3), nothing)

S = continuation(prob,
            Henderson(np0 = 4),
            CoveringPar(max_charts = 20000, 
                    max_steps = 100,
                    R0 = .1,
                    ϵ = Inf,
                    delta_angle = Inf,
                    )
                )

MPC.plotd(S; 
    plot_circle = true, 
    # draw_tangent = false, 
    draw_edges = true,
    plot_center = true,
    put_ids = true,
    ind_plot = [1,2])


