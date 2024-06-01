using Revise
cd(@__DIR__)
using Pkg
pkg"activate ."

using GLMakie
Makie.inline!(false)
Makie.inline!(true)

using BifurcationKit

using MultiParamContinuation

using Test, LinearAlgebra
const MPC = MultiParamContinuation

F(u,p) = [u[1]^2 + u[2]^2 + u[3]^2 - 1] # sphere

prob = ManifoldProblem_BK(F, 
                    [1.1,0.,0.],
                    nothing;
                    # get_tangent
                        )

# problem quand j'ajoute une charte
S = MPC.continuation(prob,
            Henderson(np0 = 5,
                        # use_curvature = true,
                        ),
            CoveringPar(max_charts = 20000, 
                    max_steps = 250,
                    verbose = 0,
                    newton_options = NewtonPar(max_iterations = 5),
                    R0 = .2,
                    ϵ = 10.15,
                    delta_angle = 10.1,
                    )
            )

MPC.plotd(S; 
    # draw_circle = true, 
    draw_tangent = true, 
    draw_edges = true,
    # plot_center = true,
    # put_ids = true,
    ind_plot = [1,3])