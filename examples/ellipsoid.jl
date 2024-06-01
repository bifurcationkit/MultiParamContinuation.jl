using Revise
cd(@__DIR__)
using Pkg
pkg"activate ."

using GLMakie
Makie.inline!(false)
Makie.inline!(true)

using MultiParamContinuation

using Test, LinearAlgebra
const MPC = MultiParamContinuation

F(u,p) = [u[1]^2 + u[2]^2 + p.R * u[3]^2 - 1] # sphere

prob = ManifoldProblem(F, 
                    [0. ,0, 1],
                    (R = 2.,);
                    # get_tangent
                        )

# problem quand j'ajoute une charte
S = MPC.continuation(prob,
            Henderson(np0 = 4,
                        θmin = 0.001,
                        # use_curvature = true,
                        ),
            CoveringPar(max_charts = 20000, 
                    max_steps = 500,
                    verbose = 0,
                    R0 = .2,
                    ϵ = 0.0125,
                    # delta_angle = 0.15,
                    )
            )

MPC.plotd(S; 
    # draw_circle = true, 
    draw_tangent = true, 
    draw_edges = true,
    # plot_center = true,
    # put_ids = true,
    ind_plot = 1:3
    )