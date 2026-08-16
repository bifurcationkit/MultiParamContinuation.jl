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

F(u,p) = [u[1]^2 + u[2]^2 + p.R * u[3]^2 - 1]

prob = ManifoldProblem(F, 
                    [0., 0, 1],
                    (R = 5.,);
                        )

S = MPC.continuation(prob,
            Henderson(np0 = 4,
                        θmin = 0.001,
                        use_curvature = true,
                        ),
            CoveringPar(max_charts = 20000, 
                    max_steps = 2000,
                    verbose = 0,
                    Rmax = 0.2,
                    R0 = 0.1,
                    ϵ = 0.005,
                    α = 1.3
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