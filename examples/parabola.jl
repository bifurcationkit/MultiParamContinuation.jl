using Revise
cd(@__DIR__)
using Pkg
pkg"activate ."

using GLMakie
Makie.inline!(false)
Makie.inline!(true)

using MultiParamContinuation

const MPC = MultiParamContinuation

F(u,p) = [u[1]^12 + u[2]^12 - u[3]]

prob = ManifoldProblem(F, zeros(3), nothing;
            finalize_solution = (u,p) -> (u[3] < 2) * (u[1]>-0.1) * (u[2]>-0.1))

S = continuation(prob,
            Henderson(
                      use_curvature = true,
                      ),
            CoveringPar(max_charts = 20000,
                    max_steps = 1000,
                    # verbose = 2,
                    newton_options = NonLinearSolveSpec(;maxiters = 5, abstol = 1e-12, reltol = 1e-10),
                    Rmax = .2,
                    R0 = 0.01,
                    ϵ = 0.005,
                    ))

MPC.plotd(S; 
    # draw_circle = true,
    draw_tangent = true, 
    draw_edges = true,
    # plot_center = true,
    # put_ids = true,
    # ind_plot = [1,3]
    )

step!(S,1000);fig = MPC.plotd(S; circle = false, draw_edges = true)

MPC.plot2d(S; 
    # draw_circle = true, 
    draw_tangent = true, 
    plot_center = true,
    # put_ids = true,
    ind_plot = [1,2])