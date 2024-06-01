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

function F(u,p) 
    z,y,x = u
    [z*x - x^3 + y/10]
end

prob = ManifoldProblem(F, [-0.,-0.,-1], nothing;
            finalize_solution = ProductSpace([-2,-0.8,-1],[2.,1,2]))

S = continuation(prob,
            Henderson(np0 = 4, 
                      θmin = 0.01,
                      use_curvature = true,
                      ),
            CoveringPar(max_charts = 20000,
                    max_steps = 2000,
                    verbose = 0,
                    newton_options = NonLinearSolveSpec(;maxiters = 6, abstol = 1e-12, reltol = 1e-10),
                    R0 = .1,
                    ϵ = 0.15,
                    delta_angle = 0.1,
                    ))

MPC.plotd(S; 
    # draw_circle = true,
    draw_tangent = true, 
    draw_edges = true,
    ind_plot = [1,3])

step!(S,500);fig = MPC.plotd(S)

MPC.plotcenters(S)

MPC.plot2d(S; 
    # draw_circle = true, 
    draw_tangent = true, 
    plot_center = true,
    # put_ids = true,
    ind_plot = [1,2])