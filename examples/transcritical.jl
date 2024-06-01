using Revise
cd(@__DIR__)
using Pkg
pkg"activate ."

using GLMakie
Makie.inline!(false)
Makie.inline!(true)

using MultiParamContinuation

const MPC = MultiParamContinuation

function F(u,p) 
    x,y,z = u
    [(z-x-.5*(y-1)*(y-1))*(x+y+z) + y/10]
end

prob = ManifoldProblem(F, [-0.6,-0.6,-0.6], nothing;
            finalize_solution = Cube(1.1))

S = continuation(prob,
            Henderson(np0 = 4, 
                      θmin = 0.1,
                    #   use_curvature = true,
                      ),
            CoveringPar(max_charts = 20000,
                    max_steps = 2000,
                    verbose = 0,
                    newton_options = NonLinearSolveSpec(;maxiters = 6, abstol = 1e-12, reltol = 1e-10),
                    R0 = .1,
                    ϵ = 0.005,
                    delta_angle = 10.15,
                    ))

MPC.plotd(S; 
    # draw_circle = true,
    draw_tangent = true, 
    # plot_center = true,
    # put_ids = true,
    ind_plot = [1,3])

step!(S, 1000);fig = MPC.plotd(S; draw_circle = false)


MPC.plot2d(S; 
    # draw_circle = true, 
    draw_tangent = true, 
    plot_center = true,
    # put_ids = true,
    ind_plot = [1,2])