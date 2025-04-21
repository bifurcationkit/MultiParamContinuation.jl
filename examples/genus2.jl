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
    [2y*(y^2-3x^2)*(1-z^2)+(x^2+y^2)^2-(9z^2-1)*(1-z^2)]
end

event_function(u, p) = u[3]

prob = ManifoldProblem(F, [0.,0,-0.1], nothing;
            event_function,
            finalize_solution = Cube(20.1))

S = MPC.continuation(prob,
            Henderson(np0 = 6, 
                      θmin = 0.01,
                      θmax = 1.5,
                      use_curvature = true,
                      ),
            CoveringPar(max_charts = 7000,
                    max_steps = 1300,
                    verbose = 0,
                    newton_options = NonLinearSolveSpec(;maxiters = 5, abstol = 1e-12),
                    R0 = .25,
                    ϵ = 0.025,
                    delta_angle = 0.1,
                    ))

MPC.plotd(S; 
    # draw_tangent = false,
    # plot_center = true,
    draw_edges = true,
    )

step!(S, 2000);fig = MPC.plotd(S; )

MPC.plot2d(S;
    # draw_circle = true, 
    # draw_tangent = true, 
    # plot_center = false,
    put_ids = true,
    ind_plot = [1,2]
    )

