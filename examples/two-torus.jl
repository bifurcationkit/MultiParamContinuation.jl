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
    [((x^2+y^2)^2−x^2+y^2)^2+z^2-0.1^2]
end

prob = ManifoldProblem(F, [0.,0,-0.1], nothing;
            finalize_solution = Cube(20.1))

S = continuation(prob,
            Henderson(np0 = 4, 
                      θmin = 0.05,
                      use_curvature = true,
                      ),
            CoveringPar(max_charts = 20000,
                    max_steps = 2000,
                    verbose = 0,
                    newton_options = NonLinearSolveSpec(;maxiters = 5, abstol = 1e-12, reltol = 1e-10),
                    R0 = .05,
                    ϵ = 0.01,
                    delta_angle = 10.01,
                    ))

MPC.plotd(S; 
    # draw_circle = true,
    draw_edges = true,
    # plot_center = true,
    )

step!(S, 3000);fig = MPC.plotd(S; draw_circle = false)


MPC.plotcenters(S)

MPC.plot2d_improved(S; 
    # draw_circle = true, 
    draw_tangent = true, 
    plot_center = false,
    # put_ids = true,
    ind_plot = [1,2])