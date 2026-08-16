using Revise
cd(@__DIR__)
using Pkg
pkg"activate ."

using GLMakie
Makie.inline!(false)
Makie.inline!(true)

using MultiParamContinuation

const MPC = MultiParamContinuation

F(u,p) = [u[2]^2 + u[3]^2 - 1]

prob = ManifoldProblem(F, [0,0,1.], nothing;
            finalize_solution = (u,p) -> -2<=u[1]<=2)

contpar = CoveringPar(max_charts = 1500, 
                                max_steps = 320,
                                verbose = 0,
                                newton_options = NonLinearSolveSpec(;maxiters = 5, abstol = 1e-12, reltol = 1e-10),
                                R0 = .3,
                                ); 
alg = Henderson(np0 = 4,
                θmin = 0.001,
                use_curvature = true
                )

S = continuation(prob,
            alg,
            contpar)

MPC.plotd(S; 
    # draw_circle = true, 
    draw_tangent = true,
    draw_edges = true,
    plot_center = true,
    # put_ids = true,
    ind_plot = 1:3
    )


MPC.plot2d(S; 
    draw_circle = true, 
    draw_tangent = true,
    draw_edges = true,
    plot_center = true,
    put_ids = true,
    # ind_plot = [1,3]
    )

step!(S,1000);fig = MPC.plotd(S, draw_edges = true)