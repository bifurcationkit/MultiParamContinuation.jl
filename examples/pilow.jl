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

F(u,p) = [u[1]^4 + u[2]^4 + u[3]^4 - u[1]^2 - u[2]^2 - u[3]^2]

prob = ManifoldProblem(F, [0,0,1.], nothing)

contpar = CoveringPar(max_charts = 1500, 
                                max_steps = 1800,
                                verbose = 0,
                                newton_options = NonLinearSolveSpec(;maxiters = 6, abstol = 1e-12),
                                R0 = .2,
                                ϵ = 0.015,
                                # delta_angle = 10.15,
                                ); 
alg = Henderson(np0 = 6,
                θmin = 0.001,
                use_curvature = true
                )

S = @time MPC.continuation(prob,
            alg,
            contpar)

f = MPC.plotd(S; 
    # draw_circle = true, 
    draw_tangent = true, 
    draw_edges = true,
    # plot_center = true,
    # put_ids = true,
    ind_plot = [1,3])

MPC.plot2d(S; 
    # draw_circle = true, 
    plot_center = true,
    put_ids = true,
    ind_plot = [2,1]
    )


step!(S,500)

###################################
# figure 
f = Figure(size = (800, 800))
ax = Axis3(f[1,1], aspect = :data, elevation = pi/4, azimuth = -pi/2)
MPC.plotd(ax, S; 
    # draw_circle = true, 
    draw_tangent = true, 
    # plot_center = true,
    # put_ids = true,
    ind_plot = [1,3])
f

ax = current_axis()
ax.azimuth = -pi/3