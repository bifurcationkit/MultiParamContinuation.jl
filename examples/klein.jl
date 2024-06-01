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
    x,y,z = u
    r = x^2+y^2+z^2
    [(r+2*y-1)*((r-2*y-1)^2-8*z^2)+16*x*z*(r-2*y-1)]
end

prob = ManifoldProblem(F, [1,1,0.], nothing)

S = @time continuation(prob,
            Henderson(np0 = 4),
            CoveringPar(max_charts = 1500, 
                    max_steps = 2000,
                    verbose = 0,
                    newton_options = NonLinearSolveSpec(;maxiters = 8),
                    R0 = .15,
                    ϵ = 0.1,
                    delta_angle = Inf,
                    ))

f = MPC.plotd(S; 
    # draw_circle = true, 
    draw_tangent = true,
    draw_edges = true,
    # plot_center = true,
    # put_ids = true,
    ind_plot = [1,3])

MPC.plot2d(S; 
draw_circle = true, 
    draw_tangent = true, 
    plot_center = true,
    put_ids = true,
    ind_plot = [2,3])


step!(S, 1000)

MPC.plotcenters(S)

###################################
# figure for the README
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