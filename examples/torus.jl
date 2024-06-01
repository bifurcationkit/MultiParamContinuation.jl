using Revise
cd(@__DIR__)
using Pkg
pkg"activate ."

using GLMakie
Makie.inline!(false)
Makie.inline!(true)

# using WGLMakie

using MultiParamContinuation

using Test, LinearAlgebra
const MPC = MultiParamContinuation

function F(u,p)
    x,y,z = u
    R2 = 0.5
    R1 = 0.8
    return [(sqrt(x^2 + y^2) - R1)^2 + z^2 - R2^2]
end

prob = ManifoldProblem(F, [-0.05,1.29,-0.05], nothing)

S = @time continuation(prob,
            Henderson(np0 = 8, use_tree = true),
            CoveringPar(max_charts = 3200, 
                    max_steps = 1000,
                    verbose = 0,
                    newton_options = NonLinearSolveSpec(;maxiters = 5, abstol = 1e-12, reltol = 1e-10),
                    R0 = .1,
                    ϵ = Inf,
                    delta_angle = Inf,
                    ))

MPC.plotd(S; 
    # draw_circle = true, 
    draw_tangent = true, 
    draw_edges = true,
    # plot_center = true,
    # put_ids = true,
    ind_plot = [1,3])

MPC.plot2d(S; 
    # draw_circle = true, 
    # plot_center = true,
    # put_ids = true,
    # ind_plot = [1,3]
    )


step!(S,1000);fig = MPC.plotd(S; draw_circle = false, ind_plot = [1,3])


###################################
# figure for the README
f = Figure(size = (800, 800))
ax = Axis3(f[1,1], aspect = :data, elevation = pi/4, azimuth = -pi/3)
MPC.plotd(ax, S; 
    # draw_circle = true, 
    draw_tangent = true, 
    # plot_center = true,
    # put_ids = true,
    ind_plot = [1,3])
f