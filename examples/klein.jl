using Revise
cd(@__DIR__)
using Pkg
pkg"activate ."

using GLMakie
Makie.inline!(false)
Makie.inline!(true)

using MultiParamContinuation

using Test, LinearAlgebra, StaticArrays
const MPC = MultiParamContinuation

function F(u, p) 
    x,y,z = u
    r = x^2+y^2+z^2
    SA[(r+2*y-1)*((r-2*y-1)^2-8*z^2)+16*x*z*(r-2*y-1)]
end

prob = ManifoldProblem(F, SA[1,1,0.], nothing)
alg = Henderson(np0 = 4, use_curvature = true)
params = CoveringPar(max_charts = 3000, 
                    max_steps = 1000,
                    verbose = 0,
                    newton_options = NonLinearSolveSpec(;maxiters = 8),
                    Rmax = .4,
                    R0 = 0.1,
                    ϵ = 0.05,
                    )
S = @time MPC.continuation(prob,
            alg,
            params
            )

f = MPC.plotd(S; 
    # draw_circle = true, 
    draw_tangent = true,
    draw_edges = true,
    # plot_center = true,
    # put_ids = true,
    ind_plot = 1:3)

MPC.plot2d(S; 
    draw_circle = true, 
    draw_tangent = true, 
    plot_center = true,
    put_ids = true,
    ind_plot = 1:3)


step!(S, 1000); MPC.plotd(S; draw_edges = true)

MPC.plotcenters(S)

###################################
begin
    f = Figure(size = (800, 800))
    ax = Axis3(f[1,1], aspect = :data, elevation = pi/4, azimuth = -pi/2)
    MPC.plotd(ax, S; 
        # draw_circle = true, 
        draw_tangent = true, 
        draw_edges = true,
        # plot_center = true,
        # put_ids = true,
        ind_plot = 1:3
    )
    f
end