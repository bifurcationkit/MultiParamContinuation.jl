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
    x,y,r,s = u
    [x*(2*x*x+y*y-r+s)+x*y,y*(x*x+2*y*y-r-s)-x*x]
end

prob = ManifoldProblem(F, [0.2,0.2,-0.,0], nothing;
            finalize_solution = ProductSpace([-5,-5,-1,-1],[5,5,1,1]),
            record_from_solution = (u,p;k...) -> [u[3],u[4],u[1]/2+u[2]/2]
)

S = continuation(prob,
            Henderson(np0 = 8, 
                      θmin = 0.05,
                    #   use_curvature = true,
                    use_tree = true,
                      ),
            CoveringPar(max_charts = 20000,
                    max_steps = 2000,
                    verbose = 0,
                    newton_options = NonLinearSolveSpec(;maxiters = 8, abstol = 1e-12, reltol = 1e-10),
                    R0 = .05,
                    ϵ = 0.02,
                    delta_angle = 10.15,
                    ))

MPC.plotcenters(S)

step!(S, 1000);MPC.plotcenters(S)

MPC.plot2d(S; 
    # draw_circle = true, 
    draw_tangent = false, 
    plot_center = true,
    # put_ids = true,
    ind_plot = [1,2]
    )