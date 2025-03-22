using Revise
cd(@__DIR__)
using Pkg
pkg"activate ."

using GLMakie
Makie.inline!(false)
Makie.inline!(true)

using MultiParamContinuation

const MPC = MultiParamContinuation

function F(u,p) # CUSP
    x,y,r,s=u
    [x*(x*x-y*y-r)-2*x*y*y-s, 2*x*x*y+y*(x*x-2*y*y-r)]
end

prob = ManifoldProblem(F, [-1,0,1,0.], nothing; 
                    record_from_solution = (u,p;k...) -> (u[3], u[4], u[1]+u[3]),
                    finalize_solution = Cube(2.),
                    )

S = @time MPC.continuation(prob,
            Henderson(np0 = 8, 
                        use_tree = true
                        ),
            CoveringPar(max_charts = 3500, 
                    max_steps = 860,
                    verbose = 0,
                    newton_options = NonLinearSolveSpec(;maxiters = 5, abstol = 1e-12, reltol = 1e-10),
                    R0 = .1,
                    # R0 = 1/sqrt(2), # for Torus
                    ϵ = Inf,
                    delta_angle = Inf,
                    )
            )

f = MPC.plotd(S; 
    # draw_circle = true, 
    draw_tangent = false, 
    plot_center = true,
    # put_ids = true,
    ind_plot = [1,3])

step!(S, 1500)