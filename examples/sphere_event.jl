using Revise
cd(@__DIR__)
using Pkg
pkg"activate ."

using GLMakie
Makie.inline!(false)
Makie.inline!(true)

using MultiParamContinuation, LinearAlgebra

const MPC = MultiParamContinuation

F(u,p) = [u[1]^2 + u[2]^2 + u[3]^2 - 1]

function get_tangent(u, par)
    # compute the normal
    n = u
    n ./= norm(u)
    u0x, u0y, u0z = n
    
    if u0x !=1
        v0x=1 -u0x*u0x
        v0y=  -u0x*u0y
        v0z=  -u0x*u0z
    elseif u0y != 1
        v0x=  -u0y*u0x
        v0y=1 -u0y*u0y
        v0z=  -u0y*u0z
        
    else
        v0x=  -u0z*u0x
        v0y=  -u0z*u0y
        v0z=1 -u0z*u0z
    end
    s = 1/sqrt(v0x*v0x+v0y*v0y+v0z*v0z)
    v0x=s*v0x
    v0y=s*v0y
    v0z=s*v0z
    
    v1x=-v0y*u0z+v0z*u0y
    v1y=-v0z*u0x+v0x*u0z
    v1z=-v0x*u0y+v0y*u0x
    T = [v0x v0y v0z;
         v1x v1y v1z]'
end

event_function(u, p) = u[3]

prob = ManifoldProblem(F, 
                    [0,0.,1.],
                    nothing;
                    get_tangent,
                    event_function
                        )

# problem quand j'ajoute une charte
S = MPC.continuation(prob,
            Henderson(np0 = 4),
            CoveringPar(max_charts = 2000, 
                    max_steps = 1000,
                    verbose = 0,
                    newton_options = NonLinearSolveSpec(;maxiters = 5, abstol = 1e-12, reltol = 1e-10),
                    R0 = .1,
                    )
            )

MPC.plotd(S; 
    draw_edges = true,
    draw_tangent = true,
    plot_center = false,
    )

for c in S.atlas
    @info c.label
end