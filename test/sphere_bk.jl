using MultiParamContinuation
using BifurcationKit

using Test, LinearAlgebra
const MPC = MultiParamContinuation

F(u,p) = [u[1]^2 + u[2]^2 + u[3]^2 - 1] # sphere

prob = ManifoldProblem_BK(F, 
                    [1,0.,0.],
                    nothing;
                        )

# problem quand j'ajoute une charte
S = MPC.continuation(prob,
            Henderson(np0 = 5,
                        use_curvature = true,
                        use_tree = true
                        ),
            CoveringPar(max_charts = 20000, 
                    max_steps = 250,
                    verbose = 0,
                    newton_options = NewtonPar(max_iterations = 5),
                    R0 = .2,
                    ϵ = 10.15,
                    delta_angle = 10.1,
                    )
            )

Sl = length(S)

# nothing is added
MPC.step!(S, 1000)
@test length(S) == Sl