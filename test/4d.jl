using MultiParamContinuation

using LinearAlgebra
const MPC = MultiParamContinuation

function F(u,p) 
    x,y,r,s = u
    [x*(2*x*x+y*y-r+s)+x*y,y*(x*x+2*y*y-r-s)-x*x]
end

prob = ManifoldProblem(F, [0.2,0.2,-0.,0], nothing;
            finalize_solution = ProductSpace([-5,-5,-1,-1],[5,5,1,1]),
            recordFromSolution = (u,p) -> [u[3],u[4],u[1]/2+u[2]/2]
)

S = continuation(prob,
            Henderson(np0 = 8, 
                      θmin = 0.05,
                      use_tree = true,
                      ),
            CoveringPar(max_charts = 2000,
                    max_steps = 4000,
                    verbose = 0,
                    newton_options = NonLinearSolveSpec(;maxiters = 8, abstol = 1e-12, reltol = 1e-10),
                    R0 = .05,
                    ϵ = 0.02,
                    delta_angle = 10.15,
                    ))

@test length(S) > 1000