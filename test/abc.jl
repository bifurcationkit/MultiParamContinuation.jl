using MultiParamContinuation
using BifurcationKit

using Test, LinearAlgebra
const MPC = MultiParamContinuation

function abc!(dz, z, p, t = 0)
    (;D, B, σ, β, α) = p
    u1, u2, u3 = z
    dz[1] = -u1 + D*(1 - u1)*exp(u3)
    dz[2] = -u2 + D*(1 - u1)*exp(u3) - D*σ*u2*exp(u3)
    dz[3] = -u3 - β*u3 + D*B*(1 - u1)*exp(u3) + D*B*α*σ*u2*exp(u3)
    dz
end

# we group the differentials together
prob_bk = BifurcationProblem(abc!, [1., 0., 0. ], (D = 0.11, B = 8., α = 1., σ = 0.04, β = 1.56), (@optic _.D), 
        record_from_solution = (x, p; k...) -> (u3 = x[3], u1 = x[1], u2 = x[2]),)

opts_br = ContinuationPar(p_max = 1.5, n_inversion = 8, nev = 3)
br = BifurcationKit.continuation(prob_bk, PALC(), opts_br; normC = norminf)

prob = MPC.ManifoldProblem_BK(
                        prob_bk, br.sol[1].x, (@optic _.D), (@optic _.β);
                        )

S_eq = @time MPC.continuation(prob,
                        Henderson(np0 = 3,
                                    θmin = 0.001,
                                    use_tree = true,
                                  ),
                        CoveringPar(max_charts = 20000,
                                max_steps = 1000,
                                verbose = 0,
                                newton_options = NewtonPar(tol = 1e-10, verbose = false),
                                R0 = .04,
                                ϵ = 0.1,
                                delta_angle = 10.15,
                                ))

# Hopf continuation
opts_cover = CoveringPar(max_charts = 1000,
    max_steps = 5,
    verbose = 1,
    newton_options = NewtonPar(tol = 1e-11, verbose = false),
    R0 = .31,
    ϵ = 0.15,
    # delta_angle = 10.15,
    )

atlas_hopf = @time MPC.continuation(deepcopy(br), 1, 
    (@optic _.α), (@optic _.β), 
    opts_cover;
    alg = Henderson(np0 = 5,
                θmin = 0.001,
                use_curvature = false,
                use_tree = true,
            ),
    )