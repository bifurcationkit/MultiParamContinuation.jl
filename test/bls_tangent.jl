using MultiParamContinuation
using BifurcationKit

using Test, LinearAlgebra
const MPC = MultiParamContinuation
const BK = BifurcationKit
const MPCExt = Base.get_extension(MultiParamContinuation, :BifurcationKitExt)

@testset "Bordered tangent with BifurcationKit BLS" begin

    @test MPCExt !== nothing

    #━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    # Matrix-free jacobian (jvp) on a 2-parameter composite problem.
    # F(u, p) = [u1 - a, u2 - b, u3 - c]  => manifold {u1 = a, u2 = b, u3 = c}, dim 2
    F(u, p) = [u[1] - p.a, u[2] - p.b, u[3] - p.c]
    prob_bk = BifurcationProblem(F, [0., 0., 0.], (a = 0., b = 0., c = 0.), (@optic _.a))

    n, k = 5, 2
    RHS = vcat(zeros(3, 2), Matrix(I, 2, 2)) # n × k

    prob = MPC.ManifoldProblem_BK(prob_bk, [0., 0., 0.], (@optic _.a), (@optic _.b);
                jacobian = BK.MatrixFree(),
                get_tangent = MPCExt.BLSBorderedTangent(BK.MatrixFreeBLS(BK.GMRESIterativeSolvers())))

    u0, par = prob.u0, prob.params
    Φ = MPC.get_tangent(prob, u0, par, RHS)
    Jfun = MPC.jacobian(prob, u0, par) # jvp
    JΦ = hcat((Jfun(Φ[:, i]) for i in 1:k)...)

    @test size(Φ) == (n, k)
    @test norm(JΦ) < 1e-6
    @test norm(Φ' * Φ - I(k)) < 1e-8

    #━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    # Dense jacobian, function-form problem, matrix based BLS
    Fs(u, p) = [u[1]^2 + u[2]^2 + u[3]^2 - 1] # sphere
    prob_s = MPC.ManifoldProblem_BK(Fs, [1., 0., 0.], nothing;
                get_tangent = MPCExt.BLSBorderedTangent(BK.MatrixBLS()))

    u0 = [1., 0., 0.]
    Φs = MPC.get_tangent(prob_s, u0, nothing, vcat(zeros(1, 2), Matrix(I, 2, 2)))
    J = MPC.jacobian(prob_s, u0, nothing)

    @test size(Φs) == (3, 2)
    @test norm(J * Φs) < 1e-8
    @test norm(Φs' * Φs - I(2)) < 1e-8
end
