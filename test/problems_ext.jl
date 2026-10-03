using MultiParamContinuation
using BifurcationKit

using Test, LinearAlgebra, ForwardDiff
const MPC = MultiParamContinuation
const BK = BifurcationKit
const MPCExt = Base.get_extension(MultiParamContinuation, :BifurcationKitExt)

@testset "BifurcationKit extension - problems" begin
    @test MPCExt !== nothing

    weights = MPC.Weight(MPC.TrivialWeight())

    #━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    @testset "constructors from a raw vector field" begin
        Fs(u, p) = [u[1]^2 + u[2]^2 + u[3]^2 - 1] # sphere

        for prob in (MPC.ManifoldProblem_BK(Fs, [1., 0., 0.], nothing),
                     MPC.ManifoldProblemBKMatrixFree(Fs, [1., 0., 0.], nothing))
            @test prob isa MPC.AbstractManifoldProblemBifurcationKit
            @test size(prob) == (3, 1)
            @test prob.u0 == [1., 0., 0.]
            @test prob.params === nothing
            @test prob.prob_cons isa MPCExt.ConstrainedProblem
            @test BK.residual(prob, prob.u0, nothing) ≈ zeros(1) atol = 1e-10
            @test MPC.jacobian(prob, prob.u0, nothing) ≈ [2 0 0] atol = 1e-8
        end

        # keyword arguments are forwarded to the manifold problem
        prob = MPC.ManifoldProblem_BK(Fs, [1., 0., 0.], nothing;
                        record_from_solution = (u, p; k...) -> u[1],
                        get_radius = (u, p) -> 0.42,
                        project = (u, p) -> [1., 0., 0.])
        @test prob.recordFromSolution([2., 0., 0.], nothing) == 2.
        @test prob.get_radius(prob.u0, nothing) == 0.42
        @test MPC._has_projection(prob)
    end

    #━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    @testset "check_dim" begin
        F3(u, p) = [u[1], u[2], u[3]] # m == n
        @test_throws AssertionError MPC.ManifoldProblem_BK(F3, zeros(3), nothing)
        prob = MPC.ManifoldProblem_BK(F3, zeros(3), nothing; check_dim = false)
        @test size(prob) == (3, 3)
    end

    #━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    @testset "constructors from a BifurcationProblem" begin
        F2(u, p) = [u[1] - p.a, u[2] - p.b, u[3] - p.c]
        par = (a = 1., b = 2., c = 3.)
        u0 = [0.5, 0.6, 0.7]
        prob_bk = BifurcationProblem(F2, u0, par, (@optic _.c))

        Z = vcat(u0, par.a, par.b)
        pb_ref = MPCExt.BifurcationProblem_2P(prob_bk, (@optic _.a), (@optic _.b))
        Jref = ForwardDiff.jacobian(z -> pb_ref(z, par), Z)
        unit5(i) = (v = zeros(5); v[i] = 1.; v)

        for (jac, is_mf) in ((nothing, false),
                             (BK.AutoDiffDense(), false),
                             (BK.AutoDiffMF(), true),
                             (BK.FiniteDifferencesMF(), true),
                             (BK.MatrixFree(), true))
            prob = MPC.ManifoldProblem_BK(prob_bk, u0, (@optic _.a), (@optic _.b); jacobian = jac)
            @test prob.u0 == Z
            @test size(prob) == (5, 3)
            @test prob.VF.VF.F isa MPCExt.BifurcationProblem_2P
            @test prob.prob_cons isa MPCExt.ConstrainedProblem
            J = MPC.jacobian(prob, prob.u0, par)
            Jmat = is_mf ? hcat((BK.apply(J, unit5(i)) for i in 1:5)...) : J
            @test Jmat ≈ Jref atol = 1e-6

            probmf = MPC.ManifoldProblemBKMatrixFree(prob_bk, u0, (@optic _.a), (@optic _.b); jacobian = jac)
            @test probmf isa MPC.ManifoldProblemBKMatrixFree
            @test probmf.u0 == Z
        end
    end

    #━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    @testset "BifurcationProblem_2P residual and jacobian" begin
        F2(u, p) = [u[1] - p.a, u[2] - p.b, u[3] - p.c]
        par = (a = 1., b = 2., c = 3.)
        prob_bk = BifurcationProblem(F2, zeros(3), par, (@optic _.c))
        pb = MPCExt.BifurcationProblem_2P(prob_bk, (@optic _.a), (@optic _.b))
        Z = [0.5, 0.6, 0.7, 1.1, 2.2]
        @test pb(Z, par) ≈ [0.5 - 1.1, 0.6 - 2.2, 0.7 - 3.0]

        J = MPCExt.jacobian(pb, Z, par)
        @test J ≈ [1 0 0 -1 0; 0 1 0 0 -1; 0 0 1 0 0] atol = 1e-6
        @test J ≈ ForwardDiff.jacobian(z -> pb(z, par), Z) atol = 1e-6
    end

    #━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    @testset "d2F of the composite problem" begin
        # the extension assumes n - m = 2 for BifurcationProblem_2P composites
        Fq(u, p) = [u[1]^2 + p.a * u[2] + p.b^2 - 1, u[2] - p.a]
        par = (a = 0., b = 0.)
        prob_bk = BifurcationProblem(Fq, [1., 0.], par, (@optic _.a))
        pb = MPCExt.BifurcationProblem_2P(prob_bk, (@optic _.a), (@optic _.b))
        Z = [1., 0., 0., 0.]
        dx1 = [0.3, 0.4, 0.6, 0.7]
        dx2 = [0.1, 0.2, 0.4, 0.5]
        d2 = MPCExt._d2F_2P(pb, Z, par, dx1, dx2)
        exact = [2 * dx1[1] * dx2[1] + (dx1[2] * dx2[3] + dx1[3] * dx2[2]) + 2 * dx1[4] * dx2[4], 0.]
        @test d2 ≈ exact atol = 1e-5

        # d2F dispatches to _d2F_2P through ManifoldProblem_BK
        prob = MPC.ManifoldProblem_BK(prob_bk, [1., 0.], (@optic _.a), (@optic _.b))
        @test MPC.d2F(prob, prob.u0, par, dx1, dx2) ≈ exact atol = 1e-5

        # for a plain vector field, the BK d2F of the underlying problem is used
        Fs(u, p) = [u[1]^2 + u[2]^2 + u[3]^2 - 1]
        prob_s = MPC.ManifoldProblem_BK(Fs, [1., 0., 0.], nothing)
        @test MPC.d2F(prob_s, [1., 0., 0.], nothing, [0., 1., 0.], [0., 1., 0.])[] ≈ 2 atol = 1e-6
        @test MPC.d2F(prob_s, [1., 0., 0.], nothing, [0., 1., 0.], [0., 0., 1.])[] ≈ 0 atol = 1e-6
    end

    #━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    @testset "ConstrainedProblem" begin
        Fs(u, p) = [u[1]^2 + u[2]^2 + u[3]^2 - 1]
        prob_bk = BifurcationProblem(Fs, [1., 0., 0.], nothing, (@optic _))
        Φ = [0. 0.; 1. 0.; 0. 1.]
        wbar = [1., 0.1, 0.1]
        cpb = MPCExt.ConstrainedProblem(prob_bk, Φ, Φ' * wbar, wbar)
        w = [1., 0.2, 0.1]
        @test BK.residual(cpb, w, nothing) ≈ vcat(Fs(w, nothing), Φ' * (w - wbar))
        J = MPCExt.jacobian(cpb, w, nothing)
        @test J ≈ vcat(BK.jacobian(prob_bk, w, nothing), Φ') atol = 1e-8
        @test BK.getu0(cpb) == prob_bk.u0
        @test BK.getparams(cpb) === nothing
    end

    #━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    @testset "dispatch glue" begin
        Fs(u, p) = [u[1]^2 + u[2]^2 + u[3]^2 - 1]
        prob = MPC.ManifoldProblem_BK(Fs, [1., 0., 0.], nothing)
        @test MPC.jacobian(prob, [1., 0., 0.], nothing) ≈ BK.jacobian(prob.VF, [1., 0., 0.], nothing)
        @test BK.residual(prob, [1., 0., 0.], nothing) == BK.residual(prob.VF, [1., 0., 0.], nothing)
        @test BK.getlens(prob) === nothing
    end

    #━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    @testset "correct_guess" begin
        F(u, p) = [u[1] - 1, u[2] - 2]
        u0 = [1.2, 1.8, 0.]
        prob = MPC.ManifoldProblem_BK(F, u0, nothing)
        cpar = CoveringPar(newton_options = BK.NewtonPar(verbose = false, max_iterations = 20))
        rhs = vcat(zeros(2, 1), Matrix(I, 1, 1))
        cache = MPC.HendersonCache(prob, cpar, Henderson(), 1.0, rhs)
        u = MPCExt.correct_guess(cache, cpar.newton_options)
        @test u ≈ [1., 2., 0.] atol = 1e-10
        @test norm(BK.residual(prob, u, nothing)) < 1e-12
    end

    #━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    @testset "project_on_M (dense)" begin
        Fs(u, p) = [u[1]^2 + u[2]^2 + u[3]^2 - 1]
        u0 = [1., 0., 0.]
        prob = MPC.ManifoldProblem_BK(Fs, u0, nothing)
        Φ = [0. 0.; 1. 0.; 0. 1.]
        chart = MPC.new_chart(u0, Φ, 0.2, MPC.init_polygonal_boundary(5, 0.2))
        guess = u0 .+ Φ * [0.1, 0.1]
        cpar = CoveringPar(newton_options = BK.NewtonPar(verbose = false, max_iterations = 20))
        wbar = copy(guess)
        u = MPC.project_on_M(prob, guess, chart, wbar, cpar, weights)
        @test u !== nothing
        @test norm(Fs(u, nothing)) < 1e-8
        @test norm(Φ' * (u - wbar)) < 1e-8
        # the constrained problem has been updated
        @test prob.prob_cons.Φ == Φ
        @test prob.prob_cons.wbar == wbar
        @test prob.prob_cons.prob.u0 == guess

        # failure to converge returns nothing
        u_fail = MPC.project_on_M(prob, [3., 3., 3.], chart, [3., 3., 3.], cpar, weights)
        @test u_fail === nothing

        # a user provided projection is used when present
        prob2 = MPC.ManifoldProblem_BK(Fs, u0, nothing; project = (u, p) -> [1., 0., 0.])
        @test MPC.project_on_M(prob2, [0.5, 0.5, 0.5], chart, [0.5, 0.5, 0.5], cpar, weights) == [1., 0., 0.]
    end

    #━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    @testset "project_on_M (matrix free)" begin
        Fs(u, p) = [u[1]^2 + u[2]^2 + u[3]^2 - 1]
        u0 = [1., 0., 0.]
        prob = MPC.ManifoldProblemBKMatrixFree(Fs, u0, nothing)
        Φ = [0. 0.; 1. 0.; 0. 1.]
        chart = MPC.new_chart(u0, Φ, 0.2, MPC.init_polygonal_boundary(5, 0.2))
        guess = u0 .+ Φ * [0.1, 0.1]
        cpar = CoveringPar(newton_options = BK.NewtonPar(verbose = false, max_iterations = 20),
                            solver_bls = BK.GMRESIterativeSolvers())
        wbar = copy(guess)
        u = MPC.project_on_M(prob, guess, chart, wbar, cpar, weights)
        @test u !== nothing
        @test norm(Fs(u, nothing)) < 1e-8
        @test norm(Φ' * (u - wbar)) < 1e-8
    end

    #━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    @testset "get_tangent with BLSBorderedTangent" begin
        F(u, p) = [u[1] - p.a, u[2] - p.b, u[3] - p.c]
        par = (a = 0., b = 0., c = 0.)
        u0 = [0., 0., 0.]
        prob_bk = BifurcationProblem(F, u0, par, (@optic _.c))
        bls = BK.MatrixFreeBLS(BK.GMRESIterativeSolvers())
        RHS = vcat(zeros(3, 2), Matrix(I, 2, 2))

        for jac in (BK.AutoDiffMF(), BK.FiniteDifferencesMF(), BK.MatrixFree())
            prob = MPC.ManifoldProblem_BK(prob_bk, u0, (@optic _.a), (@optic _.b);
                        jacobian = jac,
                        get_tangent = MPCExt.BLSBorderedTangent(bls))
            Φ = MPC.get_tangent(prob, prob.u0, par, RHS)
            J = MPC.jacobian(prob, prob.u0, par)
            JΦ = hcat((BK.apply(J, Φ[:, i]) for i in 1:2)...)
            @test size(Φ) == (5, 2)
            @test norm(JΦ) < 1e-5
            @test norm(Φ' * Φ - I(2)) < 1e-8

            # continuity of the basis when a guess Φ0 is provided
            Φ2 = MPC.get_tangent(prob, prob.u0, par, RHS, Φ)
            @test norm(Φ2' * Φ2 - I(2)) < 1e-8
            @test norm(Φ2 - Φ) < 1e-6
        end

        # unimplemented tangent algorithms throw
        Fs(u, p) = [u[1]^2 + u[2]^2 + u[3]^2 - 1]
        prob_b = MPC.ManifoldProblem_BK(Fs, [1., 0., 0.], nothing; get_tangent = MPC.BorderedTangent())
        @test_throws AssertionError MPC.get_tangent(prob_b, [1., 0., 0.], nothing, vcat(zeros(1, 2), Matrix(I, 2, 2)))
        prob_q = MPC.ManifoldProblem_BK(Fs, [1., 0., 0.], nothing; get_tangent = MPC.QRDirectTangent())
        @test_throws AssertionError MPC.get_tangent(prob_q, [1., 0., 0.], nothing, vcat(zeros(1, 2), Matrix(I, 2, 2)))
    end

    #━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    @testset "get_curvature (matrix free)" begin
        bls = BK.MatrixFreeBLS(BK.GMRESIterativeSolvers())

        # unit sphere: the curvature estimate is K = 1
        Fs(u, p) = [u[1]^2 + u[2]^2 + u[3]^2 - 1]
        u0 = [1., 0., 0.]
        Jmf = (x, p) -> dx -> ForwardDiff.derivative(t -> Fs(x .+ t .* dx, p), zero(eltype(x)))
        prob_mf = MPC.ManifoldProblemBKMatrixFree(Fs, u0, nothing;
                    J = Jmf,
                    get_tangent = MPCExt.BLSBorderedTangent(bls))
        prob_d = MPC.ManifoldProblem_BK(Fs, u0, nothing)
        Φ = [0. 0.; 1. 0.; 0. 1.]
        K_mf = MPC.get_curvature(prob_mf, u0, Φ, nothing, weights)
        K_d = MPC.get_curvature(prob_d, u0, Φ, nothing, weights)
        @test K_mf ≈ K_d rtol = 1e-6
        @test K_mf ≈ 1 rtol = 1e-6

        # 2-parameter composite with n - m = 2: cross-check the matrix-free and dense
        # implementations of the curvature estimate on the same tangent basis
        # (the extension assumes n - m = 2 for BifurcationProblem_2P composites)
        F(u, p) = [u[1]^2 + u[2]^2 - p.a, u[1] - p.b]
        par = (a = 1., b = 1.)
        prob_bk = BifurcationProblem(F, [1., 0.], par, (@optic _.a))
        prob_mf2 = MPC.ManifoldProblemBKMatrixFree(prob_bk, [1., 0.], (@optic _.a), (@optic _.b);
                    jacobian = BK.MatrixFree(),
                    get_tangent = MPCExt.BLSBorderedTangent(bls))
        prob_d2 = MPC.ManifoldProblem_BK(prob_bk, [1., 0.], (@optic _.a), (@optic _.b))
        Φ2 = [0. (1 / sqrt(6)); 1. 0.; 0. (2 / sqrt(6)); 0. (1 / sqrt(6))]
        K_mf2 = MPC.get_curvature(prob_mf2, prob_mf2.u0, Φ2, par, weights)
        K_d2 = MPC.get_curvature(prob_d2, prob_d2.u0, Φ2, par, weights)
        @test K_mf2 ≈ K_d2 rtol = 1e-4
    end
end
