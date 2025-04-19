using ForwardDiff

"""
$SIGNATURES

Make a manifold problem from a `BifurcationProblem`.
"""
function ManifoldProblem_BK(F, u0, par;
                        check_dim::Bool = true,
                        record_from_solution = (u,p;k...) -> nothing,
                        project = nothing,
                        get_radius = get_radius_default,
                        get_tangent = nothing,
                        event_function = event_default,
                        finalize_solution = finalize_default,
                        prob_cons = nothing)
    bifprob = BifurcationProblem(F, u0, par, (@optic _))
    m = length(BK.residual(bifprob.VF, u0, par))
    ManifoldProblemBK(bifprob, u0, par;
                        m,
                        check_dim,
                        record_from_solution,
                        project,
                        get_radius,
                        get_tangent,
                        event_function,
                        finalize_solution,
                        prob_cons)
end

jacobian(pb::ManifoldProblemBK, u, p) = BK.jacobian(pb.VF, u, p)
BK.residual(pb::ManifoldProblemBK, u, p) = BK.residual(pb.VF, u, p)
d2F(pb::ManifoldProblemBK, x, p, dx1, dx2) = BK.d2F(pb.VF, x, p, dx1, dx2)
BK.getlens(::ManifoldProblemBK) = nothing

function correct_guess(cache, options::BK.NewtonPar)
    prob = cache.prob
    sol = BK.solve(prob.VF, BK.Newton(), options)
    if ~BK.converged(sol)
        throw("Newton for first point did not converge!!")
    end
    return sol.u
end

##############################################################################################################
"""
$SIGNATURES

Create a BifurcationProblem with two parameter axes for which we continue the zeros.
"""
struct BifurcationProblem_2P{T1, T2, T3}
    prob::T1
    lens1::T2
    lens2::T3
end

function (pb::BifurcationProblem_2P)(Z, par)
    u = @view Z[1:end-2]
    p1 = Z[end-1]
    p2 = Z[end]
    par2 = BK._set(par, (pb.lens1, pb.lens2), (p1, p2))
    BK.residual(pb.prob, u, par2)
end

function jacobian(pb::BifurcationProblem_2P, Z, par)
    u = @view Z[1:end-2]
    p1 = Z[end-1]
    p2 = Z[end]
    par2 = BK._set(par, (pb.lens1, pb.lens2), (p1, p2))
    # @error (@which BK._get_matrix(BK.jacobian(pb.prob, u, par2)))
    J0 = BK._get_matrix(BK.jacobian(pb.prob, u, par2))
    l1 = ForwardDiff.derivative(z -> pb((@set Z[end-1] = z), par), Z[end-1])
    l2 = ForwardDiff.derivative(z -> pb((@set Z[end] = z),   par), Z[end])
    hcat(J0, l1, l2)
end

##############################################################################################################
struct _A{T1, T2, T3, T4}
    prob::T1
    Φ::T2
    Φwbar::T3
    wbar::T4
end

function (pb::_A)(w, p)
    vcat(BK.residual(pb.prob.VF, w, p), pb.Φ' * (w - pb.wbar))
end

function jacobian(pb::_A, w, p)
    J0 = BK.jacobian(pb.prob, w, p)
    vcat(J0, pb.Φ')
end

"""
$SIGNATURES

Make a manifold problem from a `BifurcationProblem` and specifying two parameter axes.
"""
function ManifoldProblem_BK(prob_bk::BK.AbstractBifurcationProblem,
                            u0, 
                            lens1, 
                            lens2;
                            check_dim::Bool = true,
                            record_from_solution = (u,p;k...) -> nothing,
                            project = nothing,
                            get_radius = get_radius_default,
                            get_tangent = nothing,
                            event_function = event_default,
                            finalize_solution = finalize_default)
    par = prob_bk.params
    m = length(BK.residual(prob_bk, prob_bk.u0, par))

    # make a bifurcation problem with two parameters axes
    pb_composite = BifurcationProblem_2P(prob_bk, lens1, lens2)
    new_u0 = vcat(u0, BK._get(par, lens1), BK._get(par, lens2))

    prob_mpc = BifurcationProblem(pb_composite, 
                    new_u0, 
                    par, 
                    (@optic _); 
                    J = (x, p) -> jacobian(pb_composite, x, p)
                    )
    𝒯 = eltype(new_u0)
    Φ = zeros(𝒯, m+2, 2)
    wbar = zeros(𝒯, m+2)
    prob_cons = _A(prob_mpc, Φ, Φ' * wbar, wbar)

    ManifoldProblemBK(
                        prob_mpc,
                        new_u0, 
                        par;
                        m,
                        check_dim,
                        record_from_solution,
                        project,
                        get_radius,
                        get_tangent,
                        event_function,
                        finalize_solution,
                        prob_cons,
                    )
end

##############################################################################################################
function project_on_M(prob, guess, chart::Chart, wbar, cpar::CoveringPar{T, <: BK.NewtonPar}) where {T}
    if _has_projection(prob)
        return project(prob, guess, prob.params)
    else
        options = cpar.newton_options
        if ~isnothing(cpar.solver_bls)
            options = @set options.linsolver = cpar.solver_bls
        end
        Φ = chart.Φ
        function f(w, p)
            vcat(BK.residual(prob.VF, w, p), Φ' * (w - wbar))
        end
        prob_bls = BifurcationProblem(f, guess, BK.getparams(prob.VF))
        sol = BK.solve(prob_bls, Newton(), options)
    end
    if BK.converged(sol)
        return sol.u
    else
        return nothing
    end
end