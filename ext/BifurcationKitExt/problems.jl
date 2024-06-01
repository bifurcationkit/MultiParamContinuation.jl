using ForwardDiff

"""
$SIGNATURES
"""
function ManifoldProblem_BK(F, u0, par;
                        check_dim::Bool = true,
                        recordFromSolution = (u,p) -> nothing,
                        project = nothing,
                        get_radius = get_radius_default,
                        get_tangent = nothing,
                        event_function = event_default,
                        finalize_solution = finalize_default,
                        prob_cons = nothing)
    bifprob = BifurcationProblem(F, u0, par, (@lens _))
    m = length(BK.residual(bifprob.VF, u0, par))
    ManifoldProblemBK(bifprob, u0, par;
                        m,
                        check_dim,
                        recordFromSolution,
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
    sol = BK.newton(prob.VF, options)
    if ~BK.converged(sol)
        throw("Newton for first point did not converge!!")
    end
    return sol.u
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
        sol = BK.newton(prob_bls, options)
    end
    if BK.converged(sol)
        return sol.u
    else
        return nothing
    end
end