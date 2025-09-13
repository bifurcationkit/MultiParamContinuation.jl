"""
$SIGNATURES

Compute the manifold of bifurcations in the parameters `(getlens(br.prob), lens2, lens3)` from the ind-th bifurcation point in `br`.
"""
function continuation(br::BK.ContResult,
                        ind::Int,
                        lens2::BK.AllOpticTypes,
                        lens3::BK.AllOpticTypes,
                        contparams::CoveringPar;
                        alg = Henderson(np0 = 3,
                                        θmin = 0.001,
                                        use_curvature = false,
                                        use_tree = false,
                                    ),
                        record_from_solution = BK.record_sol_default,
                        continuationpar_bk::ContinuationPar = br.contparams,
                        finalize_solution = (X, p; k...) -> true,
                        options_bk = (
                                jacobian_ma = BK.MinAugMatrixBased(),
                                update_minaug_every_step = 1,
                                verbosity = 3,
                            ),
                        )
    if lens3 == lens2 || BK.getlens(BK.getprob(br)) in (lens2, lens3)
        error("You must pass different lenses")
    end

    @assert options_bk.jacobian_ma == BK.MinAugMatrixBased()

    # we compute one step with continuation to get the problem and initial point
    @reset continuationpar_bk.max_steps = 1
    @reset continuationpar_bk.detect_bifurcation = 2
    br_codim2 = BK.continuation(br, ind, lens2, continuationpar_bk; options_bk...)
    # we now extract the bifurcation problem and use it for 
    prob_bk = BK.getprob(br_codim2)

    prob_mpc = ManifoldProblem_BK(
                prob_bk, 
                br_codim2.sol[1].x,
                lens2, lens3;
                record_from_solution,
                finalize_solution
                )
    continuation(prob_mpc, alg, contparams)
end