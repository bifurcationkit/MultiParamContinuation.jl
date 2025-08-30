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
                        options_cont_bk::ContinuationPar = br.contparams,
                        finalize_solution = (X,p) -> true,
                        jacobian_ma = BK.MinAugMatrixBased()
                        )
    if lens3 == lens2 || BK.getlens(BK.getprob(br)) in (lens2, lens3)
        error("You must pass different lenses")
    end
    # we compute one step with continuation to get the problem and initial point
    @reset options_cont_bk.max_steps = 1
    @reset options_cont_bk.detect_bifurcation = 2
    br_codim2 = BK.continuation(br, ind, lens2, options_cont_bk; 
                        jacobian_ma,
                        verbosity = 3)
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