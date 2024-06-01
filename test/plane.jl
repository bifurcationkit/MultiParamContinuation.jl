using MultiParamContinuation, StaticArrays

const MPC = MultiParamContinuation

F(u,p) = [u[3]]

prob = ManifoldProblem(F, SA[0.,0,0], nothing;
                        get_tangent = (u,p) -> SA[1 0; 0 1; 0 0],
                        finalize_solution = Cube(0.5),
                        project = (out,p) -> begin
                            SVector(out[1], out[2], 0)
                        end
                        )

cpar = CoveringPar(max_charts = 200_000, 
            verbose = 0,
            max_steps = 4_000,
            R0 = .1,
            ϵ = Inf,
            delta_angle = Inf,
            )

S = continuation(prob,
        Henderson(np0 = 5),
        cpar)

@test S[1] isa MPC.Chart
@test length(S) == 109

S = continuation(prob,
        Henderson(np0 = 5,
                  use_tree = true,
                  children_pre_leaf = 5,
                  ),
        cpar
        )

@test S[1] isa MPC.Chart
@test length(S) == 109


prob = ManifoldProblem(F, SA[0.,0,0], nothing;
                        get_tangent = (u,p) -> SA[1 0; 0 1; 0 0],
                        finalize_solution = ProductSpace(fill(-0.5,3), fill(0.5,3)),
                        project = (out,p) -> begin
                            SVector(out[1], out[2], 0)
                        end
                        )
S = continuation(prob,
        Henderson(np0 = 5),
        cpar)

@test S[1] isa MPC.Chart
@test length(S) == 109

show(S)
show(S[1])
