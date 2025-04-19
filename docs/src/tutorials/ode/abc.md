# 🟢 ABC problem

```@contents
Pages = ["abc.md"]
Depth = 3
```

The goal of this tutorial is to show how to compute a Manifold from a `BifurcationProblem` as function of two free parameters.

```@example TUTABC
using CairoMakie
using BifurcationKit, MultiParamContinuation

const BK = BifurcationKit
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
par_abc = (D = 0.11, B = 8., α = 1., σ = 0.04, β = 1.56)
z0 = [1., 0., 0. ]
prob_bk = BifurcationProblem(abc!, z0, par_abc, (@optic _.D), 
        record_from_solution = (x, p; k...) -> (u3 = x[3], u1 = x[1], u2 = x[2]),)

opts_br = ContinuationPar(p_max = .5, n_inversion = 8, nev = 3)
br = BK.continuation(prob_bk, PALC(), opts_br; normC = norminf)

fig, ax = BK.plot(br, plotfold=false, )
fig
```

We can also compute the stationary points as function of two free parameters:


```@example TUTABC
prob = MPC.ManifoldProblem_BK(
                        prob_bk, br.sol[1].x, (@optic _.D), (@optic _.β);
                        record_from_solution = (X, p; k...) -> begin
                            return (β = X[end], D = X[end-1], u3 = X[3])
                        end,
                        finalize_solution = (X,p) -> begin
                            D = X[end-1]
                            β = X[end]
                            keep = (0.01 <= D <= 0.5) && (1.5 <= β <= 1.65)
                            return keep
                        end,
                        )

S_eq = @time MPC.continuation(prob,
                        Henderson(np0 = 3,
                                    θmin = 0.001,
                                    # use_curvature = true,
                                    use_tree = true,
                                  ),
                        CoveringPar(max_charts = 20000,
                                max_steps = 1000,
                                verbose = 0,
                                newton_options = NewtonPar(tol = 1e-10),
                                R0 = .04,
                                ϵ = 0.1,
                                delta_angle = 10.15,
                                ))
```

You can plot the data as follows

```@example TUTABC
function plot_data(S; k...)
    fig = Figure()
    ax = Axis3(fig[1,1], zlabel = "u3", xlabel = "D", ylabel = "β", title = "$(length(S)) charts")
    plot_data!(ax, S; k...)
    fig
end

function plot_data!(ax, S; ind_col = 1, fil=x->true, cols = [real(c.index) for c in filter(x -> fil(x.u), S.atlas)])
    pts = mapreduce(c->[c.data[2], c.data[3], c.data[1]]', vcat, filter(x -> fil(x.u), S.atlas))
    cols = [real(c.data[ind_col]) for c in filter(x -> fil(x.u), S.atlas)]
    hm = scatter!(ax, pts, color = cols)
end

plot_data(S_eq)
```