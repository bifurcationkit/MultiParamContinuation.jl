# 🟢 ABC problem

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

opts_br = ContinuationPar(p_max = 1.5, n_inversion = 8, nev = 3)
br = BK.continuation(prob_bk, PALC(), opts_br; normC = norminf)
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

function plot_data!(ax, S; ind = (3,4,5),ind_col = 1, fil=x->true, cols = [real(c.index) for c in filter(x -> fil(x.u), S.atlas)])
    pts = mapreduce(c->[c.u[ind[1]], c.u[ind[2]], c.u[ind[3]]]', vcat, filter(x -> fil(x.u), S.atlas))
    hm = scatter!(ax, pts, color = cols)
end

plot_data(S_eq)
```

## Surface of Hopf points as function of 3 parameters

```@example TUTABC
opts_cover = CoveringPar(max_charts = 1000,
    max_steps = 600,
    # verbose = 1,
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
```

```@example TUTABC
fig = Figure()
ax = Axis3(fig[1,1], zlabel = "β", xlabel = "D", ylabel = "α", title = "Hopf surface $(length(atlas_hopf)) charts")
plot_data!(ax, atlas_hopf, ind = (4,5,6))
fig
```

## Surface of Periodic orbits as function of 2 parameters

We trace the curve of periodic orbits from a Hopf point: 

```@example TUTABC
argspo = (record_from_solution = (x, p; k...) -> begin
                xtt = BK.get_periodic_orbit(p.prob, x, p.p)
                return (max = maximum(xtt[3,:]), min = minimum(xtt[3,:]), period = x[end])
            end,
    plot_solution = (ax, x, p; k...) -> begin
        xtt = BK.get_periodic_orbit(p.prob, x, p.p)
        lines!(ax, xtt.t, xtt.u[1,:]; label = "u1", linewidth = 2)
        lines!(ax, xtt.t, xtt.u[2,:]; label = "u2")
        BK.plot!(get(k, :ax1, nothing), br)
    end,)

# continuation parameters
opts_po_cont = ContinuationPar(dsmax = 0.03, dsmin = 1e-4, ds = 0.0005, max_steps = 130, tol_stability = 4e-2, plot_every_step = 20)

br_po = BK.continuation(
    br, 1, opts_po_cont,
    PeriodicOrbitOCollProblem(50, 4; update_section_every_step = 1, jacobian = BK.DenseAnalyticalInplace());
    δp = 0.0001,
    argspo...,
    normC = norminf)
```

We can now compute the manifold

```@example TUTABC
using LinearAlgebra
const coll = br_po.prob

prob = MPC.ManifoldProblem_BK(
                        br_po.prob, br_po.sol[1].x, (@optic _.D), (@optic _.β),
                        record_from_solution = (X, p; k...) -> begin
                            xtt = BK.get_periodic_orbit(coll, X[1:end-2], nothing)
                            Max = maximum(xtt[3,:])
                            return (u3 = Max, D = X[end-1], β = X[end], period = X[end-2])
                        end,
                        finalize_solution = (X,p) -> begin
                            D = X[end-1]
                            β = X[end]
                            keep = (0.1 <= D <= 0.5) && (1.5 <= β <= 1.65)
                            return keep #&& norm(X) < 10
                        end,
                        )

S_po = @time MPC.continuation(prob,
                        Henderson(np0 = 6,
                                  θmin = 0.001,
                                  use_curvature = true,
                                  ),
                        CoveringPar(max_charts = 20000,
                                max_steps = 200,
                                verbose = 1,
                                newton_options = NewtonPar(tol = 1e-10, verbose = false),
                                R0 = .95,
                                ϵ = 0.4,
                                delta_angle = 4pi,
                                )
                        )
```

```@example TUTABC
fig = Figure()
ax3 = Axis3(fig[1,1], zlabel = "u3", xlabel = "D", ylabel = "β", title = "PO $(length(S_po)) charts")
pts = mapreduce(c->[c.data[1], c.data[3], c.data[4]]', vcat, S_po.atlas)
scatter!(ax3, pts, color = [c.data[2] for c in S_po.atlas])
fig
```