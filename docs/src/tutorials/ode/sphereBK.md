# 🟢 Sphere based on BifurcationKit.jl

In this tutorial, we show how to cover a sphere as solution of

$$F(u) := \|u\|^2-1 = 0$$

We use this model as a mean to introduce the basics of `MultiParamContinuation.jl` based on `BifurcationKit.jl`.

It is easy to encode the manifold using `ManifoldProblem_BK` which yields a `ManifoldProblemBK` with `BifurcationKit` internals.

```@example TUTSPHEREBK
using CairoMakie, BifurcationKit, MultiParamContinuation
const MPC = MultiParamContinuation

F(u,p) = [u[1]^2 + u[2]^2 + u[3]^2 - 1]

prob = ManifoldProblem_BK(F, 
                    [1,0.,0.],
                    nothing
                        )
```

```@example TUTSPHEREBK
show(prob)
```

We compute a covering of the manifold: 

```@example TUTSPHEREBK
S = MPC.continuation(prob,
            Henderson(np0 = 6),
            CoveringPar(max_charts = 20000, 
                    max_steps = 250,
                    # trigger the use of BifurcationKit newton solver
                    newton_options = NewtonPar(),
                    R0 = .2,
                    )
            )
show(S)
```

You can now plot the result as

```@example TUTSPHEREBK
MPC.plotd(S; 
    draw_tangent = true, 
    plot_center = false,
    draw_edges = true,
    ind_plot = [1,3])
```
