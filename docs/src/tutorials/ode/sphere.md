# 🟢 Sphere

```@contents
Pages = ["sphere.md"]
Depth = 2
```

In this tutorial, we show how to cover a sphere as solution of

$$F(u) := \|u\|^2-1 = 0$$


We use this model as a mean to introduce the basics of `MultiParamContinuation.jl`.

It is easy to encode the manifold as follows

```@example TUTSPHERE
using GLMakie, MultiParamContinuation
const MPC = MultiParamContinuation

F(u,p) = [u[1]^2 + u[2]^2 + u[3]^2 - 1]

prob = ManifoldProblem(F, 
                    [1,0.,0.],
                    nothing
                        )
```

```@example TUTSPHERE
show(prob)
```

We now compute the covering of the manifold

```@example TUTSPHERE
S = MPC.continuation(prob,
            Henderson(),
            CoveringPar(max_charts = 20000, 
                    max_steps = 250,
                    verbose = 0,
                    newton_options = NonLinearSolveSpec(;maxiters = 5),
                    R0 = .2,
                    )
            )
show(S)
```

You plot the result as follows

```@example TUTSPHERE
MPC.plotd(S; 
    draw_tangent = true, 
    plot_center = false,
    draw_edges = true,
    ind_plot = [1,3])
```
