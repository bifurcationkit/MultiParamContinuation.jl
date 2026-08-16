# 🟢 Parabola

In this tutorial, we show how to cover a parabola.

```@example TUTTORUS
using CairoMakie, MultiParamContinuation
const MPC = MultiParamContinuation

F(u,p) = [u[1]^12 + u[2]^12 - u[3]]

prob = ManifoldProblem(F, [0,0,0.], nothing;
    finalize_solution = (u,p) -> (u[3] < 2) * (u[1]>-0.1) * (u[2]>-0.1))
```

```@example TUTTORUS
show(prob)
```

We now compute a covering of the manifold 

```@example TUTTORUS
S = continuation(prob,
            Henderson(
                      use_curvature = true,
                      ),
            CoveringPar(max_charts = 10000,
                    max_steps = 2000,
                    verbose = 0,
                    newton_options = NonLinearSolveSpec(;maxiters = 5, abstol = 1e-12, reltol = 1e-10),
                    Rmax = .2, # maximal radius of validity
                    R0 = .01, # initial radius of validity
                    ϵ = 0.005, # used to estimated current radius
                    )
            )
show(S)
```

You plot the result as follows

```@example TUTTORUS
f = MPC.plotd(S; draw_edges = true,)
```
