# 🟢 Torus

In this tutorial, we show how to cover a torus.

It is easy to encode the manifold as follows

```@example TUTTORUS
using GLMakie, MultiParamContinuation
const MPC = MultiParamContinuation

function F(u,p)
    x,y,z = u
    R2 = 0.5
    R1 = 0.8
    return [(sqrt(x^2 + y^2) - R1)^2 + z^2 - R2^2]
end

prob = ManifoldProblem(F, [-0.05, 1.29, -0.05], nothing)
```

```@example TUTTORUS
show(prob)
```

We now compute a covering of the manifold 

```@example TUTTORUS
S = MPC.continuation(prob,
            Henderson(np0 = 8),
            CoveringPar(max_charts = 2000, 
                    max_steps = 16000,
                    verbose = 0,
                    R0 = .1,
                    ϵ = Inf,
                    delta_angle = Inf,
                    ))
show(S)
```

You plot the result as follows

```@example TUTTORUS
f = MPC.plotd(S; 
    draw_tangent = true, 
    draw_edges = true,
    plot_center = false,
    ind_plot = [1,3])
```
