# Manifold problems

```@contents
Pages = ["BifProblem.md"]
Depth = 3
```

The idea behind `MultiParamContinuation` is to compute immersed manifolds $\mathcal M$ in memory limited environments, the manifold being defined as (parts) of the zeros of 

$$F:\mathbb R^n\times \mathbb R^{par}\to\mathbb R^m.$$

## Generic manifold problem

[`ManifoldProblem`](@ref) is the basic / generic structure for encoding a bifurcation problem ; it holds the following fields:

- the vector field
- an initial guess
- a set of parameters

as well as user defined functions for 

- recording (`record_from_solution`) indicators about the solution when this one is too large to be saved at every continuation step.

and some other that we described below. A detailed account for the struct is [`ManifoldProblem`](@ref).

!!! tip "Tutorial"
    The [Plane](@ref plane) tutorial provides an example where these arguments are used in situation.

### Basic example

```@example TUTPROB
using MultiParamContinuation

F(u,p) = [u[3]]

prob = ManifoldProblem(F, [0,0,0.], nothing)
```

### Projection function

You can pass your own projection function which takes a point $u_{guess}\in \mathbb R^n$ and project it on $\mathcal M$. For example, in the case of the plane `F(u,p) = u[3]`, we could use:

```@example TUTPROB
prob = ManifoldProblem(F, [0,0,0.], nothing;
                        project = (u_g,p) -> begin
                        		u = copy(u_g)
                        		u[3] = 0
                        		return u
                        end,
                        )
```

### Tangent function

You can pass your own tangent function which returns an orthonormal basis of the tangent space at point $u$. For example, in the case of the plane `F(u,p) = u[3]`, we could use:

```@example TUTPROB
prob = ManifoldProblem(F, [0,0,0.], nothing;
                        get_tangent = (u,p) -> [1 0; 0 1; 0 0],
                        )
```

### Finalizing the solution

You can pass a finalizer function which returns a boolean. This boolean is used to discard or not the current chart. This can be used to restrict the manifold to a given bounding box.

You can have a look at [bounding space](@ref bouding-space) for predefined finalizer functions for simple spaces.

```@example TUTPROB
prob = ManifoldProblem(F, [0,0,0.], nothing;
                        finalize_solution = (c::Chart, p) -> true
                        )
```