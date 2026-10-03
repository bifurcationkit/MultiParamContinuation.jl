# Manifold problem for use with BifurcationKit

```@contents
Pages = ["BifProblemBK.md"]
Depth = 3
```

`MultiParamContinuation.jl` is based on newton algorithm which relies on `NonlinearSolve.jl` for the implementation. One can chose to rely on `BifurcationKit.jl` newton method.

This can be done by calling `ManifoldProblem_BK` which has the same arguments as [`ManifoldProblem`](@ref).

## Jacobian of the two-parameter problem

When two parameter axes `lens1` and `lens2` are specified (they must be different),
`ManifoldProblem_BK` builds the composite problem `Z = [u; p1; p2]`. The computation
of its jacobian can be controlled with the `jacobian` keyword:

- `jacobian = nothing` (default): the state block is the jacobian of the underlying
  `BifurcationProblem` and the two parameter blocks are computed by finite differences.
- `jacobian = BK.AutoDiffDense()`: full dense jacobian by automatic differentiation.
- `jacobian = BK.AutoDiffMF()`: matrix-free jacobian-vector product (jvp) by automatic
  differentiation. Use a matrix-free linear solver, e.g. `GMRESIterativeSolvers`.
- `jacobian = BK.FiniteDifferencesMF()`: matrix-free jacobian-vector product by finite
  differences.
- `jacobian = BK.MatrixFree()`: matrix-free jacobian-vector product. The state block is
  the jacobian of the underlying problem (which should itself be matrix-free to fully
  avoid dense allocations) and the two parameter blocks are added explicitly.

!!! warning "Two-parameter problems are 2-manifolds"
    The composite problem `Z = [u; p1; p2]` assumes `n - m = 2`, i.e. the state
    dimension equals the number of equations: `length(F(u, p)) == length(u)`.

## Matrix-free problems

When the state dimension `n` is large (e.g. a discretized PDE), it is desirable to
avoid forming dense jacobian matrices. The `ManifoldProblemBKMatrixFree` problem type
only requires jacobian-vector products (jvp) throughout the covering algorithm. Three
ingredients must be set up:

1. **Jacobian**: `jacobian = BK.MatrixFree()` (or `BK.AutoDiffMF()`,
   `BK.FiniteDifferencesMF()`).
2. **Tangent space**: the default tangent computation solves a dense bordered system,
   hence matrix-free problems must use
   `get_tangent = BLSBorderedTangent(BK.MatrixFreeBLS(BK.GMRESIterativeSolvers()))`.
   This is also required to compute the curvature with `use_curvature = true`.
3. **Projection**: the Newton projection of a new chart is matrix-free as well; it
   needs a Krylov linear solver passed via `CoveringPar(solver_bls = ...)` (or
   `newton_options.linsolver`), e.g. `solver_bls = BK.GMRESIterativeSolvers()`.

As an example, we cover the manifold defined by

$$\lVert u\rVert^2 = a,\quad u_1 = b,\quad u_2 = c,$$

which has 3 equations for 3 state unknowns (so `n - m = 2`):

```@example TUTMF
using BifurcationKit, MultiParamContinuation
const MPC = MultiParamContinuation
const MPCExt = Base.get_extension(MultiParamContinuation, :BifurcationKitExt)

F(u, p) = [u[1]^2 + u[2]^2 + u[3]^2 - p.a, u[1] - p.b, u[2] - p.c]
par = (a = 2., b = 1., c = 0.)
prob_bk = BifurcationProblem(F, [1., 0., 1.], par, (@optic _.a))

prob = MPC.ManifoldProblemBKMatrixFree(prob_bk, [1., 0., 1.], (@optic _.a), (@optic _.b);
            jacobian = BK.MatrixFree(),
            get_tangent = MPCExt.BLSBorderedTangent(BK.MatrixFreeBLS(BK.GMRESIterativeSolvers())))
show(prob)
```

We can now compute a covering of the manifold, with a curvature based radius
adaptation:

```@example TUTMF
S = MPC.continuation(prob,
            Henderson(np0 = 6, use_curvature = true),
            CoveringPar(max_charts = 20000,
                    max_steps = 250,
                    newton_options = NewtonPar(),
                    solver_bls = BK.GMRESIterativeSolvers(),
                    Rmax = 0.2,
                    R0 = 0.2,
                    ϵ = 0.01,
                    )
            )
length(S)
```

For a problem without parameters, a matrix-free jacobian can be passed directly with
the `J` keyword of `ManifoldProblemBKMatrixFree`, e.g. for the unit sphere:

```@example TUTMF
Fs(u, p) = [u[1]^2 + u[2]^2 + u[3]^2 - 1]
Jmf = (x, p) -> dx -> [2 * dot(x, dx)]

prob_s = MPC.ManifoldProblemBKMatrixFree(Fs, [1., 0., 0.], nothing;
                    J = Jmf,
                    get_tangent = MPCExt.BLSBorderedTangent(BK.MatrixFreeBLS(BK.GMRESIterativeSolvers())))
show(prob_s)
```

```@docs
MultiParamContinuation.ManifoldProblem_BK
```
