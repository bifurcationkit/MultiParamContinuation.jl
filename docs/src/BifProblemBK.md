# Manifold problem for use with BifurcationKit

```@contents
Pages = ["BifProblemBK.md"]
Depth = 3
```

`MultiParamContinuation.jl` is based on newton algorithm which relies on `NonlinearSolve.jl` for the implementation. One can chose to rely on `BifurcationKit.jl` newton method.

This can be done by calling `ManifoldProblem_BK` which has the same arguments as [`ManifoldProblem`](@ref).


```@docs
MultiParamContinuation.ManifoldProblem_BK
```