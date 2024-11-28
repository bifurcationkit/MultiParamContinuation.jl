"""
$TYPEDEF

Struct to specify the nonlinear solver and its arguments (options).

## Fields

- `nl_solver` nonlinear solver, defaults to `NewtonRaphson()`.
- `options` keyword arguments passed to the nonlinear solver. Example `abstol = 1e-12`.

## constructor

```
NonLinearSolveSpec(alg = NewtonRaphson(); kwargs...)
```
"""
struct NonLinearSolveSpec{T1, T2}
    nl_solver::T1
    options::T2
end

function NonLinearSolveSpec(alg = NewtonRaphson(); kwargs...)
    return NonLinearSolveSpec(alg, kwargs)
end

"""

$SIGNATURES

Parameters for the covering algorithm(s).

## Fields

$TYPEDFIELDS
"""
@with_kw struct CoveringPar{T, Tnl, Tbls}
    "Maximum distance between the tangent plane and the manifold."
    ϵ::T = 0.1
    "Initial radius of validity."
    R0::T = 0.1
    "[Internal] Minimal radius of polyhedra."
    Rmin::T = 0.001
    "Maximum angle difference between charts' tangent spaces in radians."
    delta_angle::T = 2π
    "Maximum number of charts."
    max_charts::UInt = 100
    "Maximum number of continuation steps. Because of mesh adaptation or failure, the number of computed charts is less or equal than `max_steps`."
    max_steps::UInt = 1000
    "Verbose mode, belongs to {0,1,2}. verbose = 0 prints nothing. verbose = 1 prints the charts, verbose = 2 print the intersection of the charts."
    verbose::Int = false
    "Newton options."
    newton_options::Tnl = NonLinearSolveSpec()
    "Bordered Linear Solver"
    solver_bls::Tbls = nothing
    "[Internal]."
    dotmin::T = 0.2
    @assert ϵ > 0 && R0 > 0 && Rmin > 0
    @assert verbose in 0:3
end