abstract type AbstractManifoldProblem end

for op in (:ManifoldProblem, :ManifoldProblemBK)
    @eval begin
    """
    $TYPEDEF

    Define a problem to perform multi-parameters continuation of a manifold defined as the zeros of F: Rⁿ → Rᵐ with n > m >= 1.

    ## Fields

    $TYPEDFIELDS
    
    ## Constructor

    ```
    ManifoldProblem(F, u0, par;
                    m = length(F(u0, par)),
                    check_dim::Bool = true,
                    recordFromSolution = (u,p) -> nothing,
                    project = nothing,
                    get_radius = get_radius_default,
                    get_tangent = nothing,
                    event_function = event_default,
                    finalize_solution = finalize_default,
                    project_for_tree = project_for_tree_default
                    )
    ```
    """
    struct $op{Tu <: AbstractVector, Tp, TVF, Trec, Tproj, Ttangent, Tradius, Tevent, Tfinalize, Tbb, Tpc} <: AbstractManifoldProblem
        "[Internal] input space dimension"
        n::Int
        "[Internal] output space dimension"
        m::Int
        "Equation representing the mapping F"
        VF::TVF
        "Guess for the initial point on the manifold"
        u0::Tu
        "Parameters passed to F"
        params::Tp
        "Record a few indicators at each chart of the manifold"
        recordFromSolution::Trec
        "Function to project a point from a tangent space to the manifold. If not provided, a newton algorithm is used. The signature is `project(u, par)` and returns a vector of solutions"
        project::Tproj
        "Compute an orthonormal basis of the tangent space at a point u on the manifold. Return a matrix of dimension n x (n-m). The signature is `get_tangent(u, par)`. If not provided, a dedicaded function is used."
        get_tangent::Ttangent
        "Get the an estimate of the curvature at a point u on the manifold. If not provided, a dedicaded function is used."
        get_radius::Tradius
        "Event function"
        event_function::Tevent
        "Finalise solution. Function to accept or not the current chart. It has signature `finalise(c::Chart, par)::Bool`"
        finalize_solution::Tfinalize
        "Function  used to project a point for the tree used to find the charts near a new point. Needs not be linear but the dimension should be at least the manifold embedding dimension."
        project_for_tree::Tbb
        "[Internal] constrained problem for projecting on manifold"
        prob_cons::Tpc

    end

    Base.size(prob::$op) = (prob.n, prob.m)
    @inline Base.eltype(prob::$op{Tu}) where Tu = eltype(Tu)
    @inline _has_projection(::$op{Tu, Tp, TVF, Trec, Tproj}) where {Tu, Tp, TVF, Trec, Tproj} = ~(Tproj == Nothing)
    @inline _has_tangent_computation(::$op{Tu, Tp, TVF, Trec, Tproj, Ttangent}) where {Tu, Tp, TVF, Trec, Tproj, Ttangent} = ~(Ttangent == Nothing)
    @inline _has_event(::$op{Tu, Tp, TVF, Trec, Tproj, Ttangent, Tradius, Tevent}) where {Tu, Tp, TVF, Trec, Tproj, Ttangent, Tradius, Tevent} = ~(Tevent == Nothing)

    function $op(F, u0, par;
                    m = length(F(u0, par)),
                    check_dim::Bool = true,
                    recordFromSolution = (u,p) -> nothing,
                    project = nothing,
                    get_radius = get_radius_default,
                    get_tangent = nothing,
                    event_function = event_default,
                    finalize_solution = finalize_default,
                    project_for_tree = project_for_tree_default,
                    prob_cons = nothing
                    )
        n = length(u0)
        dim = n-m
        if check_dim
            @assert n > m "This does not define an immersed manifold n = $n, m = $m"
        end
        $op(n,
            m,
            F,
            u0,
            par,
            recordFromSolution,
            project,
            get_tangent,
            get_radius,
            event_function,
            finalize_solution,
            project_for_tree,
            prob_cons)
        end
    end
end

@inline get_radius_default(u, p) = 1
@inline finalize_default(u,p) = true
@inline event_default(u,p) = nothing
@inline project_for_tree_default(u,p) = nothing

# empty function for BifurcationKit extension
function ManifoldProblem_BK end

function Base.show(io::IO, prob::AbstractManifoldProblem; prefix = "")
    n, m = size(prob)
    println(prefix * "$(n-m)-d Manifold Problem")
    println(prefix * "    ├─ n = ", prob.n)
    println(prefix * "    └─ m = ", prob.m)
end

function jacobian(prob::ManifoldProblem, u, p)
    return ForwardDiff.jacobian(x -> prob.VF(x, p), u)
end

function d2F(prob, u0, parms, dx1, dx2)
    d1Fad(x,p,dx1) = ForwardDiff.derivative(t -> prob.VF(x .+ t .* dx1, p), zero(eltype(dx1)))
    ForwardDiff.derivative(t -> d1Fad(u0 .+ t .* dx2, parms, dx1), zero(eltype(dx1)))
end

"""
Compute a basis for the tangent space at point `u0` on F(u, par) = 0.
If J is the jacobian dF(u0, par), it is found by solving

┌   ┐     ┌      ┐
│ J │ Φ = │  0   │
│ T │     │I(n-m)│
└   ┘     └      ┘

where T is a random matrix.
"""
function get_tangent(prob, u0, par, RHS)
    if _has_tangent_computation(prob)
        T = prob.get_tangent(u0, par)
    else
        J = jacobian(prob, u0, par)
        n, m = size(prob)
        _A = vcat(J, rand(n-m, n))
        T = _A \ RHS
        _A = vcat(J, T')
        T = _A \ RHS
        fact = qr(T)
        T = Matrix(fact.Q)
    end
end

function project(prob, u0, par)
    prob.project(u0, par)
end

function get_curvature(prob, c::Chart{Tu}, par) where {T, Tu <: AbstractVector{T}}
    u0 = c.u
    Φ = c.Φ
    n, m = size(prob)
    d = n-m
    J = jacobian(prob, u0, par)
    _A = vcat(J, Φ')
    A = zeros(T, d, n, d)
    for i in Base.OneTo(d)
        Φi = Φ[:, i]
        for j in Base.OneTo(d)
            Φj = Φ[:, j]
            d2 = d2F(prob, u0, par, Φi, Φj)
            A[i, :, j] .= _A \ vcat(d2, zeros(d))
        end
    end
    norm(A)
end