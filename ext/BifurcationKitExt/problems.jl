using ForwardDiff
record_from_solution_nothing(x, p; k...) = nothing

function correct_guess(cache, options::BK.NewtonPar)
    prob = cache.prob
    sol = BK.solve(prob.VF, BK.Newton(), options; normN = BK.norminf)
    if ~BK.converged(sol)
        throw("Newton for first point did not converge!!")
    end
    return sol.u
end

for (M, OP) in ((:ManifoldProblem_BK, :ManifoldProblemBK),
                (:ManifoldProblemBKMatrixFree, :ManifoldProblemBKMatrixFree))

    @eval begin
        """
        $SIGNATURES

        Make a manifold problem from a vector field `F`. The zeros of `F: Rⁿ → Rᵐ`
        (with `n > m`) define an `n-m`-d manifold. The mapping is wrapped into a
        `BifurcationProblem`; the same keyword arguments as in the method taking a
        `BifurcationProblem` are accepted.

        The jacobian of the underlying `BifurcationProblem` can be provided with the
        `J` keyword, either as a matrix, as a function `(x, p) -> J` or as a
        matrix-free jacobian-vector product `(x, p) -> dx -> J⋅dx`. For a matrix-free
        problem, also provide `get_tangent = BLSBorderedTangent(bls)` (see
        [`ManifoldProblemBKMatrixFree`](@ref)).
        """
        function $M(F, u0::AbstractVector, par;
                    J = nothing,
                    check_dim::Bool = true,
                    record_from_solution = record_from_solution_nothing,
                    project = nothing,
                    get_radius = get_radius_default,
                    get_tangent = nothing,
                    event_function = event_default,
                    finalize_solution = finalize_default,
                    project_for_tree = project_for_tree_default,
                    prob_cons = nothing,
                    weights = Weight(TrivialWeight()),
                    )
            prob_bk = BifurcationProblem(F, u0, par, (@optic _); J)
            if isnothing(prob_cons)
                𝒯 = eltype(u0)
                m = length(BK.residual(prob_bk, prob_bk.u0, par))
                Φ = zeros(𝒯, m+2, 2)
                wbar = zeros(𝒯, m+2)
                prob_cons = ConstrainedProblem(prob_bk, Φ, Φ' * wbar, wbar)
            end
            m = length(BK.residual(prob_bk.VF, u0, par))
            _make_manifold_problem($OP, prob_bk, u0, par, m;
                        check_dim, record_from_solution, project, get_radius,
                        get_tangent, event_function, finalize_solution, project_for_tree,
                        prob_cons, weights)
        end

        """
        $SIGNATURES

        Make a manifold problem from a `BifurcationProblem` and specifying two parameter axes.

        The unknowns are the composite vector `Z = [u; p1; p2]`, where `p1, p2` are the
        values of the two parameter axes given by the lenses `lens1, lens2`. The state
        dimension must equal the number of equations, `length(u) == length(F(u, p))`,
        i.e. `n - m = 2`.

        The jacobian of the composite problem `Z = [u; p1; p2]` can be selected with the
        `jacobian` keyword. It defaults to `nothing` (jacobian of `prob_bk` for the state
        block and finite differences for the two parameter blocks). Otherwise pass
        a `BifurcationKit` jacobian marker, e.g. `BK.AutoDiffDense()`, `BK.AutoDiffMF()`,
        `BK.MatrixFree()` or `BK.FiniteDifferencesMF()`.

        For the matrix-free markers `BK.AutoDiffMF()`, `BK.FiniteDifferencesMF()` and
        `BK.MatrixFree()`, the tangent space must be computed with a bordered linear
        solver: pass `get_tangent = BLSBorderedTangent(bls)`, e.g.
        `BLSBorderedTangent(BK.MatrixFreeBLS(BK.GMRESIterativeSolvers()))`. For
        `BK.MatrixFree()`, the state block is the jacobian of `prob_bk`, which should
        itself be matrix-free to fully avoid dense allocations.

        The metric of the embedding space can be changed with the `weights` keyword, see
        [`ManifoldProblem`](@ref).
        """
        function $M(prob_bk::BK.AbstractBifurcationProblem,
                    u0::AbstractVector,
                    lens1,
                    lens2;
                    check_dim::Bool = true,
                    project = nothing,
                    get_radius = get_radius_default,
                    get_tangent = nothing,
                    record_from_solution = record_from_solution_nothing,
                    event_function = event_default,
                    finalize_solution = finalize_default,
                    project_for_tree = project_for_tree_default,
                    weights = Weight(TrivialWeight()),
                    jacobian = nothing,
                    )
            par = BK.getparams(prob_bk)
            m = length(BK.residual(prob_bk, prob_bk.u0, par))

            # make a bifurcation problem with two parameters axes
            pb_composite = BifurcationProblem_2P(prob_bk, lens1, lens2, jacobian)
            new_u0 = vcat(u0, BK._get(par, lens1), BK._get(par, lens2))

            prob_mpc = BifurcationProblem(pb_composite,
                            new_u0,
                            par,
                            (@optic _);
                            J = (x, p) -> _jacobian_2P(pb_composite, pb_composite.jacobian, x, p)
                            )
            𝒯 = eltype(new_u0)
            Φ = zeros(𝒯, m+2, 2)
            wbar = zeros(𝒯, m+2)
            prob_cons = ConstrainedProblem(prob_mpc, Φ, Φ' * wbar, wbar)

            _make_manifold_problem($OP, prob_mpc, new_u0, par, m;
                        check_dim, record_from_solution, project, get_radius,
                        get_tangent, event_function, finalize_solution, project_for_tree,
                        prob_cons, weights)
        end
    end
end

jacobian(pb::AbstractManifoldProblemBifurcationKit, u, p) = BK.jacobian(pb.VF, u, p)
BK.residual(pb::AbstractManifoldProblemBifurcationKit, u, p) = BK.residual(pb.VF, u, p)
function d2F(prob::AbstractManifoldProblemBifurcationKit, x, p, dx1, dx2)
    pb = prob.VF.VF.F
    return pb isa BifurcationProblem_2P ? _d2F_2P(pb, x, p, dx1, dx2) :
        BK.d2F(prob.VF, x, p, dx1, dx2)
end
BK.getlens(::AbstractManifoldProblemBifurcationKit) = nothing

function _make_manifold_problem(::Type{OP}, F, u0, par, m;
                                check_dim::Bool = true,
                                record_from_solution = record_from_solution_nothing,
                                project = nothing,
                                get_radius = get_radius_default,
                                get_tangent = nothing,
                                event_function = event_default,
                                finalize_solution = finalize_default,
                                project_for_tree = project_for_tree_default,
                                prob_cons = nothing,
                                weights = Weight(TrivialWeight()),
                                ) where {OP}
    n = length(u0)
    if check_dim
        @assert n > m "This does not define an immersed manifold n = $n, m = $m"
    end
    OP(n, m, F, u0, par,
        record_from_solution, project, get_tangent, get_radius,
        event_function, finalize_solution,
        project_for_tree,
        prob_cons,
        MultiParamContinuation.update_default,
        weights)
end
#━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
"""
$TYPEDEF

Constrained problem used to project a point `w` back on the manifold. It is the mapping

    (pb::ConstrainedProblem)(w, p) = [F(w, p) ; Φ' (w - wbar)]

where `F` is the residual of the underlying bifurcation problem `prob`, `Φ` is an
orthonormal basis of the tangent space at the chart center `wbar`. The extra constraint
`Φ' (w - wbar) = 0` restricts `w` to the affine tangent plane; the projection onto the
manifold is obtained by solving `(pb::ConstrainedProblem)(w, p) = 0` for `w`.

## Fields

$TYPEDFIELDS
"""
struct ConstrainedProblem{T1, T2, T3, T4} <: BK.AbstractBifurcationProblem
    "Underlying bifurcation problem, typically a `BifurcationProblem` wrapping a `BifurcationProblem_2P`."
    prob::T1
    "Orthonormal basis of the tangent space at `wbar`, a matrix of size `n × (n-m)`."
    Φ::T2
    "Cached `Φ' * wbar`."
    Φwbar::T3
    "Center of the chart, a point on the manifold."
    wbar::T4
end

function (cpb::ConstrainedProblem)(w, p)
    # weights = get_weights(cpb.prob)
    # vcat(BK.residual(cpb.prob.VF, w, p), apply_T(weights, cpb.Φ, w - cpb.wbar))
    vcat(BK.residual(cpb.prob.VF, w, p), cpb.Φ' * (w - cpb.wbar))
end

function jacobian(cpb::ConstrainedProblem, w, p)
    J0 = BK.jacobian(cpb.prob, w, p)
    vcat(J0, cpb.Φ')
end

BK.residual(cpb::ConstrainedProblem, x, p) = cpb(x, p)
BK.jacobian(cpb::ConstrainedProblem, x, p) = jacobian(cpb, x, p)
BK.getu0(cpb::ConstrainedProblem) = cpb.prob.u0
BK.getparams(cpb::ConstrainedProblem) = cpb.prob.params
#━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
"""
$SIGNATURES

Bifurcation Problem with two parameter axes for which we continue the zeros.
"""
struct BifurcationProblem_2P{T1, T2, T3, T4}
    prob::T1
    lens1::T2
    lens2::T3
    "How to compute the jacobian of the composite problem. `nothing` (default) uses the jacobian of `prob` for the state block and automatic differentiation for the two parameter blocks. Otherwise a BifurcationKit jacobian marker such as `BK.AutoDiffDense()`, `BK.AutoDiffMF()`, `BK.MatrixFree()`, `BK.FiniteDifferences()`, ..."
    jacobian::T4
end

BifurcationProblem_2P(prob, lens1, lens2) = BifurcationProblem_2P(prob, lens1, lens2, nothing)

function (pb::BifurcationProblem_2P)(Z, par)
    u = @view Z[1:end-2]
    p1 = Z[end-1]
    p2 = Z[end]
    par2 = BK._set(par, (pb.lens1, pb.lens2), (p1, p2))
    BK.residual(pb.prob, u, par2)
end

# Parameter columns ∂F/∂p₁, ∂F/∂p₂ of the composite jacobian. The point is
# built with `@set` (through the `BifurcationProblem_2P` wrapper) instead of a `vcat`.
function _jacobian_param(pb::BifurcationProblem_2P, Z, par)
    p1, p2 = Z[end-1], Z[end]
    par2 = BK._set(par, (pb.lens1, pb.lens2), (p1, p2))
    # (ForwardDiff.derivative(z -> pb((@set Z[end-1] = z), par2), p1),
    #  ForwardDiff.derivative(z -> pb((@set Z[end] = z),   par2), p2))
    h = sqrt(eps(real(eltype(Z))))
    f0 = pb(Z, par2)
    ∂1 = (pb((@set Z[end-1] = p1 + h), par2) .- f0) ./ h
    ∂2 = (pb((@set Z[end] = p2 + h),   par2) .- f0) ./ h
    return (∂1, ∂2)
end

# Default: use the jacobian of the underlying BifurcationKit problem for the
# state block and automatic differentiation for the two parameter blocks.
function _jacobian_2P(pb::BifurcationProblem_2P, ::Nothing, Z, par)
    u = @view Z[1:end-2]
    p1 = Z[end-1]
    p2 = Z[end]
    par2 = BK._set(par, (pb.lens1, pb.lens2), (p1, p2))
    J0 = BK.jacobian(pb.prob, u, par2)
    l1, l2 = _jacobian_param(pb, Z, par)
    return hcat(J0, l1, l2)
end

# Dense jacobian obtained by automatic differentiation of the composite problem.
_jacobian_2P(pb::BifurcationProblem_2P, ::BK.AutoDiffDense, Z, par) =
    ForwardDiff.jacobian(z -> pb(z, par), Z)

# Matrix-free jacobian-vector product obtained by AD of the composite problem.
_jacobian_2P(pb::BifurcationProblem_2P, ::BK.AutoDiffMF, Z, par) =
    dx -> ForwardDiff.derivative(t -> pb(Z .+ t .* dx, par), zero(eltype(Z)))

# Matrix-free jacobian-vector product obtained by finite differences.
function _jacobian_2P(pb::BifurcationProblem_2P, ::BK.FiniteDifferencesMF, Z, par)
    h = sqrt(eps(real(eltype(Z))))
    return dx -> (pb(Z .+ h .* dx, par) .- pb(Z, par)) ./ h
end

# Matrix-free jacobian-vector product of the composite problem. The state block uses
# the matrix-free jacobian (jvp) `J0 = jacobian(pb.prob, u, par2)` of the underlying
# BifurcationKit problem; the two parameter columns are added explicitly.
# Returns `dx -> J ⋅ dx`, with `dx` either a plain vector or a `BK.BorderedArray`.
function _jacobian_2P(pb::BifurcationProblem_2P, ::BK.MatrixFree, Z, par)
    u = @view Z[1:end-2]
    p1, p2 = Z[end-1], Z[end]
    par2 = BK._set(par, (pb.lens1, pb.lens2), (p1, p2))
    l1, l2 = _jacobian_param(pb, Z, par)
    J0 = BK.jacobian(pb.prob, u, par2)
    return dx -> begin
        d = dx isa BK.BorderedArray ? dx.u : dx
        BK.apply(J0, @view(d[1:end-2])) .+ l1 .* d[end-1] .+ l2 .* d[end]
    end
end

jacobian(pb::BifurcationProblem_2P, Z, par) = _jacobian_2P(pb, pb.jacobian, Z, par)

# Second derivative of the composite residual `pb(Z, par)` along `(dx1, dx2)`.
# - state-state block: analytic `BK.d2F` of the underlying problem;
# - mixed state/parameter blocks: finite differences of the underlying state jvp;
# - parameter-parameter block: finite differences of the parameter columns.
# The jvp may be stateful (e.g. built with `updateJac!`), so it is (re)created right
# before being applied.
function _d2F_2P(pb::BifurcationProblem_2P, Z, par, dx1, dx2)
    u = @view Z[1:end-2]
    p1, p2 = Z[end-1], Z[end]
    par2 = BK._set(par, (pb.lens1, pb.lens2), (p1, p2))
    a1 = @view dx1[1:end-2]
    a2 = @view dx2[1:end-2]
    α1 = dx1[end-1:end]
    α2 = dx2[end-1:end]
    𝒯 = eltype(Z)
    h = sqrt(eps(real(𝒯)))

    # ∂²F/∂u²
    d2 = BK.d2F(pb.prob, u, par2, a1, a2)

    # mixed blocks: derivative of the state jvp w.r.t. the parameters
    S = BK.jacobian(pb.prob, u, par2)
    y1 = BK.apply(S, a1)
    y2 = BK.apply(S, a2)
    par2_α2 = BK._set(par2, (pb.lens1, pb.lens2), (p1 + h * α2[1], p2 + h * α2[2]))
    S2 = BK.jacobian(pb.prob, u, par2_α2)
    m1 = (BK.apply(S2, a1) .- y1) ./ h
    par2_α1 = BK._set(par2, (pb.lens1, pb.lens2), (p1 + h * α1[1], p2 + h * α1[2]))
    S1 = BK.jacobian(pb.prob, u, par2_α1)
    m2 = (BK.apply(S1, a2) .- y2) ./ h

    # parameter-parameter block: derivative of the parameter columns
    l1, l2 = _jacobian_param(pb, Z, par)
    w0 = l1 .* α1[1] .+ l2 .* α1[2]
    Zα2 = copy(Z)
    Zα2[end-1] += h * α2[1]
    Zα2[end] += h * α2[2]
    l1b, l2b = _jacobian_param(pb, Zα2, par)
    pp = (l1b .* α1[1] .+ l2b .* α1[2] .- w0) ./ h

    return d2 .+ m1 .+ m2 .+ pp
end
#━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
function project_on_M(prob_bk, guess, chart::Chart, wbar, cpar::CoveringPar{T, <: BK.NewtonPar}, weights) where {T}
    # @error "project_on_M - BK"
    if _has_projection(prob_bk)
        return project(prob_bk, guess, prob_bk.params)
    else
        options = cpar.newton_options
        if ~isnothing(cpar.solver_bls)
            options = @set options.linsolver = cpar.solver_bls
        end

        prob_bls = prob_bk.prob_cons
        prob_bls.Φ .= chart.Φ
        prob_bls.wbar .= wbar
        prob_bls.prob.u0 .= guess

        # f(w, p) = vcat(BK.residual(prob_bk.VF, w, p), apply_T(weights, chart.Φ, w - wbar))
        # prob_bls = BifurcationProblem(f, guess, BK.getparams(prob_bk.VF))
        normN = Base.Fix1(weighted_norm, weights)
        sol = BK.solve(prob_bls, BK.Newton(), options; normN)
    end
    if BK.converged(sol)
        return sol.u
    else
        return nothing
    end
end

"""
$SIGNATURES

Matrix-free projection on the manifold for a `ManifoldProblemBKMatrixFree`. It solves
the projection Newton system

    [F(w) ; Φ' D (w - u₀) - ω] = 0

where `D` is the diagonal metric given by the `weights` and `ω = Φ' D (wbar - u₀)`
restricts the correction to the affine tangent plane, by providing a matrix-free
jacobian-vector product `dx -> [J⋅dx ; Φ' D ⋅dx]` to `BifurcationProblem`, where
`J = jacobian(prob, w, p)` may itself be matrix-free.
The linear solver `cpar.newton_options.linsolver` (or `cpar.solver_bls` when
provided) must be a Krylov method, e.g. `GMRESIterativeSolvers()`.
"""
function project_on_M(prob::ManifoldProblemBKMatrixFree, guess, chart::Chart, wbar, cpar::CoveringPar{T, <: BK.NewtonPar}, weights) where {T}
    @error "project_on_M BK"
    if _has_projection(prob)
        return project(prob, guess, prob.params)
    end
    options = cpar.newton_options
    if ~isnothing(cpar.solver_bls)
        options = @set options.linsolver = cpar.solver_bls
    end
    Φ = chart.Φ
    u₀ = chart.u
    ########
    # θ = 0.5 # small θ favours changes in u
    # β0 = θ / (length(u₀) - 2)
    # β1 = (1 - θ)/2
    # vβ = fill(β0, length(u₀)); vβ[end-1:end] .= β1
    # D = Diagonal(vβ)
    ########
    ω = apply_T(weights, Φ, wbar - u₀) # Φ' * D * (wbar - u₀)
    function f(w, p)
        # vcat(BK.residual(prob.VF, w, p), Φ' * (w - wbar))
        vcat(BK.residual(prob.VF, w, p), 
            # Φ' * D * (w - u₀) - ω)
            # apply_T(weights, Φ, w - wbar)
            apply_T(weights, Φ, w - u₀) - ω
            )
    end
    # matrix-free jacobian of [F ; Φ'] : dx -> [J⋅dx ; Φ'⋅dx]
    Jmf = (w, p) -> let Φ=Φ
            Jp = jacobian(prob, w, p)
            return dx -> vcat(BK.apply(Jp, dx), apply_T(weights, Φ, dx))#Φ' * D * dx)
    end
    prob_bls = BifurcationProblem(f, guess, BK.getparams(prob.VF); J = Jmf)
    sol = BK.solve(prob_bls, BK.Newton(), options; normN = BK.norminf)
    if BK.converged(sol)
        return sol.u
    else
        return nothing
    end
end
#━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
for BKP in (:ManifoldProblemBK, :ManifoldProblemBKMatrixFree)
    @eval begin
        function get_tangent(prob::$BKP{Tu, Tp, TVF, Trec, Tproj, BorderedTangent}, u0, par, RHS, Φ0 = nothing) where {Tu <: AbstractVector, Tp, TVF, Trec, Tproj}
            @assert false
            return _get_tangent_bordered(prob, u0, par, RHS, Φ0)
        end

        function get_tangent(prob::$BKP{Tu, Tp, TVF, Trec, Tproj, QRDirectTangent}, u0, par, RHS, Φ0 = nothing) where {Tu <: AbstractVector, Tp, TVF, Trec, Tproj}
            @assert false
            return _get_tangent_QR(prob, u0, par, RHS, Φ0)
        end
    end
end
#━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
"""
$TYPEDEF

Marker to compute the tangent space of a `ManifoldProblemBK` by solving the bordered
system with a `BifurcationKit` bordered linear solver `bls`. It enables matrix-free
jacobians, e.g. `bls = BK.MatrixFreeBLS(BK.GMRESIterativeSolvers())`. It is also
required by [`get_curvature`](@ref) for a matrix-free problem.

Pass it as the `get_tangent` field of the manifold problem:

    ManifoldProblem_BK(...; get_tangent = BLSBorderedTangent(bls))
"""
struct BLSBorderedTangent{Tbls} <: AbstractTangentAlgorithm
    "BifurcationKit bordered linear solver used to solve the bordered system."
    bls::Tbls
end

"""
$TYPEDSIGNATURES

Compute the tangent space by solving the bordered system with the `BifurcationKit`
bordered linear solver `bls`. The jacobian `J = jacobian(prob, u0, par)` may be
matrix-free (a jacobian-vector product). It is found by solving the
bordered system

┌   ┐     ┌      ┐
│ J │ Φ = │  0   │
│ T │     │I(n-m)│
└   ┘     └      ┘

where `T` is a random matrix. The returned basis is orthonormalized in the metric
defined by the `weights` of the problem (weighted Cholesky factorization).
"""
function _get_tangent_bordered_bls(prob::AbstractManifoldProblemBifurcationKit, u0, par, RHS, bls, Φ0 = nothing)
    J = jacobian(prob, u0, par)
    n, m = size(prob)
    k = n - m
    𝒯 = eltype(u0)
    Jmat = J isa AbstractMatrix
    # matrix/operator of J[:, 1:m]
    J1 = Jmat ? view(J, :, 1:m) : (v -> BK.apply(J, vcat(v, zeros(𝒯, k))))

    # columns of J2 = J[:, m+1:n]
    pb = prob.VF.VF.F
    dₚF = if pb isa BifurcationProblem_2P
        # parameter columns ∂F/∂p₁, ∂F/∂p₂ from the composite problem
        _jacobian_param(pb, u0, par)
    else
        # k jacobian-vector products
        unit(i) = (v = zeros(𝒯, n); v[i] = one(𝒯); v)
        ntuple(i -> BK.apply(J, unit(m + i)), k)
    end

    function bordered(B0)
        B = Matrix(B0)                                   # dense border, k × n
        b = ntuple(i -> @view(B[i, 1:m]), k)             # rows of B1
        c = B[:, (m+1):n]                                # B2
        Φ = zeros(𝒯, n, k)
        local converged = true
        for i in 1:k
            rhs_bottom = zeros(𝒯, k)
            rhs_bottom[i] = one(𝒯)
            u1, u2, cv, it = BK.solve_bls_block(bls, J1, dₚF, b, c, zeros(𝒯, m), rhs_bottom)
            converged = converged & cv
            Φ[:, i] .= vcat(u1, u2)
        end
        return Φ, converged
    end

    # initial border: random, or the tangent space of a nearby chart ("tangent with
    # guess") for continuity of the basis along the manifold
    w = get_weights(get_weights(prob))
    B0 = isnothing(Φ0) ? rand(real(𝒯), k, n) :
        (w isa TrivialWeight ? Matrix(Φ0') : Matrix(Φ0' * Diagonal(w)))
    Φ, cv = bordered(B0)
    cv || return nothing
    Φ, cv = bordered(Φ') # re-border with the tangent
    cv || return nothing
    T = MultiParamContinuation.weighted_orthonormalize(Φ, get_weights(prob)) # T' * D * T = I
    return T
end

for BKP in (:ManifoldProblemBK, :ManifoldProblemBKMatrixFree)
    @eval function get_tangent(prob::$BKP{Tu, Tp, TVF, Trec, Tproj, BLSBorderedTangent{Tbls}}, u0, par, RHS, Φ0 = nothing) where {Tu <: AbstractVector, Tp, TVF, Trec, Tproj, Tbls}
        return _get_tangent_bordered_bls(prob, u0, par, RHS, prob.get_tangent.bls, Φ0)
    end
end
#━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
"""
$TYPEDSIGNATURES

Curvature of the manifold at the chart `c` for a `ManifoldProblemBKMatrixFree`.

For each pair of tangent directions `(Φi, Φj)`, the bordered system
`[J(u0); Φ' D] X = [d2F(Φi,Φj); 0]`, where `D` is the diagonal metric given by the
`weights`, is solved with the bordered linear solver
`prob.get_tangent.bls` (see [`BLSBorderedTangent`](@ref)). The entry
`cmat[i, j]` is the weighted norm of the solution.

The composite second derivative is provided by `d2F(prob, …)`: the state-state block is
analytic (underlying `BK.d2F`), the mixed and parameter blocks are finite differences.
Because the jvp may be stateful, all `d2F` are computed first and the composite jvp is
rebuilt at `u0` before the bordered solves.
"""
function get_curvature(prob::ManifoldProblemBKMatrixFree, u0::AbstractVector{T}, Φ, par, weights) where {T}
    @error "Curvature - BK"
    n, m = size(prob)
    d = n - m
    𝒯 = eltype(u0)

    bls = prob.get_tangent isa BLSBorderedTangent ? prob.get_tangent.bls :
        error("get_curvature requires `get_tangent = BLSBorderedTangent(bls)` for a matrix-free problem.")

    # 1) composite second derivatives on the upper triangle (the calls advance the stateful jvp)
    d2s = Matrix{Vector{𝒯}}(undef, d, d)
    for i in Base.OneTo(d), j in i:d
        d2s[i, j] = d2F(prob, u0, par, @view(Φ[:, i]), @view(Φ[:, j]))
    end

    # 2) rebuild the composite jvp at u0 for the bordered solves
    J = jacobian(prob, u0, par)
    Jmat = J isa AbstractMatrix
    J1 = Jmat ? view(J, :, 1:m) : (v -> BK.apply(J, vcat(v, zeros(𝒯, d))))

    pb = prob.VF.VF.F
    dₚF = if pb isa BifurcationProblem_2P
        _jacobian_param(pb, u0, par)
    else
        unit(i) = (v = zeros(𝒯, n); v[i] = one(𝒯); v)
        ntuple(i -> BK.apply(J, unit(m + i)), d)
    end

    # border Φ' D : blocs B1 = Φ[:, 1:m]' D and B2 = Φ[:, m+1:n]' D
    is_trivial = weights.weights isa TrivialWeight
    Φ = is_trivial ? Φ : LinearAlgebra.Diagonal(weights.weights) * Φ
    b = ntuple(i -> @view(Φ[1:m, i]), d)
    cb = Φ[(m+1):n, :]'

    # 3) cmat[i,j] = ‖ [J;Φ'D]⁻¹ [d²F(Φi,Φj); 0] ‖_D
    cmat = zeros(𝒯, d, d)
    for i in Base.OneTo(d), j in i:d
        rhs = zeros(𝒯, d)
        u1, u2, _, _ = BK.solve_bls_block(bls, J1, dₚF, b, cb, d2s[i, j], rhs)
        cmat[i, j] = weighted_norm(weights, vcat(u1, u2))
        cmat[j, i] = cmat[i, j]
    end
    K = LinearAlgebra.eigmax(Symmetric(cmat))
    return K
end
#━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━