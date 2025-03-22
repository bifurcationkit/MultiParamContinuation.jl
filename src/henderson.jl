abstract type AbstractCoveringAlgorithm end

"""
$TYPEDEF

Continuation algorithm from [1] computing implicitly defined 2d manifolds. 

## Fields

$TYPEDFIELDS

## Reference(s)

[1] Henderson, Michael E. “Multiple Parameter Continuation: Computing Implicitly Defined k-Manifolds.” International Journal of Bifurcation and Chaos 12, no. 03 (March 2002): 451-76. https://doi.org/10.1142/S0218127402004498.

"""
Base.@kwdef struct Henderson{T} <: AbstractCoveringAlgorithm
    "Number of initial vertices in the polyhedra on each tangent space."
    np0::Int = 4
    "Use hessian for curvature estimate."
    use_curvature::Bool = false
    "[Internal]"
    θmax::T = 1.2
    "[Internal]"
    θmin::T = 0.01
    " Use tree to find neighbors. Useful when the number of charts is large because the complexity changes from N² to N⋅log(N)."
    use_tree::Bool = false
    "Number of children per leaf in the tree. Control the depth of the tree."
    children_pre_leaf::Int = 5
end

"""
$TYPEDEF

Cache for `::Henderson` algorithm. 

## Fields

$TYPEDFIELDS
"""
mutable struct HendersonCache{T1, T2 <: CoveringPar, T3, T4, T5}
    "A manifold problem"
    const prob::T1
    "Continuation parameters"
    const contparams::T2
    "Continuation algorithm"
    const alg::T3
    "[Internal] Default value = 1"
    θ::T4
    "[Internal] cache for storing the tangent space"
    _rhs_tangent::T5
end

"""
$SIGNATURES

Main function for covering the manifold defined in `prob` using the algorithm `alg`.
"""
function continuation(prob::AbstractManifoldProblem, 
                      alg::Henderson, 
                      contparams::CoveringPar;
                      verbose = contparams.verbose,
                      θ = one(eltype(prob)))
    _verbose = verbose > 0
    # create cache for covering algorithm
    n, m = size(prob)
    dim = n - m
    cache = HendersonCache(prob, contparams, alg, θ, vcat(zeros(m, dim), I(dim)))
    chart0 = init(cache)
    n, m = size(prob)
    Ω = new_atlas(chart0, cache; dim = n - m )
    update_boundary!(Ω)
    n_steps = 1
    while length(Ω) < contparams.max_charts && 
            n_steps < contparams.max_steps
        if _verbose
            println("━"^50)
            println("─── step     = ", n_steps)
            println(" ├─ # charts = ", length(Ω))
            println(" └─ new boundary chart R = ", Ω[end].R)
            println("           ├─         u[1:3] = ", Ω[end].u[1:3])
            if ~isnothing(Ω[end].data)
                println("           ├─         data   = ", Ω[end].data)
            end
            println("           └─         id = ", Ω[end].index)
        end
        if ~step!(Ω)
            return Ω
        end
        n_steps += 1
    end
    return Ω
end 

"""
$SIGNATURES

Perform one step of the continuation algorithm.
"""
function step!(Ω::Atlas)
    alg = Ω.alg
    new_chart = generate_new_chart(Ω)
    if isnothing(new_chart)
        return false
    end
    # We first add new_chart to Ω in order to update the tree.
    # We then update the polygon set of new_chart in remove_halfspace!
    add!(Ω, new_chart)
    remove_halfspace!(Ω, new_chart)
    if ~is_on_boundary(new_chart) && alg.contparams.verbose > 0
        @warn "new chart $(new_chart.index) is not on boundary"
    end
    update_boundary!(Ω)
    return true
end

"""
$SIGNATURES

Perform n steps of the continuation algorithm.
"""
function step!(Ω::Atlas, n::Int)
    @progress for _ in Base.OneTo(n)
        step!(Ω)
    end
    Ω
end

function init(cache::HendersonCache)
    (; prob, alg, contparams) = cache
    # get first point on manifold
    if _has_projection(prob)
        u0 = project(prob, prob.u0, prob.params)
    else
        u0 = correct_guess(cache, contparams.newton_options)
    end
    # get tangent space
    T = get_tangent(prob, u0, prob.params, cache._rhs_tangent)
    R = contparams.R0
    Ps = init_polygonal_boundary(alg.np0, R)
    data = prob.recordFromSolution(u0, prob.params)
    eve = prob.event_function(u0, prob.params)
    return new_chart(u0, T, R, Ps; id = 1, data, eve)
end

function correct_guess(cache, nl_spec::NonLinearSolveSpec)
    prob = cache.prob
    probnl = NonlinearProblem(prob.VF, prob.u0, prob.params)
    sol = solve(probnl, nl_spec.nl_solver; nl_spec.options...)
    if sol.retcode != ReturnCode.Success
        throw("Newton for first point did not converge!!")
    end
    return sol.u
end

"""
$SIGNATURES

Compute the projection of guess on the manifold M
"""
function project_on_M(prob, guess, chart::Chart, wbar, cpar::CoveringPar{T, <: NonLinearSolveSpec}) where {T}
    if _has_projection(prob)
        return project(prob, guess, prob.params)
    else
        nl_spec = cpar.newton_options
        Φ = chart.Φ
        function f(w,p)
            vcat(prob.VF(w, p), Φ' * (w - wbar))
        end
        prob_bls = NonlinearProblem(f, guess, prob.params)
        sol = solve(prob_bls, nl_spec.nl_solver; nl_spec.options... )
    end

    if sol.retcode == ReturnCode.Success
        return sol.u
    else
        return nothing
    end
end

"""
Create a chart from guess after projecting it on the manifold.
"""
function _new_chart_from_guess(cache, chart, ω; 
                    R = cache.contparams.R0, 
                    id = 0)
    (;prob, contparams) = cache
    guess = chart.u .+ chart.Φ * ω
    u = project_on_M(prob, guess, chart, copy(guess), contparams)
    if isnothing(u)
        return nothing
    end
    Φ = get_tangent(prob, u, prob.params, cache._rhs_tangent)
    data = prob.recordFromSolution(u, prob.params)
    eve = prob.event_function(u, prob.params)
    if isnothing(eve)
        label = Symbol()
    else
        label = eve * chart.event_values < 0 ? :EVE : Symbol()
    end
    return new_chart(u, 
                Φ, 
                R, 
                (init_polygonal_boundary(
                        cache.alg.np0, 
                        R * 1)); 
                id, data, eve, label)
end

function get_interior_vertex_in_T(c::Chart, sₑ, ds = 1)
    uᵢ = c.u
    # get vertex inside ball
    # this vertex is null here
    if false #keep this please!
        sᵢ = c.Φ' * (uᵢ - c.u)
        # project in tangent space
        s = sₑ - sᵢ
    end
    # s = copy(sₑ)
    # get vertex at the intersection between the ray s and the Ball B(0, R)
    s = sₑ .* (c.R / norm(sₑ) * ds)
    return s
end

function generate_exterior_vertex(Ω::Atlas, clist::Vector{<: Chart})
    for chart in clist
        if is_on_boundary(chart)
            for (ind, P) in pairs(chart.P)
                if chart.inside_ball[ind] == false
                    if true#check_alphas(Ω, chart) # Henderson p. 463
                        return chart, P
                    end
                end
            end
        end
    end
    return nothing
end

function generate_new_chart(Ω::Atlas; id = length(Ω) + 1)
    Blist = get_boundary_list(Ω)
    guess = generate_exterior_vertex(Ω, Blist)
    if isnothing(guess)
        return nothing
    end
    c, sₑ = guess
    cache = Ω.alg
    contparams = cache.contparams
    verbose = contparams.verbose > 1
    (;ϵ, delta_angle) = contparams
    t = cache.θ
    (;θmin, θmax) = cache.alg
    iter = 1
    s = get_interior_vertex_in_T(c, sₑ)
    while t > θmin
        verbose && println("----> iteration = ", iter, ", t = $t, from chart = ", c.index, "\n u = ", c.u)
        ω = t .* s
        new_chart = _new_chart_from_guess(cache, c, ω; R = c.R, id)
        if isnothing(new_chart)
            @goto failed
        end
        if cache.alg.use_curvature
            K = get_curvature(cache.prob, new_chart, cache.prob.params)
            radius_estimate = sqrt(2ϵ / K)
            verbose && @error "Radius est" radius_estimate c.R K
            new_chart.R = min(c.R, radius_estimate)
        end
        # distance from guess to projected point
        dst = norm(new_chart.u .- (c.u .+ c.Φ * ω), Inf)
        # angle between tangent spaces, do not compute if delta_angle large enough
        δα = delta_angle > pi ? 0 : abs(largest_principal_angle(c.Φ, new_chart.Φ))
        if δα > delta_angle ||
                    dst > ϵ
            verbose && @error "Reduction"  δα dst cache.θ
            @goto failed
        else 
            if dst < ϵ/2 && iter == 1
                cache.θ = min(cache.θ * 1.2, θmax)
            end
            return new_chart
        end
        # alg.θ = min(1.1 * alg.θ, 1.)
        # c.R *= 1.005

        @label failed
        t *= 0.8
        iter +=1
    end
    return nothing
end

"""
$SIGNATURES

Change `charti.P` to satisfy halfspace condition. Does not mutate the second argument.

## Reference

[2] Henry, Damennick B., and Daniel J. Scheeres. “Fully Numerical Computation of Heteroclinic Connection Families in the Spatial Three-Body Problem.” Communications in Nonlinear Science and Numerical Simulation 130 (March 2024): 107780. https://doi.org/10.1016/j.cnsns.2023.107780.
"""
function remove_halfspace_first_chart!(charti::Chart, chartj::Chart)
    ui = charti.u
    Ri = charti.R
    Ti = charti.Φ

    uj = chartj.u
    Rj = chartj.R

    du = Ti' * (uj .- ui)
    Bound = Ri^2 - Rj^2 + norm(du, 2)^2
    
    testp = [2dot(s, du) < Bound for s in charti.P]
    n = length(testp)
    n₊ = sum(testp)
    n₋ = n - n₊
    # @assert n₊ > 0
    if n₊ == n || n₊ == 0
        return false
    end

    ind = 0
    while sum(testp[1:n₊]) < n₊ && ind < n
        α = testp[1]
        for ii=1:n-1
            testp[ii] = testp[ii+1]
        end
        testp[n] = α
        ind +=1
    end

    @assert sum(testp) == n₊

    _l = length(charti.P)
    rg = eachindex(charti.P) .+ ind
    Pi_circ = [charti.P[mod1(ii, _l)] for ii in rg]
    # the testp looks like [1,1,1,1,1,0,0,0,0]
    testp = ([2dot(s, du) < Bound for s in Pi_circ])

    j = findlast(testp) # exist because n₊ > 0
    newPi = [Pi_circ[k] for k = 1:j]

    # we update Pj
    dP = Pi_circ[j+1] .- Pi_circ[j]
    t = (Bound / 2 - dot(Pi_circ[j], du)) / (dot(dP, du))
    new_P_indj = @. Pi_circ[j] + dP * t
    push!(newPi, new_P_indj)

    # we update Pj+1
    j = length(Pi_circ)
    dP = Pi_circ[mod1(j+1,_l)] .- Pi_circ[mod1(j,_l)]
    t = (Bound / 2 - dot(Pi_circ[j], du)) / (dot(dP, du))
    new_P_indj = Pi_circ[j] .+ dP .* t
    push!(newPi, new_P_indj)
    charti.P = newPi

    charti.inside_ball = ([is_inside_ball(charti, P) for P in charti.P])
    update!(charti)
    return true
end

function remove_halfspace!(Ω::Atlas, c1::Chart)
    verbose = Ω.alg.contparams.verbose > 1
    int_list = intersec_list(Ω, c1)
    for id in int_list
        cΩ = Ω[id]
        @assert cΩ.index == id
        verbose && println("Merge c$(c1.index) with $(cΩ.index)")
        remove_halfspace_first_chart!(c1, cΩ)
    end
    for id in int_list
        cΩ = Ω[id]
        @assert cΩ.index == id
        verbose && println("Merge c$(cΩ.index) with $(c1.index)")
        remove_halfspace_first_chart!(cΩ, c1)
    end
end

### TODO REMOVE THIS FUNCTION
function _remove_halfspace!(Ω::Atlas, c1::Chart)
    verbose = Ω.alg.contparams.verbose > 0
    int_list = intersec_list(Ω, c1)
    for index in int_list
        cΩ = Ω[id]
        @assert cΩ.index == id
        verbose && println("Merge c$(c1.index) with $(cΩ.index)")
        remove_halfspace_first_chart!(c1, cΩ)
    end
    c1
end