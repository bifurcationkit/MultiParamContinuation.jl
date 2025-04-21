using LinearAlgebra

"""
$TYPEDEF

Structure to define a chart, that is a local system of coordinates which parametrizes a small patch of the manifold.

## Fields

$TYPEDFIELDS

## Methods

- `do_intersect(c1, c2)::Bool` whether the charts intersect

"""
mutable struct Chart{Tu, Ttg, Tr, Tp, Tin, Td, Teve, Tl}
    "Base point, solution of F(x) = 0"
    const u::Tu
    "Tangent space at u"
    const Φ::Ttg
    "Radius for region of validity (Ball)."
    R::Tr # modified by curvature
    "Polyhedron. Represents the domain of the solution manifold covered by this chart."
    P::Tp
    "Array of booleans to check if the vertices in P are inside the ball of validity"
    inside_ball::Tin
    "Is the polygon inside the validity ball?"
    interior::Bool
    "Chart index"
    const index::Int
    "Data, like from record_from_solution, or eigenvalues, etc."
    const data::Td
    "[Internal] Event values"
    const event_values::Teve
    "[Internal] Label, like for event detection"
    const label::Tl
    "[Internal] list of direct neighbors in the atlas" 
    neighbors::Vector{Int}
end

@inline is_inside_ball(c::Chart, P) = norm(P, 2) <= c.R

function is_on_boundary(c::Chart)
    c.interior = ~all(c.inside_ball)
    return c.interior
end

function new_chart(u0, Φ, Radius, P; 
                            id = 0,
                            eve = nothing,
                            data = nothing,
                            label = Symbol(),
                            neighbors = Int[])
    Chart(u0, Φ, Radius, P, 
            ([false for _ in eachindex(P)]), 
            false,
            id,
            data,
            eve,
            label,
            neighbors
            )
end

function update!(c::Chart)
    c.interior = is_on_boundary(c)
end

function init_polygonal_boundary(N, R) 
    v = R / cospi(1/N)
    return [@SVector [cospi(2*(i-1)/N) * v,
                      sinpi(2*(i-1)/N) * v] for i in Base.OneTo(N)]
end

function get_alpha(c1::Chart, c2::Chart)
    R1 = c1.R
    R2 = c2.R
    α₁₂ = (1 + (R1 - R2) * (R1 + R2)) / 2
end

"""
$TYPEDEF

Atlas of charts which represents a manifold.

## Fields

$TYPEDFIELDS

## Methods

- `add!(a::Atlas, c::Chart)`

- `new_atlas(c::Chart, alg = nothing; dim = 2)`

- `length(a::Atlas)` returns the number of charts

- `a[3]` returns the 3rd chart in the atlas `a`, see `?Chart`

"""
struct Atlas{dim, Tc, Talg, Ttree}
    "List of Charts"
    atlas::Vector{Tc}
    "[Internal] Boundary list of charts"
    BList::Vector{Tc} # Cf Henderson - 2002
    "Algorithm, for example: `Henderson()`"
    alg::Talg
    "[Internal] Tree for neighbors search"
    tree::Ttree
end
# see https://discourse.julialang.org/t/mutable-struct-vs-ref-in-an-immutable-one/92846/2
Base.length(Ω::Atlas) = length(Ω.atlas)
Base.getindex(Ω::Atlas, k::Int) = getindex(Ω.atlas, k)
Base.lastindex(Ω::Atlas) = lastindex(Ω.atlas)
@inline get_boundary_list(Ω::Atlas) = Ω.BList
@inline use_tree(Ω::Atlas) = Ω.alg.alg.use_tree

# constructor
function new_atlas(c::Tc, cache::Talg = nothing; dim::Int = 2) where {Tu, Ttg, Tc <: Chart{Tu, Ttg}, Talg}
    max_size = cache.alg.children_pre_leaf
    tree = cache.alg.use_tree ? BVHNode(length(c.u); max_size) : nothing
    Ω = Atlas{dim, Tc, Talg, typeof(tree)}([c], Vector{Tc}(), cache, tree)
    if use_tree(Ω)
        add!(Ω.tree, Ω, length(Ω))
    end
    return Ω
end

# add chart to atlas
function add!(Ω::Atlas, c::Chart)
    push!(Ω.atlas, c)
    if use_tree(Ω)
        @assert c.index == length(Ω)
        add!(Ω.tree, Ω, length(Ω))
    end
end

function check_alphas(Ω::Atlas, cᵢ::Chart)
    blist = intersec_list(Ω, cᵢ)
    for j in blist
        cΩ = Ω[j]
        @assert j == cΩ.index
        αᵢⱼ = get_alpha(cᵢ, cΩ)
        if αᵢⱼ < 0 
            return false
        end
    end
    return true
end

function update_boundary!(Ω::Atlas)
    prob = Ω.alg.prob
    BList = get_boundary_list(Ω) 
    empty!(BList)
    for c in Ω.atlas
        update!(c)
        if is_on_boundary(c) && prob.finalize_solution(c.u, prob.params)
            push!(BList, c)
        end
    end
end

"""
This is an over-estimate on the charts that intersect c. Better look at `true_intersec_list`.
"""
function intersec_list(Ω::Atlas, c::Chart, use_tree_bool::Bool = use_tree(Ω))
    Jᵢᵐ = Int[]
    if ~use_tree_bool
        for cΩ in Ω.atlas
            if do_intersect(cΩ, c)
                push!(Jᵢᵐ, cΩ.index)
            end
        end
    else
        @assert c.index <= length(Ω)
        list = neighbors(Ω.tree, Ω, c.index)
        for id in list
            if do_intersect(Ω[id], c)
                push!(Jᵢᵐ, id)
            end
        end
    end
    return Jᵢᵐ
end

"""

Function mainly for plotting.
"""
function test_P(charti::Chart, chartj::Chart)
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
    return n₊ == n || n₊ == 0
end

# function mainly for plotting
function true_intersec_list(Ω::Atlas, c::Chart, use_tree_bool::Bool = use_tree(Ω))
    inter_list = intersec_list(Ω, c, use_tree_bool)
    out = Int[]
    for id in inter_list
        if ~test_P(c, Ω[id])
            push!(out, id)
        end
    end
    out
end

function do_intersect(c1::Chart, c2::Chart)
    # !! TODO MAKE IT WORK ON GPU
    u1 = c1.u
    R1 = c1.R
    u2 = c2.u
    R2 = c2.R
    return dist2(u1, u2) <= (R1 + R2)^2
end

function Base.show(io::IO, Ω::Atlas{dim}) where {dim}
    println(io, "Surface $dim-d")
    println(io, "   ├─ # charts = ", length(Ω))
    println(io, "   └─ problem = ")
    show(io, Ω.alg.prob; prefix = " "^10)
end

function Base.show(io::IO, c::Chart)
    println(io, "Chart")
    println(io, "  ├─ id = ", c.index)
    println(io, "  ├─ R = ", c.R)
    println(io, "  ├─ on boundary = ", c.interior)
    println(io, "  ├─ u = ", c.u)
    println(io, "  ├─ inside validity ball = ", c.inside_ball)
    if ~isnothing(c.data)
        println(io, "  ├─ data = ", c.data)
    end
    if ~isnothing(c.event_values)
        println(io, "  ├─ event values = ", c.event_values)
    end
    if c.label != Symbol()
        println(io, "  ├─ event label = ", c.label)
    end
    println(io, "  └─ Polygon = ")
    display(c.P)
end