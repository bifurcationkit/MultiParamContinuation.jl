# https://jacco.ompf2.com/2022/04/13/how-to-build-a-bvh-part-1-basics/
# https://github.com/alecjacobson/computer-graphics-bounding-volume-hierarchy

# for walking through the tree
# https://docs.julialang.org/en/v1/base/numbers/#Base.digits
using StaticArrays, Statistics

struct AABB{T}
    min::T
    max::T
end

Base.show(io::IO, box::AABB) = println(io, "AABB\n  ├─ lower = ", box.min, "\n  └─ upper = ", box.max)
inside(aabb::AABB, point::AbstractVector) =  all(aabb.min[i] <= point[i] < aabb.max[i] for i in eachindex(point))
Base.in(A::AABB, B::AABB) = inside(B, A.min) && inside(B, A.max)
# intersect(A::AABB, B::AABB) = inside(A, B.min) || inside(A, B.max) # not used

"""
$SIGNATURES

Overlaps of two AABB. If one is inside the other one, it also returns true.

## Note

Touching equals overlap.
"""
function overlaps(A::AABB, B::AABB)
    result = true
    for i in eachindex(A.min)
        if (B.max[i] < A.min[i] || B.min[i] > A.max[i])
            result = false
            break
        end
    end
    return result
end

function update_bb!(A::AABB, B::AABB)
    for i in eachindex(A.min)
        A.min[i] = min(A.min[i], B.min[i])
        A.max[i] = max(A.max[i], B.max[i])
    end
    A
end

function update_bb!(A::AABB, c::Chart)
    for i in eachindex(A.min)
        A.min[i] = min(A.min[i], c.u[i] - c.R)
        A.max[i] = max(A.max[i], c.u[i] + c.R)
    end
    A
end

"""
Point of chart `S[id]` used to build the bounding boxes of the tree. It is the
projection `project_for_tree(u, p)` when provided, otherwise the raw point `u`. The
projection must be chosen so that the sphere test `dist2(weights, u1, u2) <= (R1 + R2)^2`
implies `max_i |p1[i] - p2[i]| <= R1 + R2`.
"""
function _tree_point(S::Atlas, id::Int)
    prob = S.alg.prob
    pr = prob.project_for_tree(S[id].u, prob.params)
    return isnothing(pr) ? S[id].u : pr
end

function update_bb!(A::AABB, S::Atlas, id::Int)
    c = S[id]
    p = _tree_point(S, id)
    for i in eachindex(A.min)
        A.min[i] = min(A.min[i], p[i] - c.R)
        A.max[i] = max(A.max[i], p[i] + c.R)
    end
    A
end

update_bb!(::AABB, ::Nothing) = nothing
update_bb!(::Nothing, ::AABB) = nothing

"""
Implement a Bounding Volume Hierarchy
"""
mutable struct BVHNode{T, Tids, Tsplit}
    box::AABB{T}
    parent::Union{Nothing, BVHNode{T, Tids, Tsplit}}
    left_child::Union{Nothing, BVHNode{T, Tids, Tsplit}}
    right_child::Union{Nothing, BVHNode{T, Tids, Tsplit}}
    is_leaf::Bool
    chart_ids::Tids
    max_size::Int64
    split_dim::Int64
    split_value::Tsplit
end
is_leaf(node::BVHNode) = node.is_leaf
is_root(node::BVHNode) = isnothing(node.parent)
Base.in(A::AABB, node::BVHNode) = Base.in(A, node.box)

overlaps(node::BVHNode, A::AABB) = overlaps(node.box, A)

function parent_split_dim(node)
    if isnothing(node.parent)
        return 0
    end
    return node.parent.split_dim
end

function Base.show(io::IO, node::BVHNode)
    if is_leaf(node)
        printstyled(io, "LEAF\n", bold=true, color = :blue)
    else
        printstyled(io, "NODE\n", bold=true, color = :blue)
    end
    println(io, "├─ #     = ", npoints(node), " points")
    println(io, "├─ lower = ", node.box.min)
    println(io, "├─ upper = ", node.box.max)
    if is_leaf(node)
        println(io, "└─ ids   = ", node.chart_ids)
    end
end

update_parents_aabb!(::Nothing) = nothing

function update_parents_aabb!(node::BVHNode)
    if ~isnothing(node.right_child)
        update_bb!(node.box, node.right_child.box)
    end
    if ~isnothing(node.left_child)
        update_bb!(node.box, node.left_child.box)
    end
    if ~isnothing(node.parent)
        update_parents_aabb!(node.parent)
    end
    return
end

function npoints(n::BVHNode) 
    if n.is_leaf
        return length(n.chart_ids)
    end
    return npoints(n.left_child) + npoints(n.right_child)
end

function BVHNode(dim = 3; 
                max_size::Int = 5,
                parent = nothing,
                left_child = nothing,
                right_child = nothing,
                is_leaf = true,
                chart_ids = Int[],
                T = Float64,
                split_value = Inf * one(T),
                split_dim = 0,
                )
    # sa_bd = MVector{dim, T}(1:dim...)
    sa_bd = ones(T, dim)

    BVHNode(
            AABB(
                sa_bd * Inf,
                sa_bd * (-Inf),
            ),
            parent,
            left_child,
            right_child,
            is_leaf,
            chart_ids,
            max_size,
            split_dim,
            split_value,
    )
end

function add!(node::BVHNode, S::Atlas, id::Int)
    if isnothing(node.left_child) &&
        isnothing(node.right_child)
        if npoints(node) < node.max_size && is_leaf(node)
            # node is leaf with enough space, add id to chart_ids
            push!(node.chart_ids, id)
            # update the bounding box
            update_bb!(node.box, S, id)
        elseif is_leaf(node)
            # split the node in the direction with largest spread
            node_ids = SA[node.chart_ids..., id]
            maxs = node.box.max
            mins = node.box.min
            perm = sortperm(maxs - mins, rev = true)

            split_dim = perm[1]
            # split_value = (maxs[split_dim] + mins[split_dim])/2
            split_value = median(_tree_point(S, i)[split_dim] for i in node_ids)

            node.is_leaf = false
            update_bb!(node.box, S, id)

            # update split info in current node
            node.split_dim = parent_split_dim(node) == split_dim ? perm[2] : split_dim
            node.split_value = split_value

            node.left_child = BVHNode(length(node.box.min);parent = node, split_dim = 0, max_size = node.max_size)
            node.right_child = BVHNode(length(node.box.min);parent = node, split_dim = 0, max_size = node.max_size)

            # partition of the ids. The median split may leave one side empty
            # (e.g. all points share the same coordinate along `split_dim`),
            # in which case we fall back to a half/half partition to guarantee
            # that the recursion terminates.
            ids_left  = [i for i in node_ids if _tree_point(S, i)[split_dim] <  split_value]
            ids_right = [i for i in node_ids if _tree_point(S, i)[split_dim] >= split_value]
            if isempty(ids_left) || isempty(ids_right)
                mid = cld(length(node_ids), 2)
                ids_left  = collect(node_ids[1:mid])
                ids_right = collect(node_ids[mid+1:end])
            end

            for i in ids_left
                add!(node.left_child, S, i)
            end
            for i in ids_right
                add!(node.right_child, S, i)
            end

            empty!(node.chart_ids)
        end
    else
         # we need to add the chart to the appropriate sub node in the tree
        if _tree_point(S, id)[node.split_dim] < node.split_value
            add!(node.left_child, S, id)
        else
            add!(node.right_child, S, id)
        end
    end
    update_parents_aabb!(node)
    node
end

function BVHNode(S::Atlas, dim = 3; max_size::Int = 5)
    tree = BVHNode(dim; max_size)
    for id in eachindex(S.atlas)
        add!(tree, S, id)
    end
    tree
end

"""
$SIGNATURES

Compute the nearest neighbors from `S[id]` in the BVH tree. We know that AABB(S[id]) is in the tree so it simplifies the search.
"""
function neighbors(tree::BVHNode, S::Atlas, id::Int)
    c = S[id]
    p = _tree_point(S, id)
    aabb = AABB(p .- c.R, p .+ c.R)
    list = Int[]
    _query(tree, aabb, list)
    list
end

function _query(tree::BVHNode, aabb::AABB, list)
    if overlaps(tree, aabb)
        if is_leaf(tree)
            append!(list, tree.chart_ids)
        end
        if ~isnothing(tree.right_child)
            _query(tree.right_child, aabb, list)
        end
        if ~isnothing(tree.left_child)
            _query(tree.left_child, aabb, list)
        end
    end
end