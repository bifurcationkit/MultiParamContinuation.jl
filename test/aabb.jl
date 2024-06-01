using MultiParamContinuation
const MPC = MultiParamContinuation

aabb1 = MPC.AABB([0,0], [1,1])
aabb2 = MPC.AABB([1.1,1.2], [1.2,1.5])
show(aabb1)

function _plot!(ax, box::MPC.AABB)
    side = box.max - box.min
    poly!(ax, Rect(box.min[1], box.min[2], side[1], side[2]); 
        color = :white,
        strokewidth = 2,
        alpha = 0.1,
            )
    ax
end

# f = Figure(); ax=Axis(f[1,1]);_plot!(ax, aabb1);_plot!(ax, aabb2);f
@test MPC.intersect(aabb1, aabb2) == false
@test MPC.overlaps(aabb1, aabb2) == false

tree = BVHNode(2)
@test MPC.is_leaf(tree) == true
@test MPC.is_root(tree) == true

MPC.update_bb!(tree.box, aabb1)
@test MPC.overlaps(tree, aabb2) == false