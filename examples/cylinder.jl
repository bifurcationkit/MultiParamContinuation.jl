using Revise
cd(@__DIR__)
using Pkg
pkg"activate ."

using GLMakie
Makie.inline!(false)
Makie.inline!(true)

using MultiParamContinuation

const MPC = MultiParamContinuation

F(u,p) = [u[2]^2 + u[3]^2 - 1]

prob = ManifoldProblem(F, [0,0,1.], nothing;
            finalize_solution = (u,p) -> -2<=u[1]<=2)

contpar = CoveringPar(max_charts = 1500, 
                                max_steps = 320,
                                verbose = 0,
                                newton_options = NonLinearSolveSpec(;maxiters = 5, abstol = 1e-12, reltol = 1e-10),
                                R0 = .25,
                                # ϵ = 0.1,
                                # delta_angle = 0.1,
                                ); 
alg = Henderson(np0 = 4,
                θmin = 0.001,
                # use_curvature = true
                )

S = continuation(prob,
            alg,
            contpar)

MPC.plotd(S; 
    # draw_circle = true, 
    draw_tangent = true,
    draw_edges = true,
    plot_center = true,
    # put_ids = true,
    ind_plot = [1,3]
    )


MPC.plot2d(S; 
    draw_circle = true, 
    draw_tangent = true,
    draw_edges = true,
    plot_center = true,
    put_ids = true,
    # ind_plot = [1,3]
    )

MPC.intersec_list(S, S[1])

function build_mesh(S)
    neighbors_list = Vector{Int}[]
    for c in S.atlas
        push!(neighbors_list, MPC.true_intersec_list(S, c))
    end

    # build faces 
    faces = Vector{Int}[]
    for c in S.atlas
        neighbors = neighbors_list[c.index]
        for n1 in neighbors
            for n2 in neighbors
                if n2 != n1 && n2 in neighbors_list[n1]
                    push!(faces, [c.index, n1, n2])
                end
            end
        end
    end
    centers = mapreduce(c->[c.data[2], c.data[3], c.data[1]]', vcat, S.atlas)
    @error "" size(centers)
    f = Figure()
    ax = Axis3(f[1,1], aspect = :equal, zlabel = "u3", xlabel = "D", ylabel = "β", title = "$(length(S)) charts")
    poly!(ax, centers', reduce(vcat, faces), strokewidth=1, shading=true)
    scatter!(ax, centers)
    f
end

build_mesh(S)