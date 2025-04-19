using Revise, Test, ForwardDiff, GLMakie
using BifurcationKit
const BK = BifurcationKit

using Logging: global_logger
using TerminalLoggers: TerminalLogger
global_logger(TerminalLogger())

Makie.inline!(true)
# Makie.inline!(false)
####################################################################################################
function abc!(dz, z, p, t = 0)
    (;D, B, σ, β, α) = p
    u1, u2, u3 = z
    dz[1] = -u1 + D*(1 - u1)*exp(u3)
    dz[2] = -u2 + D*(1 - u1)*exp(u3) - D*σ*u2*exp(u3)
    dz[3] = -u3 - β*u3 + D*B*(1 - u1)*exp(u3) + D*B*α*σ*u2*exp(u3)
    dz
end

# we group the differentials together
par_abc = (D = 0.11, B = 8., α = 1., σ = 0.04, β = 1.56)
z0 = [1., 0., 0. ]
prob_bk = BifurcationProblem(abc!, z0, par_abc, (@optic _.D), 
        record_from_solution = (x, p; k...) -> (u3 = x[3], u1 = x[1], u2 = x[2]),)

opts_br = ContinuationPar(p_max = .5, n_inversion = 8, nev = 3)
br = BK.continuation(prob_bk, PALC(), opts_br; normC = norminf)

BK.plot(br, plotfold=false, )
####################################################################################################
argspo = (record_from_solution = (x, p; k...) -> begin
                xtt = BK.get_periodic_orbit(p.prob, x, p.p)
                return (max = maximum(xtt[3,:]), min = minimum(xtt[3,:]), period = x[end])
            end,
    plot_solution = (ax, x, p; k...) -> begin
        xtt = BK.get_periodic_orbit(p.prob, x, p.p)
        lines!(ax, xtt.t, xtt.u[1,:]; label = "u1", linewidth = 2)
        lines!(ax, xtt.t, xtt.u[2,:]; label = "u2")
        BK.plot!(get(k, :ax1, nothing), br)
    end,)
####################################################################################################
hopfpt = get_normal_form(br, 1)

# continuation parameters
opts_po_cont = ContinuationPar(dsmax = 0.03, dsmin = 1e-4, ds = 0.0005, max_steps = 130, tol_stability = 4e-2, plot_every_step = 20)
@reset opts_po_cont.newton_options.tol = 1e-11

br_po = BK.continuation(
    br, 1, opts_po_cont,
    PeriodicOrbitOCollProblem(50, 4; update_section_every_step = 1, jacobian = BK.DenseAnalyticalInplace());
    δp = 0.0001,
    linear_algo = BK.COPBLS(),
    verbosity = 1,
    plot = true,
    argspo...,
    normC = norminf)

BK.plot(br_po, br)
####################################################################################################
function plot_data(S; k...)
    fig = Figure()
    ax = Axis3(fig[1,1], zlabel = "u3", xlabel = "D", ylabel = "β", title = "$(length(S)) charts")
    plot_data!(ax, S; k...)
    fig
end

function plot_data!(ax, S; ind_col = 1, fil=x->true, cols = [real(c.index) for c in filter(x -> fil(x.u), S.atlas)])
    pts = mapreduce(c->[c.data[2], c.data[3], c.data[1]]', vcat, filter(x -> fil(x.u), S.atlas))
    cols = [real(c.data[ind_col]) for c in filter(x -> fil(x.u), S.atlas)]
    hm = scatter!(ax, pts, color = cols)
    # cb = Colorbar(fig, hm)
    # fig[1,2] = cb
end

function build_mesh!(ax, S; fil=x->true,ind_col = 1)
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
    cols = [real(c.data[ind_col]) for c in filter(x -> fil(x.u), S.atlas)]

    poly!(ax, centers', reduce(vcat, faces), strokewidth=1, color = :blue)
    # scatter!(ax, centers)
    ax
end

function build_mesh(S; k...)
    f = Figure()
    ax = Axis3(f[1,1], aspect = :equal, zlabel = "u3", xlabel = "D", ylabel = "β", title = "$(length(S)) charts")
    build_mesh!(ax, S; k...)
    f
end
####################################################################################################
using MultiParamContinuation
const MPC = MultiParamContinuation

prob = MPC.ManifoldProblem_BK(
                        prob_bk, br.sol[1].x, (@optic _.D), (@optic _.β);
                        record_from_solution = (X, p; k...) -> begin
                            return (β = X[end], D = X[end-1], u3 = X[3])
                        end,
                        finalize_solution = (X,p) -> begin
                            D = X[end-1]
                            β = X[end]
                            keep = (0.01 <= D <= 0.5) && (1.5 <= β <= 1.65)
                            return keep
                        end,
                        )

S_eq = @time MPC.continuation(prob,
                        Henderson(np0 = 3,
                                    θmin = 0.001,
                                    # use_curvature = true,
                                    use_tree = true,
                                  ),
                        CoveringPar(max_charts = 20000,
                                max_steps = 1000,
                                verbose = 0,
                                newton_options = NewtonPar(tol = 1e-10, verbose = false),
                                R0 = .04,
                                ϵ = 0.1,
                                delta_angle = 10.15,
                                ))

MPC.plot2d(S_eq,ind_plot=4:5)

plot_data(S_eq)
build_mesh(S_eq)
