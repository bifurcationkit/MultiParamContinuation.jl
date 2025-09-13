using Makie.Colors
using Makie.GeometryBasics
using StaticArrays

@inline get_tangent(A::AbstractMatrix, i::Int) = @view A[:, i]

"""
$SIGNATURES

3d plot of the manifold.

## Optional arguments
- `draw_circle::Bool` plot the validity ball
- `draw_tangent = true` plot the tangent space
- `plot_center = true`, plot the center of the charts
- `put_ids = false` add the indices of the charts
- `draw_frame = false` plot the local frame on the tangent space
- `ind_plot = (1,2,3)` indices from `u` to plot
"""
function plotd(ax, chart::Chart; 
                ind = 1,
                angles = nothing,
                draw_circle = false,
                colorrange = nothing,
                draw_frame = false,
                draw_tangent = true,
                put_ids = false,
                plot_center = false,
                draw_edges = false,
                record_from_solution = (u, p; k...) -> u,
                ind_plot = (1,2,3),
                k...)
    u0 = chart.u#isnothing(chart.data) ? chart.u : chart.data
    T = chart.Φ
    R = chart.R

    if draw_tangent
        Φ = s -> u0 .+ T * s
        Pt = mapreduce(Φ, hcat, chart.P)
        Pt = hcat(Pt, Φ(chart.P[1]))
        _color = is_on_boundary(chart) ? HSV(40,30,60) : HSV(200, 50, 50)
        _color = chart.label == Symbol() ? _color :  :red 

        points2d = map(x -> Point2f(x[ind_plot[1]], x[ind_plot[2]]), eachcol(Pt))
        faces = Makie.GeometryBasics.earcut_triangulate([points2d])
        mesh!(ax, Pt[[ind_plot[1],ind_plot[2],ind_plot[3]], :], faces; label = nothing, alpha = 1, color = _color, colorrange)
        if draw_edges
            @views lines!(ax, Pt[ind_plot[1], :], Pt[ind_plot[2], :], Pt[ind_plot[3], :]; label = nothing, color = :black, colorrange)
        end
    end

    # plots draw_circle
    if draw_circle
        Pt = mapreduce(Φ, hcat, [R .* [cos(θ), sin(θ)] for θ in LinRange(0,2pi, 100)])
        @views lines!(ax, Pt[ind_plot[1],:], Pt[ind_plot[2],:], Pt[ind_plot[3],:]; label = nothing, color = ind, colorrange, k...)
    end

    us = get_tangent(T, 1)
    ut = get_tangent(T, 2)

    if draw_frame
        R /= 3
        lines!(ax, [u0[1], u0[1] + R * us[1]], 
                   [u0[2], u0[2] + R * us[2]], 
                   [u0[3], u0[3] + R * us[3]]; label = nothing, color = :red, k...)
        lines!(ax, [u0[1], u0[1] + R * ut[1]], 
                   [u0[2], u0[2] + R * ut[2]],
                   [u0[3], u0[3] + R * ut[3]]; label = nothing, color = :green, k...)
    end

    if plot_center ||  chart.index == colorrange[end]
        scatter!(ax, [u0[ind_plot[1]]], [u0[ind_plot[2]]], [u0[ind_plot[3]]]; label = "$ind", color = ind, colorrange)
    end
    if put_ids
        text!(ax, [u0[ind_plot[1]]], [u0[ind_plot[2]] + 0.005], [u0[ind_plot[3]] + 0.005], text = "$ind")
    end
    ax
end

function plotd(chart::Chart; k...)
    f = Figure(); ax = Axis3(f[1,1], aspect = (1, 1, 1))
    plotd(ax, chart; ind = 1, colorrange = (1,2), k...)
    f
end

function plotd(ax, Σ::Atlas; k...)
    n = length(Σ)
    for (ind, chart) in pairs(Σ.atlas)
        plotd(ax, chart; ind, colorrange = (1,n), k... )
    end
    # axislegend(ax)
end

function plotd(Σ::Atlas; size = (700,700), k...)
    f = Figure(;size)
    ax = Axis3(f[1,1], aspect = :data, title = "$(length(Σ)) charts")
    plotd(ax, Σ; k...)
    f
end

function plotcenters(Σ::Atlas; size = (700,700), k...)
    f = Figure(;size)
    if ~isnothing(Σ[1].data)
        centers = mapreduce(x -> x.data, hcat, Σ.atlas)
    else
        centers = mapreduce(x -> x.u, hcat, Σ.atlas)
    end
    ax = Axis3(f[1,1], aspect = :data, title = "$(length(Σ)) charts")
    scatter!(ax, centers[1:3, :], color = centers[1,:])
    f
end

"""
$SIGNATURES

2d plot of the manifold.

## Optional arguments
- `draw_circle::Bool` plot the validity ball
- `draw_tangent = true` plot the tangent space
- `plot_center = true`, plot the center of the charts
- `put_ids = false` add the indices of the charts
- `draw_frame = false` plot the local frame on the tangent space
- `ind_plot=1:2` indices used to plot in the 2d space
"""
function plot2d(Σ::Atlas; size = (700,700), 
                        ind_plot = 1:2,
                        plot_circle = false,
                        colorrange = nothing,
                        draw_frame = false,
                        put_ids = false,
                        plot_center = false,
                        k...)
    f = Figure(;size)
    ax = Axis(f[1,1], xlabel = "x", ylabel = "y")
    pts = Polygon[]
    plot2d_improved!(pts, Σ; ind_plot)

    n = length(Σ)
    if n == 0
        return f
    end

    bd = findall(c -> is_on_boundary(c), Σ.atlas)
    if ~isempty(bd)
        poly!(ax, pts[bd], strokecolor = :black, strokewidth = 1, colorrange = (1,n), alpha = 0.2, color = :orange)
    end

    nbd = findall(c -> !is_on_boundary(c), Σ.atlas)
    if ~isempty(nbd)
        poly!(ax, pts[nbd], strokecolor = :black, strokewidth = 1, colorrange = (1,n), alpha = 0.2, color = eachindex(nbd))
    end

    if plot_center
        centers = [Point2f(c.u[1], c.u[2]) for c in Σ.atlas]
        scatter!(ax, centers)
    end

    if put_ids
        centers = [Point2f(c.u[1], c.u[2] + 0.005) for c in Σ.atlas]
        texts = ["$(c.index)" for c in Σ.atlas]
        text!(ax, centers; text = texts)
    end

    if plot_circle
        # @assert false
        pols = ([Circle(Point2f(c.u[ind_plot]), c.R) for c in Σ.atlas])
        poly!(ax, pols, color = :white, alpha = 0.2, strokewidth = 2, strokecolor=:black)
    end
    f
end 

using LinearAlgebra

function plot2d_improved!(_polygons, Σ; ind_plot = 1:2, k...)
    tmp = copy(Σ[1].u)
    for chart in Σ.atlas

        u0 = chart.u

        T = chart.Φ
        R = chart.R

        function Φ(tmp, s)
            if u0 isa StaticArray
                tmp = u0
                tmp = tmp + T*s
            else
                tmp .= u0
                mul!(tmp, T, s, 1, 1)
            end
            tmp
        end

        pts = Point2f[Point2f(Φ(tmp, pt)[ind_plot]) for pt in chart.P]
        push!(_polygons, Polygon(pts))

    end
end

function _plotd_improved(Σ::Atlas; size = (700,700), 
                        draw_circle = false,
                        colorrange = (1,2),
                        draw_tangent = true,
                        put_ids = false,
                        plot_center = false,
                        draw_edges = false,
                        plot_circle = false,
                        k...)
    # @assert false
    f = Figure(;size)
    ax = Axis3(f[1,1], aspect = (1, 1, 1))
    pts = Point3f[]
    faces = plotd_improved!(pts, Σ; k...)
    
    n = length(Σ)
    bd = findall(c -> is_on_boundary(c), Σ.atlas)

    @error "" (faces) (pts)

    mesh!(ax, pts, faces)
    # poly!(ax, pts[bd], strokecolor = :black, strokewidth = 1, colorrange = (1,n), alpha = 0.2, color = :orange)
    
    # nbd = findall(c -> !is_on_boundary(c), Σ.atlas)
    # poly!(ax, pts[nbd], strokecolor = :black, strokewidth = 1, colorrange = (1,n), alpha = 0.2, color = eachindex(nbd))
    
    if plot_center
        centers = mapreduce(x->x.u, hcat, Σ.atlas)
        scatter!(ax, centers[1:3,:])
    end

    f
end 

function _plotd_improved!(pts, Σ; k...)
    allfaces = zeros(UInt32,1,3)
    offset = 0
    for chart in Σ.atlas
        u0 = isnothing(chart.data) ? chart.u : chart.data
        
        T = chart.Φ
        R = chart.R
        
        Φ = s -> u0 .+ T * s

        pts2d = Point2f[Φ(pt)[1:2] for pt in chart.P]
        push!(pts2d, pts2d[1])
        pts3d = Point3f[Φ(pt) for pt in chart.P]
        # push!(pts3d, pts3d[1])
        _faces = Makie.GeometryBasics.earcut_triangulate([pts2d])
        faces = mapreduce(t->[t[1].i t[2].i t[3].i], vcat, _faces)
        @error "#########" chart.id length(chart.P) offset faces
        
        # @assert false
        append!(pts, (pts3d))
        allfaces = vcat(allfaces, faces .+ offset)
        offset += (length(chart.P))
        # break
    end
    return allfaces[2:end, :]
end

##############################################################################################################

function plotkd!(ax, tree)
    plot_node!(ax, tree.head)
    ax
end

function plotkd!(ax, node::BVHNode, S::Atlas)
    n = npoints(node)
    plotkd_sub!(ax, node, S; id = 1, ncolor = n)
end

function plotkd_sub!(ax, node::BVHNode, S; id = 1, ncolor = 2)
    box = node.box
    side = box.max - box.min
    if ~is_leaf(node)
        side .+= 0.01
    else
        side .+= 0.005
    end

    if 1==1#~is_leaf(node)
        poly!(ax, Rect(box.min[1], box.min[2], side[1], side[2]); 
            color = :white,
            strokewidth = 2,
            alpha = 0.1,
            # colorrange = (1, ncolor),
            strokecolor = node.is_leaf ? :red : :black,
            # linestyle = node.is_leaf ? :dashdot : :solid,
                )
    end

    if is_leaf(node)
        for id in node.chart_ids
            u = S[id].u
            R = S[id].R
            s = 0.025
            poly!(ax, Rect(u[1]-s, u[2]-s, 2s, 2s); 
                color = :red,
                strokewidth = 1,
                alpha = 0.3,
                # colorrange = (1, ncolor),
                # strokecolor = 1,
                # linestyle = node.is_leaf ? :dashdot : :solid,
                    )
            scatter!(ax, u[1], u[2]; color = :black)
            text!(ax, u[1], u[2]; text = "$id")
        end
    end
    
        
    plotkd_sub!(ax, node.left_child, S; id = id+1, ncolor)
    plotkd_sub!(ax, node.right_child, S; id = id+1, ncolor)
end

plotkd_sub!(ax, ::Nothing, S; k...) = ax

plotkd!(ax, box::AABB) = poly!(ax, Rect(box.min[1], box.min[2], box.max[1]-box.min[1], box.max[2]-box.min[2]), color=:white, strokewidth = 2,)
