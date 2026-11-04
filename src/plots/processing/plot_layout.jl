# === DEFAULT KEYWORD ARGUMENTS ===
const PLOT_LAYOUT_HEAD_KWARGS = Dict{Symbol,Tuple{Any,String}}(
    :head_color => (nothing, "Color of the head shape outline."),
    :head_linewidth => (2, "Line width of the head shape outline."),
    :head_radius => (1, "Radius of the head shape (normalized units)."),
)

const PLOT_LAYOUT_POINT_KWARGS = Dict{Symbol,Tuple{Any,String}}(
    :point_plot => (true, "Whether to plot electrode points."),
    :point_marker => (:circle, "Marker style for electrode points."),
    :point_markersize => (12, "Size of electrode point markers."),
    :point_color => (nothing, "Color of electrode points."),
)

const PLOT_LAYOUT_LABEL_KWARGS = Dict{Symbol,Tuple{Any,String}}(
    :label_plot => (true, "Whether to plot electrode labels."),
    :label_fontsize => (20, "Font size for electrode labels."),
    :label_color => (nothing, "Color of electrode labels."),
    :label_xoffset => (0, "X-axis offset for electrode labels."),
    :label_yoffset => (0, "Y-axis offset for electrode labels."),
    :label_zoffset => (0, "Z-axis offset for electrode labels (3D only)."),
    :label_align => ((:left, :bottom), "Text alignment for electrode labels (e.g., (:center, :center))."),
)

const PLOT_LAYOUT_ROI_KWARGS = Dict{Symbol,Tuple{Any,String}}(
    :roi_border_size => (10, "Size of the border around ROI points."),
    :roi_linecolor => (nothing, "Color of ROI outline."),
    :roi_linewidth => (2, "Line width of ROI outline."),
    :roi_fill => (false, "Whether to fill the ROI area."),
    :roi_fillcolor => (:gray, "Color of ROI fill."),
    :roi_fillalpha => (0.2, "Transparency of ROI fill."),
)

const _UNIT_HEAD_OUTLINE = let
    n_circle = 100
    circle_pts = [Point2f(cos(θ), sin(θ)) for θ in range(-Float32(π), Float32(π), length = n_circle)]
    n_ear = 30
    ear_r = 1.0f0 / 7.0f0
    right_ear = [Point2f(1.0f0 + ear_r * cos(θ), ear_r * sin(θ)) for θ in range(-Float32(π) / 2, Float32(π) / 2, length = n_ear)]
    left_ear = [Point2f(-1.0f0 - ear_r * cos(θ), ear_r * sin(θ)) for θ in range(-Float32(π) / 2, Float32(π) / 2, length = n_ear)]
    nose = Point2f[(-0.1f0, 1.0f0), (0.0f0, 1.15f0), (0.1f0, 1.0f0)]
    nan_pt = Point2f(NaN32, NaN32)

    pts = Point2f[]
    sizehint!(pts, n_circle + 2 * n_ear + 6)
    append!(pts, circle_pts)
    push!(pts, nan_pt)
    append!(pts, right_ear)
    push!(pts, nan_pt)
    append!(pts, left_ear)
    push!(pts, nan_pt)
    append!(pts, nose)
    pts
end

const _ROI_CIRCLE_UNIT = [Point2{Float64}(sin(θ), cos(θ)) for θ = 0:(2π/361):2π]

@inline _layout_kw(kwargs, dict::Dict{Symbol,Tuple{Any,String}}, key::Symbol) = haskey(kwargs, key) ? kwargs[key] : dict[key][1]

@inline _cross2d(o::Point2{Float64}, a::Point2{Float64}, b::Point2{Float64})::Float64 =
    (a[1] - o[1]) * (b[2] - o[2]) - (a[2] - o[2]) * (b[1] - o[1])

function _monotone_chain_hull!(pts::Vector{Point2{Float64}})::Vector{Point2{Float64}}
    n = length(pts)
    n <= 1 && return copy(pts)
    sort!(pts; by = p -> (p[1], p[2]))

    hull = Vector{Point2{Float64}}(undef, 2 * n)
    k = 0
    @inbounds for i = 1:n
        p = pts[i]
        while k >= 2 && _cross2d(hull[k-1], hull[k], p) <= 0.0
            k -= 1
        end
        k += 1
        hull[k] = p
    end
    t = k + 1
    @inbounds for i = (n-1):-1:1
        p = pts[i]
        while k >= t && _cross2d(hull[k-1], hull[k], p) <= 0.0
            k -= 1
        end
        k += 1
        hull[k] = p
    end
    resize!(hull, max(1, k - 1))
    return hull
end

function _convex_hull_with_border(xpos::AbstractVector{<:Real}, ypos::AbstractVector{<:Real}, border_size::Real)::Vector{Point2f}
    centers = [Point2{Float64}(Float64(xpos[i]), Float64(ypos[i])) for i in eachindex(xpos)]
    base_hull = _monotone_chain_hull!(centers)
    r = abs(Float64(border_size))

    if r == 0.0
        n_b = length(base_hull)
        out = Vector{Point2f}(undef, n_b + 1)
        @inbounds for i = 1:n_b
            out[i] = Point2f(base_hull[i][1], base_hull[i][2])
        end
        out[n_b+1] = out[1]
        return out
    end

    n_unit = length(_ROI_CIRCLE_UNIT)
    pts = Vector{Point2{Float64}}(undef, length(base_hull) * n_unit)
    idx = 0
    @inbounds for c in base_hull
        cx, cy = c[1], c[2]
        for u in _ROI_CIRCLE_UNIT
            idx += 1
            pts[idx] = Point2{Float64}(muladd(r, u[1], cx), muladd(r, u[2], cy))
        end
    end

    hull = _monotone_chain_hull!(pts)
    n_h = length(hull)
    out = Vector{Point2f}(undef, n_h + 1)
    @inbounds for i = 1:n_h
        out[i] = Point2f(hull[i][1], hull[i][2])
    end
    out[n_h+1] = out[1]
    return out
end

"""
    plot_layout_2d!(fig::Figure, ax::Axis, layout::Layout; neighbors::Bool=false, correlation_matrix::Union{DataFrame,Nothing}=nothing, kwargs...) -> Axis
    plot_layout_2d(layout::Layout; neighbors::Bool=false, correlation_matrix::Union{DataFrame,Nothing}=nothing, display_plot::Bool=true, kwargs...) -> NamedTuple

Plot a 2D EEG electrode layout with customizable head outline, electrode markers, labels, and optional hover neighbor or correlation overlays.

$(_generate_kwargs_doc(PLOT_LAYOUT_HEAD_KWARGS))
$(_generate_kwargs_doc(PLOT_LAYOUT_POINT_KWARGS))
$(_generate_kwargs_doc(PLOT_LAYOUT_LABEL_KWARGS))
"""
function plot_layout_2d!(
    fig::Figure,
    ax::Axis,
    layout::Layout;
    neighbors::Bool = false,
    correlation_matrix::Union{DataFrame,Nothing} = nothing,
    kwargs...,
)
    _ensure_coordinates_2d!(layout)

    head_color = _layout_kw(kwargs, PLOT_LAYOUT_HEAD_KWARGS, :head_color)
    head_lw = _layout_kw(kwargs, PLOT_LAYOUT_HEAD_KWARGS, :head_linewidth)
    radius = _layout_kw(kwargs, PLOT_LAYOUT_HEAD_KWARGS, :head_radius)

    plot_points = _layout_kw(kwargs, PLOT_LAYOUT_POINT_KWARGS, :point_plot)
    point_marker = _layout_kw(kwargs, PLOT_LAYOUT_POINT_KWARGS, :point_marker)
    point_ms = _layout_kw(kwargs, PLOT_LAYOUT_POINT_KWARGS, :point_markersize)
    point_color = _layout_kw(kwargs, PLOT_LAYOUT_POINT_KWARGS, :point_color)

    plot_labels = _layout_kw(kwargs, PLOT_LAYOUT_LABEL_KWARGS, :label_plot)
    label_fs = _layout_kw(kwargs, PLOT_LAYOUT_LABEL_KWARGS, :label_fontsize)
    label_color = _layout_kw(kwargs, PLOT_LAYOUT_LABEL_KWARGS, :label_color)
    xoffset = _layout_kw(kwargs, PLOT_LAYOUT_LABEL_KWARGS, :label_xoffset)
    yoffset = _layout_kw(kwargs, PLOT_LAYOUT_LABEL_KWARGS, :label_yoffset)
    label_align = _layout_kw(kwargs, PLOT_LAYOUT_LABEL_KWARGS, :label_align)

    actual_head_color = _resolve_theme_linecolor(ax, head_color)
    head_pts = radius == 1 ? _UNIT_HEAD_OUTLINE : (_UNIT_HEAD_OUTLINE .* Float32(radius))
    lines!(ax, head_pts; color = actual_head_color, linewidth = head_lw)

    x2 = layout.data.x2
    y2 = layout.data.y2
    positions = Point2f.(x2, y2)

    if !isnothing(correlation_matrix)
        point_kw = (point_markersize = point_ms, point_color = point_color, point_marker = point_marker)
        label_kw = (label_color = label_color, label_fontsize = label_fs, label_align = label_align)
        _add_interactive_correlation_points!(
            fig,
            ax,
            layout.data,
            positions,
            correlation_matrix,
            point_kw,
            label_kw,
            xoffset,
            yoffset,
            plot_labels,
        )
    else
        if plot_points
            scatter!(ax, positions; marker = point_marker, markersize = point_ms, color = _resolve_theme_textcolor(ax, point_color))
        end

        if plot_labels && !isempty(layout.data.label)
            labels = String.(layout.data.label)
            label_positions =
                (xoffset == 0 && yoffset == 0) ? positions : [Point2f(x2[i] + xoffset, y2[i] + yoffset) for i in eachindex(x2)]
            text!(
                ax,
                label_positions;
                text = labels,
                fontsize = label_fs,
                color = _resolve_theme_textcolor(ax, label_color),
                align = label_align,
            )
        end

        if neighbors
            if isnothing(layout.neighbors)
                @minimal_warning "Layout has no neighbors data. Set neighbors=false or calculate neighbors first."
            else
                _add_interactive_points!(fig, ax, layout.data, layout.neighbors, positions)
            end
        end
    end

    hidedecorations!(ax)
    hidespines!(ax)
    return ax
end

plot_layout_2d!(fig, ax, dat::EegFunData; kwargs...) = plot_layout_2d!(fig, ax, dat.layout; kwargs...)

"""
    plot_layout_2d(layout::Layout; neighbors::Bool=false, correlation_matrix=nothing, display_plot::Bool=true, kwargs...) -> NamedTuple

Create a figure and plot a 2D EEG electrode layout. See [`plot_layout_2d!`](@ref) for keyword options.
"""
function plot_layout_2d(
    layout::Layout;
    neighbors::Bool = false,
    correlation_matrix::Union{DataFrame,Nothing} = nothing,
    display_plot::Bool = true,
    kwargs...,
)
    fig = _eegfun_figure()
    ax = Axis(fig[1, 1], aspect = DataAspect(), xautolimitmargin = (0.1, 0.15), yautolimitmargin = (0.1, 0.1))
    plot_layout_2d!(fig, ax, layout; neighbors = neighbors, correlation_matrix = correlation_matrix, kwargs...)
    _set_window_title("Layout 2D")
    display_plot && _display_figure(fig)
    return (fig = fig, axes = [ax])
end

plot_layout_2d(dat::EegFunData; kwargs...) = plot_layout_2d(dat.layout; kwargs...)

"""
    plot_layout(args...; kwargs...)

2D sensor layout alias for [`plot_layout_2d`](@ref).
"""
plot_layout(args...; kwargs...) = plot_layout_2d(args...; kwargs...)

"""
    plot_layout!(args...; kwargs...)

In-place 2D sensor layout alias for [`plot_layout_2d!`](@ref).
"""
plot_layout!(args...; kwargs...) = plot_layout_2d!(args...; kwargs...)

"""
    add_topo_rois!(ax::Axis, layout::Layout, rois::Vector{<:Vector{Symbol}}; kwargs...)

Draw padded convex-hull region-of-interest (ROI) outlines or fills around groups of electrodes on `ax`.

$(_generate_kwargs_doc(PLOT_LAYOUT_ROI_KWARGS))
"""
function add_topo_rois!(ax::Axis, layout::Layout, rois::Vector{<:Vector{Symbol}}; kwargs...)
    n_rois = length(rois)
    n_rois == 0 && return nothing

    roi_kwargs = _merge_plot_kwargs(PLOT_LAYOUT_ROI_KWARGS, kwargs; validate = false)
    if haskey(roi_kwargs, :roi_color) && !haskey(kwargs, :roi_linecolor)
        roi_kwargs[:roi_linecolor] = roi_kwargs[:roi_color]
    end

    merged_kwargs = Dict{Symbol,Vector}()
    for (key, value) in roi_kwargs
        if !(value isa Vector)
            merged_kwargs[key] = fill(value, n_rois)
        elseif length(value) != n_rois
            @minimal_error "Length of $(key) ($(length(value))) must match number of ROIs ($n_rois)"
        else
            merged_kwargs[key] = value
        end
    end

    labels = layout.data.label
    x2 = layout.data.x2
    y2 = layout.data.y2
    for (i, roi) in enumerate(rois)
        roi_idx = findall(in(roi), labels)
        if isempty(roi_idx)
            @minimal_warning "No electrodes found for ROI $i"
            continue
        end

        hull_pts = _convex_hull_with_border(@view(x2[roi_idx]), @view(y2[roi_idx]), merged_kwargs[:roi_border_size][i])
        linecolor_i = _resolve_theme_linecolor(ax, merged_kwargs[:roi_linecolor][i])
        lw_i = merged_kwargs[:roi_linewidth][i]

        if merged_kwargs[:roi_fill][i]
            poly!(
                ax,
                hull_pts;
                color = (merged_kwargs[:roi_fillcolor][i], merged_kwargs[:roi_fillalpha][i]),
                strokecolor = linecolor_i,
                strokewidth = lw_i,
            )
        else
            lines!(ax, hull_pts; color = linecolor_i, linewidth = lw_i)
        end
    end
    return nothing
end

"""
    plot_layout_3d!(fig::Figure, ax::Axis3, layout::Layout; neighbors::Bool=false, kwargs...) -> Axis3
    plot_layout_3d(layout::Layout; neighbors::Bool=false, display_plot::Bool=true, kwargs...) -> NamedTuple

Plot a 3D EEG electrode layout with electrode markers, labels, and optional hover neighbor connections.

$(_generate_kwargs_doc(PLOT_LAYOUT_HEAD_KWARGS))
$(_generate_kwargs_doc(PLOT_LAYOUT_POINT_KWARGS))
$(_generate_kwargs_doc(PLOT_LAYOUT_LABEL_KWARGS))
"""
function plot_layout_3d!(fig::Figure, ax::Axis3, layout::Layout; neighbors::Bool = false, kwargs...)
    _ensure_coordinates_3d!(layout)

    plot_points = _layout_kw(kwargs, PLOT_LAYOUT_POINT_KWARGS, :point_plot)
    point_marker = _layout_kw(kwargs, PLOT_LAYOUT_POINT_KWARGS, :point_marker)
    point_ms = _layout_kw(kwargs, PLOT_LAYOUT_POINT_KWARGS, :point_markersize)
    point_color = _layout_kw(kwargs, PLOT_LAYOUT_POINT_KWARGS, :point_color)

    plot_labels = _layout_kw(kwargs, PLOT_LAYOUT_LABEL_KWARGS, :label_plot)
    label_fs = _layout_kw(kwargs, PLOT_LAYOUT_LABEL_KWARGS, :label_fontsize)
    label_color = _layout_kw(kwargs, PLOT_LAYOUT_LABEL_KWARGS, :label_color)
    xoffset = _layout_kw(kwargs, PLOT_LAYOUT_LABEL_KWARGS, :label_xoffset)
    yoffset = _layout_kw(kwargs, PLOT_LAYOUT_LABEL_KWARGS, :label_yoffset)
    zoffset = _layout_kw(kwargs, PLOT_LAYOUT_LABEL_KWARGS, :label_zoffset)
    label_align = _layout_kw(kwargs, PLOT_LAYOUT_LABEL_KWARGS, :label_align)

    x3 = layout.data.x3
    y3 = layout.data.y3
    z3 = layout.data.z3

    if plot_points
        scatter!(ax, x3, y3, z3; marker = point_marker, markersize = point_ms, color = _resolve_theme_textcolor(ax, point_color))
    end

    if plot_labels && !isempty(layout.data.label)
        labels = String.(layout.data.label)
        label_positions = [Point3f(x3[i] + xoffset, y3[i] + yoffset, z3[i] + zoffset) for i in eachindex(x3)]
        text!(
            ax,
            label_positions;
            text = labels,
            fontsize = label_fs,
            color = _resolve_theme_textcolor(ax, label_color),
            align = label_align,
        )
    end

    if neighbors
        if isnothing(layout.neighbors)
            @minimal_warning "Layout has no neighbors data. Set neighbors=false or calculate neighbors first."
        else
            positions = Point3f.(x3, y3, z3)
            _add_interactive_points!(fig, ax, layout.data, layout.neighbors, positions, true)
        end
    end

    hidedecorations!(ax)
    hidespines!(ax)
    return ax
end

plot_layout_3d!(fig, ax, dat::EegFunData; kwargs...) = plot_layout_3d!(fig, ax, dat.layout; kwargs...)

"""
    plot_layout_3d(layout::Layout; neighbors::Bool=false, display_plot::Bool=true, kwargs...) -> NamedTuple

Create a figure and plot a 3D EEG electrode layout. See [`plot_layout_3d!`](@ref) for keyword options.
"""
function plot_layout_3d(layout::Layout; neighbors::Bool = false, display_plot::Bool = true, kwargs...)
    fig = _eegfun_figure()
    ax = Axis3(fig[1, 1], xautolimitmargin = (0.1, 0.15), yautolimitmargin = (0.1, 0.1), zautolimitmargin = (0.1, 0.1))
    plot_layout_3d!(fig, ax, layout; neighbors = neighbors, kwargs...)
    _set_window_title("Layout 3D")
    display_plot && _display_figure(fig)
    return (fig = fig, axes = [ax], ax = ax)
end

plot_layout_3d(dat::EegFunData; kwargs...) = plot_layout_3d(dat.layout; kwargs...)

function _add_interactive_points!(
    fig::Figure,
    ax::Union{Axis,Axis3},
    layout::DataFrame,
    neighbors::Union{Nothing,AbstractDict} = nothing,
    positions::Union{Vector{Point2f},Vector{Point3f}} = Point2f[],
    is_3d::Bool = false,
)
    isnothing(neighbors) && return nothing
    n_ch = length(layout.label)
    base_size = 15
    hover_size = 25
    sizes = Observable(fill(base_size, n_ch))

    label_to_idx = Dict{Symbol,Int}(Symbol(lbl) => i for (i, lbl) in enumerate(layout.label))
    neighbor_indices = Vector{Vector{Int}}(undef, n_ch)
    max_neighbors = 0
    for (i, label) in enumerate(layout.label)
        nb_vec = Int[]
        sym = Symbol(label)
        if haskey(neighbors, sym)
            for nb in neighbors[sym].channels
                nb_idx = get(label_to_idx, nb, 0)
                nb_idx > 0 && push!(nb_vec, nb_idx)
            end
        end
        neighbor_indices[i] = nb_vec
        length(nb_vec) > max_neighbors && (max_neighbors = length(nb_vec))
    end

    new_sizes = fill(base_size, n_ch)
    new_lines = is_3d ? Vector{Point3f}(undef, max_neighbors * 2) : Vector{Point2f}(undef, max_neighbors * 2)

    p = scatter!(ax, positions; color = _resolve_theme_textcolor(ax, nothing), markersize = sizes, inspectable = true, markerspace = :pixel)
    linesegments = Observable(is_3d ? Point3f[] : Point2f[])
    lines!(ax, linesegments, color = :gray, linewidth = 3)

    on(events(fig).mouseposition) do _
        plt, i = pick(fig)
        if plt == p && !isnothing(i) && 1 <= i <= n_ch
            fill!(new_sizes, base_size)
            new_sizes[i] = hover_size
            sizes[] = new_sizes

            hovered_pos = positions[i]
            line_count = 0
            for nb_idx in neighbor_indices[i]
                new_lines[line_count+1] = hovered_pos
                new_lines[line_count+2] = positions[nb_idx]
                line_count += 2
            end
            linesegments[] = new_lines[1:line_count]
        end
    end
end

function _add_interactive_correlation_points!(
    fig::Figure,
    ax::Axis,
    layout::DataFrame,
    positions::Vector{Point2f},
    correlation_matrix::DataFrame,
    point_kwargs,
    label_kwargs,
    xoffset,
    yoffset,
    plot_labels,
)
    n_channels = length(layout.label)
    n_channels == 0 && return nothing

    base_size = Float64(point_kwargs[:point_markersize])
    hover_size = base_size * 1.5

    sizes = Observable(fill(base_size, n_channels))
    corr_values = Observable(fill(0.0, n_channels))
    corr_texts = Observable(fill("", n_channels))

    actual_label_color = _resolve_theme_textcolor(ax, label_kwargs[:label_color])
    actual_point_color = _resolve_theme_textcolor(ax, point_kwargs[:point_color])

    channel_label_plot = if plot_labels
        label_positions =
            (xoffset == 0 && yoffset == 0) ? positions : [Point2f(layout.x2[i] + xoffset, layout.y2[i] + yoffset) for i = 1:n_channels]
        text!(
            ax,
            label_positions;
            text = String.(layout.label),
            fontsize = label_kwargs[:label_fontsize],
            color = actual_label_color,
            align = label_kwargs[:label_align],
        )
    else
        nothing
    end

    corr_text_plot = text!(
        ax,
        positions;
        text = corr_texts,
        fontsize = 16,
        color = _resolve_theme_textcolor(ax, nothing),
        align = label_kwargs[:label_align],
        visible = false,
    )

    p = scatter!(
        ax,
        positions;
        color = actual_point_color,
        markersize = sizes,
        marker = point_kwargs[:point_marker],
        inspectable = true,
        markerspace = :pixel,
    )

    corr_scatter = scatter!(
        ax,
        positions;
        color = corr_values,
        markersize = sizes,
        marker = point_kwargs[:point_marker],
        inspectable = true,
        markerspace = :pixel,
        colormap = :RdYlGn,
        colorrange = (-1.0, 1.0),
        visible = false,
    )

    row_to_idx = Dict{Symbol,Int}(Symbol(r) => k for (k, r) in enumerate(correlation_matrix[!, :row]))
    new_sizes = fill(base_size, n_channels)
    new_corr_values = fill(0.0, n_channels)
    new_corr_strings = fill("", n_channels)
    last_hovered = Ref{Int}(0)

    on(events(fig).mouseposition) do _
        plt, i = pick(fig)
        if (plt == p || plt == corr_scatter) && !isnothing(i) && 1 <= i <= n_channels
            i == last_hovered[] && return nothing
            last_hovered[] = i

            fill!(new_sizes, base_size)
            new_sizes[i] = hover_size
            sizes[] = copy(new_sizes)

            row_idx = get(row_to_idx, Symbol(layout.label[i]), 0)
            if row_idx > 0
                for (j, label) in enumerate(layout.label)
                    corr_val = hasproperty(correlation_matrix, label) ? Float64(correlation_matrix[row_idx, label]) : 0.0
                    new_corr_values[j] = corr_val
                    new_corr_strings[j] = string(round(corr_val, digits = 2))
                end
            else
                fill!(new_corr_values, 0.0)
                fill!(new_corr_strings, "")
            end

            corr_values[] = copy(new_corr_values)
            corr_texts[] = copy(new_corr_strings)

            !isnothing(channel_label_plot) && (channel_label_plot.visible = false)
            corr_text_plot.visible = true
            p.visible = false
            corr_scatter.visible = true
        elseif last_hovered[] != 0
            last_hovered[] = 0
            fill!(new_sizes, base_size)
            sizes[] = copy(new_sizes)

            !isnothing(channel_label_plot) && (channel_label_plot.visible = true)
            corr_text_plot.visible = false
            p.visible = true
            corr_scatter.visible = false
        end
        return nothing
    end
    return nothing
end
