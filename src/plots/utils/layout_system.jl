const LAYOUT_KWARGS = Dict{Symbol,Tuple{Any,String}}(
    :topo_plot_width => (0.10, "Width of individual plots in topo layout (as fraction of figure width)"),
    :topo_plot_height => (0.10, "Height of individual plots in topo layout (as fraction of figure height)"),
    :topo_margin => (0.02, "Margin between plots in topo layout"),
    :topo_scale_offset => (0.1, "Offset factor for scale plot position"),
    :topo_scale_pos => ((0.8, -0.8), "Fallback position for scale plot as (x, y) tuple"),
    :grid_dims => (nothing, "Grid dimensions as (rows, cols) tuple. If nothing, automatically calculated"),
    :grid_rowgap => (10, "Gap between rows in grid layout (in pixels)"),
    :grid_colgap => (10, "Gap between columns in grid layout (in pixels)"),
    :grid_skip_positions => (nothing, "Positions to skip in grid layout as vector of (row, col) tuples, e.g., [(2,1), (2,3)]"),
)

"""
    PlotLayout

A struct that defines how to arrange plots in a figure (`:single`, `:grid`, `:topo`, or `:custom`).
"""
struct PlotLayout
    type::Symbol
    dims::AbstractVector{Int}
    positions::Vector{Tuple{Float64,Float64}}
    channels::Vector{Symbol}
    metadata::Dict{Symbol,Any}
end

function _create_axis(fig, position; kwargs...)
    return Axis(
        position;
        xgridvisible = kwargs[:xgrid],
        ygridvisible = kwargs[:ygrid],
        xminorgridvisible = kwargs[:xminorgrid],
        yminorgridvisible = kwargs[:yminorgrid],
    )
end

function _create_topo_axis(fig, plot_width, plot_height, halign, valign; kwargs...)
    return Axis(
        fig[1, 1],
        width = Relative(plot_width),
        height = Relative(plot_height),
        halign = halign,
        valign = valign,
        xgridvisible = kwargs[:xgrid],
        ygridvisible = kwargs[:ygrid],
        xminorgridvisible = kwargs[:xminorgrid],
        yminorgridvisible = kwargs[:yminorgrid],
    )
end

function _apply_grid_spacing!(fig::Figure, plot_layout::PlotLayout)
    rowgap_val = plot_layout.metadata[:grid_rowgap]
    colgap_val = plot_layout.metadata[:grid_colgap]
    rows, cols = plot_layout.dims
    if !isnothing(rowgap_val) && rows > 1
        for row = 1:(rows-1)
            rowgap!(fig.layout, row, rowgap_val)
        end
    end
    if !isnothing(colgap_val) && cols > 1
        for col = 1:(cols-1)
            colgap!(fig.layout, col, colgap_val)
        end
    end
end

function _add_axis_and_channel!(axes::Vector{Axis}, channels::Vector{Symbol}, ax::Axis, channel::Symbol)
    push!(axes, ax)
    push!(channels, channel)
end

function _create_single_layout(channels::Vector{Symbol}; kwargs...)
    metadata = _merge_plot_kwargs(LAYOUT_KWARGS, kwargs; validate = false)
    return PlotLayout(:single, [1, 1], [(0.0, 0.0)], channels, metadata)
end

function _validate_skip_positions!(metadata::Dict{Symbol,Any})
    skip_positions = metadata[:grid_skip_positions]
    isnothing(skip_positions) && return nothing, 0
    for pos in skip_positions
        (pos isa Tuple && length(pos) == 2) || @minimal_error("grid_skip_positions must be a vector of (row, col) tuples, got: $pos")
    end
    skip_set = Set(skip_positions)
    metadata[:grid_skip_positions] = skip_set
    return skip_set, length(skip_set)
end

function _create_grid_layout(channels::Vector{Symbol}; kwargs...)
    isempty(channels) && @minimal_error("Cannot create grid layout with empty channel list")
    metadata = _merge_plot_kwargs(LAYOUT_KWARGS, kwargs; validate = false)
    skip_set, skip_count = _validate_skip_positions!(metadata)
    n_channels = length(channels)

    if !isnothing(metadata[:grid_dims])
        rows, cols = _validate_dims(metadata[:grid_dims], n_channels + skip_count)
        if !isnothing(skip_set)
            for (row, col) in skip_set
                (1 <= row <= rows && 1 <= col <= cols) || @minimal_error("Skip position ($row, $col) is outside grid bounds ($rows×$cols)")
            end
        end
        available = rows * cols - skip_count
        available < n_channels && @minimal_error(
            "Grid size ($rows×$cols) with $skip_count skipped positions has only $available available positions, but need $n_channels"
        )
    else
        rows, cols = _best_rect(n_channels + skip_count)
        if !isnothing(skip_set)
            for (row, col) in skip_set
                (1 <= row <= rows && 1 <= col <= cols) ||
                    @minimal_error("Skip position ($row, $col) is outside auto-calculated grid bounds ($rows×$cols)")
            end
        end
    end

    return PlotLayout(:grid, [rows, cols], Tuple{Float64,Float64}[], channels, metadata)
end

function _compute_topo_positions(
    labels::AbstractVector,
    xs::AbstractVector,
    ys::AbstractVector,
    channels::Vector{Symbol},
)::Vector{Tuple{Float64,Float64}}
    minx, maxx = isempty(xs) ? (0.0, 0.0) : Float64.(extrema(xs))
    miny, maxy = isempty(ys) ? (0.0, 0.0) : Float64.(extrema(ys))
    spacing = max(maxx - minx, maxy - miny) * 0.15

    positions = Vector{Tuple{Float64,Float64}}(undef, length(channels))
    new_plot_count = 0
    last_new_pos = (0.0, 0.0)

    for (ci, channel) in enumerate(channels)
        idx = findfirst(==(channel), labels)
        if !isnothing(idx)
            @inbounds positions[ci] = (Float64(xs[idx]), Float64(ys[idx]))
        else
            @minimal_warning "Channel $channel not found in layout, finding non-overlapping position"
            new_plot_count += 1
            new_pos = new_plot_count == 1 ? (maxx + spacing, maxy - spacing * 0.5) : (last_new_pos[1], last_new_pos[2] - spacing)
            positions[ci] = new_pos
            last_new_pos = new_pos
        end
    end
    return positions
end

function _create_topo_layout(layout::Layout, channels::Vector{Symbol}; kwargs...)
    topo_kwargs = _merge_plot_kwargs(LAYOUT_KWARGS, kwargs; validate = false)
    _ensure_coordinates_2d!(layout)
    positions = _compute_topo_positions(layout.data[!, :label], layout.data[!, :x2], layout.data[!, :y2], channels)
    return PlotLayout(:topo, [0, 0], positions, channels, topo_kwargs)
end

"""
    create_layout(layout_spec, channels, eeg_layout; kwargs...)

Create a `PlotLayout` object based on `layout_spec` (`:single`, `:grid`, `:topo`, or an existing `PlotLayout`).
"""
function create_layout(layout_spec::Union{Symbol,PlotLayout}, channels::Vector{Symbol}, eeg_layout::Union{Layout,Nothing}; kwargs...)
    layout_spec === :single && return _create_single_layout(channels; kwargs...)
    layout_spec === :grid && return _create_grid_layout(channels; kwargs...)
    layout_spec === :topo && return _create_topo_layout(eeg_layout, channels; kwargs...)
    layout_spec isa PlotLayout && return layout_spec
    @minimal_error "Invalid layout specification: $layout_spec. Must be :single, :grid, :topo, or a PlotLayout object."
end

"""
    create_custom_layout(positions::Vector{Tuple{Float64, Float64}}, channels::Vector{Symbol})

Create a custom layout with specific `(x, y)` positions for each channel.
"""
function create_custom_layout(positions::Vector{Tuple{Float64,Float64}}, channels::Vector{Symbol})
    length(positions) == length(channels) ||
        @minimal_error("n positions ($(length(positions))) must match n channels ($(length(channels)))")
    return PlotLayout(:custom, [0, 0], positions, channels, Dict{Symbol,Any}())
end

function _validate_dims(dims, n_items::Int)
    isnothing(dims) && return _best_rect(n_items)
    length(dims) == 2 || @minimal_error("dims must be a tuple of two positive integers (rows, cols), got: $dims")
    if dims[1] * dims[2] < n_items
        @minimal_warning "Grid size ($(dims[1])×$(dims[2])) is too small for $n_items items. Using optimal grid dimensions instead."
        return _best_rect(n_items)
    end
    return dims
end

function _best_rect(n::Int)
    n <= 0 && @minimal_error("n must be positive")
    n == 1 && return (1, 1)
    for r = isqrt(n):-1:2
        if n % r == 0
            c = n ÷ r
            c <= 3 * r && return (r, c)
            break
        end
    end
    rows = ceil(Int, sqrt(n))
    return (rows, ceil(Int, n / rows))
end

function _apply_layout!(fig::Figure, plot_layout::PlotLayout; kwargs...)
    axes = Axis[]
    channels = Symbol[]

    if plot_layout.type == :single
        ax = _create_axis(fig, fig[1, 1]; kwargs...)
        _add_axis_and_channel!(axes, channels, ax, plot_layout.channels[1])
    elseif plot_layout.type == :grid
        rows, cols = plot_layout.dims
        skip_positions = get(plot_layout.metadata, :grid_skip_positions, nothing)
        channel_idx = 1
        for row = 1:rows, col = 1:cols
            is_skipped = !isnothing(skip_positions) && (row, col) in skip_positions
            ax = _create_axis(fig, fig[row, col]; kwargs...)
            if channel_idx <= length(plot_layout.channels) && !is_skipped
                _add_axis_and_channel!(axes, channels, ax, plot_layout.channels[channel_idx])
                channel_idx += 1
            else
                hidespines!(ax)
                hidedecorations!(ax)
            end
        end
        _apply_grid_spacing!(fig, plot_layout)
    elseif plot_layout.type == :topo
        plot_width = plot_layout.metadata[:topo_plot_width]
        plot_height = plot_layout.metadata[:topo_plot_height]
        positions = plot_layout.positions
        isempty(positions) && @minimal_error("Cannot create topo layout with empty positions. Problem with layout creation!")

        minx, maxx = extrema(pos[1] for pos in positions)
        miny, maxy = extrema(pos[2] for pos in positions)
        xrange = maxx == minx ? 1.0 : maxx - minx
        yrange = maxy == miny ? 1.0 : maxy - miny

        for (channel, pos) in zip(plot_layout.channels, positions)
            ax = _create_topo_axis(fig, plot_width, plot_height, (pos[1] - minx) / xrange, (pos[2] - miny) / yrange; kwargs...)
            _add_axis_and_channel!(axes, channels, ax, channel)
        end
    elseif plot_layout.type == :custom
        for (channel, pos) in zip(plot_layout.channels, plot_layout.positions)
            ax = _create_axis(fig, fig[pos[1], pos[2]]; kwargs...)
            _add_axis_and_channel!(axes, channels, ax, channel)
        end
    end

    return axes, channels
end

function _resolve_layout_title(kwargs, fallback::String)
    get(kwargs, :show_title, true) == false && return ""
    pt = get(kwargs, :plot_title, nothing)
    return isnothing(pt) ? fallback : pt
end

function _set_grid_axis_properties!(ax::Axis, channel::Symbol, row::Int, col::Int, total_rows::Int, total_cols::Int; kwargs...)
    ax.title = _resolve_layout_title(kwargs, string(channel))
    if col != 1
        ax.ylabel = ""
        ax.yticklabelsvisible = false
    end
    if row != total_rows
        ax.xlabel = ""
        ax.xticklabelsvisible = false
    end
end

function _apply_layout_axis_properties!(axes::Vector{Axis}, plot_layout::PlotLayout; kwargs...)
    if plot_layout.type == :grid
        rows, cols = plot_layout.dims
        for (idx, ax) in enumerate(axes)
            _set_grid_axis_properties!(ax, plot_layout.channels[idx], fld(idx - 1, cols) + 1, mod(idx - 1, cols) + 1, rows, cols; kwargs...)
        end
    elseif plot_layout.type == :topo
        for (idx, ax) in enumerate(axes)
            ax.title = _resolve_layout_title(kwargs, string(plot_layout.channels[idx]))
            hidedecorations!(ax, grid = false, ticks = true, ticklabels = true)
            hidespines!(ax)
        end
    end
    for ax in axes
        _draw_inner_title_if_requested!(ax; kwargs...)
    end
end

function _apply_axis_properties!(ax::Axis; kwargs...)
    ax.title = _resolve_layout_title(kwargs, "")
    if haskey(kwargs, :plot_title_fontsize) && !isnothing(kwargs[:plot_title_fontsize])
        ax.titlesize = kwargs[:plot_title_fontsize]
    end
    ax.xlabel = kwargs[:xlabel]
    ax.ylabel = kwargs[:ylabel]

    if haskey(kwargs, :time_unit) && kwargs[:time_unit] == :ms
        ax.xtickformat = values -> [string(round(Int, v * 1000)) for v in values]
        kwargs[:xlabel] == "Time (s)" && (ax.xlabel = "Time (ms)")
    end

    !isnothing(kwargs[:xlim]) && xlims!(ax, kwargs[:xlim])
    !isnothing(kwargs[:ylim]) && ylims!(ax, kwargs[:ylim])
    ax.yreversed = kwargs[:yreversed]

    haskey(kwargs, :xticks) && !isnothing(kwargs[:xticks]) && (ax.xticks = kwargs[:xticks])
    haskey(kwargs, :yticks) && !isnothing(kwargs[:yticks]) && (ax.yticks = kwargs[:yticks])

    ax.xgridvisible = kwargs[:xgrid]
    ax.ygridvisible = kwargs[:ygrid]
    ax.xminorgridvisible = kwargs[:xminorgrid]
    ax.yminorgridvisible = kwargs[:yminorgrid]

    _draw_inner_title_if_requested!(ax; kwargs...)
end

function _draw_inner_title_if_requested!(ax::Axis; kwargs...)
    pos = get(kwargs, :plot_title_position, nothing)
    isnothing(pos) && return
    title_str = ax.title[]
    isempty(title_str) && return
    text!(
        ax,
        pos[1],
        pos[2],
        text = title_str,
        space = :relative,
        align = kwargs[:plot_title_align],
        fontsize = kwargs[:plot_title_fontsize],
        font = :bold,
    )
    ax.title = ""
end
