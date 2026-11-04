# Note: PLOT_TOPOGRAPHY_KWARGS is defined in plot_ica.jl

const _UNIT_CIRCLE_200 = Tuple{Float32,Float32}[(Float32(cos(a)), Float32(sin(a))) for a in range(0, 2π, length = 200)]

function _symmetric_colorrange(A::AbstractArray)
    max_abs = 0.0
    has_valid = false
    @inbounds for x in A
        v = Float64(x)
        if !isnan(v)
            has_valid = true
            av = abs(v)
            av > max_abs && (max_abs = av)
        end
    end
    (!has_valid || max_abs ≈ 0.0) && return (-1.0, 1.0)
    return (-max_abs, max_abs)
end

function _render_topo_surface!(
    fig::Figure,
    ax::Axis,
    channel_data::AbstractVector{<:Real},
    layout::Layout;
    data_interp = nothing,
    x_bounds = nothing,
    y_bounds = nothing,
    method::Symbol = :thin_plate,
    gridscale::Int = 75,
    colormap = nothing,
    colorrange = nothing,
    num_levels::Int = 20,
    kwargs...,
)
    if isnothing(data_interp)
        ch_f64 = channel_data isa AbstractVector{Float64} ? channel_data : Float64.(channel_data)
        data_interp, x_bounds, y_bounds = _topo_interpolate(ch_f64, layout, gridscale; method = method)
    end

    cr = !isnothing(colorrange) ? colorrange : _symmetric_colorrange(data_interp)

    co = contourf!(
        ax,
        range(x_bounds[1], x_bounds[2], length = gridscale),
        range(y_bounds[1], y_bounds[2], length = gridscale),
        data_interp,
        levels = range(cr[1], cr[2], length = num_levels);
        colormap = _resolve_theme_colormap(ax, colormap),
        extendlow = :auto,
        extendhigh = :auto,
        nan_color = :transparent,
    )
    co.colorrange = cr

    _draw_smooth_circle_mask!(ax, x_bounds, y_bounds)
    plot_layout_2d!(fig, ax, layout; kwargs...)

    return co, x_bounds, y_bounds
end

function _plot_topography!(
    fig::Figure,
    ax::Axis,
    dat::DataFrame,
    layout::Layout;
    data_interp = nothing,
    x_bounds = nothing,
    y_bounds = nothing,
    kwargs...,
)
    if !has_valid_coordinates(layout)
        @minimal_warning "Cannot create topographic plot: layout has no spatial coordinates."
        return
    end
    _ensure_coordinates_2d!(layout)

    plot_kwargs = _merge_plot_kwargs(PLOT_TOPOGRAPHY_KWARGS, kwargs)

    method = pop!(plot_kwargs, :method)
    gridscale = pop!(plot_kwargs, :gridscale)
    cb_pos = _get_colorbar_position(pop!(plot_kwargs, :colorbar_position), 1:1, 1:1)
    colorrange = pop!(plot_kwargs, :colorrange)
    colormap = pop!(plot_kwargs, :colormap)
    num_levels = pop!(plot_kwargs, :num_levels)
    show_title = pop!(plot_kwargs, :show_title, true)

    if !show_title
        ax.title = ""
    elseif !isnothing(plot_kwargs[:plot_title]) && !isempty(plot_kwargs[:plot_title])
        ax.title = plot_kwargs[:plot_title]
        ax.titlesize = plot_kwargs[:plot_title_fontsize]
    end

    channel_data = isnothing(data_interp) ? Float64[mean(dat[!, lbl]) for lbl in layout.data.label] : Float64[]
    colorbar_kwargs = pop!(plot_kwargs, :colorbar_kwargs, (;))
    colorbar_plot = pop!(plot_kwargs, :colorbar_plot)
    highlight_channels = pop!(plot_kwargs, :highlight_channels, nothing)

    co, _, _ = _render_topo_surface!(
        fig,
        ax,
        channel_data,
        layout;
        data_interp = data_interp,
        x_bounds = x_bounds,
        y_bounds = y_bounds,
        method = method,
        gridscale = gridscale,
        colormap = colormap,
        colorrange = colorrange,
        num_levels = num_levels,
        plot_kwargs...,
    )

    !isnothing(highlight_channels) && _draw_highlight_channels!(ax, highlight_channels, layout)
    colorbar_plot && Colorbar(fig[cb_pos...], co; colorbar_kwargs...)

    return (fig = fig, axes = [ax])
end

function _draw_highlight_channels!(ax::Axis, groups, layout::Layout)
    group_list = groups isa AbstractVector ? groups : [groups]
    labels_col = layout.data.label
    x2_col = layout.data.x2
    y2_col = layout.data.y2

    for group in group_list
        color = get(group, :color, :white)
        sz = get(group, :size, 8)
        marker = get(group, :marker, :circle)
        show_label = get(group, :label, true)
        label_fontsize = get(group, :fontsize, 10)

        xs = Float64[]
        ys = Float64[]
        labels = Symbol[]
        for ch in group.channels
            idx = findfirst(==(ch), labels_col)
            if !isnothing(idx)
                push!(xs, x2_col[idx])
                push!(ys, y2_col[idx])
                push!(labels, ch)
            end
        end

        if !isempty(xs)
            scatter!(ax, xs, ys, color = color, marker = marker, markersize = sz, strokewidth = 1, strokecolor = :black)
            if show_label
                text!(
                    ax,
                    Point2f.(xs, ys);
                    text = string.(labels),
                    align = (:center, :bottom),
                    fontsize = label_fontsize,
                    color = color,
                    offset = (0, sz / 2 + 2),
                )
            end
        end
    end
end

function _expand_topo_time_bins(
    dataset::SingleDataFrameEeg,
    n_topo::Int,
    channel_selection,
    sample_selection,
    interval_selection::Interval,
    time_unit::Symbol;
    include_cond_prefix::Bool = false,
)
    dat_selected =
        subset(dataset, channel_selection = channel_selection, sample_selection = sample_selection, interval_selection = interval_selection)
    t_min, t_max = extrema(time_vector(dat_selected))
    bin_edges = range(t_min, t_max, length = n_topo + 1)

    bins = typeof(dataset)[]
    sizehint!(bins, n_topo)
    for i = 1:n_topo
        bin_start, bin_stop = bin_edges[i], bin_edges[i+1]
        bin_dat = subset(
            dataset,
            channel_selection = channel_selection,
            sample_selection = x -> x.time .>= bin_start .&& x.time .<= bin_stop,
            interval_selection = interval_selection,
        )
        if hasproperty(bin_dat, :condition_name)
            range_str =
                time_unit == :ms ? @sprintf("%.0f – %.0f ms", bin_start * 1000, bin_stop * 1000) :
                @sprintf("%.3f – %.3f s", bin_start, bin_stop)
            bin_dat.condition_name = include_cond_prefix ? "$(bin_dat.condition_name) ($range_str)" : range_str
        end
        push!(bins, bin_dat)
    end
    return bins
end

"""
    plot_topography(filepath::String; input_dir = pwd(), participant_selection = participants(), kwargs...)
    plot_topography(dat::SingleDataFrameEeg; channel_selection = channels(), sample_selection = samples(), interval_selection = times(), baseline_interval = nothing, n_topo = nothing, scrubber = false, display_plot = true, interactive = true, kwargs...)
    plot_topography(dat::AbstractVector{<:SingleDataFrameEeg}; kwargs...)
    plot_topography(dat::AbstractVector{<:EpochData}; kwargs...)
    plot_topography(dat::MultiDataFrameEeg, epoch::Int; kwargs...)

Plot 2D scalp topographies for ERP, continuous, or epoched EEG data (or `.jld2` files/patterns).

Supports multi-bin time expansion via `n_topo`, interactive playback via `scrubber = true`,
and Shift+drag spatial channel selection with right-click ERP plotting.

$(_generate_kwargs_doc(PLOT_TOPOGRAPHY_KWARGS))
"""
function plot_topography(filepath::String; input_dir::String = pwd(), participant_selection = participants(), kwargs...)
    if endswith(filepath, ".jld2")
        data = read_data(filepath)
        isnothing(data) && @minimal_error "No data found in file: $filepath"
        return plot_topography(data; kwargs...)
    end
    files = _require_batch_files(filepath, input_dir, participant_selection)
    results = NamedTuple[]
    for file in files
        @info "Plotting: $file"
        data = read_data(joinpath(input_dir, file))
        isnothing(data) && continue
        push!(results, plot_topography(data; kwargs...))
    end
    return results
end

function plot_topography(
    dat::SingleDataFrameEeg;
    channel_selection = channels(),
    sample_selection = samples(),
    interval_selection::Interval = times(),
    baseline_interval::Interval = nothing,
    n_topo::Union{Int,Nothing} = nothing,
    scrubber::Bool = false,
    display_plot = true,
    interactive = true,
    kwargs...,
)
    !isnothing(baseline_interval) && (dat = baseline(dat, baseline_interval))

    if scrubber
        dat_subset =
            subset(dat, channel_selection = channel_selection, sample_selection = sample_selection, interval_selection = interval_selection)
        if nrow(dat_subset.data) > 1
            return _plot_topography_scrubber(dat_subset; display_plot = display_plot, interactive = interactive, kwargs...)
        end
    end

    if !isnothing(n_topo) && n_topo > 1
        bin_datasets = _expand_topo_time_bins(
            dat,
            n_topo,
            channel_selection,
            sample_selection,
            interval_selection,
            get(kwargs, :time_unit, :s);
            include_cond_prefix = false,
        )
        return plot_topography(
            bin_datasets;
            channel_selection = channels(),
            sample_selection = samples(),
            interval_selection = times(),
            display_plot = display_plot,
            interactive = interactive,
            kwargs...,
        )
    end

    _set_window_title(_generate_window_title(dat))
    plot_kwargs = _merge_plot_kwargs(PLOT_TOPOGRAPHY_KWARGS, kwargs)
    fontsize_kw = isnothing(plot_kwargs[:theme_fontsize]) ? (;) : (; fontsize = plot_kwargs[:theme_fontsize])
    fig = _eegfun_figure(; figure_padding = plot_kwargs[:figure_padding], fontsize_kw...)
    ax = Axis(fig[1, 1], aspect = DataAspect())

    plot_topography!(
        fig,
        ax,
        dat;
        channel_selection = channel_selection,
        sample_selection = sample_selection,
        interval_selection = interval_selection,
        kwargs...,
    )

    if interactive
        deregister_interaction!(ax, :rectanglezoom)
        _setup_topo_keyboard_handlers!(fig, ax; zoom_step = plot_kwargs[:zoom_step])
        (dat isa ErpData || dat isa EpochData) && _setup_topo_selection!(fig, ax, dat, interval_selection)
    end

    display_plot && _display_figure(fig)
    _set_window_title("Makie")
    return (fig = fig, axes = [ax])
end

"""
    plot_topography!(fig, ax, dat::SingleDataFrameEeg; kwargs...)
    plot_topography!(fig, ax, dat::MultiDataFrameEeg, epoch::Int; kwargs...)

In-place topographic plot: subset `dat` and render onto an existing `Figure` and `Axis`.
"""
function plot_topography!(
    fig,
    ax,
    dat::SingleDataFrameEeg;
    channel_selection = channels(),
    sample_selection = samples(),
    interval_selection::Interval = times(),
    data_interp = nothing,
    x_bounds = nothing,
    y_bounds = nothing,
    kwargs...,
)
    plot_kwargs = _merge_plot_kwargs(PLOT_TOPOGRAPHY_KWARGS, kwargs)
    dat_subset =
        subset(dat, channel_selection = channel_selection, sample_selection = sample_selection, interval_selection = interval_selection)
    return _plot_topography!(
        fig,
        ax,
        dat_subset.data,
        dat_subset.layout;
        data_interp = data_interp,
        x_bounds = x_bounds,
        y_bounds = y_bounds,
        plot_kwargs...,
    )
end

function plot_topography(
    dat::AbstractVector{<:SingleDataFrameEeg};
    channel_selection = channels(),
    sample_selection = samples(),
    interval_selection::Interval = times(),
    baseline_interval::Interval = nothing,
    display_plot = true,
    interactive = true,
    kwargs...,
)
    !isnothing(baseline_interval) && (dat = baseline.(dat, Ref(baseline_interval)))
    n_datasets = length(dat)
    n_datasets == 0 && @minimal_error "Cannot plot empty vector of datasets"

    kwargs_dict = _kwargs_dict(kwargs)
    colorbar_enabled = get(kwargs_dict, :colorbar_plot, true)
    user_colorbar_position = get(kwargs_dict, :colorbar_position, nothing)
    colorbar_plot_numbers = get(kwargs_dict, :colorbar_plot_numbers, [])
    dims = get(kwargs_dict, :dims, nothing)
    n_topo = pop!(kwargs_dict, :n_topo, nothing)

    if !isnothing(n_topo) && n_topo > 1
        time_unit = get(kwargs_dict, :time_unit, :s)
        expanded = eltype(dat)[]
        sizehint!(expanded, n_datasets * n_topo)
        for dataset in dat
            append!(
                expanded,
                _expand_topo_time_bins(
                    dataset,
                    n_topo,
                    channel_selection,
                    sample_selection,
                    interval_selection,
                    time_unit;
                    include_cond_prefix = true,
                ),
            )
        end

        isnothing(dims) && (kwargs_dict[:dims] = (n_datasets, n_topo))
        if !haskey(kwargs_dict, :colorbar_plot)
            if haskey(kwargs_dict, :colorbar_plot_numbers) && !isempty(kwargs_dict[:colorbar_plot_numbers])
                kwargs_dict[:colorbar_plot_numbers] = [i * n_topo for i in kwargs_dict[:colorbar_plot_numbers]]
            else
                kwargs_dict[:colorbar_plot] = false
            end
        end

        return plot_topography(
            expanded;
            channel_selection = channels(),
            sample_selection = samples(),
            interval_selection = times(),
            display_plot = display_plot,
            interactive = interactive,
            kwargs_dict...,
        )
    end

    plot_rows, plot_cols = if isnothing(dims)
        _best_rect(n_datasets)
    else
        (length(dims) != 2 || any(dims .<= 0)) && @minimal_error "Invalid dimensions: $dims. Expected [rows, cols] with positive values."
        r, c = dims
        r * c < n_datasets && @minimal_error "Grid dimensions $dims provide $(r * c) cells but need $n_datasets."
        r, c
    end

    total_cols = if colorbar_enabled && !isnothing(user_colorbar_position) && user_colorbar_position[1] > 1
        plot_cols
    elseif colorbar_enabled
        plot_cols * 2
    else
        plot_cols
    end

    _set_window_title(_generate_window_title(dat))
    plot_kwargs = _merge_plot_kwargs(PLOT_TOPOGRAPHY_KWARGS, kwargs_dict)
    fontsize_kw = isnothing(plot_kwargs[:theme_fontsize]) ? (;) : (; fontsize = plot_kwargs[:theme_fontsize])
    fig = _eegfun_figure(; figure_padding = plot_kwargs[:figure_padding], fontsize_kw...)
    axes = Axis[]

    subsets = [
        subset(
            dataset;
            channel_selection = channel_selection,
            sample_selection = sample_selection,
            interval_selection = interval_selection,
        ) for dataset in dat
    ]

    first_layout = subsets[1].layout
    method = plot_kwargs[:method]
    gridscale = plot_kwargs[:gridscale]
    can_batch =
        n_datasets > 1 &&
        has_valid_coordinates(first_layout) &&
        all(s -> has_valid_coordinates(s.layout) && s.layout.data.label == first_layout.data.label, subsets) &&
        _topo_method_supported(method)

    batched_interp = nothing
    xb_shared = nothing
    yb_shared = nothing

    if can_batch
        _ensure_coordinates_2d!(first_layout)
        labels = first_layout.data.label
        channel_mat = Matrix{Float64}(undef, length(labels), n_datasets)
        for j = 1:n_datasets
            df = subsets[j].data
            for (i, lbl) in enumerate(labels)
                channel_mat[i, j] = mean(df[!, lbl])
            end
        end
        batched_interp, xb_shared, yb_shared = _topo_interpolate(channel_mat, first_layout, gridscale; method = method)
    end

    for (idx, dataset) in enumerate(subsets)
        base_row = div(idx - 1, plot_cols) + 1
        base_col = mod1(idx, plot_cols)

        if colorbar_enabled && !isnothing(user_colorbar_position)
            cb_row_offset, cb_col_offset = user_colorbar_position
            plot_row = cb_row_offset > 1 ? ((base_row - 1) * 2 + 1) : base_row
            plot_col = cb_row_offset > 1 ? base_col : ((base_col - 1) * 2 + 1)
            colorbar_row = plot_row + (cb_row_offset - 1)
            colorbar_col = plot_col + (cb_col_offset - 1)
        elseif colorbar_enabled
            plot_row = base_row
            plot_col = (base_col - 1) * 2 + 1
            colorbar_row = plot_row
            colorbar_col = plot_col + 1
        else
            plot_row = base_row
            plot_col = base_col
            colorbar_row = plot_row
            colorbar_col = plot_col + 1
        end

        ax = Axis(fig[plot_row, plot_col], aspect = DataAspect())
        push!(axes, ax)
        if hasproperty(dataset, :condition_name) && dataset.condition_name !== ""
            ax.title = dataset.condition_name
        end

        should_show_colorbar = colorbar_enabled && (isempty(colorbar_plot_numbers) || idx in colorbar_plot_numbers)
        _plot_topography!(
            fig,
            ax,
            dataset.data,
            dataset.layout;
            data_interp = can_batch ? batched_interp[idx] : nothing,
            x_bounds = can_batch ? xb_shared : nothing,
            y_bounds = can_batch ? yb_shared : nothing,
            kwargs_dict...,
            colorbar_plot = should_show_colorbar,
            colorbar_position = (colorbar_row, colorbar_col),
        )
    end

    if colorbar_enabled && (isnothing(user_colorbar_position) || user_colorbar_position[1] <= 1)
        for col = 1:total_cols
            colsize!(fig.layout, col, isodd(col) ? Auto() : Fixed(50))
        end
    end

    if interactive
        deregister_interaction!.(axes, :rectanglezoom)
        _setup_topo_keyboard_handlers!(fig, axes; zoom_step = plot_kwargs[:zoom_step])
        if all(d -> d isa ErpData || d isa EpochData, dat)
            _setup_shared_topo_selection!(fig, dat, TopoSelectionState(axes, interval_selection))
        end
    end

    display_plot && _display_figure(fig)
    _set_window_title("Makie")
    return (fig = fig, axes = axes)
end

plot_topography(dat::AbstractVector{<:EpochData}; kwargs...) = plot_topography(average_epochs.(dat); kwargs...)

function plot_topography(dat::AbstractVector{<:EpochData}, epoch::Int; kwargs...)
    @info "Plotting epoch $(epoch) for each dataset in the vector"
    return [plot_topography(d, epoch; kwargs...) for d in dat]
end

function plot_topography!(
    fig,
    ax,
    dat::MultiDataFrameEeg,
    epoch::Int;
    channel_selection = channels(),
    sample_selection = samples(),
    interval_selection::Interval = times(),
    kwargs...,
)
    return plot_topography!(
        fig,
        ax,
        epoch_to_continuous(dat, epoch);
        channel_selection = channel_selection,
        sample_selection = sample_selection,
        interval_selection = interval_selection,
        kwargs...,
    )
end

function plot_topography(
    dat::MultiDataFrameEeg,
    epoch::Int;
    channel_selection = channels(),
    sample_selection = samples(),
    interval_selection::Interval = times(),
    scrubber::Bool = false,
    display_plot = true,
    kwargs...,
)
    cont_dat = epoch_to_continuous(dat, epoch)
    if scrubber
        return plot_topography(
            cont_dat;
            channel_selection = channel_selection,
            sample_selection = sample_selection,
            interval_selection = interval_selection,
            scrubber = true,
            display_plot = display_plot,
            kwargs...,
        )
    end

    _set_window_title(_generate_window_title(dat) * " - Epoch $epoch")
    plot_kwargs = _merge_plot_kwargs(PLOT_TOPOGRAPHY_KWARGS, kwargs)
    fontsize_kw = isnothing(plot_kwargs[:theme_fontsize]) ? (;) : (; fontsize = plot_kwargs[:theme_fontsize])
    fig = _eegfun_figure(; figure_padding = plot_kwargs[:figure_padding], fontsize_kw...)
    ax = Axis(fig[1, 1], aspect = DataAspect())

    plot_topography!(
        fig,
        ax,
        cont_dat;
        channel_selection = channel_selection,
        sample_selection = sample_selection,
        interval_selection = interval_selection,
        kwargs...,
    )

    display_plot && _display_figure(fig)
    _set_window_title("Makie")
    return (fig = fig, axes = [ax])
end

function _draw_smooth_circle_mask!(ax::Axis, x_bounds::Tuple, y_bounds::Tuple)
    x_radius = max(abs(x_bounds[1]), abs(x_bounds[2]))
    y_radius = max(abs(y_bounds[1]), abs(y_bounds[2]))
    inner_x = Float32(x_radius * 0.97)
    inner_y = Float32(y_radius * 0.97)
    outer_x = Float32(x_radius + 0.05)
    outer_y = Float32(y_radius + 0.05)

    n_points = length(_UNIT_CIRCLE_200)
    annulus = Vector{Point2f}(undef, 2 * n_points)
    @inbounds for i = 1:n_points
        ca, sa = _UNIT_CIRCLE_200[i]
        annulus[i] = Point2f(outer_x * ca, outer_y * sa)
        annulus[2*n_points-i+1] = Point2f(inner_x * ca, inner_y * sa)
    end

    if Makie.to_color(ax.backgroundcolor[]) == RGBAf(1, 1, 1, 1) &&
       Makie.to_color(Makie.root(ax.scene).backgroundcolor[]) != RGBAf(1, 1, 1, 1)
        ax.backgroundcolor[] = RGBAf(0, 0, 0, 0)
    end

    mask_color = lift(ax.backgroundcolor, Makie.root(ax.scene).backgroundcolor) do ax_bg_raw, root_bg_raw
        ax_bg = Makie.to_color(ax_bg_raw)
        root_bg = Makie.to_color(root_bg_raw)
        return (ax_bg.alpha > 0 && !(ax_bg == RGBAf(1, 1, 1, 1) && root_bg != RGBAf(1, 1, 1, 1))) ? ax_bg : root_bg
    end
    poly!(ax, annulus, color = mask_color, strokecolor = mask_color, strokewidth = 0)
    return nothing
end

# =============================================================================
# TOPOGRAPHY INTERPOLATION & OPERATOR CACHING
# =============================================================================

function _topo_plot_radius(layout::Layout)
    x2 = layout.data.x2
    y2 = layout.data.y2
    r_max = 0.0
    @inbounds for i in eachindex(x2, y2)
        r = hypot(Float64(x2[i]), Float64(y2[i]))
        r > r_max && (r_max = r)
    end
    return r_max * 1.05
end

struct SphericalSplineOperator
    W::Matrix{Float64}
    valid::Vector{Int}
    grid_scale::Int
    radius::Float64
end

const _SPLINE_CACHE = Dict{UInt,SphericalSplineOperator}()
const _SPLINE_CACHE_LOCK = ReentrantLock()

const _RBF_GRID_CACHE = Dict{Tuple{Int,Float64},Tuple{Vector{Int},Matrix{Float64}}}()
const _RBF_GRID_LOCK = ReentrantLock()

function _get_rbf_grid(grid_scale::Int, plot_radius::Float64)
    key = (grid_scale, plot_radius)
    lock(_RBF_GRID_LOCK) do
        length(_RBF_GRID_CACHE) > 32 && empty!(_RBF_GRID_CACHE)
        get!(_RBF_GRID_CACHE, key) do
            r = range(-plot_radius, plot_radius, length = grid_scale)
            valid_cells = Int[]
            sizehint!(valid_cells, grid_scale * grid_scale)
            @inbounds for j = 1:grid_scale, i = 1:grid_scale
                hypot(r[i], r[j]) <= plot_radius && push!(valid_cells, i + (j - 1) * grid_scale)
            end
            grid_points = Matrix{Float64}(undef, 2, length(valid_cells))
            @inbounds for (k, lin) in enumerate(valid_cells)
                grid_points[1, k] = r[mod1(lin, grid_scale)]
                grid_points[2, k] = r[cld(lin, grid_scale)]
            end
            return (valid_cells, grid_points)
        end
    end
end

function _build_spline_operator(layout::Layout, ch::Vector{Int}, gs::Int, m::Int, λ::Float64, nterms::Int)
    d = layout.data
    U, k = _spline_geometry(d, ch)
    n = length(ch)
    f = _g_factors(m, nterms)

    A = Matrix{Float64}(undef, n + 1, n + 1)
    A[1:n, 1:n] .= _g_from_cos!(U * U', f)
    @inbounds for i = 1:n
        A[i, i] += λ
        A[i, n+1] = 1.0
        A[n+1, i] = 1.0
    end
    A[n+1, n+1] = 0.0
    X = lu(A) \ [Matrix{Float64}(I, n, n); zeros(1, n)]

    R = _topo_plot_radius(layout)
    valid, grid_pts = _get_rbf_grid(gs, R)
    nv = length(valid)
    Ug = Matrix{Float64}(undef, nv, 3)
    @inbounds for r = 1:nv
        x = grid_pts[1, r]
        y = grid_pts[2, r]
        ρ = hypot(x, y)
        s, c = sincos(min(k * ρ, π))
        sx = ρ > 0 ? s / ρ : 0.0
        Ug[r, 1] = sx * x
        Ug[r, 2] = sx * y
        Ug[r, 3] = c
    end

    Gg = _g_from_cos!(Ug * U', f)
    W = Gg * view(X, 1:n, :)
    W .+= view(X, (n+1):(n+1), :)
    return SphericalSplineOperator(W, valid, gs, R)
end

function _get_spline_operator(layout::Layout, ch::Vector{Int}, gs::Int, m::Int, λ::Float64, nterms::Int)
    d = layout.data
    inc_vals = hasproperty(d, :inc) ? d.inc[ch] : nothing
    azi_vals = hasproperty(d, :azi) ? d.azi[ch] : nothing
    key = hash((d.label[ch], d.x2[ch], d.y2[ch], inc_vals, azi_vals, gs, m, λ, nterms))
    lock(_SPLINE_CACHE_LOCK) do
        length(_SPLINE_CACHE) > 16 && empty!(_SPLINE_CACHE)
        get!(() -> _build_spline_operator(layout, ch, gs, m, λ, nterms), _SPLINE_CACHE, key)
    end
end

function _scatter_to_grid(op::SphericalSplineOperator, v::AbstractVector)
    out = fill(NaN, op.grid_scale, op.grid_scale)
    @inbounds for (k, lin) in enumerate(op.valid)
        out[lin] = v[k]
    end
    return out
end

function _row_all_finite(dat::AbstractMatrix)
    n_rows, n_cols = size(dat)
    fin = Vector{Bool}(undef, n_rows)
    @inbounds for r = 1:n_rows
        ok = true
        for c = 1:n_cols
            if !isfinite(dat[r, c])
                ok = false
                break
            end
        end
        fin[r] = ok
    end
    return fin
end

function _data_interpolation_topo_spherical_spline(
    dat::AbstractVecOrMat{<:AbstractFloat},
    layout::Layout,
    grid_scale::Int;
    layout_idx::AbstractVector{Int} = 1:size(dat, 1),
    m::Int = 4,
    lambda::Float64 = 1e-5,
    n_legendre_terms::Int = 15,
)
    length(layout_idx) == size(dat, 1) || @minimal_error "Data has $(size(dat, 1)) channels but layout_idx has $(length(layout_idx))"
    _ensure_coordinates_2d!(layout)

    finite = dat isa AbstractVector ? isfinite.(dat) : _row_all_finite(dat)
    all_fin = all(finite)
    n_fin = all_fin ? length(finite) : count(finite)
    n_fin >= 3 || @minimal_error "Fewer than 3 finite channels for spherical spline"

    ch = all_fin ? collect(Int, layout_idx) : Int[layout_idx[i] for i in eachindex(layout_idx) if finite[i]]
    op = _get_spline_operator(layout, ch, grid_scale, m, lambda, n_legendre_terms)

    bounds = (-op.radius, op.radius)
    if dat isa AbstractVector
        d_fin = all_fin ? dat : dat[finite]
        return _scatter_to_grid(op, op.W * d_fin), bounds, bounds
    else
        d_fin = all_fin ? dat : dat[finite, :]
        V = op.W * d_fin
        return [_scatter_to_grid(op, view(V, :, c)) for c in axes(V, 2)], bounds, bounds
    end
end

function _run_rbf_topo(
    method_obj::M,
    dat::AbstractVecOrMat{<:AbstractFloat},
    pts::Matrix{Float64},
    finite::AbstractVector{Bool},
    all_fin::Bool,
    valid_cells::Vector{Int},
    grid_points::Matrix{Float64},
    grid_scale::Int,
    plot_radius::Float64,
) where {M}
    bounds = (-plot_radius, plot_radius)
    if dat isa AbstractVector
        d = (all_fin && dat isa Vector{Float64}) ? dat : Float64[dat[i] for i in eachindex(dat) if finite[i]]
        itp = ScatteredInterpolation.interpolate(method_obj, pts, d)
        eval_res = ScatteredInterpolation.evaluate(itp, grid_points)
        result = fill(NaN, grid_scale, grid_scale)
        @inbounds for (k, lin) in enumerate(valid_cells)
            result[lin] = eval_res[k]
        end
        return result, bounds, bounds
    else
        n_fin = size(pts, 2)
        n_comps = size(dat, 2)
        d = if all_fin && dat isa Matrix{Float64}
            dat
        else
            mat = Matrix{Float64}(undef, n_fin, n_comps)
            @inbounds for c = 1:n_comps
                r_out = 0
                for r = 1:size(dat, 1)
                    if finite[r]
                        r_out += 1
                        mat[r_out, c] = Float64(dat[r, c])
                    end
                end
            end
            mat
        end
        itp = ScatteredInterpolation.interpolate(method_obj, pts, d)
        eval_res = ScatteredInterpolation.evaluate(itp, grid_points)
        results = [fill(NaN, grid_scale, grid_scale) for _ = 1:n_comps]
        @inbounds for c = 1:n_comps
            res = results[c]
            for (k, lin) in enumerate(valid_cells)
                res[lin] = eval_res[k, c]
            end
        end
        return results, bounds, bounds
    end
end

function _data_interpolation_topo(
    dat::AbstractVecOrMat{<:AbstractFloat},
    layout::Layout,
    grid_scale::Int;
    method::Symbol = :thin_plate,
    layout_idx::AbstractVector{Int} = 1:size(dat, 1),
)
    length(layout_idx) == size(dat, 1) || @minimal_error "Data has $(size(dat, 1)) channels but layout_idx has $(length(layout_idx))"
    _ensure_coordinates_2d!(layout)

    finite = dat isa AbstractVector ? isfinite.(dat) : _row_all_finite(dat)
    all_fin = all(finite)
    n_fin = all_fin ? length(finite) : count(finite)
    n_fin >= 3 || @minimal_error "Fewer than 3 finite channels for topography interpolation"

    x2 = layout.data.x2
    y2 = layout.data.y2
    pts = Matrix{Float64}(undef, 2, n_fin)
    k_out = 0
    @inbounds for i in eachindex(layout_idx)
        if all_fin || finite[i]
            k_out += 1
            li = layout_idx[i]
            pts[1, k_out] = Float64(x2[li])
            pts[2, k_out] = Float64(y2[li])
        end
    end

    plot_radius = _topo_plot_radius(layout)
    valid_cells, grid_points = _get_rbf_grid(grid_scale, plot_radius)

    if method === :thin_plate
        return _run_rbf_topo(
            ScatteredInterpolation.ThinPlate(),
            dat,
            pts,
            finite,
            all_fin,
            valid_cells,
            grid_points,
            grid_scale,
            plot_radius,
        )
    elseif method === :multiquadratic
        return _run_rbf_topo(
            ScatteredInterpolation.Multiquadratic(),
            dat,
            pts,
            finite,
            all_fin,
            valid_cells,
            grid_points,
            grid_scale,
            plot_radius,
        )
    elseif method === :inverse_multiquadratic
        return _run_rbf_topo(
            ScatteredInterpolation.InverseMultiquadratic(),
            dat,
            pts,
            finite,
            all_fin,
            valid_cells,
            grid_points,
            grid_scale,
            plot_radius,
        )
    elseif method === :gaussian
        return _run_rbf_topo(
            ScatteredInterpolation.Gaussian(),
            dat,
            pts,
            finite,
            all_fin,
            valid_cells,
            grid_points,
            grid_scale,
            plot_radius,
        )
    elseif method === :inverse_quadratic
        return _run_rbf_topo(
            ScatteredInterpolation.InverseQuadratic(),
            dat,
            pts,
            finite,
            all_fin,
            valid_cells,
            grid_points,
            grid_scale,
            plot_radius,
        )
    elseif method === :polyharmonic
        return _run_rbf_topo(
            ScatteredInterpolation.Polyharmonic(3),
            dat,
            pts,
            finite,
            all_fin,
            valid_cells,
            grid_points,
            grid_scale,
            plot_radius,
        )
    elseif method === :shepard
        return _run_rbf_topo(ScatteredInterpolation.Shepard(), dat, pts, finite, all_fin, valid_cells, grid_points, grid_scale, plot_radius)
    elseif method === :nearest
        return _run_rbf_topo(
            ScatteredInterpolation.NearestNeighbor(),
            dat,
            pts,
            finite,
            all_fin,
            valid_cells,
            grid_points,
            grid_scale,
            plot_radius,
        )
    else
        @minimal_error "Unknown method: $method. Supported: $(keys(_TOPO_RBF_METHODS))"
    end
end

const _TOPO_RBF_METHODS = Dict{Symbol,Any}(
    :multiquadratic         => ScatteredInterpolation.Multiquadratic,
    :inverse_multiquadratic => ScatteredInterpolation.InverseMultiquadratic,
    :gaussian               => ScatteredInterpolation.Gaussian,
    :inverse_quadratic      => ScatteredInterpolation.InverseQuadratic,
    :thin_plate             => ScatteredInterpolation.ThinPlate,
    :polyharmonic           => () -> ScatteredInterpolation.Polyharmonic(3),
    :shepard                => ScatteredInterpolation.Shepard,
    :nearest                => ScatteredInterpolation.NearestNeighbor,
)

_topo_method_supported(m::Symbol) = m === :spherical_spline || haskey(_TOPO_RBF_METHODS, m)

function _topo_interpolate(
    dat::AbstractVecOrMat{<:AbstractFloat},
    layout::Layout,
    gridscale::Int;
    method::Symbol = :thin_plate,
    layout_idx::AbstractVector{Int} = 1:size(dat, 1),
    kwargs...,
)
    if method === :spherical_spline
        return _data_interpolation_topo_spherical_spline(dat, layout, gridscale; layout_idx = layout_idx, kwargs...)
    elseif haskey(_TOPO_RBF_METHODS, method)
        return _data_interpolation_topo(dat, layout, gridscale; method = method, layout_idx = layout_idx, kwargs...)
    else
        supported = [keys(_TOPO_RBF_METHODS)..., :spherical_spline]
        @minimal_error "Unknown interpolation method: $method. Supported: $supported"
    end
end

# === INTERACTIVITY & SELECTION ===

function _setup_topo_keyboard_handlers!(fig::Figure, axes; zoom_step::Float64 = 0.2)
    on(events(fig).keyboardbutton) do event
        if event.action == Keyboard.press
            if event.key == Keyboard.i
                _show_plot_help(:topography)
            elseif event.key == Keyboard.up
                axes isa Vector ? _topo_scale_up!.(axes, zoom_step) : _topo_scale_up!(axes, zoom_step)
            elseif event.key == Keyboard.down
                axes isa Vector ? _topo_scale_down!.(axes, zoom_step) : _topo_scale_down!(axes, zoom_step)
            end
        end
    end
end

mutable struct TopoSelectionState
    active::Observable{Bool}
    bounds::Observable{Tuple{Float64,Float64,Float64,Float64}}
    visible::Observable{Bool}
    rectangles::Vector{Vector{Makie.Poly}}
    bounds_list::Observable{Vector{Tuple{Float64,Float64,Float64,Float64}}}
    temp_rectangles::Vector{Union{Makie.Poly,Nothing}}
    selected_channels::Vector{Symbol}
    axes::Vector{Axis}
    time_interval::Interval

    function TopoSelectionState(axes::Vector{Axis}, time_interval::Interval = times())
        n_axes = length(axes)
        new(
            Observable(false),
            Observable((0.0, 0.0, 0.0, 0.0)),
            Observable(false),
            [Makie.Poly[] for _ = 1:n_axes],
            Observable(Tuple{Float64,Float64,Float64,Float64}[]),
            Union{Makie.Poly,Nothing}[nothing for _ = 1:n_axes],
            Symbol[],
            axes,
            time_interval,
        )
    end
end

function _scale_topo_levels!(ax::Axis, scale_factor::Float64)
    for plot in ax.scene.plots
        if plot isa Makie.Contourf
            current_levels = plot.levels[]
            if !isnothing(current_levels)
                level_min, level_max = extrema(current_levels)
                center = (level_min + level_max) / 2
                half_range = (level_max - level_min) * scale_factor / 2
                new_min, new_max = center - half_range, center + half_range
                plot.levels[] = range(new_min, new_max, length = length(current_levels))
                plot.colorrange[] = (new_min, new_max)
                break
            end
        end
    end
end

function _scale_topo_levels!(ax::LScene, scale_factor::Float64)
    for plot in ax.scene.plots
        if hasproperty(plot, :colorrange) && !isnothing(plot.colorrange[])
            level_min, level_max = plot.colorrange[]
            center = (level_min + level_max) / 2
            half_range = (level_max - level_min) * scale_factor / 2
            plot.colorrange[] = (center - half_range, center + half_range)
            break
        end
    end
end

_topo_scale_up!(ax, zoom_step::Float64 = 0.2) = _scale_topo_levels!(ax, 1.0 - zoom_step)
_topo_scale_down!(ax, zoom_step::Float64 = 0.2) = _scale_topo_levels!(ax, 1.0 / (1.0 - zoom_step))

function _create_rectangles_on_all_axes!(selection_state::TopoSelectionState, rect_points::Vector{Point2f}, is_temporary::Bool)
    for (idx, other_ax) in enumerate(selection_state.axes)
        rect = poly!(other_ax, rect_points, color = (:blue, 0.3), strokecolor = :black, strokewidth = 1, visible = true, overdraw = true)
        if is_temporary
            selection_state.temp_rectangles[idx] = rect
        else
            push!(selection_state.rectangles[idx], rect)
        end
    end
end

function _setup_topo_selection!(fig::Figure, ax::Axis, original_data, time_interval::Interval = times())
    _setup_shared_topo_selection!(fig, [original_data], TopoSelectionState([ax], time_interval))
end

function _setup_shared_topo_selection!(fig::Figure, datasets::Vector, shared_selection_state::TopoSelectionState)
    shift_pressed = Observable(false)

    on(events(fig).keyboardbutton) do event
        if event.key in (Keyboard.left_shift, Keyboard.right_shift)
            shift_pressed[] = event.action != Keyboard.release
        end
    end

    on(events(fig).mousebutton) do event
        mouse_pos = events(fig).mouseposition[]
        active_ax, active_dataset = _find_active_axis_with_dataset(shared_selection_state.axes, mouse_pos, datasets)
        isnothing(active_ax) && return

        if event.button == Mouse.left
            if event.action == Mouse.press
                shift_pressed[] ? _start_topo_selection!(active_ax, shared_selection_state) :
                _clear_all_topo_selections!(shared_selection_state)
            elseif event.action == Mouse.release && shared_selection_state.active[]
                _finish_topo_selection!(active_ax, shared_selection_state, active_dataset)
            end
        elseif event.button == Mouse.right && event.action == Mouse.press
            selected_channels = shared_selection_state.selected_channels
            if !isempty(selected_channels)
                datasets_erp = _convert_to_erp_for_plotting(datasets, shared_selection_state.time_interval)
                _show_topo_context_menu!(datasets_erp, selected_channels, shared_selection_state.time_interval)
            else
                @minimal_warning "No channels selected. Select channels with Shift+Left Click+Drag, then right-click to plot ERP"
            end
        end
    end

    on(events(fig).mouseposition) do pos
        if shared_selection_state.active[]
            active_ax = _find_active_axis(shared_selection_state.axes, pos)
            !isnothing(active_ax) && _update_topo_selection!(active_ax, shared_selection_state)
        end
    end
end

function _rect_points(x0::Real, y0::Real, x1::Real, y1::Real)
    return Point2f[Point2f(x0, y0), Point2f(x1, y0), Point2f(x1, y1), Point2f(x0, y1)]
end

function _start_topo_selection!(ax::Axis, selection_state::TopoSelectionState)
    selection_state.active[] = true
    selection_state.visible[] = true
    mouse_pos = mouseposition(ax)
    mx, my = Float64(mouse_pos[1]), Float64(mouse_pos[2])
    selection_state.bounds[] = (mx, my, mx, my)
    _create_rectangles_on_all_axes!(selection_state, fill(Point2f(mx, my), 4), true)
    _update_topo_selection!(ax, selection_state)
end

function _finish_topo_selection!(ax::Axis, selection_state::TopoSelectionState, original_data = nothing)
    mouse_pos = mouseposition(ax)
    start_x, start_y = selection_state.bounds[][1], selection_state.bounds[][2]
    x_min, x_max = minmax(start_x, Float64(mouse_pos[1]))
    y_min, y_max = minmax(start_y, Float64(mouse_pos[2]))
    final_bounds = (x_min, y_min, x_max, y_max)
    selection_state.bounds[] = final_bounds

    current_bounds = selection_state.bounds_list[]
    push!(current_bounds, final_bounds)
    selection_state.bounds_list[] = current_bounds

    _create_rectangles_on_all_axes!(selection_state, _rect_points(x_min, y_min, x_max, y_max), false)

    for (idx, temp_rect) in enumerate(selection_state.temp_rectangles)
        if !isnothing(temp_rect)
            delete!(temp_rect.parent, temp_rect)
            selection_state.temp_rectangles[idx] = nothing
        end
    end

    all_selected_electrodes = Symbol[]
    for (bx0, by0, bx1, by1) in selection_state.bounds_list[]
        append!(all_selected_electrodes, _find_electrodes_in_region(bx0, by0, bx1, by1, original_data))
    end
    unique_electrodes = unique(all_selected_electrodes)
    @info "$(length(selection_state.bounds_list[])) regions; Channels: $unique_electrodes"
    selection_state.selected_channels = unique_electrodes
    selection_state.active[] = false
end

function _update_topo_selection!(ax::Axis, selection_state::TopoSelectionState)
    selection_state.active[] || return
    axis_pos = mouseposition(ax)
    start_x, start_y = selection_state.bounds[][1], selection_state.bounds[][2]
    end_x, end_y = Float64(axis_pos[1]), Float64(axis_pos[2])
    selection_state.bounds[] = (start_x, start_y, end_x, end_y)
    pts = _rect_points(start_x, start_y, end_x, end_y)
    for temp_rect in selection_state.temp_rectangles
        !isnothing(temp_rect) && (temp_rect[1] = pts)
    end
end

function _clear_all_topo_selections!(selection_state::TopoSelectionState)
    for axis_rects in selection_state.rectangles
        for rect in axis_rects
            delete!(rect.parent, rect)
        end
        empty!(axis_rects)
    end
    selection_state.bounds_list[] = Tuple{Float64,Float64,Float64,Float64}[]
    empty!(selection_state.selected_channels)
    selection_state.active[] = false
    selection_state.visible[] = false
end

function _find_electrodes_in_region(x_min::Float64, y_min::Float64, x_max::Float64, y_max::Float64, original_data)
    d = original_data.layout.data
    labels = d.label
    x2 = d.x2
    y2 = d.y2
    out = Symbol[]
    @inbounds for i in eachindex(labels, x2, y2)
        if x_min <= x2[i] <= x_max && y_min <= y2[i] <= y_max
            push!(out, Symbol(labels[i]))
        end
    end
    return out
end

function _convert_to_erp_for_plotting(datasets, time_interval::Interval = times())
    datasets_vec = datasets isa Vector ? datasets : [datasets]
    return map(datasets_vec) do dat
        if dat isa ErpData
            dat
        elseif dat isa ContinuousData
            subset_data = subset(dat, interval_selection = time_interval)
            ErpData(
                subset_data.file,
                1,
                splitext(basename(subset_data.file))[1],
                subset_data.data,
                subset_data.layout,
                subset_data.sample_rate,
                subset_data.analysis_info,
                1,
            )
        else
            ErpData(
                dat.file,
                dat.condition,
                dat.condition_name,
                average_trials(dat),
                dat.layout,
                dat.sample_rate,
                dat.analysis_info,
                n_epochs(dat),
            )
        end
    end
end

function _show_topo_context_menu!(datasets, selected_channels::Vector{Symbol}, time_interval::Interval = times())
    erp_datasets = _convert_to_erp_for_plotting(datasets, time_interval)
    return _show_topo_context_menu!(length(erp_datasets) == 1 ? erp_datasets[1] : erp_datasets, selected_channels, time_interval)
end

function _show_topo_context_menu!(
    datasets::Union{ErpData,Vector{ErpData}},
    selected_channels::Vector{Symbol},
    time_interval::Interval = times(),
)
    datasets_vec = datasets isa Vector ? datasets : [datasets]
    has_multiple_conditions = length(datasets_vec) > 1
    has_multiple_channels = length(selected_channels) > 1

    plot_types, plot_configs = if has_multiple_conditions && has_multiple_channels
        (
            [
                "Separate channels, separate conditions",
                "Average channels, separate conditions",
                "Separate channels, average conditions",
                "Average channels, average conditions",
            ],
            [
                (false, false, "separate channels, separate conditions"),
                (false, true, "average channels, separate conditions"),
                (true, false, "separate channels, average conditions"),
                (true, true, "average channels, average conditions"),
            ],
        )
    elseif has_multiple_conditions
        (["Separate conditions", "Average conditions"], [(false, false, "separate conditions"), (true, false, "average conditions")])
    elseif has_multiple_channels
        (["Plot Individual Channels", "Plot Averaged Channels"], [(false, false, "individual channels"), (false, true, "averaged channels")])
    else
        (["Plot Channel"], [(false, false, "single channel")])
    end

    menu_fig = _eegfun_figure(size = (400, 200))
    for (idx, plot_type) in enumerate(plot_types)
        btn = Button(menu_fig[idx, 1], label = plot_type)
        on(btn.clicks) do _
            avg_conditions, avg_channels, msg = plot_configs[idx]
            data_to_plot = avg_conditions ? _average_conditions(datasets_vec) : datasets_vec
            @info "Plotting ERP: $msg: $selected_channels"
            plot_erp(data_to_plot; channel_selection = channels(selected_channels), average_channels = avg_channels)
        end
    end

    _display_in_screen(menu_fig; size = (400, 200))
end

function _setup_scrubber_controls!(
    fig::Figure,
    time_vec::AbstractVector,
    current_idx::Observable{Int},
    ylim_obs::Observable{Tuple{Float64,Float64}},
    interactive::Bool,
    update_sample!::Function,
    show_mean!::Function,
)
    n_samples = length(time_vec)
    t_min, t_max = extrema(time_vec)

    controls = fig[2, 1:2] = GridLayout(tellwidth = false)
    Label(controls[1, 1], "Time:", font = :bold)
    slider = Slider(controls[1, 2], range = 1:n_samples, startvalue = 1)
    Label(controls[1, 3], lift(idx -> @sprintf("%.3f s", time_vec[idx]), current_idx), font = :bold, width = 75)
    btn_prev = Button(controls[1, 4], label = "<", width = 36)
    btn_play = Button(controls[1, 5], label = "Play", width = 54)
    btn_next = Button(controls[1, 6], label = ">", width = 36)
    btn_mean = Button(controls[1, 7], label = "Mean", width = 54)

    is_playing = Observable(false)
    play_task = Ref{Union{Nothing,Task}}(nothing)

    _stop_play() = (is_playing[] = false; btn_play.label[] = "Play")

    function _start_play()
        is_playing[] = true
        btn_play.label[] = "Pause"
        play_task[] = @async begin
            t_wall_start = time()
            t_data_start = time_vec[current_idx[]]
            if t_data_start >= t_max - 0.002
                t_data_start = t_min
                set_close_to!(slider, 1)
            end
            while is_playing[] && to_value(events(fig).window_open)
                target_data_time = t_data_start + (time() - t_wall_start)
                if target_data_time > t_max
                    t_wall_start = time()
                    t_data_start = t_min
                    target_data_time = t_min
                end
                target_idx = clamp(searchsortedfirst(time_vec, target_data_time), 1, n_samples)
                target_idx != current_idx[] && set_close_to!(slider, target_idx)
                sleep(0.016)
            end
        end
    end

    _toggle_play() = is_playing[] ? _stop_play() : _start_play()

    on(slider.value) do idx
        idx != current_idx[] && update_sample!(idx)
    end
    on(btn_prev.clicks) do _
        _stop_play()
        set_close_to!(slider, max(1, current_idx[] - 1))
    end
    on(btn_next.clicks) do _
        _stop_play()
        set_close_to!(slider, min(n_samples, current_idx[] + 1))
    end
    on(btn_play.clicks) do _
        _toggle_play()
    end
    on(btn_mean.clicks) do _
        _stop_play()
        show_mean!()
    end
    on(events(fig).window_open) do is_open
        !is_open && _stop_play()
    end

    if interactive
        on(events(fig).keyboardbutton) do event
            if event.action == Keyboard.press || event.action == Keyboard.repeat
                if event.key == Keyboard.left
                    _stop_play()
                    set_close_to!(slider, max(1, current_idx[] - 1))
                elseif event.key == Keyboard.right
                    _stop_play()
                    set_close_to!(slider, min(n_samples, current_idx[] + 1))
                elseif event.key == Keyboard.space
                    _toggle_play()
                elseif event.key == Keyboard.up
                    new_clim = max(1.0, ylim_obs[][2] * 0.8)
                    ylim_obs[] = (-new_clim, new_clim)
                elseif event.key == Keyboard.down
                    new_clim = ylim_obs[][2] * 1.25
                    ylim_obs[] = (-new_clim, new_clim)
                end
            end
        end
    end
    return nothing
end

function _plot_topography_scrubber(dat::SingleDataFrameEeg; display_plot = true, interactive = true, kwargs...)
    layout = dat.layout
    if !has_valid_coordinates(layout)
        @minimal_warning "Cannot create topographic plot: layout has no spatial coordinates."
        return
    end
    _ensure_coordinates_2d!(layout)

    df = dat.data
    df_props = Set(propertynames(df))
    available_channels = Symbol[Symbol(ch) for ch in layout.data.label if Symbol(ch) in df_props]
    if isempty(available_channels)
        @minimal_warning "No layout-matched channels found in data."
        return
    end

    time_vec = df[!, :time]
    n_samples = length(time_vec)
    t_min, t_max = extrema(time_vec)

    plot_kwargs = _merge_plot_kwargs(PLOT_TOPOGRAPHY_KWARGS, kwargs)
    fontsize_kw = isnothing(plot_kwargs[:theme_fontsize]) ? (;) : (; fontsize = plot_kwargs[:theme_fontsize])
    gridscale = get(plot_kwargs, :gridscale, 65)
    method = get(plot_kwargs, :method, :thin_plate)
    colormap = get(plot_kwargs, :colormap, nothing)
    num_levels = get(plot_kwargs, :num_levels, 21)
    colorbar_plot = get(plot_kwargs, :colorbar_plot, true)
    colorbar_kwargs = pop!(plot_kwargs, :colorbar_kwargs, (;))
    user_colorrange = get(plot_kwargs, :colorrange, nothing)

    X = Matrix{Float64}(undef, length(available_channels), n_samples)
    for (i, ch) in enumerate(available_channels)
        X[i, :] .= df[!, ch]
    end
    layout_syms = Symbol.(layout.data.label)
    layout_idx = Int[findfirst(==(ch), layout_syms)::Int for ch in available_channels]

    ylim_obs = if !isnothing(user_colorrange)
        Observable{Tuple{Float64,Float64}}((Float64(user_colorrange[1]), Float64(user_colorrange[2])))
    else
        q99 = Statistics.quantile!(abs.(vec(copy(X))), 0.99)
        clim = max(round(q99, digits = 1), 1.0)
        Observable{Tuple{Float64,Float64}}((-clim, clim))
    end
    levels_obs = lift(lim -> range(lim[1], lim[2], length = num_levels), ylim_obs)

    init_mat, x_bounds, y_bounds = _topo_interpolate(@view(X[:, 1]), layout, gridscale; method = method, layout_idx = layout_idx)
    mat_obs = Observable(init_mat)

    fig = _eegfun_figure(size = (650, 720); figure_padding = plot_kwargs[:figure_padding], fontsize_kw...)
    _set_window_title("Topography Scrubber | $(_generate_window_title(dat))")

    current_idx = Observable{Int}(1)
    title_obs = Observable("Topography | Latency: $(@sprintf("%.3f s (%d ms)", time_vec[1], round(Int, time_vec[1] * 1000)))")
    ax = Axis(fig[1, 1], aspect = DataAspect(), title = title_obs, titlesize = 16)
    hidedecorations!(ax)
    hidespines!(ax)
    deregister_interaction!(ax, :rectanglezoom)
    deregister_interaction!(ax, :scrollzoom)

    co = contourf!(
        ax,
        range(x_bounds[1], x_bounds[2], length = gridscale),
        range(y_bounds[1], y_bounds[2], length = gridscale),
        mat_obs,
        levels = levels_obs,
        colormap = _resolve_theme_colormap(ax, colormap),
        extendlow = :auto,
        extendhigh = :auto,
        nan_color = :transparent,
    )

    _draw_smooth_circle_mask!(ax, x_bounds, y_bounds)
    plot_layout_2d!(fig, ax, layout; point_plot = true, label_plot = false)
    colorbar_plot && Colorbar(fig[1, 2], co; label = "μV", height = Relative(0.75), colorbar_kwargs...)

    _setup_scrubber_controls!(
        fig,
        time_vec,
        current_idx,
        ylim_obs,
        interactive,
        function (idx::Int)
            c_idx = clamp(idx, 1, n_samples)
            current_idx[] = c_idx
            t = time_vec[c_idx]
            title_obs[] = @sprintf("Topography | Latency: %.3f s (%d ms)", t, round(Int, t * 1000))
            new_mat, _, _ = _topo_interpolate(@view(X[:, c_idx]), layout, gridscale; method = method, layout_idx = layout_idx)
            mat_obs[] = new_mat
        end,
        function ()
            mean_vals = vec(mean(X, dims = 2))
            mean_mat, _, _ = _topo_interpolate(mean_vals, layout, gridscale; method = method, layout_idx = layout_idx)
            mat_obs[] = mean_mat
            title_obs[] = @sprintf("Topography | Mean over %.3f – %.3f s", t_min, t_max)
        end,
    )

    screen = display_plot ? _display_figure(fig) : nothing
    _set_window_title("Makie")
    return (fig = fig, axes = [ax], screen = screen)
end
