const PLOT_ERP_IMAGE_KWARGS = Dict{Symbol,Tuple{Any,String}}(
    # Display parameters
    :display_plot => (true, "Display the plot (true/false)"),
    :figure_title => ("", "Title drawn at the top of the entire figure canvas"),
    :figure_title_fontsize => (24, "Font size for figure title"),
    :interactive => (true, "Enable interactive features (true/false)"),
    :zoom_step => (0.2, "Fractional zoom step for arrow keys (e.g. 0.2 means 20% zoom in/out)"),
    :selection_color => (:blue, "Color for interactive selection rectangles"),
    :selection_alpha => (0.3, "Alpha (transparency) for interactive selection rectangles"),

    # Axis limits and labels
    :xlim => (nothing, "X-axis limits as (min, max) tuple. If nothing, automatically determined"),
    :ylim => (nothing, "Y-axis limits for ERP waveform as (min, max) tuple. If nothing, automatically determined"),
    :xlabel => ("Time (s)", "Label for x-axis"),
    :ylabel => ("Epoch", "Label for y-axis"),
    :yreversed => (false, "Whether to reverse the y-axis"),

    # Title
    :plot_title => (nothing, "Plot title"),
    :show_title => (true, "Show plot title (true/false)"),
    :plot_title_fontsize => (16, "Font size for plot titles"),
    :plot_title_position => (
        nothing,
        "Relative (x, y) coordinates for the plot title (e.g., (0.5, 0.95)). If provided, the title is drawn inside the axis.",
    ),
    :plot_title_align => ((:center, :top), "Alignment of the inner plot title"),

    # Image styling
    :colormap => (nothing, "Colormap for the image"),
    :colorrange => (nothing, "Color range for the image. If nothing, automatically determined"),
    :interpolate => (false, "Whether to interpolate the heatmap pixels for smoother appearance"),

    # ERP overlay
    :plot_erp => (true, "Whether to plot ERP average overlay"),
    :boxcar_average => (1, "Boxcar average window size for smoothing the ERP image (1 = no smoothing)"),

    # Colorbar parameters
    :colorbar_kwargs => (
        NamedTuple(),
        "Additional kwargs passed directly to the Makie Colorbar block (see Makie's [Colorbar documentation](https://docs.makie.org/stable/reference/blocks/colorbar/) for available attributes).",
    ),
    :colorbar_plot => (true, "Whether to show the colorbar"),
    :colorbar_position => (:right, "Position of the colorbar (:right, :left, :top, :bottom, or tuple)"),
    :colorbar_label => ("μV", "Label for the colorbar"),
    :colorbar_plot_numbers => ([], "Plot indices for which to show colorbars. Empty list shows colorbars for all plots."),

    # Grid
    :xgrid => (false, "Show x-axis grid (true/false)"),
    :ygrid => (false, "Show y-axis grid (true/false)"),
    :xminorgrid => (false, "Show x-axis minor grid (true/false)"),
    :yminorgrid => (false, "Show y-axis minor grid (true/false)"),

    # Origin lines
    :add_xy_origin => (true, "Add origin lines at x=0 and y=0 (true/false)"),

    # Layout parameters
    :layout_kwargs => (
        NamedTuple(),
        "Additional parameters to configure grid or topographical layouts (e.g., grid_dims, topo_plot_width, grid_rowgap)",
    ),
    :figure_padding => ((10, 30, 10, 10), "Padding around entire figure as (left, right, bottom, top) tuple (in pixels)"),
)

"""
    plot_erp_image(filepath::String; input_dir = pwd(), participant_selection = participants(), kwargs...)
    plot_erp_image(dat::EpochData;
                   layout::Union{Symbol, PlotLayout} = :single,
                   channel_selection = channels(),
                   channel_plot_order = nothing,
                   sample_selection = samples(),
                   interval_selection::Interval = times(),
                   kwargs...)
    plot_erp_image(data::AbstractVector{<:EpochData}; kwargs...)

Plot single-trial ERP image heatmaps across epochs for specified channels and samples.

# Arguments
- `dat::EpochData`: Epoch data structure (or `.jld2` filepath / filename pattern)
- `layout`: `:single` (default), `:grid`, `:topo`, or `PlotLayout`
- `channel_selection`: Channel selection predicate or symbol/vector (default: `channels()`)
- `sample_selection`: Sample selection predicate (default: `samples()`)
- `interval_selection`: Time interval filter (default: `times()`)

$(_generate_kwargs_doc(PLOT_ERP_IMAGE_KWARGS))

# Returns
- `NamedTuple` `(fig = fig, axes = axes)` (or vector of `NamedTuple`s for multi-condition/pattern inputs)
"""
function plot_erp_image(filepath::String; input_dir::String = pwd(), participant_selection = participants(), kwargs...)
    if endswith(filepath, ".jld2")
        data = read_data(filepath)
        isnothing(data) && @minimal_error "No data found in file: $filepath"
        return plot_erp_image(data; kwargs...)
    end

    files = _require_batch_files(filepath, input_dir, participant_selection)
    results = NamedTuple[]
    for file in files
        @info "Plotting: $file"
        data = read_data(joinpath(input_dir, file))
        isnothing(data) && continue
        push!(results, plot_erp_image(data; kwargs...))
    end
    return results
end

function plot_erp_image(
    dat::EpochData;
    layout::Union{Symbol,PlotLayout} = :single,
    channel_selection = channels(),
    channel_plot_order::Union{Nothing,AbstractVector} = nothing,
    sample_selection = samples(),
    interval_selection::Interval = times(),
    kwargs...,
)
    plot_kwargs = _merge_plot_kwargs(PLOT_ERP_IMAGE_KWARGS, kwargs)

    # Avoid deep-copying all EpochData DataFrames when no sample/interval subsetting is requested
    dat_subset = if _combine_interval_sample(interval_selection, sample_selection) isa _AllSamples
        dat
    else
        subset(dat; channel_selection = channel_selection, sample_selection = sample_selection, interval_selection = interval_selection)
    end

    _set_window_title(_generate_window_title(dat_subset))

    avail_channels = unique(vcat(channel_labels(dat_subset), extra_labels(dat_subset)))
    selected_channels = get_selected_channels(dat_subset, channel_selection; include_meta = false, include_extra = true)
    all_plot_channels = _order_plot_channels(selected_channels, avail_channels, channel_plot_order)
    isempty(all_plot_channels) && @minimal_error("No channels available for plotting after filtering")

    if isnothing(plot_kwargs[:plot_title]) && layout == :single
        plot_kwargs[:plot_title] = length(all_plot_channels) == 1 ? string(all_plot_channels[1]) : "$(_print_vector(all_plot_channels))"
    end

    colorbar_enabled = plot_kwargs[:colorbar_plot]
    colorbar_plot_numbers = plot_kwargs[:colorbar_plot_numbers]
    layout_kwargs = pop!(plot_kwargs, :layout_kwargs, (;))
    colorbar_kwargs = pop!(plot_kwargs, :colorbar_kwargs, (;))

    fig = _eegfun_figure()
    plot_layout = create_layout(layout, all_plot_channels, dat_subset.layout; layout_kwargs...)

    if layout == :grid && colorbar_enabled
        axes = Axis[]
        channels_vec = Symbol[]
        _, cols = plot_layout.dims
        for (idx, channel) in enumerate(plot_layout.channels)
            sg = GridLayout(fig[div(idx-1, cols)+1, mod1(idx, cols)])
            push!(axes, Axis(sg[1, 1]))
            push!(channels_vec, channel)
        end
    else
        axes, channels_vec = _apply_layout!(fig, plot_layout; plot_kwargs...)
    end

    if plot_kwargs[:interactive]
        for ax in axes
            deregister_interaction!(ax, :rectanglezoom)
        end
    end

    first_df = first(dat_subset.data)
    time_vec = first_df[!, :time]
    t_limits = !isnothing(plot_kwargs[:xlim]) ? plot_kwargs[:xlim] : (first(time_vec), last(time_vec))
    n_eps = n_epochs(dat_subset)
    epoch_range = 1:n_eps
    cb_pos = _get_colorbar_position(plot_kwargs[:colorbar_position], 1:1, 1:1)
    single_col_indices = plot_layout.type == :single ? [DataFrames.columnindex(first_df, ch) for ch in all_plot_channels] : Int[]

    single_raw_img = nothing
    heatmaps = Heatmap[]
    for (idx, (ax, channel)) in enumerate(zip(axes, channels_vec))
        img = if plot_layout.type == :single
            raw = _extract_erp_image_single(dat_subset.data, single_col_indices)
            single_raw_img = raw
            plot_kwargs[:boxcar_average] > 1 ? _apply_boxcar_average(raw, plot_kwargs[:boxcar_average]) : raw
        else
            raw = _extract_erp_image_channel(dat_subset.data, DataFrames.columnindex(first_df, channel))
            plot_kwargs[:boxcar_average] > 1 ? _apply_boxcar_average!(raw, plot_kwargs[:boxcar_average]) : raw
        end

        resolved_cmap = _resolve_theme_colormap(ax, plot_kwargs[:colormap])
        cr = if !isnothing(plot_kwargs[:colorrange])
            plot_kwargs[:colorrange]
        elseif _is_diverging_colormap(resolved_cmap)
            _symmetric_colorrange(extrema(img)...)
        else
            dmin, dmax = extrema(img)
            dmin == dmax ? (dmin - 1.0, dmax + 1.0) : (dmin, dmax)
        end

        hm = heatmap!(ax, time_vec, epoch_range, img, colormap = resolved_cmap, colorrange = cr, interpolate = plot_kwargs[:interpolate])
        push!(heatmaps, hm)

        if plot_layout.type == :single
            ax.xlabel = plot_kwargs[:plot_erp] ? "" : plot_kwargs[:xlabel]
            ax.ylabel = plot_kwargs[:ylabel]
            plot_kwargs[:plot_erp] && (ax.xticklabelsvisible = false)
        elseif plot_layout.type == :topo && plot_kwargs[:show_title]
            ax.title = string(channel)
        end

        ax.yreversed = plot_kwargs[:yreversed]
        xlims!(ax, t_limits)
        ylims!(ax, (1, n_eps))

        if colorbar_enabled
            if plot_layout.type == :single
                Colorbar(fig[cb_pos...], hm; colorbar_kwargs...)
            elseif plot_layout.type == :grid && (isempty(colorbar_plot_numbers) || idx in colorbar_plot_numbers)
                sg = ax.layoutobservables.gridcontent[].parent
                Colorbar(sg[cb_pos...], hm; colorbar_kwargs...)
            end
        end
    end

    # Apply axis properties to heatmap axes without waveform ylim (ylim applies only to ax_erp)
    heatmap_axis_kwargs = copy(plot_kwargs)
    heatmap_axis_kwargs[:ylim] = nothing
    _apply_axis_properties!.(axes; heatmap_axis_kwargs...)
    _apply_layout_axis_properties!(axes, plot_layout; heatmap_axis_kwargs...)

    for ax in axes
        ylims!(ax, (1, n_eps))
        if plot_layout.type == :topo
            ax.xticksvisible = true
            ax.yticksvisible = true
            ax.xticklabelsvisible = false
            ax.yticklabelsvisible = false
        end
        _set_origin_lines!(ax; add_xy_origin = plot_kwargs[:add_xy_origin])
    end

    if plot_layout.type == :topo
        scale_pos = get(plot_kwargs, :layout_topo_scale_pos, get(plot_layout.metadata, :topo_scale_pos, (0.8, -0.8)))
        scale_ax = Axis(
            fig[1, 1],
            width = Relative(plot_layout.metadata[:topo_plot_width]),
            height = Relative(plot_layout.metadata[:topo_plot_height]),
            halign = scale_pos[1],
            valign = scale_pos[2],
        )
        push!(axes, scale_ax)
        scale_ax.title = ""
        xlims!(scale_ax, t_limits)
        ylims!(scale_ax, (1, n_eps))
        _set_axis_properties!(
            scale_ax;
            xlim = t_limits,
            ylim = (1, n_eps),
            xlabel = plot_kwargs[:xlabel],
            ylabel = plot_kwargs[:ylabel],
            yreversed = plot_kwargs[:yreversed],
        )
        _set_axis_grid!(
            scale_ax;
            xgrid = plot_kwargs[:xgrid],
            ygrid = plot_kwargs[:ygrid],
            xminorgrid = plot_kwargs[:xminorgrid],
            yminorgrid = plot_kwargs[:yminorgrid],
        )
        scale_ax.bottomspinevisible = true
        scale_ax.topspinevisible = true
        scale_ax.leftspinevisible = true
        scale_ax.rightspinevisible = true
        _set_origin_lines!(scale_ax; add_xy_origin = get(plot_kwargs, :add_xy_origin, true))

        if colorbar_enabled && !isempty(heatmaps)
            Colorbar(
                fig[1, 1],
                heatmaps[1];
                label = plot_kwargs[:colorbar_label],
                halign = scale_pos[1] + 0.015,
                valign = scale_pos[2],
                tellheight = false,
                tellwidth = false,
                colorbar_kwargs...,
            )
        end
    end

    if plot_kwargs[:interactive]
        _setup_erp_image_interactivity!(fig, axes, heatmaps; zoom_step = plot_kwargs[:zoom_step])
        selection_state =
            SharedSelectionState(axes; selection_color = plot_kwargs[:selection_color], selection_alpha = plot_kwargs[:selection_alpha])
        _setup_unified_selection!(fig, axes, selection_state, dat_subset)
        if plot_layout.type in (:topo, :grid)
            _setup_channel_selection_events!(fig, selection_state, plot_layout, dat_subset, axes)
        end
    end

    if plot_kwargs[:plot_erp] && plot_layout.type == :single && !isnothing(single_raw_img)
        ax_erp = Axis(fig[2, 1])
        push!(axes, ax_erp)

        avg_data = _row_means(single_raw_img)
        lines!(ax_erp, time_vec, avg_data, color = :black, linewidth = 1)

        ax_erp.xlabel = plot_kwargs[:xlabel]
        ax_erp.ylabel = "Amplitude (μV)"
        ax_erp.xgridvisible = false
        ax_erp.ygridvisible = false
        ax_erp.xminorgridvisible = false
        ax_erp.yminorgridvisible = false
        ylims!(ax_erp, !isnothing(plot_kwargs[:ylim]) ? plot_kwargs[:ylim] : extrema(avg_data))

        vlines!(ax_erp, [0], color = :gray, linewidth = 0.5)
        hlines!(ax_erp, [0], color = :gray, linewidth = 0.5)
        xlims!(ax_erp, !isnothing(plot_kwargs[:xlim]) ? plot_kwargs[:xlim] : axes[1].xaxis.attributes.limits[])

        rowsize!(fig.layout, 1, Relative(2 / 3))
        rowsize!(fig.layout, 2, Relative(1 / 3))
    end

    length(axes) > 1 && linkxaxes!(axes...)
    plot_kwargs[:display_plot] && _display_figure(fig)

    _set_window_title("Makie")
    return (fig = fig, axes = axes)
end

function plot_erp_image(data::AbstractVector{<:EpochData}; kwargs...)
    results = plot_erp_image.(data; kwargs...)
    return length(results) == 1 ? only(results) : results
end

function _setup_erp_image_interactivity!(fig::Figure, axes::Vector{Axis}, heatmaps::Vector; zoom_step::Float64 = 0.2)
    _setup_keyboard_tracking(fig)
    active_axis = Ref(axes[1])

    on(events(fig).mouseposition) do mouse_pos
        for ax in axes
            vp = ax.scene.viewport[]
            if vp.origin[1] <= mouse_pos[1] <= vp.origin[1] + vp.widths[1] && vp.origin[2] <= mouse_pos[2] <= vp.origin[2] + vp.widths[2]
                active_axis[] = ax
                break
            end
        end
    end

    on(events(fig).keyboardbutton) do event
        if event.action == Keyboard.press && event.key == Keyboard.i
            _show_plot_help(:erp_image)
        elseif event.action in (Keyboard.press, Keyboard.repeat)
            action = if event.key == Keyboard.up
                :y_more
            elseif event.key == Keyboard.down
                :y_less
            elseif event.key == Keyboard.left
                :x_less
            elseif event.key == Keyboard.right
                :x_more
            else
                :none
            end
            action != :none && _handle_erp_image_navigation!(axes, heatmaps, action, active_axis[], zoom_step)
        end
    end
end

function _handle_erp_image_navigation!(axes::Vector{Axis}, heatmaps::Vector, action::Symbol, active_axis::Axis, zoom_step::Float64)
    if action in (:x_less, :x_more)
        (action == :x_less ? _xless! : _xmore!)(active_axis, zoom_step)
    elseif active_axis == axes[1]
        factor = action == :y_more ? (1.0 / (1.0 - zoom_step)) : (1.0 - zoom_step)
        for hm in heatmaps
            hm.colorrange[] = hm.colorrange[] .* factor
        end
    else
        (action == :y_more ? _ymore! : _yless!)(active_axis, zoom_step)
    end
end

function _row_means(mat::Matrix{Float64})
    n_samples, n_eps = size(mat)
    out = zeros(Float64, n_samples)
    inv_n = 1.0 / n_eps
    @inbounds for ep = 1:n_eps
        @simd for s = 1:n_samples
            out[s] += mat[s, ep]
        end
    end
    out .*= inv_n
    return out
end

# Smooths an (n_samples × n_epochs) image matrix across epochs (columns) using a running column sum
function _apply_boxcar_average(img::Matrix{Float64}, window_size::Int)
    return _apply_boxcar_average!(Matrix{Float64}(undef, size(img)...), img, window_size)
end

function _apply_boxcar_average!(img::Matrix{Float64}, window_size::Int)
    return _apply_boxcar_average!(Matrix{Float64}(undef, size(img)...), img, window_size)
end

function _apply_boxcar_average!(out::Matrix{Float64}, img::Matrix{Float64}, window_size::Int)
    n_samples, n_eps = size(img)
    half_window = window_size ÷ 2
    acc = zeros(Float64, n_samples)

    win_start = 1
    win_end = min(n_eps, 1 + half_window)
    @inbounds for ep = win_start:win_end
        @simd for s = 1:n_samples
            acc[s] += img[s, ep]
        end
    end
    inv_w = 1.0 / (win_end - win_start + 1)
    @inbounds @simd for s = 1:n_samples
        out[s, 1] = acc[s] * inv_w
    end

    @inbounds for ep = 2:n_eps
        new_start = max(1, ep - half_window)
        new_end = min(n_eps, ep + half_window)
        if new_start > win_start
            @simd for s = 1:n_samples
                acc[s] -= img[s, win_start]
            end
        end
        if new_end > win_end
            @simd for s = 1:n_samples
                acc[s] += img[s, new_end]
            end
        end
        win_start = new_start
        win_end = new_end
        inv_w = 1.0 / (win_end - win_start + 1)
        @simd for s = 1:n_samples
            out[s, ep] = acc[s] * inv_w
        end
    end
    return out
end

@inline function _copy_col_f64!(img::Matrix{Float64}, ep::Int, col::AbstractVector, n_samples::Int)
    if col isa Vector{Float64}
        @inbounds @simd for s = 1:n_samples
            img[s, ep] = col[s]
        end
    else
        @inbounds @simd for s = 1:n_samples
            img[s, ep] = Float64(col[s])
        end
    end
end

@inline function _add_col_f64!(img::Matrix{Float64}, ep::Int, col::AbstractVector, n_samples::Int)
    if col isa Vector{Float64}
        @inbounds @simd for s = 1:n_samples
            img[s, ep] += col[s]
        end
    else
        @inbounds @simd for s = 1:n_samples
            img[s, ep] += Float64(col[s])
        end
    end
end

# Returns (n_samples × n_epochs) directly for contiguous column-major extraction and zero-transpose heatmap!
function _extract_erp_image_single(epochs_data::Vector{DataFrame}, col_indices::Vector{Int})
    n_eps = length(epochs_data)
    n_samples = nrow(first(epochs_data))
    img = Matrix{Float64}(undef, n_samples, n_eps)
    n_cols = length(col_indices)
    inv_cols = 1.0 / n_cols

    for (ep, df) in enumerate(epochs_data)
        _copy_col_f64!(img, ep, df[!, col_indices[1]], n_samples)
        for c = 2:n_cols
            _add_col_f64!(img, ep, df[!, col_indices[c]], n_samples)
        end
        if n_cols > 1
            @inbounds @simd for s = 1:n_samples
                img[s, ep] *= inv_cols
            end
        end
    end
    return img
end

function _extract_erp_image_channel(epochs_data::Vector{DataFrame}, col_idx::Int)
    n_eps = length(epochs_data)
    n_samples = nrow(first(epochs_data))
    img = Matrix{Float64}(undef, n_samples, n_eps)
    for (ep, df) in enumerate(epochs_data)
        _copy_col_f64!(img, ep, df[!, col_idx], n_samples)
    end
    return img
end
