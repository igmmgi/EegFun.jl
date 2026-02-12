"""
Topographic visualization of statistical test results across time windows.
Creates a grid of topographic maps showing t-statistics or difference amplitudes
with significant channels highlighted — inspired by FieldTrip's ft_clusterplot.
"""


"""
    plot_topo_stats(result::StatsResult;
                     n_topos::Int = 10,
                     time_range::Union{Nothing, Tuple{Real, Real}} = nothing,
                     topo_data::Symbol = :tvalues,
                     highlight_significant::Bool = true,
                     highlight_color = :white,
                     highlight_marker::Symbol = :circle,
                     highlight_size::Real = 8,
                     highlight_threshold::Real = 0.5,
                     kwargs...)

Plot a grid of topographic maps showing statistical results across time windows.

Each panel shows the spatial distribution of t-statistics (or difference amplitudes)
averaged within a time window, with significant channels highlighted as markers.

# Arguments
- `result::StatsResult`: Results from `analytic_test` or `permutation_test`
- `n_topos::Int`: Number of topographic panels (default: 10)
- `time_range::Union{Nothing, Tuple{Real, Real}}`: Time window to display in seconds.
  If `nothing`, uses the full time range from the result
- `topo_data::Symbol`: What to display on the topographic maps:
  - `:tvalues` (default) — t-statistics from the statistical test
  - `:difference` — difference wave amplitude (condition A − condition B)
- `highlight_significant::Bool`: Whether to overlay markers on significant channels (default: true)
- `highlight_color`: Color for significance markers (default: `:white`)
- `highlight_marker::Symbol`: Marker symbol for significant channels (default: `:circle`).
  Options: `:circle`, `:cross`, `:diamond`, `:star5`, `:xcross`, `:utriangle`, `:dtriangle`, etc.
- `highlight_size::Real`: Size of significance markers (default: 8)
- `highlight_threshold::Real`: Proportion of time points within a window that must be significant
  to mark a channel (default: 0.5). `0.0` = any single time point (union), `0.5` = majority,
  `1.0` = all time points (intersection)

Additional keyword arguments from PLOT_TOPOGRAPHY_KWARGS are supported:
$(_generate_kwargs_doc(PLOT_TOPOGRAPHY_KWARGS))

# Returns
- Named tuple `(fig, axes)`

# Examples
```julia
# Basic usage with analytic test
result = analytic_test(prepared, correction_method=:no)
plot_topo_stats(result)

# Show difference wave amplitudes instead of t-values
plot_topo_stats(result, topo_data=:difference)

# Focus on a specific time window with more panels
plot_topo_stats(result, time_range=(0.1, 0.4), n_topos=15)

# Custom marker style for significant channels
plot_topo_stats(result, highlight_marker=:star5, highlight_color=:yellow, highlight_size=12)
```
"""
function plot_topo_stats(
    result::StatsResult;
    n_topos::Int = 10,
    time_range::Union{Nothing,Tuple{Real,Real}} = nothing,
    topo_data::Symbol = :tvalues,
    highlight_significant::Bool = true,
    highlight_color = :white,
    highlight_marker::Symbol = :circle,
    highlight_size::Real = 8,
    highlight_threshold::Real = 0.5,
    kwargs...,
)
    # Validate topo_data
    topo_data in (:tvalues, :difference) || error("topo_data must be :tvalues or :difference, got :$topo_data")

    # Merge user kwargs with shared topography defaults
    plot_kwargs = _merge_plot_kwargs(PLOT_TOPOGRAPHY_KWARGS, kwargs)

    # Extract topography-specific parameters
    method = pop!(plot_kwargs, :method)
    gridscale = pop!(plot_kwargs, :gridscale)
    colormap = pop!(plot_kwargs, :colormap)
    ylim = pop!(plot_kwargs, :ylim)
    display_plot = pop!(plot_kwargs, :display_plot)
    figure_title = pop!(plot_kwargs, :figure_title)
    num_levels = pop!(plot_kwargs, :num_levels)

    # Pop ICA-specific keys not used here
    pop!(plot_kwargs, :interactive, nothing)
    pop!(plot_kwargs, :use_global_scale, nothing)
    pop!(plot_kwargs, :component_selection, nothing)

    # Override defaults for stats topo context
    plot_kwargs[:label_plot] = get(kwargs, :label_plot, false)   # Labels off by default for grid
    plot_kwargs[:point_plot] = get(kwargs, :point_plot, false)   # Points off by default for grid

    # Get time points and determine range
    all_time_points = result.time_points
    if isnothing(time_range)
        t_start, t_end = first(all_time_points), last(all_time_points)
    else
        t_start, t_end = time_range
    end

    # Find indices within the time range
    time_mask = t_start .<= all_time_points .<= t_end
    time_indices = findall(time_mask)
    if isempty(time_indices)
        error("No time points found in range ($t_start, $t_end). Data range: $(first(all_time_points)) to $(last(all_time_points))")
    end

    # Divide time indices into n_topos equal bins
    n_topos = min(n_topos, length(time_indices))
    bins = _partition_indices(time_indices, n_topos)

    # Get layout and electrodes
    layout = result.data[1].layout
    electrodes = result.electrodes

    # Validate layout
    if !has_valid_coordinates(layout)
        error("Cannot create topographic plot: layout has no spatial coordinates.")
    end

    # Ensure coordinates
    _ensure_coordinates_2d!(layout)
    _ensure_coordinates_3d!(layout)

    # Compute data for each bin
    topo_values = Vector{Vector{Float64}}(undef, n_topos)
    sig_masks = Vector{BitVector}(undef, n_topos)
    bin_time_labels = Vector{String}(undef, n_topos)

    for (i, bin_indices) in enumerate(bins)
        bin_time_start = all_time_points[first(bin_indices)]
        bin_time_end = all_time_points[last(bin_indices)]
        bin_time_labels[i] = @sprintf("%.0f – %.0f ms", bin_time_start * 1000, bin_time_end * 1000)

        if topo_data == :tvalues
            topo_values[i] = vec(mean(result.stat_matrix.t[:, bin_indices], dims = 2))
        else  # :difference
            cond_A = result.data[1]
            cond_B = result.data[2]
            diff_per_electrode = Float64[]
            for ch_sym in electrodes
                if ch_sym in propertynames(cond_A.data) && ch_sym in propertynames(cond_B.data)
                    a_vals = cond_A.data[bin_indices, ch_sym]
                    b_vals = cond_B.data[bin_indices, ch_sym]
                    push!(diff_per_electrode, mean(a_vals .- b_vals))
                else
                    push!(diff_per_electrode, 0.0)
                end
            end
            topo_values[i] = diff_per_electrode
        end

        # Significance: channel is significant if proportion of significant time points >= threshold
        n_bin_points = length(bin_indices)
        sig_combined = result.masks.positive[:, bin_indices] .| result.masks.negative[:, bin_indices]
        sig_proportion = vec(sum(sig_combined, dims = 2)) ./ n_bin_points
        sig_masks[i] = sig_proportion .>= highlight_threshold
    end

    # Compute symmetric color limits across all panels for consistent scaling
    all_vals = vcat(topo_values...)
    max_abs = maximum(abs, all_vals)
    if max_abs ≈ 0.0
        max_abs = 1.0
    end
    # Use user-provided ylim or compute symmetric limits
    if isnothing(ylim)
        ylim = (-max_abs, max_abs)
    end

    # Grid layout (dims comes from PLOT_TOPOGRAPHY_KWARGS)
    dims = pop!(plot_kwargs, :dims, nothing)
    if isnothing(dims)
        n_rows, n_cols = _best_rect(n_topos)
    else
        n_rows, n_cols = dims
        if n_rows * n_cols < n_topos
            error("Grid dimensions ($n_rows × $n_cols) provide $(n_rows * n_cols) cells but need $n_topos")
        end
    end

    # Create figure
    test_type = isa(result, AnalyticResult) ? "Analytic" : "Permutation"
    data_label = topo_data == :tvalues ? "t-statistic" : "Difference (μV)"
    fig_title = figure_title == "Topography Plot" ? "$test_type Test — $data_label" : figure_title

    fig = Figure(size = (200 * n_cols + 100, 200 * n_rows + 50))
    if plot_kwargs[:show_title]
        Label(fig[0, 1:n_cols], fig_title, fontsize = 18, font = :bold)
    end

    # Supported interpolation methods
    supported_methods =
        [:multiquadratic, :inverse_multiquadratic, :gaussian, :inverse_quadratic, :thin_plate, :polyharmonic, :shepard, :nearest]

    axes = Axis[]
    layout_labels = layout.data.label

    for i = 1:n_topos
        row = div(i - 1, n_cols) + 1
        col = mod1(i, n_cols)

        ax = Axis(fig[row, col], aspect = DataAspect(), title = bin_time_labels[i], titlesize = plot_kwargs[:title_fontsize])
        push!(axes, ax)

        # Map electrode values to layout order
        channel_data = topo_values[i]
        layout_values = Float64[]
        for lbl in layout_labels
            ch_idx = findfirst(==(lbl), electrodes)
            push!(layout_values, ch_idx !== nothing ? channel_data[ch_idx] : 0.0)
        end

        # Interpolate
        if method ∈ supported_methods
            data_interp, x_bounds, y_bounds = _data_interpolation_topo(layout_values, layout, gridscale; method = method)
        elseif method == :spherical_spline
            data_interp = _data_interpolation_topo_spherical_spline(layout_values, layout, gridscale)
            x_coords = layout.data.x2
            y_coords = layout.data.y2
            max_radius = maximum(sqrt.(x_coords .^ 2 .+ y_coords .^ 2))
            margin = max_radius * 0.05
            plot_radius = max_radius + margin
            x_bounds = (-plot_radius, plot_radius)
            y_bounds = (-plot_radius, plot_radius)
        else
            error("Unknown interpolation method: $method")
        end

        # Render topography
        co = contourf!(
            ax,
            range(x_bounds[1], x_bounds[2], length = gridscale),
            range(y_bounds[1], y_bounds[2], length = gridscale),
            data_interp,
            levels = range(ylim[1], ylim[2], gridscale);
            extendlow = :auto,
            extendhigh = :auto,
            colormap = colormap,
            nan_color = :transparent,
        )
        co.colorrange = ylim

        # Circle mask and head shape (uses shared kwargs)
        _draw_smooth_circle_mask!(ax, x_bounds, y_bounds)
        plot_layout_2d!(
            fig,
            ax,
            layout;
            point_plot = plot_kwargs[:point_plot],
            point_marker = plot_kwargs[:point_marker],
            point_markersize = plot_kwargs[:point_markersize],
            point_color = plot_kwargs[:point_color],
            label_plot = plot_kwargs[:label_plot],
            label_fontsize = plot_kwargs[:label_fontsize],
            label_color = plot_kwargs[:label_color],
            label_xoffset = plot_kwargs[:label_xoffset],
            label_yoffset = plot_kwargs[:label_yoffset],
            head_color = plot_kwargs[:head_color],
            head_linewidth = plot_kwargs[:head_linewidth],
            head_radius = plot_kwargs[:head_radius],
        )

        # Highlight significant channels
        if highlight_significant && any(sig_masks[i])
            sig_channel_indices = findall(sig_masks[i])
            sig_x = Float64[]
            sig_y = Float64[]
            for ch_idx in sig_channel_indices
                ch_sym = electrodes[ch_idx]
                layout_idx = findfirst(==(ch_sym), layout_labels)
                if layout_idx !== nothing
                    push!(sig_x, layout.data.x2[layout_idx])
                    push!(sig_y, layout.data.y2[layout_idx])
                end
            end
            if !isempty(sig_x)
                scatter!(
                    ax,
                    sig_x,
                    sig_y,
                    color = highlight_color,
                    marker = highlight_marker,
                    markersize = highlight_size,
                    strokewidth = 1,
                    strokecolor = :black,
                )
            end
        end

        hidedecorations!(ax)
    end

    # Add shared colorbar (uses shared colorbar kwargs)
    cb_label = topo_data == :tvalues ? "t-statistic" : "Difference (μV)"
    colorbar_kwargs = _extract_colorbar_kwargs!(plot_kwargs)
    # Remove keys we set explicitly to avoid conflicts
    pop!(colorbar_kwargs, :colorrange, nothing)
    pop!(colorbar_kwargs, :label, nothing)
    Colorbar(fig[1:n_rows, n_cols+1]; colorbar_kwargs..., colormap = colormap, colorrange = ylim, label = cb_label)

    if display_plot
        _display_figure(fig)
    end

    return (fig = fig, axes = axes)
end


"""
    _partition_indices(indices::Vector{Int}, n::Int) -> Vector{UnitRange{Int}}

Divide a vector of sorted indices into `n` approximately equal contiguous bins.
Returns vector of index ranges into the original `indices` vector.
"""
function _partition_indices(indices::Vector{Int}, n::Int)
    len = length(indices)
    bins = Vector{UnitRange{Int}}(undef, n)
    base_size = div(len, n)
    remainder = mod(len, n)

    start = 1
    for i = 1:n
        bin_size = base_size + (i <= remainder ? 1 : 0)
        stop = start + bin_size - 1
        bins[i] = indices[start]:indices[stop]
        start = stop + 1
    end

    return bins
end
