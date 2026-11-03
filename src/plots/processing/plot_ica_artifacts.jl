"""
    _plot_component_metric!(ax, metrics_df, metric_col::Symbol, identified_comps::Vector{Int}; kwargs...)

Helper function to plot generic component metrics and highlight identified components.
"""
function _plot_component_metric!(
    ax,
    metrics_df,
    metric_col::Symbol,
    identified_comps::Vector{Int};
    threshold = nothing,
    threshold_style = :dash,
    threshold_color = :gray,
    threshold_label = nothing,
    add_legend = false,
    label_all = "All Components",
    label_highlight = "Identified Components",
)
    # Plot all components
    scatter!(ax, metrics_df.Component, metrics_df[!, metric_col], color = :gray, markersize = 12, label = add_legend ? label_all : nothing)

    # Add reference lines
    if !isnothing(threshold)
        t_label = add_legend && !isnothing(threshold_label) ? threshold_label : nothing
        hlines!(ax, threshold, color = threshold_color, linestyle = threshold_style, label = t_label)
    end

    # Highlight identified components
    if !isempty(identified_comps)
        highlight_df = metrics_df[in.(metrics_df.Component, Ref(identified_comps)), :]
        scatter!(
            ax,
            highlight_df.Component,
            highlight_df[!, metric_col],
            color = :red,
            markersize = 15,
            label = add_legend ? label_highlight : nothing,
        )

        # Add labels
        for row in eachrow(highlight_df)
            text!(
                ax,
                row.Component,
                row[metric_col],
                text = string(row.Component),
                color = :red,
                align = (:center, :bottom),
                fontsize = 14,
                offset = (0, 5),
            )
        end
    end

    if add_legend
        axislegend(ax, position = (1.0, 1.0))
    end
end

"""
    plot_eog_component_features(identified_comps::Dict, metrics_df::DataFrame; kwargs...)

Plot z-scores of EOG correlations from the metrics DataFrame and highlight identified components.

Uses the results from `identify_eye_components`.

# Arguments
- `identified_comps::Dict`: Dictionary returned by `identify_eye_components` (containing `:vEOG`, `:hEOG`).
- `metrics_df::DataFrame`: DataFrame returned by `identify_eye_components`. Expected to have columns `:vEOG_zscore`, `:hEOG_zscore`, and `:Component`.

# Keyword Arguments
- `z_threshold::Float64`: Z-score threshold to draw lines on the plot (default: 3.0).
- `display_plot::Bool`: Whether to display the plot (default: true).

# Returns
- Named tuple `(fig, axes)` where `axes = (ax_v, ax_h)`.
"""
function plot_eog_component_features(identified_comps::Dict, metrics_df::DataFrame; kwargs...)
    # Merge user kwargs with defaults
    plot_kwargs = _merge_plot_kwargs(PLOT_ICA_QUALITY_KWARGS, _normalize_plot_kwargs(kwargs))

    # Extract commonly used values
    z_threshold = plot_kwargs[:z_threshold]
    display_plot = plot_kwargs[:display_plot]
    show_title = plot_kwargs[:show_title]
    thresh_color = plot_kwargs[:zscore_threshold_color]
    thresh_style = plot_kwargs[:zscore_threshold_style]

    # Check if data is empty
    if isempty(metrics_df) || isempty(metrics_df.vEOG_zscore) || isempty(metrics_df.hEOG_zscore)
        @minimal_warning "Could not plot eye component features, input DataFrame or z-scores are empty."
        return (fig = _eegfun_figure(), axes = Axis[])
    end

    fig = _eegfun_figure()

    # Plot vEOG Correlation Z-Scores
    ax_v = Axis(fig[1, 1], xlabel = "Component Number", ylabel = "Z-Score", title = show_title ? "vEOG Correlation Z-Scores" : "")
    _plot_component_metric!(
        ax_v,
        metrics_df,
        :vEOG_zscore,
        identified_comps[:vEOG];
        threshold = [z_threshold, -z_threshold],
        threshold_color = thresh_color,
        threshold_style = thresh_style,
    )

    # Plot hEOG Correlation Z-Scores
    ax_h = Axis(fig[1, 2], xlabel = "Component Number", ylabel = "Z-Score", title = show_title ? "hEOG Correlation Z-Scores" : "")
    _plot_component_metric!(
        ax_h,
        metrics_df,
        :hEOG_zscore,
        identified_comps[:hEOG];
        threshold = [z_threshold, -z_threshold],
        threshold_color = thresh_color,
        threshold_style = thresh_style,
    )

    if display_plot
        _display_figure(fig)
    end

    return (fig = fig, axes = [ax_v, ax_h])
end

"""
    plot_spatial_kurtosis_components(kurtosis_comps::Vector{Int}, metrics_df::DataFrame; kwargs...)

Plot spatial kurtosis z-scores for all ICA components and highlight those exceeding the threshold.

# Arguments
- `kurtosis_comps::Vector{Int}`: Vector of component indices identified as having high spatial kurtosis.
- `metrics_df::DataFrame`: DataFrame containing spatial kurtosis metrics. Expected to have columns `:Component` and `:z_spatial_kurtosis`.

# Keyword Arguments
- `z_threshold::Float64`: Z-score threshold for the reference line (default: 3.0).
- `display_plot::Bool`: Whether to display the plot (default: true).

# Returns
- Named tuple `(fig = fig, axes = [ax])`.
"""
function plot_spatial_kurtosis_components(kurtosis_comps::Vector{Int}, metrics_df::DataFrame; kwargs...)
    # Merge user kwargs with defaults
    plot_kwargs = _merge_plot_kwargs(PLOT_ICA_QUALITY_KWARGS, _normalize_plot_kwargs(kwargs))

    # Extract commonly used values
    z_threshold = plot_kwargs[:z_threshold]
    display_plot = plot_kwargs[:display_plot]
    show_title = plot_kwargs[:show_title]
    default_title = "Component Spatial Kurtosis Z-Scores"
    title_str = !show_title ? "" : something(plot_kwargs[:plot_title], default_title)

    # Create figure
    fig = _eegfun_figure()

    # Plot spatial kurtosis z-scores
    ax = Axis(fig[1, 1], xlabel = "Component", ylabel = "Spatial Kurtosis Z-Score", title = title_str)

    _plot_component_metric!(
        ax,
        metrics_df,
        :z_spatial_kurtosis,
        kurtosis_comps;
        threshold = [z_threshold],
        threshold_color = plot_kwargs[:zscore_threshold_color],
        threshold_style = plot_kwargs[:zscore_threshold_style],
    )
    hlines!(ax, [0.0], color = :gray, linestyle = :dot)

    if display_plot
        _display_figure(fig)
    end

    return (fig = fig, axes = [ax])
end

"""
    plot_line_noise_components(line_noise_comps::Vector{Int}, metrics_df::DataFrame; kwargs...)

Plot spectral metrics used for line noise component identification.

# Arguments
- `line_noise_comps::Vector{Int}`: Vector of component indices identified as line noise.
- `metrics_df::DataFrame`: DataFrame containing line noise metrics. Expected to have columns `:Component`, `:power_ratio_zscore`, and `:harmonic_ratio`.

# Keyword Arguments
- `z_threshold::Float64`: Z-score threshold for the reference line (default: 3.0).
- `min_harmonic_power::Real`: Minimum harmonic power reference line (default: 0.5).
- `display_plot::Bool`: Whether to display the plot (default: true).

# Returns
- Named tuple `(fig = fig, axes = [ax1, ax2])`.
"""
function plot_line_noise_components(line_noise_comps::Vector{Int}, metrics_df::DataFrame; kwargs...)
    # Merge user kwargs with defaults
    plot_kwargs = _merge_plot_kwargs(PLOT_ICA_QUALITY_KWARGS, _normalize_plot_kwargs(kwargs))

    # Extract commonly used values
    z_threshold = plot_kwargs[:z_threshold]
    min_harmonic_power = plot_kwargs[:min_harmonic_power]
    display_plot = plot_kwargs[:display_plot]
    show_title = plot_kwargs[:show_title]
    show_legend = plot_kwargs[:legend]

    # Create figure with two subplots
    fig = _eegfun_figure(size = (1000, 400))

    # Plot 1: Power Ratio Z-Scores
    ax1 = Axis(
        fig[1, 1],
        xlabel = "Component",
        ylabel = "Power Ratio Z-Score",
        title = show_title ? "Line Frequency Power Ratio Z-Scores" : "",
    )
    _plot_component_metric!(
        ax1,
        metrics_df,
        :power_ratio_zscore,
        line_noise_comps;
        threshold = [z_threshold],
        threshold_color = plot_kwargs[:zscore_threshold_color],
        threshold_style = plot_kwargs[:zscore_threshold_style],
        threshold_label = "Threshold",
        add_legend = show_legend,
        label_highlight = "Line Noise Components",
    )

    # Plot 2: Harmonic Ratios
    ax2 = Axis(fig[1, 2], xlabel = "Component", ylabel = "Power Ratio", title = show_title ? "Harmonic Power Ratios" : "")
    _plot_component_metric!(
        ax2,
        metrics_df,
        :harmonic_ratio,
        line_noise_comps;
        threshold = [min_harmonic_power],
        threshold_color = plot_kwargs[:threshold_color],
        threshold_style = plot_kwargs[:threshold_style],
        threshold_label = "Min Harmonic Power",
        add_legend = show_legend,
        label_highlight = "Line Noise Components",
    )

    if display_plot
        _display_figure(fig)
    end

    return (fig = fig, axes = [ax1, ax2])
end

# ECG Component Identification 
"""
    _find_peaks(data::AbstractVector; min_prominence_std::Real=2.0, window_size::Int=1)

Enhanced peak finder. Finds indices where data point is greater than its
neighbors within a window (positive peaks) OR less than its neighbors within a window (negative peaks)
and exceeds the threshold. Returns indices of peaks.

# Arguments
- `data::AbstractVector`: Input data vector
- `min_prominence_std::Real`: Minimum prominence threshold in standard deviations (default: 2.0)
- `window_size::Int`: Number of samples to look on each side for comparison (default: 1)
"""
function _find_peaks(data::AbstractVector; min_prominence_std::Real = 2.0, window_size::Int = 1)
    if length(data) < 2 * window_size + 1
        return Int[]
    end

    peaks = Int[]
    mean_val = mean(data)
    std_val = std(data)
    # Handle zero std deviation case
    threshold_pos = (std_val ≈ 0) ? mean_val : mean_val + min_prominence_std * std_val
    threshold_neg = (std_val ≈ 0) ? mean_val : mean_val - min_prominence_std * std_val

    @inbounds for i = (window_size+1):(length(data)-window_size)
        val = data[i]

        # Positive peaks (greater than all neighbors in window and above threshold)
        if val > threshold_pos
            is_peak = true
            for j = (i-window_size):(i-1)
                if val <= data[j]
                    is_peak = false
                    break
                end
            end
            if is_peak
                for j = (i+1):(i+window_size)
                    if val <= data[j]
                        is_peak = false
                        break
                    end
                end
            end
            is_peak && push!(peaks, i)
            # Negative peaks (less than all neighbors in window and below threshold)
        elseif val < threshold_neg
            is_trough = true
            for j = (i-window_size):(i-1)
                if val >= data[j]
                    is_trough = false
                    break
                end
            end
            if is_trough
                for j = (i+1):(i+window_size)
                    if val >= data[j]
                        is_trough = false
                        break
                    end
                end
            end
            is_trough && push!(peaks, i)
        end
    end
    return peaks
end

"""
    plot_ecg_component_features(identified_comps::Vector{Int64}, metrics_df::DataFrame)

Create a simplified visualization of ECG component detection metrics.

# Arguments
- `identified_comps::Vector{Int64}`: Vector of component indices identified as ECG artifacts
- `metrics_df::DataFrame`: DataFrame with component metrics

# Returns
- `fig::Figure`: The Makie Figure containing the plot
"""
function plot_ecg_component_features(identified_comps::AbstractVector{<:Integer}, metrics_df::DataFrame; kwargs...)
    # Merge user kwargs with defaults
    plot_kwargs = _merge_plot_kwargs(PLOT_ICA_QUALITY_KWARGS, _normalize_plot_kwargs(kwargs))

    # Extract commonly used values
    display_plot = plot_kwargs[:display_plot]
    show_title = plot_kwargs[:show_title]
    show_legend = plot_kwargs[:legend]
    thresh_color = plot_kwargs[:zscore_threshold_color]
    thresh_style = plot_kwargs[:zscore_threshold_style]

    # Create figure with two panels
    fig = _eegfun_figure(size = (1000, 600))

    # Copy DataFrame so we never mutate caller's input
    df = copy(metrics_df)
    heart_rates = [isnan(ibi) || ibi <= 0 ? NaN : 60.0 / ibi for ibi in df.mean_ibi_s]
    df[!, :heart_rate_bpm] = heart_rates

    # Left panel: Heart Rate vs Peak Ratio
    ax1 =
        Axis(fig[1, 1], xlabel = "Heart Rate (BPM)", ylabel = "Peak Ratio (valid/total)", title = show_title ? "ECG Detection Metrics" : "")

    # Right panel: Heart Rate vs IBI Regularity (std)
    ax2 = Axis(fig[1, 2], xlabel = "BPM", ylabel = "IBI Std Dev (seconds)", title = show_title ? "Heart Rate Regularity" : "")

    # Plot non-ECG components
    non_ecg_idx = setdiff(1:nrow(df), identified_comps)
    non_ecg_df = df[non_ecg_idx, :]

    # Filter out NaNs for plotting
    valid_non_ecg = findall(.!isnan.(non_ecg_df.heart_rate_bpm) .& .!isnan.(non_ecg_df.peak_ratio))
    if !isempty(valid_non_ecg)
        scatter!(
            ax1,
            non_ecg_df.heart_rate_bpm[valid_non_ecg],
            non_ecg_df.peak_ratio[valid_non_ecg],
            color = :gray,
            markersize = 8,
            label = "Non-ECG",
        )
    end

    # Filter valid points for second plot
    valid_non_ecg2 = findall(.!isnan.(non_ecg_df.heart_rate_bpm) .& .!isnan.(non_ecg_df.std_ibi_s))
    if !isempty(valid_non_ecg2)
        scatter!(
            ax2,
            non_ecg_df.heart_rate_bpm[valid_non_ecg2],
            non_ecg_df.std_ibi_s[valid_non_ecg2],
            color = :gray,
            markersize = 8,
            label = "Non-ECG",
        )
    end

    # Plot ECG components
    ecg_df = df[in.(df.Component, Ref(identified_comps)), :]

    # Filter out NaNs
    valid_ecg = findall(.!isnan.(ecg_df.heart_rate_bpm) .& .!isnan.(ecg_df.peak_ratio))
    if !isempty(valid_ecg)
        scatter!(
            ax1,
            ecg_df.heart_rate_bpm[valid_ecg],
            ecg_df.peak_ratio[valid_ecg],
            color = :red,
            markersize = 12,
            marker = :diamond,
            label = "ECG",
        )

        # Add component labels
        for i in valid_ecg
            text!(
                ax1,
                ecg_df.heart_rate_bpm[i],
                ecg_df.peak_ratio[i],
                text = string(ecg_df.Component[i]),
                align = (:center, :bottom),
                offset = (0, 3),
                fontsize = 12,
            )
        end
    end

    # Plot ECG components in second panel
    valid_ecg2 = findall(.!isnan.(ecg_df.heart_rate_bpm) .& .!isnan.(ecg_df.std_ibi_s))
    if !isempty(valid_ecg2)
        scatter!(
            ax2,
            ecg_df.heart_rate_bpm[valid_ecg2],
            ecg_df.std_ibi_s[valid_ecg2],
            color = :red,
            markersize = 12,
            marker = :diamond,
            label = "ECG",
        )

        # Add component labels
        for i in valid_ecg2
            text!(
                ax2,
                ecg_df.heart_rate_bpm[i],
                ecg_df.std_ibi_s[i],
                text = string(ecg_df.Component[i]),
                align = (:center, :bottom),
                offset = (0, 3),
                fontsize = 12,
            )
        end
    end

    # Add reference ranges using pre-defined values
    # Normal heart rate range (typical values)
    vlines!(ax1, [60, 100], color = (:green, 0.5), linestyle = thresh_style, label = "Normal HR Range")
    vlines!(ax2, [60, 100], color = (:green, 0.5), linestyle = thresh_style)

    # Get threshold values from actual data when possible
    min_peak_ratio = plot_kwargs[:min_peak_ratio]
    max_std = 0.12

    if !isempty(identified_comps)
        ecg_identified = ecg_df[in.(ecg_df.Component, Ref(identified_comps)), :]
        if !isempty(ecg_identified)
            min_peak_ratio = minimum(ecg_identified.peak_ratio)
            max_std = maximum(ecg_identified.std_ibi_s)
        end
    end

    # Add threshold lines
    hlines!(ax1, [min_peak_ratio], color = (thresh_color, 0.5), linestyle = thresh_style, label = "Min Peak Ratio")
    hlines!(ax2, [max_std], color = (thresh_color, 0.5), linestyle = thresh_style, label = "Max StdDev")

    # Add legends
    if show_legend
        axislegend(ax1, position = :rt)
        axislegend(ax2, position = :rt)
    end

    # Calculate BPM range from actual data
    min_hr = 40   # Default minimum heart rate
    max_hr = 120  # Default maximum heart rate
    if !isempty(ecg_df) && any(.!isnan.(ecg_df.heart_rate_bpm))
        # Use the actual range from identified components
        valid_hrs = filter(!isnan, ecg_df.heart_rate_bpm)
        if !isempty(valid_hrs)
            min_hr = floor(Int, minimum(valid_hrs))
            max_hr = ceil(Int, maximum(valid_hrs))
        end
    end

    # Display summary text
    Label(
        fig[2, 1:2],
        "Found $(length(identified_comps)) ECG components: $(join(identified_comps, ", "))\n" *
        "Criteria: Heart rate $min_hr-$max_hr BPM, StdDev ≤ $(round(max_std, digits=3))s, Peak Ratio ≥ $(round(min_peak_ratio, digits=2))",
        fontsize = 14,
        tellwidth = false,
    )

    if display_plot
        _display_figure(fig)
    end

    return (fig = fig, axes = [ax1, ax2])
end

"""
    plot_artifact_components(ica::InfoIca, artifacts::ArtifactComponents; kwargs...)

Plot topoplots for all artifact components organized by type (vEOG, hEOG, ECG, Line Noise, Channel Noise).

This is a convenience function that takes the output of `combine_artifact_components` and creates
a comprehensive visualization showing all identified artifact components with clear labels.

# Arguments
- `ica::InfoIca`: The ICA result object
- `artifacts::ArtifactComponents`: The artifact components structure from `combine_artifact_components`

# Keyword Arguments
All keyword arguments from `plot_topography` are supported, including:
- `method::Symbol`: Interpolation method (see `plot_topography` for supported methods)
- `gridscale::Int`: Grid resolution for interpolation
- `colormap`: Colormap for the topography
- `display_plot::Bool`: Whether to display the plot

# Returns
- `Figure`: The Makie Figure containing all topoplots

# Examples
```julia
# Identify artifact components
eog_comps, _ = identify_eog_components(dat, ica)
ecg_comps, _ = identify_ecg_components(dat, ica)
line_noise_comps, _ = identify_line_noise_components(dat, ica)
channel_noise_comps, _ = identify_spatial_kurtosis_components(ica)

# Combine them
artifacts = combine_artifact_components(eog_comps, ecg_comps, line_noise_comps, channel_noise_comps)

# Plot all artifact components
fig = plot_artifact_components(ica, artifacts)
```
"""
function plot_artifact_components(ica::InfoIca, artifacts::ArtifactComponents; kwargs...)
    # ICA weights are dimensionless — override the default "μV" label
    norm_kw = _normalize_plot_kwargs(kwargs, _TOPO_KWARG_ALIASES)
    merged_kwargs = Dict{Symbol,Any}(:colorbar_label => "a.u.", pairs(norm_kw)...)

    # Merge user kwargs with defaults
    plot_kwargs = _merge_plot_kwargs(PLOT_TOPOGRAPHY_KWARGS, merged_kwargs)

    # Ensure coordinates are 2d and 3d
    _ensure_coordinates_2d!(ica.layout)
    _ensure_coordinates_3d!(ica.layout)

    # Extract commonly used kwargs
    method = pop!(plot_kwargs, :method)
    gridscale = pop!(plot_kwargs, :gridscale)
    colormap = pop!(plot_kwargs, :colormap)
    colorrange = pop!(plot_kwargs, :colorrange)
    display_plot = pop!(plot_kwargs, :display_plot)
    num_levels = pop!(plot_kwargs, :num_levels)
    show_title = pop!(plot_kwargs, :show_title)

    # Extract head shape parameters
    head_color = pop!(plot_kwargs, :head_color)
    head_linewidth = pop!(plot_kwargs, :head_linewidth)
    head_radius = pop!(plot_kwargs, :head_radius)

    # Extract electrode plotting parameters
    point_plot = pop!(plot_kwargs, :point_plot)
    label_plot = pop!(plot_kwargs, :label_plot)

    # Get all component types and their components
    component_data = [
        ("vEOG", artifacts.eog[:vEOG]),
        ("hEOG", artifacts.eog[:hEOG]),
        ("ECG", artifacts.ecg),
        ("Line Noise", artifacts.line_noise),
        ("Channel Noise", artifacts.channel_noise),
    ]

    # Filter out empty component lists
    component_types = [name for (name, comps) in component_data if !isempty(comps)]
    component_lists = [comps for (name, comps) in component_data if !isempty(comps)]

    # If no components, return empty figure
    if isempty(component_types)
        @minimal_warning "No artifact components found to plot"
        fig = _eegfun_figure()
        Label(fig[1, 1], "No artifact components found", fontsize = 16)
        return (fig = fig, axes = Axis[])
    end

    # Calculate total number of components to plot
    total_comps = sum(length(comps) for comps in component_lists)
    n_cols = min(4, total_comps)  # Max 4 columns
    n_rows = ceil(Int, total_comps / n_cols)

    # Create figure
    fig = _eegfun_figure(size = (n_cols * 200, n_rows * 200))

    # Plot each component individually
    axes = Axis[]
    plot_idx = 1
    for (comp_type, comps) in zip(component_types, component_lists)
        for comp_idx in comps
            row = ((plot_idx - 1) ÷ n_cols) + 1
            col = ((plot_idx - 1) % n_cols) + 1

            # Create subplot for this component
            ax = Axis(fig[row, col], aspect = DataAspect())
            ax.title = show_title ? "$comp_type $comp_idx" : ""
            ax.titlesize = 12
            push!(axes, ax)

            # Create topoplot data
            topo_data, x_bounds, y_bounds = _prepare_ica_topo_data(ica, comp_idx, method, gridscale)
            levels = if !isnothing(colorrange)
                _calculate_topo_levels(
                    topo_data;
                    use_global_scale = true,
                    global_min = colorrange[1],
                    global_max = colorrange[2],
                    num_levels = num_levels,
                )
            else
                _calculate_topo_levels(topo_data; num_levels = num_levels)
            end

            # Use the generic topo plotting function
            _plot_topo_on_axis!(
                ax,
                fig,
                topo_data,
                ica.layout,
                levels;
                x_bounds = x_bounds,
                y_bounds = y_bounds,
                gridscale = gridscale,
                colormap = colormap,
                head_color = head_color,
                head_linewidth = head_linewidth,
                head_radius = head_radius,
                point_plot = point_plot,
                label_plot = label_plot,
            )
            hidedecorations!(ax, grid = false)

            plot_idx += 1
        end
    end

    # Add overall title
    if show_title
        header_title = something(plot_kwargs[:plot_title], "ICA Artifact Components")
        Label(fig[0, :], header_title, fontsize = 18, font = :bold)
    end

    # Display plot if requested
    if display_plot
        _display_figure(fig)
    end

    return (fig = fig, axes = axes)
end
