function _plot_component_metric!(
    ax::Axis,
    metrics_df::DataFrame,
    metric_col::Symbol,
    identified_comps::AbstractVector{<:Integer};
    threshold = nothing,
    threshold_style = :dash,
    threshold_color = :gray,
    threshold_label = nothing,
    add_legend::Bool = false,
    label_all::String = "All Components",
    label_highlight::String = "Identified Components",
)
    comps_col = metrics_df[!, :Component]
    vals_col = metrics_df[!, metric_col]
    scatter!(ax, comps_col, vals_col, color = :gray, markersize = 12, label = add_legend ? label_all : nothing)

    if !isnothing(threshold)
        t_label = add_legend && !isnothing(threshold_label) ? threshold_label : nothing
        hlines!(ax, threshold, color = threshold_color, linestyle = threshold_style, label = t_label)
    end

    if !isempty(identified_comps)
        id_set = Set{Int}(identified_comps)
        hi_idx = findall(c -> Int(c) in id_set, comps_col)
        if !isempty(hi_idx)
            hi_comps = comps_col[hi_idx]
            hi_vals = vals_col[hi_idx]
            scatter!(ax, hi_comps, hi_vals, color = :red, markersize = 15, label = add_legend ? label_highlight : nothing)
            text!(ax, hi_comps, hi_vals; text = string.(hi_comps), color = :red, align = (:center, :bottom), fontsize = 14, offset = (0, 5))
        end
    end

    add_legend && axislegend(ax, position = (1.0, 1.0))
end

"""
    plot_eog_component_features(identified_comps::Dict, metrics_df::DataFrame; kwargs...) -> NamedTuple

Plot vEOG and hEOG correlation z-scores across ICA components and highlight identified ocular components.

$(_generate_kwargs_doc(PLOT_ICA_QUALITY_KWARGS))
"""
function plot_eog_component_features(identified_comps::Dict, metrics_df::DataFrame; kwargs...)
    plot_kwargs = _merge_plot_kwargs(PLOT_ICA_QUALITY_KWARGS, _normalize_plot_kwargs(kwargs))
    z_threshold = plot_kwargs[:z_threshold]
    show_title = plot_kwargs[:show_title]
    thresh_color = plot_kwargs[:zscore_threshold_color]
    thresh_style = plot_kwargs[:zscore_threshold_style]

    if isempty(metrics_df) || isempty(metrics_df.vEOG_zscore) || isempty(metrics_df.hEOG_zscore)
        @minimal_warning "Could not plot eye component features, input DataFrame or z-scores are empty."
        return (fig = _eegfun_figure(), axes = Axis[])
    end

    fig = _eegfun_figure()
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

    plot_kwargs[:display_plot] && _display_figure(fig)
    return (fig = fig, axes = [ax_v, ax_h])
end

"""
    plot_spatial_kurtosis_components(kurtosis_comps::AbstractVector{<:Integer}, metrics_df::DataFrame; kwargs...) -> NamedTuple

Plot spatial kurtosis z-scores across ICA components and highlight focal channel-noise components.

$(_generate_kwargs_doc(PLOT_ICA_QUALITY_KWARGS))
"""
function plot_spatial_kurtosis_components(kurtosis_comps::AbstractVector{<:Integer}, metrics_df::DataFrame; kwargs...)
    plot_kwargs = _merge_plot_kwargs(PLOT_ICA_QUALITY_KWARGS, _normalize_plot_kwargs(kwargs))
    title_str = !plot_kwargs[:show_title] ? "" : something(plot_kwargs[:plot_title], "Component Spatial Kurtosis Z-Scores")

    fig = _eegfun_figure()
    ax = Axis(fig[1, 1], xlabel = "Component", ylabel = "Spatial Kurtosis Z-Score", title = title_str)
    _plot_component_metric!(
        ax,
        metrics_df,
        :z_spatial_kurtosis,
        kurtosis_comps;
        threshold = [plot_kwargs[:z_threshold]],
        threshold_color = plot_kwargs[:zscore_threshold_color],
        threshold_style = plot_kwargs[:zscore_threshold_style],
    )
    hlines!(ax, [0.0], color = :gray, linestyle = :dot)

    plot_kwargs[:display_plot] && _display_figure(fig)
    return (fig = fig, axes = [ax])
end

"""
    plot_line_noise_components(line_noise_comps::AbstractVector{<:Integer}, metrics_df::DataFrame; kwargs...) -> NamedTuple

Plot line-frequency power ratio z-scores and harmonic power ratios across ICA components.

$(_generate_kwargs_doc(PLOT_ICA_QUALITY_KWARGS))
"""
function plot_line_noise_components(line_noise_comps::AbstractVector{<:Integer}, metrics_df::DataFrame; kwargs...)
    plot_kwargs = _merge_plot_kwargs(PLOT_ICA_QUALITY_KWARGS, _normalize_plot_kwargs(kwargs))
    show_title = plot_kwargs[:show_title]
    show_legend = plot_kwargs[:legend]

    fig = _eegfun_figure(size = (1000, 400))
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
        threshold = [plot_kwargs[:z_threshold]],
        threshold_color = plot_kwargs[:zscore_threshold_color],
        threshold_style = plot_kwargs[:zscore_threshold_style],
        threshold_label = "Threshold",
        add_legend = show_legend,
        label_highlight = "Line Noise Components",
    )

    ax2 = Axis(fig[1, 2], xlabel = "Component", ylabel = "Power Ratio", title = show_title ? "Harmonic Power Ratios" : "")
    _plot_component_metric!(
        ax2,
        metrics_df,
        :harmonic_ratio,
        line_noise_comps;
        threshold = [plot_kwargs[:min_harmonic_power]],
        threshold_color = plot_kwargs[:threshold_color],
        threshold_style = plot_kwargs[:threshold_style],
        threshold_label = "Min Harmonic Power",
        add_legend = show_legend,
        label_highlight = "Line Noise Components",
    )

    plot_kwargs[:display_plot] && _display_figure(fig)
    return (fig = fig, axes = [ax1, ax2])
end

function _scatter_ecg_panel!(ax::Axis, hr_bpm::Vector{Float64}, y_col::AbstractVector, comps_col::AbstractVector, is_ecg::Set{Int})
    non_x = Float64[]
    non_y = Float64[]
    ecg_x = Float64[]
    ecg_y = Float64[]
    ecg_lbls = String[]

    @inbounds for i in eachindex(hr_bpm, y_col, comps_col)
        x = hr_bpm[i]
        y = Float64(y_col[i])
        (isnan(x) || isnan(y)) && continue
        c = Int(comps_col[i])
        if c in is_ecg
            push!(ecg_x, x)
            push!(ecg_y, y)
            push!(ecg_lbls, string(c))
        else
            push!(non_x, x)
            push!(non_y, y)
        end
    end

    !isempty(non_x) && scatter!(ax, non_x, non_y, color = :gray, markersize = 8, label = "Non-ECG")
    if !isempty(ecg_x)
        scatter!(ax, ecg_x, ecg_y, color = :red, markersize = 12, marker = :diamond, label = "ECG")
        text!(ax, ecg_x, ecg_y; text = ecg_lbls, align = (:center, :bottom), offset = (0, 3), fontsize = 12)
    end
end

"""
    plot_ecg_component_features(identified_comps::AbstractVector{<:Integer}, metrics_df::DataFrame; kwargs...) -> NamedTuple

Plot heart-rate vs peak-ratio and inter-beat-interval regularity metrics for ECG ICA component identification.

$(_generate_kwargs_doc(PLOT_ICA_QUALITY_KWARGS))
"""
function plot_ecg_component_features(identified_comps::AbstractVector{<:Integer}, metrics_df::DataFrame; kwargs...)
    plot_kwargs = _merge_plot_kwargs(PLOT_ICA_QUALITY_KWARGS, _normalize_plot_kwargs(kwargs))
    show_title = plot_kwargs[:show_title]
    show_legend = plot_kwargs[:legend]
    thresh_color = plot_kwargs[:zscore_threshold_color]
    thresh_style = plot_kwargs[:zscore_threshold_style]

    fig = _eegfun_figure(size = (1000, 600))
    ax1 =
        Axis(fig[1, 1], xlabel = "Heart Rate (BPM)", ylabel = "Peak Ratio (valid/total)", title = show_title ? "ECG Detection Metrics" : "")
    ax2 = Axis(fig[1, 2], xlabel = "BPM", ylabel = "IBI Std Dev (seconds)", title = show_title ? "Heart Rate Regularity" : "")

    mean_ibi = metrics_df[!, :mean_ibi_s]
    std_ibi = metrics_df[!, :std_ibi_s]
    peak_ratio = metrics_df[!, :peak_ratio]
    comps_col = metrics_df[!, :Component]

    n = length(mean_ibi)
    hr_bpm = Vector{Float64}(undef, n)
    @inbounds for i = 1:n
        ibi = Float64(mean_ibi[i])
        hr_bpm[i] = (isnan(ibi) || ibi <= 0.0) ? NaN : 60.0 / ibi
    end

    is_ecg = Set{Int}(identified_comps)
    _scatter_ecg_panel!(ax1, hr_bpm, peak_ratio, comps_col, is_ecg)
    _scatter_ecg_panel!(ax2, hr_bpm, std_ibi, comps_col, is_ecg)

    vlines!(ax1, [60, 100], color = (:green, 0.5), linestyle = thresh_style, label = "Normal HR Range")
    vlines!(ax2, [60, 100], color = (:green, 0.5), linestyle = thresh_style)

    min_peak_ratio = Float64(plot_kwargs[:min_peak_ratio])
    max_std = 0.12
    min_hr = 40
    max_hr = 120

    if !isempty(is_ecg)
        pr_min = Inf
        std_max = -Inf
        hr_min = Inf
        hr_max = -Inf
        has_ecg_rows = false
        has_hr = false
        @inbounds for i = 1:n
            if Int(comps_col[i]) in is_ecg
                has_ecg_rows = true
                pr = Float64(peak_ratio[i])
                sd = Float64(std_ibi[i])
                pr < pr_min && (pr_min = pr)
                sd > std_max && (std_max = sd)
                hr = hr_bpm[i]
                if !isnan(hr)
                    hr < hr_min && (hr_min = hr)
                    hr > hr_max && (hr_max = hr)
                    has_hr = true
                end
            end
        end
        if has_ecg_rows
            min_peak_ratio = pr_min
            max_std = std_max
        end
        if has_hr
            min_hr = floor(Int, hr_min)
            max_hr = ceil(Int, hr_max)
        end
    end

    hlines!(ax1, [min_peak_ratio], color = (thresh_color, 0.5), linestyle = thresh_style, label = "Min Peak Ratio")
    hlines!(ax2, [max_std], color = (thresh_color, 0.5), linestyle = thresh_style, label = "Max StdDev")

    if show_legend
        axislegend(ax1, position = :rt)
        axislegend(ax2, position = :rt)
    end

    Label(
        fig[2, 1:2],
        "Found $(length(identified_comps)) ECG components: $(join(identified_comps, ", "))\n" *
        "Criteria: Heart rate $min_hr-$max_hr BPM, StdDev ≤ $(round(max_std, digits=3))s, Peak Ratio ≥ $(round(min_peak_ratio, digits=2))",
        fontsize = 14,
        tellwidth = false,
    )

    plot_kwargs[:display_plot] && _display_figure(fig)
    return (fig = fig, axes = [ax1, ax2])
end

"""
    plot_artifact_components(ica::InfoIca, artifacts::ArtifactComponents; kwargs...) -> NamedTuple

Plot scalp topographies for all identified ICA artifact components grouped by category (`vEOG`, `hEOG`, `ECG`, `Line Noise`, `Channel Noise`).

$(_generate_kwargs_doc(PLOT_TOPOGRAPHY_KWARGS))
"""
function plot_artifact_components(ica::InfoIca, artifacts::ArtifactComponents; kwargs...)
    norm_kw = _normalize_plot_kwargs(kwargs, _TOPO_KWARG_ALIASES)
    merged_kwargs = Dict{Symbol,Any}(:colorbar_label => "a.u.", pairs(norm_kw)...)
    plot_kwargs = _merge_plot_kwargs(PLOT_TOPOGRAPHY_KWARGS, merged_kwargs)

    _ensure_coordinates_2d!(ica.layout)
    _ensure_coordinates_3d!(ica.layout)

    method = plot_kwargs[:method]
    gridscale = plot_kwargs[:gridscale]
    colormap = plot_kwargs[:colormap]
    colorrange = plot_kwargs[:colorrange]
    display_plot = plot_kwargs[:display_plot]
    num_levels = plot_kwargs[:num_levels]
    show_title = plot_kwargs[:show_title]

    categories = (
        ("vEOG", artifacts.eog[:vEOG]),
        ("hEOG", artifacts.eog[:hEOG]),
        ("ECG", artifacts.ecg),
        ("Line Noise", artifacts.line_noise),
        ("Channel Noise", artifacts.channel_noise),
    )

    entries = Tuple{String,Int}[]
    for (name, comps) in categories, c in comps
        push!(entries, (name, c))
    end

    if isempty(entries)
        @minimal_warning "No artifact components found to plot"
        fig = _eegfun_figure()
        Label(fig[1, 1], "No artifact components found", fontsize = 16)
        return (fig = fig, axes = Axis[])
    end

    total_comps = length(entries)
    n_cols = min(4, total_comps)
    n_rows = cld(total_comps, n_cols)
    fig = _eegfun_figure(size = (n_cols * 200, n_rows * 200))

    unique_comps = unique([c for (_, c) in entries])
    all_topo_data, x_bounds, y_bounds = _prepare_ica_topo_data(ica, unique_comps, method, gridscale)
    comp_to_data = Dict{Int,Matrix{Float64}}(unique_comps[k] => all_topo_data[k] for k in eachindex(unique_comps))

    axes = Vector{Axis}(undef, total_comps)
    for (plot_idx, (comp_type, comp_idx)) in enumerate(entries)
        row, col = divrem(plot_idx - 1, n_cols) .+ (1, 1)
        ax = Axis(fig[row, col], aspect = DataAspect(), title = show_title ? "$comp_type $comp_idx" : "", titlesize = 12)
        axes[plot_idx] = ax

        topo_data = comp_to_data[comp_idx]
        levels =
            !isnothing(colorrange) ?
            _calculate_topo_levels(
                topo_data;
                use_global_scale = true,
                global_min = colorrange[1],
                global_max = colorrange[2],
                num_levels = num_levels,
            ) : _calculate_topo_levels(topo_data; num_levels = num_levels)

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
            head_color = plot_kwargs[:head_color],
            head_linewidth = plot_kwargs[:head_linewidth],
            head_radius = plot_kwargs[:head_radius],
            point_plot = plot_kwargs[:point_plot],
            label_plot = plot_kwargs[:label_plot],
        )
        hidedecorations!(ax, grid = false)
    end

    if show_title
        Label(fig[0, :], something(plot_kwargs[:plot_title], "ICA Artifact Components"), fontsize = 18, font = :bold)
    end

    display_plot && _display_figure(fig)
    return (fig = fig, axes = axes)
end
