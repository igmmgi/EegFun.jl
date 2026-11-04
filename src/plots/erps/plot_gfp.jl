const PLOT_GFP_KWARGS = Dict{Symbol,Tuple{Any,String}}(
    # Display parameters
    :display_plot => (true, "Whether to display the plot"),

    # Axis limits and labels
    :xlim => (nothing, "X-axis limits as (min, max) tuple. If nothing, automatically determined"),
    :ylim => (nothing, "Y-axis limits as (min, max) tuple. If nothing, automatically determined"),
    :xlabel => ("Time (s)", "Label for x-axis"),
    :ylabel => (nothing, "Label for y-axis (auto-determined based on normalize flag)"),

    # Title
    :plot_title => (nothing, "Plot title"),
    :show_title => (true, "Whether to show the plot title"),
    :plot_title_fontsize => (16, "Font size for plot titles"),
    :plot_title_position => (
        nothing,
        "Relative (x, y) coordinates for the plot title (e.g., (0.5, 0.95)). If provided, the title is drawn inside the axis.",
    ),
    :plot_title_align => ((:center, :top), "Alignment of the inner plot title"),

    # Legend
    :legend => (true, "Whether to show the legend when plotting multiple conditions"),

    # Line styling
    :linewidth => (2, "Line width for GFP/dissimilarity traces"),
    :color => (nothing, "Color for traces"),
    :linestyle => (nothing, "Line style for traces"),

    # Plot configuration
    :show_erp_traces => (false, "Whether to show individual ERP channel traces in top panel"),
    :show_dissimilarity => (false, "Whether to include Global Dissimilarity panel"),

    # Grid
    :xgrid => (false, "Whether to show x-axis grid"),
    :ygrid => (false, "Whether to show y-axis grid"),
    :xminorgrid => (false, "Whether to show x-axis minor grid"),
    :yminorgrid => (false, "Whether to show y-axis minor grid"),

    # Origin lines
    :add_x_origin => (true, "Whether to add vertical line at x=0"),
)

function _is_gfp_df_normalized(df::DataFrame, normalize::Union{Nothing,Bool})
    !isnothing(normalize) && return normalize
    if hasproperty(df, :normalized) && !isempty(df.normalized)
        return Bool(first(df.normalized))
    end
    meta = DataFrames.metadata(df)
    if haskey(meta, "normalized")
        val = meta["normalized"]
        return val isa Bool ? val : lowercase(string(val)) == "true"
    end
    return false
end

function _copy_gfp_trace_segment!(
    x_vec::Vector{Float64},
    y_vec::Vector{Float64},
    pos::Int,
    time_vec::AbstractVector,
    ch_col::AbstractVector,
    n_samples::Int,
)
    @inbounds @simd for s = 1:n_samples
        x_vec[pos+s-1] = Float64(time_vec[s])
        y_vec[pos+s-1] = Float64(ch_col[s])
    end
    return nothing
end

function _concat_gfp_erp_traces(df::DataFrame, selected_channels)
    time_vec = df[!, :time]
    n_samples = length(time_vec)
    n_chs = length(selected_channels)
    total_len = n_chs * n_samples + (n_chs - 1)
    x_vec = Vector{Float64}(undef, total_len)
    y_vec = Vector{Float64}(undef, total_len)
    pos = 1
    for (ch_i, ch) in enumerate(selected_channels)
        col_idx = ch isa Int ? ch : DataFrames.columnindex(df, ch)
        _copy_gfp_trace_segment!(x_vec, y_vec, pos, time_vec, df[!, col_idx], n_samples)
        pos += n_samples
        if ch_i != n_chs
            x_vec[pos] = NaN
            y_vec[pos] = NaN
            pos += 1
        end
    end
    return x_vec, y_vec
end

function _plot_gfp_panels!(
    fig::Figure,
    axes::Vector{Axis},
    dfs::AbstractVector{<:DataFrame},
    is_normalized::Bool,
    show_dissimilarity::Bool,
    plot_kwargs::Dict{Symbol,Any};
    start_panel::Int = 1,
    has_erp_traces::Bool = false,
    add_cond_labels::Bool = length(dfs) > 1,
)
    ax_gfp = Axis(fig[start_panel, 1])
    push!(axes, ax_gfp)

    ylabel_gfp = !isnothing(plot_kwargs[:ylabel]) ? plot_kwargs[:ylabel] : (is_normalized ? "GFP (%)" : "GFP (μV)")
    resolved_color = _resolve_theme_linecolor(ax_gfp, plot_kwargs[:color])
    resolved_linestyle = isnothing(plot_kwargs[:linestyle]) ? :solid : plot_kwargs[:linestyle]
    colors = (length(dfs) == 1 || !isnothing(plot_kwargs[:color])) ? [resolved_color] : Makie.wong_colors()
    lw = plot_kwargs[:linewidth]

    for (i, df) in enumerate(dfs)
        c = colors[mod1(i, length(colors))]
        lines!(
            ax_gfp,
            df[!, :time],
            df[!, :gfp];
            color = c,
            linewidth = lw,
            linestyle = resolved_linestyle,
            label = add_cond_labels ? "Condition $i" : nothing,
        )
    end

    if plot_kwargs[:show_title]
        ax_gfp.title = (!has_erp_traces && !isnothing(plot_kwargs[:plot_title])) ? plot_kwargs[:plot_title] : "Global Field Power"
    end

    _set_axis_properties!(
        ax_gfp;
        xlim = plot_kwargs[:xlim],
        ylim = plot_kwargs[:ylim],
        xlabel = show_dissimilarity ? "" : plot_kwargs[:xlabel],
        ylabel = ylabel_gfp,
    )
    _set_axis_grid!(
        ax_gfp;
        xgrid = plot_kwargs[:xgrid],
        ygrid = plot_kwargs[:ygrid],
        xminorgrid = plot_kwargs[:xminorgrid],
        yminorgrid = plot_kwargs[:yminorgrid],
    )
    plot_kwargs[:add_x_origin] && vlines!(ax_gfp, [0.0], color = :black, linewidth = 1, linestyle = :dash)
    length(dfs) > 1 && plot_kwargs[:legend] && axislegend(ax_gfp, position = :rt)

    if show_dissimilarity
        ax_diss = Axis(fig[start_panel+1, 1])
        push!(axes, ax_diss)
        for (i, df) in enumerate(dfs)
            c = colors[mod1(i, length(colors))]
            lines!(ax_diss, df[!, :time], df[!, :dissimilarity]; color = c, linewidth = lw, linestyle = resolved_linestyle)
        end
        plot_kwargs[:show_title] && (ax_diss.title = "Global Dissimilarity")
        _set_axis_properties!(
            ax_diss;
            xlim = plot_kwargs[:xlim],
            ylim = plot_kwargs[:ylim],
            xlabel = plot_kwargs[:xlabel],
            ylabel = is_normalized ? "Dissimilarity (%)" : "Dissimilarity",
        )
        _set_axis_grid!(
            ax_diss;
            xgrid = plot_kwargs[:xgrid],
            ygrid = plot_kwargs[:ygrid],
            xminorgrid = plot_kwargs[:xminorgrid],
            yminorgrid = plot_kwargs[:yminorgrid],
        )
        plot_kwargs[:add_x_origin] && vlines!(ax_diss, [0.0], color = :black, linewidth = 1, linestyle = :dash)
    end
    return nothing
end

"""
    plot_gfp(dat::ErpData; channel_selection = channels(), normalize = true, kwargs...)
    plot_gfp(datasets::AbstractVector{<:ErpData}; condition_selection = conditions(), channel_selection = channels(), normalize = true, kwargs...)
    plot_gfp(gfp_data::Union{DataFrame,AbstractVector{<:DataFrame}}; normalize = nothing, kwargs...)
    plot_gfp(filepath::String; input_dir = pwd(), participant_selection = participants(), kwargs...)

Plot Global Field Power (GFP) and optionally ERP channel traces (`show_erp_traces = true`)
and Global Dissimilarity (`show_dissimilarity = true`).

# Arguments
- `dat` / `datasets` / `gfp_data` / `filepath`: `ErpData`, vector of `ErpData`, precomputed GFP `DataFrame`(s), or `.jld2` file/pattern.
- `channel_selection`: Channel predicate for selecting channels (default: `channels()`).
- `condition_selection`: Condition filter for `Vector{ErpData}` (default: `conditions()`).
- `normalize`: Normalize GFP to 0–100% (`true` by default for `ErpData`; inferred from metadata when `nothing` for `DataFrame`s).

$(_generate_kwargs_doc(PLOT_GFP_KWARGS))
"""
function plot_gfp(filepath::String; input_dir::String = pwd(), participant_selection = participants(), kwargs...)
    if endswith(filepath, ".jld2")
        data = read_data(filepath)
        isnothing(data) && @minimal_error "No data found in file: $filepath"
        return plot_gfp(data; kwargs...)
    end
    files = _require_batch_files(filepath, input_dir, participant_selection)
    results = NamedTuple[]
    for file in files
        @info "Plotting: $file"
        data = read_data(joinpath(input_dir, file))
        isnothing(data) && continue
        push!(results, plot_gfp(data; kwargs...))
    end
    return results
end

function plot_gfp(dat::ErpData; channel_selection = channels(), normalize::Bool = true, kwargs...)
    return plot_gfp([dat]; channel_selection = channel_selection, normalize = normalize, kwargs...)
end

function plot_gfp(
    datasets::AbstractVector{<:ErpData};
    condition_selection = conditions(),
    channel_selection = channels(),
    normalize::Bool = true,
    kwargs...,
)
    datasets_filtered = datasets[get_selected_conditions(datasets, condition_selection)]
    _set_window_title(_generate_window_title(datasets_filtered))
    plot_kwargs = _merge_plot_kwargs(PLOT_GFP_KWARGS, kwargs)

    show_dissimilarity = plot_kwargs[:show_dissimilarity]
    show_erp_traces = plot_kwargs[:show_erp_traces]
    results = if show_dissimilarity
        [gfp_and_dissimilarity(dat; channel_selection = channel_selection, normalize = normalize) for dat in datasets_filtered]
    else
        [gfp(dat; channel_selection = channel_selection, normalize = normalize) for dat in datasets_filtered]
    end

    n_panels = 1 + Int(show_erp_traces) + Int(show_dissimilarity)
    fig = _eegfun_figure(size = (800, 200 * n_panels))
    axes = Axis[]
    panel_idx = 1

    if show_erp_traces
        ax_erp = Axis(fig[panel_idx, 1])
        push!(axes, ax_erp)
        for dat in datasets_filtered
            selected_channels = get_selected_channels(dat, channel_selection, include_meta = false, include_extra = false)
            isempty(selected_channels) && continue
            x_vec, y_vec = _concat_gfp_erp_traces(dat.data, selected_channels)
            lines!(ax_erp, x_vec, y_vec, color = (:gray, 0.5), linewidth = 0.5)
        end
        ax_erp.xlabel = ""
        ax_erp.ylabel = "Amplitude (μV)"
        plot_kwargs[:show_title] && (ax_erp.title = isnothing(plot_kwargs[:plot_title]) ? "EEG Channels" : plot_kwargs[:plot_title])
        plot_kwargs[:add_x_origin] && vlines!(ax_erp, [0.0], color = :black, linewidth = 1, linestyle = :dash)
        panel_idx += 1
    end

    _plot_gfp_panels!(
        fig,
        axes,
        results,
        normalize,
        show_dissimilarity,
        plot_kwargs;
        start_panel = panel_idx,
        has_erp_traces = show_erp_traces,
        add_cond_labels = length(datasets_filtered) > 1,
    )

    plot_kwargs[:display_plot] && _display_figure(fig)
    _set_window_title("Makie")
    return (fig = fig, axes = axes)
end

function plot_gfp(gfp_data::DataFrame; normalize::Union{Nothing,Bool} = nothing, kwargs...)
    plot_kwargs = _merge_plot_kwargs(PLOT_GFP_KWARGS, kwargs)
    hasproperty(gfp_data, :time) || @minimal_error("DataFrame must have :time column")
    hasproperty(gfp_data, :gfp) || @minimal_error("DataFrame must have :gfp column")

    show_dissimilarity = plot_kwargs[:show_dissimilarity] && hasproperty(gfp_data, :dissimilarity)
    fig = _eegfun_figure(size = (800, 300 * (1 + Int(show_dissimilarity))))
    axes = Axis[]

    _plot_gfp_panels!(
        fig,
        axes,
        [gfp_data],
        _is_gfp_df_normalized(gfp_data, normalize),
        show_dissimilarity,
        plot_kwargs;
        add_cond_labels = false,
    )

    _set_window_title("GFP")
    plot_kwargs[:display_plot] && _display_figure(fig)
    _set_window_title("Makie")
    return (fig = fig, axes = axes)
end

function plot_gfp(gfp_data::AbstractVector{<:DataFrame}; normalize::Union{Nothing,Bool} = nothing, kwargs...)
    plot_kwargs = _merge_plot_kwargs(PLOT_GFP_KWARGS, kwargs)
    for (i, df) in enumerate(gfp_data)
        hasproperty(df, :time) || @minimal_error("DataFrame $i must have :time column")
        hasproperty(df, :gfp) || @minimal_error("DataFrame $i must have :gfp column")
    end

    show_dissimilarity = plot_kwargs[:show_dissimilarity] && all(df -> hasproperty(df, :dissimilarity), gfp_data)
    fig = _eegfun_figure(size = (800, 300 * (1 + Int(show_dissimilarity))))
    axes = Axis[]
    is_normalized = isempty(gfp_data) ? false : _is_gfp_df_normalized(gfp_data[1], normalize)

    _plot_gfp_panels!(fig, axes, gfp_data, is_normalized, show_dissimilarity, plot_kwargs; add_cond_labels = true)

    _set_window_title("GFP")
    plot_kwargs[:display_plot] && _display_figure(fig)
    _set_window_title("Makie")
    return (fig = fig, axes = axes)
end
