"""
Plotting functions for Global Field Power (GFP) and Global Dissimilarity.
"""

# === DEFAULT KEYWORD ARGUMENTS ===
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

"""
    plot_gfp(dat::ErpData;
             channel_selection::Function = channels(),
             normalize::Bool = true,
             kwargs...)

Plot Global Field Power (GFP) for ERP data.

# Arguments
- `dat::ErpData`: ERP data structure
- `channel_selection::Function`: Channel predicate for selecting channels (default: all channels)
- `normalize::Bool`: If true, normalize GFP to 0-100% (default: true)
- `kwargs`: Additional keyword arguments

$(_generate_kwargs_doc(PLOT_GFP_KWARGS))

# Examples
```julia
using EegFun, JLD2

# Load and plot GFP
erp_data = load("participant_1_erps.jld2", "erps")[1]
plot_gfp(erp_data)

# Plot with individual channel traces
plot_gfp(erp_data, show_erp_traces = true)

# Plot with Global Dissimilarity
plot_gfp(erp_data, show_dissimilarity = true)

# Plot all three panels
plot_gfp(erp_data, show_erp_traces = true, show_dissimilarity = true)

# Custom styling
plot_gfp(erp_data, 
         color = :blue,
         linewidth = 3,
         xlim = (-0.2, 0.8))
```

# Plot Layout
By default, shows only the GFP trace. With `show_erp_traces=true`, adds a top panel
with all channel traces. With `show_dissimilarity=true`, adds a bottom panel with
Global Dissimilarity.

    plot_gfp(filepath::String; input_dir=pwd(), participant_selection=participants(), kwargs...)

Load ERP data and plot GFP. Accepts either a `.jld2` filepath or a pattern
to discover and plot all matching files (one plot per file).

# Examples
```julia
plot_gfp("erps_final.jld2")
plot_gfp("erps_final")
```
"""
function plot_gfp(filepath::String; input_dir::String = pwd(), participant_selection = participants(), kwargs...)
    if endswith(filepath, ".jld2")
        data = read_data(filepath)
        isnothing(data) && @minimal_error "No data found in file: $filepath"
        return plot_gfp(data; kwargs...)
    else
        files = _find_batch_files(filepath, input_dir, participant_selection)
        isempty(files) && @minimal_error "No files matching pattern '$filepath' in $input_dir"

        results = NamedTuple[]
        for file in sort(files, by = _natural_sort_key)
            file_path = joinpath(input_dir, file)
            @info "Plotting: $file"
            data = read_data(file_path)
            isnothing(data) && continue
            result = plot_gfp(data; kwargs...)
            push!(results, result)
        end
        return results
    end
end

function _is_gfp_df_normalized(df::DataFrame, normalize::Union{Nothing,Bool})
    !isnothing(normalize) && return normalize
    if hasproperty(df, :normalized) && !isempty(df.normalized)
        return Bool(first(df.normalized))
    end
    meta = DataFrames.metadata(df)
    for key in ("normalized", :normalized)
        if haskey(meta, string(key))
            val = meta[string(key)]
            return val isa Bool ? val : lowercase(string(val)) == "true"
        end
    end
    return false
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

    # Apply condition_selection first
    datasets_filtered = datasets[get_selected_conditions(datasets, condition_selection)]

    # Generate window title from datasets
    title_str = _generate_window_title(datasets_filtered)
    _set_window_title(title_str)

    # Merge user kwargs and default kwargs
    plot_kwargs = _merge_plot_kwargs(PLOT_GFP_KWARGS, _normalize_plot_kwargs(kwargs))

    # Calculate GFP (and optionally dissimilarity) for all datasets
    show_dissimilarity = plot_kwargs[:show_dissimilarity]
    show_erp_traces = plot_kwargs[:show_erp_traces]

    if show_dissimilarity
        # Calculate both GFP and dissimilarity
        results = [gfp_and_dissimilarity(dat; channel_selection = channel_selection, normalize = normalize) for dat in datasets_filtered]
    else
        # Calculate only GFP
        results = [gfp(dat; channel_selection = channel_selection, normalize = normalize) for dat in datasets_filtered]
    end

    # Determine number of panels
    n_panels = 1 + (show_erp_traces ? 1 : 0) + (show_dissimilarity ? 1 : 0)

    # Create figure
    fig = _eegfun_figure(size = (800, 200 * n_panels))
    axes = Axis[]

    panel_idx = 1

    # Panel 1: ERP traces (if requested)
    if show_erp_traces
        ax_erp = Axis(fig[panel_idx, 1])
        push!(axes, ax_erp)

        # Vectorized background ERP traces: single NaN-separated lines! call per dataset
        for dat in datasets_filtered
            selected_channels = get_selected_channels(dat, channel_selection, include_meta = false, include_extra = false)
            isempty(selected_channels) && continue

            time_vec = dat.data.time
            n_samples = length(time_vec)
            n_chs = length(selected_channels)
            total_len = n_chs * n_samples + (n_chs - 1)

            x_vec = Vector{Float64}(undef, total_len)
            y_vec = Vector{Float64}(undef, total_len)

            pos = 1
            for (ch_i, ch) in enumerate(selected_channels)
                col_idx = ch isa Int ? ch : DataFrames.columnindex(dat.data, ch)
                ch_data = dat.data[!, col_idx]::Vector{Float64}
                @simd for s = 1:n_samples
                    x_vec[pos+s-1] = time_vec[s]
                    y_vec[pos+s-1] = ch_data[s]
                end
                pos += n_samples
                if ch_i != n_chs
                    x_vec[pos] = NaN
                    y_vec[pos] = NaN
                    pos += 1
                end
            end
            lines!(ax_erp, x_vec, y_vec, color = (:gray, 0.5), linewidth = 0.5)
        end

        ax_erp.xlabel = ""
        ax_erp.ylabel = "Amplitude (μV)"
        if plot_kwargs[:show_title]
            ax_erp.title = isnothing(plot_kwargs[:plot_title]) ? "EEG Channels" : plot_kwargs[:plot_title]
        end

        if plot_kwargs[:add_x_origin]
            vlines!(ax_erp, [0.0], color = :black, linewidth = 1, linestyle = :dash)
        end

        panel_idx += 1
    end

    # Panel 2 (or 1): GFP
    ax_gfp = Axis(fig[panel_idx, 1])
    push!(axes, ax_gfp)

    # Determine y-label
    ylabel_gfp = if !isnothing(plot_kwargs[:ylabel])
        plot_kwargs[:ylabel]
    else
        normalize ? "GFP (%)" : "GFP (μV)"
    end

    resolved_color = _resolve_theme_linecolor(ax_gfp, plot_kwargs[:color])
    resolved_linestyle = isnothing(plot_kwargs[:linestyle]) ? :solid : plot_kwargs[:linestyle]

    # Plot GFP for each dataset
    colors = length(datasets_filtered) == 1 ? [resolved_color] : (!isnothing(plot_kwargs[:color]) ? [resolved_color] : Makie.wong_colors())
    for (i, result) in enumerate(results)
        color = colors[mod1(i, length(colors))]
        lines!(
            ax_gfp,
            result.time,
            result.gfp,
            color = color,
            linewidth = plot_kwargs[:linewidth],
            linestyle = resolved_linestyle,
            label = length(datasets_filtered) > 1 ? "Condition $i" : nothing,
        )
    end

    # Apply styling and limits
    if plot_kwargs[:show_title]
        if !show_erp_traces
            ax_gfp.title = isnothing(plot_kwargs[:plot_title]) ? "Global Field Power" : plot_kwargs[:plot_title]
        else
            ax_gfp.title = "Global Field Power"
        end
    end

    # Set axis limits and labels using the shared function
    _set_axis_properties!(
        ax_gfp;
        xlim = plot_kwargs[:xlim],
        ylim = plot_kwargs[:ylim],
        xlabel = show_dissimilarity ? "" : plot_kwargs[:xlabel],
        ylabel = ylabel_gfp,
    )

    # Set grid using the shared function
    _set_axis_grid!(
        ax_gfp;
        xgrid = plot_kwargs[:xgrid],
        ygrid = plot_kwargs[:ygrid],
        xminorgrid = plot_kwargs[:xminorgrid],
        yminorgrid = plot_kwargs[:yminorgrid],
    )

    if plot_kwargs[:add_x_origin]
        vlines!(ax_gfp, [0.0], color = :black, linewidth = 1, linestyle = :dash)
    end

    # Add legend if multiple datasets
    if length(datasets_filtered) > 1 && plot_kwargs[:legend]
        axislegend(ax_gfp, position = :rt)
    end

    panel_idx += 1

    # Panel 3: Global Dissimilarity (if requested)
    if show_dissimilarity
        ax_diss = Axis(fig[panel_idx, 1])
        push!(axes, ax_diss)

        # Determine y-label
        ylabel_diss = normalize ? "Dissimilarity (%)" : "Dissimilarity"

        # Plot dissimilarity for each dataset
        for (i, result) in enumerate(results)
            color = colors[mod1(i, length(colors))]
            lines!(
                ax_diss,
                result.time,
                result.dissimilarity,
                color = color,
                linewidth = plot_kwargs[:linewidth],
                linestyle = resolved_linestyle,
            )
        end

        if plot_kwargs[:show_title]
            ax_diss.title = "Global Dissimilarity"
        end

        # Set axis limits and labels using the shared function
        _set_axis_properties!(
            ax_diss;
            xlim = plot_kwargs[:xlim],
            ylim = plot_kwargs[:ylim],
            xlabel = plot_kwargs[:xlabel],
            ylabel = ylabel_diss,
        )

        # Set grid using the shared function
        _set_axis_grid!(
            ax_diss;
            xgrid = plot_kwargs[:xgrid],
            ygrid = plot_kwargs[:ygrid],
            xminorgrid = plot_kwargs[:xminorgrid],
            yminorgrid = plot_kwargs[:yminorgrid],
        )

        if plot_kwargs[:add_x_origin]
            vlines!(ax_diss, [0.0], color = :black, linewidth = 1, linestyle = :dash)
        end
    end

    # Display the plot if requested
    if plot_kwargs[:display_plot]
        _display_figure(fig)
    end

    _set_window_title("Makie")
    return (fig = fig, axes = axes)
end

function plot_gfp(gfp_data::DataFrame; normalize::Union{Nothing,Bool} = nothing, kwargs...)

    # Merge user kwargs and default kwargs
    plot_kwargs = _merge_plot_kwargs(PLOT_GFP_KWARGS, _normalize_plot_kwargs(kwargs))

    # Check required columns
    if !hasproperty(gfp_data, :time)
        @minimal_error("DataFrame must have :time column")
    end
    if !hasproperty(gfp_data, :gfp)
        @minimal_error("DataFrame must have :gfp column")
    end

    # Determine if dissimilarity should be plotted
    show_dissimilarity = plot_kwargs[:show_dissimilarity] && hasproperty(gfp_data, :dissimilarity)

    # Determine number of panels
    n_panels = 1 + (show_dissimilarity ? 1 : 0)

    # Create figure
    fig = _eegfun_figure(size = (800, 300 * n_panels))
    axes = Axis[]

    panel_idx = 1

    # GFP panel
    ax_gfp = Axis(fig[panel_idx, 1])
    push!(axes, ax_gfp)

    is_normalized = _is_gfp_df_normalized(gfp_data, normalize)
    ylabel_gfp = if !isnothing(plot_kwargs[:ylabel])
        plot_kwargs[:ylabel]
    else
        is_normalized ? "GFP (%)" : "GFP (μV)"
    end

    resolved_color = _resolve_theme_linecolor(ax_gfp, plot_kwargs[:color])
    resolved_linestyle = isnothing(plot_kwargs[:linestyle]) ? :solid : plot_kwargs[:linestyle]

    lines!(ax_gfp, gfp_data.time, gfp_data.gfp, color = resolved_color, linewidth = plot_kwargs[:linewidth], linestyle = resolved_linestyle)

    if plot_kwargs[:show_title]
        ax_gfp.title = isnothing(plot_kwargs[:plot_title]) ? "Global Field Power" : plot_kwargs[:plot_title]
    end

    # Set axis limits and labels using the shared function
    _set_axis_properties!(
        ax_gfp;
        xlim = plot_kwargs[:xlim],
        ylim = plot_kwargs[:ylim],
        xlabel = show_dissimilarity ? "" : plot_kwargs[:xlabel],
        ylabel = ylabel_gfp,
    )

    # Set grid using the shared function
    _set_axis_grid!(
        ax_gfp;
        xgrid = plot_kwargs[:xgrid],
        ygrid = plot_kwargs[:ygrid],
        xminorgrid = plot_kwargs[:xminorgrid],
        yminorgrid = plot_kwargs[:yminorgrid],
    )

    if plot_kwargs[:add_x_origin]
        vlines!(ax_gfp, [0.0], color = :black, linewidth = 1, linestyle = :dash)
    end

    panel_idx += 1

    # Dissimilarity panel (if requested and available)
    if show_dissimilarity
        ax_diss = Axis(fig[panel_idx, 1])
        push!(axes, ax_diss)

        ylabel_diss = is_normalized ? "Dissimilarity (%)" : "Dissimilarity"

        lines!(
            ax_diss,
            gfp_data.time,
            gfp_data.dissimilarity,
            color = resolved_color,
            linewidth = plot_kwargs[:linewidth],
            linestyle = resolved_linestyle,
        )

        if plot_kwargs[:show_title]
            ax_diss.title = "Global Dissimilarity"
        end

        # Set axis limits and labels using the shared function
        _set_axis_properties!(
            ax_diss;
            xlim = plot_kwargs[:xlim],
            ylim = plot_kwargs[:ylim],
            xlabel = plot_kwargs[:xlabel],
            ylabel = ylabel_diss,
        )

        # Set grid using the shared function
        _set_axis_grid!(
            ax_diss;
            xgrid = plot_kwargs[:xgrid],
            ygrid = plot_kwargs[:ygrid],
            xminorgrid = plot_kwargs[:xminorgrid],
            yminorgrid = plot_kwargs[:yminorgrid],
        )

        if plot_kwargs[:add_x_origin]
            vlines!(ax_diss, [0.0], color = :black, linewidth = 1, linestyle = :dash)
        end
    end

    _set_window_title("GFP")
    if plot_kwargs[:display_plot]
        _display_figure(fig)
    end
    _set_window_title("Makie")

    return (fig = fig, axes = axes)
end

function plot_gfp(gfp_data::AbstractVector{<:DataFrame}; normalize::Union{Nothing,Bool} = nothing, kwargs...)

    # Merge user kwargs and default kwargs
    plot_kwargs = _merge_plot_kwargs(PLOT_GFP_KWARGS, _normalize_plot_kwargs(kwargs))

    # Check required columns in all DataFrames
    for (i, df) in enumerate(gfp_data)
        if !hasproperty(df, :time)
            @minimal_error("DataFrame $i must have :time column")
        end
        if !hasproperty(df, :gfp)
            @minimal_error("DataFrame $i must have :gfp column")
        end
    end

    # Check if dissimilarity is available in all datasets
    show_dissimilarity = plot_kwargs[:show_dissimilarity] && all(hasproperty(df, :dissimilarity) for df in gfp_data)

    # Determine number of panels
    n_panels = 1 + (show_dissimilarity ? 1 : 0)

    # Create figure
    fig = _eegfun_figure(size = (800, 300 * n_panels))
    axes = Axis[]

    panel_idx = 1

    # GFP panel
    ax_gfp = Axis(fig[panel_idx, 1])
    push!(axes, ax_gfp)

    is_normalized = isempty(gfp_data) ? false : _is_gfp_df_normalized(gfp_data[1], normalize)
    ylabel_gfp = if !isnothing(plot_kwargs[:ylabel])
        plot_kwargs[:ylabel]
    else
        is_normalized ? "GFP (%)" : "GFP (μV)"
    end

    resolved_color = _resolve_theme_linecolor(ax_gfp, plot_kwargs[:color])
    resolved_linestyle = isnothing(plot_kwargs[:linestyle]) ? :solid : plot_kwargs[:linestyle]

    # Plot GFP for each dataset
    colors = length(gfp_data) == 1 ? [resolved_color] : (!isnothing(plot_kwargs[:color]) ? [resolved_color] : Makie.wong_colors())
    for (i, df) in enumerate(gfp_data)
        color = colors[mod1(i, length(colors))]
        lines!(
            ax_gfp,
            df.time,
            df.gfp,
            color = color,
            linewidth = plot_kwargs[:linewidth],
            linestyle = resolved_linestyle,
            label = "Condition $i",
        )
    end

    if plot_kwargs[:show_title]
        ax_gfp.title = isnothing(plot_kwargs[:plot_title]) ? "Global Field Power" : plot_kwargs[:plot_title]
    end

    # Set axis limits and labels using the shared function
    _set_axis_properties!(
        ax_gfp;
        xlim = plot_kwargs[:xlim],
        ylim = plot_kwargs[:ylim],
        xlabel = show_dissimilarity ? "" : plot_kwargs[:xlabel],
        ylabel = ylabel_gfp,
    )

    # Set grid using the shared function
    _set_axis_grid!(
        ax_gfp;
        xgrid = plot_kwargs[:xgrid],
        ygrid = plot_kwargs[:ygrid],
        xminorgrid = plot_kwargs[:xminorgrid],
        yminorgrid = plot_kwargs[:yminorgrid],
    )

    if plot_kwargs[:add_x_origin]
        vlines!(ax_gfp, [0.0], color = :black, linewidth = 1, linestyle = :dash)
    end

    # Add legend
    if length(gfp_data) > 1 && plot_kwargs[:legend]
        axislegend(ax_gfp, position = :rt)
    end

    panel_idx += 1

    # Dissimilarity panel (if requested and available)
    if show_dissimilarity
        ax_diss = Axis(fig[panel_idx, 1])
        push!(axes, ax_diss)

        ylabel_diss = is_normalized ? "Dissimilarity (%)" : "Dissimilarity"

        for (i, df) in enumerate(gfp_data)
            color = colors[mod1(i, length(colors))]
            lines!(ax_diss, df.time, df.dissimilarity, color = color, linewidth = plot_kwargs[:linewidth], linestyle = resolved_linestyle)
        end

        if plot_kwargs[:show_title]
            ax_diss.title = "Global Dissimilarity"
        end

        # Set axis limits and labels using the shared function
        _set_axis_properties!(
            ax_diss;
            xlim = plot_kwargs[:xlim],
            ylim = plot_kwargs[:ylim],
            xlabel = plot_kwargs[:xlabel],
            ylabel = ylabel_diss,
        )

        # Set grid using the shared function
        _set_axis_grid!(
            ax_diss;
            xgrid = plot_kwargs[:xgrid],
            ygrid = plot_kwargs[:ygrid],
            xminorgrid = plot_kwargs[:xminorgrid],
            yminorgrid = plot_kwargs[:yminorgrid],
        )

        if plot_kwargs[:add_x_origin]
            vlines!(ax_diss, [0.0], color = :black, linewidth = 1, linestyle = :dash)
        end
    end

    _set_window_title("GFP")
    if plot_kwargs[:display_plot]
        _display_figure(fig)
    end
    _set_window_title("Makie")

    return (fig = fig, axes = axes)
end
