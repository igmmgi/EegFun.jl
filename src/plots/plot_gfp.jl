"""
Plotting functions for Global Field Power (GFP) and Global Dissimilarity.
"""

# =============================================================================
# DEFAULT KEYWORD ARGUMENTS
# =============================================================================
const PLOT_GFP_KWARGS = Dict{Symbol,Tuple{Any,String}}(
    # Display parameters
    :display_plot => (true, "Whether to display the plot"),

    # Axis limits and labels
    :xlim => (nothing, "X-axis limits as (min, max) tuple. If nothing, automatically determined"),
    :ylim => (nothing, "Y-axis limits as (min, max) tuple. If nothing, automatically determined"),
    :xlabel => ("Time (s)", "Label for x-axis"),
    :ylabel => (nothing, "Label for y-axis (auto-determined based on normalize flag)"),

    # Title
    :title => ("", "Plot title"),
    :show_title => (true, "Whether to show the title"),

    # Line styling
    :linewidth => (2, "Line width for GFP/dissimilarity traces"),
    :color => (:black, "Color for traces"),
    :linestyle => (:solid, "Line style for traces"),

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

$(generate_kwargs_doc(PLOT_GFP_KWARGS))

# Examples
```julia
using eegfun, JLD2

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
"""
function plot_gfp(dat::ErpData; channel_selection::Function = channels(), normalize::Bool = true, kwargs...)
    return plot_gfp([dat]; channel_selection = channel_selection, normalize = normalize, kwargs...)
end


"""
    plot_gfp(datasets::Vector{ErpData};
             channel_selection::Function = channels(),
             normalize::Bool = true,
             kwargs...)

Plot Global Field Power for multiple ERP datasets (e.g., different conditions).

Each dataset is plotted with a different color/line style on the same axes.
"""
function plot_gfp(
    datasets::Vector{ErpData};
    channel_selection::Function = channels(),
    normalize::Bool = true,
    kwargs...,
)

    # Merge user kwargs and default kwargs
    plot_kwargs = _merge_plot_kwargs(PLOT_GFP_KWARGS, kwargs)

    # Calculate GFP (and optionally dissimilarity) for all datasets
    show_dissimilarity = plot_kwargs[:show_dissimilarity]
    show_erp_traces = plot_kwargs[:show_erp_traces]

    if show_dissimilarity
        # Calculate both GFP and dissimilarity
        results = [
            gfp_and_dissimilarity(dat; channel_selection = channel_selection, normalize = normalize) for dat in datasets
        ]
    else
        # Calculate only GFP
        results = [gfp(dat; channel_selection = channel_selection, normalize = normalize) for dat in datasets]
    end

    # Determine number of panels
    n_panels = 1 + (show_erp_traces ? 1 : 0) + (show_dissimilarity ? 1 : 0)

    # Create figure
    fig = Figure(size = (800, 200 * n_panels))

    panel_idx = 1

    # Panel 1: ERP traces (if requested)
    if show_erp_traces
        ax_erp = Axis(fig[panel_idx, 1])

        # Plot all channel traces for each dataset
        for (dataset_idx, dat) in enumerate(datasets)
            # Get selected channels
            selected_channels =
                get_selected_channels(dat, channel_selection, include_meta = false, include_extra = false)

            # Plot each channel
            for ch in selected_channels
                lines!(ax_erp, dat.data.time, dat.data[!, ch], color = (:gray, 0.5), linewidth = 0.5)
            end
        end

        ax_erp.xlabel = ""
        ax_erp.ylabel = "Amplitude (μV)"
        ax_erp.title =
            plot_kwargs[:show_title] ? (isempty(plot_kwargs[:title]) ? "EEG Channels" : plot_kwargs[:title]) : ""

        if plot_kwargs[:add_x_origin]
            vlines!(ax_erp, [0.0], color = :black, linewidth = 1, linestyle = :dash)
        end

        panel_idx += 1
    end

    # Panel 2 (or 1): GFP
    ax_gfp = Axis(fig[panel_idx, 1])

    # Determine y-label
    ylabel_gfp = if plot_kwargs[:ylabel] !== nothing
        plot_kwargs[:ylabel]
    else
        normalize ? "GFP (%)" : "GFP (μV)"
    end

    # Plot GFP for each dataset
    colors = length(datasets) == 1 ? [plot_kwargs[:color]] : Makie.wong_colors()
    for (i, result) in enumerate(results)
        color = colors[mod1(i, length(colors))]
        lines!(
            ax_gfp,
            result.time,
            result.gfp,
            color = color,
            linewidth = plot_kwargs[:linewidth],
            linestyle = plot_kwargs[:linestyle],
            label = length(datasets) > 1 ? "Condition $i" : nothing,
        )
    end

    # Apply styling
    ax_gfp.xlabel = show_dissimilarity ? "" : plot_kwargs[:xlabel]
    ax_gfp.ylabel = ylabel_gfp
    if !show_erp_traces
        ax_gfp.title =
            plot_kwargs[:show_title] ? (isempty(plot_kwargs[:title]) ? "Global Field Power" : plot_kwargs[:title]) : ""
    else
        ax_gfp.title = "Global Field Power"
    end

    # Set axis limits using the shared function
    _setup_axis_limits!(ax_gfp; xlim = plot_kwargs[:xlim], ylim = plot_kwargs[:ylim])

    # Set grid using the shared function
    _setup_axis_grid!(ax_gfp; 
                     xgrid = plot_kwargs[:xgrid], 
                     ygrid = plot_kwargs[:ygrid],
                     xminorgrid = plot_kwargs[:xminorgrid], 
                     yminorgrid = plot_kwargs[:yminorgrid])

    if plot_kwargs[:add_x_origin]
        vlines!(ax_gfp, [0.0], color = :black, linewidth = 1, linestyle = :dash)
    end

    # Add legend if multiple datasets
    if length(datasets) > 1
        axislegend(ax_gfp, position = :rt)
    end

    panel_idx += 1

    # Panel 3: Global Dissimilarity (if requested)
    if show_dissimilarity
        ax_diss = Axis(fig[panel_idx, 1])

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
                linestyle = plot_kwargs[:linestyle],
            )
        end

        ax_diss.xlabel = plot_kwargs[:xlabel]
        ax_diss.ylabel = ylabel_diss
        ax_diss.title = "Global Dissimilarity"

        # Set axis limits using the shared function
        _setup_axis_limits!(ax_diss; xlim = plot_kwargs[:xlim], ylim = plot_kwargs[:ylim])

        # Set grid using the shared function
        _setup_axis_grid!(ax_diss; 
                         xgrid = plot_kwargs[:xgrid], 
                         ygrid = plot_kwargs[:ygrid],
                         xminorgrid = plot_kwargs[:xminorgrid], 
                         yminorgrid = plot_kwargs[:yminorgrid])

        if plot_kwargs[:add_x_origin]
            vlines!(ax_diss, [0.0], color = :black, linewidth = 1, linestyle = :dash)
        end
    end

    # Display the plot if requested
    if plot_kwargs[:display_plot]
        display(fig)
    end

    return fig
end


"""
    plot_gfp(gfp_data::DataFrame; kwargs...)

Plot GFP from pre-computed GFP data.

This version accepts a DataFrame that was returned by the `gfp()` function,
allowing you to plot pre-computed results without recalculating.

# Arguments
- `gfp_data::DataFrame`: DataFrame containing `:time` and `:gfp` columns (from `gfp()` function)
- `kwargs`: Additional keyword arguments (same as main `plot_gfp`)

# Examples
```julia
# Calculate GFP separately
gfp_result = gfp(erp_data, normalize = true)

# Plot later
plot_gfp(gfp_result)

# Plot with custom styling
plot_gfp(gfp_result, color = :red, linewidth = 3)
```
"""
function plot_gfp(gfp_data::DataFrame; kwargs...)

    # Merge user kwargs and default kwargs
    plot_kwargs = _merge_plot_kwargs(PLOT_GFP_KWARGS, kwargs)

    # Check required columns
    if !hasproperty(gfp_data, :time)
        @minimal_error_throw("DataFrame must have :time column")
    end
    if !hasproperty(gfp_data, :gfp)
        @minimal_error_throw("DataFrame must have :gfp column")
    end

    # Determine if dissimilarity should be plotted
    show_dissimilarity = plot_kwargs[:show_dissimilarity] && hasproperty(gfp_data, :dissimilarity)

    # Determine number of panels
    n_panels = 1 + (show_dissimilarity ? 1 : 0)

    # Create figure
    fig = Figure(size = (800, 300 * n_panels))

    panel_idx = 1

    # GFP panel
    ax_gfp = Axis(fig[panel_idx, 1])

    # Determine if data is normalized (simple heuristic: check if values are 0-100)
    is_normalized = all(0 .<= gfp_data.gfp .<= 100)
    ylabel_gfp = if plot_kwargs[:ylabel] !== nothing
        plot_kwargs[:ylabel]
    else
        is_normalized ? "GFP (%)" : "GFP (μV)"
    end

    lines!(
        ax_gfp,
        gfp_data.time,
        gfp_data.gfp,
        color = plot_kwargs[:color],
        linewidth = plot_kwargs[:linewidth],
        linestyle = plot_kwargs[:linestyle],
    )

    ax_gfp.xlabel = show_dissimilarity ? "" : plot_kwargs[:xlabel]
    ax_gfp.ylabel = ylabel_gfp
    ax_gfp.title =
        plot_kwargs[:show_title] ? (isempty(plot_kwargs[:title]) ? "Global Field Power" : plot_kwargs[:title]) : ""

    # Set axis limits using the shared function
    _setup_axis_limits!(ax_gfp; xlim = plot_kwargs[:xlim], ylim = plot_kwargs[:ylim])

    # Set grid using the shared function
    _setup_axis_grid!(ax_gfp; 
                     xgrid = plot_kwargs[:xgrid], 
                     ygrid = plot_kwargs[:ygrid],
                     xminorgrid = plot_kwargs[:xminorgrid], 
                     yminorgrid = plot_kwargs[:yminorgrid])

    if plot_kwargs[:add_x_origin]
        vlines!(ax_gfp, [0.0], color = :black, linewidth = 1, linestyle = :dash)
    end

    panel_idx += 1

    # Dissimilarity panel (if requested and available)
    if show_dissimilarity
        ax_diss = Axis(fig[panel_idx, 1])

        # Determine if dissimilarity is normalized
        is_diss_normalized = all(0 .<= gfp_data.dissimilarity .<= 100)
        ylabel_diss = is_diss_normalized ? "Dissimilarity (%)" : "Dissimilarity"

        lines!(
            ax_diss,
            gfp_data.time,
            gfp_data.dissimilarity,
            color = plot_kwargs[:color],
            linewidth = plot_kwargs[:linewidth],
            linestyle = plot_kwargs[:linestyle],
        )

        ax_diss.xlabel = plot_kwargs[:xlabel]
        ax_diss.ylabel = ylabel_diss
        ax_diss.title = "Global Dissimilarity"

        # Set axis limits using the shared function
        _setup_axis_limits!(ax_diss; xlim = plot_kwargs[:xlim], ylim = plot_kwargs[:ylim])

        # Set grid using the shared function
        _setup_axis_grid!(ax_diss; 
                         xgrid = plot_kwargs[:xgrid], 
                         ygrid = plot_kwargs[:ygrid],
                         xminorgrid = plot_kwargs[:xminorgrid], 
                         yminorgrid = plot_kwargs[:yminorgrid])

        if plot_kwargs[:add_x_origin]
            vlines!(ax_diss, [0.0], color = :black, linewidth = 1, linestyle = :dash)
        end
    end

    # Display the plot if requested
    if plot_kwargs[:display_plot]
        display(fig)
    end

    return fig
end


"""
    plot_gfp(gfp_data::Vector{DataFrame}; kwargs...)

Plot GFP from multiple pre-computed GFP datasets.

# Examples
```julia
# Calculate GFP for multiple conditions
erps = load("participant_1_erps.jld2", "erps")
gfp_results = gfp.(erps, normalize = true)

# Plot all conditions together
plot_gfp(gfp_results)
```
"""
function plot_gfp(gfp_data::Vector{DataFrame}; kwargs...)

    # Merge user kwargs and default kwargs
    plot_kwargs = _merge_plot_kwargs(PLOT_GFP_KWARGS, kwargs)

    # Check required columns in all DataFrames
    for (i, df) in enumerate(gfp_data)
        if !hasproperty(df, :time)
            @minimal_error_throw("DataFrame $i must have :time column")
        end
        if !hasproperty(df, :gfp)
            @minimal_error_throw("DataFrame $i must have :gfp column")
        end
    end

    # Check if dissimilarity is available in all datasets
    show_dissimilarity = plot_kwargs[:show_dissimilarity] && all(hasproperty(df, :dissimilarity) for df in gfp_data)

    # Determine number of panels
    n_panels = 1 + (show_dissimilarity ? 1 : 0)

    # Create figure
    fig = Figure(size = (800, 300 * n_panels))

    panel_idx = 1

    # GFP panel
    ax_gfp = Axis(fig[panel_idx, 1])

    # Determine if data is normalized (check first dataset)
    is_normalized = all(0 .<= gfp_data[1].gfp .<= 100)
    ylabel_gfp = if plot_kwargs[:ylabel] !== nothing
        plot_kwargs[:ylabel]
    else
        is_normalized ? "GFP (%)" : "GFP (μV)"
    end

    # Plot GFP for each dataset
    colors = length(gfp_data) == 1 ? [plot_kwargs[:color]] : Makie.wong_colors()
    for (i, df) in enumerate(gfp_data)
        color = colors[mod1(i, length(colors))]
        lines!(
            ax_gfp,
            df.time,
            df.gfp,
            color = color,
            linewidth = plot_kwargs[:linewidth],
            linestyle = plot_kwargs[:linestyle],
            label = "Condition $i",
        )
    end

    ax_gfp.xlabel = show_dissimilarity ? "" : plot_kwargs[:xlabel]
    ax_gfp.ylabel = ylabel_gfp
    ax_gfp.title =
        plot_kwargs[:show_title] ? (isempty(plot_kwargs[:title]) ? "Global Field Power" : plot_kwargs[:title]) : ""

    # Set axis limits using the shared function
    _setup_axis_limits!(ax_gfp; xlim = plot_kwargs[:xlim], ylim = plot_kwargs[:ylim])

    # Set grid using the shared function
    _setup_axis_grid!(ax_gfp; 
                     xgrid = plot_kwargs[:xgrid], 
                     ygrid = plot_kwargs[:ygrid],
                     xminorgrid = plot_kwargs[:xminorgrid], 
                     yminorgrid = plot_kwargs[:yminorgrid])

    if plot_kwargs[:add_x_origin]
        vlines!(ax_gfp, [0.0], color = :black, linewidth = 1, linestyle = :dash)
    end

    # Add legend
    if length(gfp_data) > 1
        axislegend(ax_gfp, position = :rt)
    end

    panel_idx += 1

    # Dissimilarity panel (if requested and available)
    if show_dissimilarity
        ax_diss = Axis(fig[panel_idx, 1])

        is_diss_normalized = all(0 .<= gfp_data[1].dissimilarity .<= 100)
        ylabel_diss = is_diss_normalized ? "Dissimilarity (%)" : "Dissimilarity"

        for (i, df) in enumerate(gfp_data)
            color = colors[mod1(i, length(colors))]
            lines!(
                ax_diss,
                df.time,
                df.dissimilarity,
                color = color,
                linewidth = plot_kwargs[:linewidth],
                linestyle = plot_kwargs[:linestyle],
            )
        end

        ax_diss.xlabel = plot_kwargs[:xlabel]
        ax_diss.ylabel = ylabel_diss
        ax_diss.title = "Global Dissimilarity"

        # Set axis limits using the shared function
        _setup_axis_limits!(ax_diss; xlim = plot_kwargs[:xlim], ylim = plot_kwargs[:ylim])

        # Set grid using the shared function
        _setup_axis_grid!(ax_diss; 
                         xgrid = plot_kwargs[:xgrid], 
                         ygrid = plot_kwargs[:ygrid],
                         xminorgrid = plot_kwargs[:xminorgrid], 
                         yminorgrid = plot_kwargs[:yminorgrid])

        if plot_kwargs[:add_x_origin]
            vlines!(ax_diss, [0.0], color = :black, linewidth = 1, linestyle = :dash)
        end
    end

    # Display the plot if requested
    if plot_kwargs[:display_plot]
        display(fig)
    end

    return fig
end
