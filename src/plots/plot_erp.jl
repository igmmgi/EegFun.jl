# plot_erp: Unified ERP plotting with layout system support

# =============================================================================
# DEFAULT KEYWORD ARGUMENTS
# =============================================================================

const DEFAULT_ERP_KWARGS = Dict(
    :xlim => nothing,
    :ylim => nothing,
    :xlabel => "Time (S)",
    :ylabel => "mV",
    :linewidth => 1,
    :color => :black,
    :linestyle => :solid,
    :colormap => :jet,
    :yreversed => false,
    :average_channels => false,
    :legend => true,
    :legend_label => "",
)

# =============================================================================

"""
    plot_erp(dat::ErpData; 
             layout::Union{Symbol, PlotLayout, Vector{Int}} = :single,
             channel_selection::Function = channels(),
             sample_selection::Function = samples(),
             kwargs...)

Create ERP plots with flexible layout options.

# Arguments
- `dat::ErpData`: ERP data structure
- `layout`: Layout specification:
  - `:single` (default): Single plot with all channels
  - `:grid`: Auto-calculated grid layout
  - `:topo`: Topographic layout based on channel positions
  - `PlotLayout`: Custom layout object
  - `Vector{Int}`: Custom grid dimensions [rows, cols]
  
- `channel_selection::Function`: Function that returns boolean vector for channel filtering
- `sample_selection::Function`: Function that returns boolean vector for sample filtering
- `kwargs`: Additional keyword arguments

# Keyword Arguments
- `xlim`: X-axis limits (default: auto-calculated)
- `ylim`: Y-axis limits (default: auto-calculated)
- `title`: Plot title (default: auto-generated)
- `xlabel`: X-axis label (default: "Time (S)")
- `ylabel`: Y-axis label (default: "mV")
- `linewidth`: Line width (default: 2)
- `color`: Line color for single channel (default: :black)
- `linestyle`: Line style (default: :solid)
- `colormap`: Color map for multiple channels (default: :jet)
- `yreversed`: Whether to reverse Y-axis (default: false)
- `average_channels`: Whether to average channels (default: false)
- `legend`: Whether to show legend (default: true)
- `legend_label`: Legend label prefix (default: "")
- `hidedecorations`: Whether to hide axis decorations in grid/topo layouts (default: false)
- `theme_fontsize`: Font size for theme (default: 24)
- `plot_width`: Plot width for topo layout (default: 0.10)
- `plot_height`: Plot height for topo layout (default: 0.10)
- `margin`: Margin between plots for topo layout (default: 0.02)

# Examples
```julia
# Single plot with all channels
plot_erp(dat)

# Grid layout
plot_erp(dat, layout = :grid)

# Custom grid dimensions
plot_erp(dat, layout = [3, 4])

# Topographic layout
plot_erp(dat, layout = :topo)

# Custom layout object
layout = create_grid_layout(channels(dat), rows = 2, cols = 3)
plot_erp(dat, layout = layout)
```
"""
function plot_erp(dat::ErpData; 
                 layout::Union{Symbol, PlotLayout, Vector{Int}} = :single,
                 channel_selection::Function = channels(),
                 sample_selection::Function = samples(),
                 kwargs...)
    return plot_erp([dat]; layout=layout, channel_selection=channel_selection, 
                    sample_selection=sample_selection, kwargs...)
end

"""
    plot_erp(datasets::Vector{ErpData}; 
             layout::Union{Symbol, PlotLayout, Vector{Int}} = :single,
             channel_selection::Function = channels(), 
             sample_selection::Function = samples(), 
             kwargs...)

Plot multiple ERP datasets on the same axis (e.g., conditions).
"""
function plot_erp(datasets::Vector{ErpData}; 
                 layout::Union{Symbol, PlotLayout, Vector{Int}} = :single,
                 channel_selection::Function = channels(), 
                 sample_selection::Function = samples(), 
                 kwargs...)
    
    dat_subset = subset(
        datasets[1];
        channel_selection = channel_selection,
        sample_selection = sample_selection,
        include_extra = true,
    )
    
    selected_channels = channel_labels(dat_subset)
    extra_channels = extra_labels(dat_subset)
    all_plot_channels = vcat(selected_channels, extra_channels)
    
    # Create figure and apply layout system
    fig = Figure()
    plot_layout = create_layout(layout, all_plot_channels, dat_subset.layout)

    # Use defaults with title override
    default_kwargs = copy(DEFAULT_ERP_KWARGS)
    default_kwargs[:title] = nothing
    plot_kwargs = merge(default_kwargs, kwargs)

    # Apply layout to create axes and get channel assignments
    axes, channel_assignments = apply_layout!(fig, plot_layout; kwargs...)
    
    # Now do the actual plotting for each axis
    for (ax, (idx, channel)) in zip(axes, channel_assignments)
        # For grid and topo layouts, disable legend by default
        if plot_layout.type == :grid || plot_layout.type == :topo
            layout_kwargs = Dict{Symbol, Any}(kwargs)
            layout_kwargs[:legend] = false
            _plot_erp!(ax, datasets, [channel]; layout_kwargs...)
        else
            # For single layout, plot all channels on this axis
            if plot_layout.type == :single
                _plot_erp!(ax, datasets, all_plot_channels; kwargs...)
            else
                _plot_erp!(ax, datasets, [channel]; kwargs...)
            end
        end
    end
    
    # Apply common axis properties (but preserve grid-specific axis cleanup)
    for ax in axes
        _apply_axis_properties!(ax; kwargs...)
    end
    
    # For grid layouts, ensure axis labels are properly cleaned up
    if plot_layout.type == :grid
        for (idx, ax) in enumerate(axes)
            row = fld(idx-1, plot_layout.cols) + 1
            col = mod(idx-1, plot_layout.cols) + 1
            
            # Re-apply grid axis properties to ensure they're not overridden
            _set_grid_axis_properties!(ax, plot_layout, plot_layout.channels[idx], row, col, plot_layout.rows, plot_layout.cols; kwargs...)
        end
    end
    
    # For topo layouts, remove axis labels and ticks, and add scale plot
    if plot_layout.type == :topo
        for ax in axes
            # Remove axis labels and ticks for topographic plots
            ax.xlabel = ""
            ax.ylabel = ""
            hidedecorations!(ax, grid = false, ticks = true, ticklabels = true)
            # Hide all axis spines for clean topographic appearance
            hidespines!(ax)
        end
    end
    
    # Link axes for consistent navigation
    if length(axes) > 1
        linkaxes!(axes...)
    end
    
    return fig, axes
end


# ===== INTERNAL PLOTTING FUNCTIONS =====
"""
    _plot_erp!(ax::Axis, datasets::Vector{ErpData}, channels::Vector{Symbol}; kwargs...)

Internal function to plot ERP data on an axis.
Handles both single and multiple datasets.
Note: datasets should already be subset based on channel_selection and sample_selection.
"""
function _plot_erp!(ax::Axis, datasets::Vector{ErpData}, channels::Vector{Symbol}; kwargs...)
    
    # Use defaults
    kwargs = merge(DEFAULT_ERP_KWARGS, kwargs)
    
    # Handle individual channel plotting
    # Styling for multiple datasets and channels
    linestyles = [:solid, :dot, :dash, :dashdot, :dashdotdot]
    dataset_colors = Makie.cgrad(kwargs[:colormap], length(datasets), categorical = true)
    channel_colors = Makie.cgrad(kwargs[:colormap], length(channels), categorical = true)
    
    # Plot each dataset for ALL channels in this subplot
    for (idx, dat) in enumerate(datasets)
        # Set styling for this condition
        if length(datasets) > 1
            linestyle = linestyles[(idx-1)%length(linestyles)+1]
            dataset_color = dataset_colors[idx]
        else
            linestyle = kwargs[:linestyle]
            dataset_color = kwargs[:color]
        end
        
        # Plot ALL channels for this dataset
        for (ch_idx, channel) in enumerate(channels)
            # Create unique label for this condition + channel
            if length(datasets) > 1
                label = string("Cond: ", idx, " ", channel)
            else
                label = string(channel)
            end
            
            # Use different color for each channel within each condition
            if length(channels) > 1
                color = channel_colors[ch_idx]
            else
                color = dataset_color
            end
            
            lines!(
                ax,
                dat.data[!, :time],
                dat.data[!, channel],
                color = color,
                linewidth = kwargs[:linewidth],
                linestyle = linestyle,
                label = label,
            )
        end
    end
    
    # Set title to show all channels
    ax.title = length(channels) == 1 ? string(channels[1]) : "$(print_vector(channels))"
    
    # Add zero lines
    vlines!(ax, [0], color = :black, linewidth = 0.5)
    hlines!(ax, [0], color = :black, linewidth = 0.5)
    
    # Set axis labels
    ax.xlabel = kwargs[:xlabel]
    ax.ylabel = kwargs[:ylabel]
    
    # Apply limits if provided
    if !isnothing(kwargs[:xlim])
        xlims!(ax, kwargs[:xlim])
    end
    
    if !isnothing(kwargs[:ylim])
        ylims!(ax, kwargs[:ylim])
    end
    
    # Show legend if requested
    if kwargs[:legend]
        axislegend(ax, framevisible = false, position = :lt)
    end
    
    return ax
end


