# plot_erp: Unified ERP plotting with layout system support

# =============================================================================
# DEFAULT KEYWORD ARGUMENTS
# =============================================================================

const DEFAULT_ERP_KWARGS = Dict(
    :xlim => nothing,
    :ylim => nothing,
    :xlabel => "Time (S)",
    :ylabel => "mV",
    :title => "",
    :show_title => true,
    :linewidth => 1,
    :color => :black,
    :linestyle => :solid,
    :colormap => :jet,
    :yreversed => false,
    :average_channels => false,
    :legend => true,
    :legend_label => "",
    :xgrid => false,
    :ygrid => false,
    :xminorgrid => false,
    :yminorgrid => false,
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

    # Merge user kwargs and default kwargs
    default_kwargs = copy(DEFAULT_ERP_KWARGS)
    plot_kwargs = merge(default_kwargs, kwargs)
    
    # data subsetting
    dat_subset = subset(
        datasets[1];
        channel_selection = channel_selection,
        sample_selection = sample_selection,
        include_extra = true,
    )

    if plot_kwargs[:average_channels]
        original_channels = channel_labels(dat_subset) # keep for better default plot title
        dat_subset = channel_average(dat_subset, reduce = true)
    end
    
    selected_channels = channel_labels(dat_subset)
    extra_channels = extra_labels(dat_subset)
    all_plot_channels = vcat(selected_channels, extra_channels)

    # set default plot title
    if plot_kwargs[:show_title] && plot_kwargs[:title] == ""
        plot_kwargs[:title] = length(all_plot_channels) == 1 ? string(all_plot_channels[1]) : "$(print_vector(all_plot_channels))"
        if plot_kwargs[:average_channels]
            plot_kwargs[:title] = "Avg: $(print_vector(original_channels))"
        end
    end
    
    # Create figure and apply layout system
    fig = Figure()
    plot_layout = create_layout(layout, all_plot_channels, dat_subset.layout)

    # Apply layout to create axes and get channel assignments
    axes, channel_assignments = apply_layout!(fig, plot_layout; plot_kwargs...)
    
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
                _plot_erp!(ax, [dat_subset], all_plot_channels; plot_kwargs...)
            else
                _plot_erp!(ax, [dat_subset], [channel]; plot_kwargs...)
            end
        end
    end
    
    # Apply our axis stuff
    _apply_axis_properties!.(axes; plot_kwargs...)
    _apply_layout_axis_properties!(axes, plot_layout; plot_kwargs...) # slightly different for grid and topo layouts
    
    # Link axes for consistent navigation
    length(axes) > 1 && linkaxes!(axes...)
    
    return fig, axes
end


"""
    _plot_erp!(ax::Axis, datasets::Vector{ErpData}, channels::Vector{Symbol}; kwargs...)

Internal function to plot ERP data on an axis.
Handles both single and multiple datasets.
Note: datasets should already be subset based on channel_selection and sample_selection.
"""
function _plot_erp!(ax::Axis, datasets::Vector{ErpData}, channels::Vector{Symbol}; kwargs...)
    
    # Use defaults + overrides
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
    
    # Add zero lines with tick marks
    vlines!(ax, [0], color = :black, linewidth = 0.5)
    hlines!(ax, [0], color = :black, linewidth = 0.5)
    
    # Show legend if requested and there are multiple channels or datasets
    if kwargs[:legend] && (length(channels) > 1 || length(datasets) > 1)
        axislegend(ax, framevisible = false, position = :lt)
    end
    
    return ax
end


