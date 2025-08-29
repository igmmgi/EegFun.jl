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
    :add_xy_origin => true,
    :interactive => true,  # New: enable/disable keyboard interactivity
)

# =============================================================================
# INTERACTIVITY CONSTANTS
# =============================================================================

const ERP_KEYBOARD_ACTIONS = Dict(
    Keyboard.up => :up,
    Keyboard.down => :down,
    Keyboard.left => :left,
    Keyboard.right => :right
)

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
- `interactive`: Whether to enable keyboard interactivity (default: true)

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

# Disable interactivity
plot_erp(dat, interactive = false)
```

# Interactive Controls
When `interactive = true` (default):
- **Up Arrow**: Zoom in on Y-axis (compress Y limits)
- **Down Arrow**: Zoom out on Y-axis (expand Y limits)
- **Left Arrow**: Zoom in on X-axis (compress time range)
- **Right Arrow**: Zoom out on X-axis (expand time range)
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
        datasets;
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

    # set default plot title only for single layouts
    # For grid/topo layouts, we want individual channel names, not a global title
    if plot_kwargs[:show_title] && plot_kwargs[:title] == "" && layout == :single
        plot_kwargs[:title] = length(all_plot_channels) == 1 ? string(all_plot_channels[1]) : "$(print_vector(all_plot_channels))"
        if plot_kwargs[:average_channels]
            plot_kwargs[:title] = "Avg: $(print_vector(original_channels))"
        end
    end
    
    # Create figure and apply layout system
    fig = Figure()
    plot_layout = create_layout(layout, all_plot_channels, first(dat_subset).layout)
    axes, channels = apply_layout!(fig, plot_layout; plot_kwargs...)
    
    # Now do the actual plotting for each axis
    for (ax, channel) in zip(axes, channels)
        @info "plot_erp: plotting channel: $channel, plot_layout.type: $(plot_layout.type)"
        channels_to_plot = plot_layout.type == :single ? all_plot_channels : [channel]
        _plot_erp!(ax, dat_subset, channels_to_plot; plot_kwargs...)
    end
    
    # Apply our axis stuff
    _apply_axis_properties!.(axes; plot_kwargs...)
    _apply_layout_axis_properties!(axes, plot_layout; plot_kwargs...) # slightly different for grid and topo layouts
    
    # Link axes for consistent navigation
    length(axes) > 1 && linkaxes!(axes...)
    
    # Add keyboard interactivity if enabled
    if plot_kwargs[:interactive]
        _setup_erp_interactivity!(fig, axes)
    end
    
    return fig, axes
end

"""
    _setup_erp_interactivity!(fig::Figure, axes::Vector{Axis})

Set up keyboard interactivity for ERP plots.
"""
function _setup_erp_interactivity!(fig::Figure, axes::Vector{Axis})
    # Handle keyboard events
    on(events(fig).keyboardbutton) do event
        if event.action in (Keyboard.press, Keyboard.repeat) && haskey(ERP_KEYBOARD_ACTIONS, event.key)
            action = ERP_KEYBOARD_ACTIONS[event.key]
            _handle_erp_navigation!(axes, action)
        end
    end
end

"""
    _handle_erp_navigation!(axes::Vector{Axis}, action::Symbol)

Handle navigation actions for ERP plots.
"""
function _handle_erp_navigation!(axes::Vector{Axis}, action::Symbol)
    # Only zoom the first axis - the linkaxes! will handle synchronizing all others
    ax = first(axes)
    if action == :up
        ymore!(ax)
    elseif action == :down
        yless!(ax)
    elseif action == :left
        xless!(ax)
    elseif action == :right
        xmore!(ax)
    end
end

"""
    _erp_zoom_in_y!(ax::Axis)

Zoom in on Y-axis by compressing the limits (zoom in on waveforms).
"""
function ymore!(ax::Axis)
    ylims!(ax, ax.yaxis.attributes.limits[] .* 0.9)
end

"""
    _erp_zoom_out_y!(ax::Axis)

Zoom out on Y-axis by expanding the limits (zoom out from waveforms).
"""
function yless!(ax::Axis)
    ylims!(ax, ax.yaxis.attributes.limits[] .* 1.1)
end

"""
    xmore!(ax::Axis)

Zoom in on X-axis by compressing the limits (zoom in on time range).
"""
function xmore!(ax::Axis)
    xlims!(ax, ax.xaxis.attributes.limits[] .* 0.9)
end

"""
    xless!(ax::Axis)

Zoom out on X-axis by expanding the limits (zoom out from time range).
"""
function xless!(ax::Axis)
    xlims!(ax, ax.xaxis.attributes.limits[] .* 1.1)
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
            linestyle = linestyles[(idx-1)%length(linestyles)+1] # wrap probably not great but better than crash
            dataset_color = dataset_colors[idx]
        else
            linestyle = kwargs[:linestyle]
            dataset_color = kwargs[:color]
        end
        
        # Plot ALL channels for this dataset
        for (ch_idx, channel) in enumerate(channels)

            # labels/ colours
            label = length(datasets) > 1 ?  string("Cond: ", idx, " ", channel) : string(channel)
            color = length(channels) > 1 ? channel_colors[ch_idx] : dataset_color
            
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
    if kwargs[:add_xy_origin]
        vlines!(ax, [0], color = :black, linewidth = 0.5)
        hlines!(ax, [0], color = :black, linewidth = 0.5)
    end
    
    # Show legend if requested and there are multiple channels or datasets
    if kwargs[:legend] && (length(channels) > 1 || length(datasets) > 1)
        axislegend(ax, framevisible = false, position = :lt)
    end
    
    return ax
end
