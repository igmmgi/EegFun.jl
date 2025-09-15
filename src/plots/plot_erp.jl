# =============================================================================
# DEFAULT KEYWORD ARGUMENTS
# =============================================================================
const PLOT_ERP_KWARGS = Dict{Symbol,Tuple{Any,String}}(
    # Display parameters
    :display_plot => (true, "Whether to display the plot"),
    
    # Axis limits and labels
    :xlim => (nothing, "X-axis limits as (min, max) tuple. If nothing, automatically determined"),
    :ylim => (nothing, "Y-axis limits as (min, max) tuple. If nothing, automatically determined"),
    :xlabel => ("Time (S)", "Label for x-axis"),
    :ylabel => ("μV", "Label for y-axis"),
    
    # Title
    :title => ("", "Plot title"),
    :show_title => (true, "Whether to show the title"),
    
    # Line styling
    :linewidth => (1, "Line width for ERP traces"),
    :color => (:black, "Color for ERP traces"),
    :linestyle => (:solid, "Line style for ERP traces"),
    :colormap => (:jet, "Colormap for multi-channel plots"),
    
    # Plot configuration
    :yreversed => (false, "Whether to reverse the y-axis"),
    :average_channels => (false, "Whether to average across channels"),
    :interactive => (true, "Whether to enable interactive features"),
    
    # Legend
    :legend => (true, "Whether to show the legend"),
    :legend_label => ("", "Custom label for the legend"),
    
    # Grid
    :xgrid => (false, "Whether to show x-axis grid"),
    :ygrid => (false, "Whether to show y-axis grid"),
    :xminorgrid => (false, "Whether to show x-axis minor grid"),
    :yminorgrid => (false, "Whether to show y-axis minor grid"),
    
    # Origin lines
    :add_xy_origin => (true, "Whether to add origin lines at x=0 and y=0"),
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

$(generate_kwargs_doc(PLOT_ERP_KWARGS))

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

# Mouse Selection
- **Shift + Left Click + Drag**: Select a time region (blue rectangle)
- **Right Click on Selection**: Open context menu with plot options
  - Topoplot (multiquadratic)
  - Topoplot (spherical_spline)
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
    plot_kwargs = _merge_plot_kwargs(PLOT_ERP_KWARGS, kwargs)
    
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
        channels_to_plot = plot_layout.type == :single ? all_plot_channels : [channel]
        @info "plot_erp ($layout): $(print_vector(channels_to_plot))"
        _plot_erp!(ax, dat_subset, channels_to_plot; plot_kwargs...)
    end
    
    # Apply our axis stuff
    _apply_axis_properties!.(axes; plot_kwargs...)
    _apply_layout_axis_properties!(axes, plot_layout; plot_kwargs...) # slightly different for grid and topo layouts
    
    # Link axes for consistent navigation
    length(axes) > 1 && linkaxes!(axes...)
    
    # Add keyboard interactivity if enabled
    if plot_kwargs[:interactive]
        _setup_shared_interactivity!(fig, axes)
        
        # Disable default interactions that conflict with our custom selection
        # Need to disable on ALL axes for grid layouts to work properly
        for ax in axes
            deregister_interaction!(ax, :rectanglezoom)
        end
        
        # Set up selection system for all axes (will work with linked axes)
        # Pass the ORIGINAL datasets (not dat_subset) so we can subset by time for topo plots
        selection_state = SharedSelectionState(axes)
        
        # Set up selection system that works for all layouts
        # Use figure-level events to avoid conflicts with multiple axis handlers
        _setup_unified_selection!(fig, axes, selection_state, datasets, plot_layout, _handle_erp_right_click!)
        
        # Set up channel selection events for topo and grid layouts
        if plot_layout.type == :topo
            _setup_channel_selection_events!(fig, selection_state, plot_layout, datasets, axes, :topo)
        elseif plot_layout.type == :grid
            _setup_channel_selection_events!(fig, selection_state, plot_layout, datasets, axes, :grid)
        end
    end
    
    if plot_kwargs[:display_plot]
        display_figure(fig)
    end
    
    return fig, axes
end

# Unified selection setup is now handled by _setup_unified_selection! in shared_interactivity.jl

function _handle_erp_right_click!(selection_state, mouse_x, data)
    if selection_state.visible[] && _is_within_selection(selection_state, mouse_x)
        _show_erp_context_menu!(selection_state, data)
    end
end

function _show_erp_context_menu!(selection_state, data)
    # Create the menu figure
    menu_fig = Figure()
    plot_types = ["Topoplot (multiquadratic)", "Topoplot (spherical_spline)"]

    menu_buttons = [Button(menu_fig[idx, 1], label = plot_type) for (idx, plot_type) in enumerate(plot_types)]

    for btn in menu_buttons
        on(btn.clicks) do n
            original_data, x_min, x_max = _subset_erp_selected_data(selection_state, data)
            
            # Create time-based sample selection for the topo plot
            time_sample_selection = x -> (x.time .>= x_min) .& (x.time .<= x_max)
            
            if btn.label[] == "Topoplot (multiquadratic)"
                plot_topography(original_data, sample_selection = time_sample_selection, method = :multiquadratic)
            elseif btn.label[] == "Topoplot (spherical_spline)"
                plot_topography(original_data, sample_selection = time_sample_selection, method = :spherical_spline)
            end
        end
    end

    new_screen = getfield(Main, :GLMakie).Screen()
    display(new_screen, menu_fig)
end

function _subset_erp_selected_data(selection_state, data)
    x_min, x_max = minmax(selection_state.bounds[]...)
    
    # Return the original data and time bounds instead of subsetting
    # This preserves all electrodes for topo plots
    return (data, x_min, x_max)
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
