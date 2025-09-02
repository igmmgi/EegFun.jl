# plot_erp: Unified ERP plotting with layout system support

# Shared interactivity functions are now included in src/eegfun.jl

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
# ERP-SPECIFIC SELECTION STATE
# =============================================================================

# Use the shared selection state
const ErpSelectionState = SharedSelectionState

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
        _setup_shared_interactivity!(fig, axes)
        
        # Disable default interactions that conflict with our custom selection
        # Need to disable on ALL axes for grid layouts to work properly
        for ax in axes
            deregister_interaction!(ax, :rectanglezoom)
        end
        
        # Set up selection system for all axes (will work with linked axes)
        # Pass the ORIGINAL datasets (not dat_subset) so we can subset by time for topo plots
        selection_state = ErpSelectionState(axes)
        
        # Set up selection system that works for all layouts
        # Use figure-level events to avoid conflicts with multiple axis handlers
        _setup_erp_selection_unified!(fig, axes, selection_state, datasets, plot_layout)
    end
    
    return fig, axes
end

"""
    _setup_erp_selection_unified!(fig::Figure, axes::Vector{Axis}, selection_state::ErpSelectionState, data, plot_layout::PlotLayout)

Set up unified mouse selection and context menu for ERP plots that works across all layouts.
Uses figure-level events to avoid conflicts with multiple axis handlers.
"""
function _setup_erp_selection_unified!(fig::Figure, axes::Vector{Axis}, selection_state::ErpSelectionState, data, plot_layout::PlotLayout)
    # Track if Shift and Ctrl are currently pressed
    shift_pressed = Ref(false)
    ctrl_pressed = Ref(false)

    # Use figure-level keyboard events for Shift and Ctrl tracking
    on(events(fig).keyboardbutton) do key_event
        if key_event.key == Keyboard.left_shift
            shift_pressed[] = key_event.action == Keyboard.press
        elseif key_event.key == Keyboard.left_control
            ctrl_pressed[] = key_event.action == Keyboard.press
        end
    end

    # Use figure-level mouse events to avoid conflicts
    on(events(fig).mousebutton) do event
        # Find which axis the mouse is over
        mouse_pos = events(fig).mouseposition[]
        active_ax = nothing
        
        for ax in axes
            if _is_mouse_in_axis(ax, mouse_pos)
                active_ax = ax
                break
            end
        end
        
        if isnothing(active_ax)
            return
        end

        mouse_x = mouseposition(active_ax)[1]

        if event.button == Mouse.left
            if event.action == Mouse.press
                if shift_pressed[] && _is_within_selection(selection_state, mouse_x)
                    _clear_shared_selection!(selection_state)
                elseif shift_pressed[]
                    _start_shared_selection!(active_ax, selection_state, mouse_x)
                end
            elseif event.action == Mouse.release && selection_state.active[]
                _finish_shared_selection!(active_ax, selection_state, mouse_x)
            end
        elseif event.button == Mouse.right && event.action == Mouse.press
            _handle_erp_right_click!(selection_state, mouse_x, data)
        end
    end

    # Update selection rectangle while dragging using figure-level mouse position
    on(events(fig).mouseposition) do _
        if selection_state.active[]
            # Update time selection rectangle
            mouse_pos = events(fig).mouseposition[]
            active_ax = nothing
            
            for ax in axes
                if _is_mouse_in_axis(ax, mouse_pos)
                    active_ax = ax
                    break
                end
            end
            
            if !isnothing(active_ax)
                world_pos = mouseposition(active_ax)[1]
                _update_shared_selection!(active_ax, selection_state, selection_state.bounds[][1], world_pos)
            end
        end
    end
    
    # Add separate figure-level event handlers for channel selection in topo and grid layouts
    # This provides Ctrl+Click+Drag for multiple channel selections, Left Click to clear all, Right Click for info
    if _is_topo_layout(plot_layout)
        _setup_topo_channel_selection_events!(fig, selection_state, plot_layout, data, axes)
    elseif _is_grid_layout(plot_layout)
        _setup_grid_channel_selection_events!(fig, selection_state, plot_layout, data, axes)
    end
end

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

"""
    _setup_topo_channel_selection_events!(fig::Figure, selection_state::ErpSelectionState, plot_layout::PlotLayout, data, axes::Vector{Axis})

Set up separate figure-level event handlers for channel selection in topo layouts.

# Mouse Controls:
- **Ctrl + Left Click + Drag**: Select channels (draw selection rectangle)
- **Multiple selections**: Each Ctrl+Click+Drag creates a new selection area
- **Left Click (without Ctrl)**: Clear all channel selections and selection areas
- **Right Click**: Print info TODO (for future functionality)
"""
function _setup_topo_channel_selection_events!(fig::Figure, selection_state::ErpSelectionState, plot_layout::PlotLayout, data, axes::Vector{Axis})
    channel_selection_active = Ref(false)  # Separate state for channel selection
    shift_pressed = Ref(false)
    ctrl_pressed = Ref(false)

    # Use figure-level keyboard events for Shift and Ctrl tracking
    on(events(fig).keyboardbutton) do key_event
        if key_event.key == Keyboard.left_shift
            shift_pressed[] = key_event.action == Keyboard.press
        elseif key_event.key == Keyboard.left_control
            ctrl_pressed[] = key_event.action == Keyboard.press
        end
        if key_event.action == Keyboard.release
            channel_selection_active[] = false
        end
    end   

    # Handle mouse events for channel selection
    on(events(fig).mousebutton) do event
        # Ctrl + Left Click + Drag: Select channels (existing functionality)
        if event.button == Mouse.left && ctrl_pressed[]
            if event.action == Mouse.press
                _start_figure_channel_selection!(fig, selection_state, plot_layout, data, channel_selection_active)
            elseif event.action == Mouse.release && channel_selection_active[]
                _finish_figure_channel_selection!(fig, selection_state, plot_layout, data, channel_selection_active, axes)
            end
        end
        
        # Left Click (without Ctrl): Clear all channel selections
        if event.button == Mouse.left && event.action == Mouse.press && !ctrl_pressed[] && !shift_pressed[]
            _clear_all_shared_channel_selections!(fig, selection_state)
        end
        
        # Right Click: Print info TODO (for future functionality)
        if event.button == Mouse.right && event.action == Mouse.press
            @info "TODO: implement right click functionality for ERP topo"
            println("TODO: implement right click functionality for ERP topo")
        end
    end
    
    # Handle mouse movement for updating channel selection rectangle
    on(events(fig).mouseposition) do _
        if channel_selection_active[]
            _update_figure_channel_selection!(fig, selection_state, plot_layout, data)
        end
    end
end

"""
    _setup_grid_channel_selection_events!(fig::Figure, selection_state::ErpSelectionState, plot_layout::PlotLayout, data, axes::Vector{Axis})

Set up separate figure-level event handlers for channel selection in grid layouts.
This provides the same functionality as topo layouts but for grid arrangements.
"""
function _setup_grid_channel_selection_events!(fig::Figure, selection_state::ErpSelectionState, plot_layout::PlotLayout, data, axes::Vector{Axis})
    # Track Ctrl key state for channel selection
    ctrl_pressed = Ref(false)
    channel_selection_active = Ref(false)  # Separate state for channel selection
    shift_pressed = Ref(false)
    
    # Use figure-level keyboard events for Shift and Ctrl tracking
    on(events(fig).keyboardbutton) do key_event
        if key_event.key == Keyboard.left_shift
            shift_pressed[] = key_event.action == Keyboard.press
        elseif key_event.key == Keyboard.left_control
            ctrl_pressed[] = key_event.action == Keyboard.press
        end
        if key_event.action == Keyboard.release
            channel_selection_active[] = false
        end
    end   

    # Handle mouse events for channel selection
    on(events(fig).mousebutton) do event
        # Ctrl + Left Click + Drag: Select channels (existing functionality)
        if event.button == Mouse.left && ctrl_pressed[]
            if event.action == Mouse.press
                _start_figure_channel_selection!(fig, selection_state, plot_layout, data, channel_selection_active)
            elseif event.action == Mouse.release && channel_selection_active[]
                _finish_figure_channel_selection!(fig, selection_state, plot_layout, data, channel_selection_active, axes)
            end
        end
        
        # Left Click (without Ctrl): Clear all channel selections
        if event.button == Mouse.left && event.action == Mouse.press && !ctrl_pressed[] && !shift_pressed[]
            _clear_all_shared_channel_selections!(fig, selection_state)
        end
        
        # Right Click: Print info TODO (for future functionality)
        if event.button == Mouse.right && event.action == Mouse.press
            @info "TODO: implement right click functionality for ERP grid"
            println("TODO: implement right click functionality for ERP grid")
        end
    end
    
    # Handle mouse movement for updating channel selection rectangle
    on(events(fig).mouseposition) do _
        if channel_selection_active[]
            _update_figure_channel_selection!(fig, selection_state, plot_layout, data)
        end
    end
end