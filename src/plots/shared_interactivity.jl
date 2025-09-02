# Shared interactivity functions for plot_erp and plot_epochs
# This file provides common interactive features to avoid code duplication

# =============================================================================
# SHARED CONSTANTS
# =============================================================================

const SHARED_KEYBOARD_ACTIONS = Dict(
    Keyboard.up => :up,
    Keyboard.down => :down,
    Keyboard.left => :left,
    Keyboard.right => :right
)

# =============================================================================
# SHARED SELECTION STATE
# =============================================================================

"""
    SharedSelectionState

Generic selection state that can be used by both plot_erp and plot_epochs.
"""
mutable struct SharedSelectionState
    active::Observable{Bool}
    bounds::Observable{Tuple{Float64,Float64}}
    visible::Observable{Bool}
    rectangles::Vector{Makie.Poly}  # Store rectangles for all axes (time selection)
    channel_rectangles::Vector{Makie.Poly}  # Store channel selection rectangles
    selection_rectangles::Vector{Makie.Poly}  # Store multiple selection rectangles
    selection_bounds::Vector{Tuple{Float64,Float64,Float64,Float64}}  # Store bounds for each selection
    current_selection_idx::Union{Int, Nothing}  # Index of currently active selection
    
    function SharedSelectionState(axes::Vector{Axis})
        rectangles = Makie.Poly[]
        for ax in axes
            initial_points = [Point2f(0.0, 0.0)]
            poly_element = poly!(ax, initial_points, color = (:blue, 0.3), visible = false)
            push!(rectangles, poly_element)
        end
        new(Observable(false), Observable((0.0, 0.0)), Observable(false), rectangles, Makie.Poly[], Makie.Poly[], Tuple{Float64,Float64,Float64,Float64}[], nothing)
    end
end

# =============================================================================
# SHARED INTERACTIVITY FUNCTIONS
# =============================================================================

"""
    _setup_shared_interactivity!(fig::Figure, axes::Vector{Axis}, keyboard_actions::Dict)

Set up keyboard interactivity for plots.
"""
function _setup_shared_interactivity!(fig::Figure, axes::Vector{Axis}, keyboard_actions::Dict = SHARED_KEYBOARD_ACTIONS)
    # Handle keyboard events
    on(events(fig).keyboardbutton) do event
        if event.action in (Keyboard.press, Keyboard.repeat) && haskey(keyboard_actions, event.key)
            action = keyboard_actions[event.key]
            _handle_shared_navigation!(axes, action)
        end
    end
end

"""
    _handle_shared_navigation!(axes::Vector{Axis}, action::Symbol)

Handle navigation actions for plots.
"""
function _handle_shared_navigation!(axes::Vector{Axis}, action::Symbol)
    for ax in axes
        current_xlims = xlims(ax)
        current_ylims = ylims(ax)
        
        if action == :up
            # Zoom in on Y-axis (compress Y limits)
            y_center = (current_ylims[1] + current_ylims[2]) / 2
            y_range = current_ylims[2] - current_ylims[1]
            new_y_range = y_range * 0.8
            new_ylims = (y_center - new_y_range/2, y_center + new_y_range/2)
            ylims!(ax, new_ylims)
        elseif action == :down
            # Zoom out on Y-axis (expand Y limits)
            y_center = (current_ylims[1] + current_ylims[2]) / 2
            y_range = current_ylims[2] - current_ylims[1]
            new_y_range = y_range * 1.25
            new_ylims = (y_center - new_y_range/2, y_center + new_y_range/2)
            ylims!(ax, new_ylims)
        elseif action == :left
            # Zoom in on X-axis (compress time range)
            x_center = (current_xlims[1] + current_xlims[2]) / 2
            x_range = current_xlims[2] - current_xlims[1]
            new_x_range = x_range * 0.8
            new_xlims = (x_center - new_x_range/2, x_center + new_x_range/2)
            xlims!(ax, new_xlims)
        elseif action == :right
            # Zoom out on X-axis (expand time range)
            x_center = (current_xlims[1] + current_xlims[2]) / 2
            x_range = current_xlims[2] - current_xlims[1]
            new_x_range = x_range * 1.25
            new_xlims = (x_center - new_x_range/2, x_center + new_x_range/2)
            xlims!(ax, new_xlims)
        end
    end
end

# =============================================================================
# SHARED SELECTION FUNCTIONS
# =============================================================================

"""
    _is_mouse_in_axis(ax::Axis, mouse_pos)

Check if mouse position is within the axis bounds.
"""
function _is_mouse_in_axis(ax::Axis, mouse_pos)
    ax_pos = ax.scene.viewport[]
    ax_x = ax_pos.origin[1]
    ax_y = ax_pos.origin[2]
    ax_w = ax_pos.widths[1]
    ax_h = ax_pos.widths[2]
    
    return (ax_x <= mouse_pos[1] <= ax_x + ax_w) && (ax_y <= mouse_pos[2] <= ax_y + ax_h)
end

"""
    _is_within_selection(selection_state::SharedSelectionState, mouse_x)

Check if mouse position is within an existing selection.
"""
function _is_within_selection(selection_state::SharedSelectionState, mouse_x)
    if !selection_state.visible[]
        return false
    end
    
    x_min, x_max = minmax(selection_state.bounds[]...)
    return x_min <= mouse_x <= x_max
end

"""
    _start_shared_selection!(ax::Axis, selection_state::SharedSelectionState, mouse_x)

Start a new time selection.
"""
function _start_shared_selection!(ax::Axis, selection_state::SharedSelectionState, mouse_x)
    selection_state.active[] = true
    selection_state.bounds[] = (mouse_x, mouse_x)
    selection_state.visible[] = false
    _update_shared_selection!(ax, selection_state, mouse_x, mouse_x)
end

"""
    _finish_shared_selection!(ax::Axis, selection_state::SharedSelectionState, mouse_x)

Finish the current time selection.
"""
function _finish_shared_selection!(ax::Axis, selection_state::SharedSelectionState, mouse_x)
    selection_state.active[] = false
    selection_state.visible[] = true
    selection_state.bounds[] = (selection_state.bounds[][1], mouse_x)
    _update_shared_selection!(ax, selection_state, selection_state.bounds[][1], mouse_x)
    # Make all rectangles visible
    for rect in selection_state.rectangles
        rect.visible[] = true
    end
end

"""
    _update_shared_selection!(ax::Axis, selection_state::SharedSelectionState, x1, x2)

Update the selection rectangle while dragging.
"""
function _update_shared_selection!(ax::Axis, selection_state::SharedSelectionState, x1, x2)
    # Update all rectangles across all axes
    for (i, rect) in enumerate(selection_state.rectangles)
        # Get the y-limits for the corresponding axis
        if i <= length(selection_state.rectangles)
            # Use fixed y-range for consistency across all subplots
            rect[1] = Point2f[
                Point2f(Float64(x1), Float64(-1000)),  # Use fixed y-range for consistency
                Point2f(Float64(x2), Float64(-1000)),
                Point2f(Float64(x2), Float64(1000)),
                Point2f(Float64(x1), Float64(1000)),
            ]
            rect.visible[] = true
        end
    end
end

"""
    _clear_shared_selection!(selection_state::SharedSelectionState)

Clear all time selections.
"""
function _clear_shared_selection!(selection_state::SharedSelectionState)
    # Clear all rectangles
    for rect in selection_state.rectangles
        rect[1] = [Point2f(0.0, 0.0)]
        rect.visible[] = false
    end
    selection_state.bounds[] = (0.0, 0.0)
    selection_state.visible[] = false
end

# =============================================================================
# SHARED CHANNEL SELECTION FUNCTIONS
# =============================================================================

"""
    _start_figure_channel_selection!(fig::Figure, selection_state::SharedSelectionState, plot_layout, data, channel_selection_active)

Start channel selection with Ctrl + drag.
"""
function _start_figure_channel_selection!(fig::Figure, selection_state::SharedSelectionState, plot_layout, data, channel_selection_active)
    if !_is_topo_layout(plot_layout) && !_is_grid_layout(plot_layout)
        return
    end
    selection_start = events(fig).mouseposition[]
    selection_state.bounds[] = (selection_start[1], selection_start[2])
    channel_selection_active[] = true
    selection_state.current_selection_idx = length(selection_state.selection_rectangles) + 1
end

"""
    _update_figure_channel_selection!(fig::Figure, selection_state::SharedSelectionState, plot_layout, data)

Update channel selection rectangle while dragging.
"""
function _update_figure_channel_selection!(fig::Figure, selection_state::SharedSelectionState, plot_layout, data)
    if !_is_topo_layout(plot_layout) && !_is_grid_layout(plot_layout)
        return
    end
    current_pos = events(fig).mouseposition[]
    start_screen = selection_state.bounds[]
    fig_size = size(fig.scene)
    start_norm = (start_screen[1] / fig_size[1], start_screen[2] / fig_size[2])
    end_norm = (current_pos[1] / fig_size[1], current_pos[2] / fig_size[2])
    x1, x2 = minmax(start_norm[1], end_norm[1])
    y1, y2 = minmax(start_norm[2], end_norm[2])
    current_idx = selection_state.current_selection_idx
    if isnothing(current_idx) || current_idx > length(selection_state.selection_rectangles)
        rect_points = [Point2f(x1, y1), Point2f(x2, y1), Point2f(x2, y2), Point2f(x1, y2)]
        rect = poly!(fig.scene, rect_points, color = (:blue, 0.3), strokecolor = :red, strokewidth = 2, overdraw = true, space = :relative)
        push!(selection_state.selection_rectangles, rect)
        push!(selection_state.selection_bounds, (x1, y1, x2, y2))
    else
        rect = selection_state.selection_rectangles[current_idx]
        rect_points = [Point2f(x1, y1), Point2f(x2, y1), Point2f(x2, y2), Point2f(x1, y2)]
        rect[1] = rect_points
        selection_state.selection_bounds[current_idx] = (x1, y1, x2, y2)
    end
end

"""
    _finish_figure_channel_selection!(fig::Figure, selection_state::SharedSelectionState, plot_layout, data, channel_selection_active, axes::Vector{Axis})

Finish channel selection and highlight selected channels.
"""
function _finish_figure_channel_selection!(fig::Figure, selection_state::SharedSelectionState, plot_layout, data, channel_selection_active, axes::Vector{Axis})
    if (!_is_topo_layout(plot_layout) && !_is_grid_layout(plot_layout)) || !channel_selection_active[]
        return
    end
    selection_end = events(fig).mouseposition[]
    start_screen = selection_state.bounds[]
    fig_size = size(fig.scene)
    start_norm = (start_screen[1] / fig_size[1], start_screen[2] / fig_size[2])
    end_norm = (selection_end[1] / fig_size[1], selection_end[2] / fig_size[2])
    x1, x2 = minmax(start_norm[1], end_norm[1])
    y1, y2 = minmax(start_norm[2], end_norm[2])
    x1 = max(0.0, min(1.0, x1))
    x2 = max(0.0, min(1.0, x2))
    y1 = max(0.0, min(1.0, y1))
    y2 = max(0.0, min(1.0, y2))
    # Handle both PlotLayout struct and Dict
    channels = isa(plot_layout, PlotLayout) ? plot_layout.channels : get(plot_layout, :channels, Symbol[])
    axes_rects = _get_axes_rectangles(axes, channels, fig)
    all_selected_channels = Symbol[]
    all_overlapping_axes_rects = []
    for (selection_idx, selection_bounds) in enumerate(selection_state.selection_bounds)
        for (axis_rect, channels) in axes_rects
            if _rectangles_overlap(selection_bounds, axis_rect)
                append!(all_selected_channels, channels)
                push!(all_overlapping_axes_rects, (axis_rect, channels))
            end
        end
    end
    current_selection_bounds = (x1, y1, x2, y2)
    for (axis_rect, channels) in axes_rects
        if _rectangles_overlap(current_selection_bounds, axis_rect)
            append!(all_selected_channels, channels)
            push!(all_overlapping_axes_rects, (axis_rect, channels))
        end
    end
    unique_channels = unique(all_selected_channels)
    if !isempty(selection_state.channel_rectangles)
        for rect in selection_state.channel_rectangles
            delete!(rect.parent, rect)
        end
        empty!(selection_state.channel_rectangles)
    end
    if !isempty(all_overlapping_axes_rects)
        new_rectangles = _draw_channel_rectangles!(fig, all_overlapping_axes_rects)
        append!(selection_state.channel_rectangles, new_rectangles)
    end
    if !isempty(unique_channels)
        println("Total selected channels across all regions: $unique_channels")
        println("Number of selection areas: $(length(selection_state.selection_rectangles))")
    else
        println("No channels selected")
    end
    channel_selection_active[] = false
end

"""
    _clear_all_shared_channel_selections!(fig::Figure, selection_state::SharedSelectionState)

Clear all channel selections.
"""
function _clear_all_shared_channel_selections!(fig::Figure, selection_state::SharedSelectionState)
    for rect in selection_state.selection_rectangles
        delete!(rect.parent, rect)
    end
    empty!(selection_state.selection_rectangles)
    for rect in selection_state.channel_rectangles
        delete!(rect.parent, rect)
    end
    empty!(selection_state.channel_rectangles)
    empty!(selection_state.selection_bounds)
    println("Cleared all channel selections")
end

"""
    _draw_channel_rectangles!(fig::Figure, all_overlapping_axes_rects)

Draw highlight rectangles for selected channels.
"""
function _draw_channel_rectangles!(fig::Figure, all_overlapping_axes_rects)
    rectangles = []
    for (axis_rect, channels) in all_overlapping_axes_rects
        rect = poly!(fig.scene, [Point2f(axis_rect[1], axis_rect[2]), Point2f(axis_rect[3], axis_rect[2]), Point2f(axis_rect[3], axis_rect[4]), Point2f(axis_rect[1], axis_rect[4])], color = (:blue, 0.6), strokecolor = :darkblue, strokewidth = 2, overdraw = true, space = :relative)
        push!(rectangles, rect)
    end
    return rectangles
end

"""
    _get_axes_rectangles(axes::Vector{Axis}, channels::Vector{Symbol}, fig::Figure)

Get normalized rectangle coordinates for all axes.
"""
function _get_axes_rectangles(axes::Vector{Axis}, channels::Vector{Symbol}, fig::Figure)
    axes_rects = []
    fig_size = size(fig.scene)
    for (idx, ax) in enumerate(axes)
        if idx <= length(channels)
            # Use viewport instead of deprecated px_area
            ax_pos = ax.scene.viewport[]
            ax_x = ax_pos.origin[1]
            ax_y = ax_pos.origin[2]
            ax_w = ax_pos.widths[1]
            ax_h = ax_pos.widths[2]
            norm_x1 = ax_x / fig_size[1]
            norm_y1 = ax_y / fig_size[2]
            norm_x2 = (ax_x + ax_w) / fig_size[1]
            norm_y2 = (ax_y + ax_h) / fig_size[2]
            push!(axes_rects, ((norm_x1, norm_y1, norm_x2, norm_y2), [channels[idx]]))
        end
    end
    return axes_rects
end

# =============================================================================
# SHARED UTILITY FUNCTIONS
# =============================================================================

"""
    _is_topo_layout(plot_layout)

Check if the layout is topographic.
"""
function _is_topo_layout(plot_layout)
    if isa(plot_layout, PlotLayout)
        return plot_layout.type == :topo
    elseif isa(plot_layout, Dict)
        return get(plot_layout, :type, :single) == :topo
    else
        return false
    end
end

"""
    _is_grid_layout(plot_layout)

Check if the layout is grid-based.
"""
function _is_grid_layout(plot_layout)
    if isa(plot_layout, PlotLayout)
        return plot_layout.type == :grid
    elseif isa(plot_layout, Dict)
        return get(plot_layout, :type, :single) == :grid
    else
        return false
    end
end

"""
    _rectangles_overlap(rect1, rect2)

Check if two rectangles overlap.
"""
function _rectangles_overlap(rect1, rect2)
    x1_min, y1_min, x1_max, y1_max = rect1
    x2_min, y2_min, x2_max, y2_max = rect2
    return !(x1_max < x2_min || x2_max < x1_min || y1_max < y2_min || y2_max < y1_min)
end


