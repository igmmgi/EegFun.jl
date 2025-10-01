# Shared interactivity functions for plot_erp and plot_epochs

# =============================================================================
# SHARED CONSTANTS
# =============================================================================
const SHARED_KEYBOARD_ACTIONS =
    Dict(Keyboard.up => :up, Keyboard.down => :down, Keyboard.left => :left, Keyboard.right => :right)

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
    rectangles::Vector{Makie.Poly}
    channel_rectangles::Vector{Makie.Poly}
    selection_rectangles::Vector{Makie.Poly}
    selection_bounds::Vector{Tuple{Float64,Float64,Float64,Float64}}
    current_selection_idx::Union{Int,Nothing}

    function SharedSelectionState(axes::Vector{Axis})
        # Create time selection rectangles for each axis
        rectangles = [
            poly!(
                ax,
                [Point2f(0.0, 0.0)],
                color = (:blue, 0.3),
                strokecolor = :darkblue,
                strokewidth = 1,
                visible = false,
                overdraw = true,
                space = :data,
            ) for ax in axes
        ]

        new(
            Observable(false),
            Observable((0.0, 0.0)),
            Observable(false),
            rectangles,
            Makie.Poly[],
            Makie.Poly[],
            Tuple{Float64,Float64,Float64,Float64}[],
            nothing,
        )
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
    on(events(fig).keyboardbutton) do event
        if event.action in (Keyboard.press, Keyboard.repeat) && haskey(keyboard_actions, event.key)
            action = keyboard_actions[event.key]
            _handle_shared_navigation!(axes, action)
        end
    end
end

"""
    _handle_shared_navigation!(axes::Vector{Axis}, action::Symbol)

Handle navigation actions for plots using arrow keys.
"""
function _handle_shared_navigation!(axes::Vector{Axis}, action::Symbol)
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
    ymore!(ax::Axis)

Zoom in on Y-axis by compressing the limits (zoom in on waveforms).
"""
function ymore!(ax::Axis)
    ylims!(ax, ax.yaxis.attributes.limits[] .* 0.8)
end

"""
    yless!(ax::Axis)

Zoom out on Y-axis by expanding the limits (zoom out from waveforms).
"""
function yless!(ax::Axis)
    ylims!(ax, ax.yaxis.attributes.limits[] .* 1.25)
end

"""
    xmore!(ax::Axis)

Zoom in on X-axis by compressing the limits (zoom in on time range).
"""
function xmore!(ax::Axis)
    xlims!(ax, ax.xaxis.attributes.limits[] .* 0.8)
end

"""
    xless!(ax::Axis)

Zoom out on X-axis by expanding the limits (zoom out from time range).
"""
function xless!(ax::Axis)
    xlims!(ax, ax.xaxis.attributes.limits[] .* 1.25)
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
            # Get the actual y-limits of the axis for proper rectangle positioning
            y_lims = ax.yaxis.attributes.limits[]
            y_min = y_lims[1]
            y_max = y_lims[2]

            rect[1] = Point2f[
                Point2f(Float64(x1), Float64(y_min)),
                Point2f(Float64(x2), Float64(y_min)),
                Point2f(Float64(x2), Float64(y_max)),
                Point2f(Float64(x1), Float64(y_max)),
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
function _start_figure_channel_selection!(
    fig::Figure,
    selection_state::SharedSelectionState,
    plot_layout,
    data,
    channel_selection_active,
)
    if plot_layout.type != :topo && plot_layout.type != :grid
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
    if plot_layout.type != :topo && plot_layout.type != :grid
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
        rect = poly!(
            fig.scene,
            rect_points,
            color = (:blue, 0.3),
            strokecolor = :red,
            strokewidth = 2,
            overdraw = true,
            space = :relative,
        )
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
function _finish_figure_channel_selection!(
    fig::Figure,
    selection_state::SharedSelectionState,
    plot_layout,
    data,
    channel_selection_active,
    axes::Vector{Axis},
)
    if (plot_layout.type != :topo && plot_layout.type != :grid) || !channel_selection_active[]
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
        rect = poly!(
            fig.scene,
            [
                Point2f(axis_rect[1], axis_rect[2]),
                Point2f(axis_rect[3], axis_rect[2]),
                Point2f(axis_rect[3], axis_rect[4]),
                Point2f(axis_rect[1], axis_rect[4]),
            ],
            color = (:blue, 0.6),
            strokecolor = :darkblue,
            strokewidth = 2,
            overdraw = true,
            space = :relative,
        )
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
# SHARED UNIFIED SELECTION SETUP
# =============================================================================

"""
    _find_active_axis(axes::Vector{Axis}, mouse_pos)

Find which axis the mouse is currently over.
"""
function _find_active_axis(axes::Vector{Axis}, mouse_pos)
    for ax in axes
        if _is_mouse_in_axis(ax, mouse_pos)
            return ax
        end
    end
    return nothing
end

"""
    _setup_keyboard_tracking(fig::Figure)

Set up keyboard event tracking for Shift and Ctrl keys.
Returns Refs for tracking key states.
"""
function _setup_keyboard_tracking(fig::Figure)
    shift_pressed = Ref(false)
    ctrl_pressed = Ref(false)

    on(events(fig).keyboardbutton) do key_event
        if key_event.key == Keyboard.left_shift
            shift_pressed[] = key_event.action == Keyboard.press
        elseif key_event.key == Keyboard.left_control
            ctrl_pressed[] = key_event.action == Keyboard.press
        end
    end

    return shift_pressed, ctrl_pressed
end

"""
    _handle_mouse_button_events(fig::Figure, axes::Vector{Axis}, selection_state::SharedSelectionState, data, right_click_handler)

Handle mouse button events for selection.
"""
function _handle_mouse_button_events(
    fig::Figure,
    axes::Vector{Axis},
    selection_state::SharedSelectionState,
    data,
    right_click_handler,
    shift_pressed,
)
    on(events(fig).mousebutton) do event
        mouse_pos = events(fig).mouseposition[]
        active_ax = _find_active_axis(axes, mouse_pos)

        isnothing(active_ax) && return

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
            if !isnothing(right_click_handler)
                right_click_handler(selection_state, mouse_x, data)
            else
                @info "TODO: implement right click functionality"
            end
        end
    end
end

"""
    _handle_mouse_movement(fig::Figure, axes::Vector{Axis}, selection_state::SharedSelectionState)

Handle mouse movement for updating selection rectangles.
"""
function _handle_mouse_movement(fig::Figure, axes::Vector{Axis}, selection_state::SharedSelectionState)
    on(events(fig).mouseposition) do _
        if selection_state.active[]
            mouse_pos = events(fig).mouseposition[]
            active_ax = _find_active_axis(axes, mouse_pos)

            if !isnothing(active_ax)
                world_pos = mouseposition(active_ax)[1]
                _update_shared_selection!(active_ax, selection_state, selection_state.bounds[][1], world_pos)
            end
        end
    end
end

"""
    _setup_unified_selection!(fig::Figure, axes::Vector{Axis}, selection_state::SharedSelectionState, data, plot_layout, right_click_handler=nothing)

Set up unified mouse selection that works across all layouts.
Uses figure-level events to avoid conflicts with multiple axis handlers.
"""
function _setup_unified_selection!(
    fig::Figure,
    axes::Vector{Axis},
    selection_state::SharedSelectionState,
    data,
    plot_layout,
    right_click_handler = nothing,
)
    shift_pressed, ctrl_pressed = _setup_keyboard_tracking(fig)
    _handle_mouse_button_events(fig, axes, selection_state, data, right_click_handler, shift_pressed)
    _handle_mouse_movement(fig, axes, selection_state)
end

# =============================================================================
# SHARED CHANNEL SELECTION EVENT FUNCTIONS
# =============================================================================

"""
    _setup_channel_selection_events!(fig::Figure, selection_state::SharedSelectionState, plot_layout, data, axes::Vector{Axis}, layout_type::Symbol)

Set up figure-level event handlers for channel selection in topo and grid layouts.
"""
function _setup_channel_selection_events!(
    fig::Figure,
    selection_state::SharedSelectionState,
    plot_layout,
    data,
    axes::Vector{Axis},
    layout_type::Symbol,
)
    channel_selection_active = Ref(false)
    shift_pressed = Ref(false)
    ctrl_pressed = Ref(false)

    # Keyboard event tracking
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

    # Mouse events
    on(events(fig).mousebutton) do event
        if event.button == Mouse.left && ctrl_pressed[]
            if event.action == Mouse.press
                _start_figure_channel_selection!(fig, selection_state, plot_layout, data, channel_selection_active)
            elseif event.action == Mouse.release && channel_selection_active[]
                _finish_figure_channel_selection!(
                    fig,
                    selection_state,
                    plot_layout,
                    data,
                    channel_selection_active,
                    axes,
                )
            end
        end

        if event.button == Mouse.left && event.action == Mouse.press && !ctrl_pressed[] && !shift_pressed[]
            _clear_all_shared_channel_selections!(fig, selection_state)
        end

        if event.button == Mouse.right && event.action == Mouse.press
            @info "TODO: implement right click functionality for $layout_type"
        end
    end

    # Mouse movement
    on(events(fig).mouseposition) do _
        if channel_selection_active[]
            _update_figure_channel_selection!(fig, selection_state, plot_layout, data)
        end
    end
end

# =============================================================================
# SHARED UTILITY FUNCTIONS
# =============================================================================

"""
    _rectangles_overlap(rect1, rect2)

Check if two rectangles overlap.
"""
function _rectangles_overlap(rect1, rect2)
    x1_min, y1_min, x1_max, y1_max = rect1
    x2_min, y2_min, x2_max, y2_max = rect2
    return !(x1_max < x2_min || x2_max < x1_min || y1_max < y2_min || y2_max < y1_min)
end
