# Essential helper functions for channel selection in plot_epochs

function _start_figure_channel_selection!(fig::Figure, selection_state::EpochsSelectionState, plot_layout, data, channel_selection_active)
    if !_is_topo_layout(plot_layout) && !_is_grid_layout(plot_layout)
        return
    end
    selection_start = events(fig).mouseposition[]
    selection_state.bounds[] = (selection_start[1], selection_start[2])
    channel_selection_active[] = true
    selection_state.current_selection_idx = length(selection_state.selection_rectangles) + 1
end

function _update_figure_channel_selection!(fig::Figure, selection_state::EpochsSelectionState, plot_layout, data)
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

function _finish_figure_channel_selection!(fig::Figure, selection_state::EpochsSelectionState, plot_layout, data, channel_selection_active, axes::Vector{Axis})
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
    axes_rects = _get_axes_rectangles(axes, plot_layout[:channels], fig)
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

function _clear_all_epochs_channel_selections!(fig::Figure, selection_state::EpochsSelectionState)
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

function _draw_channel_rectangles!(fig::Figure, all_overlapping_axes_rects)
    rectangles = []
    for (axis_rect, channels) in all_overlapping_axes_rects
        rect = poly!(fig.scene, [Point2f(axis_rect[1], axis_rect[2]), Point2f(axis_rect[3], axis_rect[2]), Point2f(axis_rect[3], axis_rect[4]), Point2f(axis_rect[1], axis_rect[4])], color = (:blue, 0.6), strokecolor = :darkblue, strokewidth = 2, overdraw = true, space = :relative)
        push!(rectangles, rect)
    end
    return rectangles
end

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

_is_topo_layout(plot_layout) = get(plot_layout, :type, :single) == :topo
_is_grid_layout(plot_layout) = get(plot_layout, :type, :single) == :grid

function _rectangles_overlap(rect1, rect2)
    x1_min, y1_min, x1_max, y1_max = rect1
    x2_min, y2_min, x2_max, y2_max = rect2
    return !(x1_max < x2_min || x2_max < x1_min || y1_max < y2_min || y2_max < y1_min)
end
