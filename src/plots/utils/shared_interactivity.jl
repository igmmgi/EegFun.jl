const SHARED_KEYBOARD_ACTIONS = Dict(Keyboard.up => :up, Keyboard.down => :down, Keyboard.left => :left, Keyboard.right => :right)

mutable struct SharedSelectionState
    active::Observable{Bool}
    bounds::Observable{Tuple{Float64,Float64}}
    visible::Observable{Bool}
    rectangles::Vector{Makie.Poly}
    channel_rectangles::Vector{Makie.Poly}
    selection_rectangles::Vector{Makie.Poly}
    selection_bounds::Vector{Tuple{Float64,Float64,Float64,Float64}}
    current_selection_idx::Union{Int,Nothing}
    selection_color::Symbol
    selection_alpha::Float64

    function SharedSelectionState(axes::Vector{Axis}; selection_color::Symbol = :blue, selection_alpha::Float64 = 0.3)
        rectangles = [
            poly!(
                ax,
                [Point2f(0.0, 0.0)],
                color = (selection_color, selection_alpha),
                strokecolor = selection_color,
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
            selection_color,
            selection_alpha,
        )
    end
end

function _setup_shared_interactivity!(
    fig::Figure,
    axes::Vector{Axis},
    keyboard_actions::Dict = SHARED_KEYBOARD_ACTIONS;
    zoom_step::Float64 = 0.2,
)
    on(events(fig).keyboardbutton) do event
        if event.action in (Keyboard.press, Keyboard.repeat) && haskey(keyboard_actions, event.key)
            _handle_shared_navigation!(axes, keyboard_actions[event.key], zoom_step)
        end
    end
end

function _setup_shared_interactivity!(
    fig::Figure,
    axes::Vector{Axis},
    plot_type::Symbol,
    keyboard_actions::Dict = SHARED_KEYBOARD_ACTIONS;
    zoom_step::Float64 = 0.2,
)
    _setup_shared_interactivity!(fig, axes, keyboard_actions; zoom_step = zoom_step)
    _setup_help_interaction!(fig, plot_type)
end

function _handle_shared_navigation!(axes::Vector{Axis}, action::Symbol, zoom_step::Float64)
    ax = first(axes)
    action == :up && _ymore!(ax, zoom_step)
    action == :down && _yless!(ax, zoom_step)
    action == :left && _xless!(ax, zoom_step)
    action == :right && _xmore!(ax, zoom_step)
end

_ymore!(ax::Axis, zoom_step::Float64) = ylims!(ax, ax.yaxis.attributes.limits[] .* (1.0 - zoom_step))
_yless!(ax::Axis, zoom_step::Float64) = ylims!(ax, ax.yaxis.attributes.limits[] .* (1.0 / (1.0 - zoom_step)))
_xmore!(ax::Axis, zoom_step::Float64) = xlims!(ax, ax.xaxis.attributes.limits[] .* (1.0 - zoom_step))
_xless!(ax::Axis, zoom_step::Float64) = xlims!(ax, ax.xaxis.attributes.limits[] .* (1.0 / (1.0 - zoom_step)))

function _is_mouse_in_axis(ax::Axis, mouse_pos)
    ax_pos = ax.scene.viewport[]
    ax_x, ax_y = ax_pos.origin[1], ax_pos.origin[2]
    ax_w, ax_h = ax_pos.widths[1], ax_pos.widths[2]
    return (ax_x <= mouse_pos[1] <= ax_x + ax_w) && (ax_y <= mouse_pos[2] <= ax_y + ax_h)
end

function _is_within_selection(selection_state::SharedSelectionState, mouse_x)
    selection_state.visible[] || return false
    x_min, x_max = minmax(selection_state.bounds[]...)
    return x_min <= mouse_x <= x_max
end

function _start_shared_selection!(ax::Axis, selection_state::SharedSelectionState, mouse_x)
    selection_state.active[] = true
    selection_state.bounds[] = (mouse_x, mouse_x)
    selection_state.visible[] = false
    _update_shared_selection!(ax, selection_state, mouse_x, mouse_x)
end

function _finish_shared_selection!(ax::Axis, selection_state::SharedSelectionState, mouse_x)
    selection_state.active[] = false
    selection_state.visible[] = true
    x1 = selection_state.bounds[][1]
    selection_state.bounds[] = (x1, mouse_x)
    _update_shared_selection!(ax, selection_state, x1, mouse_x)
    for rect in selection_state.rectangles
        rect.visible[] = true
    end
end

function _update_shared_selection!(ax::Axis, selection_state::SharedSelectionState, x1, x2)
    y_min, y_max = ax.yaxis.attributes.limits[]
    pts = Point2f[
        Point2f(Float64(x1), Float64(y_min)),
        Point2f(Float64(x2), Float64(y_min)),
        Point2f(Float64(x2), Float64(y_max)),
        Point2f(Float64(x1), Float64(y_max)),
    ]
    for rect in selection_state.rectangles
        rect[1] = pts
        rect.visible[] = true
    end
end

function _clear_shared_selection!(selection_state::SharedSelectionState)
    for rect in selection_state.rectangles
        rect[1] = [Point2f(0.0, 0.0)]
        rect.visible[] = false
    end
    selection_state.bounds[] = (0.0, 0.0)
    selection_state.visible[] = false
end

function _norm_screen_rect(fig::Figure, p1, p2; clamp_unit::Bool = false)
    w, h = size(fig.scene)
    x1, x2 = minmax(p1[1] / w, p2[1] / w)
    y1, y2 = minmax(p1[2] / h, p2[2] / h)
    clamp_unit && return (clamp(x1, 0.0, 1.0), clamp(y1, 0.0, 1.0), clamp(x2, 0.0, 1.0), clamp(y2, 0.0, 1.0))
    return (x1, y1, x2, y2)
end

function _start_figure_channel_selection!(fig::Figure, selection_state::SharedSelectionState, plot_layout, channel_selection_active)
    plot_layout.type in (:topo, :grid) || return
    selection_start = events(fig).mouseposition[]
    selection_state.bounds[] = (selection_start[1], selection_start[2])
    channel_selection_active[] = true
    selection_state.current_selection_idx = length(selection_state.selection_rectangles) + 1
end

function _update_figure_channel_selection!(fig::Figure, selection_state::SharedSelectionState, plot_layout)
    plot_layout.type in (:topo, :grid) || return
    x1, y1, x2, y2 = _norm_screen_rect(fig, selection_state.bounds[], events(fig).mouseposition[])
    rect_points = [Point2f(x1, y1), Point2f(x2, y1), Point2f(x2, y2), Point2f(x1, y2)]
    current_idx = selection_state.current_selection_idx
    if isnothing(current_idx) || current_idx > length(selection_state.selection_rectangles)
        rect = poly!(
            fig.scene,
            rect_points,
            color = (selection_state.selection_color, selection_state.selection_alpha),
            strokecolor = :red,
            strokewidth = 2,
            overdraw = true,
            space = :relative,
        )
        push!(selection_state.selection_rectangles, rect)
        push!(selection_state.selection_bounds, (x1, y1, x2, y2))
    else
        rect = selection_state.selection_rectangles[current_idx]
        rect[1] = rect_points
        selection_state.selection_bounds[current_idx] = (x1, y1, x2, y2)
    end
end

function _finish_figure_channel_selection!(
    fig::Figure,
    selection_state::SharedSelectionState,
    plot_layout,
    channel_selection_active,
    axes::Vector{Axis},
)
    (plot_layout.type in (:topo, :grid) && channel_selection_active[]) || return
    current_bounds = _norm_screen_rect(fig, selection_state.bounds[], events(fig).mouseposition[]; clamp_unit = true)
    current_idx = selection_state.current_selection_idx
    if isnothing(current_idx) || current_idx > length(selection_state.selection_bounds)
        push!(selection_state.selection_bounds, current_bounds)
    else
        selection_state.selection_bounds[current_idx] = current_bounds
    end
    channels = plot_layout isa PlotLayout ? plot_layout.channels : get(plot_layout, :channels, Symbol[])
    axes_rects = _get_axes_rectangles(axes, channels, fig)

    all_selected_channels = Symbol[]
    all_overlapping_axes_rects = Tuple{NTuple{4,Float64},Vector{Symbol}}[]
    seen_rects = Set{NTuple{4,Float64}}()

    for selection_bounds in selection_state.selection_bounds
        for (axis_rect, chs) in axes_rects
            if _rectangles_overlap(selection_bounds, axis_rect)
                append!(all_selected_channels, chs)
                if axis_rect ∉ seen_rects
                    push!(seen_rects, axis_rect)
                    push!(all_overlapping_axes_rects, (axis_rect, chs))
                end
            end
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
        append!(
            selection_state.channel_rectangles,
            _draw_channel_rectangles!(fig, all_overlapping_axes_rects, selection_state.selection_color, selection_state.selection_alpha),
        )
    end
    if !isempty(unique_channels)
        @info "Total selected channels across all regions: $unique_channels"
        @info "Number of selection areas: $(length(selection_state.selection_rectangles))"
    else
        @info "No channels selected"
    end
    channel_selection_active[] = false
end

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
    @info "Cleared all channel selections"
end

function _draw_channel_rectangles!(fig::Figure, all_overlapping_axes_rects, selection_color::Symbol = :blue, selection_alpha::Float64 = 0.6)
    rectangles = Makie.Poly[]
    alpha_val = min(1.0, selection_alpha * 2.0)
    for (axis_rect, _) in all_overlapping_axes_rects
        rect = poly!(
            fig.scene,
            [
                Point2f(axis_rect[1], axis_rect[2]),
                Point2f(axis_rect[3], axis_rect[2]),
                Point2f(axis_rect[3], axis_rect[4]),
                Point2f(axis_rect[1], axis_rect[4]),
            ],
            color = (selection_color, alpha_val),
            strokecolor = selection_color,
            strokewidth = 2,
            overdraw = true,
            space = :relative,
        )
        push!(rectangles, rect)
    end
    return rectangles
end

function _get_axes_rectangles(axes::Vector{Axis}, channels::Vector{Symbol}, fig::Figure)
    axes_rects = Tuple{NTuple{4,Float64},Vector{Symbol}}[]
    w, h = size(fig.scene)
    for (idx, ax) in enumerate(axes)
        if idx <= length(channels)
            ax_pos = ax.scene.viewport[]
            ax_x, ax_y = ax_pos.origin[1], ax_pos.origin[2]
            ax_w, ax_h = ax_pos.widths[1], ax_pos.widths[2]
            push!(axes_rects, ((ax_x / w, ax_y / h, (ax_x + ax_w) / w, (ax_y + ax_h) / h), [channels[idx]]))
        end
    end
    return axes_rects
end

function _find_active_axis(axes::Vector{Axis}, mouse_pos)
    for ax in axes
        _is_mouse_in_axis(ax, mouse_pos) && return ax
    end
    return nothing
end

function _find_active_axis_with_dataset(axes::Vector{Axis}, mouse_pos, datasets::Vector)
    for (idx, ax) in enumerate(axes)
        if _is_mouse_in_axis(ax, mouse_pos)
            return (ax, idx <= length(datasets) ? datasets[idx] : nothing)
        end
    end
    return (nothing, nothing)
end

function _setup_keyboard_tracking(fig::Figure)
    shift_pressed = Ref(false)
    ctrl_pressed = Ref(false)
    on(events(fig).keyboardbutton) do key_event
        if key_event.key in (Keyboard.left_shift, Keyboard.right_shift)
            shift_pressed[] = key_event.action in (Keyboard.press, Keyboard.repeat)
        elseif key_event.key in (Keyboard.left_control, Keyboard.right_control)
            ctrl_pressed[] = key_event.action in (Keyboard.press, Keyboard.repeat)
        end
    end
    return shift_pressed, ctrl_pressed
end

function _handle_mouse_button_events(
    fig::Figure,
    axes::Vector{Axis},
    selection_state::SharedSelectionState,
    data,
    right_click_handler,
    shift_pressed,
)
    on(events(fig).mousebutton) do event
        active_ax = _find_active_axis(axes, events(fig).mouseposition[])
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
        elseif event.button == Mouse.right && event.action == Mouse.press && !isnothing(right_click_handler)
            right_click_handler(selection_state, mouse_x, data)
        end
    end
end

function _handle_mouse_movement(fig::Figure, axes::Vector{Axis}, selection_state::SharedSelectionState)
    on(events(fig).mouseposition) do _
        if selection_state.active[]
            active_ax = _find_active_axis(axes, events(fig).mouseposition[])
            if !isnothing(active_ax)
                _update_shared_selection!(active_ax, selection_state, selection_state.bounds[][1], mouseposition(active_ax)[1])
            end
        end
    end
end

function _setup_unified_selection!(
    fig::Figure,
    axes::Vector{Axis},
    selection_state::SharedSelectionState,
    data,
    right_click_handler = nothing,
)
    shift_pressed, _ = _setup_keyboard_tracking(fig)
    _handle_mouse_button_events(fig, axes, selection_state, data, right_click_handler, shift_pressed)
    _handle_mouse_movement(fig, axes, selection_state)
end

function _setup_channel_selection_events!(
    fig::Figure,
    selection_state::SharedSelectionState,
    plot_layout,
    data,
    axes::Vector{Axis};
    channel_right_click_handler = nothing,
)
    channel_selection_active = Ref(false)
    shift_pressed, ctrl_pressed = _setup_keyboard_tracking(fig)

    on(events(fig).keyboardbutton) do key_event
        if key_event.action == Keyboard.release && key_event.key in (Keyboard.left_control, Keyboard.right_control)
            channel_selection_active[] = false
        end
    end

    on(events(fig).mousebutton) do event
        if event.button == Mouse.left && ctrl_pressed[]
            if event.action == Mouse.press
                _start_figure_channel_selection!(fig, selection_state, plot_layout, channel_selection_active)
            elseif event.action == Mouse.release && channel_selection_active[]
                _finish_figure_channel_selection!(fig, selection_state, plot_layout, channel_selection_active, axes)
            end
        elseif event.button == Mouse.left && event.action == Mouse.press && !ctrl_pressed[] && !shift_pressed[]
            _clear_all_shared_channel_selections!(fig, selection_state)
        elseif event.button == Mouse.right &&
               event.action == Mouse.press &&
               !isnothing(channel_right_click_handler) &&
               !isempty(selection_state.selection_bounds)
            selected_chs = _resolve_selected_channels(selection_state, plot_layout, axes, fig)
            !isempty(selected_chs) && channel_right_click_handler(selected_chs, data)
        end
    end

    on(events(fig).mouseposition) do _
        channel_selection_active[] && _update_figure_channel_selection!(fig, selection_state, plot_layout)
    end
end

_is_shift_held(fig::Figure) = Keyboard.left_shift in events(fig).keyboardstate || Keyboard.right_shift in events(fig).keyboardstate

function _get_mouse_position(ax::Axis)
    pos = mouseposition(ax.scene)
    return (!isnothing(pos[1]) && !isnothing(pos[2])) ? (pos[1], pos[2]) : nothing
end

@noinline _col_val_f64(col::AbstractVector, idx::Int)::Float64 = Float64(@inbounds col[idx])
@inline _closest_time_idx_typed(time_points::AbstractVector, mouse_time::Float64)::Int =
    isempty(time_points) ? 0 : find_closest_time_index(time_points, mouse_time)

function _find_closest_channel(mouse_time, mouse_amp, current_epoch, selected_channels)
    tcol = current_epoch[!, :time]
    mt = Float64(mouse_time)
    closest_time_idx = _closest_time_idx_typed(tcol, mt)
    closest_time_idx == 0 && return nothing, Inf

    lookup = DataFrames.index(current_epoch).lookup
    cols = DataFrames._columns(current_epoch)
    min_distance = Inf
    closest_channel = nothing
    mouse_amp_f64 = Float64(mouse_amp)

    for ch in selected_channels
        col_idx = get(lookup, ch, 0)
        if col_idx != 0
            col = cols[col_idx]
            val = col isa Vector{Float64} ? (@inbounds col[closest_time_idx]) : _col_val_f64(col, closest_time_idx)
            distance = abs(val - mouse_amp_f64)
            if distance < min_distance
                min_distance = distance
                closest_channel = ch
            end
        end
    end
    return closest_channel, min_distance
end

function _toggle_channel_selection(channel, selected_channels_set)
    if channel in selected_channels_set
        delete!(selected_channels_set, channel)
        @info "Deselected channel: $channel"
    else
        push!(selected_channels_set, channel)
        @info "Selected channel: $channel"
    end
end

function _rectangles_overlap(rect1, rect2)
    x1_min, y1_min, x1_max, y1_max = rect1
    x2_min, y2_min, x2_max, y2_max = rect2
    return !(x1_max < x2_min || x2_max < x1_min || y1_max < y2_min || y2_max < y1_min)
end

function _resolve_selected_channels(selection_state::SharedSelectionState, plot_layout, axes::Vector{Axis}, fig::Figure)
    channels_list = plot_layout isa PlotLayout ? plot_layout.channels : get(plot_layout, :channels, Symbol[])
    isempty(channels_list) && return Symbol[]
    axes_rects = _get_axes_rectangles(axes, channels_list, fig)
    selected = Symbol[]
    for bounds in selection_state.selection_bounds
        for (axis_rect, chs) in axes_rects
            _rectangles_overlap(bounds, axis_rect) && append!(selected, chs)
        end
    end
    return unique(selected)
end
