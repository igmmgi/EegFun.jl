using GLMakie

##################################################################
# Data Browser: Continuous Data
struct Marker
    data::Any
    line::Any
    text::Any
end

@kwdef mutable struct IcaState
    removed_activations::Union{Nothing,Any} = nothing
    components_to_remove::Union{Nothing,Int} = nothing
    components_removed::Union{Nothing,Int} = nothing
end

mutable struct ExtraChannelVis
    channel::Union{Nothing,Makie.Lines,Makie.PolyElement}
    label::Union{Nothing,Makie.Text}
end

@kwdef mutable struct SelectionState
    active::Observable{Bool} = Observable(false)
    bounds::Observable{Tuple{Float64,Float64}} = Observable((0.0, 0.0))
    visible::Observable{Bool} = Observable(false)
    rectangle::Any  # Makie.Poly - needs ax to initialize

    function SelectionState(ax; kwargs...)
        new(Observable(false), Observable((0.0, 0.0)), Observable(false), poly!(ax, Point2f[], color = (:blue, 0.3)))
    end

end

@kwdef mutable struct FilterState
    active::Observable{NamedTuple{(:hp, :lp),Tuple{Bool,Bool}}} = Observable((hp = false, lp = false))
    hp_freq::Observable{Float64} = Observable(0.1)
    lp_freq::Observable{Float64} = Observable(40.0)
end

struct ToggleConfig
    label::String
    action::Function
end

@kwdef mutable struct ViewState
    xrange::Observable{UnitRange{Int64}} = Observable(1:5000)
    yrange::Observable{UnitRange{Int64}} = Observable(-1500:1500)
    offset::Vector{Float64}
    crit_val::Observable{Float64} = Observable(0.0)
    butterfly::Observable{Bool} = Observable(false)

    function ViewState(n_channels::Int; kwargs...)
        offset = n_channels > 1 ? LinRange(1500 * 0.9, -1500 * 0.9, n_channels + 2)[2:end-1] : zeros(n_channels)
        new(Observable(1:5000), Observable(-1500:1500), offset, Observable(0.0), Observable(false))
    end

end

@kwdef mutable struct ChannelState

    labels::Vector{Symbol}
    visible::Vector{Bool}
    data_labels::Dict{Symbol,Makie.Text} = Dict{Symbol,Makie.Text}()
    data_lines::Dict{Symbol,Makie.Lines} = Dict{Symbol,Makie.Lines}()

    function ChannelState(channel_labels::Vector{Symbol}; kwargs...)
        new(channel_labels, fill(true, length(channel_labels)), Dict{Symbol,Makie.Text}(), Dict{Symbol,Makie.Lines}())
    end

end

@kwdef mutable struct ContinuousDataState
    current::Observable{DataFrame}
    original::EegData
    filter_state::FilterState = FilterState()

    function ContinuousDataState(data::EegData; kwargs...)
        new(Observable(copy(data.data)), data, FilterState())
    end
end

@kwdef mutable struct EpochedDataState
    current::Vector{Observable{DataFrame}}
    original::EegData
    filter_state::FilterState = FilterState()
    current_epoch::Observable{Int} = Observable(1)

    function EpochedDataState(data::EegData; kwargs...)
        new([Observable(copy(df)) for df in data.data], data, FilterState(), Observable(1))
    end
end


@kwdef mutable struct ContinuousDataBrowserState
    view::ViewState
    channels::ChannelState
    data::ContinuousDataState
    selection::SelectionState
    markers::Vector{Marker} = Vector{Marker}()
    ica_state::Union{Nothing,IcaState} = nothing
    extra_channel_vis::Union{Nothing,ExtraChannelVis} = nothing
end

@kwdef mutable struct EpochedDataBrowserState
    view::ViewState
    channels::ChannelState
    data::EpochedDataState
    selection::SelectionState
    markers::Vector{Marker} = Vector{Marker}()
    ica_state::Union{Nothing,IcaState} = nothing
    extra_channel_vis::Union{Nothing,ExtraChannelVis} = nothing
end



function clear_axes!(ax, datas)
    [delete!(ax, value) for data in datas for (key, value) in data]
end

function set_axes!(ax, state::ContinuousDataBrowserState)
    @lift xlims!(
        ax,
        $(state.data.current).time[$(state.view.xrange)[1]],
        $(state.data.current).time[$(state.view.xrange)[end]],
    )
    @lift ylims!(ax, $(state.view.yrange)[1], $(state.view.yrange)[end])
end

function set_axes!(ax, state::EpochedDataBrowserState)
    @lift xlims!(ax, $(state.data.current[1]).time[1], $(state.data.current[1]).time[end])
    @lift ylims!(ax, $(state.view.yrange)[1], $(state.view.yrange)[end])
end


function yless!(ax, state)
    (state.view.yrange.val[1] + 100 >= 0 || state.view.yrange.val[end] - 100 <= 0) && return
    state.view.yrange.val = state.view.yrange.val[1]+100:state.view.yrange.val[end]-100
    ylims!(ax, state.view.yrange.val[1], state.view.yrange.val[end])
end

function ymore!(ax, state)
    state.view.yrange.val = state.view.yrange.val[1]-100:state.view.yrange.val[end]+100
    ylims!(ax, state.view.yrange.val[1], state.view.yrange.val[end])
end

function xback!(ax, state::ContinuousDataBrowserState)
    state.view.xrange.val[1] - 200 < 1 && return
    state.view.xrange.val = state.view.xrange.val .- 200
    xlims!(
        ax,
        state.data.current[].time[state.view.xrange.val[1]],
        state.data.current[].time[state.view.xrange.val[end]],
    )
end

function xforward!(ax, state::ContinuousDataBrowserState)
    state.view.xrange.val[1] + 200 > nrow(state.data.current[]) && return
    state.view.xrange.val = state.view.xrange.val .+ 200
    xlims!(
        ax,
        state.data.current[].time[state.view.xrange.val[1]],
        state.data.current[].time[state.view.xrange.val[end]],
    )
end

function step_epoch_backward(ax, state::EpochedDataBrowserState)
    println("backward")
    clear_axes!(ax, [state.channels.data_lines, state.channels.data_labels])
    state.data.current_epoch[] = max(1, state.data.current_epoch[] - 1)
    #ax.title = "Epoch $(trial.val)/$(length(dat.data))"
    #menu_trial.i_selected[] = trial
    #update_extreme_spans!()
    #update_markers!(markers)
    draw(ax, state)
end

function step_epoch_forward(ax, state::EpochedDataBrowserState)
    println("forward")
    clear_axes!(ax, [state.channels.data_lines, state.channels.data_labels])
    state.data.current_epoch[] = min(n_epochs(state.data.original), state.data.current_epoch[] + 1)
    #trial[] = min(length(dat.data), trial.val[1] + 1)
    #ax.title = "Epoch $(trial.val)/$(length(dat.data))"
    #menu_trial.i_selected[] = trial
    #update_extreme_spans!()
    #update_markers!(markers)
    #draw()
    draw(ax, state)
end



function xforward!(ax, state::EpochedDataBrowserState)
    state.view.xrange.val[1] + 200 > nrow(state.data.current[state.data.current_epoch[]][]) && return
    state.view.xrange.val = state.view.xrange.val .+ 200
    xlims!(
        ax,
        state.data.current[state.data.current_epoch[]][].time[1],
        state.data.current[state.data.current_epoch[]][].time[end],
    )
end




function add_marker!(markers, ax, data, col; label = nothing, trial = nothing)
    if isnothing(trial)
        marker_data = data[findall(x -> x != 0, data[!, col]), [:time, col]]
    else
        marker_data = data[trial][findall(x -> x != 0, data[trial][!, col]), [:time, col]]
    end
    if isnothing(label)
        label = string.(marker_data[!, col])
    else
        label = repeat([label], nrow(marker_data))
    end
    push!(
        markers,
        Marker(
            marker_data,
            vlines!(marker_data.time, color = :grey, linewidth = 1, visible = false),
            text!(
                label,
                position = [(x, ax.yaxis.attributes.limits[][2] * 0.98) for x in marker_data.time],
                space = :data,
                align = (:center, :center),
                fontsize = 22,
                visible = false,
            ),
        ),
    )
end

function plot_vertical_lines!(ax, marker, active)
    marker.line.visible = active
    marker.text.visible = active
    marker.text.position = [(x, ax.yaxis.attributes.limits[][2] * 0.98) for x in marker.data.time] # incase y changed
end

function init_markers(ax, state::ContinuousDataBrowserState)
    markers = Marker[]
    add_marker!(markers, ax, state.data.current[], :triggers)
    if all(in.(["is_vEOG", "is_hEOG"], Ref(names(state.data.current[]))))
        add_marker!(markers, ax, state.data.current[], :is_vEOG, label = "v")
        add_marker!(markers, ax, state.data.current[], :is_hEOG, label = "h")
    end
    return markers
end

function init_markers(ax, state::EpochedDataBrowserState)
    markers = Marker[]
    add_marker!(markers, ax, state.data.current[state.data.current_epoch[]][], :triggers)
    if all(in.(["is_vEOG", "is_hEOG"], Ref(names(state.data.current[state.data.current_epoch[]][]))))
        add_marker!(markers, ax, state.data.current[state.data.current_epoch[]][], :is_vEOG, label = "v")
        add_marker!(markers, ax, state.data.current[state.data.current_epoch[]][], :is_hEOG, label = "h")
    end
    return markers
end

notify_data_update(state::ContinuousDataState) = notify(state.current)
notify_data_update(state::EpochedDataState) = notify(state.current[state.current_epoch[]])

# Function to plot the butterfly plot
function butterfly_plot!(state)
    state.view.butterfly[] = !state.view.butterfly[]
    if state.view.butterfly[]
        state.view.offset = zeros(length(state.channels.labels))
    else
        update_channel_offsets!(state)
    end
    for (_, label) in state.channels.data_labels
        label.visible = !state.view.butterfly[]
    end
    notify_data_update(state.data)
end

function update_channel_offsets!(state)
    nchannels = count(state.channels.visible)
    if nchannels > 1 && !state.view.butterfly[]
        state.view.offset[state.channels.visible] .=
            LinRange((state.view.yrange[][end] * 0.9), state.view.yrange[][1] * 0.9, nchannels + 2)[2:end-1]
    else
        state.view.offset[state.channels.visible] .= zeros(nchannels)
    end
end

function reset_to_original!(state::ContinuousDataState)
    state.current[] = copy(state.original.data)
end

function reset_to_original!(state::EpochedDataState)
    for (i, df) in enumerate(state.original.data)
        state.current[i][] = copy(df)
    end
end

function apply_filter!(state::ContinuousDataBrowserState, filter_type, freq)
    state.current[] = filter_data(
        state.current[],
        state.channels.labels,
        String(filter_type),
        "iir",
        freq,
        sample_rate(state.original),
        order = filter_type == :hp ? 1 : 3,
    )
end

function apply_filter!(state::EpochedDataBrowserState, filter_type, freq)
    for df in state.data.current
        df[] = filter_data(
            df[],
            state.channels.labels,
            String(filter_type),
            "iir",
            freq,
            sample_rate(state.data.original),
            order = filter_type == :hp ? 1 : 3,
        )
    end
end

function apply_filters!(state)
    # Reset to original if no filters active
    if !state.data.filter_state.active[].hp && !state.data.filter_state.active[].lp
        reset_to_original!(state.data)
        return
    end

    # Start with fresh data if changing filter configuration
    if state.data.filter_state.active[].hp != state.data.filter_state.active[].lp
        reset_to_original!(state.data)
    end

    # Apply active filters
    for (filter_type, freq) in zip([:hp, :lp], [state.data.filter_state.hp_freq[], state.data.filter_state.lp_freq[]])
        if state.data.filter_state.active[][filter_type]
            apply_filter!(state, filter_type, freq)
        end
    end

    # # Apply active filters
    # for (filter_type, freq) in zip([:hp, :lp], [state.data.filter_state.hp_freq[], state.data.filter_state.lp_freq[]])
    #     if state.data.filter_state.active[][filter_type]
    #         state.data.current[] = filter_data(
    #             state.data.current[],
    #             state.channels.labels,
    #             String(filter_type),
    #             "iir",
    #             freq,
    #             sample_rate(state.data.original),
    #             order = filter_type == :hp ? 1 : 3,
    #         )
    #     end
    # end

    notify_data_update(state.data)

end

function apply_hp_filter!(state)
    current_state = state.data.filter_state.active[]
    state.data.filter_state.active[] = (hp = !current_state.hp, lp = current_state.lp)
    apply_filters!(state)
end

function apply_lp_filter!(state)
    current_state = state.data.filter_state.active[]
    state.data.filter_state.active[] = (hp = current_state.hp, lp = !current_state.lp)
    apply_filters!(state)
end

function update_x_region_selection!(ax, state, x1, x2)
    ylims = ax.limits[][2]
    state.selection.rectangle[1] = Point2f[
        Point2f(Float64(x1), Float64(ylims[1])),
        Point2f(Float64(x2), Float64(ylims[1])),
        Point2f(Float64(x2), Float64(ylims[2])),
        Point2f(Float64(x1), Float64(ylims[2])),
    ]
end

function clear_x_region_selection!(state)
    state.selection.rectangle[1] = Point2f[]
    state.selection.bounds[] = (0.0, 0.0)
    state.selection.visible[] = false
end

function get_x_region_data(state)
    x_min, x_max = minmax(state.selection.bounds[]...)
    time_mask = (x_min .<= state.data.current[].time .<= x_max)
    selected_data = state.data.current[][time_mask, :]
    println("Selected data: $(round(x_min, digits = 2)) to $(round(x_max, digits = 2)) S, size $(size(selected_data))")
    return selected_data
end

function has_column(state::ContinuousDataState, col::String)
    col in names(state.current[])
end

function has_column(state::EpochedDataState, col::String)
    col in names(state.current[state.current_epoch[]][])
end

function create_toggles(fig, ax, state)

    configs = [
        ToggleConfig("Butterfly Plot", (active) -> butterfly_plot!(state)),
        ToggleConfig("Trigger", (active) -> plot_vertical_lines!(ax, state.markers[1], active)),
    ]

    # Add EOG toggles if available
    if has_column(state.data, "is_vEOG")
        push!(configs, ToggleConfig("vEOG", (active) -> plot_vertical_lines!(ax, state.markers[2], active)))
    end
    if has_column(state.data, "is_hEOG")
        push!(configs, ToggleConfig("hEOG", (active) -> plot_vertical_lines!(ax, state.markers[3], active)))
    end


    # Add filter toggles if available
    if state.data.original.analysis_info.hp_filter == 0.0
        push!(configs, ToggleConfig("HP-Filter On/Off", (_) -> apply_hp_filter!(state)))
    end
    if state.data.original.analysis_info.lp_filter == 0.0
        push!(configs, ToggleConfig("LP-Filter On/Off", (_) -> apply_lp_filter!(state)))
    end

    # Create toggles
    toggles = [(config.label, Toggle(fig), config.action) for config in configs]

    # Setup observers
    for toggle in toggles
        on(toggle[2].active) do active
            toggle[3](active)
        end
    end

    # Return as grid layout
    return permutedims(reduce(hcat, [[t[2], Label(fig, t[1], fontsize = 22, halign = :left)] for t in toggles]))

end

function show_menu(state)
    menu_fig = Figure()
    plot_types = ["Topoplot", "Plot2", "Plot3"]  # Could be expanded based on data type

    menu_buttons = [Button(menu_fig[idx, 1], label = plot_type) for (idx, plot_type) in enumerate(plot_types)]

    for btn in menu_buttons
        on(btn.clicks) do n
            selected_data = get_x_region_data(state)
            if btn.label[] == "Topoplot"
                plot_topoplot(selected_data, state.data_original.layout)
            elseif btn.label[] == "Plot2"
                println("Plot2: TODO")
            elseif btn.label[] == "Plot3"
                println("Plot3: TODO")
            end
        end
    end

    display(GLMakie.Screen(), menu_fig)
end

# Mouse events
function is_mouse_in_axis(ax, pos)
    bbox = ax.layoutobservables.computedbbox[]
    return bbox.origin[1] <= pos[1] <= (bbox.origin[1] + bbox.widths[1]) &&
           bbox.origin[2] <= pos[2] <= (bbox.origin[2] + bbox.widths[2])
end

function is_within_selection(state, mouse_x)
    bounds = state.selection.bounds[]
    return mouse_x >= min(bounds[1], bounds[2]) && mouse_x <= max(bounds[1], bounds[2])
end

# Selection management functions
function start_selection!(ax, state, mouse_x)
    state.selection.active[] = true
    state.selection.bounds[] = (mouse_x, mouse_x)
    update_x_region_selection!(ax, state, mouse_x, mouse_x)
end

function finish_selection!(ax, state, mouse_x)
    state.selection.active[] = false
    state.selection.visible[] = true
    state.selection.bounds[] = (state.selection.bounds[][1], mouse_x)
    update_x_region_selection!(ax, state, state.selection.bounds[][1], mouse_x)
end

# Main event handlers
function handle_mouse_events!(ax, state)
    on(events(ax).mousebutton) do event
        pos = events(ax).mouseposition[]
        if !is_mouse_in_axis(ax, pos)
            return
        end

        mouse_x = mouseposition(ax)[1]

        if event.button == Mouse.left
            handle_left_click!(ax, state, event, mouse_x)
        elseif event.button == Mouse.right && event.action == Mouse.press
            handle_right_click!(ax, state, mouse_x)
        end
    end

    # Update selection rectangle while dragging
    on(events(ax).mouseposition) do _
        if state.selection.active[]
            world_pos = mouseposition(ax)[1]
            update_x_region_selection!(ax, state, state.selection.bounds[][1], world_pos)
        end
    end
end

function handle_left_click!(ax, state, event, mouse_x)
    if event.action == Mouse.press
        if state.selection.visible[] && is_within_selection(state, mouse_x)
            clear_x_region_selection!(state)
        else
            start_selection!(ax, state, mouse_x)
        end
    elseif event.action == Mouse.release && state.selection.active[]
        finish_selection!(ax, state, mouse_x)
    end
end

function handle_right_click!(ax, state, mouse_x)
    if state.selection.visible[] && is_within_selection(state, mouse_x)
        show_menu(state)
    end
end

# Define type-specific keyboard actions
const KEYBOARD_ACTIONS =
    Dict(Keyboard.left => :left, Keyboard.right => :right, Keyboard.up => :up, Keyboard.down => :down)



function handle_keyboard_events!(fig, ax, state)
    on(events(fig).keyboardbutton) do event
        actions = get_keyboard_actions(state)
        if event.action in (Keyboard.press, Keyboard.repeat) && haskey(KEYBOARD_ACTIONS, event.key)
            action = KEYBOARD_ACTIONS[event.key]
            if state.selection.visible[]
                handle_selection_movement!(ax, state, action)
            else
                println("Navigation: $action")
                handle_navigation!(ax, state, action)
            end
        end
    end
end

function get_time_bounds(dat::ContinuousDataState)
    (dat.current[].time[1], dat.current[].time[end])
end

function get_time_bounds(dat::EpochedDataState)
    (dat.current[dat.current_epoch[]][].time[1], dat.current[dat.current_epoch[]][].time[end])
end

function handle_selection_movement!(ax, state, action::Symbol)
    width = state.selection.bounds[][2] - state.selection.bounds[][1]
    time_start, time_end = get_time_bounds(state.data)
    if action == :left
        new_start = max(time_start, state.selection.bounds[][1] - width / 5)
    elseif action == :right
        new_start = min(time_end - width, state.selection.bounds[][1] + width / 5)
    else
        return
    end
    state.selection.bounds[] = (new_start, new_start + width)
    update_x_region_selection!(ax, state, state.selection.bounds[][1], state.selection.bounds[][2])
end

function handle_navigation!(ax, state::ContinuousDataBrowserState, action::Symbol)
    if action == :left
        xback!(ax, state)
    elseif action == :right
        xforward!(ax, state)
    elseif action == :up
        ymore!(ax, state)
    elseif action == :down
        yless!(ax, state)
    end
end

function handle_navigation!(ax, state::EpochedDataBrowserState, action::Symbol)
    if action == :left
        step_epoch_backward(ax, state)
    elseif action == :right
        step_epoch_forward(ax, state)
    elseif action == :up
        ymore!(ax, state)
    elseif action == :down
        yless!(ax, state)
    end
end



function create_menu(fig, options, default, label; kwargs...)
    menu = Menu(fig, options = options, default = default, direction = :down, fontsize = 18, width = 200, kwargs...)
    return hcat(menu, Label(fig, label, fontsize = 22, halign = :left))
end

function create_labels_menu(fig, ax, state)
    options = vcat(["All", "Left", "Right", "Central"], state.channels.labels)
    menu = create_menu(fig, options, "All", "Labels")

    on(menu[1].selection) do s
        if s == "All"
            state.channels.visible .= true
        elseif s == "Left"
            state.channels.visible .= occursin.(r"\d*[13579]$", String.(state.channels.labels))
        elseif s == "Right"
            state.channels.visible .= occursin.(r"\d*[24680]$", String.(state.channels.labels))
        elseif s == "Central"
            state.channels.visible .= occursin.(r"z$", String.(state.channels.labels))
        else
            state.channels.visible .= (state.channels.labels .== Symbol(s))
        end

        clear_axes!(ax, [state.channels.data_lines, state.channels.data_labels])
        update_channel_offsets!(state)
        draw(ax, state)

    end

    return menu

end

function create_reference_menu(fig, state, dat)

    options = vcat([:none, :avg, :mastoid], state.channels.labels)
    menu = create_menu(fig, options, String(dat.analysis_info.reference), "Reference")

    on(menu[1].selection) do s
        s == :none && return
        rereference!(state.data.current[], channels(dat), resolve_reference(dat, s))
        notify(state.data.current)
    end

    return menu

end

# Update create_ica_menu to use the new struct
function create_ica_menu(fig, ax, state, ica)

    options = vcat(["None"], ica.ica_label)
    menu = create_menu(fig, options, "None", "ICA Components")

    on(menu[1].selection) do s
        clear_axes!(ax, [state.channels.channel_data_original, state.channels.channel_data_labels])
        state.ica_state.components_to_remove = extract_int(String(s))

        if !isnothing(state.ica_state.components_to_remove)
            if !isnothing(state.ica_state.components_removed)
                state.current[] = restore_original_data(
                    state.data.current[],
                    ica,
                    [state.ica_state.components_removed],
                    state.ica_state.removed_activations,
                )
            end
            state.data.current[], state.ica_state.removed_activations =
                remove_ica_components(state.data.current[], ica, [state.ica_state.components_to_remove])
            state.ica_state.components_removed = state.ica_state.components_to_remove
        else
            state.data.current[], state.ica_state.removed_activations = restore_original_data(
                state.data.current[],
                ica,
                [state.ica_state.components_removed],
                state.ica_state.removed_activations,
            )
        end

        notify(state.current)
        draw(ax, state)
    end

    return menu

end

function create_sliders(fig, state::ContinuousDataBrowserState, dat)

    sliders = []

    # Extreme value slider
    slider_extreme = Slider(fig[1, 2], range = 0:10:200, startvalue = 0, width = 100)
    on(slider_extreme.value) do x
        state.view.crit_val[] = x
    end
    push!(sliders, hcat(slider_extreme, Label(fig, @lift("Extreme: $($(slider_extreme.value)) μV"), fontsize = 22)))

    # HP filter slider
    if dat.analysis_info.hp_filter == 0.0
        slider_hp = Slider(fig[1, 2], range = 0.1:0.1:2, startvalue = 0.5, width = 100)
        on(slider_hp.value) do val
            state.data.filter_state.hp_freq[] = val
        end
        push!(sliders, hcat(slider_hp, Label(fig, @lift("HP-Filter: $($(slider_hp.value)) Hz"), fontsize = 22)))
    end

    # LP filter slider
    if dat.analysis_info.lp_filter == 0.0
        slider_lp = Slider(fig[1, 2], range = 5:5:60, startvalue = 20, width = 100)
        on(slider_lp.value) do val
            state.data.filter_state.lp_freq[] = val
        end
        push!(sliders, hcat(slider_lp, Label(fig, @lift("LP-Filter: $($(slider_lp.value)) Hz"), fontsize = 22)))
    end

    slider_range = Slider(fig[3, 1], range = 100:50:30000, startvalue = state.view.xrange[][end], snap = true)
    on(slider_range.value) do x
        new_range = slider_x.value.val:min(nrow(state.data.current[]), x + slider_x.value.val)
        if length(new_range) > 1
            state.view.xrange[] = new_range
        end
    end

    slider_x = Slider(fig[2, 1], range = 1:50:nrow(state.data.current[]), startvalue = 1, snap = true)
    on(slider_x.value) do x
        new_range = x:min(nrow(state.data.current[]), (x + slider_range.value.val) - 1)
        if length(new_range) > 1
            state.view.xrange[] = new_range
        end
    end

    return sliders

end

function create_sliders(fig, state::EpochedDataBrowserState, dat)

    sliders = []

    # Extreme value slider
    slider_extreme = Slider(fig[1, 2], range = 0:10:200, startvalue = 0, width = 100)
    on(slider_extreme.value) do x
        state.view.crit_val[] = x
    end
    push!(sliders, hcat(slider_extreme, Label(fig, @lift("Extreme: $($(slider_extreme.value)) μV"), fontsize = 22)))

    # HP filter slider
    if dat.analysis_info.hp_filter == 0.0
        slider_hp = Slider(fig[1, 2], range = 0.1:0.1:2, startvalue = 0.5, width = 100)
        on(slider_hp.value) do val
            state.data.filter_state.hp_freq[] = val
        end
        push!(sliders, hcat(slider_hp, Label(fig, @lift("HP-Filter: $($(slider_hp.value)) Hz"), fontsize = 22)))
    end

    # LP filter slider
    if dat.analysis_info.lp_filter == 0.0
        slider_lp = Slider(fig[1, 2], range = 5:5:60, startvalue = 20, width = 100)
        on(slider_lp.value) do val
            state.data.filter_state.lp_freq[] = val
        end
        push!(sliders, hcat(slider_lp, Label(fig, @lift("LP-Filter: $($(slider_lp.value)) Hz"), fontsize = 22)))
    end

    # slider_range = Slider(fig[3, 1], range = 100:50:30000, startvalue = state.view.xrange[][end], snap = true)
    # on(slider_range.value) do x
    #     new_range = slider_x.value.val:min(nrow(state.data.current[]), x + slider_x.value.val)
    #     if length(new_range) > 1
    #         state.view.xrange[] = new_range
    #     end
    # end

    # slider_x = Slider(fig[2, 1], range = 1:50:nrow(state.data.current[]), startvalue = 1, snap = true)
    # on(slider_x.value) do x
    #     new_range = x:min(nrow(state.data.current[]), (x + slider_range.value.val) - 1)
    #     if length(new_range) > 1
    #         state.view.xrange[] = new_range
    #     end
    # end

    return sliders

end



function create_extra_channel_menu(fig, ax, state, dat)

    menu = Menu(
        fig,
        options = [:none; extra_channels(dat)],
        default = "none",
        direction = :down,
        fontsize = 18,
        width = 200,
    )

    on(menu.selection) do s
        if !isnothing(state.channels.extra_channel_vis)
            !isnothing(state.channels.extra_channel_vis.channel) &&
                delete!(ax, state.channels.extra_channel_vis.channel)
            !isnothing(state.channels.extra_channel_vis.label) && delete!(ax, state.channels.extra_channel_vis.label)
        end

        s == :none && return

        if eltype(state.data.current[][!, s]) == Bool
            highlight_data = @views splitgroups(findall(state.data.current[][!, s]))
            region_offset = all(iszero, highlight_data[2] .- highlight_data[1]) ? 5 : 0

            state.channels.extra_channel_vis = ExtraChannelVis(
                vspan!(
                    ax,
                    state.data.current[][highlight_data[1], :time],
                    state.data.current[][highlight_data[2].+region_offset, :time],
                    color = "LightGrey",
                    alpha = 0.5,
                    visible = true,
                ),
                nothing,
            )
        else
            state.channels.extra_channel_vis = ExtraChannelVis(
                lines!(
                    ax,
                    @lift($(state.data.current).time),
                    @lift(begin
                        df = $(state.data.current)
                        current_offset = state.channels.offset[end] + mean(diff(state.channels.offset))
                        df[!, s] .+ current_offset
                    end),
                    color = :black,
                    linewidth = 2,
                ),
                text!(
                    ax,
                    @lift($(state.data.current).time[$(state.view.xrange)[1]]),
                    @lift(begin
                        df = $(state.data.current)
                        current_offset = state.channels.offset[end] + mean(diff(state.channels.offset))
                        df[!, s][$(state.view.xrange)[1]] .+ current_offset
                    end),
                    text = String(s),
                    align = (:left, :center),
                    fontsize = 18,
                ),
            )
        end
    end

    return hcat(menu, Label(fig, "Extra Channels", fontsize = 22, halign = :left))

end

function build_grid_components!(
    fig,
    dat,
    state,
    toggles,
    labels_menu,
    reference_menu,
    ica_menu = nothing,
    extra_menu = nothing,
)
    # Start with common components
    grid_components = Matrix{Any}[]  # Explicitly type as Matrix{Any}
    push!(grid_components, toggles[:, 1:2])

    # Add common controls that exist
    append!(grid_components, create_sliders(fig, state, dat))
    push!(grid_components, labels_menu)
    push!(grid_components, reference_menu)

    # Only add optional menus if they exist
    !isnothing(ica_menu) && push!(grid_components, ica_menu)
    !isnothing(extra_menu) && push!(grid_components, extra_menu)

    fig[1, 2] = grid!(vcat(grid_components...), tellheight = false)
    colsize!(fig.layout, 2, Relative(1 / 6))

end

# Do the actual drawing
function draw(ax, state::ContinuousDataBrowserState)
    for (idx, col) in enumerate(state.channels.visible)
        if state.channels.visible[idx]  # Only plot if channel is visible
            col = state.channels.labels[idx]
            state.channels.data_lines[col] = lines!(
                ax,
                @lift($(state.data.current).time),
                @lift($(state.data.current)[!, $col] .+ state.view.offset[idx]),
                color = @lift(abs.($(state.data.current)[!, col]) .>= $(state.view.crit_val)),
                colormap = [:darkgrey, :darkgrey, :red],
                linewidth = 2,
            )
            if !state.view.butterfly[]
                state.channels.data_labels[col] = text!(
                    ax,
                    @lift($(state.data.current)[$(state.view.xrange), :time][1]),
                    @lift($(state.data.current)[$(state.view.xrange), $col][1] .+ state.view.offset[idx]),
                    text = String(col),
                    align = (:left, :center),
                    fontsize = 18,
                )
            end
        end
    end
end

function draw(ax, state::EpochedDataBrowserState)
    current_epoch = state.data.current_epoch[]  # Get current epoch value
    current_data = state.data.current[current_epoch]  # Get current DataFrame Observable

    # Ensure view range is within data bounds
    n_samples = @lift(nrow($(current_data)))
    state.view.xrange[] = 1:min(state.view.xrange[][end], n_samples[])

    for (idx, col) in enumerate(state.channels.visible)
        if state.channels.visible[idx]  # Only plot if channel is visible
            col = state.channels.labels[idx]

            state.channels.data_lines[col] = lines!(
                ax,
                @lift($(current_data).time),
                @lift($(current_data)[!, $col] .+ state.view.offset[idx]),
                color = @lift(abs.($(current_data)[!, $col]) .>= $(state.view.crit_val)),
                colormap = [:darkgrey, :darkgrey, :red],
                linewidth = 2,
            )
            if !state.view.butterfly[]
                state.channels.data_labels[col] = text!(
                    ax,
                    @lift($(current_data).time[1]),
                    @lift($(current_data)[!, $col][1] .+ state.view.offset[idx]),
                    text = String(col),
                    align = (:left, :center),
                    fontsize = 18,
                )
            end
        end
    end
end








function plot_databrowser(dat::ContinuousData, channel_labels::Vector{Symbol}, ica::Union{InfoIca,Nothing} = nothing)

    # Setup figure and axis first
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel = "Time (S)", ylabel = "Amplitude (mV)")

    state = ContinuousDataBrowserState(
        view = ViewState(length(channel_labels)),
        channels = ChannelState(channel_labels),
        data = ContinuousDataState(dat),
        selection = SelectionState(ax),
        ica_state = !isnothing(ica) ? IcaState(nothing, nothing, nothing) : nothing,
    )

    # Controls
    deregister_interaction!(ax, :rectanglezoom) # interactions(ax)

    # axes
    set_axes!(ax, state)

    # Mouse events
    handle_mouse_events!(ax, state)
    handle_keyboard_events!(fig, ax, state)

    toggles = create_toggles(fig, ax, state)

    # Vertical line markers
    state.markers = init_markers(ax, state)

    # Create menus
    labels_menu = create_labels_menu(fig, ax, state)
    reference_menu = create_reference_menu(fig, state, dat)
    ica_menu = !isnothing(ica) ? create_ica_menu(fig, ax, state, ica) : nothing
    extra_menu = create_extra_channel_menu(fig, ax, state, dat)

    # GUI Control Panel
    build_grid_components!(fig, dat, state, toggles, labels_menu, reference_menu, ica_menu, extra_menu)

    # plot theme adjustments
    update_theme!(Theme(fontsize = 24))

    hideydecorations!(ax, label = true)
    draw(ax, state)
    display(fig)
    # DataInspector(fig)

end

# Convenience methods
plot_databrowser(dat::EegData) = plot_databrowser(dat, dat.layout.label)
plot_databrowser(dat::EegData, ica::InfoIca) = plot_databrowser(dat, dat.layout.label, ica)
plot_databrowser(dat::EegData, channel_label::Symbol) = plot_databrowser(dat, [channel_label])
plot_databrowser(dat::EegData, channel_label::Symbol, ica::InfoIca) = plot_databrowser(dat, [channel_label], ica)


function plot_databrowser(dat::EpochData, channel_labels::Vector{Symbol}, ica::Union{InfoIca,Nothing} = nothing)

    # Setup figure and axis first
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel = "Time (S)", ylabel = "Amplitude (mV)", title = "Epoch 1/$(n_epochs(dat))")

    state = EpochedDataBrowserState(
        view = ViewState(length(channel_labels)),
        channels = ChannelState(channel_labels),
        data = EpochedDataState(dat),
        selection = SelectionState(ax),
        ica_state = !isnothing(ica) ? IcaState(nothing, nothing, nothing) : nothing,
    )

    # Controls
    deregister_interaction!(ax, :rectanglezoom) # interactions(ax)

    # axes
    set_axes!(ax, state)

    # Mouse events
    handle_mouse_events!(ax, state)
    handle_keyboard_events!(fig, ax, state)

    toggles = create_toggles(fig, ax, state)

    # Vertical line markers
    state.markers = init_markers(ax, state)

    # Create menus
    labels_menu = create_labels_menu(fig, ax, state)
    reference_menu = create_reference_menu(fig, state, dat)
    ica_menu = !isnothing(ica) ? create_ica_menu(fig, ax, state, ica) : nothing
    extra_menu = create_extra_channel_menu(fig, ax, state, dat)

    # GUI Control Panel
    build_grid_components!(fig, dat, state, toggles, labels_menu, reference_menu, ica_menu, extra_menu)

    # plot theme adjustments
    update_theme!(Theme(fontsize = 24))

    hideydecorations!(ax, label = true)
    draw(ax, state)
    display(fig)
    # DataInspector(fig)

end

# Convenience methods
plot_databrowser(dat::EpochData) = plot_databrowser(dat, dat.layout.label)
plot_databrowser(dat::EpochData, ica::InfoIca) = plot_databrowser(dat, dat.layout.label, ica)
plot_databrowser(dat::EpochData, channel_label::Symbol) = plot_databrowser(dat, [channel_label])
plot_databrowser(dat::EpochData, channel_label::Symbol, ica::InfoIca) = plot_databrowser(dat, [channel_label], ica)





###########################################################
# function plot_databrowser(dat::EpochData, channel_labels::Vector{Symbol}, ica::Union{InfoIca,Nothing} = nothing)
 
# 

# 
# 
#     function update_markers!(markers)
#         for marker in markers
#             delete!(ax, marker.line)
#             delete!(ax, marker.text)
#         end
#         empty!(markers)
#         add_marker!(markers, ax, data, :triggers, trial = trial.val)
#         if ("is_vEOG" in names(dat.data[trial.val]) && "is_hEOG" in names(dat.data[trial.val]))
#             add_marker!(markers, ax, data, :is_vEOG, trial = trial.val, label = "v")
#             add_marker!(markers, ax, data, :is_hEOG, trial = trial.val, label = "h")
#         end
#     end
# 
# 
#     function toggle_button_group(fig, labels)
#         # Pre-define all possible toggles with their functions
#         toggle_configs = [
#             ("Butterfly Plot", butterfly_plot),
#             ("Trigger", plot_lines, "triggers"),
#             ("vEOG", plot_lines, "is_vEOG"),
#             ("hEOG", plot_lines, "is_hEOG"),
#             ("LP-Filter On/Off", apply_lp_filter),
#         ]
# 
#         # Filter toggles based on available labels
#         toggles = [
#             (config[1:2]..., Toggle(fig, active = false)) for
#             config in toggle_configs if length(config) == 2 || config[3] in labels
#         ]
# 
#         # Create grid layout
#         grid = [[toggle[3], Label(fig, toggle[1], fontsize = 22, halign = :left), toggle[2]] for toggle in toggles]
# 
#         return permutedims(reduce(hcat, grid))  # Return as a proper matrix
#     end
# 

# 
#     # Makie Figure
#     fig = Figure()
#     ax = Axis(fig[1, 1])
# 
#     # controls
#     # interactions(ax)
#     # deregister_interaction!(ax, :rectanglezoom)
#     # deregister_interaction!(ax, :dragpan)
#     # deregister_interaction!(ax, :scrollzoom)
#     # deregister_interaction!(ax, :limitreset)
# 
#     # data to plot
#     data = deepcopy(dat.data)
#     data_filtered = nothing
# 
#     channel_data_original = Dict()
#     channel_data_labels = Dict()
# 
#     # default xrange/yrange
#     xlimit = nrow(dat.data[1])
#     xrange = Observable(1:xlimit)
#     trial = Observable(1) # first trial
#     yrange = Observable(-1500:1500)
#     nchannels = length(channel_labels)
#     channel_labels_original = channel_labels
# 
#     if nchannels > 1
#         offset = LinRange((yrange.val[end] * 0.9), yrange.val[1] * 0.9, nchannels + 2)[2:end-1]
#     else # just centre
#         offset = zeros(length(channel_labels))
#     end
# 
# 
#     xlims!(ax, data[1].time[xrange.val[1]], data[1].time[xrange.val[end]])
#     ylims!(ax, yrange.val[1], yrange.val[end])
#     ax.title = "Epoch $(trial.val)/$(length(dat.data))"
#     ax.xlabel = "Time (S)"
#     ax.ylabel = "Amplitude (mV)"
# 
#     # toggle buttons for showing events (triggers, vEOG/hEOG, extreme values ...)
#     toggles = toggle_button_group(fig, names(data[trial.val]))
#     for t = 1:length(toggles[:, 1])
#         if toggles[t, 2].text.val ∈ ["Trigger", "vEOG", "hEOG"]
#             on(toggles[t, 1].active) do _
#                 toggles[t, 3](ax, state.markers[t-1], toggles[t, 1].active.val)
#             end
#         else
#             on(toggles[t, 1].active) do _
#                 toggles[t, 3](toggles[t, 1].active.val)
#             end
#         end
#     end
# 

# 
#     menu_trial = hcat(
#         Menu(fig, options = 1:length(data), default = 1, direction = :down, fontsize = 18),
#         Label(fig, "Epoch", fontsize = 22, halign = :left),
#     )
#     on(menu_trial[1].selection) do s
#         clear_axes(ax, [channel_data_original, channel_data_labels])
#         trial[] = s
#         ax.title = "Epoch $(trial.val)/$(length(dat.data))"
#         update_extreme_spans!()
#         update_markers!(markers)
#         draw()
#     end
# 
#     slider_extreme   = Slider(fig[1, 2], range = 0:5:100, startvalue = 200, width = 100)
#     slider_lp_filter = Slider(fig[1, 2], range = 5:5:60, startvalue = 20, width = 100)
# 
#     state.crit_val = lift(slider_extreme.value) do x
#         state.crit_val[] = x
#     end
# 

# 
#     # position GUI controls
#     fig[1, 2] = grid!(
#         vcat(
#             toggles[:, 1:2],
#             hcat(
#                 slider_lp_filter,
#                 Label(fig, @lift("LP-Filter: $($(slider_lp_filter.value)) Hz"), fontsize = 22, halign = :left),
#             ),
#             hcat(slider_extreme, Label(fig, @lift("Extreme: $($(slider_extreme.value)) μV"), fontsize = 22)),
#             menu,
#             menu_trial,
#         ),
#         tellheight = false,
#     )
#     colsize!(fig.layout, 2, Relative(1 / 6))
# 
# 
#     ################### vertical line markers ###############################
#     # Vertical line markers
#     markers = []
#     update_markers!(markers)
# 
#     extreme_spans = []
#     function update_extreme_spans!()
#         if length(extreme_spans) > 0
#             delete!(ax, extreme_spans[1])
#             extreme_spans = []
#         end
#         tmp = findall(x -> x != 0, data[trial.val][!, :].is_extreme)
#         if length(tmp) > 0
#             extreme = splitgroups(tmp)
#             if length(extreme) > 0
#                 # TODO: hard coded Toggle index!!!
#                 push!(
#                     extreme_spans,
#                     vspan!(
#                         ax,
#                         data[trial.val][extreme[1], :time],
#                         data[trial.val][extreme[2], :time],
#                         color = "LightGrey",
#                         alpha = 0.5,
#                         visible = toggles[5, 1].active.val,
#                     ),
#                 )
#             end
#         end
#     end
# 
#     #################### Extreme Values ###############################
#     #if ("is_extreme" in names(data[1]))
#     #    update_extreme_spans!()
#     #end
