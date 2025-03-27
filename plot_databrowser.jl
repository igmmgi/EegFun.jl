using GLMakie

##################################################################
# Data Browser: Continuous Data
struct Marker
    data::Any
    line::Any
    text::Any
end

mutable struct IcaState
    removed_activations::Union{Nothing,Any}
    components_to_remove::Union{Nothing,Int}
    components_removed::Union{Nothing,Int}
end

mutable struct ExtraChannelVis
    channel::Union{Nothing,Any}  #  vspan or lines
    label::Union{Nothing,Any}    # label for non-bool channels
end

mutable struct SelectionState
    active::Observable{Bool}
    bounds::Observable{Tuple{Float64,Float64}}
    visible::Observable{Bool}
    rectangle::Any  # Makie.Poly
end

mutable struct FilterState
    active::Observable{NamedTuple{(:hp, :lp),Tuple{Bool,Bool}}}
    hp_freq::Observable{Float64}
    lp_freq::Observable{Float64}
end

struct ToggleConfig
    label::String
    action::Function
    # marker_index::Union{Nothing, Int}  # Remove default value
end

mutable struct DataBrowserState
    # Data and display state
    xrange::Observable{UnitRange{Int64}}
    yrange::Observable{UnitRange{Int64}}
    filter_state::FilterState
    offset::Vector{Float64}
    crit_val::Observable{Float64}

    # Channel information
    channel_labels::Vector{Symbol}  # Current channel labels

    # Plot elements
    channel_data_labels::Dict{Symbol,Makie.Text}
    channel_data_original::Dict{Symbol,Makie.Lines}

    # Data
    data_obs::Observable{DataFrame}  # Observable data
    data_original::EegData  # Original data

    # Selection state
    selection::SelectionState

    markers::Vector{Marker}
    ica_state::Union{Nothing,IcaState}
    extra_channel_vis::Union{Nothing,ExtraChannelVis}

    channel_visible::Vector{Bool}  # Boolean mask for which channels to show
end

function clear_axes!(ax, datas)
    [delete!(ax, value) for data in datas for (key, value) in data]
end

function set_axes!(ax, state)
    @lift xlims!(ax, $(state.data_obs).time[$(state.xrange)[1]], $(state.data_obs).time[$(state.xrange)[end]])
    @lift ylims!(ax, $(state.yrange)[1], $(state.yrange)[end])
end

function yless!(ax, state)
    (state.yrange.val[1] + 100 >= 0 || state.yrange.val[end] - 100 <= 0) && return
    state.yrange.val = state.yrange.val[1]+100:state.yrange.val[end]-100
    ylims!(ax, state.yrange.val[1], state.yrange.val[end])
end

function ymore!(ax, state)
    state.yrange.val = state.yrange.val[1]-100:state.yrange.val[end]+100
    ylims!(ax, state.yrange.val[1], state.yrange.val[end])
end

function xback!(ax, state)
    state.xrange.val[1] - 200 < 1 && return
    state.xrange.val = state.xrange.val .- 200
    xlims!(ax, state.data_obs[].time[state.xrange.val[1]], state.data_obs[].time[state.xrange.val[end]])
end

function xforward!(ax, state)
    state.xrange.val[1] + 200 > nrow(state.data_obs[]) && return
    state.xrange.val = state.xrange.val .+ 200
    xlims!(ax, state.data_obs[].time[state.xrange.val[1]], state.data_obs[].time[state.xrange.val[end]])
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

function init_markers(ax, state)
    markers = Marker[]
    add_marker!(markers, ax, state.data_obs[], :triggers)
    if all(in.(["is_vEOG", "is_hEOG"], Ref(names(state.data_obs[]))))
        add_marker!(markers, ax, state.data_obs[], :is_vEOG, label = "v")
        add_marker!(markers, ax, state.data_obs[], :is_hEOG, label = "h")
    end
    return markers
end

# Function to plot the butterfly plot
function butterfly_plot!(active::Bool, state::DataBrowserState)
    if active
        state.offset = LinRange(0, 0, length(state.channel_labels))
    else
        state.offset = LinRange((state.yrange.val[end] * 0.9), state.yrange.val[1] * 0.9, length(state.channel_labels))
    end
    for (_, label) in state.channel_data_labels
        label.visible = !active
    end
    notify(state.data_obs)
end

function update_channel_offsets!(state)
    n_visible = count(state.channel_visible)  # Count number of true values
    state.offset[state.channel_visible] = calculate_channel_offsets(n_visible, state.yrange[])
end

function apply_filters!(state::DataBrowserState)

    # If turning off both filters, reset to original
    if !state.filter_state.active[].hp && !state.filter_state.active[].lp
        state.data_obs[] = copy(state.data_original.data)
        return
    end

    # If turning off one filter, start fresh and apply remaining filter
    if state.filter_state.active[].hp != state.filter_state.active[].lp  # Only one filter active
        state.data_obs[] = copy(state.data_original.data)
        if state.filter_state.active[].hp
            state.data_obs[] = filter_data(
                state.data_obs[],
                state.channel_labels,
                "hp",
                "iir",
                state.filter_state.hp_freq[],
                sample_rate(state.data_original),
                order = 1,
            )
        else  # lp is active
            state.data_obs[] = filter_data(
                state.data_obs[],
                state.channel_labels,
                "lp",
                "iir",
                state.filter_state.lp_freq[],
                sample_rate(state.data_original),
                order = 3,
            )
        end
    else  # Both filters active - only apply the new filter to current data
        if state.filter_state.active[].lp  # LP was just turned on
            state.data_obs[] = filter_data(
                state.data_obs[],
                state.channel_labels,
                "lp",
                "iir",
                state.filter_state.lp_freq[],
                sample_rate(state.data_original),
                order = 3,
            )
        else  # HP was just turned on
            state.data_obs[] = filter_data(
                state.data_obs[],
                state.channel_labels,
                "hp",
                "iir",
                state.filter_state.hp_freq[],
                sample_rate(state.data_original),
                order = 1,
            )
        end
    end

    notify(state.data_obs)

end

function apply_hp_filter!(state)
    current_state = state.filter_state.active[]
    state.filter_state.active[] = (hp = !current_state.hp, lp = current_state.lp)
    apply_filters!(state)
end

function apply_lp_filter!(state)
    current_state = state.filter_state.active[]
    state.filter_state.active[] = (hp = current_state.hp, lp = !current_state.lp)
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
    time_mask = (x_min .<= state.data_obs[].time .<= x_max)
    selected_data = state.data_obs[][time_mask, :]
    println("Selected data: $(round(x_min, digits = 2)) to $(round(x_max, digits = 2)) S, size $(size(selected_data))")
    return selected_data
end

function calculate_channel_offsets(nchannels::Int, yrange::UnitRange{Int64})
    if nchannels > 1
        return LinRange((yrange[end] * 0.9), yrange[1] * 0.9, nchannels + 2)[2:end-1]
    else
        return zeros(nchannels)
    end
end

function create_toggles(fig, ax, state)

    configs = [
        ToggleConfig("Butterfly Plot", (active) -> butterfly_plot!(active, state)),
        ToggleConfig("Trigger", (active) -> plot_vertical_lines!(ax, state.markers[1], active))
    ]

    # Add EOG toggles if available
    if "is_vEOG" in names(state.data_obs[])
        push!(configs, ToggleConfig("vEOG", (active) -> plot_vertical_lines!(ax, state.markers[2], active)))
    end
    if "is_hEOG" in names(state.data_obs[])
        push!(configs, ToggleConfig("hEOG", (active) -> plot_vertical_lines!(ax, state.markers[3], active)))
    end
    
    # Add filter toggles if available
    if state.data_original.analysis_info.hp_filter == 0.0
        push!(configs, ToggleConfig("HP-Filter On/Off", (_) -> apply_hp_filter!(state)))
    end
    if state.data_original.analysis_info.lp_filter == 0.0
        push!(configs, ToggleConfig("LP-Filter On/Off", (_) -> apply_lp_filter!(state)))
    end

    # Create toggles
    toggles = [(
        config.label,
        Toggle(fig),
        config.action
    ) for config in configs]

    # Setup observers
    for toggle in toggles
        on(toggle[2].active) do active
            toggle[3](active)
        end
    end

    # Return as grid layout
    return permutedims(reduce(hcat, [[t[2], Label(fig, t[1], fontsize = 22, halign = :left)] for t in toggles]))

end


function show_menu()

    menu_fig = Figure()
    menu_buttons = [
        Button(menu_fig[idx, 1], label = plot_type) for
        (idx, plot_type) in enumerate(["Topoplot", "Topoplot", "Topoplot"])
    ]
    for btn in menu_buttons
        on(btn.clicks) do n
            selected_data = get_x_region_data(state)
            if btn.label[] == "Topoplot"
                plot_topoplot(selected_data, dat.layout)
            end
        end
    end

    display(GLMakie.Screen(), menu_fig)

end

function handle_mouse_events!(ax, state)

    on(events(ax).mousebutton) do event
        # Check if mouse is within axis bounds
        pos = events(ax).mouseposition[]
        bbox = ax.layoutobservables.computedbbox[]

        # Check if position is within bounding box
        if bbox.origin[1] <= pos[1] <= (bbox.origin[1] + bbox.widths[1]) &&
           bbox.origin[2] <= pos[2] <= (bbox.origin[2] + bbox.widths[2])
            mouse_x = mouseposition(ax)[1]

            if event.button == Mouse.left
                if event.action == Mouse.press
                    if state.selection.visible[] &&
                       mouse_x >= min(state.selection.bounds[][1], state.selection.bounds[][2]) &&
                       mouse_x <= max(state.selection.bounds[][1], state.selection.bounds[][2])
                        clear_x_region_selection!(state)
                    else
                        state.selection.active[] = true
                        state.selection.bounds[] = (mouse_x, mouse_x)
                        update_x_region_selection!(ax, state, mouse_x, mouse_x)
                    end
                elseif event.action == Mouse.release && state.selection.active[]
                    state.selection.active[] = false
                    state.selection.visible[] = true
                    state.selection.bounds[] = (state.selection.bounds[][1], mouse_x)
                    update_x_region_selection!(ax, state, state.selection.bounds[][1], mouse_x)
                end
            elseif event.button == Mouse.right && event.action == Mouse.press
                if state.selection.visible[] &&
                   mouse_x >= min(state.selection.bounds[][1], state.selection.bounds[][2]) &&
                   mouse_x <= max(state.selection.bounds[][1], state.selection.bounds[][2])
                    show_menu()
                end
            end
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

# keyboard events
function handle_keyboard_events!(fig, ax, state)

    on(events(fig).keyboardbutton) do event
        if event.action in (Keyboard.press, Keyboard.repeat)
            if state.selection.visible[]
                # Move selection region
                if event.key == Keyboard.left
                    width = state.selection.bounds[][2] - state.selection.bounds[][1]
                    new_start = max(state.data_obs[].time[1], state.selection.bounds[][1] - width / 5)
                    state.selection.bounds[] = (new_start, new_start + width)
                    update_x_region_selection!(ax, state, state.selection.bounds[][1], state.selection.bounds[][2])
                elseif event.key == Keyboard.right
                    width = state.selection.bounds[][2] - state.selection.bounds[][1]
                    new_start = min(state.data_obs[].time[end] - width, state.selection.bounds[][1] + width / 5)
                    state.selection.bounds[] = (new_start, new_start + width)
                    update_x_region_selection!(ax, state, state.selection.bounds[][1], state.selection.bounds[][2])
                end
            else
                event.key == Keyboard.left && xback!(ax, state)
                event.key == Keyboard.right && xforward!(ax, state)
            end
            event.key == Keyboard.down && yless!(ax, state)
            event.key == Keyboard.up && ymore!(ax, state)
        end
    end

end


function create_labels_menu(fig, ax, state)

    menu = Menu(
        fig,
        options = vcat(["All", "Left", "Right", "Central"], state.channel_labels),
        default = "All",
        direction = :down,
        fontsize = 18,
        width = 200,
    )

    on(menu.selection) do s
        if s == "All"
            state.channel_visible .= true
        elseif s == "Left"
            state.channel_visible .= occursin.(r"\d*[13579]$", String.(state.channel_labels))
        elseif s == "Right"
            state.channel_visible .= occursin.(r"\d*[24680]$", String.(state.channel_labels))
        elseif s == "Central"
            state.channel_visible .= occursin.(r"z$", String.(state.channel_labels))
        else
            state.channel_visible .= (state.channel_labels .== Symbol(s))
        end

        clear_axes!(ax, [state.channel_data_original, state.channel_data_labels])
        update_channel_offsets!(state)
        draw(ax, state)

    end

    return hcat(menu, Label(fig, "Labels", fontsize = 22, halign = :left))

end

function create_reference_menu(fig, state, dat)

    menu = Menu(
        fig,
        options = vcat([:none, :avg, :mastoid], state.channel_labels),
        default = String(dat.analysis_info.reference),
        direction = :down,
        fontsize = 18,
        width = 200,
    )

    on(menu.selection) do s
        s == :none && return
        rereference!(state.data_obs[], channels(dat), resolve_reference(dat, s))
        notify(state.data_obs)
    end

    return hcat(menu, Label(fig, "Reference", fontsize = 22, halign = :left))

end


# Update create_ica_menu to use the new struct
function create_ica_menu(fig, ax, state, ica)
    menu = Menu(fig, options = vcat(["None"], ica.ica_label), default = "None", direction = :down, fontsize = 18)

    on(menu.selection) do s
        clear_axes!(ax, [state.channel_data_original, state.channel_data_labels])
        state.ica_state.components_to_remove = extract_int(String(s))

        if !isnothing(state.ica_state.components_to_remove)
            if !isnothing(state.ica_state.components_removed)
                state.data_obs[] = restore_original_data(
                    state.data_obs[],
                    ica,
                    [state.ica_state.components_removed],
                    state.ica_state.removed_activations,
                )
            end
            state.data_obs[], state.ica_state.removed_activations =
                remove_ica_components(state.data_obs[], ica, [state.ica_state.components_to_remove])
            state.ica_state.components_removed = state.ica_state.components_to_remove
        else
            state.data_obs[] = restore_original_data(
                state.data_obs[],
                ica,
                [state.ica_state.components_removed],
                state.ica_state.removed_activations,
            )
        end

        notify(state.data_obs)
        draw(ax, state)
    end

    return hcat(menu, Label(fig, "ICA Components", fontsize = 22, halign = :left))

end

function create_sliders(fig, state, dat)

    sliders = []

    # Extreme value slider
    slider_extreme = Slider(fig[1, 2], range = 0:10:200, startvalue = 0, width = 100)
    on(slider_extreme.value) do x
        state.crit_val[] = x
    end
    push!(sliders, hcat(slider_extreme, Label(fig, @lift("Extreme: $($(slider_extreme.value)) μV"), fontsize = 22)))

    # HP filter slider
    if dat.analysis_info.hp_filter == 0.0
        slider_hp = Slider(fig[1, 2], range = 0.1:0.1:2, startvalue = 0.5, width = 100)
        on(slider_hp.value) do val
            state.filter_state.hp_freq[] = val
        end
        push!(sliders, hcat(slider_hp, Label(fig, @lift("HP-Filter: $($(slider_hp.value)) Hz"), fontsize = 22)))
    end

    # LP filter slider
    if dat.analysis_info.lp_filter == 0.0
        slider_lp = Slider(fig[1, 2], range = 5:5:60, startvalue = 20, width = 100)
        on(slider_lp.value) do val
            state.filter_state.lp_freq[] = val
        end
        push!(sliders, hcat(slider_lp, Label(fig, @lift("LP-Filter: $($(slider_lp.value)) Hz"), fontsize = 22)))
    end

    slider_range = Slider(fig[3, 1], range = 100:50:30000, startvalue = state.xrange[][end], snap = true)
    on(slider_range.value) do x
        new_range = slider_x.value.val:min(nrow(state.data_obs[]), x + slider_x.value.val)
        if length(new_range) > 1
            state.xrange[] = new_range
        end
    end

    slider_x = Slider(fig[2, 1], range = 1:50:nrow(state.data_obs[]), startvalue = 1, snap = true)
    on(slider_x.value) do x
        new_range = x:min(nrow(state.data_obs[]), (x + slider_range.value.val) - 1)
        if length(new_range) > 1
            state.xrange[] = new_range
        end
    end

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
        if !isnothing(state.extra_channel_vis)
            !isnothing(state.extra_channel_vis.channel) && delete!(ax, state.extra_channel_vis.channel)
            !isnothing(state.extra_channel_vis.label) && delete!(ax, state.extra_channel_vis.label)
        end

        s == :none && return

        if eltype(state.data_obs[][!, s]) == Bool
            highlight_data = @views splitgroups(findall(state.data_obs[][!, s]))
            region_offset = all(iszero, highlight_data[2] .- highlight_data[1]) ? 5 : 0

            state.extra_channel_vis = ExtraChannelVis(
                vspan!(
                    ax,
                    state.data_obs[][highlight_data[1], :time],
                    state.data_obs[][highlight_data[2].+region_offset, :time],
                    color = "LightGrey",
                    alpha = 0.5,
                    visible = true,
                ),
                nothing,
            )
        else
            state.extra_channel_vis = ExtraChannelVis(
                lines!(
                    ax,
                    @lift($(state.data_obs).time),
                    @lift(begin
                        df = $(state.data_obs)
                        current_offset = $(state.offset)[end] + mean(diff($(state.offset)))
                        df[!, s] .+ current_offset
                    end),
                    color = :black,
                    linewidth = 2,
                ),
                text!(
                    ax,
                    @lift($(state.data_obs).time[$(state.xrange)[1]]),
                    @lift(begin
                        df = $(state.data_obs)
                        current_offset = $(state.offset)[end] + mean(diff($(state.offset)))
                        df[!, s][$(state.xrange)[1]] .+ current_offset
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
    # Only add menu_ica if it exists
    if !isnothing(ica_menu)
        push!(grid_components, ica_menu)
    end
    push!(grid_components, extra_menu)

    fig[1, 2] = grid!(vcat(grid_components...), tellheight = false)
    colsize!(fig.layout, 2, Relative(1 / 6))

end

# Do the actuall drawing
function draw(ax, state; plot_labels = true)
    for (idx, col) in enumerate(state.channel_visible)
        if state.channel_visible[idx]  # Only plot if channel is visible
            col = state.channel_labels[idx]
            state.channel_data_original[col] = lines!(
                ax,
                @lift($(state.data_obs).time),
                @lift($(state.data_obs)[!, $col] .+ state.offset[idx]),
                color = @lift(abs.($(state.data_obs)[!, col]) .>= $(state.crit_val)),
                colormap = [:darkgrey, :darkgrey, :red],
                linewidth = 2,
            )
            if plot_labels
                state.channel_data_labels[col] = text!(
                    ax,
                    @lift($(state.data_obs)[$(state.xrange), :time][1]),
                    @lift($(state.data_obs)[$(state.xrange), $col][1] .+ state.offset[idx]),
                    text = String(col),
                    align = (:left, :center),
                    fontsize = 18,
                )
            end
        end
    end
end

function plot_databrowser(dat::ContinuousData, channel_labels::Vector{Symbol}, ica::Union{InfoIca, Nothing} = nothing)

    # Setup figure and axis first
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel = "Time (S)", ylabel = "Amplitude (mV)")

    state = DataBrowserState(
        Observable(1:5000),           # xrange
        Observable(-1500:1500),       # yrange
        FilterState(                  # New filter state struct
            Observable((hp = false, lp = false)),
            Observable(0.1),          # hp_freq
            Observable(40.0),          # lp_freq
        ),
        length(channel_labels) > 1 ? LinRange((1500 * 0.9), -1500 * 0.9, length(channel_labels) + 2)[2:end-1] :
        Observable(zeros(length(channel_labels))),  # offset
        Observable(0.0),
        channel_labels,
        Dict{Symbol,Makie.Text}(),
        Dict{Symbol,Makie.Lines}(),
        Observable(copy(dat.data)),
        dat,
        SelectionState(
            Observable(false),
            Observable((0.0, 0.0)),
            Observable(false),
            poly!(ax, Point2f[], color = (:blue, 0.3)),
        ),
        Vector{Marker}(),  # Initialize empty markers vecto
        !isnothing(ica) ? IcaState(nothing, nothing, nothing) : nothing,
        nothing,
        fill(true, length(channel_labels)),  # Initially show all channels
    )

    set_axes!(ax, state)

    # Controls
    deregister_interaction!(ax, :rectanglezoom) # interactions(ax)

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
    draw(ax, state; plot_labels = true)
    display(fig)
    # DataInspector(fig)

end

# Convenience methods
plot_databrowser(dat::ContinuousData) = plot_databrowser(dat, dat.layout.label)
plot_databrowser(dat::ContinuousData, ica::InfoIca) = plot_databrowser(dat, dat.layout.label, ica)
plot_databrowser(dat::ContinuousData, channel_label::Symbol) = plot_databrowser(dat, [channel_label])
plot_databrowser(dat::ContinuousData, channel_label::Symbol, ica::InfoIca) = plot_databrowser(dat, [channel_label], ica)




# ###########################################################
# function plot_databrowser(dat::EpochData, channel_labels::Vector{Symbol}, ica::Union{InfoIca,Nothing} = nothing)
# 
#     function butterfly_plot(active::Bool)
#         if active
#             # Set all offsets to 0 for butterfly plot
#             offset[] = LinRange(0, 0, length(channel_labels))
#             # Hide labels
#             for (_, label) in channel_data_labels
#                 label.visible = false
#             end
#         else
#             # Reset to original spacing
#             offset[] = LinRange((yrange.val[end] * 0.9), yrange.val[1] * 0.9, length(channel_labels))
#             # Show labels
#             for (_, label) in channel_data_labels
#                 label.visible = true
#             end
#         end
#     end
# 
# 
#     function apply_lp_filter(active)
#         clear_axes(ax, [channel_data_original, channel_data_labels])
#         if active
#             data = filter_data(data, channel_labels, "lp", slider_lp_filter.value.val, 6, dat.sample_rate)
#         else
#             data = copy(dat.data)
#         end
#         draw(plot_labels = true)
#     end
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
#     function step_epoch_forward()
#         clear_axes(ax, [channel_data_original, channel_data_labels])
#         trial[] = min(length(dat.data), trial.val[1] + 1)
#         ax.title = "Epoch $(trial.val)/$(length(dat.data))"
#         menu_trial.i_selected[] = trial
#         update_extreme_spans!()
#         update_markers!(markers)
#         draw()
#     end
# 
#     function step_epoch_backward()
#         clear_axes(ax, [channel_data_original, channel_data_labels])
#         trial[] = max(1, trial.val[1] - 1)
#         ax.title = "Epoch $(trial.val)/$(length(dat.data))"
#         menu_trial.i_selected[] = trial
#         update_extreme_spans!()
#         update_markers!(markers)
#         draw()
#     end
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
#     # menu for electrode/channel selection
#     labels_menu = hcat(
#         Menu(
#             fig,
#             options = vcat(["All", "Left", "Right", "Central"], channel_labels_original),
#             default = "All",
#             direction = :down,
#             fontsize = 18,
#             width = 200,  # Match width
#         ),
#         Label(fig, "Labels", fontsize = 22, halign = :left),
#     )
#     on(labels_menu[1].selection) do s
#         channel_labels = [s]
#         if s == "All"
#             channel_labels = channel_labels_original
#         elseif s == "Left"
#             channel_labels = channel_labels_original[findall(occursin.(r"\d*[13579]$", channel_labels_original))]
#         elseif s == "Right"
#             channel_labels = channel_labels_original[findall(occursin.(r"\d*[24680]$", channel_labels_original))]
#         elseif s == "Central"
#             channel_labels = channel_labels_original[findall(occursin.(r"z$", channel_labels_original))]
#         end
# 
#         nchannels = length(channel_labels)
# 
#         clear_axes(ax, [channel_data_original, channel_data_labels])
# 
#         data = copy(dat.data)
#         if nchannels > 1
#             offset = LinRange(yrange.val[end] * 0.9, yrange.val[1] * 0.9, nchannels + 2)[2:end-1]
#         else # just centre
#             offset = zeros(nchannels)
#         end
#         draw(plot_labels = true)
# 
# 
#     end
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
#     # keyboard events
#     on(events(fig).keyboardbutton) do event
#         if event.action in (Keyboard.press,)
#             event.key == Keyboard.left && step_epoch_backward()
#             event.key == Keyboard.right && step_epoch_forward()
#             event.key == Keyboard.down && yless!(ax, yrange)
#             event.key == Keyboard.up && ymore!(ax, yrange)
#         end
#         # TODO: what is best here?
#         # return Consume()
#         # return Consume(false)
#     end
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
# 
#     function draw(; plot_labels = true)
#         # Create one plot per channel and update its data
#         for (idx, col) in enumerate(channel_labels)
#             # Initial empty plot with label
#             line = series!(
#                 ax,
#                 Point2f[],  # Will be updated by the lift
#                 color = :darkgrey,
#                 linewidth = 2,
#                 label = col,  # Add label to the line
#                 markersize = 0,  # Hide markers
#             )
#             channel_data_original[col] = line
# 
#             # Update plot data when range changes
#             lift(xrange, state.crit_val) do current_range, cv
#                 points = [Point2f(data[i, :time], data[i, col] + offset[idx]) for i in current_range]
#                 line[1][] = points
#             end
#         end
# 
#         # Update xlims separately
#         lift(xrange) do current_range
#             xlims!(ax, data[current_range[1], :time], data[current_range[end], :time])
#         end
#     end
# 
#     # plot theme adjustments
#     fontsize_theme = Theme(fontsize = 24)
#     update_theme!(fontsize_theme)
# 
#     hideydecorations!(ax, label = true)
#     draw(plot_labels = true)
#     display(fig)
#     # DataInspector(fig)
# 
# end
# 
# # Convenience methods
# plot_databrowser(dat::EpochData) = plot_databrowser(dat, dat.layout.label)
# plot_databrowser(dat::EpochData, ica::InfoIca) = plot_databrowser(dat, dat.layout.label, ica)
# plot_databrowser(dat::EpochData, channel_label::Symbol) = plot_databrowser(dat, [channel_label])
# plot_databrowser(dat::EpochData, channel_label::Symbol, ica::InfoIca) = plot_databrowser(dat, [channel_label], ica)

