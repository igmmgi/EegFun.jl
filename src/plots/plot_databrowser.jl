
# Base type for data states
abstract type AbstractDataState end

# Data Browser: Continuous Data
mutable struct Marker
    data::Any
    line::Any
    text::Any
    name::Symbol
    visible::Bool
end

mutable struct IcaState
    removed_activations::Union{Nothing,Any}
    components_to_remove::Union{Nothing,Vector{Int}}
    components_removed::Union{Nothing,Vector{Int}}
    IcaState() = new(nothing, nothing, nothing)
end

mutable struct ExtraChannelVis
    visualization::Union{Nothing,Makie.Lines,Makie.PolyElement,Any}
    label::Union{Nothing,Makie.Text}
end

mutable struct SelectionState
    active::Observable{Bool}
    bounds::Observable{Tuple{Float64,Float64}}
    visible::Observable{Bool}
    rectangle::Makie.Poly
    function SelectionState(ax)
        initial_points = [Point2f(0.0, 0.0)]
        poly_element = poly!(ax, initial_points, color = (:blue, 0.3), visible = false)
        new(Observable(false), Observable((0.0, 0.0)), Observable(false), poly_element)
    end
end

mutable struct FilterState
    active::Observable{NamedTuple{(:hp, :lp),Tuple{Bool,Bool}}}
    hp_freq::Observable{Float64}
    lp_freq::Observable{Float64}
    FilterState() = new(Observable((hp = false, lp = false)), Observable(0.1), Observable(40.0))
end

struct ToggleConfig
    label::String
    action::Function
end

mutable struct ViewState
    xrange::Observable{UnitRange{Int64}}
    yrange::Observable{UnitRange{Int64}}
    offset::Vector{Float64}
    crit_val::Observable{Float64}
    butterfly::Observable{Bool}
    function ViewState(n_channels::Int)
        # TODO: could this be done better?
        offset = n_channels > 1 ? LinRange(1500 * 0.9, -1500 * 0.9, n_channels + 2)[2:(end-1)] : zeros(n_channels)
        new(Observable(1:5000), Observable(-1500:1500), offset, Observable(0.0), Observable(false))
    end
end

mutable struct ChannelState
    labels::Vector{Symbol}
    visible::Vector{Bool}
    selected::Vector{Bool}
    data_labels::Dict{Symbol,Makie.Text}
    data_lines::Dict{Symbol,Union{Makie.Lines,Makie.PolyElement,Any}}
    function ChannelState(channel_labels::Vector{Symbol})
        new(
            channel_labels,
            fill(true, length(channel_labels)),
            fill(false, length(channel_labels)),
            Dict{Symbol,Makie.Text}(),
            Dict{Symbol,Union{Makie.Lines,Makie.PolyElement,Any}}(),
        )
    end
end

# Concrete data state implementations
mutable struct ContinuousDataState <: AbstractDataState
    current::Observable{EegData}
    original::EegData
    filter_state::FilterState
    function ContinuousDataState(data::EegData)
        new(Observable(copy(data)), data, FilterState())
    end
end

mutable struct EpochedDataState <: AbstractDataState
    current::Observable{EegData}
    original::EegData
    filter_state::FilterState
    current_epoch::Observable{Int}
    function EpochedDataState(data::EegData)
        new(Observable(copy(data)), data, FilterState(), Observable(1))
    end
end

mutable struct ExtraChannelInfo
    channel::Union{Nothing,Symbol}
    visible::Bool
    data_lines::Dict{Symbol,Union{Makie.Lines,Makie.PolyElement,Any}}
    data_labels::Dict{Symbol,Makie.Text}
    ExtraChannelInfo() =
        new(nothing, false, Dict{Symbol,Union{Makie.Lines,Makie.PolyElement,Any}}(), Dict{Symbol,Makie.Text}())
end


mutable struct DataBrowserState{T<:AbstractDataState}
    view::ViewState
    channels::ChannelState
    data::T
    selection::SelectionState
    markers::Vector{Marker}
    ica_state::Union{Nothing,IcaState}
    extra_channel::ExtraChannelInfo
    reference_state::Symbol

    # Constructor
    function DataBrowserState{T}(;
        view::ViewState,
        channels::ChannelState,
        data::T,
        selection::SelectionState,
        ica_state::Union{Nothing,IcaState} = nothing,
        extra_channel::ExtraChannelInfo = ExtraChannelInfo(),
    ) where {T<:AbstractDataState}
        return new{T}(
            view,
            channels,
            data,
            selection,
            Vector{Marker}(),
            ica_state,
            extra_channel,
            data.original.analysis_info.reference,
        )
    end
end

# Type aliases for code readability and simplified type annotations
const ContinuousDataBrowserState = DataBrowserState{ContinuousDataState}
const EpochedDataBrowserState = DataBrowserState{EpochedDataState}

# Single function using multiple dispatch for the data state creation
function create_browser_state(dat::T, channel_labels, ax, ica) where {T<:EegData}
    state_type = data_state_type(T)
    return DataBrowserState{state_type}(
        view = ViewState(length(channel_labels)),
        channels = ChannelState(channel_labels),
        data = state_type(dat),  # This directly calls the constructor!
        selection = SelectionState(ax),
        ica_state = !isnothing(ica) ? IcaState() : nothing,
    )
end

# Type mapping
data_state_type(::Type{ContinuousData}) = ContinuousDataState
data_state_type(::Type{EpochData}) = EpochedDataState

# Helper functions for common data access/resetting/updating
get_current_data(state::ContinuousDataState) = state.current[].data
get_current_data(state::EpochedDataState) = state.current[].data[state.current_epoch[]]
get_time_bounds(dat::ContinuousDataState) = (dat.current[].data.time[1], dat.current[].data.time[end])
get_time_bounds(dat::EpochedDataState) =
    (dat.current[].data[dat.current_epoch[]].time[1], dat.current[].data[dat.current_epoch[]].time[end])
has_column(state::ContinuousDataState, col::String) = col in names(state.current[].data)
has_column(state::EpochedDataState, col::String) = col in names(state.current[].data[state.current_epoch[]])

notify_data_update(state::AbstractDataState) = notify(state.current)

# This single function works for BOTH types
function reset_to_original!(state::AbstractDataState)
    state.current[] = copy(state.original)
end

############
# UI
############
function setup_ui_base(fig, ax, state, dat, ica = nothing)
    # Controls
    deregister_interaction!(ax, :rectanglezoom)

    # Set axes
    set_axes!(ax, state)

    # Mouse and keyboard events
    handle_mouse_events!(ax, state)
    handle_keyboard_events!(fig, ax, state)

    # Create toggles
    toggles = create_toggles(fig, ax, state)

    # Vertical line markers
    state.markers = init_markers(ax, state)

    # Create standard menus
    labels_menu = create_labels_menu(fig, ax, state)
    reference_menu = create_reference_menu(fig, state, dat)

    # Create optional menus
    ica_menu = nothing
    if !isnothing(ica) && !isnothing(state.ica_state)
        ica_menu = create_ica_menu(fig, ax, state, ica)
    end

    # Extra channel menu
    extra_menu = create_extra_channel_menu(fig, ax, state, dat)

    return (toggles, labels_menu, reference_menu, ica_menu, extra_menu)
end

# Unified setup_ui method using multiple dispatch for the epoch menu
function setup_ui(fig, ax, state::DataBrowserState{<:AbstractDataState}, dat, ica = nothing)
    # Get common UI elements
    toggles, labels_menu, reference_menu, ica_menu, extra_menu = setup_ui_base(fig, ax, state, dat, ica)

    # Get type-specific epoch menu (or nothing)
    epoch_menu = get_epoch_menu(fig, ax, state)

    # Build the grid components
    build_grid_components!(fig, dat, state, toggles, labels_menu, reference_menu, ica_menu, extra_menu, epoch_menu)

    # Apply theme
    update_theme!(Theme(fontsize = 18))
    hideydecorations!(ax, label = true)

    return state
end

function create_toggles(fig, ax, state)
    configs = [
        ToggleConfig("Butterfly Plot", (active) -> butterfly_plot!(ax, state)),
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

function create_menu(fig, options, default, label; kwargs...)
    menu = Menu(
        fig,
        options = options,
        default = default,
        direction = :down,
        fontsize = 18,
        width = Auto(),
        tellwidth = false,
        kwargs...,
    )
    return hcat(menu, Label(fig, label, fontsize = 22, halign = :left, tellwidth = false))
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
        state.reference_state = s
        s == :none && return
        rereference!(state.data, s)
        notify_data_update(state.data)
    end

    return menu
end

# Update create_ica_menu to use the new struct
function create_ica_menu(fig, ax, state, ica)

    options = vcat(["None"], ica.ica_label)
    menu = create_menu(fig, options, "None", "ICA Components")

    on(menu[1].selection) do s
        clear_axes!(ax, [state.channels.data_lines, state.channels.data_labels])
        component_to_remove_int = extract_int(String(s))
        state.ica_state.components_to_remove = isnothing(component_to_remove_int) ? nothing : [component_to_remove_int]

        if !isnothing(state.ica_state.components_to_remove)
            # Restore previous state if necessary before applying new removal
            if !isnothing(state.ica_state.components_removed)
                apply_ica_restore!(
                    state.data,
                    ica,
                    state.ica_state.components_removed,
                    state.ica_state.removed_activations,
                )
            end
            # Apply removal using the helper function
            state.ica_state.removed_activations =
                apply_ica_removal!(state.data, ica, state.ica_state.components_to_remove)
            state.ica_state.components_removed = state.ica_state.components_to_remove
        else # Selected "None"
            # Restore previous state if components were removed
            if !isnothing(state.ica_state.components_removed)
                apply_ica_restore!(
                    state.data,
                    ica,
                    state.ica_state.components_removed,
                    state.ica_state.removed_activations,
                )
                # Reset ICA state tracking
                state.ica_state.removed_activations = nothing
                state.ica_state.components_removed = nothing
            end
        end

        notify_data_update(state.data)
        draw(ax, state)
    end

    return menu

end

function create_epoch_menu(fig, ax, state)
    menu = create_menu(fig, 1:n_epochs(state.data.original), state.data.current_epoch[], "Epoch")

    on(menu[1].selection) do s
        clear_axes!(ax, [state.channels.data_lines, state.channels.data_labels])
        state.data.current_epoch[] = s
        ax.title = "Epoch $(s)/$(n_epochs(state.data.original))"
        update_markers!(ax, state)
        draw(ax, state)
        draw_extra_channel!(ax, state)
    end

    return menu
end

function show_additional_menu(state)

    # Create the menu figure
    menu_fig = Figure()
    plot_types = ["Topoplot (multiquadratic)", "Topoplot (spherical_spline)", "Spectrum"]

    menu_buttons = [Button(menu_fig[idx, 1], label = plot_type) for (idx, plot_type) in enumerate(plot_types)]

    for btn in menu_buttons
        on(btn.clicks) do n
            selected_data = subset_selected_data(state)
            if btn.label[] == "Topoplot (multiquadratic)"
                plot_topography(selected_data, method = :multiquadratic)
            elseif btn.label[] == "Topoplot (spherical_spline)"
                plot_topography(selected_data, method = :spherical_spline)
            elseif btn.label[] == "Spectrum"
                plot_channel_spectrum(selected_data)
            end
        end
    end

    new_screen = getfield(Main, :GLMakie).Screen()
    display(new_screen, menu_fig)
end

# Create common sliders for both continuous and epoched data
function create_common_sliders(fig, state, dat)
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

    return sliders
end

# Main create_sliders function for ContinuousDataBrowserState
function create_sliders(fig, state::ContinuousDataBrowserState, dat)
    sliders = create_common_sliders(fig, state, dat)

    # Add navigation sliders specific to continuous data
    slider_range = Slider(fig[3, 1], range = 100:50:30000, startvalue = state.view.xrange[][end], snap = true)
    slider_x = Slider(fig[2, 1], range = 1:50:nrow(state.data.current[].data), startvalue = 1, snap = true)

    on(slider_range.value) do x
        new_range = slider_x.value.val:min(nrow(state.data.current[].data), x+slider_x.value.val)
        if length(new_range) > 1
            state.view.xrange[] = new_range
        end
    end

    on(slider_x.value) do x
        new_range = x:min(nrow(state.data.current[].data), (x+slider_range.value.val)-1)
        if length(new_range) > 1
            state.view.xrange[] = new_range
        end
    end

    return sliders
end

# Create_sliders function for EpochedDataBrowserState (only common sliders)
function create_sliders(fig, state::EpochedDataBrowserState, dat)
    return create_common_sliders(fig, state, dat)
end

function create_extra_channel_menu(fig, ax, state, dat)
    menu =
        Menu(fig, options = [:none; extra_labels(dat)], default = "none", direction = :down, fontsize = 18, width = 200)

    on(menu.selection) do s
        state.extra_channel.channel = s == :none ? nothing : s
        state.extra_channel.visible = s != :none
        draw_extra_channel!(ax, state)
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
    epoch_menu = nothing,
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
    !isnothing(epoch_menu) && push!(grid_components, epoch_menu)

    # Use a Grid with auto-sizing for better responsiveness
    control_panel = grid!(vcat(grid_components...), tellheight = false)

    # Make control panel responsive with automatic sizing
    fig[1, 2] = control_panel

    # Use relative sizing for the control panel column
    colsize!(fig.layout, 2, Relative(0.25))
    rowsize!(fig.layout, 1, Relative(1.0))

end

# Get epoch menu based on state type
get_epoch_menu(fig, ax, state::ContinuousDataBrowserState) = nothing
get_epoch_menu(fig, ax, state::EpochedDataBrowserState) = create_epoch_menu(fig, ax, state)

############
# Navigation
############
const KEYBOARD_ACTIONS =
    Dict(Keyboard.left => :left, Keyboard.right => :right, Keyboard.up => :up, Keyboard.down => :down)

function handle_navigation!(ax, state::DataBrowserState{<:AbstractDataState}, action::Symbol)
    if action == :up
        ymore!(ax, state)
    elseif action == :down
        yless!(ax, state)
    elseif action == :left
        _handle_left_navigation(ax, state, state.data)
    elseif action == :right
        _handle_right_navigation(ax, state, state.data)
    end
end

# Type-specific left/right navigation
function _handle_left_navigation(ax, state, data::ContinuousDataState)
    xback!(ax, state)
end

function _handle_left_navigation(ax, state, data::EpochedDataState)
    step_epoch_backward(ax, state)
end

function _handle_right_navigation(ax, state, data::ContinuousDataState)
    xforward!(ax, state)
end

function _handle_right_navigation(ax, state, data::EpochedDataState)
    step_epoch_forward(ax, state)
end

function xback!(ax, state::ContinuousDataBrowserState)
    state.view.xrange.val[1] - 200 < 1 && return
    state.view.xrange[] = state.view.xrange.val .- 200
    xlims!(
        ax,
        state.data.current[].data.time[state.view.xrange.val[1]],
        state.data.current[].data.time[state.view.xrange.val[end]],
    )
end

function xforward!(ax, state::ContinuousDataBrowserState)
    state.view.xrange.val[1] + 200 > nrow(state.data.current[].data) && return
    state.view.xrange[] = state.view.xrange.val .+ 200
    xlims!(
        ax,
        state.data.current[].data.time[state.view.xrange.val[1]],
        state.data.current[].data.time[state.view.xrange.val[end]],
    )
end

function step_epoch_backward(ax, state::EpochedDataBrowserState)
    clear_axes!(ax, [state.channels.data_lines, state.channels.data_labels])
    state.data.current_epoch[] = max(1, state.data.current_epoch[] - 1)
    ax.title = "Epoch $(state.data.current_epoch[])/$(n_epochs(state.data.original))"
    update_markers!(ax, state)
    draw(ax, state)
    draw_extra_channel!(ax, state)
end

function step_epoch_forward(ax, state::EpochedDataBrowserState)
    clear_axes!(ax, [state.channels.data_lines, state.channels.data_labels])
    state.data.current_epoch[] = min(n_epochs(state.data.original), state.data.current_epoch[] + 1)
    ax.title = "Epoch $(state.data.current_epoch[])/$(n_epochs(state.data.original))"
    update_markers!(ax, state)
    draw(ax, state)
    draw_extra_channel!(ax, state)
end

function yless!(ax, state)
    (state.view.yrange.val[1] + 100 >= 0 || state.view.yrange.val[end] - 100 <= 0) && return
    state.view.yrange[] = (state.view.yrange.val[1]+100):(state.view.yrange.val[end]-100)
    ylims!(ax, state.view.yrange.val[1], state.view.yrange.val[end])
end

function ymore!(ax, state)
    state.view.yrange[] = (state.view.yrange.val[1]-100):(state.view.yrange.val[end]+100)
    ylims!(ax, state.view.yrange.val[1], state.view.yrange.val[end])
end

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
    # Make sure the polygon is visible when selection is visible
    state.selection.rectangle.visible[] = true
end

function handle_mouse_events!(ax, state)
    # Track if Shift and Ctrl are currently pressed
    shift_pressed = Ref(false)
    ctrl_pressed = Ref(false)

    # Listen for keyboard events to track Shift and Ctrl state
    on(events(ax).keyboardbutton) do key_event
        if key_event.key == Keyboard.left_shift
            shift_pressed[] = key_event.action == Keyboard.press
        elseif key_event.key == Keyboard.left_control
            ctrl_pressed[] = key_event.action == Keyboard.press
        end
    end

    on(events(ax).mousebutton) do event
        pos = events(ax).mouseposition[]
        if !is_mouse_in_axis(ax, pos)
            return
        end

        mouse_x = mouseposition(ax)[1]

        if event.button == Mouse.left
            if event.action == Mouse.press
                if ctrl_pressed[]
                    # Ctrl+Left press: Check for channel selection (immediate response)
                    mouse_y = mouseposition(ax)[2]
                    clicked_channel = find_closest_channel(ax, state, mouse_x, mouse_y)
                    if !isnothing(clicked_channel)
                        toggle_channel_visibility!(ax, state, clicked_channel)
                        return
                    end
                elseif shift_pressed[]
                    # Shift+Left press: Start time selection
                    handle_left_click!(ax, state, event, mouse_x)
                end
            elseif event.action == Mouse.release
                if shift_pressed[]
                    # Shift+Left release: Finish time selection
                    handle_left_click!(ax, state, event, mouse_x)
                end
            end
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
        show_additional_menu(state)
    end
end



# Helper function to find the closest channel to a click
function find_closest_channel(ax, state, mouse_x, mouse_y)
    current_data = get_current_data(state.data)

    min_distance = Inf
    closest_channel = nothing

    for (idx, visible) in enumerate(state.channels.visible)
        if visible
            col = state.channels.labels[idx]

            # Find the data point closest to the clicked x position
            time_data = current_data.time
            x_idx = argmin(abs.(time_data .- mouse_x))

            # Calculate the y position of this channel at the clicked x position
            channel_y = current_data[x_idx, col] + state.view.offset[idx]

            # Calculate distance to the clicked point
            distance = abs(mouse_y - channel_y)

            # If this is the closest channel so far and within reasonable distance
            if distance < min_distance && distance < 50  # 50 pixel tolerance
                min_distance = distance
                closest_channel = idx
            end
        end
    end

    return closest_channel
end

# Helper function to toggle channel selection
function toggle_channel_visibility!(ax, state, channel_idx)
    # Toggle the selection of the clicked channel
    state.channels.selected[channel_idx] = !state.channels.selected[channel_idx]

    # Immediate redraw for responsive feedback
    clear_axes!(ax, [state.channels.data_lines, state.channels.data_labels])
    draw(ax, state)
end

function handle_keyboard_events!(fig, ax, state)
    on(events(fig).keyboardbutton) do event
        if event.action in (Keyboard.press, Keyboard.repeat) && haskey(KEYBOARD_ACTIONS, event.key)
            action = KEYBOARD_ACTIONS[event.key]
            if state.selection.visible[]
                handle_selection_movement!(ax, state, action)
            else
                handle_navigation!(ax, state, action)
            end
        end
    end
end

function handle_selection_movement!(ax, state, action::Symbol)
    if action in (:left, :right)
        _handle_selection_movement_impl(ax, state, action)
    elseif action in (:up, :down)
        handle_navigation!(ax, state, action)
    end
end

function _handle_selection_movement_impl(ax, state::ContinuousDataBrowserState, action::Symbol)
    width = state.selection.bounds[][2] - state.selection.bounds[][1]
    time_start, time_end = get_time_bounds(state.data)
    if action == :left
        new_start = max(time_start, state.selection.bounds[][1] - width / 5)
    else  # :right
        new_start = min(time_end - width, state.selection.bounds[][1] + width / 5)
    end
    state.selection.bounds[] = (new_start, new_start + width)
    update_x_region_selection!(ax, state, state.selection.bounds[][1], state.selection.bounds[][2])
end

function _handle_selection_movement_impl(ax, state::EpochedDataBrowserState, action::Symbol)
    width = state.selection.bounds[][2] - state.selection.bounds[][1]
    current_epoch = state.data.current_epoch[]
    current_data = state.data.current[].data[current_epoch]
    time_start, time_end = current_data.time[1], current_data.time[end]

    if action == :left
        new_start = max(time_start, state.selection.bounds[][1] - width / 5)
    else  # :right
        new_start = min(time_end - width, state.selection.bounds[][1] + width / 5)
    end
    state.selection.bounds[] = (new_start, new_start + width)
    update_x_region_selection!(ax, state, state.selection.bounds[][1], state.selection.bounds[][2])
end

function update_x_region_selection!(ax, state, x1, x2)
    ylims = ax.limits[][2]
    state.selection.rectangle[1] = Point2f[
        Point2f(Float64(x1), Float64(ylims[1])),
        Point2f(Float64(x2), Float64(ylims[1])),
        Point2f(Float64(x2), Float64(ylims[2])),
        Point2f(Float64(x1), Float64(ylims[2])),
    ]
    # Ensure the polygon is visible when it has valid points
    state.selection.rectangle.visible[] = true
end

function clear_x_region_selection!(state)
    # Set to a single point instead of empty vector to avoid CairoMakie issues
    state.selection.rectangle[1] = [Point2f(0.0, 0.0)]
    state.selection.bounds[] = (0.0, 0.0)
    state.selection.visible[] = false
    # Make sure the polygon is hidden when selection is cleared
    state.selection.rectangle.visible[] = false
end

function subset_selected_data(state::ContinuousDataBrowserState)
    x_min, x_max = minmax(state.selection.bounds[]...)
    selected_channels = state.channels.labels[state.channels.visible]
    return subset(
        state.data.current[],
        sample_selection = x -> (x.time .>= x_min) .& (x.time .<= x_max),
        channel_selection = channels(selected_channels),
    )
end


function get_x_region_data(state::EpochedDataBrowserState)
    x_min, x_max = minmax(state.selection.bounds[]...)
    current_epoch = state.data.current_epoch[]
    current_data = state.data.current[].data[current_epoch]
    time_mask = (x_min .<= current_data.time .<= x_max)
    selected_data = current_data[time_mask, :]
    @info "Selected data: $(round(x_min, digits = 2)) to $(round(x_min, digits = 2)) S, size $(size(selected_data))"
    return selected_data
end



############
# Filtering
############
function apply_filter!(state::DataBrowserState{T}, filter_type, freq) where {T<:AbstractDataState}
    filter_data!(
        state.data.current[],
        String(filter_type),
        "iir",
        freq,
        order = filter_type == :hp ? 1 : 3,
        channel_selection = (channels) -> [ch in state.channels.labels for ch in channels],
    )
end

function apply_filters!(state)
    # Reset to original if no filters active
    if !state.data.filter_state.active[].hp && !state.data.filter_state.active[].lp
        reset_to_original!(state.data)
        rereference!(state.data, state.reference_state)
        notify_data_update(state.data)
        return
    end

    # Start with fresh data if changing filter configuration
    if state.data.filter_state.active[].hp != state.data.filter_state.active[].lp
        reset_to_original!(state.data)
        rereference!(state.data, state.reference_state)
    end

    # Apply active filters
    for (filter_type, freq) in zip([:hp, :lp], [state.data.filter_state.hp_freq[], state.data.filter_state.lp_freq[]])
        if state.data.filter_state.active[][filter_type]
            apply_filter!(state, filter_type, freq)
        end
    end

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

########################
# Reference
########################
function rereference!(state::AbstractDataState, ref)
    rereference!(state.current[], ref, channels())
end

########################
# ICA Component Application Helpers
########################

# Apply ICA component removal based on state type
function apply_ica_removal!(state::ContinuousDataState, ica::InfoIca, components_to_remove::Vector{Int})
    # Create new EegData with ICA components removed
    ica_data = copy(state.current[])
    ica_data.data, activations = remove_ica_components(state.current[].data, ica, components_to_remove)
    state.current[] = ica_data
    return activations
end

function apply_ica_removal!(state::EpochedDataState, ica::InfoIca, components_to_remove::Vector{Int})
    # Create new EegData with ICA components removed
    ica_data = copy(state.current[])
    first_epoch_activations = nothing
    for (i, epoch_df) in enumerate(ica_data.data)
        ica_data.data[i], activations = remove_ica_components(epoch_df, ica, components_to_remove)
        if i == 1
            first_epoch_activations = activations
        end
    end
    state.current[] = ica_data
    return first_epoch_activations
end

# Apply ICA component restoration based on state type
function apply_ica_restore!(
    state::ContinuousDataState,
    ica::InfoIca,
    components_removed::Vector{Int},
    removed_activations,
)
    if isnothing(removed_activations)
        @warn "Cannot restore ICA components: No previous activations stored."
        return
    end
    # Create new EegData with restored data
    restored_data = copy(state.current[])
    restored_data.data = restore_original_data(state.current[].data, ica, components_removed, removed_activations)
    state.current[] = restored_data
end

function apply_ica_restore!(state::EpochedDataState, ica::InfoIca, components_removed::Vector{Int}, removed_activations)
    if isnothing(removed_activations)
        @warn "Cannot restore ICA components: No previous activations stored."
        return
    end
    # Create new EegData with restored data
    restored_data = copy(state.current[])
    for (i, epoch_df) in enumerate(restored_data.data)
        restored_data.data[i] = restore_original_data(epoch_df, ica, components_removed, removed_activations)
    end
    state.current[] = restored_data
end


########################
# Drawing
########################
function add_marker!(markers, ax, data, col; label = nothing, trial = nothing, visible = false)
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
            vlines!(marker_data.time, color = :grey, linewidth = 1, visible = visible),
            text!(
                label,
                position = [(x, ax.yaxis.attributes.limits[][2] * 0.98) for x in marker_data.time],
                space = :data,
                align = (:center, :center),
                fontsize = 22,
                visible = visible,
            ),
            col,
            visible,
        ),
    )
end

function update_channel_offsets!(state)
    nchannels = count(state.channels.visible)
    if nchannels > 1 && !state.view.butterfly[]
        state.view.offset[state.channels.visible] .=
            LinRange((state.view.yrange[][end] * 0.9), state.view.yrange[][1] * 0.9, nchannels + 2)[2:(end-1)]
    else
        state.view.offset[state.channels.visible] .= zeros(nchannels)
    end
end

function clear_axes!(ax, datas)
    [delete!(ax, value) for data in datas for (key, value) in data]
    [empty!(data) for data in datas]
end

# Generic set_axes! function that handles both types
function set_axes!(ax, state::DataBrowserState{<:AbstractDataState})
    # Set y limits for both types
    @lift ylims!(ax, $(state.view.yrange)[1], $(state.view.yrange)[end])

    # Use type-specific method for setting x limits
    set_x_limits!(ax, state, state.data)
end

# Type-specific x limit setting for continuous data
function set_x_limits!(ax, state, data::ContinuousDataState)
    @lift xlims!(
        ax,
        $(data.current).data.time[$(state.view.xrange)[1]],
        $(data.current).data.time[$(state.view.xrange)[end]],
    )
end

# Type-specific x limit setting for epoched data
function set_x_limits!(ax, state, data::EpochedDataState)
    @lift xlims!(ax, $(data.current).data[1].time[1], $(data.current).data[1].time[end])
end

# Common marker initialization
function init_markers(ax, state; marker_visible = Dict{Symbol,Bool}())
    markers = Marker[]
    data = get_current_data(state.data)
    if has_column(state.data, "triggers")
        add_marker!(markers, ax, data, :triggers, visible = get(marker_visible, :triggers, false))
    end
    if has_column(state.data, "is_vEOG")
        add_marker!(markers, ax, data, :is_vEOG, label = "v", visible = get(marker_visible, :is_vEOG, false))
    end
    if has_column(state.data, "is_hEOG")
        add_marker!(markers, ax, data, :is_hEOG, label = "h", visible = get(marker_visible, :is_hEOG, false))
    end

    return markers
end

function update_markers!(ax, state)
    marker_visible = Dict{Symbol,Bool}()
    for marker in state.markers
        marker_visible[marker.name] = marker.visible
        delete!(ax, marker.line)
        delete!(ax, marker.text)
    end
    empty!(state.markers)
    state.markers = init_markers(ax, state; marker_visible = marker_visible)
end

function butterfly_plot!(ax, state)
    state.view.butterfly[] = !state.view.butterfly[]
    clear_axes!(ax, [state.channels.data_lines, state.channels.data_labels])
    update_channel_offsets!(state)
    draw(ax, state)
end

# Single function with data access abstraction
function draw(ax, state::DataBrowserState{<:AbstractDataState})
    # Get data access functions based on type
    get_data, get_time, get_label_y = get_data_accessors(state.data)

    # Pre-compute shared observables
    visible_time_obs = @lift(get_time($(state.data.current), $(state.view.xrange)))
    time_start_obs = @lift(get_time($(state.data.current), $(state.view.xrange)[1:1])[1])

    @sync for (idx, visible) in enumerate(state.channels.visible)
        col = state.channels.labels[idx]
        if visible
            is_selected = state.channels.selected[idx]

            # Line properties
            line_color =
                is_selected ? :black :
                @lift(abs.(get_data($(state.data.current), $(state.view.xrange), $col)) .>= $(state.view.crit_val))
            line_colormap = is_selected ? [:black] : [:darkgrey, :darkgrey, :red]
            line_width = is_selected ? 4 : 2

            # Channel data
            channel_data_obs = @lift(get_data($(state.data.current), $(state.view.xrange), $col))
            channel_data_with_offset = @lift($(channel_data_obs) .+ state.view.offset[idx])

            # Update or create line
            update_or_create_line!(
                state.channels.data_lines,
                col,
                ax,
                visible_time_obs,
                channel_data_with_offset,
                line_color,
                line_colormap,
                line_width,
            )

            # Handle labels
            if !state.view.butterfly[]
                label_y_obs = @lift(get_label_y($(state.data.current), $col, state.view.offset[idx]))
                update_or_create_label!(state.channels.data_labels, col, ax, time_start_obs, label_y_obs, is_selected)
            else
                hide_channel_label!(state.channels.data_labels, col)
            end
        else
            hide_channel_objects!(state.channels, col)
        end
    end
end

# Data accessor functions
function get_data_accessors(data::ContinuousDataState)
    get_data = (current, range, col) -> current.data[range, col]
    get_time = (current, range) -> current.data.time[range]
    get_label_y = (current, col, offset) -> current.data[1, col] .+ offset
    return get_data, get_time, get_label_y
end

function get_data_accessors(data::EpochedDataState)
    get_data = (current, range, col) -> current.data[current.current_epoch[]][range, col]
    get_time = (current, range) -> current.data[current.current_epoch[]].time[range]
    get_label_y = (current, col, offset) -> current.data[current.current_epoch[]][!, col][1] .+ offset
    return get_data, get_time, get_label_y
end

# Helper functions for line/label management
function update_or_create_line!(data_lines, col, ax, x_obs, y_obs, color, colormap, linewidth)
    if haskey(data_lines, col)
        data_lines[col].x[] = x_obs[]
        data_lines[col].y[] = y_obs[]
        data_lines[col].color[] = color
        data_lines[col].colormap[] = colormap
        data_lines[col].linewidth[] = linewidth
        show!(data_lines[col])
    else
        data_lines[col] = lines!(ax, x_obs, y_obs, color = color, colormap = colormap, linewidth = linewidth)
    end
end

function update_or_create_label!(data_labels, col, ax, x_obs, y_obs, is_selected)
    if haskey(data_labels, col)
        data_labels[col].position[] = Point2f(x_obs[], y_obs[])
        data_labels[col].color[] = is_selected ? :red : :black
        show!(data_labels[col])
    else
        data_labels[col] = text!(
            ax,
            x_obs,
            y_obs,
            text = String(col),
            align = (:left, :center),
            fontsize = 18,
            color = is_selected ? :red : :black,
        )
    end
end

function hide_channel_label!(data_labels, col)
    if haskey(data_labels, col)
        hide!(data_labels[col])
    end
end

function hide_channel_objects!(channels, col)
    if haskey(channels.data_lines, col)
        hide!(channels.data_lines[col])
    end
    if haskey(channels.data_labels, col)
        hide!(channels.data_labels[col])
    end
end

# Single function with data access abstraction
function draw_extra_channel!(ax, state::DataBrowserState{<:AbstractDataState})
    clear_axes!(ax, [state.extra_channel.data_lines, state.extra_channel.data_labels])

    if state.extra_channel.visible && !isnothing(state.extra_channel.channel)
        current_offset = state.view.offset[end] + mean(diff(state.view.offset))
        channel = state.extra_channel.channel

        # Get data access functions based on type
        get_data, get_time, get_label_y = get_data_accessors(state.data)

        if eltype(get_data(state.data.current[], 1:1, channel)) == Bool
            # Boolean data - create highlights
            highlight_data = @views splitgroups(findall(get_data(state.data.current[], :, channel)))
            region_offset = all(iszero, highlight_data[2] .- highlight_data[1]) ? 5 : 0
            state.extra_channel.data_lines[channel] = vspan!(
                ax,
                get_time(state.data.current[], highlight_data[1]),
                get_time(state.data.current[], highlight_data[2] .+ region_offset),
                color = :Red,
                alpha = 0.5,
                visible = true,
            )
        else
            # Regular data - create line and label
            state.extra_channel.data_lines[channel] = lines!(
                ax,
                @lift(get_time($(state.data.current), :)),
                @lift(get_data($(state.data.current), :, $channel) .+ $current_offset),
                color = :black,
                linewidth = 2,
            )
            state.extra_channel.data_labels[channel] = text!(
                ax,
                @lift(get_time($(state.data.current), $(state.view.xrange)[1:1])[1]),
                @lift(get_label_y($(state.data.current), $channel, $current_offset)),
                text = String(channel),
                align = (:left, :center),
                fontsize = 18,
            )
        end
    end
end

# Helper functions for getting type-specific information
get_title(dat::EpochData) = "Epoch 1/$(n_epochs(dat))"
get_title(dat::ContinuousData) = ""

function plot_databrowser(dat::EegData, ica = nothing)

    @info "plot_databrowser: ..."

    # Check if CairoMakie is being used and warn about interactivity
    if string(Makie.current_backend()) == "CairoMakie"
        @minimal_warning "CairoMakie detected. For full interactivity in plot_databrowser, use GLMakie."
    end

    # Common fig/ax/state/ui setup
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel = "Time (S)", ylabel = "Amplitude (μV)", title = get_title(dat))
    state = create_browser_state(dat, dat.layout.data.label, ax, ica)
    setup_ui(fig, ax, state, dat, ica)

    # Render and return
    draw(ax, state)
    draw_extra_channel!(ax, state)
    display(fig)
    return fig, ax
end

function plot_vertical_lines!(ax, marker, active)
    marker.line.visible = active
    marker.text.visible = active
    marker.visible = active
    marker.text.position = [(x, ax.yaxis.attributes.limits[][2] * 0.98) for x in marker.data.time] # incase y changed
end



# Subset-based plot_databrowser methods
"""
    plot_databrowser(dat::Union{ContinuousData,EpochData}; 
                    channel_selection::Function = channels(), 
                    sample_selection::Function = samples(),
                    epoch_selection::Function = epochs(),
                    ica = nothing)

Plot a subset of EEG data using the databrowser interface.

# Arguments
- `dat`: EEG data object (ContinuousData or EpochData)
- `channel_selection`: Function that returns channel labels to include (default: all channels)
- `sample_selection`: Function that returns sample mask to include (default: all samples)
- `epoch_selection`: Function that returns epoch mask to include (default: all epochs, only for EpochData)
- `ica`: Optional ICA result object for component visualization

# Returns
- Figure and Axis objects for the databrowser plot

# Examples
```julia
# Plot subset by channels only
plot_databrowser(dat, channel_selection = channels([:Fp1, :Fp2]))

# Plot subset by samples only
plot_databrowser(dat, sample_selection = x -> x.sample .< 1000)

# Plot subset by epochs only (EpochData only)
plot_databrowser(dat, epoch_selection = epochs(1:10))

# Plot subset with all filters
plot_databrowser(dat, 
    channel_selection = channels([:Fp1, :Fp2]), 
    sample_selection = x -> x.sample .< 1000,
    epoch_selection = epochs([1, 3, 5]))

# Plot subset with ICA
plot_databrowser(dat, 
    channel_selection = channels([:Fp1, :Fp2]), 
    ica = ica_result)
```
"""
# TODO: cannot dispatch on kwargs so need another function name: is there a better way?
function plot_databrowser_subset(
    dat::ContinuousData;
    channel_selection::Function = channels(),
    sample_selection::Function = samples(),
    ica = nothing,
)
    # Create data subset and plot
    filtered_dat = subset(dat, channel_selection = channel_selection, sample_selection = sample_selection)
    return plot_databrowser(filtered_dat, ica)
end

function plot_databrowser_subset(
    dat::EpochData;
    channel_selection::Function = channels(),
    sample_selection::Function = samples(),
    epoch_selection::Function = epochs(),
    ica = nothing,
)
    # Create data subset and plot
    filtered_dat = subset(
        dat,
        channel_selection = channel_selection,
        sample_selection = sample_selection,
        epoch_selection = epoch_selection,
    )
    return plot_databrowser(filtered_dat, ica)
end
