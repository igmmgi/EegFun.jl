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
    # Add drag detection state
    drag_start_pos::Union{Nothing,Tuple{Float64,Float64}}
    drag_start_time::Union{Nothing,Float64}
    is_dragging::Bool
    drag_threshold::Float64  # Minimum distance to consider it a drag
    click_timeout::Float64   # Maximum time for a click (seconds)
    
    function SelectionState(ax)
        # Initialize with a single point to avoid empty vector issues with CairoMakie
        initial_points = [Point2f(0.0, 0.0)]
        poly_element = poly!(ax, initial_points, color = (:blue, 0.3), visible = false)
        new(Observable(false), Observable((0.0, 0.0)), Observable(false), poly_element, 
            nothing, nothing, false, 3.0, 0.15)  # 3 pixel threshold, 150ms timeout
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
    display_scale::Observable{Float64}  # New display scale factor
    function ViewState(n_channels::Int)
        offset = n_channels > 1 ? LinRange(1500 * 0.9, -1500 * 0.9, n_channels + 2)[2:end-1] : zeros(n_channels)
        new(Observable(1:5000), Observable(-1500:1500), offset, Observable(0.0), Observable(false), Observable(1.0))
    end
end

mutable struct ChannelState
    labels::Vector{Symbol}
    visible::Vector{Bool}
    data_labels::Dict{Symbol,Makie.Text}
    data_lines::Dict{Symbol,Union{Makie.Lines,Makie.PolyElement,Any}}
    # Cache for computed Observables to avoid recreation
    cached_observables::Dict{Symbol,Any}
    # Track highlighted channels
    highlighted::Observable{Symbol}
    function ChannelState(channel_labels::Vector{Symbol})
        new(
            channel_labels,
            fill(true, length(channel_labels)),
            Dict{Symbol,Makie.Text}(),
            Dict{Symbol,Union{Makie.Lines,Makie.PolyElement,Any}}(),
            Dict{Symbol,Any}(),  # Cache for computed Observables
            Observable{Symbol}(:none),  # Track highlighted channel
        )
    end
end

# Concrete data state implementations
mutable struct ContinuousDataState <: AbstractDataState
    current::Observable{DataFrame}
    original::EegData
    filter_state::FilterState
    function ContinuousDataState(data::EegData)
        new(Observable(copy(data.data)), data, FilterState())
    end
end

mutable struct EpochedDataState <: AbstractDataState
    current::Vector{Observable{DataFrame}}
    original::EegData
    filter_state::FilterState
    current_epoch::Observable{Int}
    function EpochedDataState(data::EegData)
        new([Observable(copy(df)) for df in data.data], data, FilterState(), Observable(1))
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

# Single unified browser state with a type parameter
mutable struct DataBrowserState{T<:AbstractDataState}
    view::ViewState
    channels::ChannelState
    data::T
    selection::SelectionState
    markers::Vector{Marker}
    ica_state::Union{Nothing,IcaState}
    extra_channel::ExtraChannelInfo
    reference_state::Symbol

    # Unified constructor
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

# Create the appropriate state type
function create_browser_state(dat::ContinuousData, channel_labels, ax, ica)
    return DataBrowserState{ContinuousDataState}(
        view = ViewState(length(channel_labels)),
        channels = ChannelState(channel_labels),
        data = ContinuousDataState(dat),
        selection = SelectionState(ax),
        ica_state = !isnothing(ica) ? IcaState() : nothing,
    )
end

function create_browser_state(dat::EpochData, channel_labels, ax, ica)
    return DataBrowserState{EpochedDataState}(
        view = ViewState(length(channel_labels)),
        channels = ChannelState(channel_labels),
        data = EpochedDataState(dat),
        selection = SelectionState(ax),
        ica_state = !isnothing(ica) ? IcaState() : nothing,
    )
end

# Helper functions for Observable caching
function get_or_create_cached_observable!(channel_state::ChannelState, key::Symbol, create_func::Function)
    """
    Get an Observable from cache or create it if it doesn't exist.
    
    # Arguments
    - `channel_state::ChannelState`: The channel state containing the cache
    - `key::Symbol`: The cache key for this Observable
    - `create_func::Function`: Function to create the Observable if not cached
    
    # Returns
    - The cached or newly created Observable
    """
    if haskey(channel_state.cached_observables, key)
        return channel_state.cached_observables[key]
    else
        obs = create_func()
        channel_state.cached_observables[key] = obs
        return obs
    end
end

function clear_observable_cache!(channel_state::ChannelState)
    """
    Clear the Observable cache when data changes significantly.
    """
    empty!(channel_state.cached_observables)
end

function clear_observable_cache!(state::DataBrowserState)
    """
    Clear the Observable cache for the entire browser state.
    """
    clear_observable_cache!(state.channels)
end

# Helper functions for common data access/resetting/updating
get_current_data(state::ContinuousDataState) = state.current[]
get_current_data(state::EpochedDataState) = state.current[state.current_epoch[]][]
get_time_bounds(dat::ContinuousDataState) = (dat.current[].time[1], dat.current[].time[end])
get_time_bounds(dat::EpochedDataState) =
    (dat.current[dat.current_epoch[]][].time[1], dat.current[dat.current_epoch[]][].time[end])
has_column(state::ContinuousDataState, col::String) = col in names(state.current[])
has_column(state::EpochedDataState, col::String) = col in names(state.current[state.current_epoch[]][])

notify_data_update(state::ContinuousDataState) = notify(state.current)
notify_data_update(state::EpochedDataState) = notify(state.current[state.current_epoch[]])

function reset_to_original!(state::ContinuousDataState)
    state.current[] = copy(state.original.data)
end

function reset_to_original!(state::EpochedDataState)
    for (i, df) in enumerate(state.original.data)
        state.current[i][] = copy(df)
    end
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
    handle_mouse_events!(fig, ax, state)
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
        # Clear cache when changing channel visibility
        clear_observable_cache!(state.channels)
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
        rereference!(state.data, dat, s)
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
        # Clear cache when applying ICA changes
        clear_observable_cache!(state.channels)
        component_to_remove_int = extract_int(String(s))
        state.ica_state.components_to_remove = isnothing(component_to_remove_int) ? nothing : [component_to_remove_int]

        if !isnothing(state.ica_state.components_to_remove)
            # Restore previous state if necessary before applying new removal
            if !isnothing(state.ica_state.components_removed)
                 apply_ica_restore!(
                    state.data, 
                    ica, 
                    state.ica_state.components_removed, 
                    state.ica_state.removed_activations
                 )
            end
            # Apply removal using the helper function
            state.ica_state.removed_activations = apply_ica_removal!(
                state.data, 
                ica, 
                state.ica_state.components_to_remove
            )
            state.ica_state.components_removed = state.ica_state.components_to_remove
        else # Selected "None"
             # Restore previous state if components were removed
             if !isnothing(state.ica_state.components_removed)
                 apply_ica_restore!(
                    state.data, 
                    ica, 
                    state.ica_state.components_removed, 
                    state.ica_state.removed_activations
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
        # Clear cache when changing epochs
        clear_observable_cache!(state.channels)
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
    plot_types = ["Topoplot", "Spectrum", "Plot3"]

    menu_buttons = [Button(menu_fig[idx, 1], label = plot_type) for (idx, plot_type) in enumerate(plot_types)]

    for btn in menu_buttons
        on(btn.clicks) do n
            selected_data = get_x_region_data(state)
            if btn.label[] == "Topoplot"
                plot_topoplot(selected_data, state.data.original.layout)
            elseif btn.label[] == "Spectrum"
                selected_channel = state.channels.labels[state.channels.visible]
                plot_channel_spectrum(selected_data, channel_selection = channels(selected_channel))
            elseif btn.label[] == "Plot3"
                println("Plot3: TODO")
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
    slider_x = Slider(fig[2, 1], range = 1:50:nrow(state.data.current[]), startvalue = 1, snap = true)
    
    on(slider_range.value) do x
        new_range = slider_x.value.val:min(nrow(state.data.current[]), x + slider_x.value.val)
        if length(new_range) > 1
            state.view.xrange[] = new_range
        end
    end

    on(slider_x.value) do x
        new_range = x:min(nrow(state.data.current[]), (x + slider_range.value.val) - 1)
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
    menu = Menu(
        fig,
        options = [:none; extra_channels(dat)],
        default = "none",
        direction = :down,
        fontsize = 18,
        width = 200,
    )

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
    current_range = state.view.xrange[]
    current_range[1] - 200 < 1 && return
    state.view.xrange[] = current_range .- 200
    xlims!(
        ax,
        state.data.current[].time[state.view.xrange[][1]],
        state.data.current[].time[state.view.xrange[][end]],
    )
end

function xforward!(ax, state::ContinuousDataBrowserState)
    current_range = state.view.xrange[]
    current_range[1] + 200 > nrow(state.data.current[]) && return
    state.view.xrange[] = current_range .+ 200
    xlims!(
        ax,
        state.data.current[].time[state.view.xrange[][1]],
        state.data.current[].time[state.view.xrange[][end]],
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
    # Zoom out by increasing display scale
    current_scale = state.view.display_scale[]
    state.view.display_scale[] = current_scale * 1.2
end

function ymore!(ax, state)
    # Zoom in by decreasing display scale
    current_scale = state.view.display_scale[]
    state.view.display_scale[] = current_scale * 0.8
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
    
    # Only make selection visible if it's actually a meaningful selection (not just a click)
    if abs(mouse_x - state.selection.bounds[][1]) > 0.01  # Small threshold to detect actual drag
        state.selection.visible[] = true
        state.selection.bounds[] = (state.selection.bounds[][1], mouse_x)
        update_x_region_selection!(ax, state, state.selection.bounds[][1], mouse_x)
        # Make sure the polygon is visible when selection is visible
        state.selection.rectangle.visible[] = true
    else
        # Just a click - clear the selection to ensure arrow keys work
        clear_x_region_selection!(state)
    end
end

function handle_mouse_events!(fig, ax, state)
    # Handle regular mouse events (selection, navigation)
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

    # Handle line clicks for highlighting (immediate response on press)
    on(events(fig).mousebutton) do event
        if event.button == Mouse.left && event.action == Mouse.press
            # Try to handle line click immediately
            handle_line_click!(fig, ax, state, event)
        end
    end

    # Update selection rectangle while dragging
    on(events(ax).mouseposition) do _
        if state.selection.active[]
            world_pos = mouseposition(ax)[1]
            update_x_region_selection!(ax, state, state.selection.bounds[][1], world_pos)
            
            # Check if we've started dragging
            if state.selection.drag_start_pos !== nothing && !state.selection.is_dragging
                current_pos = (world_pos, mouseposition(ax)[2])
                distance = sqrt((current_pos[1] - state.selection.drag_start_pos[1])^2 + 
                              (current_pos[2] - state.selection.drag_start_pos[2])^2)
                
                if distance > state.selection.drag_threshold
                    state.selection.is_dragging = true
                end
            end
        end
    end
end

# Add line highlighting functionality
function handle_line_click!(fig, ax, state, event)
    plt = pick(fig)
    if plt !== nothing && plt[1] isa Lines
        # Check if the clicked line corresponds to a channel
        for (channel, line) in state.channels.data_lines
            if plt[1] == line
                # Toggle highlighting: if already highlighted, clear it; otherwise highlight
                if state.channels.highlighted[] == channel
                    state.channels.highlighted[] = :none
                else
                    state.channels.highlighted[] = channel
                end
                # No need to redraw - the Observables will update automatically
                break
            end
        end
    else
        # Clicked on empty space - clear highlighting to ensure arrow keys work
        state.channels.highlighted[] = :none
    end
end

function handle_left_click!(ax, state, event, mouse_x)
    if event.action == Mouse.press
        # Check if we clicked on a line first
        plt = pick(ax.scene)
        clicked_on_line = false
        if plt !== nothing && plt[1] isa Lines
            for (channel, line) in state.channels.data_lines
                if plt[1] == line
                    clicked_on_line = true
                    break
                end
            end
        end
        
        # Only start selection if we didn't click on a line
        if !clicked_on_line
            # Record start position and time for drag detection
            state.selection.drag_start_pos = (mouse_x, mouseposition(ax)[2])
            state.selection.drag_start_time = time()
            state.selection.is_dragging = false
            
            if state.selection.visible[] && is_within_selection(state, mouse_x)
                clear_x_region_selection!(state)
            else
                start_selection!(ax, state, mouse_x)
            end
        end
    elseif event.action == Mouse.release && state.selection.active[]
        finish_selection!(ax, state, mouse_x)
        
        # Reset drag detection state
        state.selection.drag_start_pos = nothing
        state.selection.drag_start_time = nothing
        state.selection.is_dragging = false
    end
end

function handle_right_click!(ax, state, mouse_x)
    if state.selection.visible[] && is_within_selection(state, mouse_x)
        show_additional_menu(state)
    end
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
    current_data = state.data.current[current_epoch][]
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

function get_x_region_data(state::ContinuousDataBrowserState)
    x_min, x_max = minmax(state.selection.bounds[]...)
    time_mask = (x_min .<= state.data.current[].time .<= x_max)
    selected_data = state.data.current[][time_mask, :]
    println("Selected data: $(round(x_min, digits = 2)) to $(round(x_max, digits = 2)) S, size $(size(selected_data))")
    return selected_data
end

function get_x_region_data(state::EpochedDataBrowserState)
    x_min, x_max = minmax(state.selection.bounds[]...)
    current_epoch = state.data.current_epoch[]
    current_data = state.data.current[current_epoch][]
    time_mask = (x_min .<= current_data.time .<= x_max)
    selected_data = current_data[time_mask, :]
    println("Selected data: $(round(x_min, digits = 2)) to $(round(x_max, digits = 2)) S, size $(size(selected_data))")
    return selected_data
end



############
# Filtering
############
function apply_filter!(state::ContinuousDataBrowserState, filter_type, freq)
    state.data.current[] = filter_data(
        state.data.current[],
        state.channels.labels,
        String(filter_type),
        "iir",
        freq,
        sample_rate(state.data.original),
        order = filter_type == :hp ? 1 : 3,
    )
end

function apply_filter!(state::EpochedDataBrowserState, filter_type, freq)
    for (i, df) in enumerate(state.data.current)
        df[] = filter_data(
            df[],
            state.channels.labels,
            String(filter_type),
            "iir",
            freq,
            sample_rate(state.data.original),
            order = filter_type == :hp ? 1 : 3,
        )
        notify(df)  # Notify each Observable individually
    end
end

function apply_filters!(state)
    # Reset to original if no filters active
    if !state.data.filter_state.active[].hp && !state.data.filter_state.active[].lp
        reset_to_original!(state.data)
        rereference!(state.data, state.data.original, state.reference_state)
        notify_data_update(state.data)
        return
    end

    # Start with fresh data if changing filter configuration
    if state.data.filter_state.active[].hp != state.data.filter_state.active[].lp
        reset_to_original!(state.data)
        rereference!(state.data, state.data.original, state.reference_state)
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
function rereference!(state::ContinuousDataState, dat, ref)
    rereference!(state.current[], channels(dat), resolve_reference(dat, ref))
end

function rereference!(state::EpochedDataState, dat, ref)
    for df in state.current
        rereference!(df[], channels(dat), resolve_reference(dat, ref))
    end
end

########################
# ICA Component Application Helpers
########################

# Apply ICA component removal based on state type
function apply_ica_removal!(state::ContinuousDataState, ica::InfoIca, components_to_remove::Vector{Int})
    df_new, activations = remove_ica_components(state.current[], ica, components_to_remove)
    state.current[] = df_new # Update observable
    return activations
end

function apply_ica_removal!(state::EpochedDataState, ica::InfoIca, components_to_remove::Vector{Int})
    first_epoch_activations = nothing
    for (i, df_obs) in enumerate(state.current)
        df_new, activations = remove_ica_components(df_obs[], ica, components_to_remove)
        df_obs[] = df_new # Update observable
        if i == 1
            first_epoch_activations = activations
        end
    end
    return first_epoch_activations
end

# Apply ICA component restoration based on state type
function apply_ica_restore!(state::ContinuousDataState, ica::InfoIca, components_removed::Vector{Int}, removed_activations)
    if isnothing(removed_activations)
         @warn "Cannot restore ICA components: No previous activations stored."
         return
    end
    df_new = restore_original_data(state.current[], ica, components_removed, removed_activations)
    state.current[] = df_new # Update observable
end

function apply_ica_restore!(state::EpochedDataState, ica::InfoIca, components_removed::Vector{Int}, removed_activations)
    if isnothing(removed_activations)
         @warn "Cannot restore ICA components: No previous activations stored."
         return
    end
    for df_obs in state.current
        df_new = restore_original_data(df_obs[], ica, components_removed, removed_activations)
        df_obs[] = df_new # Update observable
    end
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
            LinRange((state.view.yrange[][end] * 0.9), state.view.yrange[][1] * 0.9, nchannels + 2)[2:end-1]
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
    # Create computed Observables for limits to reduce nesting
    y_range = state.view.yrange
    y_min_obs = @lift($(y_range)[1])
    y_max_obs = @lift($(y_range)[end])
    
    # Set y limits using computed Observables
    @lift ylims!(ax, $(y_min_obs), $(y_max_obs))
    
    # Use type-specific method for setting x limits
    set_x_limits!(ax, state, state.data)
end

# Type-specific x limit setting for continuous data
function set_x_limits!(ax, state, data::ContinuousDataState)
    # Create computed Observables for x limits to reduce nesting
    time_data = @lift($(data.current).time)
    x_start_obs = @lift($(time_data)[$(state.view.xrange)[1]])
    x_end_obs = @lift($(time_data)[$(state.view.xrange)[end]])
    
    @lift xlims!(ax, $(x_start_obs), $(x_end_obs))
end

# Type-specific x limit setting for epoched data
function set_x_limits!(ax, state, data::EpochedDataState)
    # Create computed Observable for x limits
    x_start_obs = @lift($(data.current[1]).time[1])
    x_end_obs = @lift($(data.current[1]).time[end])
    
    @lift xlims!(ax, $(x_start_obs), $(x_end_obs))
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
    # Clear the Observable cache when switching modes to ensure fresh Observables
    clear_observable_cache!(state.channels)
    update_channel_offsets!(state)
    draw(ax, state)
end

function draw(ax, state::DataBrowserState{<:AbstractDataState})
    _draw_implementation(ax, state, state.data)
end

# Type-specific implementations
function _draw_implementation(ax, state, data::ContinuousDataState)
    for (idx, visible) in enumerate(state.channels.visible)
        if visible
            col = state.channels.labels[idx]
            
            # Use cached Observables or create them if they don't exist
            current_data = data.current
            time_obs = get_or_create_cached_observable!(state.channels, :time, () -> @lift(@views $(data.current).time))
            
            # Create base data Observable without offset - use @views to avoid copying
            base_data_obs = get_or_create_cached_observable!(state.channels, Symbol("base_data_$(col)"), 
                () -> @lift(@views $(data.current)[!, $col]))
            
            # Create optimized data Observable with offset
            data_obs = create_optimized_data_observable(base_data_obs, col, state.view.offset[idx], state.view.display_scale)
            
            # Create optimized color Observable
            color_obs = create_optimized_color_observable(base_data_obs, state.view.crit_val)
            
            # Create reactive line style based on highlighting
            linewidth_obs = @lift($(state.channels.highlighted) == $col ? 4 : 2)
            alpha_obs = @lift($(state.channels.highlighted) == $col ? 1.0 : 0.7)
            
            state.channels.data_lines[col] = lines!(
                ax,
                time_obs,
                data_obs,
                color = color_obs,
                colormap = [:darkgrey, :darkgrey, :red],
                linewidth = linewidth_obs,
                alpha = alpha_obs,
                inspectable = true,  # Make lines clickable
            )
            
            if !state.view.butterfly[]
                # Create separate Observables for label position to reduce nesting
                label_x_obs = get_or_create_cached_observable!(state.channels, Symbol("label_x_$(col)"), 
                    () -> @lift($(data.current).time[$(state.view.xrange)[1]]))
                
                # Create intermediate Observable for the data point
                data_point = get_or_create_cached_observable!(state.channels, Symbol("data_point_$(col)"), 
                    () -> @lift($(data.current)[$(state.view.xrange)[1], $col]))
                label_y_obs = @lift($(data_point) .+ state.view.offset[idx])
                
                # Create reactive label style based on highlighting
                label_fontsize_obs = @lift($(state.channels.highlighted) == $col ? 22 : 18)
                label_color_obs = @lift($(state.channels.highlighted) == $col ? :red : :black)
                
                state.channels.data_labels[col] = text!(
                    ax,
                    label_x_obs,
                    label_y_obs,
                    text = String(col),
                    align = (:left, :center),
                    fontsize = label_fontsize_obs,
                    color = label_color_obs,
                )
            end
        end
    end
end

function _draw_implementation(ax, state, data::EpochedDataState)
    current_epoch = data.current_epoch[]  # Get current epoch value
    current_data = data.current[current_epoch]  # Get current DataFrame Observable

    # Pre-compute nrow outside Observable for better performance
    n_samples = nrow(current_data[])
    state.view.xrange[] = 1:min(state.view.xrange[][end], n_samples)

    for (idx, visible) in enumerate(state.channels.visible)
        if visible  # Only plot if channel is visible
            col = state.channels.labels[idx]

            # Use cached Observables or create them if they don't exist
            time_obs = get_or_create_cached_observable!(state.channels, Symbol("time_epoch_$(current_epoch)"), 
                () -> @lift($(current_data).time))
            
            # Create base data Observable without offset
            base_data_obs = get_or_create_cached_observable!(state.channels, Symbol("base_data_$(col)_epoch_$(current_epoch)"), 
                () -> @lift($(current_data)[!, $col]))
            
            # Create optimized data Observable with offset
            data_obs = create_optimized_data_observable(base_data_obs, col, state.view.offset[idx], state.view.display_scale)
            
            # Create optimized color Observable
            color_obs = create_optimized_color_observable(base_data_obs, state.view.crit_val)

            # Determine line style based on highlighting
            is_highlighted = state.channels.highlighted[] == col
            linewidth = is_highlighted ? 4 : 2
            alpha = is_highlighted ? 1.0 : 0.7

            state.channels.data_lines[col] = lines!(
                ax,
                time_obs,
                data_obs,
                color = color_obs,
                colormap = [:darkgrey, :darkgrey, :red],
                linewidth = linewidth,
                alpha = alpha,
                inspectable = true,  # Make lines clickable
            )
            
            if !state.view.butterfly[]
                # Create separate Observables for label position to reduce nesting
                label_x_obs = get_or_create_cached_observable!(state.channels, Symbol("label_x_$(col)_epoch_$(current_epoch)"), 
                    () -> @lift($(current_data).time[1]))
                
                # Create intermediate Observable for the data point
                data_point = get_or_create_cached_observable!(state.channels, Symbol("data_point_$(col)_epoch_$(current_epoch)"), 
                    () -> @lift($(current_data)[!, $col][1]))
                label_y_obs = @lift($(data_point) .+ state.view.offset[idx])
                
                # Create reactive label style based on highlighting
                label_fontsize_obs = @lift($(state.channels.highlighted) == $col ? 22 : 18)
                label_color_obs = @lift($(state.channels.highlighted) == $col ? :red : :black)
                
                state.channels.data_labels[col] = text!(
                    ax,
                    label_x_obs,
                    label_y_obs,
                    text = String(col),
                    align = (:left, :center),
                    fontsize = label_fontsize_obs,
                    color = label_color_obs,
                )
            end
        end
    end
end

# Generic draw_extra_channel function that dispatches based on data type
function draw_extra_channel!(ax, state::DataBrowserState{<:AbstractDataState})
    clear_axes!(ax, [state.extra_channel.data_lines, state.extra_channel.data_labels])
    _draw_extra_channel_implementation(ax, state, state.data)
end

# Type-specific implementations
function _draw_extra_channel_implementation(ax, state, data::ContinuousDataState)
    if state.extra_channel.visible && !isnothing(state.extra_channel.channel)
        current_offset = state.view.offset[end] + mean(diff(state.view.offset))
        channel = state.extra_channel.channel  # Get the channel symbol
        if eltype(data.current[][!, channel]) == Bool
            highlight_data = @views splitgroups(findall(data.current[][!, channel]))
            region_offset = all(iszero, highlight_data[2] .- highlight_data[1]) ? 5 : 0
            state.extra_channel.data_lines[channel] = vspan!(
                ax,
                data.current[][highlight_data[1], :time],
                data.current[][highlight_data[2].+region_offset, :time],
                color = :Red,
                alpha = 0.5,
                visible = true,
            )
        else
            # Use optimized Observable creation
            current_data = data.current
            time_obs = @lift($(data.current).time)
            channel_data_obs = create_optimized_data_observable(
                @lift($(data.current)[!, $channel]), 
                channel, 
                current_offset, 
                state.view.display_scale
            )
            
            state.extra_channel.data_lines[channel] = lines!(
                ax,
                time_obs,
                channel_data_obs,
                color = :black,
                linewidth = 2,
            )
            
            # Simplified label position Observable
            label_x_obs = @lift($(data.current).time[$(state.view.xrange)[1]])
            
            # Create intermediate Observable for the channel data point
            channel_data_point = @lift($(data.current)[!, $channel][$(state.view.xrange)[1]])
            label_y_obs = create_optimized_data_observable(
                channel_data_point, 
                channel, 
                current_offset, 
                state.view.display_scale
            )
            
            state.extra_channel.data_labels[channel] = text!(
                ax,
                label_x_obs,
                label_y_obs,
                text = String(channel),
                align = (:left, :center),
                fontsize = 18,
            )
        end
    end
end

function _draw_extra_channel_implementation(ax, state, data::EpochedDataState)
    if state.extra_channel.visible && !isnothing(state.extra_channel.channel)
        current_offset = state.view.offset[end] + mean(diff(state.view.offset))
        channel = state.extra_channel.channel  # Get the channel symbol
        current_epoch = state.data.current_epoch[]  # Get current epoch value
        current_data = state.data.current[current_epoch]  # Get current DataFrame Observable

        if eltype(current_data[][!, channel]) == Bool
            highlight_data = @views splitgroups(findall(current_data[][!, channel]))
            region_offset = all(iszero, highlight_data[2] .- highlight_data[1]) ? 5 : 0
            state.extra_channel.data_lines[channel] = vspan!(
                ax,
                current_data[][highlight_data[1], :time],
                current_data[][highlight_data[2].+region_offset, :time],
                color = :Red,
                alpha = 0.5,
                visible = true,
            )
        else
            # Use optimized Observable creation
            time_obs = @lift($(current_data).time)
            channel_data_obs = create_optimized_data_observable(
                @lift($(current_data)[!, $channel]), 
                channel, 
                current_offset, 
                state.view.display_scale
            )
            
            state.extra_channel.data_lines[channel] = lines!(
                ax,
                time_obs,
                channel_data_obs,
                color = :black,
                linewidth = 2,
            )
            
            # Simplified label position Observable
            label_x_obs = @lift($(current_data).time[1])
            
            # Create intermediate Observable for the channel data point
            channel_data_point = @lift($(current_data)[!, $channel][1])
            label_y_obs = create_optimized_data_observable(
                channel_data_point, 
                channel, 
                current_offset, 
                state.view.display_scale
            )
            
            state.extra_channel.data_labels[channel] = text!(
                ax,
                label_x_obs,
                label_y_obs,
                text = String(channel),
                align = (:left, :center),
                fontsize = 18,
            )
        end
    end
end

# Optimized zoom functions with minimal reactive overhead
function zoom_in!(state::DataBrowserState)
    """Zoom in by increasing display scale with minimal reactive overhead"""
    current_scale = state.view.display_scale[]
    new_scale = min(current_scale * 1.5, 10.0)  # Cap at 10x
    if new_scale != current_scale
        state.view.display_scale[] = new_scale
    end
end

function zoom_out!(state::DataBrowserState)
    """Zoom out by decreasing display scale with minimal reactive overhead"""
    current_scale = state.view.display_scale[]
    new_scale = max(current_scale / 1.5, 0.1)  # Floor at 0.1x
    if new_scale != current_scale
        state.view.display_scale[] = new_scale
    end
end

# Simple data Observable creation
function create_optimized_data_observable(data_obs::Observable, channel::Symbol, offset::Float64, display_scale_obs::Observable)
    """
    Create a simple data Observable that scales data by display scale.
    """
    # Create intermediate Observables to reduce nesting
    scaled_data = @lift($(data_obs) .* $(display_scale_obs))
    return @lift($(scaled_data) .+ offset)
end

# Optimized color Observable creation
function create_optimized_color_observable(data_obs::Observable, crit_val_obs::Observable)
    """
    Create an optimized color Observable with minimal reactive overhead.
    """
    # Create intermediate Observable for absolute values
    abs_data = @lift(abs.($(data_obs)))
    return @lift($(abs_data) .>= $(crit_val_obs))
end

# Helper function to clear cache when epoch changes
function clear_cache_on_epoch_change!(state::DataBrowserState{EpochedDataState})
    """
    Clear the Observable cache when the epoch changes for epoched data.
    This ensures fresh Observables are created for the new epoch.
    """
    clear_observable_cache!(state.channels)
end

# Helper functions for getting type-specific information
get_title(dat::EpochData) = "Epoch 1/$(n_epochs(dat))"
get_title(dat::ContinuousData) = ""

function plot_databrowser(dat::Union{ContinuousData,EpochData}, channel_labels::Vector{Symbol}, ica = nothing)

    # Check if CairoMakie is being used and warn about interactivity
    if string(Makie.current_backend()) == "CairoMakie"
        @minimal_warning "CairoMakie detected. For full interactivity in plot_databrowser, use GLMakie."
    end

    # Common setup
    fig = Figure(size = (1200, 800), fontsize = 18)

    # Create axis with appropriate title using multiple dispatch
    ax = Axis(fig[1, 1], xlabel = "Time (S)", ylabel = "Amplitude (μV)", title = get_title(dat))

    # Create appropriate state based on data type
    state = create_browser_state(dat, channel_labels, ax, ica)

    # Common UI setup
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

# Convenience methods
plot_databrowser(dat::Union{ContinuousData,EpochData}) = plot_databrowser(dat, dat.layout.label)
plot_databrowser(dat::Union{ContinuousData,EpochData}, ica::InfoIca) = plot_databrowser(dat, dat.layout.label, ica)
plot_databrowser(dat::Union{ContinuousData,EpochData}, channel_label::Symbol) = plot_databrowser(dat, [channel_label])
plot_databrowser(dat::Union{ContinuousData,EpochData}, channel_label::Symbol, ica::InfoIca) =
    plot_databrowser(dat, [channel_label], ica)
