const PLOT_DATABROWSER_KWARGS = Dict{Symbol,Tuple{Any,String}}(

    # Figure and layout
    :figure_padding => ((50, 30, 50, 50), "Figure padding as (left, right, bottom, top)"),

    # Axis styling
    :xlabel => ("Time (S)", "X-axis label"),
    :ylabel => ("Amplitude (μV)", "Y-axis label"),

    # UI styling
    :theme_fontsize => (nothing, "Font size for theme and UI elements"),

    # Line styling
    :channel_line_width => (1, "Line width for channel lines"),
    :channel_line_alpha => (1, "Transparency for channel lines"),
    :selected_channel_color => (:black, "Color for selected channels"),
    :unselected_channel_color => (:darkgrey, "Color for unselected channels"),
    :cleaned_channel_color => (:green, "Color for cleaned/corrected channels (Data - ICA)"),
    :show_cleaned_ica => (false, "Whether to show cleaned/corrected channels (Data - ICA) by default"),
    :channel_offset_scale => (1500, "Scale factor for channel vertical offset"),
    :channel_offset_margin => (0.9, "Margin factor for channel offset range"),

    # Navigation
    :scroll_step => (200, "Number of samples to scroll per arrow key press"),
    :zoom_step => (100, "Y-axis zoom increment per arrow key press"),

    # Selection styling
    :selection_color => (:blue, "Color for selection rectangle"),
    :selection_alpha => (0.1, "Transparency for selection rectangle"),

    # Filter parameters
    :default_hp_freq => (0.1, "Default high-pass filter frequency in Hz"),
    :default_lp_freq => (40.0, "Default low-pass filter frequency in Hz"),

    # Scale indicator
    :show_scale_indicator => (true, "Show scale indicator bar"),
    :scale_indicator_value => (100.0, "Scale indicator value in μV"),
    :scale_indicator_position => ((0.92, 0.96), "Scale indicator position as (x, y) in axis coordinates (0-1)"),

    # Display
    :display_plot => (true, "Whether to display the plot"),
)

# Base type for data states
"""Base type for databrowser data states."""
abstract type AbstractDataState end

# Data Browser: Continuous Data
"""Vertical marker line (trigger, EOG) on the databrowser axis."""
mutable struct Marker
    data::Any
    line::Any
    text::Any
    name::Symbol
    visible::Bool
    observers::Vector{Any}
    Marker(data, line, text, name, visible, observers = Any[]) = new(data, line, text, name, visible, observers)
end


"""Track user-selected time regions in the databrowser."""
mutable struct SelectionState
    active::Observable{Bool}
    bounds::Observable{Tuple{Float64,Float64}}
    visible::Observable{Bool}
    rectangle::Makie.Poly
    selected_regions::Observable{Vector{Tuple{Float64,Float64}}}  # Store multiple regions
    region_plots::Vector{Makie.Poly}  # Store plot objects for each region
    function SelectionState(ax, plot_kwargs)
        poly_element =
            poly!(ax, [Point2f(0.0, 0.0)], color = (plot_kwargs[:selection_color], plot_kwargs[:selection_alpha]), visible = false)
        new(Observable(false), Observable((0.0, 0.0)), Observable(false), poly_element, Observable([]), [])
    end
end

"""Observable high-pass / low-pass filter toggle state."""
mutable struct FilterState
    hp_active::Observable{Bool}
    lp_active::Observable{Bool}
    hp_freq::Observable{Float64}
    lp_freq::Observable{Float64}
    hp_method::Observable{String}
    lp_method::Observable{String}
    hp_order::Observable{Int}
    lp_order::Observable{Int}
    hp_func::Observable{String}
    lp_func::Observable{String}
    FilterState(plot_kwargs) = new(
        Observable(false),
        Observable(false),
        Observable(plot_kwargs[:default_hp_freq]),
        Observable(plot_kwargs[:default_lp_freq]),
        Observable("iir"),
        Observable("iir"),
        Observable(2),
        Observable(2),
        Observable("filtfilt"),
        Observable("filtfilt"),
    )
end

"""Label + action pair for a databrowser toggle button."""
struct ToggleConfig
    label::String
    action::Function
end

"""Databrowser view settings: x/y range, offsets, butterfly mode, scale."""
mutable struct ViewState
    xrange::Observable{UnitRange{Int64}}
    yrange::Observable{UnitRange{Int64}}
    offset::Vector{Float64}
    crit_val::Observable{Float64}
    butterfly::Observable{Bool}
    amplitude_scale::Observable{Float64}
    show_original_ica::Observable{Bool}
    show_subtracted_ica::Observable{Bool}
    function ViewState(n_channels::Int, n_samples::Int, plot_kwargs)
        offset_scale = plot_kwargs[:channel_offset_scale]
        offset_margin = plot_kwargs[:channel_offset_margin]
        offset =
            n_channels > 1 ? LinRange(offset_scale * offset_margin, -offset_scale * offset_margin, n_channels + 2)[2:(end-1)] :
            zeros(n_channels)
        new(
            Observable(1:n_samples),
            Observable(-1500:1500),
            offset,
            Observable(0.0),
            Observable(false),
            Observable(1.0),
            Observable(false),
            Observable(false),
        )
    end
end

"""Track channel visibility, selection, labels, and line objects."""
mutable struct ChannelState
    labels::Vector{Symbol}
    visible::Vector{Bool}
    selected::Vector{Bool}
    data_labels::Dict{Symbol,Makie.Text}
    data_lines::Dict{Symbol,Any}
    original_lines::Dict{Symbol,Any}    # ICA ghost: pre-removal signal
    subtracted_lines::Dict{Symbol,Any}  # ICA ghost: removed artefact
    draw_observers::Vector{Any}  # observables/observers created by `_draw`; disconnected on redraw
    function ChannelState(channel_labels::Vector{Symbol})
        new(
            channel_labels,
            fill(true, length(channel_labels)),
            fill(false, length(channel_labels)),
            Dict{Symbol,Makie.Text}(),
            Dict{Symbol,Any}(),
            Dict{Symbol,Any}(),
            Dict{Symbol,Any}(),
            Any[],
        )
    end
end

# Concrete data state implementations
"""Databrowser state wrapping continuous EEG data."""
mutable struct ContinuousDataState <: AbstractDataState
    current::Observable{EegData}
    original::EegData
    filter_state::FilterState
    function ContinuousDataState(data::EegData, plot_kwargs)
        new(Observable(copy(data)), data, FilterState(plot_kwargs))
    end
end

"""Databrowser state wrapping epoched EEG data."""
mutable struct EpochedDataState <: AbstractDataState
    current::Observable{EegData}
    original::EegData
    filter_state::FilterState
    current_epoch::Observable{Int}
    function EpochedDataState(data::EegData, plot_kwargs)
        new(Observable(copy(data)), data, FilterState(plot_kwargs), Observable(1))
    end
end

"""
    ExtraChannelInfo

A mutable struct to store additional channel information for the data browser.
Contains names, types, data dictionaries, and temporal alignment metadata.
"""
mutable struct ExtraChannelInfo
    channels::Vector{Symbol}
    data_lines::Dict{Symbol,Any}
    data_labels::Dict{Symbol,Makie.Text}
    fig::Union{Nothing,Figure}          # owning figure for legend placement
    legend::Union{Nothing,Any}          # current Legend object (replaced on redraw)
    draw_observers::Vector{Any}
    ExtraChannelInfo() = new(Symbol[], Dict{Symbol,Any}(), Dict{Symbol,Makie.Text}(), nothing, nothing, Any[])
end

"""
Visual state tracker for dynamic Independent Component Analysis (ICA) artifact removal.

`activation_cache` holds, per epoch (key `1` for continuous data), the activations of the removed
components (samples × removed components) for a sample range around the visible window, so scrolling
does not recompute them. It is invalidated whenever `removed_components` differs from `activation_key`.
"""
mutable struct IcaVisualState
    original::Union{Nothing,InfoIca}
    removed_components::Vector{Int}
    is_active::Bool
    activation_cache::Dict{Int,Tuple{UnitRange{Int},Matrix{Float64}}}
    activation_key::Vector{Int}

    function IcaVisualState(ica::Union{Nothing,InfoIca})
        new(ica, Int[], !isnothing(ica), Dict{Int,Tuple{UnitRange{Int},Matrix{Float64}}}(), Int[])
    end
end

"""Top-level composite state for the interactive data browser."""
mutable struct DataBrowserState{T<:AbstractDataState}
    view::ViewState
    channels::ChannelState
    data::T
    selection::SelectionState
    markers::Vector{Marker}
    marker_visibility::Dict{Symbol,Bool}
    ica::IcaVisualState
    extra_channel::ExtraChannelInfo
    reference_state::Symbol
    channel_repair_history::Vector{Tuple{Vector{Symbol},Symbol,Any}}  # (channels, method, original_data) - stack for multiple undos
    analysis_settings::Observable{AnalysisSettings}  # Observable analysis settings
    plot_kwargs::Dict{Symbol,Any}

    # Constructor
    function DataBrowserState{T}(;
        view::ViewState,
        channels::ChannelState,
        data::T,
        selection::SelectionState,
        ica_original::Union{Nothing,InfoIca} = nothing,
        extra_channel::ExtraChannelInfo = ExtraChannelInfo(),
        plot_kwargs::Dict{Symbol,Any} = Dict{Symbol,Any}(),
    ) where {T<:AbstractDataState}
        return new{T}(
            view,
            channels,
            data,
            selection,
            Marker[],
            Dict{Symbol,Bool}(),
            IcaVisualState(ica_original),
            extra_channel,
            data.original.analysis_info.reference,
            [],  # empty repair history
            Observable(AnalysisSettings()),  # Initialize with empty settings
            plot_kwargs,
        )
    end
end

# Type aliases 
const ContinuousDataBrowserState = DataBrowserState{ContinuousDataState}
const EpochedDataBrowserState = DataBrowserState{EpochedDataState}

# Analysis settings functions
"""
    _update_analysis_settings!(state::DataBrowserState)

Update the analysis settings observable with current state.
"""
function _update_analysis_settings!(state::DataBrowserState)
    # Get current filter settings
    hp_freq = state.data.filter_state.hp_active[] ? state.data.filter_state.hp_freq[] : 0.0
    lp_freq = state.data.filter_state.lp_active[] ? state.data.filter_state.lp_freq[] : 0.0

    # Get repaired channels and their method
    repaired_channels = Symbol[]
    for (channels, _, _) in state.channel_repair_history
        append!(repaired_channels, channels)
    end
    # If multiple repairs occurred, record the most recently applied method
    repair_method = isempty(state.channel_repair_history) ? :none : last(state.channel_repair_history)[2]

    # Get selected regions
    selected_regions = copy(state.selection.selected_regions[])

    # Update the observable
    state.analysis_settings[] = AnalysisSettings(
        hp_freq,
        lp_freq,
        state.reference_state,
        repaired_channels,
        repair_method,
        selected_regions,
        copy(state.ica.removed_components),
    )
end

"""Mark rows in `df` whose `:time` falls within any of `selected_regions`."""
function _mark_regions!(df::DataFrame, selected_regions::Vector{Tuple{Float64,Float64}})
    mask = falses(nrow(df))
    for (t0, t1) in selected_regions
        mask .|= (df.time .>= t0) .& (df.time .<= t1)
    end
    df[!, :selected_region] = mask
end

"""Add a boolean `:selected_region` column marking user-selected time regions."""
_add_selected_regions!(dat::EegData, regions::Vector{Tuple{Float64,Float64}}) = _mark_regions!(dat.data, regions)

function _add_selected_regions!(dat::EpochData, regions::Vector{Tuple{Float64,Float64}})
    foreach(df -> _mark_regions!(df, regions), dat.data)
end

# Data browser state creation
"""Construct a `DataBrowserState` from data, layout, axis, and optional ICA."""
function _create_browser_state(dat::T, channel_labels, ax, ica, plot_kwargs) where {T<:EegData}
    if !haskey(plot_kwargs, :show_cleaned_ica) || !(plot_kwargs[:show_cleaned_ica] isa Observable)
        val = get(plot_kwargs, :show_cleaned_ica, false)
        plot_kwargs[:show_cleaned_ica] = Observable{Bool}(val isa Observable ? val[] : val)
    end
    state_type = _data_state_type(T)
    initial_window = _get_initial_window_size(dat)
    return DataBrowserState{state_type}(
        view = ViewState(length(channel_labels), initial_window, plot_kwargs),
        channels = ChannelState(channel_labels),
        data = state_type(dat, plot_kwargs),  # Pass kwargs to data state constructor
        selection = SelectionState(ax, plot_kwargs),
        ica_original = ica,
        plot_kwargs = plot_kwargs,
    )
end

"""Return the initial x-range window size for the data browser."""
_get_initial_window_size(dat::ContinuousData) = min(5000, n_samples(dat))
_get_initial_window_size(dat::ErpData) = n_samples(dat)
_get_initial_window_size(dat::EpochData) = nrow(dat.data[1])

"""Map an `EegData` subtype to its corresponding `AbstractDataState` type."""
_data_state_type(::Type{ContinuousData}) = ContinuousDataState
_data_state_type(::Type{EpochData}) = EpochedDataState
_data_state_type(::Type{ErpData}) = ContinuousDataState

"""Return the current DataFrame for the active data state."""
_get_current_data(state::ContinuousDataState) = state.current[].data
_get_current_data(state::EpochedDataState) = state.current[].data[state.current_epoch[]]
"""Check whether a column exists in the current data."""
_has_column(state::ContinuousDataState, col::Symbol) = col in propertynames(state.current[].data)
_has_column(state::EpochedDataState, col::Symbol) = col in propertynames(state.current[].data[state.current_epoch[]])

######
# UI #
######
"""Create core UI widgets (toggles, menus, markers) shared by all data types."""
function _setup_ui_base(fig, ax, state, dat, ica = nothing)

    deregister_interaction!(ax, :rectanglezoom)
    _set_axes!(ax, state)

    # Mouse and keyboard events
    _handle_mouse_events!(ax, state)
    _handle_keyboard_events!(fig, ax, state)

    # Create toggles/markers/menus
    state.markers = _init_markers(ax, state)
    toggles = _create_toggles(fig, ax, state, dat)
    labels_menu = _create_labels_menu(fig, ax, state)
    reference_menu = _create_reference_menu(fig, state, dat)

    # Create optional menus (ica/extra channels)
    ica_menu = nothing
    if !isnothing(ica) && state.ica.is_active
        ica_menu = _create_ica_menu(fig, ax, state, ica)
    end

    extra_menu = nothing
    extra_labels_result = extra_labels(state.data.original)
    if !isnothing(extra_labels_result) && !isempty(extra_labels_result)
        extra_menu = _create_extra_channel_menu(fig, ax, state, dat)
    end

    return (toggles, labels_menu, reference_menu, ica_menu, extra_menu)
end

# Unified setup_ui method using multiple dispatch for the epoch menu
"""Assemble all UI controls into the figure grid."""
function _setup_ui(fig, ax, state::DataBrowserState{<:AbstractDataState}, dat, ica = nothing, plot_kwargs = nothing)
    # Get common UI elements
    toggles, labels_menu, reference_menu, ica_menu, extra_menu = _setup_ui_base(fig, ax, state, dat, ica)

    # Get type-specific epoch menu (or nothing)
    epoch_menu = _get_epoch_menu(fig, ax, state)

    # Build the grid components
    _build_grid_components!(fig, dat, state, toggles, labels_menu, reference_menu, ica_menu, extra_menu, epoch_menu)

    # Apply theme
    hideydecorations!(ax, label = true)

    return state
end

"""Build toggle buttons for butterfly, markers, and filter controls."""
function _create_toggles(fig, ax, state, dat)
    configs = [ToggleConfig("Butterfly Plot", (active) -> _butterfly_plot!(ax, state))]

    # Add marker toggles based on configuration
    marker_toggle_configs = [(:trigger, "Trigger"), (:is_vEOG, "vEOG"), (:is_hEOG, "hEOG")]

    for (marker_symbol, toggle_label) in marker_toggle_configs
        if _has_column(state.data, marker_symbol)
            state.marker_visibility[marker_symbol] = false
            push!(configs, ToggleConfig(toggle_label, (active) -> _toggle_marker_visibility!(ax, state, marker_symbol, active)))
        end
    end

    # Add filter toggles if available
    if state.data.original.analysis_info.hp_filter == 0.0
        push!(configs, ToggleConfig("HP-Filter", (active) -> _apply_hp_filter!(state, active)))
    end
    if state.data.original.analysis_info.lp_filter == 0.0
        push!(configs, ToggleConfig("LP-Filter", (active) -> _apply_lp_filter!(state, active)))
    end

    # Create toggles
    toggles_grid = []
    for config in configs
        tog = Toggle(fig)
        on(tog.active) do active
            config.action(active)
        end

        if config.label == "Trigger"
            btn = Button(fig, label = "Trigger", halign = :left)
            on(btn.clicks) do _
                _show_trigger_menu(state, ax, :trigger)
            end
            push!(toggles_grid, [tog, btn])
        elseif config.label == "HP-Filter" || config.label == "LP-Filter"
            btn = Button(fig, label = config.label, halign = :left)
            filter_type = config.label == "HP-Filter" ? :hp : :lp
            on(btn.clicks) do _
                _show_single_filter_menu(state, dat, filter_type)
            end
            push!(toggles_grid, [tog, btn])
        else
            lbl = Label(fig, config.label, fontsize = 22, halign = :left)
            push!(toggles_grid, [tog, lbl])
        end
    end

    # Return as grid layout
    return permutedims(reduce(hcat, toggles_grid))

end

"""Toggle visibility for a marker type by name and track visibility in state."""
function _toggle_marker_visibility!(ax, state, marker_symbol::Symbol, active::Bool)
    state.marker_visibility[marker_symbol] = active
    marker_index = findfirst(m -> m.name == marker_symbol, state.markers)
    if !isnothing(marker_index)
        _plot_vertical_lines!(ax, state.markers[marker_index], active)
    end
end

"""Popup: select which triggers to display."""
function _show_trigger_menu(state, ax, marker_symbol)
    trigger_col = if state.data isa ContinuousDataState
        state.data.original.data[!, marker_symbol]
    else
        reduce(vcat, (df[!, marker_symbol] for df in state.data.original.data))
    end
    unique_triggers = sort(unique(filter(!iszero, trigger_col)))

    if !haskey(state.plot_kwargs, :active_triggers)
        state.plot_kwargs[:active_triggers] = Set(unique_triggers)
    end
    active_set = state.plot_kwargs[:active_triggers]

    menu_fig = _eegfun_figure(size = (600, 600))
    Label(menu_fig[1, 1], "Select Triggers to Display", fontsize = 18, halign = :center)

    scroll_area = menu_fig[2, 1] = GridLayout()
    cbs = _build_checkbox_grid!(scroll_area, string.(unique_triggers), 5, (trig, _) -> (parse(Float64, trig) in active_set))

    groups = [("All", ch -> true), ("None", ch -> false)]
    _add_group_buttons!(menu_fig, 3, cbs, string.(unique_triggers), groups)

    action_area = menu_fig[4, 1] = GridLayout()
    btn_apply = Button(action_area[1, 1], label = "Apply", width = 200)
    on(btn_apply.clicks) do _
        empty!(active_set)
        for (i, cb) in enumerate(cbs)
            if cb.checked[]
                push!(active_set, unique_triggers[i])
            end
        end
        # Force redraw with new selected triggers
        _update_markers!(ax, state)

        # Ensure the toggle matches the new state (if any triggers are selected, turn the main toggle ON)
        # Note: In the future, we could sync the main Toggle UI state, but for now just updating markers works.
    end

    _display_in_screen(menu_fig)
end

"""Popup: select filter cutoff frequency for a specific filter."""
function _show_single_filter_menu(state, dat, filter_type::Symbol)
    menu_fig = _eegfun_figure(size = (450, 350))

    title_str = filter_type == :hp ? "Highpass Filter Settings" : "Lowpass Filter Settings"
    Label(menu_fig[1, 1], title_str, fontsize = 18, halign = :center, font = :bold)

    grid = menu_fig[2, 1] = GridLayout(valign = :top)
    row = 1

    if filter_type == :hp
        range_freq = 0.1:0.1:2
        freq_field = :hp_freq
        method_field = :hp_method
        order_field = :hp_order
        func_field = :hp_func
        default_val = 0.5
    else
        range_freq = 5:5:60
        freq_field = :lp_freq
        method_field = :lp_method
        order_field = :lp_order
        func_field = :lp_func
        default_val = 20.0
    end

    fs = state.data.filter_state
    if getfield(fs, freq_field)[] == 0.0
        getfield(fs, freq_field)[] = default_val
    end

    # Cutoff
    Label(grid[row, 1], "Cutoff (Hz):", halign = :left)
    slider_freq = Slider(grid[row, 2], range = range_freq, startvalue = getfield(fs, freq_field)[])
    Label(grid[row, 3], @lift(string(round($(slider_freq.value), digits = 1))), halign = :left)
    on(slider_freq.value) do val

        getfield(fs, freq_field)[] = val
    end
    row += 1

    # Order
    Label(grid[row, 1], "Order:", halign = :left)
    slider_order = Slider(grid[row, 2], range = 1:10, startvalue = getfield(fs, order_field)[])
    Label(grid[row, 3], @lift(string($(slider_order.value))), halign = :left)
    on(slider_order.value) do val

        getfield(fs, order_field)[] = val
    end
    row += 1

    # Method
    Label(grid[row, 1], "Method:", halign = :left)
    method_labels = ["Butterworth (IIR)", "FIR (Hamming)"]
    method_values = ["iir", "fir"]
    current_method = getfield(fs, method_field)[]
    menu_method = Menu(
        grid[row, 2],
        options = zip(method_labels, method_values),
        default = current_method == "iir" ? "Butterworth (IIR)" : "FIR (Hamming)",
        width = 220,
    )
    on(menu_method.selection) do val

        getfield(fs, method_field)[] = val
    end
    row += 1

    # Function
    Label(grid[row, 1], "Function:", halign = :left)
    func_labels = ["filtfilt (zero-phase)", "filt (one-pass)"]
    func_values = ["filtfilt", "filt"]
    current_func = getfield(fs, func_field)[]
    menu_func = Menu(
        grid[row, 2],
        options = zip(func_labels, func_values),
        default = current_func == "filtfilt" ? "filtfilt (zero-phase)" : "filt (one-pass)",
        width = 220,
    )
    on(menu_func.selection) do val

        getfield(fs, func_field)[] = val
    end
    row += 1

    rowgap!(grid, 10)
    colgap!(grid, 15)

    action_area = menu_fig[3, 1] = GridLayout()
    btn_apply = Button(action_area[1, 1], label = "Apply", width = 150)
    on(btn_apply.clicks) do _
        _apply_filters!(state)
    end

    _display_in_screen(menu_fig)
end


"""Popup: select which channels are visible."""
function _show_labels_menu(state, ax)
    all_channels = state.channels.labels
    menu_fig = _eegfun_figure(size = (900, 700))
    Label(menu_fig[1, 1], "Select Channels to Display", fontsize = 18, halign = :center)
    scroll_area = menu_fig[2, 1] = GridLayout()
    cbs = _build_checkbox_grid!(scroll_area, all_channels, 9, (ch, i) -> state.channels.visible[i])

    groups = [
        ("All", ch -> true),
        ("None", ch -> false),
        ("Left", ch -> occursin(r"\d*[13579]$", String(ch))),
        ("Right", ch -> occursin(r"\d*[24680]$", String(ch))),
        ("Central", ch -> occursin(r"z$", String(ch))),
    ]
    _add_group_buttons!(menu_fig, 3, cbs, all_channels, groups)

    action_area = menu_fig[4, 1] = GridLayout()
    btn_apply = Button(action_area[1, 1], label = "Apply", width = 200)
    on(btn_apply.clicks) do _
        for (i, cb) in enumerate(cbs)
            state.channels.visible[i] = cb.checked[]
        end
        _clear_and_save_limits!(ax, state)
        _update_channel_offsets!(state)
        _draw(ax, state)
        if !isempty(state.extra_channel.channels)
            _draw_extra_channel!(ax, state)
        end
    end

    _display_in_screen(menu_fig)
end

"""Create the channel-label selection button for the control panel."""
_create_labels_menu(fig, ax, state) = _popup_button(fig, "Labels", () -> _show_labels_menu(state, ax))

"""Popup: select reference channels and apply re-referencing."""
function _show_reference_menu(state, dat)
    all_channels = state.channels.labels
    menu_fig = _eegfun_figure(size = (900, 700))
    Label(menu_fig[1, 1], "Select Reference Channels", fontsize = 18, halign = :center)

    # Determine currently-checked channels from stored state
    current_syms = if state.reference_state == :avg
        all_channels  # all checked when avg
    elseif state.reference_state == :mastoid
        [:M1, :M2]
    elseif state.reference_state == :none
        Symbol[]
    else
        Symbol.(split(string(state.reference_state), "_"))
    end

    scroll_area = menu_fig[2, 1] = GridLayout()
    cbs = _build_checkbox_grid!(scroll_area, all_channels, 9, (ch, _) -> ch in current_syms)

    groups = [
        ("All (Avg)", ch -> true),
        ("None", ch -> false),
        ("Left", ch -> occursin(r"\d*[13579]$", String(ch))),
        ("Right", ch -> occursin(r"\d*[24680]$", String(ch))),
        ("Central", ch -> occursin(r"z$", String(ch))),
        ("Mastoids", ch -> ch == :M1 || ch == :M2),
    ]
    _add_group_buttons!(menu_fig, 3, cbs, all_channels, groups)

    action_area = menu_fig[4, 1] = GridLayout()
    btn_apply = Button(action_area[1, 1], label = "Apply", width = 200)
    on(btn_apply.clicks) do _
        selected = [all_channels[i] for (i, cb) in enumerate(cbs) if cb.checked[]]
        opt = if isempty(selected)
            :none
        elseif length(selected) == length(all_channels)
            :avg
        elseif length(selected) == 2 && :M1 in selected && :M2 in selected
            :mastoid
        else
            selected
        end
        state.reference_state = Symbol(opt isa Symbol ? string(opt) : join(opt, "_"))
        _rebuild_current!(state)
    end

    _display_in_screen(menu_fig)
end

"""Create the re-reference button for the control panel."""
_create_reference_menu(fig, state, dat) = _popup_button(fig, "Reference", () -> _show_reference_menu(state, dat))

"""Apply re-referencing to `cur` according to `ref_state`."""
function _apply_reference!(cur::EegData, ref_state::Symbol)
    ref_state == :none && return cur
    ref = if ref_state == :avg
        :avg
    elseif ref_state == :mastoid
        :mastoid
    else
        syms = Symbol.(split(string(ref_state), "_"))
        length(syms) == 1 ? syms[1] : syms
    end
    rereference!(cur, ref)
    return cur
end

"""Apply the isolated artifact of a single ICA component linearly to a Continuous dataset."""
function _apply_single_ica_removal!(cur::ContinuousData, orig::ContinuousData, ica::InfoIca, comp::Int)
    dat = cur.data
    cols = ica.layout.data.label

    unmix_vec = ica.unmixing[comp, :]
    mix_vec = ica.mixing[:, comp]
    n_samples = nrow(dat)

    # Calculate 1D activation series for the component
    activation = zeros(Float64, n_samples)
    for (i, col_sym) in enumerate(cols)
        norm_ch = (orig.data[!, col_sym] .- ica.mean[i]) ./ ica.scale
        activation .+= unmix_vec[i] .* norm_ch
    end

    # Dynamically inject the isolated footprint back into sensor space
    for (i, col_sym) in enumerate(cols)
        artifact_ch = mix_vec[i] .* ica.scale .* activation
        dat[!, col_sym] .-= artifact_ch
    end
end

"""Apply the isolated artifact of a single ICA component linearly to an Epoched dataset."""
function _apply_single_ica_removal!(cur::EpochData, orig::EpochData, ica::InfoIca, comp::Int)
    cols = ica.layout.data.label
    unmix_vec = ica.unmixing[comp, :]
    mix_vec = ica.mixing[:, comp]

    for ep = 1:length(cur.data)
        dat_ep = cur.data[ep]
        orig_ep = orig.data[ep]
        n_samples = nrow(dat_ep)

        activation = zeros(Float64, n_samples)
        for (i, col_sym) in enumerate(cols)
            norm_ch = (orig_ep[!, col_sym] .- ica.mean[i]) ./ ica.scale
            activation .+= unmix_vec[i] .* norm_ch
        end

        for (i, col_sym) in enumerate(cols)
            artifact_ch = mix_vec[i] .* ica.scale .* activation
            dat_ep[!, col_sym] .-= artifact_ch
        end
    end
end

"""Apply all active ICA component removals to `cur`."""
function _apply_all_ica_removals!(cur::EegData, state)
    if state.ica.is_active && !isempty(state.ica.removed_components)
        for comp in state.ica.removed_components
            _apply_single_ica_removal!(cur, state.data.original, state.ica.original, comp)
        end
    end
    return cur
end

"""Apply active high-pass and low-pass filters to `cur`."""
function _apply_active_filters!(cur::EegData, state)
    fs = state.data.filter_state
    if fs.hp_active[]
        sym_method = fs.hp_method[] isa Symbol ? fs.hp_method[] : Symbol(fs.hp_method[])
        sym_func = fs.hp_func[] isa Symbol ? fs.hp_func[] : Symbol(fs.hp_func[])
        highpass_filter!(cur, fs.hp_freq[]; filter_method = sym_method, order = fs.hp_order[], filter_func = sym_func)
    end
    if fs.lp_active[]
        sym_method = fs.lp_method[] isa Symbol ? fs.lp_method[] : Symbol(fs.lp_method[])
        sym_func = fs.lp_func[] isa Symbol ? fs.lp_func[] : Symbol(fs.lp_func[])
        lowpass_filter!(cur, fs.lp_freq[]; filter_method = sym_method, order = fs.lp_order[], filter_func = sym_func)
    end
    return cur
end

"""
Rebuild `state.data.current` from `state.data.original` in the canonical sequence:
repairs -> reference -> ICA removals -> active filters.
Updates `state.data.current[]` with a single notification and synchronizes analysis settings.
"""
function _rebuild_current!(state)
    cur = copy(state.data.original)
    for (chs, method, _) in state.channel_repair_history
        repair_channels!(cur, chs; method = method)
    end
    _apply_reference!(cur, state.reference_state)
    _apply_all_ica_removals!(cur, state)
    _apply_active_filters!(cur, state)
    state.data.current[] = cur
    _update_analysis_settings!(state)
    return cur
end

"""Toggle an ICA component for removal and rebuild dataset."""
function _toggle_ica_component!(state, component_id::Union{Nothing,Int})
    if isnothing(component_id)
        empty!(state.ica.removed_components)
        _rebuild_current!(state)
        return
    end

    is_removing = !(component_id in state.ica.removed_components)
    if is_removing
        push!(state.ica.removed_components, component_id)
    else
        filter!(x -> x != component_id, state.ica.removed_components)
    end

    _rebuild_current!(state)
end

"""Popup: toggle ICA component removals and display options."""
function _show_ica_menu(state, ax, ica)
    all_components = ica.ica_label
    menu_fig = _eegfun_figure(size = (800, 600))
    Label(menu_fig[1, 1], "Select ICA Components to Remove", fontsize = 18, halign = :center)

    scroll_area = menu_fig[2, 1] = GridLayout()
    cbs = _build_checkbox_grid!(scroll_area, all_components, 9, (comp, _) -> begin
        n = _extract_int(String(comp))
        !isnothing(n) && n in state.ica.removed_components
    end)

    _add_group_buttons!(menu_fig, 3, cbs, all_components, [("None (Reset)", _ -> false)])

    action_area = menu_fig[4, 1] = GridLayout()
    btn_apply = Button(action_area[1, 1], label = "Apply", width = 200)
    on(btn_apply.clicks) do _
        _clear_and_save_limits!(ax, state)
        empty!(state.ica.removed_components)
        for (i, cb) in enumerate(cbs)
            if cb.checked[]
                n = _extract_int(String(all_components[i]))
                !isnothing(n) && push!(state.ica.removed_components, n)
            end
        end
        _rebuild_current!(state)
        _draw(ax, state)
    end

    # Display options row
    display_area = menu_fig[5, 1] = GridLayout()
    Label(display_area[1, 1:6], "Display Options", fontsize = 16, halign = :center)
    display_items = [
        (state.view.show_original_ica, "Data", :grey),
        (state.plot_kwargs[:show_cleaned_ica], "Data - ICA", get(state.plot_kwargs, :cleaned_channel_color, :green)),
        (state.view.show_subtracted_ica, "ICA Activation", :red),
    ]
    for (j, (obs, lbl, colr)) in enumerate(display_items)
        col = (j - 1) * 2 + 1
        cb = Checkbox(display_area[2, col], checked = obs[], width = 16, height = 16)
        Label(display_area[2, col+1], lbl, halign = :left, color = colr)
        on(cb.checked) do active
            obs[] = active
        end
        colgap!(display_area, col, 5)
        j < length(display_items) && colgap!(display_area, col + 1, 30)
    end

    _display_in_screen(menu_fig)
end

"""Create the ICA components button for the control panel."""
function _create_ica_menu(fig, ax, state, ica)
    if !haskey(state.plot_kwargs, :show_cleaned_ica) || !(state.plot_kwargs[:show_cleaned_ica] isa Observable)
        state.plot_kwargs[:show_cleaned_ica] = Observable(false)
    end
    state.view.show_original_ica[] = true
    state.view.show_subtracted_ica[] = false
    return _popup_button(fig, "ICA Components", () -> _show_ica_menu(state, ax, ica))
end

"""Create the epoch navigation slider and label."""
function _create_epoch_menu(fig, ax, state)
    slider_epoch = Slider(fig[2, 1], range = 1:n_epochs(state.data.original), startvalue = state.data.current_epoch[], snap = true)
    label = Label(
        fig,
        @lift("Epoch: $($(slider_epoch.value))/$(n_epochs(state.data.original))"),
        fontsize = 22,
        halign = :left,
        tellwidth = false,
    )

    # Flag to prevent circular updates
    updating_from_keyboard = Ref(false)

    # Handle slider input
    on(slider_epoch.value) do epoch_num
        if !updating_from_keyboard[]
            state.data.current_epoch[] = epoch_num
            ax.title = "Epoch $(epoch_num)/$(n_epochs(state.data.original))"
            _update_markers!(ax, state)
            if !isempty(state.extra_channel.channels)
                _draw_extra_channel!(ax, state)
            end
        end
    end

    # Make slider observe current_epoch changes (for left/right key navigation)
    on(state.data.current_epoch) do epoch_num
        updating_from_keyboard[] = true
        set_close_to!(slider_epoch, epoch_num)
        updating_from_keyboard[] = false
    end

    return hcat(slider_epoch, label)
end

"""Find trigger onset events within the selected region for continuous data."""
function _find_region_triggers(state::ContinuousDataBrowserState, clicked_region_idx = nothing)
    bounds = _get_selected_region_bounds(state, clicked_region_idx)
    isnothing(bounds) && return NamedTuple{(:sample, :code, :time),Tuple{Int,Int,Float64}}[]
    x_min, x_max = bounds

    df = state.data.current[].data
    !hasproperty(df, :trigger) && return NamedTuple{(:sample, :code, :time),Tuple{Int,Int,Float64}}[]

    active_triggers = get(state.plot_kwargs, :active_triggers, nothing)
    return [
        (sample = i, code = Int(df.trigger[i]), time = Float64(df.time[i])) for i = 1:nrow(df) if df.time[i] >= x_min &&
            df.time[i] <= x_max &&
            df.trigger[i] != 0 &&
            (isnothing(active_triggers) || df.trigger[i] in active_triggers) &&
            (i == 1 || df.trigger[i] != df.trigger[i-1])
    ]
end

_find_region_triggers(state, clicked_region_idx = nothing) = NamedTuple{(:sample, :code, :time),Tuple{Int,Int,Float64}}[]

"""Extract and plot a single epoch around a trigger event within the selected region bounds."""
function _plot_single_trigger_epoch(state::ContinuousDataBrowserState, ev, x_min::Real, x_max::Real; kwargs...)
    curr_dat = state.data.current[]
    idx_min = searchsortedfirst(curr_dat.data.time, x_min)
    idx_max = searchsortedlast(curr_dat.data.time, x_max)
    t_trig = ev.time

    cols_to_keep = setdiff(propertynames(curr_dat.data), [:condition, :condition_name, :file])
    epoch_df = copy(curr_dat.data[idx_min:idx_max, cols_to_keep])
    epoch_df.time .-= t_trig
    if !hasproperty(epoch_df, :epoch)
        insertcols!(epoch_df, 1, :epoch => 1)
    end

    single_epoch = EpochData(
        curr_dat.file,
        1,
        "Trigger_$(ev.code)",
        [epoch_df],
        copy(curr_dat.layout),
        curr_dat.sample_rate,
        copy(curr_dat.analysis_info),
    )

    selected_channels = state.channels.labels[state.channels.visible]
    return plot_epochs(single_epoch; channel_selection = channels(selected_channels), kwargs...)
end

"""Extract and plot all epochs matching a trigger code using the selected region offset interval."""
function _plot_all_trigger_epochs(state::ContinuousDataBrowserState, code::Integer, x_min::Real, x_max::Real, t_trig::Real; kwargs...)
    curr_dat = state.data.current[]
    interval = (x_min - t_trig, x_max - t_trig)
    epoch_cond = EpochCondition(name = "Trigger_$code", trigger_sequences = [[code]])
    try
        all_epochs = extract_epochs(curr_dat, 1, epoch_cond, interval)
        selected_channels = state.channels.labels[state.channels.visible]
        return plot_epochs(all_epochs; channel_selection = channels(selected_channels), kwargs...)
    catch e
        @minimal_warning "Could not extract all epochs for trigger $code: $e"
        return nothing
    end
end

"""Open a popup menu with Topoplot / Spectrum / Epoch / Region actions."""
function _show_additional_menu(state, clicked_region_idx = nothing)

    # Check if layout has valid coordinates for topographic plots
    has_coords = has_valid_coordinates(state.data.original.layout)

    menu_items = Tuple{String,Function}[]

    if has_coords
        push!(menu_items, ("Topoplot", () -> begin
            selected_data = _subset_selected_data(state, clicked_region_idx)
            !isnothing(selected_data) && plot_topography(selected_data; scrubber = true)
        end))
        push!(menu_items, ("Topoplot 3D", () -> begin
            selected_data = _subset_selected_data(state, clicked_region_idx)
            !isnothing(selected_data) && plot_topography_3d(selected_data; scrubber = true)
        end))
    end
    push!(menu_items, ("Spectrum", () -> begin
        selected_data = _subset_selected_data(state, clicked_region_idx)
        !isnothing(selected_data) && plot_channel_spectrum(selected_data)
    end))

    # Check for trigger codes in the selected region for continuous data
    # Only offer epoch plotting if the trigger toggle is active and the region contains an actual trigger
    trigger_active = get(state.marker_visibility, :trigger, false)
    if state.data isa ContinuousDataState && trigger_active
        trigger_events = _find_region_triggers(state, clicked_region_idx)
        bounds = _get_selected_region_bounds(state, clicked_region_idx)
        if !isempty(trigger_events) && !isnothing(bounds)
            x_min, x_max = bounds
            unique_codes = unique([ev.code for ev in trigger_events])

            # Specific epoch action(s)
            if length(trigger_events) == 1
                ev = trigger_events[1]
                push!(menu_items, ("Plot Epoch (Trigger $(ev.code))", () -> _plot_single_trigger_epoch(state, ev, x_min, x_max)))
            else
                for (i, ev) in enumerate(trigger_events)
                    lbl =
                        length(unique_codes) == 1 ? "Plot Epoch $i (Trigger $(ev.code))" :
                        "Plot Epoch (Trigger $(ev.code) at $(round(ev.time, digits = 2))s)"
                    push!(menu_items, (lbl, () -> _plot_single_trigger_epoch(state, ev, x_min, x_max)))
                end
            end

            # All matching epochs action(s)
            for code in unique_codes
                ev = trigger_events[findfirst(e -> e.code == code, trigger_events)]
                push!(menu_items, ("Plot All Epochs (Trigger $code)", () -> _plot_all_trigger_epochs(state, code, x_min, x_max, ev.time)))
            end
        end
    end

    push!(menu_items, ("Get Selected Regions", () -> begin
        selected_regions_bool = _get_selected_regions_bool(state)
        @info "Selected regions boolean vector:"
        @info "Length: $(length(selected_regions_bool))"
        @info "Number of selected samples: $(sum(selected_regions_bool))"
        @info "Selected regions: $(state.selection.selected_regions[])"
    end))

    menu_fig = _eegfun_figure(size = (340, max(160, 50 * length(menu_items))))
    for (idx, (label_str, action_fn)) in enumerate(menu_items)
        btn = Button(menu_fig[idx, 1], label = label_str, tellwidth = false)
        on(btn.clicks) do _
            action_fn()
        end
    end

    _display_in_screen(menu_fig; size = (340, max(160, 50 * length(menu_items))))
end

"""Open the interactive channel repair window with checkboxes and method selection."""
function _channel_repair_menu(state, selected_channels, ax)
    all_channels = state.channels.labels

    menu_fig = _eegfun_figure(size = (900, 700))
    Label(menu_fig[1, 1], "Channel Repair Interface", fontsize = 18, halign = :center)

    scroll_area = menu_fig[2, 1] = GridLayout()
    channel_checkboxes = _build_checkbox_grid!(
        scroll_area,
        all_channels,
        9,
        (ch, _) -> begin
            is_repaired = _is_channel_repaired(state, ch)
            # Pre-select channels passed in, show repair status via label colour
            ch in selected_channels
        end;
        fontsize = 14,
        label_fn = ch -> _is_channel_repaired(state, ch) ? "$(ch) ✓" : string(ch),
    )

    # Add repair method selection
    Label(menu_fig[3, 1], "Repair Method:", fontsize = 14)
    method_area = menu_fig[4, 1] = GridLayout()
    method_buttons = [
        Button(method_area[1, 1], label = "Neighbor Interpolation", width = 200),
        Button(method_area[1, 2], label = "Spherical Spline", width = 200),
    ]

    # Add action buttons
    action_area = menu_fig[5, 1] = GridLayout()
    action_buttons =
        [Button(action_area[1, 1], label = "Apply Repair", width = 200), Button(action_area[1, 2], label = "Undo Last Repair", width = 200)]

    # Method selection via shared radio-button helper
    selected_method = Observable(:neighbor_interpolation)
    _radio_buttons!(method_buttons, selected_method, [:neighbor_interpolation, :spherical_spline])

    # Apply repair
    on(action_buttons[1].clicks) do n
        channels_to_repair = all_channels[findall(cb -> cb.checked[], channel_checkboxes)]
        if !isempty(channels_to_repair)
            _repair_selected_channels!(state, channels_to_repair, selected_method[], ax)
        else
            @info "No channels selected for repair"
        end
    end

    # Undo last repair
    on(action_buttons[2].clicks) do n
        if !isempty(state.channel_repair_history)
            _undo_last_repair!(state, ax)
        else
            @info "No repairs to undo"
        end
    end

    _display_in_screen(menu_fig)
end

"""Interpolate the selected channels, store in repair history, and rebuild."""
function _repair_selected_channels!(state, selected_channels, method, ax)
    # Check if any of these channels have already been repaired
    already_repaired = Set{Symbol}()
    for (repaired_channels, _, _) in state.channel_repair_history
        for ch in repaired_channels
            if ch in selected_channels
                push!(already_repaired, ch)
            end
        end
    end

    if !isempty(already_repaired)
        @info "Channels $(join(string.(collect(already_repaired)), ", ")) have already been repaired. Please undo first or select different channels."
        return
    end

    # Store repair in history
    push!(state.channel_repair_history, (selected_channels, method, nothing))

    _rebuild_current!(state)

    # Clear and redraw the plot
    _clear_and_save_limits!(ax, state)
    _draw(ax, state)

    total_repairs = length(state.channel_repair_history)
    @info "Successfully repaired channels: $(join(string.(selected_channels), ", ")) using $method"
    @info "Total repairs in history: $total_repairs"
end

"""Pop the last channel repair from history and rebuild."""
function _undo_last_repair!(state, ax)
    if isempty(state.channel_repair_history)
        @info "No repairs to undo"
        return
    end

    channels, method, _ = pop!(state.channel_repair_history)

    _rebuild_current!(state)

    # Clear and redraw the plot
    _clear_and_save_limits!(ax, state)
    _draw(ax, state)

    remaining_repairs = length(state.channel_repair_history)
    @info "Undid repair of channels: $(join(string.(channels), ", ")) (was $method)"
    @info "Remaining repairs in history: $remaining_repairs"
end

"""Create extreme-value slider. Returns a list of `hcat` rows for stacking."""
function _create_extreme_slider(fig, state)
    slider_extreme = Slider(fig, range = 0:5:100, startvalue = 0, width = 100)
    on(slider_extreme.value) do x
        state.view.crit_val[] = x
    end
    return [hcat(slider_extreme, Label(fig, @lift(string("Extreme: ", $(slider_extreme.value), " μV")), fontsize = 22))]
end

"""Create sliders for continuous data — includes x-range and position."""
function _create_sliders(fig, state::ContinuousDataBrowserState, dat)
    sliders = state.ica.is_active ? [] : _create_extreme_slider(fig, state)

    # Add navigation sliders specific to continuous data
    slider_range = Slider(fig[3, 1], range = 100:50:30000, startvalue = length(state.view.xrange[]), snap = true)
    slider_x = Slider(fig[2, 1], range = 1:50:nrow(state.data.current[].data), startvalue = first(state.view.xrange[]), snap = true)

    updating_from_key = Ref(false)

    on(slider_range.value) do w
        updating_from_key[] && return
        start = slider_x.value[]
        new_range = start:min(nrow(state.data.current[].data), start+w-1)
        if length(new_range) > 1
            state.view.xrange[] = new_range
        end
    end

    on(slider_x.value) do start
        updating_from_key[] && return
        w = slider_range.value[]
        new_range = start:min(nrow(state.data.current[].data), start+w-1)
        if length(new_range) > 1
            state.view.xrange[] = new_range
        end
    end

    on(state.view.xrange) do rng
        updating_from_key[] = true
        set_close_to!(slider_x, first(rng))
        set_close_to!(slider_range, length(rng))
        updating_from_key[] = false
    end

    return sliders
end

"""Create sliders for epoched data (extreme slider only)."""
function _create_sliders(fig, state::EpochedDataBrowserState, dat)
    return state.ica.is_active ? [] : _create_extreme_slider(fig, state)
end

"""Popup: select which extra channels to overlay."""
function _show_extra_channel_menu(state, ax, dat)
    all_extras = extra_labels(dat)
    (isnothing(all_extras) || isempty(all_extras)) && return

    menu_fig = _eegfun_figure(size = (900, 700))
    Label(menu_fig[1, 1], "Select Extra Channels", fontsize = 18, halign = :center)

    scroll_area = menu_fig[2, 1] = GridLayout()
    cbs = _build_checkbox_grid!(scroll_area, all_extras, 9, (ch, _) -> ch in state.extra_channel.channels)

    action_area = menu_fig[3, 1] = GridLayout()
    btn_apply = Button(action_area[1, 1], label = "Apply", width = 150)
    on(btn_apply.clicks) do _
        empty!(state.extra_channel.channels)
        for (i, cb) in enumerate(cbs)
            cb.checked[] && push!(state.extra_channel.channels, all_extras[i])
        end
        _draw_extra_channel!(ax, state)
    end

    _display_in_screen(menu_fig)
end

"""Create the extra-channels button for the control panel."""
_create_extra_channel_menu(fig, ax, state, dat) = _popup_button(fig, "Extra Channels", () -> _show_extra_channel_menu(state, ax, dat))

"""Wrap a widget so it can be stacked with `vcat` in `grid!`."""
_as_grid_row(fig, w::Button) = hcat(w, Label(fig, "", tellwidth = false))
_as_grid_row(_, w) = w  # already a matrix (toggles, sliders, hcat menus)

function _build_grid_components!(
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
    grid_components = Matrix{Any}[]
    push!(grid_components, toggles[:, 1:2])

    append!(grid_components, _create_sliders(fig, state, dat))
    push!(grid_components, _as_grid_row(fig, labels_menu))
    push!(grid_components, _as_grid_row(fig, reference_menu))

    !isnothing(ica_menu) && push!(grid_components, _as_grid_row(fig, ica_menu))
    !isnothing(extra_menu) && push!(grid_components, _as_grid_row(fig, extra_menu))
    !isnothing(epoch_menu) && push!(grid_components, _as_grid_row(fig, epoch_menu))

    control_panel = grid!(vcat(grid_components...), tellheight = false)
    fig[1, 2] = control_panel

    colsize!(fig.layout, 2, Relative(0.25))
    rowsize!(fig.layout, 1, Relative(1.0))

end

"""Return an epoch menu for epoched data, `nothing` for continuous."""
_get_epoch_menu(fig, ax, state::ContinuousDataBrowserState) = nothing
_get_epoch_menu(fig, ax, state::EpochedDataBrowserState) = _create_epoch_menu(fig, ax, state)

############
# Navigation
############
const KEYBOARD_ACTIONS = Dict(Keyboard.left => :left, Keyboard.right => :right, Keyboard.up => :up, Keyboard.down => :down)

"""Dispatch arrow-key actions to pan / zoom / epoch-step."""
function _handle_navigation!(ax, state::DataBrowserState{<:AbstractDataState}, action::Symbol)
    if action == :up
        _ymore!(ax, state)
    elseif action == :down
        _yless!(ax, state)
    elseif action == :left
        _handle_left_navigation(ax, state)
    elseif action == :right
        _handle_right_navigation(ax, state)
    end
end

"""Left/right arrow: pan for continuous, epoch-step for epoched."""
_handle_left_navigation(ax, state::ContinuousDataBrowserState) = _xback!(ax, state)
_handle_left_navigation(ax, state::EpochedDataBrowserState) = _step_epoch_backward(ax, state)
_handle_right_navigation(ax, state::ContinuousDataBrowserState) = _xforward!(ax, state)
_handle_right_navigation(ax, state::EpochedDataBrowserState) = _step_epoch_forward(ax, state)

"""Scroll the x-range backward."""
function _xback!(ax, state::ContinuousDataBrowserState)
    step = state.plot_kwargs[:scroll_step]
    current_range = state.view.xrange.val
    current_range[1] <= 1 && return
    shift = min(step, current_range[1] - 1)
    state.view.xrange[] = current_range .- shift
    _restore_axis_limits!(ax, state)
end

"""Scroll the x-range forward."""
function _xforward!(ax, state::ContinuousDataBrowserState)
    step = state.plot_kwargs[:scroll_step]
    current_range = state.view.xrange.val
    max_idx = nrow(state.data.current[].data)
    current_range[end] >= max_idx && return
    shift = min(step, max_idx - current_range[end])
    state.view.xrange[] = current_range .+ shift
    _restore_axis_limits!(ax, state)
end

"""Step one epoch backward / forward."""
_step_epoch_backward(ax, state::EpochedDataBrowserState) = _step_epoch!(ax, state, -1)
_step_epoch_forward(ax, state::EpochedDataBrowserState) = _step_epoch!(ax, state, 1)

"""Change the current epoch by `direction` (+1 or -1) and redraw."""
function _step_epoch!(ax, state::EpochedDataBrowserState, direction::Int)
    current = state.data.current_epoch[]
    total = n_epochs(state.data.original)
    new_epoch = clamp(current + direction, 1, total)
    new_epoch == current && return
    state.data.current_epoch[] = new_epoch
    ax.title = "Epoch $(new_epoch)/$total"
    _update_markers!(ax, state)
    if !isempty(state.extra_channel.channels)
        _draw_extra_channel!(ax, state)
    end
end

"""Zoom y-axis out / in."""
_yless!(ax, state) = _yzoom!(ax, state, 1.2)
_ymore!(ax, state) = _yzoom!(ax, state, 0.8)

"""Scale the y-range or amplitude by `factor`."""
function _yzoom!(ax, state, factor::Float64)
    zoom_step = state.plot_kwargs[:zoom_step]
    if state.view.butterfly[]
        # In butterfly mode: adjust y-range
        y_min, y_max = state.view.yrange.val[1], state.view.yrange.val[end]
        if factor > 1.0  # Zoom in (yless)
            (y_min + zoom_step >= 0 || y_max - zoom_step <= 0) && return
            state.view.yrange[] = (y_min+zoom_step):(y_max-zoom_step)
        else  # Zoom out (ymore)
            state.view.yrange[] = (y_min-zoom_step):(y_max+zoom_step)
        end
        ylims!(ax, state.view.yrange.val[1], state.view.yrange.val[end])
    else # In non-butterfly mode: adjust amplitude scale
        state.view.amplitude_scale[] = state.view.amplitude_scale[] * factor
    end
end

"""Return `true` if the mouse position falls inside the browser axis bbox."""
function _is_mouse_in_browser_axis(ax, pos)
    bbox = ax.layoutobservables.computedbbox[]
    return bbox.origin[1] <= pos[1] <= (bbox.origin[1] + bbox.widths[1]) && bbox.origin[2] <= pos[2] <= (bbox.origin[2] + bbox.widths[2])
end

"""Return the index of the selected region containing `mouse_x`, or `nothing`."""
function _find_clicked_region(state, mouse_x)
    # Check if click is within any of the selected regions
    regions = state.selection.selected_regions[]
    for (i, (start_time, end_time)) in enumerate(regions)
        if mouse_x >= min(start_time, end_time) && mouse_x <= max(start_time, end_time)
            return i
        end
    end
    return nothing
end

"""Delete region at `region_idx` from the selection list and remove its plot."""
function _remove_region_from_selection!(ax, state, region_idx)
    # Remove the region from the list
    regions = state.selection.selected_regions[]
    if 1 <= region_idx <= length(regions)
        # Remove the region plot
        plot_to_remove = state.selection.region_plots[region_idx]
        delete!(ax.scene, plot_to_remove)
        deleteat!(state.selection.region_plots, region_idx)

        # Remove the region from the list
        deleteat!(regions, region_idx)
        state.selection.selected_regions[] = regions
        _update_analysis_settings!(state)
    end
end

"""Begin a new time-region selection at `mouse_x`."""
function _start_selection!(ax, state, mouse_x)
    state.selection.active[] = true
    state.selection.bounds[] = (mouse_x, mouse_x)
    _update_x_region_selection!(ax, state, mouse_x, mouse_x)
end

"""Finalise the current drag selection and add it to the region list."""
function _finish_selection!(ax, state, mouse_x)
    t1 = state.selection.bounds[][1]
    t2 = mouse_x
    state.selection.active[] = false
    _clear_x_region_selection!(state)

    if abs(t2 - t1) < 1e-4
        return
    end

    state.selection.bounds[] = (t1, t2)
    _add_region_to_selection!(ax, state, t1, t2)
end

"""Wire up mouse button and drag events for selection and channel clicking."""
function _handle_mouse_events!(ax, state)
    # Track if Shift and Ctrl are currently pressed
    shift_pressed = Ref(false)
    ctrl_pressed = Ref(false)

    on(events(ax).keyboardbutton) do key_event
        if key_event.key == Keyboard.left_shift
            shift_pressed[] = key_event.action == Keyboard.press
        elseif key_event.key == Keyboard.left_control
            ctrl_pressed[] = key_event.action == Keyboard.press
        end
    end

    on(events(ax).mousebutton) do event
        pos = events(ax).mouseposition[]
        if !_is_mouse_in_browser_axis(ax, pos)
            return
        end

        mouse_x = mouseposition(ax)[1]

        if event.button == Mouse.left
            if event.action == Mouse.press
                if ctrl_pressed[]
                    mouse_y = mouseposition(ax)[2]
                    clicked_channel = _find_closest_browser_channel(ax, state, mouse_x, mouse_y)
                    if !isnothing(clicked_channel)
                        _toggle_channel_visibility!(ax, state, clicked_channel)
                        return
                    end
                elseif shift_pressed[]
                    _handle_left_click!(ax, state, event, mouse_x)
                end
            elseif event.action == Mouse.release
                if state.selection.active[] || shift_pressed[]
                    _handle_left_click!(ax, state, event, mouse_x)
                end
            end
        elseif event.button == Mouse.right && event.action == Mouse.press
            _handle_right_click!(ax, state, mouse_x)
        end
    end

    # Update selection rectangle while dragging
    on(events(ax).mouseposition) do _
        if state.selection.active[]
            world_pos = mouseposition(ax)[1]
            _update_x_region_selection!(ax, state, state.selection.bounds[][1], world_pos)
        end
    end
end

"""Handle left-button press/release for region creation or removal."""
function _handle_left_click!(ax, state, event, mouse_x)
    if event.action == Mouse.press
        # Check if click is within any existing selected region
        clicked_region_idx = _find_clicked_region(state, mouse_x)
        if !isnothing(clicked_region_idx)
            # Remove the clicked region
            _remove_region_from_selection!(ax, state, clicked_region_idx)
        else
            # Start a new selection
            _start_selection!(ax, state, mouse_x)
        end
    elseif event.action == Mouse.release && state.selection.active[]
        _finish_selection!(ax, state, mouse_x)
    end
end

"""Handle right-click: show context menu if inside a selected region."""
function _handle_right_click!(ax, state, mouse_x)
    # Check if right-click is within any selected region
    clicked_region_idx = _find_clicked_region(state, mouse_x)
    if !isnothing(clicked_region_idx)
        _show_additional_menu(state, clicked_region_idx)
    end
    # Right-click outside regions does nothing (use 'r' key for channel repair)
end

"""Open the channel repair menu (continuous data only)."""
function _show_channel_repair_menu(state::DataBrowserState{<:ContinuousDataState}, ax)
    if !has_valid_coordinates(state.data.original.layout)
        @minimal_warning "Channel repair requires electrode layout coordinates. Please load a layout file (e.g. using EegFun.read_layout)."
        return
    end
    _channel_repair_menu(state, _get_selected_channels(state), ax)
end

function _show_channel_repair_menu(state::DataBrowserState{<:EpochedDataState}, ax)
    @info "Channel repair is only available for continuous data"
end

"""Return currently Ctrl-selected channel labels."""
function _get_selected_channels(state)
    selected_indices = findall(state.channels.selected)
    return state.channels.labels[selected_indices]
end

"""Find the channel line closest to the given mouse coordinates."""
function _find_closest_browser_channel(ax, state, mouse_x, mouse_y)
    current_data = _get_current_data(state.data)
    tolerance = 10  # data units tolerance for snap detection

    min_distance = Inf
    closest_channel = nothing

    # Find the data point closest to the clicked x position 
    time_data = current_data.time
    x_idx = searchsortedfirst(time_data, mouse_x)
    if x_idx > 1 && x_idx <= length(time_data) && abs(time_data[x_idx-1] - mouse_x) < abs(time_data[x_idx] - mouse_x)
        x_idx -= 1
    end
    x_idx = clamp(x_idx, 1, length(time_data))

    for (idx, visible) in enumerate(state.channels.visible)
        if visible
            col = state.channels.labels[idx]

            # Calculate the y position of this channel at the clicked x position
            channel_y = current_data[x_idx, col] * state.view.amplitude_scale[] + state.view.offset[idx]

            # Calculate distance to the clicked point
            distance = abs(mouse_y - channel_y)

            # Track closest channel overall
            if distance < min_distance
                min_distance = distance
                closest_channel = idx
            end
        end
    end

    # Only return if within tolerance
    return min_distance <= tolerance ? closest_channel : nothing
end

"""Toggle the selected/highlighted state of one channel instantly."""
function _toggle_channel_selection!(ax, state, channel_idx)
    col = state.channels.labels[channel_idx]
    is_selected = state.channels.selected[channel_idx] = !state.channels.selected[channel_idx]
    is_repaired = _is_channel_repaired(state, col)

    if haskey(state.channels.data_lines, col) && haskey(state.channels.data_labels, col)
        line_handle = state.channels.data_lines[col]
        label_handle = state.channels.data_labels[col]

        if is_repaired || is_selected
            sel_color = state.plot_kwargs[:selected_channel_color]
            line_handle.colormap[] = [sel_color, sel_color, sel_color]
            line_handle.linewidth[] = state.plot_kwargs[:channel_line_width] * 2
            label_handle.color[] = :red
        else
            is_cleaned = state.ica.is_active && haskey(state.plot_kwargs, :show_cleaned_ica) && state.plot_kwargs[:show_cleaned_ica][]
            base_color = is_cleaned ? get(state.plot_kwargs, :cleaned_channel_color, :green) : state.plot_kwargs[:unselected_channel_color]
            line_handle.colormap[] = [base_color, base_color, :red]
            line_handle.linewidth[] = state.plot_kwargs[:channel_line_width]
            label_handle.color[] = :black
        end
    else
        # Fallback if handles do not exist yet
        _clear_and_save_limits!(ax, state)
        _draw(ax, state)
    end
end
const _toggle_channel_visibility! = _toggle_channel_selection!

"""Register keyboard shortcuts (arrows, i=help, r=repair, c=clear)."""
function _handle_keyboard_events!(fig, ax, state)
    on(events(fig).keyboardbutton) do event
        if event.action == Keyboard.press && event.key == Keyboard.i
            # Show help for databrowser
            _show_plot_help(:databrowser)
        elseif event.action == Keyboard.press && event.key == Keyboard.r
            # Open channel repair menu
            _show_channel_repair_menu(state, ax)
        elseif event.action == Keyboard.press && event.key == Keyboard.c
            # Clear all selected regions
            _clear_all_selected_regions!(ax, state)
        elseif event.action in (Keyboard.press, Keyboard.repeat) && haskey(KEYBOARD_ACTIONS, event.key)
            action = KEYBOARD_ACTIONS[event.key]
            _handle_navigation!(ax, state, action)
        end
    end
end


"""Redraw the temporary selection rectangle between x1 and x2."""
function _update_x_region_selection!(ax, state, x1, x2)
    limits_y = ax.limits[][2]
    ylims = isnothing(limits_y) ? (Float64(state.view.yrange[][1]), Float64(state.view.yrange[][end])) : limits_y
    state.selection.rectangle[1] = Point2f[
        Point2f(Float64(x1), Float64(ylims[1])),
        Point2f(Float64(x2), Float64(ylims[1])),
        Point2f(Float64(x2), Float64(ylims[2])),
        Point2f(Float64(x1), Float64(ylims[2])),
    ]
    state.selection.rectangle.visible[] = true
end

"""Store a finalised region, draw its permanent overlay, and update settings."""
function _add_region_to_selection!(ax, state, x1, x2)
    # Ensure x1 <= x2
    if x1 > x2
        x1, x2 = x2, x1
    end

    # Add to selected regions
    current_regions = state.selection.selected_regions[]
    new_region = (x1, x2)
    push!(current_regions, new_region)
    state.selection.selected_regions[] = current_regions

    # Create a permanent region plot
    limits_y = ax.limits[][2]
    ylims = isnothing(limits_y) ? (Float64(state.view.yrange[][1]), Float64(state.view.yrange[][end])) : limits_y
    region_points = Point2f[
        Point2f(Float64(x1), Float64(ylims[1])),
        Point2f(Float64(x2), Float64(ylims[1])),
        Point2f(Float64(x2), Float64(ylims[2])),
        Point2f(Float64(x1), Float64(ylims[2])),
    ]
    region_plot = poly!(
        ax,
        region_points,
        color = (state.plot_kwargs[:selection_color], state.plot_kwargs[:selection_alpha]),
        strokecolor = :transparent,
    )
    push!(state.selection.region_plots, region_plot)

    # Update analysis settings
    _update_analysis_settings!(state)

end

"""Reset the temporary selection rectangle to invisible."""
function _clear_x_region_selection!(state)
    # Set to a single point instead of empty vector to avoid CairoMakie issues
    state.selection.rectangle[1] = [Point2f(0.0, 0.0)]
    state.selection.bounds[] = (0.0, 0.0)
    state.selection.visible[] = false
    state.selection.rectangle.visible[] = false
end

"""Delete all permanent region overlays and clear the region list."""
function _clear_all_selected_regions!(ax, state)
    # Clear all region plots
    for plot in state.selection.region_plots
        delete!(ax.scene, plot)
    end
    empty!(state.selection.region_plots)

    # Clear the selected regions list
    state.selection.selected_regions[] = Tuple{Float64,Float64}[]
    _update_analysis_settings!(state)
end

"""
    _get_selected_regions_bool(state::DataBrowserState) -> Vector{Bool}

Returns a boolean vector indicating which samples are within the selected regions.
The vector has the same length as the total number of samples in the data.
"""
function _get_selected_regions_bool(state::DataBrowserState)
    current_data = _get_current_data(state.data)
    total_samples = nrow(current_data)
    time_data = current_data.time
    bool_vector = falses(total_samples)

    for (start_time, end_time) in state.selection.selected_regions[]
        # Find the bounding sample indices in O(log N) without allocations
        start_idx = searchsortedfirst(time_data, start_time)
        end_idx = searchsortedlast(time_data, end_time)

        # Ensure indices are within bounds and start <= end
        start_idx = max(1, min(start_idx, total_samples))
        end_idx = max(1, min(end_idx, total_samples))
        if start_idx > end_idx
            start_idx, end_idx = end_idx, start_idx
        end

        # Mark the region as selected
        bool_vector[start_idx:end_idx] .= true
    end

    return bool_vector
end

"""
    _get_selected_regions_info(state::DataBrowserState) -> NamedTuple

Returns detailed information about the selected regions including:
- `bool_vector`: Boolean vector of selected samples
- `regions`: List of (start_time, end_time) tuples
- `n_samples`: Number of selected samples
- `n_regions`: Number of selected regions
"""
function _get_selected_regions_info(state::DataBrowserState)
    bool_vector = _get_selected_regions_bool(state)
    regions = state.selection.selected_regions[]

    return (bool_vector = bool_vector, regions = regions, n_samples = sum(bool_vector), n_regions = length(regions))
end

"""Return (x_min, x_max) tuple for the clicked or most recent selected region."""
function _get_selected_region_bounds(state, clicked_region_idx = nothing)
    if !isnothing(clicked_region_idx) && 1 <= clicked_region_idx <= length(state.selection.selected_regions[])
        r = state.selection.selected_regions[][clicked_region_idx]
        return (min(r[1], r[2]), max(r[1], r[2]))
    elseif !isempty(state.selection.selected_regions[])
        r = state.selection.selected_regions[][end]
        return (min(r[1], r[2]), max(r[1], r[2]))
    else
        return nothing
    end
end

"""Return a subset of the current continuous data within the selected region."""
function _subset_selected_data(state::ContinuousDataBrowserState, clicked_region_idx = nothing)
    bounds = _get_selected_region_bounds(state, clicked_region_idx)
    if isnothing(bounds)
        @minimal_warning "No region selected"
        return nothing
    end
    x_min, x_max = bounds

    selected_channels = state.channels.labels[state.channels.visible]
    return subset(
        state.data.current[],
        sample_selection = x -> (x.time .>= x_min) .& (x.time .<= x_max),
        channel_selection = channels(selected_channels),
    )
end

"""Return a subset of the current epoch data within the selected region."""
function _subset_selected_data(state::EpochedDataBrowserState, clicked_region_idx = nothing)
    bounds = _get_selected_region_bounds(state, clicked_region_idx)
    if isnothing(bounds)
        @minimal_warning "No region selected"
        return nothing
    end
    x_min, x_max = bounds

    selected_channels = state.channels.labels[state.channels.visible]
    current_epoch = state.data.current_epoch[]

    # Create a sample_selection function that works with epoch DataFrame
    # get_selected_samples for MultiDataFrameEeg expects a function that takes the first epoch DataFrame
    sample_selection = epoch_df -> begin
        time_mask = (epoch_df.time .>= x_min) .& (epoch_df.time .<= x_max)
        return time_mask
    end

    epoch_data = subset(
        state.data.current[],
        sample_selection = sample_selection,
        channel_selection = channels(selected_channels),
        epoch_selection = epochs([current_epoch]),
    )

    # Check if we have data before converting
    if isempty(epoch_data)
        @minimal_warning "Selected region contains no data"
        return nothing
    end

    return epoch_to_continuous(epoch_data, 1)
end


############
# Filtering
############
############
# Filtering
############
"""Rebuild current dataset from original, applying repairs, reference, ICA, and active filters."""
_apply_filters!(state) = _rebuild_current!(state)

"""Toggle or set the high-pass filter and reapply."""
function _apply_hp_filter!(state, active::Bool = !state.data.filter_state.hp_active[])
    state.data.filter_state.hp_active[] = active
    _apply_filters!(state)
end

"""Toggle or set the low-pass filter and reapply."""
function _apply_lp_filter!(state, active::Bool = !state.data.filter_state.lp_active[])
    state.data.filter_state.lp_active[] = active
    _apply_filters!(state)
end

########################
# Drawing
########################

"""Add a vertical marker (trigger/EOG) to the marker list, window-scoped for performance."""
function _add_marker!(markers, ax, state, data, col; label = nothing, trial = nothing, visible = false, active_values = nothing)
    if isnothing(trial)
        # More efficient: filter directly without findall
        mask = data[!, col] .!= 0
        if !isnothing(active_values)
            mask .&= map(v -> v in active_values, data[!, col])
        end
        marker_data = data[mask, [:time, col]]
    else
        mask = data[trial][!, col] .!= 0
        if !isnothing(active_values)
            mask .&= map(v -> v in active_values, data[trial][!, col])
        end
        marker_data = data[trial][mask, [:time, col]]
    end

    # if no markers, return
    nrow(marker_data) == 0 && return

    label_vec = isnothing(label) ? string.(marker_data[!, col]) : repeat([label], nrow(marker_data))

    # Reactive atomic window-scoped marker time, label, & position subviews
    marker_obs = @lift(begin
        if state isa EpochedDataBrowserState
            idxs = 1:nrow(marker_data)
        else
            win = $(state.view.xrange)
            t_start = data.time[win[1]]
            t_end = data.time[win[end]]
            times = marker_data.time
            i1 = searchsortedfirst(times, t_start)
            i2 = searchsortedlast(times, t_end)
            idxs = i1:i2
        end
        times = marker_data.time[idxs]
        lbls = label_vec[idxs]
        y_top = Float32($(ax.yaxis.attributes.limits)[2] * 0.98)
        pos = Point2f[Point2f(x, y_top) for x in times]
        (times, lbls, pos)
    end)

    visible_times = @lift($marker_obs[1])
    visible_labels = @lift($marker_obs[2])
    visible_positions = @lift($marker_obs[3])

    push!(
        markers,
        Marker(
            marker_data,
            vlines!(ax, visible_times, color = :grey, linewidth = 1, visible = visible),
            text!(
                ax,
                visible_labels,
                position = visible_positions,
                space = :data,
                align = (:center, :center),
                fontsize = 22,
                visible = visible,
            ),
            col,
            visible,
            Any[marker_obs, visible_times, visible_labels, visible_positions],
        ),
    )
end

"""Disconnect an observer or observable from its upstream sources."""
function _disconnect_observer(o)
    if o isa Observable
        Makie.Observables.clear(o)
    elseif o isa AbstractVector
        foreach(_disconnect_observer, o)
    else
        Makie.off(o)
    end
end

"""Return the four line/label dict collections for clearing."""
_all_line_dicts(state) =
    [state.channels.data_lines, state.channels.data_labels, state.channels.original_lines, state.channels.subtracted_lines]

"""Check whether a channel has been repaired in the current session."""
_is_channel_repaired(state, col) = any(col in chs for (chs, _, _) in state.channel_repair_history)

"""Recompute channel vertical offsets for the current visible channels."""
function _update_channel_offsets!(state)
    nchannels = count(state.channels.visible)
    margin = state.plot_kwargs[:channel_offset_margin]
    if nchannels > 1 && !state.view.butterfly[]
        y_max = state.view.yrange[][end] * margin
        y_min = state.view.yrange[][1] * margin
        step = (y_min - y_max) / (nchannels - 1)
        visible_indices = findall(state.channels.visible)
        for (i, idx) in enumerate(visible_indices)
            state.view.offset[idx] = y_max + (i - 1) * step
        end
    else
        state.view.offset[state.channels.visible] .= 0.0
    end
end

"""Delete all plotted lines/labels from the given dict collections."""
function _clear_axes!(ax, datas)
    for data in datas
        for value in values(data)
            delete!(ax, value)
        end
        empty!(data)
    end
end

"""Set reactive x/y limits on the axis from view state."""
function _set_axes!(ax, state::DataBrowserState{<:AbstractDataState})
    @lift ylims!(ax, $(state.view.yrange)[1], $(state.view.yrange)[end])
    _set_x_limits!(ax, state, state.data)
end

"""Bind x-limits to the continuous data time range."""
function _set_x_limits!(ax, state, data::ContinuousDataState)
    @lift xlims!(ax, data.current[].data.time[$(state.view.xrange)[1]], data.current[].data.time[$(state.view.xrange)[end]])
end

"""Bind x-limits to the epoch time range."""
function _set_x_limits!(ax, state, data::EpochedDataState)
    @lift xlims!(ax, $(data.current).data[1].time[1], $(data.current).data[1].time[end])
end

"""Create initial marker objects for triggers and EOG channels."""
function _init_markers(ax, state)
    markers = Marker[]
    data = _get_current_data(state.data)

    # Define marker configurations
    marker_configs = [(:trigger, nothing), (:is_vEOG, "v"), (:is_hEOG, "h")]

    # Add markers based on configuration
    for (symbol, label) in marker_configs
        if _has_column(state.data, symbol)
            active_vals = symbol == :trigger ? get(state.plot_kwargs, :active_triggers, nothing) : nothing
            _add_marker!(
                markers,
                ax,
                state,
                data,
                symbol,
                label = label,
                visible = get(state.marker_visibility, symbol, false),
                active_values = active_vals,
            )
        end
    end

    return markers
end

"""Re-create markers for the current epoch (preserving visibility)."""
function _update_markers!(ax, state)
    for marker in state.markers
        for o in marker.observers
            _disconnect_observer(o)
        end
        delete!(ax, marker.line)
        delete!(ax, marker.text)
    end
    empty!(state.markers)
    state.markers = _init_markers(ax, state)
end

"""Toggle butterfly mode on/off and redraw channels."""
function _butterfly_plot!(ax, state)
    state.view.butterfly[] = !state.view.butterfly[]
    _clear_and_save_limits!(ax, state)
    _update_channel_offsets!(state)
    _draw(ax, state)
    if !isempty(state.extra_channel.channels)
        _draw_extra_channel!(ax, state)
    end
end

# ICA ghost-line helpers
"""Fill `pts` in place with `(x[i], y[i] * scale + offset)`; `pts`, `x` and `y` must have the same length."""
function _window_points!(pts::Vector{Point2{Float64}}, x, y, scale, offset)
    @inbounds for i in eachindex(pts)
        pts[i] = Point2{Float64}(x[i], y[i] * scale + offset)
    end
    return pts
end

"""Fill `pts` in place with NaN y-values (an invisible line segment) at the x-positions `x`."""
function _nan_points!(pts::Vector{Point2{Float64}}, x)
    @inbounds for i in eachindex(pts)
        pts[i] = Point2{Float64}(x[i], NaN)
    end
    return pts
end

"""
    _points_observable(fill_points!, xrange, deps...)

Observable holding one `Vector{Point2{Float64}}` per line (x and y update together, so Makie sees a single
notification per change). The vector is reused and mutated in place; it is only reallocated when the window
length changes. `fill_points!(buf, resized)` fills `buf` and returns `false` if nothing changed (no notification).
The observer functions it registers are appended to `tracker` (see `_disconnect_draw_observers!`).
"""
function _points_observable(fill_points!::F, tracker::Vector{Any}, xrange::Observable, deps::Observable...) where {F<:Function}
    obs = Observable(Vector{Point2{Float64}}(undef, length(xrange[])))
    resized = Ref(true)
    function refresh(_...)
        n = length(xrange[])
        buf = obs[]
        if length(buf) != n
            buf = Vector{Point2{Float64}}(undef, n)
            resized[] = true
        end
        if fill_points!(buf, resized[])
            resized[] = false
            obs[] = buf   # same array mutated in place: Observables always notify on assignment
        end
        return nothing
    end
    refresh()
    push!(tracker, onany(refresh, xrange, deps...))   # tracked so a redraw can disconnect them
    return obs
end

"""Raw (unfiltered) original series for `col` from the browser's original data."""
_ica_original_series(data::ContinuousDataState, ep, col) = data.original.data[!, col]
_ica_original_series(data::EpochedDataState, ep, col) = data.original.data[ep][!, col]

"""Accumulate `u * ((x[i] - m) / scale)` for `i in rng` into column `k` of `acts` (type-stable kernel)."""
function _accumulate_activation!(
    acts::Matrix{Float64},
    k::Int,
    x::AbstractVector{<:Real},
    rng::UnitRange{Int},
    u::Float64,
    m::Float64,
    scale::Float64,
)
    @inbounds for (j, i) in enumerate(rng)
        acts[j, k] += u * ((x[i] - m) / scale)
    end
    return acts
end

"""Activations (samples × removed components) of the removed ICA components for sample range `rng`."""
function _compute_ica_activations(state, ep, rng::UnitRange{Int})
    ica = state.ica.original::InfoIca
    comps = state.ica.removed_components
    acts = zeros(Float64, length(rng), length(comps))
    for (k, comp) in enumerate(comps), (ci, col_sym) in enumerate(ica.layout.data.label)
        x = _ica_original_series(state.data, ep, col_sym)
        _accumulate_activation!(acts, k, x, rng, ica.unmixing[comp, ci], ica.mean[ci], ica.scale)
    end
    return acts
end

"""
Return `(range, acts)` where `acts` holds the removed-component activations covering `win`.

Results are cached per epoch and extended to a margin of one window width on each side, so
consecutive scroll steps reuse them. The cache is dropped when the removed components change.
"""
function _ica_activations!(state, ep, win::UnitRange{Int})
    ica = state.ica
    if ica.activation_key != ica.removed_components
        empty!(ica.activation_cache)
        ica.activation_key = copy(ica.removed_components)
    end
    entry = get(ica.activation_cache, ep, nothing)
    if isnothing(entry) || first(win) < first(entry[1]) || last(win) > last(entry[1])
        n_total = length(_ica_original_series(state.data, ep, first(ica.original.layout.data.label)))
        w = length(win)
        rng = max(1, first(win)-w):min(n_total, last(win)+w)
        entry = (rng, _compute_ica_activations(state, ep, rng))
        ica.activation_cache[ep] = entry
    end
    return entry::Tuple{UnitRange{Int},Matrix{Float64}}
end

"""Add `mix_weight * scale * acts[:, k]` (rows of `win` offset by `off`) onto `restored` in place."""
function _add_component!(restored::Vector{Float64}, acts::Matrix{Float64}, k::Int, win::UnitRange{Int}, off::Int, coeff::Float64)
    @inbounds for (j, i) in enumerate(win)
        restored[j] += coeff * acts[i-off, k]
    end
    return restored
end

"""Current-window values with the removed components added back (the 'original' ICA signal) for one channel."""
function _ica_restored_window(state, ep, win::UnitRange{Int}, current_raw, ch_row)
    restored = collect(Float64, current_raw)
    comps = state.ica.removed_components
    isempty(comps) && return restored
    rng, acts = _ica_activations!(state, ep, win)
    ica = state.ica.original::InfoIca
    off = first(rng) - 1
    for (k, comp) in enumerate(comps)
        # (mix_weight * scale) first, then times the activation: same operation order as the original code
        _add_component!(restored, acts, k, win, off, ica.mixing[ch_row, comp] * ica.scale)
    end
    return restored
end

# Single function with data access abstraction
"""Draw all visible channel lines and labels on the axis."""
function _draw(ax, state::DataBrowserState{<:AbstractDataState})
    # Get data access functions based on type
    get_data, get_time, get_label_y = _get_data_accessors(state.data)

    epoch_obs = state.data isa EpochedDataState ? state.data.current_epoch : Observable(1)

    # Every derived observable / observer created below is registered in `state.channels.draw_observers`,
    # so `_clear_and_save_limits!` can disconnect them from the shared inputs (xrange, current, ...) on redraw.
    track(x) = (push!(state.channels.draw_observers, x); x)

    # Pre-compute shared observables
    visible_time_obs = track(@lift(get_time($(state.data.current), $(epoch_obs), $(state.view.xrange))))
    time_start_obs   = track(@lift($(visible_time_obs)[1]))  # Fix 6: derived from visible_time_obs

    ch_idx_map = state.ica.is_active ? Dict(lbl => i for (i, lbl) in enumerate(state.ica.original.layout.data.label)) : Dict{Symbol,Int}()

    for (idx, visible) in enumerate(state.channels.visible)
        col = state.channels.labels[idx]
        if visible
            is_selected = state.channels.selected[idx]

            ch_row = get(ch_idx_map, col, 0)

            # Channel data (compute once)
            channel_data_obs = track(
                @lift(
                    begin
                        cur = $(state.data.current)
                        ep = $(epoch_obs)
                        win = $(state.view.xrange)
                        raw = get_data(cur, ep, win, col)
                        if state.ica.is_active &&
                           ch_row > 0 &&
                           !($(state.plot_kwargs[:show_cleaned_ica])) &&
                           $(state.view.show_original_ica)
                            _ica_restored_window(state, ep, win, raw, ch_row)
                        else
                            collect(Float64, raw)
                        end
                    end
                )
            )

            # One Point2 observable per line, updated in place (no per-step allocation of the point array)
            channel_points_obs = _points_observable(
                state.channels.draw_observers,
                state.view.xrange,
                state.data.current,
                epoch_obs,
                state.view.amplitude_scale,
                state.plot_kwargs[:show_cleaned_ica],
                state.view.show_original_ica,
            ) do buf, _
                win = state.view.xrange[]
                ep = epoch_obs[]
                cur = state.data.current[]
                data_series =
                    if state.ica.is_active && ch_row > 0 && !state.plot_kwargs[:show_cleaned_ica][] && state.view.show_original_ica[]
                        _ica_restored_window(state, ep, win, get_data(cur, ep, win, col), ch_row)
                    else
                        get_data(cur, ep, win, col)
                    end
                _window_points!(buf, get_time(cur, ep, win), data_series, state.view.amplitude_scale[], state.view.offset[idx])
                return true
            end

            # Ghost lines data (activations are only computed while the respective ghost is visible;
            # while hidden the buffer is set to NaN once and then left untouched)
            if state.ica.is_active && ch_row > 0
                original_nan = Ref(false)
                original_points_obs = _points_observable(
                    state.channels.draw_observers,
                    state.view.xrange,
                    state.data.current,
                    epoch_obs,
                    state.view.amplitude_scale,
                    state.view.show_original_ica,
                ) do buf, resized
                    win = state.view.xrange[]
                    ep = epoch_obs[]
                    cur = state.data.current[]
                    if state.view.show_original_ica[]
                        restored = _ica_restored_window(state, ep, win, get_data(cur, ep, win, col), ch_row)
                        _window_points!(buf, get_time(cur, ep, win), restored, state.view.amplitude_scale[], state.view.offset[idx])
                        original_nan[] = false
                        return true
                    elseif resized || !original_nan[]
                        _nan_points!(buf, get_time(cur, ep, win))
                        original_nan[] = true
                        return true
                    end
                    return false
                end

                subtracted_nan = Ref(false)
                subtracted_points_obs = _points_observable(
                    state.channels.draw_observers,
                    state.view.xrange,
                    state.data.current,
                    epoch_obs,
                    state.view.amplitude_scale,
                    state.view.show_subtracted_ica,
                ) do buf, resized
                    win = state.view.xrange[]
                    ep = epoch_obs[]
                    cur = state.data.current[]
                    if state.view.show_subtracted_ica[]
                        current_raw = get_data(cur, ep, win, col)
                        restored = _ica_restored_window(state, ep, win, current_raw, ch_row)
                        restored .-= current_raw
                        _window_points!(buf, get_time(cur, ep, win), restored, state.view.amplitude_scale[], state.view.offset[idx])
                        subtracted_nan[] = false
                        return true
                    elseif resized || !subtracted_nan[]
                        _nan_points!(buf, get_time(cur, ep, win))
                        subtracted_nan[] = true
                        return true
                    end
                    return false
                end

                ghost_orig_visible = track(@lift($(state.view.show_original_ica) && $(state.plot_kwargs[:show_cleaned_ica])))
                _create_line!(
                    state.channels.original_lines,
                    col,
                    ax,
                    original_points_obs,
                    :grey,
                    [:grey],
                    1, # linewidth
                    0.5, # clear alpha
                    visible = ghost_orig_visible,
                )

                _create_line!(
                    state.channels.subtracted_lines,
                    col,
                    ax,
                    subtracted_points_obs,
                    :red,
                    [:red],
                    1, # linewidth
                    0.5, # clear alpha
                    visible = state.view.show_subtracted_ica,
                )
            end

            is_repaired = _is_channel_repaired(state, col)

            # Line properties (reuse channel_data_obs for efficiency)
            # Keep line_color as Vector{Bool} observable at all times so GLMakie shader pipeline type is fixed
            line_color = track(
                @lift(
                    ($(state.view.crit_val) > 0.0) ? (abs.($(channel_data_obs)) .>= $(state.view.crit_val)) :
                    fill(false, length($(channel_data_obs)))
                )
            )

            if is_repaired || is_selected
                sel_color = state.plot_kwargs[:selected_channel_color]
                line_colormap = [sel_color, sel_color, sel_color]
                line_width = state.plot_kwargs[:channel_line_width] * 2  # Make repaired/selected channels thicker
            else
                unsel_color = state.plot_kwargs[:unselected_channel_color]
                cleaned_color = get(state.plot_kwargs, :cleaned_channel_color, :green)
                line_colormap = track(
                    @lift(
                        state.ica.is_active && $(state.plot_kwargs[:show_cleaned_ica]) ? [cleaned_color, cleaned_color, :red] :
                        [unsel_color, unsel_color, :red]
                    )
                )
                line_width = state.plot_kwargs[:channel_line_width]
            end

            data_lines_visible =
                track(@lift(!state.ica.is_active || $(state.plot_kwargs[:show_cleaned_ica]) || $(state.view.show_original_ica)))

            # Update or create line
            _create_line!(
                state.channels.data_lines,
                col,
                ax,
                channel_points_obs,
                line_color,
                line_colormap,
                line_width,
                state.plot_kwargs[:channel_line_alpha];
                visible = data_lines_visible,
            )

            # Handle labels
            if !state.view.butterfly[]
                label_y_obs = track(
                    @lift(
                        get_label_y(
                            $(state.data.current),
                            $(epoch_obs),
                            $col,
                            state.view.offset[idx],
                            $(state.view.xrange)[1],
                            $(state.view.amplitude_scale),
                        )
                    )
                )
                _create_label!(state.channels.data_labels, col, ax, time_start_obs, label_y_obs, is_selected)
            else
                _hide_channel_label!(state.channels.data_labels, col)
            end
        else
            _hide_channel_objects!(state.channels, col)
        end
    end

    # Restore axis limits after clearing and redrawing (prevents Makie auto-limits from resetting zoom)
    _restore_axis_limits!(ax, state)
end

"""
    _disconnect_draw_observers!(state)

Disconnect every observable/observer registered by `_draw` from its shared inputs (`xrange`, `current`, ...).
Without this, each redraw (butterfly toggle, ICA apply, hiding a channel, ...) leaves the previous draw's
observers attached, and every later scroll step pays for them again.
"""
function _disconnect_draw_observers!(state)
    for o in state.channels.draw_observers
        _disconnect_observer(o)
    end
    empty!(state.channels.draw_observers)
    return nothing
end

"""Save current axis limits and clear all channel lines/labels."""
function _clear_and_save_limits!(ax, state)
    state.plot_kwargs[:_saved_xlims] = ax.xaxis.attributes.limits[]
    state.plot_kwargs[:_saved_ylims] = ax.yaxis.attributes.limits[]
    _clear_axes!(ax, _all_line_dicts(state))
    _disconnect_draw_observers!(state)
end

"""Restore axis limits from saved values (preserves user zoom across redraws)."""
function _restore_axis_limits!(ax, state)
    if haskey(state.plot_kwargs, :_saved_xlims)
        xlims!(ax, state.plot_kwargs[:_saved_xlims]...)
        ylims!(ax, state.plot_kwargs[:_saved_ylims]...)
        delete!(state.plot_kwargs, :_saved_xlims)
        delete!(state.plot_kwargs, :_saved_ylims)
    end
end

# Data accessor functions
"""Return `(get_data, get_time, get_label_y)` closures for continuous data."""
function _get_data_accessors(data::ContinuousDataState) # data used for dispatch
    get_data = (current, ep, range, col) -> @view current.data[!, col][range]
    get_time = (current, ep, range) -> @view current.data.time[range]
    get_label_y = (current, ep, col, offset, sample_idx, scale) -> current.data[sample_idx, col] * scale + offset
    return get_data, get_time, get_label_y
end

"""Return `(get_data, get_time, get_label_y)` closures for epoched data."""
function _get_data_accessors(state::EpochedDataState)
    get_data = (current, ep, range, col) -> @view current.data[ep][!, col][range]
    get_time = (current, ep, range) -> @view current.data[ep].time[range]
    get_label_y = (current, ep, col, offset, sample_idx, scale) -> current.data[ep][!, col][sample_idx] * scale + offset
    return get_data, get_time, get_label_y
end

"""Plot a channel line (from one `Point2` observable) on the axis and store it."""
function _create_line!(data_lines, col, ax, points_obs, color, colormap, linewidth, alpha; visible = true)
    data_lines[col] = lines!(
        ax,
        points_obs,
        color = color,
        colormap = colormap,
        linewidth = linewidth,
        alpha = alpha,
        visible = visible,
        colorrange = (0, 1),
    )
end

"""Add a channel-name text label on the axis."""
function _create_label!(data_labels, col, ax, x_obs, y_obs, is_selected)
    data_labels[col] =
        text!(ax, x_obs, y_obs, text = String(col), align = (:left, :center), fontsize = 18, color = is_selected ? :red : :black)
end

"""Hide a single channel label."""
function _hide_channel_label!(data_labels, col)
    haskey(data_labels, col) && (data_labels[col].visible = false)
end

"""Hide both line and label for a channel."""
function _hide_channel_objects!(channels, col)
    haskey(channels.data_lines, col) && (channels.data_lines[col].visible = false)
    haskey(channels.data_labels, col) && (channels.data_labels[col].visible = false)
    haskey(channels.original_lines, col) && (channels.original_lines[col].visible = false)
    haskey(channels.subtracted_lines, col) && (channels.subtracted_lines[col].visible = false)
end

# Single function with data access abstraction
"""Draw the extra-channel overlay (line or boolean highlight)."""
function _draw_extra_channel!(ax, state::DataBrowserState{<:AbstractDataState})
    for o in state.extra_channel.draw_observers
        _disconnect_observer(o)
    end
    empty!(state.extra_channel.draw_observers)

    _clear_axes!(ax, [state.extra_channel.data_lines, state.extra_channel.data_labels])

    # Remove old legend if present
    if !isnothing(state.extra_channel.legend)
        delete!(state.extra_channel.legend)
        state.extra_channel.legend = nothing
    end

    if !isempty(state.extra_channel.channels)
        margin = state.plot_kwargs[:channel_offset_margin]
        y_max = state.view.yrange[][end] * margin
        y_min = state.view.yrange[][1] * margin
        n_vis = count(state.channels.visible)
        default_step = n_vis > 1 ? (y_min - y_max) / (n_vis - 1) : -50.0

        step = 0.0
        if !state.view.butterfly[] && length(state.view.offset) > 1
            step = state.view.offset[end] - state.view.offset[end-1]
        end
        if abs(step) < 1e-4
            step = default_step
        end
        base_offset = state.view.butterfly[] ? (y_min + step) : state.view.offset[end]

        # Get data access functions based on type
        get_data, get_time, get_label_y = _get_data_accessors(state.data)

        colors = Makie.wong_colors()

        # Accumulators for the legend
        legend_elements = Any[]
        legend_labels   = String[]

        epoch_obs = state.data isa EpochedDataState ? state.data.current_epoch : Observable(1)

        for (idx, channel) in enumerate(state.extra_channel.channels)
            channel_offset = base_offset + idx * step
            line_color = :black

            if eltype(get_data(state.data.current[], epoch_obs[], 1:1, channel)) == Bool # Boolean data - vspan + legend entry only
                highlight_data = @views _splitgroups(findall(get_data(state.data.current[], epoch_obs[], :, channel)))
                if !isempty(highlight_data[1])
                    # Widen single-sample regions so they're visible (per-region, not global)
                    end_indices = copy(highlight_data[2])
                    n_samples = length(get_time(state.data.current[], epoch_obs[], :))
                    for i in eachindex(end_indices)
                        if end_indices[i] == highlight_data[1][i]
                            end_indices[i] = min(n_samples, end_indices[i] + 5)
                        end
                    end
                    state.extra_channel.data_lines[channel] = vspan!(
                        ax,
                        get_time(state.data.current[], epoch_obs[], highlight_data[1]),
                        get_time(state.data.current[], epoch_obs[], end_indices),
                        color = line_color,
                        alpha = 0.5,
                        visible = true,
                    )
                end
                # Boolean channels identified via the legend only — no in-plot text
                push!(legend_elements, PolyElement(color = (line_color, 0.5), strokecolor = :transparent))
                push!(legend_labels, String(channel))
            else # Regular channel - in-plot label is sufficient, skip legend
                time_obs = @lift(get_time($(state.data.current), $(epoch_obs), :))
                data_obs =
                    @lift(get_data($(state.data.current), $(epoch_obs), :, $channel) .* $(state.view.amplitude_scale) .+ $channel_offset)
                push!(state.extra_channel.draw_observers, time_obs)
                push!(state.extra_channel.draw_observers, data_obs)

                state.extra_channel.data_lines[channel] = lines!(ax, time_obs, data_obs, color = :black, linewidth = 2)

                lbl_x_obs = @lift(get_time($(state.data.current), $(epoch_obs), $(state.view.xrange)[1:1])[1])
                lbl_y_obs = @lift(
                    get_label_y(
                        $(state.data.current),
                        $(epoch_obs),
                        $channel,
                        $channel_offset,
                        $(state.view.xrange)[1],
                        $(state.view.amplitude_scale),
                    )
                )
                push!(state.extra_channel.draw_observers, lbl_x_obs)
                push!(state.extra_channel.draw_observers, lbl_y_obs)

                state.extra_channel.data_labels[channel] =
                    text!(ax, lbl_x_obs, lbl_y_obs, text = String(channel), align = (:left, :center), color = :black, fontsize = 18)
            end
        end

        # Build / rebuild the legend in the figure if we have a figure reference
        if !isnothing(state.extra_channel.fig) && !isempty(legend_labels)
            state.extra_channel.legend = Legend(
                state.extra_channel.fig[1, 1],
                legend_elements,
                legend_labels;
                tellwidth = false,
                tellheight = false,
                halign = :right,
                valign = :bottom,
                framevisible = true,
                padding = (8.0f0, 8.0f0, 8.0f0, 8.0f0),
                margin = (4.0f0, 4.0f0, 4.0f0, 4.0f0),
            )
        end
    end
end

"""Return a descriptive window title based on the data type."""
_get_title(dat::EpochData) = "Epoch 1/$(n_epochs(dat))"
_get_title(dat::ContinuousData) = ""
_get_title(dat::ErpData) = "Epoch Average (n=$(n_epochs(dat)))"

"""Toggle visibility of a vertical marker."""
function _plot_vertical_lines!(ax, marker, active)
    marker.line.visible = active
    marker.text.visible = active
    marker.visible = active
end

"""
    _add_scale_indicator!(ax, state, plot_kwargs)

Add a scale indicator bar to the plot showing the amplitude scale.
"""
function _add_scale_indicator!(ax, state, plot_kwargs)

    scale_value = plot_kwargs[:scale_indicator_value]

    # Use an observable for position so it can be updated by dragging
    if !(plot_kwargs[:scale_indicator_position] isa Observable)
        plot_kwargs[:scale_indicator_position] = Observable{Tuple{Float64,Float64}}(Float64.(plot_kwargs[:scale_indicator_position]))
    end
    pos_obs = plot_kwargs[:scale_indicator_position]

    # Get axis limits observables
    xlims_obs = ax.xaxis.attributes.limits
    ylims_obs = ax.yaxis.attributes.limits

    # Calculate position in data coordinates
    x_pos = @lift($xlims_obs[1] + ($xlims_obs[2] - $xlims_obs[1]) * $pos_obs[1])
    y_top = @lift($ylims_obs[1] + ($ylims_obs[2] - $ylims_obs[1]) * $pos_obs[2])
    y_bottom = @lift($y_top - scale_value * $(state.view.amplitude_scale))

    tick_length = @lift(($xlims_obs[2] - $xlims_obs[1]) * 0.005)
    tick_left = @lift($x_pos - $tick_length)
    tick_right = @lift($x_pos + $tick_length)

    # Use linesegments! to draw all parts of the scale indicator in a single draw call
    segments = @lift([
        Point2f($x_pos, $y_bottom),
        Point2f($x_pos, $y_top),           # Vertical line
        Point2f($tick_left, $y_bottom),
        Point2f($tick_right, $y_bottom), # Bottom tick
        Point2f($tick_left, $y_top),
        Point2f($tick_right, $y_top),        # Top tick
    ])
    linesegments!(ax, segments, color = :black, linewidth = 1.5, overdraw = true)

    # Add label (formatted as integer if whole number, e.g. "100 μV")
    val_str = isinteger(scale_value) ? "$(Int(scale_value)) μV" : "$scale_value μV"
    label_x = @lift($x_pos + ($xlims_obs[2] - $xlims_obs[1]) * 0.01)
    text!(
        ax,
        @lift(Point2f($label_x, $y_top)),
        text = val_str,
        align = (:left, :center),
        fontsize = 14,
        color = :black,
        space = :data,
        overdraw = true,
    )

    # Enable mouse dragging
    is_dragging = Observable(false)

    on(events(ax).mousebutton, priority = 10) do event
        if event.button == Mouse.left
            if event.action == Mouse.press
                mpos = mouseposition(ax)
                # Check if mouse is near the scale indicator (hit box)
                x_tol = (xlims_obs[][2] - xlims_obs[][1]) * 0.05
                y_tol = (ylims_obs[][2] - ylims_obs[][1]) * 0.05
                if abs(mpos[1] - x_pos[]) < x_tol && mpos[2] >= (y_bottom[] - y_tol) && mpos[2] <= (y_top[] + y_tol)
                    is_dragging[] = true
                    return Consume(true)
                end
            elseif event.action == Mouse.release
                if is_dragging[]
                    is_dragging[] = false
                    return Consume(true)
                end
            end
        end
        return Consume(false)
    end

    on(events(ax).mouseposition, priority = 10) do _
        if is_dragging[]
            mpos = mouseposition(ax)
            new_rel_x = (mpos[1] - xlims_obs[][1]) / (xlims_obs[][2] - xlims_obs[][1])
            new_rel_y = (mpos[2] - ylims_obs[][1]) / (ylims_obs[][2] - ylims_obs[][1])
            pos_obs[] = (clamp(new_rel_x, 0.0, 1.0), clamp(new_rel_y, 0.0, 1.0))
            return Consume(true)
        end
        return Consume(false)
    end
end


"""Launch the interactive data browser for a single EEG dataset."""
function plot_databrowser(dat::EegData, ica = nothing; screen = nothing, kwargs...)

    # Check if CairoMakie is being used and warn about lack of interactivity
    if string(Makie.current_backend()) == "CairoMakie"
        @minimal_warning "CairoMakie detected. For full interactivity in plot_databrowser, use GLMakie."
    end

    # Generate window title from dataset
    title_str = _generate_window_title(dat)
    _set_window_title(title_str)

    try
        # Merge user kwargs with defaults
        plot_kwargs = _merge_plot_kwargs(PLOT_DATABROWSER_KWARGS, kwargs)

        # Common fig/ax/state/ui setup
        fontsize_kw = isnothing(plot_kwargs[:theme_fontsize]) ? (;) : (; fontsize = plot_kwargs[:theme_fontsize])
        fig = _eegfun_figure(; figure_padding = plot_kwargs[:figure_padding], fontsize_kw...)
        ax = Axis(fig[1, 1], xlabel = plot_kwargs[:xlabel], ylabel = plot_kwargs[:ylabel], title = _get_title(dat))

        state = _create_browser_state(dat, dat.layout.data.label, ax, ica, plot_kwargs)
        state.extra_channel.fig = fig  # store figure reference for legend placement
        _setup_ui(fig, ax, state, dat, ica, plot_kwargs)

        _draw(ax, state)
        _draw_extra_channel!(ax, state)

        # Add scale indicator
        if plot_kwargs[:show_scale_indicator]
            _add_scale_indicator!(ax, state, plot_kwargs)
        end

        # Display on the provided screen if given, otherwise use default display
        screen_obj = nothing
        if plot_kwargs[:display_plot]
            if !isnothing(screen)
                display(screen, fig)
                screen_obj = screen
            else
                screen_obj = _display_figure(fig)
            end
        end

        return (fig = fig, ax = ax, analysis_settings = state.analysis_settings, screen = screen_obj)
    finally
        _set_window_title("Makie")
    end
end

"""
    plot_databrowser(; launcher::Bool = false, screen = nothing, kwargs...)

Launch the interactive continuous EEG data browser.

When called with no arguments, automatically loads the bundled example raw EEG dataset
(`example1.bdf`) along with its full 72-channel layout (`biosemi72.csv`) for an immediate demonstration.

Pass `launcher = true` to open the drag-and-drop file selector window instead.

# Examples
```julia
# Quick demo with bundled raw dataset and layout:
plot_databrowser()

# Open file drag-and-drop launcher:
plot_databrowser(launcher = true)
```
"""
function plot_databrowser(; launcher::Bool = false, screen = nothing, kwargs...)
    if launcher
        return _plot_databrowser_launcher(; screen = screen, kwargs...)
    end
    raw_path = example_path("data/bdf/example1.bdf")
    layout_path = example_path("layouts/biosemi/biosemi72.csv")
    dat = read_raw_data(raw_path)
    lay = read_layout(layout_path)
    polar_to_cartesian_xy!(lay)
    dat = create_eegfun_data(dat, lay)
    dat.file = basename(raw_path)
    return plot_databrowser(dat; screen = screen, kwargs...)
end

"""
    _plot_databrowser_launcher(; screen = nothing, kwargs...)

Open a file-selector window that lets you drag-and-drop EEG files or click to load the bundled example data.
"""
function _plot_databrowser_launcher(; screen = nothing, kwargs...)
    shared_layout = Ref{Union{Nothing,Layout}}(nothing)
    layout_status = Observable{String}("Layout: NA")

    fig = _eegfun_figure(size = (420, 260))
    Label(fig[1, 1], "EegFun Data Browser", fontsize = 22, halign = :center)
    Label(fig[2, 1], "Drop EEG file(s) here", fontsize = 14, halign = :center, color = :gray50)
    status_label = Label(fig[3, 1], "", fontsize = 12, halign = :center, color = :gray40)
    btn_layout = Button(fig[4, 1], label = "Select Layout File", width = 180, height = 30, fontsize = 14)
    Label(fig[5, 1], layout_status, fontsize = 12, color = :gray30)
    btn_example = Button(fig[6, 1], label = "Load Example (example1.bdf)", width = 230, height = 30, fontsize = 14)
    n_loaded = Ref(0)

    function _load_layout(filepath::AbstractString)
        try
            lay = read_layout(filepath)
            polar_to_cartesian_xy!(lay)
            shared_layout[] = lay
            layout_status[] = "Layout: $(basename(filepath))"
            status_label.text = "Loaded layout successfully"
        catch e
            @error "Failed to load layout" exception = (e, catch_backtrace())
            status_label.text = "Layout Error: $(sprint(showerror, e))"
        end
    end

    function _load_file(filepath::AbstractString)
        filepath = _sanitize_dropped_path(filepath)
        if !isfile(filepath)
            status_label.text = "File not found: $(basename(filepath))"
            return
        end

        if lowercase(splitext(filepath)[2]) == ".csv"
            _load_layout(filepath)
            return
        end

        status_label.text = "Loading $(basename(filepath))…"
        try
            dat_eeg, ica_data = _read_any_eeg_file(filepath, shared_layout[])
            dat_eeg.file = basename(filepath)
            _set_window_title(_generate_window_title(dat_eeg))
            plot_databrowser(dat_eeg, ica_data; screen = _create_screen(), kwargs...)
            n_loaded[] += 1
            status_label.text = n_loaded[] == 1 ? "Opened: $(basename(filepath))" : "Opened $(n_loaded[]) files"
        catch e
            @error "Failed to load file" exception = (e, catch_backtrace())
            status_label.text = "Error: $(sprint(showerror, e))"
        end
    end

    on(btn_layout.clicks) do _
        filepath = NativeFileDialog.pick_file(pwd())
        if filepath != "" && isfile(filepath)
            _load_layout(filepath)
        end
    end

    on(btn_example.clicks) do _
        status_label.text = "Loading example1.bdf…"
        try
            raw_path = example_path("data/bdf/example1.bdf")
            layout_path = example_path("layouts/biosemi/biosemi72.csv")
            dat = read_raw_data(raw_path)
            lay = read_layout(layout_path)
            polar_to_cartesian_xy!(lay)
            dat = create_eegfun_data(dat, lay)
            dat.file = basename(raw_path)
            _set_window_title(_generate_window_title(dat))
            plot_databrowser(dat; screen = _create_screen(), kwargs...)
            n_loaded[] += 1
            status_label.text = "Opened example1.bdf"
        catch e
            @error "Failed to load example data" exception = (e, catch_backtrace())
            status_label.text = "Error: $(sprint(showerror, e))"
        end
    end

    on(events(fig).dropped_files) do files
        for f in files
            _load_file(f)
        end
    end

    screen_obj = _display_in_screen(fig)
    _set_window_title("EegFun Data Browser")
    return (fig = fig, screen = screen_obj)
end

"""Decode percent-encoded bytes into a UTF-8 string."""
function _percent_decode(s::AbstractString)
    bytes = UInt8[]
    i = 1
    n = ncodeunits(s)
    while i <= n
        c = s[i]
        if c == '%' && i + 2 <= n && isxdigit(s[i+1]) && isxdigit(s[i+2])
            b = parse(UInt8, s[(i+1):(i+2)], base = 16)
            push!(bytes, b)
            i += 3
        else
            for b in codeunits(string(c))
                push!(bytes, b)
            end
            i = nextind(s, i)
        end
    end
    return String(bytes)
end

"""Sanitize a filepath from drag-and-drop (strip URI scheme, URL-decode)."""
function _sanitize_dropped_path(filepath::AbstractString)
    s = strip(filepath, ['\r', '\n', '"', '\''])
    if startswith(s, "file://")
        s = replace(s, r"^file://localhost" => "")
        s = replace(s, r"^file://" => "")
        s = _percent_decode(s)
        if Sys.iswindows() && occursin(r"^/[A-Za-z]:", s)
            s = s[2:end]
        end
    end
    return String(s)
end

"""Load any supported EEG file and return `(dat, ica)` — ica may be `nothing`."""
function _read_any_eeg_file(filepath::String, layout::Union{Nothing,Layout} = nothing)
    ext = lowercase(splitext(filepath)[2])
    ica_data = nothing

    dat_eeg = if ext == ".jld2"
        read_data(filepath)
    elseif ext == ".set"
        res = read_eeglab(filepath)
        if res isa Tuple
            dat_eeg_tmp, ica_data = res
            dat_eeg_tmp
        else
            res
        end
    else
        raw = read_raw_data(filepath)
        if isnothing(layout)
            create_eegfun_data(raw)
        else
            create_eegfun_data(raw, layout)
        end
    end

    isnothing(dat_eeg) && @minimal_error("Unsupported file format or failed to read: $filepath")
    return (dat_eeg, ica_data)
end

"""Load data from a `.jld2` file or pattern and open the data browser."""
function plot_databrowser(
    filename::String,
    ica = nothing;
    input_dir::String = pwd(),
    participant_selection = participants(),
    screen = nothing,
    kwargs...,
)
    if endswith(filename, ".jld2")
        dat_eeg = read_data(filename)
        if !isnothing(ica)
            ica = read_data(ica)
        end
        return plot_databrowser(dat_eeg, ica; screen = screen, kwargs...)
    else
        files = _find_batch_files(filename, input_dir, participant_selection)
        isempty(files) && @minimal_error "No files matching pattern '$filename' in $input_dir"

        results = NamedTuple[]
        for file in sort(files, by = _natural_sort_key)
            file_path = joinpath(input_dir, file)
            @info "Browsing: $file"
            dat_eeg = read_data(file_path)
            isnothing(dat_eeg) && continue
            result = plot_databrowser(dat_eeg, nothing; screen = screen, kwargs...)
            push!(results, result)
        end
        return results
    end
end

"""Open a separate data browser window for each dataset in the vector."""
function plot_databrowser(data::Vector{<:EegData}, ica = nothing; screen = nothing, kwargs...)
    @info "Vector of $(length(data)) datasets provided — opening a browser for each"
    results = NamedTuple[]
    display_plot = get(kwargs, :display_plot, true)
    for dat in data
        _set_window_title(_generate_window_title(dat))  # must be before Screen()
        s = if !isnothing(screen)
            screen
        elseif display_plot
            _create_screen()
        else
            nothing
        end
        res = plot_databrowser(dat, ica; screen = s, kwargs...)
        push!(results, res)
    end
    return results
end
