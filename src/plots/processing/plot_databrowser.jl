const PLOT_DATABROWSER_KWARGS = Dict{Symbol,Tuple{Any,String}}(

    # Figure and layout
    :figure_padding => ((50, 0, 50, 50), "Figure padding as (left, right, bottom, top)"),

    # Axis styling
    :xlabel => ("Time (S)", "X-axis label"),
    :ylabel => ("Amplitude (μV)", "Y-axis label"),

    # UI styling
    :theme_fontsize => (18, "Font size for theme and UI elements"),

    # Line styling
    :channel_line_width => (1, "Line width for channel lines"),
    :channel_line_alpha => (1, "Transparency for channel lines"),
    :selected_channel_color => (:black, "Color for selected channels"),
    :unselected_channel_color => (:darkgrey, "Color for unselected channels"),
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
        poly_element = poly!(ax, [Point2f(0.0, 0.0)], color = (plot_kwargs[:selection_color], plot_kwargs[:selection_alpha]), visible = false)
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
        Observable(false), Observable(false),
        Observable(plot_kwargs[:default_hp_freq]), Observable(plot_kwargs[:default_lp_freq]),
        Observable("iir"), Observable("iir"),
        Observable(2), Observable(2),
        Observable("filtfilt"), Observable("filtfilt")
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
    individually_selected::Vector{Symbol}
    data_labels::Dict{Symbol,Makie.Text}
    data_lines::Dict{Symbol,Union{Makie.Lines,Makie.PolyElement,Any}}
    original_lines::Dict{Symbol,Union{Makie.Lines,Makie.PolyElement,Any}}    # ICA ghost: pre-removal signal
    subtracted_lines::Dict{Symbol,Union{Makie.Lines,Makie.PolyElement,Any}}  # ICA ghost: removed artefact
    function ChannelState(channel_labels::Vector{Symbol})
        new(
            channel_labels,
            fill(true, length(channel_labels)),
            fill(false, length(channel_labels)),
            Symbol[],
            Dict{Symbol,Makie.Text}(),
            Dict{Symbol,Union{Makie.Lines,Makie.PolyElement,Any}}(),
            Dict{Symbol,Union{Makie.Lines,Makie.PolyElement,Any}}(),
            Dict{Symbol,Union{Makie.Lines,Makie.PolyElement,Any}}(),
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

"""Track optional extra channel overlays (e.g. EOG, EMG)."""
mutable struct ExtraChannelInfo
    channels::Vector{Symbol}
    data_lines::Dict{Symbol,Union{Makie.Lines,Makie.PolyElement,Any}}
    data_labels::Dict{Symbol,Makie.Text}
    fig::Union{Nothing,Figure}          # owning figure for legend placement
    legend::Union{Nothing,Any}          # current Legend object (replaced on redraw)
    ExtraChannelInfo() = new(Symbol[], Dict(), Dict(), nothing, nothing)
end

"""Visual state tracker for dynamic Independent Component Analysis (ICA) artifact removal."""
mutable struct IcaVisualState
    original::Union{Nothing,InfoIca}
    removed_components::Vector{Int}
    is_active::Bool

    function IcaVisualState(ica::Union{Nothing,InfoIca})
        new(ica, Int[], !isnothing(ica))
    end
end

"""Top-level composite state for the interactive data browser."""
mutable struct DataBrowserState{T<:AbstractDataState}
    view::ViewState
    channels::ChannelState
    data::T
    selection::SelectionState
    markers::Vector{Marker}
    ica::IcaVisualState
    extra_channel::ExtraChannelInfo
    reference_state::Symbol
    channel_repair_history::Vector{Tuple{Vector{Symbol},Symbol,Matrix{Float64}}}  # (channels, method, original_data) - stack for multiple undos
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
    selected_regions = state.selection.selected_regions[]

    # Update the observable
    state.analysis_settings[] = AnalysisSettings(
        hp_freq,
        lp_freq,
        state.reference_state,
        repaired_channels,
        repair_method,
        selected_regions,
        state.ica.removed_components,
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
_get_initial_window_size(dat::ContinuousData) = min(5000, nrow(dat.data))
_get_initial_window_size(dat::ErpData) = nrow(dat.data)
_get_initial_window_size(dat::EpochData) = nrow(dat.data[1])

"""Map an `EegData` subtype to its corresponding `AbstractDataState` type."""
_data_state_type(::Type{ContinuousData}) = ContinuousDataState
_data_state_type(::Type{EpochData}) = EpochedDataState
_data_state_type(::Type{ErpData}) = ContinuousDataState

"""Return the current DataFrame for the active data state."""
_get_current_data(state::ContinuousDataState) = state.current[].data
_get_current_data(state::EpochedDataState) = state.current[].data[state.current_epoch[]]
"""Return `(t_start, t_end)` of the current data."""
_get_time_bounds(dat::ContinuousDataState) = (dat.current[].data.time[1], dat.current[].data.time[end])
_get_time_bounds(dat::EpochedDataState) =
    (dat.current[].data[dat.current_epoch[]].time[1], dat.current[].data[dat.current_epoch[]].time[end])
"""Check whether a column exists in the current data."""
_has_column(state::ContinuousDataState, col::Symbol) = col in propertynames(state.current[].data)
_has_column(state::EpochedDataState, col::Symbol) = col in propertynames(state.current[].data[state.current_epoch[]])
"""Trigger the current-data observable so listeners redraw."""
_notify_data_update(state::AbstractDataState) = notify(state.current)
"""Replace current data with a fresh copy of the original."""
_reset_to_original!(state::AbstractDataState) = state.current[] = copy(state.original)

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
    update_theme!(Theme(fontsize = plot_kwargs[:theme_fontsize]))
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
            marker_index = findfirst(m -> m.name == marker_symbol, state.markers)
            if !isnothing(marker_index)
                push!(configs, ToggleConfig(toggle_label, (active) -> _plot_vertical_lines!(ax, state.markers[marker_index], active)))
            end
        end
    end

    # Add filter toggles if available
    if state.data.original.analysis_info.hp_filter == 0.0
        push!(configs, ToggleConfig("HP-Filter", (_) -> _apply_hp_filter!(state)))
    end
    if state.data.original.analysis_info.lp_filter == 0.0
        push!(configs, ToggleConfig("LP-Filter", (_) -> _apply_lp_filter!(state)))
    end

    # Create toggles
    toggles_grid = []
    for config in configs
        tog = Toggle(fig)
        on(tog.active) do active
            config.action(active)
        end
        
        if config.label == "Trigger"
            btn = Button(fig, label="Trigger", halign=:left)
            on(btn.clicks) do _
                _show_trigger_menu(state, ax, :trigger)
            end
            push!(toggles_grid, [tog, btn])
        elseif config.label == "HP-Filter" || config.label == "LP-Filter"
            btn = Button(fig, label=config.label, halign=:left)
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

"""Popup: select which triggers to display."""
function _show_trigger_menu(state, ax, marker_symbol)
    trigger_col = if state.data isa ContinuousDataState
        state.data.original.data[!, marker_symbol]
    else
        vcat([df[!, marker_symbol] for df in state.data.original.data]...)
    end
    unique_triggers = sort(unique(filter(!iszero, trigger_col)))

    if !haskey(state.plot_kwargs, :active_triggers)
        state.plot_kwargs[:active_triggers] = Set(unique_triggers)
    end
    active_set = state.plot_kwargs[:active_triggers]

    menu_fig = Figure(size = (600, 600))
    Label(menu_fig[1, 1], "Select Triggers to Display", fontsize=18, halign=:center)

    scroll_area = menu_fig[2, 1] = GridLayout()
    cbs = _build_checkbox_grid!(scroll_area, string.(unique_triggers), 5, (trig, _) -> (parse(Float64, trig) in active_set))

    groups = [
        ("All", ch -> true),
        ("None", ch -> false),
    ]
    _add_group_buttons!(menu_fig, 3, cbs, string.(unique_triggers), groups)

    action_area = menu_fig[4, 1] = GridLayout()
    btn_apply = Button(action_area[1, 1], label="Apply", width=200)
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

    display(GLMakie.Screen(), menu_fig)
end
"""Popup: select filter cutoff frequency for a specific filter."""
function _show_single_filter_menu(state, dat, filter_type::Symbol)
    menu_fig = Figure(size = (450, 350))
    
    title_str = filter_type == :hp ? "Highpass Filter Settings" : "Lowpass Filter Settings"
    Label(menu_fig[1, 1], title_str, fontsize = 18, halign = :center, font=:bold)

    grid = menu_fig[2, 1] = GridLayout(valign=:top)
    row = 1
    
    if filter_type == :hp
        range_freq = 0.1:0.1:2
        freq_field = :hp_freq; method_field = :hp_method; order_field = :hp_order; func_field = :hp_func
        default_val = 0.5
    else
        range_freq = 5:5:60
        freq_field = :lp_freq; method_field = :lp_method; order_field = :lp_order; func_field = :lp_func
        default_val = 20.0
    end
    
    fs = state.data.filter_state
    if getfield(fs, freq_field)[] == 0.0
        getfield(fs, freq_field)[] = default_val
    end

    # Cutoff
    Label(grid[row, 1], "Cutoff (Hz):", halign = :left)
    slider_freq = Slider(grid[row, 2], range = range_freq, startvalue = getfield(fs, freq_field)[])
    Label(grid[row, 3], @lift(string(round($(slider_freq.value), digits=1))), halign = :left)
    on(slider_freq.value) do val; getfield(fs, freq_field)[] = val; end
    row += 1

    # Order
    Label(grid[row, 1], "Order:", halign = :left)
    slider_order = Slider(grid[row, 2], range = 1:10, startvalue = getfield(fs, order_field)[])
    Label(grid[row, 3], @lift(string($(slider_order.value))), halign = :left)
    on(slider_order.value) do val; getfield(fs, order_field)[] = val; end
    row += 1

    # Method
    Label(grid[row, 1], "Method:", halign = :left)
    method_labels = ["Butterworth (IIR)", "FIR (Hamming)"]
    method_values = ["iir", "fir"]
    current_method = getfield(fs, method_field)[]
    menu_method = Menu(grid[row, 2], options = zip(method_labels, method_values), default = current_method == "iir" ? "Butterworth (IIR)" : "FIR (Hamming)", width = 220)
    on(menu_method.selection) do val; getfield(fs, method_field)[] = val; end
    row += 1

    # Function
    Label(grid[row, 1], "Function:", halign = :left)
    func_labels = ["filtfilt (zero-phase)", "filt (one-pass)"]
    func_values = ["filtfilt", "filt"]
    current_func = getfield(fs, func_field)[]
    menu_func = Menu(grid[row, 2], options = zip(func_labels, func_values), default = current_func == "filtfilt" ? "filtfilt (zero-phase)" : "filt (one-pass)", width = 220)
    on(menu_func.selection) do val; getfield(fs, func_field)[] = val; end
    row += 1

    rowgap!(grid, 10)
    colgap!(grid, 15)

    action_area = menu_fig[3, 1] = GridLayout()
    btn_apply = Button(action_area[1, 1], label = "Apply", width = 150)
    on(btn_apply.clicks) do _
        _apply_filters!(state)
    end

    display(GLMakie.Screen(), menu_fig)
end


"""Popup: select which channels are visible."""
function _show_labels_menu(state, ax)
    all_channels = state.channels.labels
    menu_fig = Figure(size = (900, 700))
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
        state.channels.individually_selected = Symbol[]
        for (i, cb) in enumerate(cbs)
            state.channels.visible[i] = cb.checked[]
        end
        _clear_and_save_limits!(ax, state)
        _update_channel_offsets!(state)
        _draw(ax, state)
    end

    display(GLMakie.Screen(), menu_fig)
end

"""Create the channel-label selection button for the control panel."""
_create_labels_menu(fig, ax, state) = _popup_button(fig, "Labels", () -> _show_labels_menu(state, ax))

"""Popup: select reference channels and apply re-referencing."""
function _show_reference_menu(state, dat)
    all_channels = state.channels.labels
    menu_fig = Figure(size = (900, 700))
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
        if opt != :none
            _rereference!(state.data, opt)
            _notify_data_update(state.data)
            _update_analysis_settings!(state)
        end
    end

    display(GLMakie.Screen(), menu_fig)
end

"""Create the re-reference button for the control panel."""
_create_reference_menu(fig, state, dat) = _popup_button(fig, "Reference", () -> _show_reference_menu(state, dat))

"""Toggle an ICA component for removal and apply subtraction algebra to the dataset."""
function _toggle_ica_component!(state, component_id::Union{Nothing,Int})
    if isnothing(component_id)
        # Fast path reset
        empty!(state.ica.removed_components)
        _reset_to_original!(state.data)
        _reapply_active_filters!(state)
        return
    end

    # Determine direction of toggle
    is_removing = !(component_id in state.ica.removed_components)

    if is_removing
        push!(state.ica.removed_components, component_id)
    else
        filter!(x -> x != component_id, state.ica.removed_components)
    end

    _apply_incremental_ica_update!(state.data, state.ica.original, component_id, is_removing)
end

"""Re-apply all currently removed ICA components to the dataset."""
function _reapply_all_ica_removals!(state)
    if state.ica.is_active && !isempty(state.ica.removed_components)
        for comp in state.ica.removed_components
            _apply_incremental_ica_update!(state.data, state.ica.original, comp, true)
        end
    end
end

"""Re-apply active filters without resetting data (used after ICA reset to preserve filter state)."""
function _reapply_active_filters!(state)
    fs = state.data.filter_state
    fs.hp_active[] && _apply_filter!(state, :hp, fs.hp_freq[], fs.hp_method[], fs.hp_order[], fs.hp_func[])
    fs.lp_active[] && _apply_filter!(state, :lp, fs.lp_freq[], fs.lp_method[], fs.lp_order[], fs.lp_func[])
end

"""Apply the isolated artifact of a single ICA component linearly to a Continuous dataset."""
function _apply_incremental_ica_update!(state::ContinuousDataState, ica::InfoIca, comp::Int, subtract::Bool)
    dat = state.current[].data
    cols = ica.layout.data.label

    unmix_vec = ica.unmixing[comp, :]
    mix_vec = ica.mixing[:, comp]
    n_samples = nrow(dat)

    # Calculate 1D activation series for the component
    activation = zeros(Float64, n_samples)
    for (i, col_sym) in enumerate(cols)
        norm_ch = (state.original.data[!, col_sym] .- ica.mean[i]) ./ ica.scale
        activation .+= unmix_vec[i] .* norm_ch
    end

    # Dynamically inject the isolated footprint back into current sensor space
    for (i, col_sym) in enumerate(cols)
        artifact_ch = mix_vec[i] .* ica.scale .* activation
        if subtract
            dat[!, col_sym] .-= artifact_ch
        else
            dat[!, col_sym] .+= artifact_ch
        end
    end
end

"""Apply the isolated artifact of a single ICA component linearly to an Epoched dataset."""
function _apply_incremental_ica_update!(state::EpochedDataState, ica::InfoIca, comp::Int, subtract::Bool)
    cols = ica.layout.data.label
    unmix_vec = ica.unmixing[comp, :]
    mix_vec = ica.mixing[:, comp]

    for ep = 1:length(state.current[].data)
        dat_ep = state.current[].data[ep]
        orig_ep = state.original.data[ep]
        n_samples = nrow(dat_ep)

        activation = zeros(Float64, n_samples)
        for (i, col_sym) in enumerate(cols)
            norm_ch = (orig_ep[!, col_sym] .- ica.mean[i]) ./ ica.scale
            activation .+= unmix_vec[i] .* norm_ch
        end

        for (i, col_sym) in enumerate(cols)
            artifact_ch = mix_vec[i] .* ica.scale .* activation
            if subtract
                dat_ep[!, col_sym] .-= artifact_ch
            else
                dat_ep[!, col_sym] .+= artifact_ch
            end
        end
    end
end

"""Popup: toggle ICA component removals and display options."""
function _show_ica_menu(state, ax, ica)
    all_components = ica.ica_label
    menu_fig = Figure(size = (800, 600))
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
        _toggle_ica_component!(state, nothing)  # reset
        for (i, cb) in enumerate(cbs)
            if cb.checked[]
                n = _extract_int(String(all_components[i]))
                !isnothing(n) && _toggle_ica_component!(state, n)
            end
        end
        _notify_data_update(state.data)
        _draw(ax, state)
        _update_analysis_settings!(state)
    end

    # Display options row
    display_area = menu_fig[5, 1] = GridLayout()
    Label(display_area[1, 1:6], "Display Options", fontsize = 16, halign = :center)
    display_items = [
        (state.view.show_original_ica, "Data"),
        (state.plot_kwargs[:show_cleaned_ica], "Data - ICA"),
        (state.view.show_subtracted_ica, "ICA Activation"),
    ]
    for (j, (obs, lbl)) in enumerate(display_items)
        col = (j - 1) * 2 + 1
        cb = Checkbox(display_area[2, col], checked = obs[], width = 16, height = 16)
        Label(display_area[2, col+1], lbl, halign = :left)
        on(cb.checked) do active
            obs[] = active
        end
        colgap!(display_area, col, 5)
        j < length(display_items) && colgap!(display_area, col + 1, 30)
    end

    display(GLMakie.Screen(), menu_fig)
end

"""Create the ICA components button for the control panel."""
function _create_ica_menu(fig, ax, state, ica)
    if !haskey(state.plot_kwargs, :show_cleaned_ica)
        state.plot_kwargs[:show_cleaned_ica] = Observable(false)
        state.view.show_original_ica[] = true
    end
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
            _clear_and_save_limits!(ax, state)
            state.data.current_epoch[] = epoch_num
            ax.title = "Epoch $(epoch_num)/$(n_epochs(state.data.original))"
            _update_markers!(ax, state)
            _draw(ax, state)
            _draw_extra_channel!(ax, state)
        end
    end

    # Make slider observe current_epoch changes (for left/right key navigation)
    on(state.data.current_epoch) do epoch_num
        updating_from_keyboard[] = true
        slider_epoch.value[] = epoch_num
        updating_from_keyboard[] = false
    end

    return hcat(slider_epoch, label)
end

"""Open a popup menu with Topoplot / Spectrum / Region actions."""
function _show_additional_menu(state, clicked_region_idx = nothing)

    # Create the menu figure
    # TODO: why does new window not always take this size?
    menu_fig = Figure(size = (300, 300))
    plot_types = ["Topoplot", "Spectrum", "Get Selected Regions"]

    menu_buttons = [Button(menu_fig[idx, 1], label = plot_type) for (idx, plot_type) in enumerate(plot_types)]

    for btn in menu_buttons
        on(btn.clicks) do _
            if btn.label[] == "Get Selected Regions"
                # Get the boolean vector of selected regions
                selected_regions_bool = _get_selected_regions_bool(state)
                @info "Selected regions boolean vector:"
                @info "Length: $(length(selected_regions_bool))"
                @info "Number of selected samples: $(sum(selected_regions_bool))"
                @info "Selected regions: $(state.selection.selected_regions[])"
            else
                selected_data = _subset_selected_data(state, clicked_region_idx)
                isnothing(selected_data) && return # No data available, just return
                if btn.label[] == "Topoplot"
                    plot_topography(selected_data)
                elseif btn.label[] == "Spectrum"
                    plot_channel_spectrum(selected_data)
                end
            end
        end
    end

    new_screen = GLMakie.Screen(size = (300, 300))
    display(new_screen, menu_fig)
end

"""Open the interactive channel repair window with checkboxes and method selection."""
function _channel_repair_menu(state, selected_channels, ax)
    all_channels = state.channels.labels

    menu_fig = Figure(size = (900, 700))
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

    new_screen = GLMakie.Screen()
    display(new_screen, menu_fig)
end

"""Interpolate the selected channels, store originals in repair history, and redraw."""
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

    # Store original data before repair
    original_data = copy(_get_channel_data_matrix(state.data.current[], selected_channels))

    # Perform the repair
    repair_channels!(state.data.current[], selected_channels, method = method)

    # Store repair in history
    push!(state.channel_repair_history, (selected_channels, method, original_data))

    # Notify that data has been updated
    _notify_data_update(state.data)

    # Update analysis settings
    _update_analysis_settings!(state)

    # Clear and redraw the plot
    _clear_and_save_limits!(ax, state)
    _draw(ax, state)

    total_repairs = length(state.channel_repair_history)
    @info "Successfully repaired channels: $(join(string.(selected_channels), ", ")) using $method"
    @info "Total repairs in history: $total_repairs"
end

"""Pop the last channel repair from history and restore original data."""
function _undo_last_repair!(state, ax)

    if isempty(state.channel_repair_history)
        @info "No repairs to undo"
        return
    end

    # Get the last repair
    channels, method, original_data = pop!(state.channel_repair_history)

    # Restore original data
    _restore_channel_data!(state.data.current[], channels, original_data)

    # Notify that data has been updated
    _notify_data_update(state.data)

    # Clear and redraw the plot
    _clear_and_save_limits!(ax, state)
    _draw(ax, state)

    remaining_repairs = length(state.channel_repair_history)
    @info "Undid repair of channels: $(join(string.(channels), ", ")) (was $method)"
    @info "Remaining repairs in history: $remaining_repairs"
end


"""Extract a matrix of selected channel columns from EEG data."""
_get_channel_data_matrix(data::EegData, channels) = Matrix(data.data[:, channels])

"""Write previously saved channel data back into the dataset."""
_restore_channel_data!(data::EegData, channels, original_data) = data.data[:, channels] = original_data


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
    sliders = _create_extreme_slider(fig, state)

    # Add navigation sliders specific to continuous data
    slider_range = Slider(fig[3, 1], range = 100:50:30000, startvalue = state.view.xrange[][end], snap = true)
    slider_x = Slider(fig[2, 1], range = 1:50:nrow(state.data.current[].data), startvalue = 1, snap = true)

    on(slider_range.value) do x
        new_range = slider_x.value.val:min(nrow(state.data.current[].data), x + slider_x.value.val)
        if length(new_range) > 1
            state.view.xrange[] = new_range
        end
    end

    on(slider_x.value) do x
        new_range = x:min(nrow(state.data.current[].data), (x + slider_range.value.val) - 1)
        if length(new_range) > 1
            state.view.xrange[] = new_range
        end
    end

    return sliders
end

"""Create sliders for epoched data (extreme slider only)."""
function _create_sliders(fig, state::EpochedDataBrowserState, dat)
    return _create_extreme_slider(fig, state)
end

"""Popup: select which extra channels to overlay."""
function _show_extra_channel_menu(state, ax, dat)
    all_extras = extra_labels(dat)
    isnothing(all_extras) || isempty(all_extras) && return

    menu_fig = Figure(size = (900, 700))
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

    display(GLMakie.Screen(), menu_fig)
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
    _clear_and_save_limits!(ax, state)
    current = state.data.current_epoch[]
    total = n_epochs(state.data.original)
    state.data.current_epoch[] = clamp(current + direction, 1, total)
    ax.title = "Epoch $(state.data.current_epoch[])/$total"
    _update_markers!(ax, state)
    _draw(ax, state)
    _draw_extra_channel!(ax, state)
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
    state.selection.active[] = false
    state.selection.visible[] = true
    state.selection.bounds[] = (state.selection.bounds[][1], mouse_x)
    _update_x_region_selection!(ax, state, state.selection.bounds[][1], mouse_x)
    state.selection.rectangle.visible[] = true

    # Add this selection to the list of selected regions
    _add_region_to_selection!(ax, state, state.selection.bounds[][1], mouse_x)

    # Clear the temporary selection rectangle after adding to permanent regions
    _clear_x_region_selection!(state)
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
                if shift_pressed[]
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
    tolerance = 10  # pixel tolerance for snap detection

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
            channel_y = current_data[x_idx, col] + state.view.offset[idx]

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

"""Toggle the selected/highlighted state of one channel and redraw."""
function _toggle_channel_visibility!(ax, state, channel_idx)
    # Toggle the selection of the clicked channel
    state.channels.selected[channel_idx] = !state.channels.selected[channel_idx]

    # Immediate redraw for responsive feedback
    _clear_and_save_limits!(ax, state)
    _draw(ax, state)
end

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
            if state.selection.visible[]
                _handle_selection_movement!(ax, state, action)
            else
                _handle_navigation!(ax, state, action)
            end
        end
    end
end

"""Move or zoom a selection region with arrow keys."""
function _handle_selection_movement!(ax, state, action::Symbol)
    if action in (:left, :right)
        _handle_selection_movement_impl(ax, state, action)
    elseif action in (:up, :down)
        _handle_navigation!(ax, state, action)
    end
end

"""Shift the selected region left or right by 1/5 of its width."""
function _handle_selection_movement_impl(ax, state::DataBrowserState{<:AbstractDataState}, action::Symbol)
    width = state.selection.bounds[][2] - state.selection.bounds[][1]
    time_start, time_end = _get_time_bounds(state.data)
    if action == :left
        new_start = max(time_start, state.selection.bounds[][1] - width / 5)
    else  # :right
        new_start = min(time_end - width, state.selection.bounds[][1] + width / 5)
    end
    state.selection.bounds[] = (new_start, new_start + width)
    _update_x_region_selection!(ax, state, state.selection.bounds[][1], state.selection.bounds[][2])
end

"""Redraw the temporary selection rectangle between x1 and x2."""
function _update_x_region_selection!(ax, state, x1, x2)
    ylims = ax.limits[][2]
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
    ylims = ax.limits[][2]
    region_points = Point2f[
        Point2f(Float64(x1), Float64(ylims[1])),
        Point2f(Float64(x2), Float64(ylims[1])),
        Point2f(Float64(x2), Float64(ylims[2])),
        Point2f(Float64(x1), Float64(ylims[2])),
    ]
    region_plot = poly!(ax, region_points, color = (state.plot_kwargs[:selection_color], state.plot_kwargs[:selection_alpha]), strokecolor = :transparent)
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

"""Return a subset of the current continuous data within the selected region."""
function _subset_selected_data(state::ContinuousDataBrowserState, clicked_region_idx = nothing)
    # Use the clicked region if specified, otherwise use the most recent region, or fall back to bounds
    if !isnothing(clicked_region_idx) && 1 <= clicked_region_idx <= length(state.selection.selected_regions[])
        # Use the specific clicked region
        x_min, x_max = state.selection.selected_regions[][clicked_region_idx]
    elseif !isempty(state.selection.selected_regions[])
        # Use the last (most recent) selected region
        x_min, x_max = state.selection.selected_regions[][end]
    else
        # Fall back to the old bounds format for backward compatibility
        x_min, x_max = minmax(state.selection.bounds[]...)
    end

    selected_channels = state.channels.labels[state.channels.visible]
    return subset(
        state.data.current[],
        sample_selection = x -> (x.time .>= x_min) .& (x.time .<= x_max),
        channel_selection = channels(selected_channels),
    )
end

"""Return a subset of the current epoch data within the selected region."""
function _subset_selected_data(state::EpochedDataBrowserState, clicked_region_idx = nothing)
    # Get the selected region
    if !isnothing(clicked_region_idx) && 1 <= clicked_region_idx <= length(state.selection.selected_regions[])
        x_min, x_max = state.selection.selected_regions[][clicked_region_idx]
    elseif !isempty(state.selection.selected_regions[])
        x_min, x_max = state.selection.selected_regions[][end]
    else
        @minimal_warning "No region selected"
        return nothing
    end

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
    if isempty(epoch_data.data) || isempty(epoch_data.data[1])
        @minimal_warning "Selected region contains no data"
        return nothing
    end

    return epoch_to_continuous(epoch_data, 1)
end


############
# Filtering
############
"""Apply a single HP or LP filter to the current data in-place."""
function _apply_filter!(state::DataBrowserState{T}, filter_type, freq, method, order, func) where {T<:AbstractDataState}
    # Get the current data, apply filter, then update the observable
    current_data = state.data.current[]
    if filter_type == :hp
        highpass_filter!(current_data, freq; filter_method=method, order=order, filter_func=func)
    elseif filter_type == :lp
        lowpass_filter!(current_data, freq; filter_method=method, order=order, filter_func=func)
    end
    state.data.current[] = current_data  # Explicitly update the observable
end

"""Reset to original, re-reference, then re-apply all active filters and ICA removals."""
function _apply_filters!(state)
    # Reset to original if no filters active
    if !state.data.filter_state.hp_active[] && !state.data.filter_state.lp_active[]
        _reset_to_original!(state.data)
        _rereference!(state.data, state.reference_state)
        _reapply_all_ica_removals!(state)
        _notify_data_update(state.data)
        return
    end

    # Always start with fresh data when applying filters to ensure clean filtering
    _reset_to_original!(state.data)
    _rereference!(state.data, state.reference_state)

    # Apply active filters
    _reapply_active_filters!(state)

    # Re-apply ICA removals that were wiped by _reset_to_original!
    _reapply_all_ica_removals!(state)

    _notify_data_update(state.data)

end

"""Toggle the high-pass filter on/off and reapply."""
function _apply_hp_filter!(state)
    state.data.filter_state.hp_active[] = !state.data.filter_state.hp_active[]
    _apply_filters!(state)
    _update_analysis_settings!(state)
end

"""Toggle the low-pass filter on/off and reapply."""
function _apply_lp_filter!(state)
    state.data.filter_state.lp_active[] = !state.data.filter_state.lp_active[]
    _apply_filters!(state)
    _update_analysis_settings!(state)
end

########################
# Reference
########################
"""Re-reference the current data to the given reference."""
function _rereference!(state::AbstractDataState, ref)
    rereference!(state.current[], ref)
    _notify_data_update(state)  # Notify that data has been updated
end

########################
# Drawing
########################
"""Add a vertical marker (trigger/EOG) to the marker list."""
function _add_marker!(markers, ax, data, col; label = nothing, trial = nothing, visible = false, active_values = nothing)
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

    label = isnothing(label) ? string.(marker_data[!, col]) : repeat([label], nrow(marker_data))

    push!(
        markers,
        Marker(
            marker_data,
            vlines!(ax, marker_data.time, color = :grey, linewidth = 1, visible = visible),
            text!(
                ax,
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
function _init_markers(ax, state; marker_visible = Dict{Symbol,Bool}())
    markers = Marker[]
    data = _get_current_data(state.data)

    # Define marker configurations
    marker_configs = [(:trigger, nothing), (:is_vEOG, "v"), (:is_hEOG, "h")]

    # Add markers based on configuration
    for (symbol, label) in marker_configs
        if _has_column(state.data, symbol)
            active_vals = symbol == :trigger ? get(state.plot_kwargs, :active_triggers, nothing) : nothing
            _add_marker!(markers, ax, data, symbol, label = label, visible = get(marker_visible, symbol, false), active_values = active_vals)
        end
    end

    return markers
end

"""Re-create markers for the current epoch (preserving visibility)."""
function _update_markers!(ax, state)
    marker_visible = Dict{Symbol,Bool}()
    for marker in state.markers
        marker_visible[marker.name] = marker.visible
        delete!(ax, marker.line)
        delete!(ax, marker.text)
    end
    empty!(state.markers)
    state.markers = _init_markers(ax, state; marker_visible = marker_visible)
end

"""Toggle butterfly mode on/off and redraw channels."""
function _butterfly_plot!(ax, state)
    state.view.butterfly[] = !state.view.butterfly[]
    _clear_and_save_limits!(ax, state)
    _update_channel_offsets!(state)
    _draw(ax, state)
end

# Single function with data access abstraction
"""Draw all visible channel lines and labels on the axis."""
function _draw(ax, state::DataBrowserState{<:AbstractDataState})
    # Get data access functions based on type
    get_data, get_time, get_label_y = _get_data_accessors(state.data)

    # Pre-compute shared observables
    visible_time_obs = @lift(get_time($(state.data.current), $(state.view.xrange)))
    time_start_obs   = @lift($(visible_time_obs)[1])  # Fix 6: derived from visible_time_obs

    ch_idx_map = state.ica.is_active ? Dict(lbl => i for (i, lbl) in enumerate(state.ica.original.layout.data.label)) : Dict{Symbol,Int}()

    ica_activations_obs = if state.ica.is_active
        @lift begin
            win   = $(state.view.xrange)
            comps = state.ica.removed_components
            n_win = length(win)
            if isempty(comps)
                zeros(Float64, 0, n_win)
            else
                acts = zeros(Float64, length(comps), n_win)
                for (k, comp) in enumerate(comps)
                    unmix_vec = state.ica.original.unmixing[comp, :]
                    for (ci, col_sym) in enumerate(state.ica.original.layout.data.label)
                        norm_ch = (get_data(state.data.original, win, col_sym) .- state.ica.original.mean[ci]) ./ state.ica.original.scale
                        @views acts[k, :] .+= unmix_vec[ci] .* norm_ch
                    end
                end
                acts
            end
        end
    else
        Observable(zeros(Float64, 0, 0))
    end

    for (idx, visible) in enumerate(state.channels.visible)
        col = state.channels.labels[idx]
        if visible
            is_selected = state.channels.selected[idx]

            # Channel data (compute once)
            channel_data_obs = @lift(get_data($(state.data.current), $(state.view.xrange), col))
            channel_data_with_offset = @lift($(channel_data_obs) .* $(state.view.amplitude_scale) .+ state.view.offset[idx])

            # Ghost lines data
            if state.ica.is_active
                ch_row = ch_idx_map[col]

                original_data_with_offset = @lift(if $(state.view.show_original_ica)
                    current_raw  = get_data($(state.data.current), $(state.view.xrange), col)
                    original_raw = copy(current_raw)
                    acts         = $(ica_activations_obs)
                    if size(acts, 2) == length(current_raw)
                        for (k, comp) in enumerate(state.ica.removed_components)
                            if k <= size(acts, 1)
                                mix_weight = state.ica.original.mixing[ch_row, comp]
                                @views original_raw .+= mix_weight .* state.ica.original.scale .* acts[k, :]
                            end
                        end
                    end
                    original_raw .* $(state.view.amplitude_scale) .+ state.view.offset[idx]
                else
                    fill(NaN, length($(state.view.xrange)))
                end)

                subtracted_data_obs = @lift(if $(state.view.show_subtracted_ica)
                    current_raw  = get_data($(state.data.current), $(state.view.xrange), col)
                    original_raw = copy(current_raw)
                    acts         = $(ica_activations_obs)
                    if size(acts, 2) == length(current_raw)
                        for (k, comp) in enumerate(state.ica.removed_components)
                            if k <= size(acts, 1)
                                mix_weight = state.ica.original.mixing[ch_row, comp]
                                @views original_raw .+= mix_weight .* state.ica.original.scale .* acts[k, :]
                            end
                        end
                    end
                    (original_raw .- current_raw) .* $(state.view.amplitude_scale) .+ state.view.offset[idx]
                else
                    fill(NaN, length($(state.view.xrange)))
                end)

                _create_line!(
                    state.channels.original_lines,
                    col,
                    ax,
                    visible_time_obs,
                    original_data_with_offset,
                    :grey,
                    [:grey],
                    1, # linewidth
                    0.5, # clear alpha
                    visible = state.view.show_original_ica,
                )

                _create_line!(
                    state.channels.subtracted_lines,
                    col,
                    ax,
                    visible_time_obs,
                    subtracted_data_obs,
                    :red,
                    [:red],
                    1, # linewidth
                    0.5, # clear alpha
                    visible = state.view.show_subtracted_ica,
                )
            end

            is_repaired = _is_channel_repaired(state, col)

            # Line properties (reuse channel_data_obs for efficiency)
            if is_repaired || is_selected
                # Repaired channels get black color and thicker lines
                line_color = state.plot_kwargs[:selected_channel_color]
                line_colormap = [:black]
                line_width = state.plot_kwargs[:channel_line_width] * 2  # Make repaired channels thicker
            else
                # Normal channels
                line_color = @lift(abs.($(channel_data_obs)) .>= $(state.view.crit_val))

                if !state.ica.is_active
                    line_colormap = [state.plot_kwargs[:unselected_channel_color], state.plot_kwargs[:unselected_channel_color], :red]
                else
                    base_color = state.plot_kwargs[:unselected_channel_color]
                    line_colormap = @lift($(state.view.show_original_ica) ? [:green, :green, :red] : [base_color, base_color, :red])
                end

                line_width = state.plot_kwargs[:channel_line_width]
            end

            # Update or create line
            _create_line!(
                state.channels.data_lines,
                col,
                ax,
                visible_time_obs,
                channel_data_with_offset,
                line_color,
                line_colormap,
                line_width,
                state.plot_kwargs[:channel_line_alpha];
                visible = get(state.plot_kwargs, :show_cleaned_ica, Observable(true)),
            )

            # Handle labels
            if !state.view.butterfly[]
                label_y_obs = @lift(get_label_y($(state.data.current), $col, state.view.offset[idx]))
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

"""Save current axis limits and clear all channel lines/labels."""
function _clear_and_save_limits!(ax, state)
    state.plot_kwargs[:_saved_xlims] = ax.xaxis.attributes.limits[]
    state.plot_kwargs[:_saved_ylims] = ax.yaxis.attributes.limits[]
    _clear_axes!(ax, _all_line_dicts(state))
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
    get_data = (current, range, col) -> current.data[range, col]
    get_time = (current, range) -> current.data.time[range]
    get_label_y = (current, col, offset) -> current.data[1, col] .+ offset
    return get_data, get_time, get_label_y
end

"""Return `(get_data, get_time, get_label_y)` closures for epoched data."""
function _get_data_accessors(state::EpochedDataState)
    get_data = (current, range, col) -> current.data[state.current_epoch[]][range, col]
    get_time = (current, range) -> current.data[state.current_epoch[]].time[range]
    get_label_y = (current, col, offset) -> current.data[state.current_epoch[]][!, col][1] .+ offset
    return get_data, get_time, get_label_y
end

"""Plot a channel line on the axis and store it."""
function _create_line!(data_lines, col, ax, x_obs, y_obs, color, colormap, linewidth, alpha; visible = true)
    data_lines[col] = lines!(ax, x_obs, y_obs, color = color, colormap = colormap, linewidth = linewidth, alpha = alpha, visible = visible)
end

"""Add a channel-name text label on the axis."""
function _create_label!(data_labels, col, ax, x_obs, y_obs, is_selected)
    data_labels[col] =
        text!(ax, x_obs, y_obs, text = String(col), align = (:left, :center), fontsize = 18, color = is_selected ? :red : :black)
end

"""Hide a single channel label."""
function _hide_channel_label!(data_labels, col)
    haskey(data_labels, col) && hide!(data_labels[col])
end

"""Hide both line and label for a channel."""
function _hide_channel_objects!(channels, col)
    haskey(channels.data_lines, col) && hide!(channels.data_lines[col])
    haskey(channels.data_labels, col) && hide!(channels.data_labels[col])
    haskey(channels.original_lines, col) && hide!(channels.original_lines[col])
    haskey(channels.subtracted_lines, col) && hide!(channels.subtracted_lines[col])
end

# Single function with data access abstraction
"""Draw the extra-channel overlay (line or boolean highlight)."""
function _draw_extra_channel!(ax, state::DataBrowserState{<:AbstractDataState})
    _clear_axes!(ax, [state.extra_channel.data_lines, state.extra_channel.data_labels])

    # Remove old legend if present
    if !isnothing(state.extra_channel.legend)
        delete!(state.extra_channel.legend)
        state.extra_channel.legend = nothing
    end

    if !isempty(state.extra_channel.channels)
        if length(state.view.offset) > 1
            current_offset = state.view.offset[end] + (state.view.offset[end] - state.view.offset[end-1])
        else
            current_offset = state.view.offset[end] + 100.0  # Default spacing
        end

        # Get data access functions based on type
        get_data, get_time, get_label_y = _get_data_accessors(state.data)

        colors = Makie.wong_colors()

        # Accumulators for the legend
        legend_elements = Any[]
        legend_labels   = String[]

        for (idx, channel) in enumerate(state.extra_channel.channels)
            # Calculate specific offset for this channel to prevent overlapping
            channel_offset = current_offset
            if length(state.view.offset) > 1
                channel_offset = state.view.offset[end] + idx * (state.view.offset[end] - state.view.offset[end-1])
            else
                channel_offset = state.view.offset[end] + idx * 100.0
            end

            line_color = colors[mod1(idx, length(colors))]

            if eltype(get_data(state.data.current[], 1:1, channel)) == Bool # Boolean data - vspan + legend entry only
                highlight_data = @views _splitgroups(findall(get_data(state.data.current[], :, channel)))
                if !isempty(highlight_data[1])
                    # Widen single-sample regions so they're visible (per-region, not global)
                    end_indices = copy(highlight_data[2])
                    for i in eachindex(end_indices)
                        if end_indices[i] == highlight_data[1][i]
                            end_indices[i] += 5
                        end
                    end
                    state.extra_channel.data_lines[channel] = vspan!(
                        ax,
                        get_time(state.data.current[], highlight_data[1]),
                        get_time(state.data.current[], end_indices),
                        color = line_color,
                        alpha = 0.5,
                        visible = true,
                    )
                end
                # Boolean channels identified via the legend only — no in-plot text
                push!(legend_elements, PolyElement(color = (line_color, 0.5), strokecolor = :transparent))
                push!(legend_labels, String(channel))
            else # Regular channel - in-plot label is sufficient, skip legend
                state.extra_channel.data_lines[channel] = lines!(
                    ax,
                    @lift(get_time($(state.data.current), :)),
                    @lift(get_data($(state.data.current), :, $channel) .* $(state.view.amplitude_scale) .+ $channel_offset),
                    color = line_color,
                    linewidth = 2,
                )
                state.extra_channel.data_labels[channel] = text!(
                    ax,
                    @lift(get_time($(state.data.current), $(state.view.xrange)[1:1])[1]),
                    @lift(get_label_y($(state.data.current), $channel, $channel_offset)),
                    text = String(channel),
                    align = (:left, :center),
                    color = line_color,
                    fontsize = 18,
                )
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

"""Toggle visibility of a vertical marker and update its label positions."""
function _plot_vertical_lines!(ax, marker, active)
    marker.line.visible = active
    marker.text.visible = active
    marker.visible = active
    marker.text.position = [(x, ax.yaxis.attributes.limits[][2] * 0.98) for x in marker.data.time] # incase y changed
end

"""
    _add_scale_indicator!(ax, state, plot_kwargs)

Add a scale indicator bar to the plot showing the amplitude scale.
"""
function _add_scale_indicator!(ax, state, plot_kwargs)

    scale_value = plot_kwargs[:scale_indicator_value]
    pos = plot_kwargs[:scale_indicator_position]

    # Get axis limits observables (they're already observables)
    xlims_obs = ax.xaxis.attributes.limits
    ylims_obs = ax.yaxis.attributes.limits

    # Calculate position in data coordinates
    # pos[1] is x position (0 = left, 1 = right), pos[2] is y position (0 = bottom, 1 = top)
    x_pos = @lift($xlims_obs[1] + ($xlims_obs[2] - $xlims_obs[1]) * pos[1])
    # Position scale bar at the top of the plot
    y_top = @lift($ylims_obs[1] + ($ylims_obs[2] - $ylims_obs[1]) * pos[2])
    y_bottom = @lift($y_top - scale_value * $(state.view.amplitude_scale))

    # Draw vertical line
    lines!(ax, @lift([$x_pos, $x_pos]), @lift([$y_bottom, $y_top]), color = :black, linewidth = 1)

    # Draw horizontal tick marks
    tick_length = @lift(($xlims_obs[2] - $xlims_obs[1]) * 0.005)
    tick_left = @lift($x_pos - $tick_length)
    tick_right = @lift($x_pos + $tick_length)
    lines!(ax, @lift([$tick_left, $tick_right]), @lift([$y_bottom, $y_bottom]), color = :black, linewidth = 1)
    lines!(ax, @lift([$tick_left, $tick_right]), @lift([$y_top, $y_top]), color = :black, linewidth = 1)

    # Add label
    label_x = @lift($x_pos + ($xlims_obs[2] - $xlims_obs[1]) * 0.01)
    text!(
        ax,
        @lift(Point2f($label_x, $y_top)),
        text = "$(round(scale_value, digits=0)) μV",
        align = (:left, :center),
        fontsize = 14,
        color = :black,
        space = :data,
    )

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

    # Merge user kwargs with defaults
    plot_kwargs = _merge_plot_kwargs(PLOT_DATABROWSER_KWARGS, kwargs)

    # Common fig/ax/state/ui setup
    fig = Figure(figure_padding = plot_kwargs[:figure_padding])
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
    if plot_kwargs[:display_plot]
        if !isnothing(screen)
            display(screen, fig)
        else
            display(fig)
        end
    end

    _set_window_title("Makie")
    # Return the observable analysis settings
    return (fig = fig, ax = ax, analysis_settings = state.analysis_settings)
end

"""Load data from a `.jld2` file or pattern and open the data browser."""
function plot_databrowser(
    filename::String,
    ica = nothing;
    input_dir::String = pwd(),
    participant_selection::Function = participants(),
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
    for dat in data
        _set_window_title(_generate_window_title(dat))  # must be before Screen()
        s = GLMakie.Screen()
        plot_databrowser(dat, ica; screen = s, kwargs...)
    end
end
