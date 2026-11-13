const PLOT_DATABROWSER_KWARGS = Dict{Symbol,Tuple{Any,String}}(
    # Figure and layout
    :figure_size => ((1300, 800), "Figure size as (width, height)"),
    :figure_padding => ((50, 30, 50, 50), "Figure padding as (left, right, bottom, top)"),

    # Axis styling
    :xlabel => ("Time (s)", "X-axis label"),
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
    :remove_dc => (false, "Subtract each channel's window mean (DC offset) for display"),
    :channel_offset_scale => (1500, "Scale factor for channel vertical offset"),
    :channel_offset_margin => (0.9, "Margin factor for channel offset range"),

    # Navigation
    :window_duration => (10.0, "Initial window duration in seconds for continuous data"),
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

abstract type AbstractDataState end

mutable struct Marker
    data::Any
    line::Any
    text::Any
    name::Symbol
    visible::Bool
    observers::Vector{Any}
    Marker(data, line, text, name, visible, observers = Any[]) = new(data, line, text, name, visible, observers)
end

mutable struct SelectionState
    active::Observable{Bool}
    bounds::Observable{Tuple{Float64,Float64}}
    visible::Observable{Bool}
    rectangle::Makie.Poly
    selected_regions::Observable{Vector{Tuple{Float64,Float64}}}
    region_plots::Vector{Makie.Poly}
    function SelectionState(ax, plot_kwargs)
        poly_element =
            poly!(ax, [Point2f(0.0, 0.0)], color = (plot_kwargs[:selection_color], plot_kwargs[:selection_alpha]), visible = false)
        new(Observable(false), Observable((0.0, 0.0)), Observable(false), poly_element, Observable(Tuple{Float64,Float64}[]), Makie.Poly[])
    end
end

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
    remove_dc::Observable{Bool}
    amplitude_scale::Observable{Float64}
    show_original_ica::Observable{Bool}
    show_subtracted_ica::Observable{Bool}
    function ViewState(n_channels::Int, n_samples::Int, plot_kwargs)
        offset_scale = plot_kwargs[:channel_offset_scale]
        offset_margin = plot_kwargs[:channel_offset_margin]
        offset =
            n_channels > 1 ? collect(LinRange(offset_scale * offset_margin, -offset_scale * offset_margin, n_channels + 2)[2:(end-1)]) :
            zeros(n_channels)
        new(
            Observable(1:n_samples),
            Observable(-1500:1500),
            offset,
            Observable(0.0),
            Observable(false),
            Observable(Bool(get(plot_kwargs, :remove_dc, false))),
            Observable(1.0),
            Observable(false),
            Observable(false),
        )
    end
end

mutable struct ChannelState
    labels::Vector{Symbol}
    visible::Vector{Bool}
    selected::Vector{Bool}
    data_labels::Dict{Symbol,Makie.Text}
    data_lines::Dict{Symbol,Any}
    original_lines::Dict{Symbol,Any}
    subtracted_lines::Dict{Symbol,Any}
    draw_observers::Vector{Any} # Disconnected on redraw to avoid stacking listeners
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

mutable struct ContinuousDataState <: AbstractDataState
    current::Observable{EegData}
    original::EegData
    filter_state::FilterState
    ContinuousDataState(data::EegData, plot_kwargs) = new(Observable(copy(data)), data, FilterState(plot_kwargs))
end

mutable struct EpochedDataState <: AbstractDataState
    current::Observable{EegData}
    original::EegData
    filter_state::FilterState
    current_epoch::Observable{Int}
    EpochedDataState(data::EegData, plot_kwargs) = new(Observable(copy(data)), data, FilterState(plot_kwargs), Observable(1))
end

mutable struct ExtraChannelInfo
    channels::Vector{Symbol}
    data_lines::Dict{Symbol,Any}
    data_labels::Dict{Symbol,Makie.Text}
    fig::Union{Nothing,Figure}
    legend::Union{Nothing,Any}
    draw_observers::Vector{Any}
    ExtraChannelInfo() = new(Symbol[], Dict{Symbol,Any}(), Dict{Symbol,Makie.Text}(), nothing, nothing, Any[])
end

# `activation_cache` holds removed-component activations (samples × removed components) per epoch
# with one window-width margin on each side so scrolling reuses cached activations.
mutable struct IcaVisualState
    original::Union{Nothing,InfoIca}
    removed_components::Vector{Int}
    is_active::Bool
    activation_cache::Dict{Int,Tuple{UnitRange{Int},Matrix{Float64}}}
    activation_key::Vector{Int}
    IcaVisualState(ica::Union{Nothing,InfoIca}) = new(ica, Int[], !isnothing(ica), Dict{Int,Tuple{UnitRange{Int},Matrix{Float64}}}(), Int[])
end

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
    custom_reference_channels::Vector{Symbol}
    channel_repair_history::Vector{Tuple{Vector{Symbol},Symbol,Any}}
    analysis_settings::Observable{AnalysisSettings}
    plot_kwargs::Dict{Symbol,Any}

    function DataBrowserState{T}(;
        view::ViewState,
        channels::ChannelState,
        data::T,
        selection::SelectionState,
        ica_original::Union{Nothing,InfoIca} = nothing,
        extra_channel::ExtraChannelInfo = ExtraChannelInfo(),
        plot_kwargs::Dict{Symbol,Any} = Dict{Symbol,Any}(),
    ) where {T<:AbstractDataState}
        init_ref = data.original.analysis_info.reference
        init_custom = (init_ref ∉ (:none, :avg, :mastoid) && init_ref in channels.labels) ? Symbol[init_ref] : Symbol[]
        return new{T}(
            view,
            channels,
            data,
            selection,
            Marker[],
            Dict{Symbol,Bool}(),
            IcaVisualState(ica_original),
            extra_channel,
            init_ref,
            init_custom,
            Tuple{Vector{Symbol},Symbol,Any}[],
            Observable(AnalysisSettings()),
            plot_kwargs,
        )
    end
end

const ContinuousDataBrowserState = DataBrowserState{ContinuousDataState}
const EpochedDataBrowserState = DataBrowserState{EpochedDataState}

function _update_analysis_settings!(state::DataBrowserState)
    fs = state.data.filter_state
    hp_freq = fs.hp_active[] ? fs.hp_freq[] : 0.0
    lp_freq = fs.lp_active[] ? fs.lp_freq[] : 0.0

    repaired_channels = Symbol[]
    for (channels, _, _) in state.channel_repair_history
        append!(repaired_channels, channels)
    end
    repair_method = isempty(state.channel_repair_history) ? :none : last(state.channel_repair_history)[2]

    state.analysis_settings[] = AnalysisSettings(
        hp_freq,
        lp_freq,
        state.reference_state,
        repaired_channels,
        repair_method,
        copy(state.selection.selected_regions[]),
        copy(state.ica.removed_components),
    )
end

function _mark_regions_col!(mask::BitVector, times::AbstractVector{<:Real}, selected_regions::Vector{Tuple{Float64,Float64}})
    @inbounds for (t0, t1) in selected_regions
        for i in eachindex(times, mask)
            t = times[i]
            (t >= t0 && t <= t1) && (mask[i] = true)
        end
    end
    return mask
end

function _mark_regions!(df::DataFrame, selected_regions::Vector{Tuple{Float64,Float64}})
    mask = falses(nrow(df))
    !isempty(selected_regions) && _mark_regions_col!(mask, df.time, selected_regions)
    df[!, :selected_region] = mask
end

_add_selected_regions!(dat::EegData, regions::Vector{Tuple{Float64,Float64}}) = _mark_regions!(dat.data, regions)
_add_selected_regions!(dat::EpochData, regions::Vector{Tuple{Float64,Float64}}) = foreach(df -> _mark_regions!(df, regions), dat.data)

function _create_browser_state(dat::T, channel_labels, ax, ica, plot_kwargs) where {T<:EegData}
    if !haskey(plot_kwargs, :show_cleaned_ica) || !(plot_kwargs[:show_cleaned_ica] isa Observable)
        val = get(plot_kwargs, :show_cleaned_ica, false)
        plot_kwargs[:show_cleaned_ica] = Observable{Bool}(val isa Observable ? val[] : val)
    end
    state_type = _data_state_type(T)
    return DataBrowserState{state_type}(
        view = ViewState(length(channel_labels), _get_initial_window_size(dat, plot_kwargs), plot_kwargs),
        channels = ChannelState(channel_labels),
        data = state_type(dat, plot_kwargs),
        selection = SelectionState(ax, plot_kwargs),
        ica_original = ica,
        plot_kwargs = plot_kwargs,
    )
end

function _get_initial_window_size(dat::ContinuousData, plot_kwargs = Dict{Symbol,Any}())
    dur_sec = Float64(get(plot_kwargs, :window_duration, 10.0))
    sr = dat.sample_rate > 0 ? Float64(dat.sample_rate) : 500.0
    return clamp(round(Int, dur_sec * sr), min(2, n_samples(dat)), n_samples(dat))
end
_get_initial_window_size(dat::ErpData, plot_kwargs = Dict{Symbol,Any}()) = n_samples(dat)
_get_initial_window_size(dat::EpochData, plot_kwargs = Dict{Symbol,Any}()) = nrow(dat.data[1])

_data_state_type(::Type{ContinuousData}) = ContinuousDataState
_data_state_type(::Type{EpochData}) = EpochedDataState
_data_state_type(::Type{ErpData}) = ContinuousDataState

_get_current_data(state::ContinuousDataState) = state.current[].data
_get_current_data(state::EpochedDataState) = state.current[].data[state.current_epoch[]]
_has_column(state::ContinuousDataState, col::Symbol) = hasproperty(state.current[].data, col)
_has_column(state::EpochedDataState, col::Symbol) = hasproperty(state.current[].data[state.current_epoch[]], col)

# =============================================================================
# UI
# =============================================================================

function _setup_ui_base(fig, ax, state, dat, ica = nothing)
    deregister_interaction!(ax, :rectanglezoom)
    _set_axes!(ax, state)
    _handle_mouse_events!(ax, state)
    _handle_keyboard_events!(fig, ax, state)

    state.markers = _init_markers(ax, state)
    toggles = _create_toggles(fig, ax, state, dat)
    labels_menu = _create_labels_menu(fig, ax, state)
    reference_menu = _create_reference_menu(fig, state, dat)

    ica_menu = (!isnothing(ica) && state.ica.is_active) ? _create_ica_menu(fig, ax, state, ica) : nothing
    extra_lbls = extra_labels(state.data.original)
    extra_menu = (!isnothing(extra_lbls) && !isempty(extra_lbls)) ? _create_extra_channel_menu(fig, ax, state, dat) : nothing

    return (toggles, labels_menu, reference_menu, ica_menu, extra_menu)
end

function _setup_ui(fig, ax, state::DataBrowserState{<:AbstractDataState}, dat, ica = nothing, plot_kwargs = nothing)
    toggles, labels_menu, reference_menu, ica_menu, extra_menu = _setup_ui_base(fig, ax, state, dat, ica)
    epoch_menu = _get_epoch_menu(fig, ax, state)
    _build_grid_components!(fig, dat, state, toggles, labels_menu, reference_menu, ica_menu, extra_menu, epoch_menu)
    hideydecorations!(ax, label = true)
    return state
end

function _create_toggles(fig, ax, state, dat)
    configs = [
        ToggleConfig("Remove DC", active -> (state.view.remove_dc[] = active)),
        ToggleConfig("Butterfly Plot", _ -> _butterfly_plot!(ax, state)),
    ]

    for (marker_symbol, toggle_label) in ((:trigger, "Trigger"), (:is_vEOG, "vEOG"), (:is_hEOG, "hEOG"))
        if _has_column(state.data, marker_symbol)
            state.marker_visibility[marker_symbol] = false
            push!(configs, ToggleConfig(toggle_label, active -> _toggle_marker_visibility!(ax, state, marker_symbol, active)))
        end
    end

    if state.data.original.analysis_info.hp_filter == 0.0
        push!(configs, ToggleConfig("HP-Filter", active -> _apply_hp_filter!(state, active)))
    end
    if state.data.original.analysis_info.lp_filter == 0.0
        push!(configs, ToggleConfig("LP-Filter", active -> _apply_lp_filter!(state, active)))
    end

    toggles_grid = Vector{Any}[]
    for config in configs
        init_active = config.label == "Remove DC" ? state.view.remove_dc[] : false
        tog = Toggle(fig, active = init_active)
        on(config.action, tog.active)

        if config.label == "Trigger"
            btn = Button(fig, label = "Trigger", halign = :left)
            on(_ -> _show_trigger_menu(state, ax, :trigger), btn.clicks)
            push!(toggles_grid, Any[tog, btn])
        elseif config.label == "HP-Filter" || config.label == "LP-Filter"
            btn = Button(fig, label = config.label, halign = :left)
            filter_type = config.label == "HP-Filter" ? :hp : :lp
            on(_ -> _show_single_filter_menu(state, dat, filter_type), btn.clicks)
            push!(toggles_grid, Any[tog, btn])
        else
            push!(toggles_grid, Any[tog, Label(fig, config.label, fontsize = 14, halign = :left)])
        end
    end

    return permutedims(reduce(hcat, toggles_grid))
end

function _toggle_marker_visibility!(ax, state, marker_symbol::Symbol, active::Bool)
    state.marker_visibility[marker_symbol] = active
    marker_index = findfirst(m -> m.name == marker_symbol, state.markers)
    !isnothing(marker_index) && _plot_vertical_lines!(ax, state.markers[marker_index], active)
end

function _show_trigger_menu(state, ax, marker_symbol)
    trigger_col = if state.data isa ContinuousDataState
        state.data.original.data[!, marker_symbol]
    else
        reduce(vcat, (df[!, marker_symbol] for df in state.data.original.data))
    end
    unique_triggers = sort(unique(filter(!iszero, trigger_col)))

    active_set = get!(state.plot_kwargs, :active_triggers) do
        Set(unique_triggers)
    end

    menu_fig = _eegfun_figure(size = (600, 600))
    Label(menu_fig[1, 1], "Select Triggers to Display", fontsize = 18, halign = :center)

    scroll_area = menu_fig[2, 1] = GridLayout()
    trig_strs = string.(unique_triggers)
    cbs = _build_checkbox_grid!(scroll_area, trig_strs, 5, (trig, _) -> (parse(Float64, trig) in active_set))
    _add_group_buttons!(menu_fig, 3, cbs, trig_strs, [("All", _ -> true), ("None", _ -> false)])

    action_area = menu_fig[4, 1] = GridLayout()
    btn_apply = Button(action_area[1, 1], label = "Apply", width = 200)
    on(btn_apply.clicks) do _
        empty!(active_set)
        for (i, cb) in enumerate(cbs)
            cb.checked[] && push!(active_set, unique_triggers[i])
        end
        _update_markers!(ax, state)
    end

    _display_in_screen(menu_fig)
end

function _show_single_filter_menu(state, dat, filter_type::Symbol)
    menu_fig = _eegfun_figure(size = (450, 350))
    title_str = filter_type == :hp ? "Highpass Filter Settings" : "Lowpass Filter Settings"
    Label(menu_fig[1, 1], title_str, fontsize = 18, halign = :center, font = :bold)

    grid = menu_fig[2, 1] = GridLayout(valign = :top)
    if filter_type == :hp
        range_freq, freq_field, method_field, order_field, func_field, default_val =
            0.1:0.1:2, :hp_freq, :hp_method, :hp_order, :hp_func, 0.5
    else
        range_freq, freq_field, method_field, order_field, func_field, default_val = 5:5:60, :lp_freq, :lp_method, :lp_order, :lp_func, 20.0
    end

    fs = state.data.filter_state
    getfield(fs, freq_field)[] == 0.0 && (getfield(fs, freq_field)[] = default_val)

    Label(grid[1, 1], "Cutoff (Hz):", halign = :left)
    slider_freq = Slider(grid[1, 2], range = range_freq, startvalue = getfield(fs, freq_field)[])
    Label(grid[1, 3], @lift(string(round($(slider_freq.value), digits = 1))), halign = :left)
    on(val -> (getfield(fs, freq_field)[] = val), slider_freq.value)

    Label(grid[2, 1], "Order:", halign = :left)
    slider_order = Slider(grid[2, 2], range = 1:10, startvalue = getfield(fs, order_field)[])
    Label(grid[2, 3], @lift(string($(slider_order.value))), halign = :left)
    on(val -> (getfield(fs, order_field)[] = val), slider_order.value)

    Label(grid[3, 1], "Method:", halign = :left)
    current_method = getfield(fs, method_field)[]
    menu_method = Menu(
        grid[3, 2],
        options = zip(["Butterworth (IIR)", "FIR (Hamming)"], ["iir", "fir"]),
        default = current_method == "iir" ? "Butterworth (IIR)" : "FIR (Hamming)",
        width = 220,
    )
    on(val -> (getfield(fs, method_field)[] = val), menu_method.selection)

    Label(grid[4, 1], "Function:", halign = :left)
    current_func = getfield(fs, func_field)[]
    func_labels = ["filtfilt (2-pass zero-phase)", "filt_zerophase (1-pass FIR)", "filt (1-pass causal)"]
    func_vals = ["filtfilt", "filt_zerophase", "filt"]
    default_func_label = current_func == "filt_zerophase" ? func_labels[2] : (current_func == "filt" ? func_labels[3] : func_labels[1])
    menu_func = Menu(grid[4, 2], options = zip(func_labels, func_vals), default = default_func_label, width = 220)
    on(val -> (getfield(fs, func_field)[] = val), menu_func.selection)

    rowgap!(grid, 10)
    colgap!(grid, 15)

    action_area = menu_fig[3, 1] = GridLayout()
    btn_apply = Button(action_area[1, 1], label = "Apply", width = 150)
    on(_ -> _apply_filters!(state), btn_apply.clicks)

    _display_in_screen(menu_fig)
end

const _CHANNEL_HEMISPHERE_GROUPS = Tuple{String,Function}[
    ("Left", ch -> occursin(r"\d*[13579]$", String(ch))),
    ("Right", ch -> occursin(r"\d*[24680]$", String(ch))),
    ("Central", ch -> occursin(r"z$", String(ch))),
]

function _show_labels_menu(state, ax)
    all_channels = state.channels.labels
    menu_fig = _eegfun_figure(size = (900, 700))
    Label(menu_fig[1, 1], "Select Channels to Display", fontsize = 18, halign = :center)
    scroll_area = menu_fig[2, 1] = GridLayout()
    cbs = _build_checkbox_grid!(scroll_area, all_channels, 9, (_, i) -> state.channels.visible[i])

    groups = vcat(Tuple{String,Function}[("All", _ -> true), ("None", _ -> false)], _CHANNEL_HEMISPHERE_GROUPS)
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
        !isempty(state.extra_channel.channels) && _draw_extra_channel!(ax, state)
    end

    _display_in_screen(menu_fig)
end

_create_labels_menu(fig, ax, state) = _popup_button(fig, "Labels", () -> _show_labels_menu(state, ax))

function _show_reference_menu(state, dat)
    all_channels = state.channels.labels
    menu_fig = _eegfun_figure(size = (900, 700))
    Label(menu_fig[1, 1], "Select Reference Channels", fontsize = 18, halign = :center)

    current_syms = if state.reference_state == :avg
        all_channels
    elseif state.reference_state == :mastoid
        [:M1, :M2]
    elseif state.reference_state == :none
        Symbol[]
    elseif !isempty(state.custom_reference_channels)
        state.custom_reference_channels
    elseif state.reference_state in all_channels
        Symbol[state.reference_state]
    else
        Symbol[]
    end

    scroll_area = menu_fig[2, 1] = GridLayout()
    cbs = _build_checkbox_grid!(scroll_area, all_channels, 9, (ch, _) -> ch in current_syms)

    groups = vcat(
        Tuple{String,Function}[("All (Avg)", _ -> true), ("None", _ -> false)],
        _CHANNEL_HEMISPHERE_GROUPS,
        Tuple{String,Function}[("Mastoids", ch -> ch == :M1 || ch == :M2)],
    )
    _add_group_buttons!(menu_fig, 3, cbs, all_channels, groups)

    action_area = menu_fig[4, 1] = GridLayout()
    btn_apply = Button(action_area[1, 1], label = "Apply", width = 200)
    on(btn_apply.clicks) do _
        selected = [all_channels[i] for (i, cb) in enumerate(cbs) if cb.checked[]]
        if isempty(selected)
            state.reference_state = :none
            empty!(state.custom_reference_channels)
        elseif length(selected) == length(all_channels)
            state.reference_state = :avg
            empty!(state.custom_reference_channels)
        elseif length(selected) == 2 && :M1 in selected && :M2 in selected
            state.reference_state = :mastoid
            empty!(state.custom_reference_channels)
        elseif length(selected) == 1
            state.reference_state = selected[1]
            state.custom_reference_channels = copy(selected)
        else
            state.reference_state = :custom
            state.custom_reference_channels = copy(selected)
        end
        _rebuild_current!(state)
    end

    _display_in_screen(menu_fig)
end

_create_reference_menu(fig, state, dat) = _popup_button(fig, "Reference", () -> _show_reference_menu(state, dat))

function _apply_reference!(cur::EegData, ref_state::Symbol, custom_channels::Vector{Symbol} = Symbol[])
    ref_state == :none && return cur
    ref = if ref_state == :avg || ref_state == :mastoid
        ref_state
    elseif !isempty(custom_channels)
        length(custom_channels) == 1 ? custom_channels[1] : custom_channels
    else
        ref_state
    end
    rereference!(cur, ref)
    return cur
end

function _accum_ica_comp_1d!(act::Vector{Float64}, col::AbstractVector{<:Real}, u::Float64, m::Float64, scale::Float64)
    @inbounds for s in eachindex(act, col)
        act[s] += u * ((col[s] - m) / scale)
    end
    return nothing
end

function _sub_ica_comp_1d!(col::AbstractVector{<:AbstractFloat}, act::Vector{Float64}, w::Float64)
    @inbounds for s in eachindex(col, act)
        col[s] -= w * act[s]
    end
    return nothing
end

function _apply_df_ica_removal!(cur_df::DataFrame, orig_df::DataFrame, ica::InfoIca, comp::Int, act::Vector{Float64})
    n = nrow(cur_df)
    length(act) == n || resize!(act, n)
    fill!(act, 0.0)
    cols = ica.layout.data.label
    scale = ica.scale
    @inbounds for i in eachindex(cols)
        col_sym = cols[i]
        _accum_ica_comp_1d!(act, orig_df[!, col_sym], ica.unmixing[comp, i], ica.mean[i], scale)
    end
    @inbounds for i in eachindex(cols)
        col_sym = cols[i]
        _sub_ica_comp_1d!(cur_df[!, col_sym], act, ica.mixing[i, comp] * scale)
    end
    return nothing
end

function _apply_single_ica_removal!(cur::SingleDataFrameEeg, orig::SingleDataFrameEeg, ica::InfoIca, comp::Int)
    act = Vector{Float64}(undef, nrow(cur.data))
    _apply_df_ica_removal!(cur.data, orig.data, ica, comp, act)
    return nothing
end

function _apply_single_ica_removal!(cur::EpochData, orig::EpochData, ica::InfoIca, comp::Int)
    isempty(cur.data) && return nothing
    act = Vector{Float64}(undef, nrow(cur.data[1]))
    for ep in eachindex(cur.data)
        _apply_df_ica_removal!(cur.data[ep], orig.data[ep], ica, comp, act)
    end
    return nothing
end

function _apply_all_ica_removals!(cur::SingleDataFrameEeg, state)
    if state.ica.is_active && !isempty(state.ica.removed_components)
        ica = state.ica.original::InfoIca
        orig_df = state.data.original.data
        act = Vector{Float64}(undef, nrow(cur.data))
        for comp in state.ica.removed_components
            _apply_df_ica_removal!(cur.data, orig_df, ica, comp, act)
        end
    end
    return cur
end

function _apply_all_ica_removals!(cur::EpochData, state)
    if state.ica.is_active && !isempty(state.ica.removed_components) && !isempty(cur.data)
        ica = state.ica.original::InfoIca
        orig_eps = state.data.original.data
        act = Vector{Float64}(undef, nrow(cur.data[1]))
        for comp in state.ica.removed_components
            for ep in eachindex(cur.data)
                _apply_df_ica_removal!(cur.data[ep], orig_eps[ep], ica, comp, act)
            end
        end
    end
    return cur
end

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

# Canonical rebuild order: repairs -> reference -> ICA removals -> active filters.
function _rebuild_current!(state)
    cur = copy(state.data.original)
    for (chs, method, _) in state.channel_repair_history
        repair_channels!(cur, chs; method = method)
    end
    _apply_reference!(cur, state.reference_state, state.custom_reference_channels)
    _apply_all_ica_removals!(cur, state)
    _apply_active_filters!(cur, state)
    state.data.current[] = cur
    _update_analysis_settings!(state)
    return cur
end

function _toggle_ica_component!(state, component_id::Union{Nothing,Int})
    if isnothing(component_id)
        empty!(state.ica.removed_components)
    elseif component_id in state.ica.removed_components
        filter!(!=(component_id), state.ica.removed_components)
    else
        push!(state.ica.removed_components, component_id)
    end
    _rebuild_current!(state)
    return nothing
end

function _ica_correlatable_channels(dat::EegData)::Vector{Symbol}
    df = dat isa SingleDataFrameEeg ? dat.data : (isempty(dat.data) ? DataFrame() : dat.data[1])
    nrow(df) == 0 && return Symbol[]
    is_numeric_col(c::Symbol) = hasproperty(df, c) && (eltype(df[!, c]) <: Real) && !(eltype(df[!, c]) <: Bool)
    extras = Symbol[c for c in extra_labels(dat) if is_numeric_col(c)]
    scalp = Symbol[c for c in channel_labels(dat) if is_numeric_col(c) && c ∉ extras]
    return vcat(extras, scalp)
end

function _ica_channel_vector(dat::SingleDataFrameEeg, col::Symbol)
    hasproperty(dat.data, col) || return nothing
    v = dat.data[!, col]
    (eltype(v) <: Real && !(eltype(v) <: Bool)) || return nothing
    return Float64.(v)
end

function _ica_channel_vector(dat::EpochData, col::Symbol)
    (isempty(dat.data) || !hasproperty(dat.data[1], col)) && return nothing
    v1 = dat.data[1][!, col]
    (eltype(v1) <: Real && !(eltype(v1) <: Bool)) || return nothing
    return reduce(vcat, (Float64.(df[!, col]) for df in dat.data))
end

function _ica_activation_matrix(dat::EegData, ica::InfoIca)::Matrix{Float64}
    cols = ica.layout.data.label
    X = if dat isa SingleDataFrameEeg
        _pack_df_to_matrix(dat.data, cols)
    elseif dat isa EpochData && !isempty(dat.data)
        n_tot = sum(nrow(df) for df in dat.data)
        _pack_dfs_to_matrix(dat.data, cols, n_tot)
    else
        return Matrix{Float64}(undef, 0, length(ica.ica_label))
    end
    size(X, 1) <= 1 && return Matrix{Float64}(undef, 0, length(ica.ica_label))
    inv_scale = 1.0 / ica.scale
    @inbounds for j in axes(X, 2)
        mj = ica.mean[j]
        @simd for i in axes(X, 1)
            X[i, j] = (X[i, j] - mj) * inv_scale
        end
    end
    return X * transpose(ica.unmixing)
end

function _ica_component_metrics(
    dat::EegData,
    ica::InfoIca;
    channel::Union{Nothing,Symbol} = nothing,
    acts::Union{Nothing,Matrix{Float64}} = nothing,
)
    n_comps = length(ica.ica_label)
    var_total = sum(ica.variance)
    variance_pct = var_total > 0 ? (Float64.(ica.variance) ./ var_total .* 100.0) : zeros(Float64, n_comps)
    channel_corr = fill(NaN, n_comps)

    if !isnothing(channel)
        ch_vec = _ica_channel_vector(dat, channel)
        if !isnothing(ch_vec)
            act_mat = isnothing(acts) ? _ica_activation_matrix(dat, ica) : acts
            if size(act_mat, 1) == length(ch_vec) && length(ch_vec) > 1 && std(ch_vec) > eps(Float64)
                for c = 1:n_comps
                    sc = std(view(act_mat, :, c))
                    channel_corr[c] = sc > eps(Float64) ? abs(cor(view(act_mat, :, c), ch_vec)) : 0.0
                end
            end
        end
    end

    return (
        variance_pct = variance_pct,
        var_rank = sortperm(variance_pct, rev = true),
        channel = channel,
        channel_corr = channel_corr,
        corr_rank = all(isnan, channel_corr) ? collect(1:n_comps) : sortperm(replace(channel_corr, NaN => -1.0), rev = true),
    )
end

function _format_ica_component_label(comp_label, comp_idx::Int, metrics)
    var_str = @sprintf("%.1f%%", metrics.variance_pct[comp_idx])
    if !isnan(metrics.channel_corr[comp_idx])
        return @sprintf("%s (%s | r=%.2f)", string(comp_label), var_str, metrics.channel_corr[comp_idx])
    else
        return @sprintf("%s (%s)", string(comp_label), var_str)
    end
end

function _ica_ranking_summary_text(ica::InfoIca, metrics; top_k::Int = 3)
    k = min(top_k, length(ica.ica_label))
    k == 0 && return ""
    var_parts = [@sprintf("%s (%.1f%%)", string(ica.ica_label[c]), metrics.variance_pct[c]) for c in metrics.var_rank[1:k]]
    sections = ["Top Variance: " * join(var_parts, ", ")]
    if !isnothing(metrics.channel) && !all(isnan, metrics.channel_corr)
        c_parts = [@sprintf("%s (%.2f)", string(ica.ica_label[c]), metrics.channel_corr[c]) for c in metrics.corr_rank[1:k]]
        push!(sections, "Top |r| with $(metrics.channel): " * join(c_parts, ", "))
    end
    return join(sections, "   |   ")
end

function _show_ica_menu(state, ax, ica; display_plot::Bool = true)
    all_components = ica.ica_label
    n_comps = length(all_components)
    orig_dat = state.data.original
    avail_channels = _ica_correlatable_channels(orig_dat)
    acts_cache = Ref{Union{Nothing,Matrix{Float64}}}(nothing)
    metrics_ref = Ref(_ica_component_metrics(orig_dat, ica))
    n_cols = 5

    menu_fig = _eegfun_figure(size = (960, 680))
    header_area = menu_fig[1, 1] = GridLayout()
    Label(header_area[1, 1], "Select ICA Components to Remove", fontsize = 18, halign = :center)
    summary_obs = Observable(_ica_ranking_summary_text(ica, metrics_ref[]))
    Label(header_area[2, 1], summary_obs, fontsize = 13, color = :gray30, halign = :center)

    ctrl_area = header_area[3, 1] = GridLayout()
    Label(ctrl_area[1, 1], "Correlate with:", fontsize = 13, halign = :right)
    menu_labels = vcat(["None"], string.(avail_channels))
    menu_vals = Union{Nothing,Symbol}[nothing; avail_channels]
    ch_menu = Menu(ctrl_area[1, 2], options = zip(menu_labels, menu_vals), default = "None", width = 160)

    Label(ctrl_area[1, 3], "Order by:", fontsize = 13, halign = :right)
    sort_options = Tuple{String,Symbol}[("IC # (Variance %)", :id), ("Channel |r|", :corr)]
    sort_btns = [Button(ctrl_area[1, 3+j], label = lbl, height = 26, fontsize = 12) for (j, (lbl, _)) in enumerate(sort_options)]
    sort_mode = Observable(:id)
    colgap!(ctrl_area, 10)

    selected_set = Set{Int}(state.ica.removed_components)
    scroll_area = menu_fig[2, 1] = GridLayout()
    current_order = collect(1:n_comps)
    updating_grid = Ref(false)
    cbs = Checkbox[]
    lbls = Label[]

    for k = 1:n_comps
        row = ((k - 1) ÷ n_cols) + 1
        col_base = ((k - 1) % n_cols) * 2
        c_idx = current_order[k]
        cb = Checkbox(scroll_area[row, col_base+1], checked = c_idx in selected_set, width = 16, height = 16)
        lbl = Label(
            scroll_area[row, col_base+2],
            _format_ica_component_label(all_components[c_idx], c_idx, metrics_ref[]),
            fontsize = 13,
            halign = :left,
        )
        push!(cbs, cb)
        push!(lbls, lbl)
        slot = k
        on(cb.checked) do active
            updating_grid[] && return
            comp_id = current_order[slot]
            if active
                push!(selected_set, comp_id)
            else
                delete!(selected_set, comp_id)
            end
        end
    end
    actual_cols = min(n_cols, n_comps)
    for c = 1:actual_cols
        col_base = (c - 1) * 2
        colgap!(scroll_area, col_base + 1, 5)
        c < actual_cols && colgap!(scroll_area, col_base + 2, 20)
    end
    rowgap!(scroll_area, 5)

    function _refresh_ica_grid!()
        m = metrics_ref[]
        summary_obs[] = _ica_ranking_summary_text(ica, m)
        order = (sort_mode[] == :corr && !all(isnan, m.channel_corr)) ? m.corr_rank : (1:n_comps)
        updating_grid[] = true
        for k = 1:n_comps
            c_idx = order[k]
            current_order[k] = c_idx
            lbls[k].text = _format_ica_component_label(all_components[c_idx], c_idx, m)
            cbs[k].checked[] = c_idx in selected_set
        end
        updating_grid[] = false
    end

    _radio_buttons!(sort_btns, sort_mode, [opt[2] for opt in sort_options])
    on(_ -> _refresh_ica_grid!(), sort_mode)

    on(ch_menu.selection) do sel_ch
        if !isnothing(sel_ch) && isnothing(acts_cache[])
            acts_cache[] = _ica_activation_matrix(orig_dat, ica)
        end
        metrics_ref[] = _ica_component_metrics(orig_dat, ica; channel = sel_ch, acts = acts_cache[])
        if isnothing(sel_ch)
            sort_mode[] == :id ? _refresh_ica_grid!() : (sort_mode[] = :id)
        else
            sort_mode[] == :corr ? _refresh_ica_grid!() : (sort_mode[] = :corr)
        end
        for (b, v) in zip(sort_btns, [opt[2] for opt in sort_options])
            b.buttoncolor[] = v == sort_mode[] ? :lightblue : :white
        end
    end

    group_area = menu_fig[3, 1] = GridLayout()
    btn_reset = Button(group_area[1, 1], label = "None (Reset)")
    on(btn_reset.clicks) do _
        empty!(selected_set)
        updating_grid[] = true
        for cb in cbs
            cb.checked[] = false
        end
        updating_grid[] = false
    end
    btn_high_corr = Button(group_area[1, 2], label = "Select |r| ≥ 0.5")
    on(btn_high_corr.clicks) do _
        m = metrics_ref[]
        all(isnan, m.channel_corr) && return
        empty!(selected_set)
        for c = 1:n_comps
            !isnan(m.channel_corr[c]) && m.channel_corr[c] >= 0.5 && push!(selected_set, c)
        end
        updating_grid[] = true
        for (k, cb) in enumerate(cbs)
            cb.checked[] = current_order[k] in selected_set
        end
        updating_grid[] = false
    end

    action_area = menu_fig[4, 1] = GridLayout()
    btn_apply = Button(action_area[1, 1], label = "Apply", width = 200)
    on(btn_apply.clicks) do _
        _clear_and_save_limits!(ax, state)
        empty!(state.ica.removed_components)
        append!(state.ica.removed_components, sort!(collect(selected_set)))
        _rebuild_current!(state)
        _draw(ax, state)
    end

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
        on(active -> (obs[] = active), cb.checked)
        colgap!(display_area, col, 5)
        j < length(display_items) && colgap!(display_area, col + 1, 30)
    end

    display_plot && _display_in_screen(menu_fig; size = (960, 680))
    return menu_fig
end

function _create_ica_menu(fig, ax, state, ica)
    if !haskey(state.plot_kwargs, :show_cleaned_ica) || !(state.plot_kwargs[:show_cleaned_ica] isa Observable)
        state.plot_kwargs[:show_cleaned_ica] = Observable(false)
    end
    state.view.show_original_ica[] = true
    state.view.show_subtracted_ica[] = false
    return _popup_button(fig, "ICA Components", () -> _show_ica_menu(state, ax, ica))
end

function _create_epoch_menu(fig, ax, state)
    n_ep = n_epochs(state.data.original)
    slider_epoch = Slider(fig[2, 1], range = 1:n_ep, startvalue = state.data.current_epoch[], snap = true)
    label = Label(fig, @lift("Epoch: $($(slider_epoch.value))/$n_ep"), fontsize = 22, halign = :left, tellwidth = false)

    updating_from_keyboard = Ref(false)

    on(slider_epoch.value) do epoch_num
        if !updating_from_keyboard[]
            state.data.current_epoch[] = epoch_num
            ax.title = "Epoch $(epoch_num)/$n_ep"
            _update_markers!(ax, state)
            !isempty(state.extra_channel.channels) && _draw_extra_channel!(ax, state)
        end
    end

    on(state.data.current_epoch) do epoch_num
        updating_from_keyboard[] = true
        set_close_to!(slider_epoch, epoch_num)
        updating_from_keyboard[] = false
    end

    return hcat(slider_epoch, label)
end

const _RegionTriggerEvent = NamedTuple{(:sample, :code, :time),Tuple{Int,Int,Float64}}

function _scan_region_triggers(trig_col::AbstractVector, time_col::AbstractVector{<:Real}, x_min::Float64, x_max::Float64, active_triggers)
    events = _RegionTriggerEvent[]
    @inbounds for i in eachindex(trig_col, time_col)
        t = Float64(time_col[i])
        (t < x_min || t > x_max) && continue
        tr = trig_col[i]
        tr == 0 && continue
        (!isnothing(active_triggers) && !(tr in active_triggers)) && continue
        (i == 1 || tr != trig_col[i-1]) && push!(events, (sample = i, code = Int(tr), time = t))
    end
    return events
end

function _find_region_triggers(state::ContinuousDataBrowserState, clicked_region_idx = nothing)
    bounds = _get_selected_region_bounds(state, clicked_region_idx)
    isnothing(bounds) && return _RegionTriggerEvent[]
    df = state.data.current[].data
    !hasproperty(df, :trigger) && return _RegionTriggerEvent[]
    return _scan_region_triggers(df.trigger, df.time, bounds[1], bounds[2], get(state.plot_kwargs, :active_triggers, nothing))
end

_find_region_triggers(state, clicked_region_idx = nothing) = _RegionTriggerEvent[]

function _plot_all_trigger_epochs(state::ContinuousDataBrowserState, code::Integer, x_min::Real, x_max::Real, t_trig::Real; kwargs...)
    epoch_cond = EpochCondition(name = "Trigger_$code", trigger_sequences = [[code]])
    try
        all_epochs = extract_epochs(state.data.current[], 1, epoch_cond, (x_min - t_trig, x_max - t_trig))
        return plot_epochs(all_epochs; channel_selection = channels(_get_browser_target_channels(state)), kwargs...)
    catch e
        @minimal_warning "Could not extract all epochs for trigger $code: $e"
        return nothing
    end
end

function _plot_trigger_erp(state::ContinuousDataBrowserState, code::Integer, x_min::Real, x_max::Real, t_trig::Real; kwargs...)
    epoch_cond = EpochCondition(name = "Trigger_$code", trigger_sequences = [[code]])
    try
        all_epochs = extract_epochs(state.data.current[], 1, epoch_cond, (x_min - t_trig, x_max - t_trig))
        erp = average_epochs(all_epochs)
        return plot_erp(erp; channel_selection = channels(_get_browser_target_channels(state)), kwargs...)
    catch e
        @minimal_warning "Could not compute ERP for trigger $code: $e"
        return nothing
    end
end

function _show_additional_menu(state, clicked_region_idx = nothing)
    menu_items = Tuple{String,Function}[]

    if has_valid_coordinates(state.data.original.layout)
        push!(
            menu_items,
            ("Topoplot", () -> begin
                selected_data = _subset_selected_data(state, clicked_region_idx; all_visible_channels = true)
                !isnothing(selected_data) && plot_topography(selected_data; scrubber = true)
            end),
        )
        push!(
            menu_items,
            ("Topoplot 3D", () -> begin
                selected_data = _subset_selected_data(state, clicked_region_idx; all_visible_channels = true)
                !isnothing(selected_data) && plot_topography_3d(selected_data; scrubber = true)
            end),
        )
    end
    push!(menu_items, ("Spectrum", () -> begin
        selected_data = _subset_selected_data(state, clicked_region_idx)
        !isnothing(selected_data) && plot_channel_spectrum(selected_data)
    end))

    if state.data isa ContinuousDataState && get(state.marker_visibility, :trigger, false)
        trigger_events = _find_region_triggers(state, clicked_region_idx)
        bounds = _get_selected_region_bounds(state, clicked_region_idx)
        if !isempty(trigger_events) && !isnothing(bounds)
            x_min, x_max = bounds
            for code in unique(ev.code for ev in trigger_events)
                ev = trigger_events[findfirst(e -> e.code == code, trigger_events)]
                push!(menu_items, ("Plot Epochs (Trigger $code)", () -> _plot_all_trigger_epochs(state, code, x_min, x_max, ev.time)))
                push!(menu_items, ("Plot ERP (Trigger $code)", () -> _plot_trigger_erp(state, code, x_min, x_max, ev.time)))
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

    fig_size = (340, max(160, 50 * length(menu_items)))
    menu_fig = _eegfun_figure(size = fig_size)
    for (idx, (label_str, action_fn)) in enumerate(menu_items)
        btn = Button(menu_fig[idx, 1], label = label_str, tellwidth = false)
        on(_ -> action_fn(), btn.clicks)
    end

    _display_in_screen(menu_fig; size = fig_size)
end

function _channel_repair_menu(state, selected_channels, ax)
    all_channels = state.channels.labels
    menu_fig = _eegfun_figure(size = (900, 700))
    Label(menu_fig[1, 1], "Channel Repair Interface", fontsize = 18, halign = :center)

    scroll_area = menu_fig[2, 1] = GridLayout()
    channel_checkboxes = _build_checkbox_grid!(
        scroll_area,
        all_channels,
        9,
        (ch, _) -> ch in selected_channels;
        fontsize = 14,
        label_fn = ch -> _is_channel_repaired(state, ch) ? "$(ch) ✓" : string(ch),
    )

    Label(menu_fig[3, 1], "Repair Method:", fontsize = 14)
    method_area = menu_fig[4, 1] = GridLayout()
    method_buttons = [
        Button(method_area[1, 1], label = "Neighbor Interpolation", width = 200),
        Button(method_area[1, 2], label = "Spherical Spline", width = 200),
    ]

    action_area = menu_fig[5, 1] = GridLayout()
    action_buttons =
        [Button(action_area[1, 1], label = "Apply Repair", width = 200), Button(action_area[1, 2], label = "Undo Last Repair", width = 200)]

    selected_method = Observable(:neighbor_interpolation)
    _radio_buttons!(method_buttons, selected_method, [:neighbor_interpolation, :spherical_spline])

    on(action_buttons[1].clicks) do _
        channels_to_repair = all_channels[findall(cb -> cb.checked[], channel_checkboxes)]
        if !isempty(channels_to_repair)
            _repair_selected_channels!(state, channels_to_repair, selected_method[], ax)
        else
            @info "No channels selected for repair"
        end
    end

    on(action_buttons[2].clicks) do _
        if !isempty(state.channel_repair_history)
            _undo_last_repair!(state, ax)
        else
            @info "No repairs to undo"
        end
    end

    _display_in_screen(menu_fig)
end

function _repair_selected_channels!(state, selected_channels, method, ax)
    already_repaired = Set{Symbol}()
    for (repaired_channels, _, _) in state.channel_repair_history, ch in repaired_channels
        ch in selected_channels && push!(already_repaired, ch)
    end

    if !isempty(already_repaired)
        @info "Channels $(join(string.(collect(already_repaired)), ", ")) have already been repaired. Please undo first or select different channels."
        return
    end

    push!(state.channel_repair_history, (selected_channels, method, nothing))
    _rebuild_current!(state)
    _clear_and_save_limits!(ax, state)
    _draw(ax, state)

    @info "Successfully repaired channels: $(join(string.(selected_channels), ", ")) using $method"
    @info "Total repairs in history: $(length(state.channel_repair_history))"
end

function _undo_last_repair!(state, ax)
    if isempty(state.channel_repair_history)
        @info "No repairs to undo"
        return
    end

    channels, method, _ = pop!(state.channel_repair_history)
    _rebuild_current!(state)
    _clear_and_save_limits!(ax, state)
    _draw(ax, state)

    @info "Undid repair of channels: $(join(string.(channels), ", ")) (was $method)"
    @info "Remaining repairs in history: $(length(state.channel_repair_history))"
end

function _create_extreme_slider(fig, state)
    slider_extreme = Slider(fig, range = 0:5:100, startvalue = 0, tellwidth = false)
    on(x -> (state.view.crit_val[] = x), slider_extreme.value)
    lbl = Label(fig, @lift(string("Extreme: ", $(slider_extreme.value), " μV")), halign = :left, fontsize = 14)
    return (lbl, slider_extreme)
end

function _set_window_duration!(state::ContinuousDataBrowserState, duration_sec::Real)
    n_total = nrow(state.data.current[].data)
    n_total <= 1 && return
    sr = state.data.original.sample_rate > 0 ? Float64(state.data.original.sample_rate) : 1.0
    w = clamp(round(Int, Float64(duration_sec) * sr), 2, n_total)
    start_idx = clamp(first(state.view.xrange[]), 1, max(1, n_total - w + 1))
    state.view.xrange[] = start_idx:(start_idx+w-1)
end

function _create_sliders(fig, state::ContinuousDataBrowserState, dat)
    extreme_widgets = state.ica.is_active ? nothing : _create_extreme_slider(fig, state)

    sr = dat.sample_rate > 0 ? Float64(dat.sample_rate) : 1.0
    n_total = nrow(state.data.current[].data)
    total_sec = n_total / sr

    sliders_grid = fig[2, 1] = GridLayout(tellwidth = true, tellheight = true)

    Label(sliders_grid[1, 1], "Position:", fontsize = 13, halign = :left, width = 62)
    step_samples = max(1, round(Int, 0.1 * sr))
    slider_x = Slider(sliders_grid[1, 2], range = 1:step_samples:n_total, startvalue = first(state.view.xrange[]), snap = true)
    pos_lbl_text = @lift begin
        rng = $(state.view.xrange)
        tvec = state.data.current[].data.time
        t0 = round(tvec[clamp(first(rng), 1, length(tvec))], digits = 1)
        t1 = round(tvec[clamp(last(rng), 1, length(tvec))], digits = 1)
        t_end = round(tvec[end], digits = 1)
        "$(t0)–$(t1) / $(t_end) s"
    end
    Label(sliders_grid[1, 3], pos_lbl_text, fontsize = 13, halign = :left, tellwidth = false)

    Label(sliders_grid[2, 1], "Window:", fontsize = 13, halign = :left, width = 62)
    min_sec = round(max(0.1, min(1.0, total_sec)), digits = 1)
    max_sec = round(max(min_sec, min(60.0, total_sec)), digits = 1)
    step_sec = max_sec <= 5.0 ? 0.1 : 0.5
    sec_range = min_sec < max_sec ? (min_sec:step_sec:max_sec) : (min_sec:0.1:(min_sec+0.1))
    init_sec = clamp(round(length(state.view.xrange[]) / sr, digits = 1), first(sec_range), last(sec_range))

    slider_range = Slider(sliders_grid[2, 2], range = sec_range, startvalue = init_sec, snap = true)
    win_right = sliders_grid[2, 3] = GridLayout(halign = :left)
    win_lbl_text = @lift "$(round(length($(state.view.xrange)) / sr, digits = 1)) s"
    Label(win_right[1, 1], win_lbl_text, fontsize = 13, halign = :left, width = 48)

    preset_secs = [s for s in (1.0, 5.0, 10.0, 30.0) if s <= max_sec + 1e-6]
    for (j, psec) in enumerate(preset_secs)
        btn = Button(win_right[1, 1+j], label = "$(Int(psec))s", width = 38, height = 24, fontsize = 12)
        on(_ -> _set_window_duration!(state, psec), btn.clicks)
    end
    colgap!(win_right, 6)
    colgap!(sliders_grid, 10)
    rowgap!(sliders_grid, 6)

    updating_from_key = Ref(false)

    function _update_window(start::Int, dur_sec::Real)
        n_cur = nrow(state.data.current[].data)
        n_cur <= 1 && return
        w = clamp(round(Int, Float64(dur_sec) * sr), 2, n_cur)
        s = clamp(start, 1, max(1, n_cur - w + 1))
        state.view.xrange[] = s:(s+w-1)
    end

    on(dur_sec -> (updating_from_key[] || _update_window(slider_x.value[], dur_sec)), slider_range.value)
    on(start -> (updating_from_key[] || _update_window(start, slider_range.value[])), slider_x.value)

    on(state.view.xrange) do rng
        updating_from_key[] = true
        set_close_to!(slider_x, first(rng))
        set_close_to!(slider_range, length(rng) / sr)
        updating_from_key[] = false
    end

    return extreme_widgets
end

_create_sliders(fig, state::EpochedDataBrowserState, dat) = state.ica.is_active ? nothing : _create_extreme_slider(fig, state)

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

_create_extra_channel_menu(fig, ax, state, dat) = _popup_button(fig, "Extra Channels", () -> _show_extra_channel_menu(state, ax, dat))

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
    control_panel = GridLayout(fig[1, 2], tellwidth = true, tellheight = false)
    colsize!(control_panel, 1, Fixed(42))
    rowgap!(control_panel, 8)

    cur_row = 1
    for i = 1:size(toggles, 1)
        control_panel[cur_row, 1] = toggles[i, 1]
        control_panel[cur_row, 2] = toggles[i, 2]
        cur_row += 1
    end

    extreme_widgets = _create_sliders(fig, state, dat)
    if !isnothing(extreme_widgets)
        lbl_ext, slider_ext = extreme_widgets
        control_panel[cur_row, 1:2] = lbl_ext
        control_panel[cur_row+1, 1:2] = slider_ext
        cur_row += 2
    end

    for widget in (labels_menu, reference_menu, ica_menu, extra_menu, epoch_menu)
        if !isnothing(widget)
            control_panel[cur_row, 1:2] = widget
            cur_row += 1
        end
    end

    colsize!(fig.layout, 2, Fixed(220))
    colgap!(fig.layout, 15)
    rowsize!(fig.layout, 1, Relative(1.0))
end

_get_epoch_menu(fig, ax, state::ContinuousDataBrowserState) = nothing
_get_epoch_menu(fig, ax, state::EpochedDataBrowserState) = _create_epoch_menu(fig, ax, state)

# =============================================================================
# Navigation & Interaction
# =============================================================================

const KEYBOARD_ACTIONS = Dict(
    Keyboard.left => :left,
    Keyboard.right => :right,
    Keyboard.up => :up,
    Keyboard.down => :down,
    Keyboard.page_up => :page_back,
    Keyboard.page_down => :page_forward,
)

function _handle_navigation!(ax, state::DataBrowserState{<:AbstractDataState}, action::Symbol)
    if action == :up
        _ymore!(ax, state)
    elseif action == :down
        _yless!(ax, state)
    elseif action == :left
        _handle_left_navigation(ax, state)
    elseif action == :right
        _handle_right_navigation(ax, state)
    elseif action == :page_back
        _handle_page_back(ax, state)
    elseif action == :page_forward
        _handle_page_forward(ax, state)
    end
end

_handle_left_navigation(ax, state::ContinuousDataBrowserState) = _xback!(ax, state)
_handle_left_navigation(ax, state::EpochedDataBrowserState) = _step_epoch_backward(ax, state)
_handle_right_navigation(ax, state::ContinuousDataBrowserState) = _xforward!(ax, state)
_handle_right_navigation(ax, state::EpochedDataBrowserState) = _step_epoch_forward(ax, state)
_handle_page_back(ax, state::ContinuousDataBrowserState) = _xback!(ax, state, length(state.view.xrange.val))
_handle_page_back(ax, state::EpochedDataBrowserState) = _step_epoch!(ax, state, -5)
_handle_page_forward(ax, state::ContinuousDataBrowserState) = _xforward!(ax, state, length(state.view.xrange.val))
_handle_page_forward(ax, state::EpochedDataBrowserState) = _step_epoch!(ax, state, 5)

function _xback!(ax, state::ContinuousDataBrowserState, step::Int = state.plot_kwargs[:scroll_step])
    current_range = state.view.xrange.val
    current_range[1] <= 1 && return
    shift = min(step, current_range[1] - 1)
    state.view.xrange[] = current_range .- shift
    _restore_axis_limits!(ax, state)
end

function _xforward!(ax, state::ContinuousDataBrowserState, step::Int = state.plot_kwargs[:scroll_step])
    current_range = state.view.xrange.val
    max_idx = nrow(state.data.current[].data)
    current_range[end] >= max_idx && return
    shift = min(step, max_idx - current_range[end])
    state.view.xrange[] = current_range .+ shift
    _restore_axis_limits!(ax, state)
end

_step_epoch_backward(ax, state::EpochedDataBrowserState) = _step_epoch!(ax, state, -1)
_step_epoch_forward(ax, state::EpochedDataBrowserState) = _step_epoch!(ax, state, 1)

function _step_epoch!(ax, state::EpochedDataBrowserState, direction::Int)
    current = state.data.current_epoch[]
    total = n_epochs(state.data.original)
    new_epoch = clamp(current + direction, 1, total)
    new_epoch == current && return
    state.data.current_epoch[] = new_epoch
    ax.title = "Epoch $(new_epoch)/$total"
    _update_markers!(ax, state)
    !isempty(state.extra_channel.channels) && _draw_extra_channel!(ax, state)
end

_yless!(ax, state) = _yzoom!(ax, state, 1.2)
_ymore!(ax, state) = _yzoom!(ax, state, 0.8)

function _yzoom!(ax, state, factor::Float64)
    zoom_step = state.plot_kwargs[:zoom_step]
    if state.view.butterfly[]
        y_min, y_max = state.view.yrange.val[1], state.view.yrange.val[end]
        if factor > 1.0
            (y_min + zoom_step >= 0 || y_max - zoom_step <= 0) && return
            state.view.yrange[] = (y_min+zoom_step):(y_max-zoom_step)
        else
            state.view.yrange[] = (y_min-zoom_step):(y_max+zoom_step)
        end
        ylims!(ax, state.view.yrange.val[1], state.view.yrange.val[end])
    else
        state.view.amplitude_scale[] = state.view.amplitude_scale[] * factor
    end
end

function _is_mouse_in_browser_axis(ax, pos)
    bbox = ax.layoutobservables.computedbbox[]
    return bbox.origin[1] <= pos[1] <= (bbox.origin[1] + bbox.widths[1]) && bbox.origin[2] <= pos[2] <= (bbox.origin[2] + bbox.widths[2])
end

function _find_clicked_region(state, mouse_x)
    for (i, (start_time, end_time)) in enumerate(state.selection.selected_regions[])
        if min(start_time, end_time) <= mouse_x <= max(start_time, end_time)
            return i
        end
    end
    return nothing
end

function _remove_region_from_selection!(ax, state, region_idx)
    regions = state.selection.selected_regions[]
    if 1 <= region_idx <= length(regions)
        delete!(ax.scene, state.selection.region_plots[region_idx])
        deleteat!(state.selection.region_plots, region_idx)
        deleteat!(regions, region_idx)
        state.selection.selected_regions[] = regions
        _update_analysis_settings!(state)
    end
end

function _start_selection!(ax, state, mouse_x)
    state.selection.active[] = true
    state.selection.bounds[] = (mouse_x, mouse_x)
    _update_x_region_selection!(ax, state, mouse_x, mouse_x)
end

function _finish_selection!(ax, state, mouse_x)
    t1 = state.selection.bounds[][1]
    t2 = mouse_x
    state.selection.active[] = false
    _clear_x_region_selection!(state)
    abs(t2 - t1) < 1e-4 && return
    state.selection.bounds[] = (t1, t2)
    _add_region_to_selection!(ax, state, t1, t2)
end

function _handle_mouse_events!(ax, state)
    shift_pressed = Ref(false)
    ctrl_pressed = Ref(false)

    on(events(ax).keyboardbutton) do key_event
        if key_event.key in (Keyboard.left_shift, Keyboard.right_shift)
            shift_pressed[] = (key_event.action != Keyboard.release)
        elseif key_event.key in (Keyboard.left_control, Keyboard.right_control)
            ctrl_pressed[] = (key_event.action != Keyboard.release)
        end
    end

    on(events(ax).mousebutton) do event
        _is_mouse_in_browser_axis(ax, events(ax).mouseposition[]) || return
        mouse_x = mouseposition(ax)[1]

        if event.button == Mouse.left
            if event.action == Mouse.press
                if ctrl_pressed[]
                    clicked_channel = _find_closest_browser_channel(ax, state, mouse_x, mouseposition(ax)[2])
                    if !isnothing(clicked_channel)
                        _toggle_channel_selection!(ax, state, clicked_channel)
                        return
                    end
                elseif shift_pressed[]
                    _handle_left_click!(ax, state, event, mouse_x)
                else
                    clicked_region_idx = _find_clicked_region(state, mouse_x)
                    if !isnothing(clicked_region_idx)
                        _remove_region_from_selection!(ax, state, clicked_region_idx)
                        return
                    end
                end
            elseif event.action == Mouse.release
                (state.selection.active[] || shift_pressed[]) && _handle_left_click!(ax, state, event, mouse_x)
            end
        elseif event.button == Mouse.right && event.action == Mouse.press
            _handle_right_click!(ax, state, mouse_x)
        end
    end

    on(events(ax).mouseposition) do _
        if state.selection.active[]
            _update_x_region_selection!(ax, state, state.selection.bounds[][1], mouseposition(ax)[1])
        end
    end
end

function _handle_left_click!(ax, state, event, mouse_x)
    if event.action == Mouse.press
        clicked_region_idx = _find_clicked_region(state, mouse_x)
        if !isnothing(clicked_region_idx)
            _remove_region_from_selection!(ax, state, clicked_region_idx)
        else
            _start_selection!(ax, state, mouse_x)
        end
    elseif event.action == Mouse.release && state.selection.active[]
        _finish_selection!(ax, state, mouse_x)
    end
end

function _handle_right_click!(ax, state, mouse_x)
    clicked_region_idx = _find_clicked_region(state, mouse_x)
    !isnothing(clicked_region_idx) && _show_additional_menu(state, clicked_region_idx)
end

function _show_channel_repair_menu(state::DataBrowserState{<:ContinuousDataState}, ax)
    if !has_valid_coordinates(state.data.original.layout)
        @minimal_warning "Channel repair requires electrode layout coordinates. Please load a layout file (e.g. using EegFun.read_layout)."
        return
    end
    _channel_repair_menu(state, _get_selected_channels(state), ax)
end

_show_channel_repair_menu(::DataBrowserState{<:EpochedDataState}, ax) = @info "Channel repair is only available for continuous data"

_get_selected_channels(state) = state.channels.labels[state.channels.selected]

_get_browser_target_channels(state) =
    any(state.channels.selected) ? state.channels.labels[state.channels.selected] : state.channels.labels[state.channels.visible]

function _find_closest_browser_channel(ax, state, mouse_x, mouse_y)
    current_data = _get_current_data(state.data)
    time_data = current_data.time
    x_idx = searchsortedfirst(time_data, mouse_x)
    if x_idx > 1 && x_idx <= length(time_data) && abs(time_data[x_idx-1] - mouse_x) < abs(time_data[x_idx] - mouse_x)
        x_idx -= 1
    end
    x_idx = clamp(x_idx, 1, length(time_data))

    min_distance = Inf
    closest_channel = nothing
    scale = state.view.amplitude_scale[]
    rem_dc = state.view.remove_dc[]
    win = state.view.xrange[]
    for (idx, visible) in enumerate(state.channels.visible)
        visible || continue
        col = state.channels.labels[idx]
        dc = rem_dc ? mean(@view current_data[!, col][win]) : 0.0
        channel_y = (current_data[x_idx, col] - dc) * scale + state.view.offset[idx]
        distance = abs(mouse_y - channel_y)
        if distance < min_distance
            min_distance = distance
            closest_channel = idx
        end
    end

    visible_indices = findall(state.channels.visible)
    tol = if length(visible_indices) >= 2
        step = abs(state.view.offset[visible_indices[end]] - state.view.offset[visible_indices[1]]) / (length(visible_indices) - 1)
        max(step * 0.9, 50.0)
    else
        100.0
    end

    return min_distance <= tol ? closest_channel : nothing
end

function _apply_channel_visual_state!(ax, state)
    sel_color = get(state.plot_kwargs, :selected_channel_color, :black)
    normal_lw = state.plot_kwargs[:channel_line_width]
    selected_lw = normal_lw * 2
    base_unselected_color = state.plot_kwargs[:unselected_channel_color]
    is_cleaned = state.ica.is_active && haskey(state.plot_kwargs, :show_cleaned_ica) && state.plot_kwargs[:show_cleaned_ica][]
    base_color = is_cleaned ? get(state.plot_kwargs, :cleaned_channel_color, :green) : base_unselected_color

    for (idx, col) in enumerate(state.channels.labels)
        !state.channels.visible[idx] && continue
        if haskey(state.channels.data_lines, col) && haskey(state.channels.data_labels, col)
            line_handle = state.channels.data_lines[col]
            label_handle = state.channels.data_labels[col]
            if state.channels.selected[idx] || _is_channel_repaired(state, col)
                line_handle.colormap[] = [sel_color, sel_color, sel_color]
                line_handle.linewidth[] = selected_lw
                label_handle.color[] = :red
            else
                line_handle.colormap[] = [base_color, base_color, :red]
                line_handle.linewidth[] = normal_lw
                label_handle.color[] = :black
            end
        end
    end
end

function _select_single_channel!(ax, state, target_col::Symbol)
    state.channels.selected .= false
    idx = findfirst(==(target_col), state.channels.labels)
    !isnothing(idx) && (state.channels.selected[idx] = true)
    _apply_channel_visual_state!(ax, state)
end

function _clear_channel_selection!(ax, state)
    state.channels.selected .= false
    _apply_channel_visual_state!(ax, state)
end

function _toggle_channel_selection!(ax, state, channel_idx)
    state.channels.selected[channel_idx] = !state.channels.selected[channel_idx]
    _apply_channel_visual_state!(ax, state)
end

function _handle_keyboard_events!(fig, ax, state)
    on(events(fig).keyboardbutton) do event
        if event.action == Keyboard.press && event.key == Keyboard.i
            _show_plot_help(:databrowser)
        elseif event.action == Keyboard.press && event.key == Keyboard.r
            _show_channel_repair_menu(state, ax)
        elseif event.action == Keyboard.press && event.key == Keyboard.c
            _clear_all_selected_regions!(ax, state)
            _clear_channel_selection!(ax, state)
        elseif event.action in (Keyboard.press, Keyboard.repeat) && haskey(KEYBOARD_ACTIONS, event.key)
            _handle_navigation!(ax, state, KEYBOARD_ACTIONS[event.key])
        end
    end
end

function _region_rect_points(ax, state, x1, x2)
    limits_y = ax.limits[][2]
    y1, y2 = isnothing(limits_y) ? (Float64(state.view.yrange[][1]), Float64(state.view.yrange[][end])) : limits_y
    return Point2f[Point2f(x1, y1), Point2f(x2, y1), Point2f(x2, y2), Point2f(x1, y2)]
end

function _update_x_region_selection!(ax, state, x1, x2)
    state.selection.rectangle[1] = _region_rect_points(ax, state, x1, x2)
    state.selection.rectangle.visible[] = true
end

function _add_region_to_selection!(ax, state, x1, x2)
    x1 > x2 && ((x1, x2) = (x2, x1))
    current_regions = state.selection.selected_regions[]
    push!(current_regions, (x1, x2))
    state.selection.selected_regions[] = current_regions

    region_plot = poly!(
        ax,
        _region_rect_points(ax, state, x1, x2),
        color = (state.plot_kwargs[:selection_color], state.plot_kwargs[:selection_alpha]),
        strokecolor = :transparent,
    )
    push!(state.selection.region_plots, region_plot)
    _update_analysis_settings!(state)
end

function _clear_x_region_selection!(state)
    # Single point avoids CairoMakie empty-polygon issues
    state.selection.rectangle[1] = [Point2f(0.0, 0.0)]
    state.selection.bounds[] = (0.0, 0.0)
    state.selection.visible[] = false
    state.selection.rectangle.visible[] = false
end

function _clear_all_selected_regions!(ax, state)
    for plot in state.selection.region_plots
        delete!(ax.scene, plot)
    end
    empty!(state.selection.region_plots)
    state.selection.selected_regions[] = Tuple{Float64,Float64}[]
    _update_analysis_settings!(state)
end

function _get_selected_regions_bool(state::DataBrowserState)
    current_data = _get_current_data(state.data)
    total_samples = nrow(current_data)
    time_data = current_data.time
    bool_vector = falses(total_samples)

    for (start_time, end_time) in state.selection.selected_regions[]
        start_idx = clamp(searchsortedfirst(time_data, start_time), 1, total_samples)
        end_idx = clamp(searchsortedlast(time_data, end_time), 1, total_samples)
        start_idx <= end_idx && (bool_vector[start_idx:end_idx] .= true)
    end

    return bool_vector
end

function _get_selected_regions_info(state::DataBrowserState)
    bool_vector = _get_selected_regions_bool(state)
    regions = state.selection.selected_regions[]
    return (bool_vector = bool_vector, regions = regions, n_samples = sum(bool_vector), n_regions = length(regions))
end

function _get_selected_region_bounds(state, clicked_region_idx = nothing)
    regions = state.selection.selected_regions[]
    if !isnothing(clicked_region_idx) && 1 <= clicked_region_idx <= length(regions)
        r = regions[clicked_region_idx]
        return (min(r[1], r[2]), max(r[1], r[2]))
    elseif !isempty(regions)
        r = regions[end]
        return (min(r[1], r[2]), max(r[1], r[2]))
    end
    return nothing
end

function _subset_selected_data(state::ContinuousDataBrowserState, clicked_region_idx = nothing; all_visible_channels::Bool = false)
    bounds = _get_selected_region_bounds(state, clicked_region_idx)
    if isnothing(bounds)
        @minimal_warning "No region selected"
        return nothing
    end
    x_min, x_max = bounds
    selected_channels = all_visible_channels ? state.channels.labels[state.channels.visible] : _get_browser_target_channels(state)
    return subset(
        state.data.current[],
        sample_selection = x -> (x.time .>= x_min) .& (x.time .<= x_max),
        channel_selection = channels(selected_channels),
    )
end

function _subset_selected_data(state::EpochedDataBrowserState, clicked_region_idx = nothing; all_visible_channels::Bool = false)
    bounds = _get_selected_region_bounds(state, clicked_region_idx)
    if isnothing(bounds)
        @minimal_warning "No region selected"
        return nothing
    end
    x_min, x_max = bounds
    selected_channels = all_visible_channels ? state.channels.labels[state.channels.visible] : _get_browser_target_channels(state)

    epoch_data = subset(
        state.data.current[],
        sample_selection = epoch_df -> (epoch_df.time .>= x_min) .& (epoch_df.time .<= x_max),
        channel_selection = channels(selected_channels),
        epoch_selection = epochs([state.data.current_epoch[]]),
    )
    if isempty(epoch_data)
        @minimal_warning "Selected region contains no data"
        return nothing
    end
    return epoch_to_continuous(epoch_data, 1)
end

# =============================================================================
# Filtering
# =============================================================================

_apply_filters!(state) = _rebuild_current!(state)

function _apply_hp_filter!(state, active::Bool = !state.data.filter_state.hp_active[])
    state.data.filter_state.hp_active[] = active
    _apply_filters!(state)
end

function _apply_lp_filter!(state, active::Bool = !state.data.filter_state.lp_active[])
    state.data.filter_state.lp_active[] = active
    _apply_filters!(state)
end

# =============================================================================
# Drawing
# =============================================================================

"""Add a vertical marker (trigger/EOG) to the marker list, window-scoped for performance."""
function _add_marker!(markers, ax, state, data, col; label = nothing, trial = nothing, visible = false, active_values = nothing)
    src_df = isnothing(trial) ? data : data[trial]
    col_vec = src_df[!, col]
    mask = col_vec .!= 0
    if !isnothing(active_values)
        @inbounds for i in eachindex(mask, col_vec)
            mask[i] && (mask[i] = col_vec[i] in active_values)
        end
    end
    marker_data = src_df[mask, [:time, col]]
    nrow(marker_data) == 0 && return

    label_vec = isnothing(label) ? string.(marker_data[!, col]) : fill(label, nrow(marker_data))

    marker_obs = @lift(begin
        if state isa EpochedDataBrowserState
            idxs = 1:nrow(marker_data)
        else
            win = $(state.view.xrange)
            t_start = data.time[win[1]]
            t_end = data.time[win[end]]
            times = marker_data.time
            idxs = searchsortedfirst(times, t_start):searchsortedlast(times, t_end)
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

function _disconnect_observer(o)
    if o isa Observable
        Makie.Observables.clear(o)
    elseif o isa AbstractVector
        foreach(_disconnect_observer, o)
    else
        Makie.off(o)
    end
end

_all_line_dicts(state) =
    (state.channels.data_lines, state.channels.data_labels, state.channels.original_lines, state.channels.subtracted_lines)

_is_channel_repaired(state, col) = any(col in chs for (chs, _, _) in state.channel_repair_history)

function _update_channel_offsets!(state)
    nchannels = count(state.channels.visible)
    margin = state.plot_kwargs[:channel_offset_margin]
    if nchannels > 1 && !state.view.butterfly[]
        y_max = state.view.yrange[][end] * margin
        y_min = state.view.yrange[][1] * margin
        step = (y_min - y_max) / (nchannels - 1)
        i = 0
        @inbounds for idx in eachindex(state.channels.visible)
            if state.channels.visible[idx]
                state.view.offset[idx] = y_max + i * step
                i += 1
            end
        end
    else
        state.view.offset[state.channels.visible] .= 0.0
    end
end

function _clear_axes!(ax, datas)
    for data in datas
        for value in values(data)
            delete!(ax, value)
        end
        empty!(data)
    end
end

function _set_axes!(ax, state::DataBrowserState{<:AbstractDataState})
    @lift ylims!(ax, $(state.view.yrange)[1], $(state.view.yrange)[end])
    _set_x_limits!(ax, state, state.data)
end

function _set_x_limits!(ax, state, data::ContinuousDataState)
    @lift xlims!(ax, data.current[].data.time[$(state.view.xrange)[1]], data.current[].data.time[$(state.view.xrange)[end]])
end

function _set_x_limits!(ax, state, data::EpochedDataState)
    @lift xlims!(ax, $(data.current).data[1].time[1], $(data.current).data[1].time[end])
end

function _init_markers(ax, state)
    markers = Marker[]
    data = _get_current_data(state.data)
    for (symbol, label) in ((:trigger, nothing), (:is_vEOG, "v"), (:is_hEOG, "h"))
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

function _update_markers!(ax, state)
    for marker in state.markers
        foreach(_disconnect_observer, marker.observers)
        delete!(ax, marker.line)
        delete!(ax, marker.text)
    end
    empty!(state.markers)
    state.markers = _init_markers(ax, state)
end

function _butterfly_plot!(ax, state)
    state.view.butterfly[] = !state.view.butterfly[]
    _clear_and_save_limits!(ax, state)
    _update_channel_offsets!(state)
    _draw(ax, state)
    !isempty(state.extra_channel.channels) && _draw_extra_channel!(ax, state)
end

function _window_points!(pts::Vector{Point2{Float64}}, x, y, scale, offset, remove_dc::Bool = false)
    dc = remove_dc && !isempty(y) ? mean(y) : 0.0
    @inbounds for i in eachindex(pts)
        pts[i] = Point2{Float64}(x[i], (y[i] - dc) * scale + offset)
    end
    return pts
end

function _nan_points!(pts::Vector{Point2{Float64}}, x)
    @inbounds for i in eachindex(pts)
        pts[i] = Point2{Float64}(x[i], NaN)
    end
    return pts
end

# Reuses a single `Vector{Point2{Float64}}` per line in place; only reallocates when window length changes.
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
            obs[] = buf
        end
        return nothing
    end
    refresh()
    push!(tracker, onany(refresh, xrange, deps...))
    return obs
end

_ica_original_series(data::ContinuousDataState, ep, col) = data.original.data[!, col]
_ica_original_series(data::EpochedDataState, ep, col) = data.original.data[ep][!, col]

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

function _fill_ica_points!(
    buf::Vector{Point2{Float64}},
    x,
    raw,
    state,
    ep::Int,
    win::UnitRange{Int},
    ch_row::Int,
    scale::Float64,
    offset::Float64,
    subtract_raw::Bool,
    remove_dc::Bool = false,
)
    comps = state.ica.removed_components
    if isempty(comps)
        if subtract_raw
            @inbounds for j in eachindex(buf)
                buf[j] = Point2{Float64}(x[j], offset)
            end
        else
            _window_points!(buf, x, raw, scale, offset, remove_dc)
        end
        return buf
    end
    rng, acts = _ica_activations!(state, ep, win)
    ica = state.ica.original::InfoIca
    win_off = first(win) - first(rng)
    n_comps = length(comps)
    val_sum = 0.0
    @inbounds for j in eachindex(buf)
        val = Float64(raw[j])
        r = j + win_off
        for k = 1:n_comps
            comp = comps[k]
            val += (ica.mixing[ch_row, comp] * ica.scale) * acts[r, k]
        end
        subtract_raw && (val -= raw[j])
        val_sum += val
        buf[j] = Point2{Float64}(x[j], val)
    end
    dc = (remove_dc && !subtract_raw && !isempty(buf)) ? (val_sum / length(buf)) : 0.0
    @inbounds for j in eachindex(buf)
        p = buf[j]
        buf[j] = Point2{Float64}(p[1], (p[2] - dc) * scale + offset)
    end
    return buf
end

# Keeps `color` as `Observable{Vector{Bool}}` for GLMakie shader type stability while skipping
# allocations/notifications on constant-length scrolls when `crit_val <= 0`.
function _color_mask_observable(state, epoch_obs, get_data, col::Symbol, ch_row::Int)
    xrange = state.view.xrange
    crit_obs = state.view.crit_val
    remove_dc_obs = state.view.remove_dc
    obs = Observable(fill(false, length(xrange[])))
    all_false = Ref(true)
    function refresh(_...)
        win = xrange[]
        n = length(win)
        crit = crit_obs[]
        buf = obs[]
        resized = length(buf) != n
        resized && (buf = Vector{Bool}(undef, n))
        if crit > 0.0
            cur = state.data.current[]
            ep = epoch_obs[]
            raw = get_data(cur, ep, win, col)
            rem_dc = remove_dc_obs[]
            if state.ica.is_active &&
               ch_row > 0 &&
               !state.plot_kwargs[:show_cleaned_ica][] &&
               state.view.show_original_ica[] &&
               !isempty(state.ica.removed_components)
                rng, acts = _ica_activations!(state, ep, win)
                ica = state.ica.original::InfoIca
                comps = state.ica.removed_components
                win_off = first(win) - first(rng)
                dc = 0.0
                if rem_dc && n > 0
                    s = 0.0
                    @inbounds for j = 1:n
                        val = Float64(raw[j])
                        r = j + win_off
                        for k in eachindex(comps)
                            comp = comps[k]
                            val += (ica.mixing[ch_row, comp] * ica.scale) * acts[r, k]
                        end
                        s += val
                    end
                    dc = s / n
                end
                @inbounds for j = 1:n
                    val = Float64(raw[j])
                    r = j + win_off
                    for k in eachindex(comps)
                        comp = comps[k]
                        val += (ica.mixing[ch_row, comp] * ica.scale) * acts[r, k]
                    end
                    buf[j] = abs(val - dc) >= crit
                end
            else
                dc = (rem_dc && !isempty(raw)) ? mean(raw) : 0.0
                @inbounds for j = 1:n
                    buf[j] = abs(raw[j] - dc) >= crit
                end
            end
            all_false[] = false
            obs[] = buf
        elseif resized || !all_false[]
            fill!(buf, false)
            all_false[] = true
            obs[] = buf
        end
        return nothing
    end
    refresh()
    push!(
        state.channels.draw_observers,
        onany(
            refresh,
            xrange,
            crit_obs,
            remove_dc_obs,
            state.data.current,
            epoch_obs,
            state.plot_kwargs[:show_cleaned_ica],
            state.view.show_original_ica,
        ),
    )
    return obs
end

function _draw(ax, state::DataBrowserState{<:AbstractDataState})
    get_data, get_time, get_label_y = _get_data_accessors(state.data)
    epoch_obs = state.data isa EpochedDataState ? state.data.current_epoch : Observable(1)
    track(x) = (push!(state.channels.draw_observers, x); x)

    visible_time_obs = track(@lift(get_time($(state.data.current), $(epoch_obs), $(state.view.xrange))))
    time_start_obs = track(@lift($(visible_time_obs)[1]))

    ch_idx_map = state.ica.is_active ? Dict(lbl => i for (i, lbl) in enumerate(state.ica.original.layout.data.label)) : Dict{Symbol,Int}()
    ghost_orig_visible =
        state.ica.is_active ? track(@lift($(state.view.show_original_ica) && $(state.plot_kwargs[:show_cleaned_ica]))) : nothing
    data_lines_visible =
        state.ica.is_active ? track(@lift($(state.plot_kwargs[:show_cleaned_ica]) || $(state.view.show_original_ica))) : true

    any_selected = any(state.channels.selected)
    unsel_color = state.plot_kwargs[:unselected_channel_color]
    cleaned_color = get(state.plot_kwargs, :cleaned_channel_color, :green)
    sel_color = state.plot_kwargs[:selected_channel_color]
    normal_lw = state.plot_kwargs[:channel_line_width]
    line_alpha = state.plot_kwargs[:channel_line_alpha]

    for (idx, visible) in enumerate(state.channels.visible)
        col = state.channels.labels[idx]
        if visible
            is_selected = state.channels.selected[idx]
            ch_row = get(ch_idx_map, col, 0)

            line_color = _color_mask_observable(state, epoch_obs, get_data, col, ch_row)

            channel_points_obs = _points_observable(
                state.channels.draw_observers,
                state.view.xrange,
                state.data.current,
                epoch_obs,
                state.view.amplitude_scale,
                state.view.remove_dc,
                state.plot_kwargs[:show_cleaned_ica],
                state.view.show_original_ica,
            ) do buf, _
                win = state.view.xrange[]
                ep = epoch_obs[]
                cur = state.data.current[]
                raw = get_data(cur, ep, win, col)
                t_view = get_time(cur, ep, win)
                scale = state.view.amplitude_scale[]
                off = state.view.offset[idx]
                rem_dc = state.view.remove_dc[]
                if state.ica.is_active && ch_row > 0 && !state.plot_kwargs[:show_cleaned_ica][] && state.view.show_original_ica[]
                    _fill_ica_points!(buf, t_view, raw, state, ep, win, ch_row, scale, off, false, rem_dc)
                else
                    _window_points!(buf, t_view, raw, scale, off, rem_dc)
                end
                return true
            end

            if state.ica.is_active && ch_row > 0
                original_nan = Ref(false)
                original_points_obs = _points_observable(
                    state.channels.draw_observers,
                    state.view.xrange,
                    state.data.current,
                    epoch_obs,
                    state.view.amplitude_scale,
                    state.view.remove_dc,
                    state.view.show_original_ica,
                ) do buf, resized
                    win = state.view.xrange[]
                    ep = epoch_obs[]
                    cur = state.data.current[]
                    if state.view.show_original_ica[]
                        _fill_ica_points!(
                            buf,
                            get_time(cur, ep, win),
                            get_data(cur, ep, win, col),
                            state,
                            ep,
                            win,
                            ch_row,
                            state.view.amplitude_scale[],
                            state.view.offset[idx],
                            false,
                            state.view.remove_dc[],
                        )
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
                        _fill_ica_points!(
                            buf,
                            get_time(cur, ep, win),
                            get_data(cur, ep, win, col),
                            state,
                            ep,
                            win,
                            ch_row,
                            state.view.amplitude_scale[],
                            state.view.offset[idx],
                            true,
                        )
                        subtracted_nan[] = false
                        return true
                    elseif resized || !subtracted_nan[]
                        _nan_points!(buf, get_time(cur, ep, win))
                        subtracted_nan[] = true
                        return true
                    end
                    return false
                end

                _create_line!(
                    state.channels.original_lines,
                    col,
                    ax,
                    original_points_obs,
                    :grey,
                    [:grey],
                    1,
                    0.5;
                    visible = ghost_orig_visible,
                )
                _create_line!(
                    state.channels.subtracted_lines,
                    col,
                    ax,
                    subtracted_points_obs,
                    :red,
                    [:red],
                    1,
                    0.5;
                    visible = state.view.show_subtracted_ica,
                )
            end

            if is_selected || _is_channel_repaired(state, col)
                line_colormap = [sel_color, sel_color, sel_color]
                line_width = normal_lw * 2
            elseif state.ica.is_active
                line_colormap = track(
                    @lift(
                        $(state.plot_kwargs[:show_cleaned_ica]) ? [cleaned_color, cleaned_color, :red] : [unsel_color, unsel_color, :red]
                    ),
                )
                line_width = normal_lw
            else
                line_colormap = [unsel_color, unsel_color, :red]
                line_width = normal_lw
            end

            _create_line!(
                state.channels.data_lines,
                col,
                ax,
                channel_points_obs,
                line_color,
                line_colormap,
                line_width,
                line_alpha;
                visible = data_lines_visible,
            )

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
                            $(state.view.remove_dc) ? mean(get_data($(state.data.current), $(epoch_obs), $(state.view.xrange), $col)) : 0.0,
                        )
                    ),
                )
                _create_label!(state.channels.data_labels, col, ax, time_start_obs, label_y_obs, is_selected, any_selected)
            else
                _hide_channel_label!(state.channels.data_labels, col)
            end
        else
            _hide_channel_objects!(state.channels, col)
        end
    end

    _restore_axis_limits!(ax, state)
end

function _disconnect_draw_observers!(state)
    foreach(_disconnect_observer, state.channels.draw_observers)
    empty!(state.channels.draw_observers)
    return nothing
end

function _clear_and_save_limits!(ax, state)
    state.plot_kwargs[:_saved_xlims] = ax.xaxis.attributes.limits[]
    state.plot_kwargs[:_saved_ylims] = ax.yaxis.attributes.limits[]
    _clear_axes!(ax, _all_line_dicts(state))
    _disconnect_draw_observers!(state)
end

function _restore_axis_limits!(ax, state)
    if haskey(state.plot_kwargs, :_saved_xlims)
        xlims!(ax, state.plot_kwargs[:_saved_xlims]...)
        ylims!(ax, state.plot_kwargs[:_saved_ylims]...)
        delete!(state.plot_kwargs, :_saved_xlims)
        delete!(state.plot_kwargs, :_saved_ylims)
    end
end

function _get_data_accessors(::ContinuousDataState)
    get_data = (current, ep, range, col) -> @view current.data[!, col][range]
    get_time = (current, ep, range) -> @view current.data.time[range]
    get_label_y = (current, ep, col, offset, sample_idx, scale, dc = 0.0) -> (current.data[sample_idx, col] - dc) * scale + offset
    return get_data, get_time, get_label_y
end

function _get_data_accessors(::EpochedDataState)
    get_data = (current, ep, range, col) -> @view current.data[ep][!, col][range]
    get_time = (current, ep, range) -> @view current.data[ep].time[range]
    get_label_y = (current, ep, col, offset, sample_idx, scale, dc = 0.0) -> (current.data[ep][!, col][sample_idx] - dc) * scale + offset
    return get_data, get_time, get_label_y
end

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

function _create_label!(data_labels, col, ax, x_obs, y_obs, is_selected, any_selected = false)
    label_color = is_selected ? :red : (any_selected ? :gray50 : :black)
    data_labels[col] = text!(ax, x_obs, y_obs, text = String(col), align = (:left, :center), fontsize = 18, color = label_color)
end

_hide_channel_label!(data_labels, col) = haskey(data_labels, col) && (data_labels[col].visible = false)

function _hide_channel_objects!(channels, col)
    haskey(channels.data_lines, col) && (channels.data_lines[col].visible = false)
    haskey(channels.data_labels, col) && (channels.data_labels[col].visible = false)
    haskey(channels.original_lines, col) && (channels.original_lines[col].visible = false)
    haskey(channels.subtracted_lines, col) && (channels.subtracted_lines[col].visible = false)
end

function _draw_extra_channel!(ax, state::DataBrowserState{<:AbstractDataState})
    foreach(_disconnect_observer, state.extra_channel.draw_observers)
    empty!(state.extra_channel.draw_observers)
    _clear_axes!(ax, (state.extra_channel.data_lines, state.extra_channel.data_labels))

    if !isnothing(state.extra_channel.legend)
        delete!(state.extra_channel.legend)
        state.extra_channel.legend = nothing
    end

    isempty(state.extra_channel.channels) && return

    margin = state.plot_kwargs[:channel_offset_margin]
    y_max = state.view.yrange[][end] * margin
    y_min = state.view.yrange[][1] * margin
    n_vis = count(state.channels.visible)
    default_step = n_vis > 1 ? (y_min - y_max) / (n_vis - 1) : -50.0

    step = (!state.view.butterfly[] && length(state.view.offset) > 1) ? (state.view.offset[end] - state.view.offset[end-1]) : 0.0
    abs(step) < 1e-4 && (step = default_step)
    base_offset = state.view.butterfly[] ? (y_min + step) : state.view.offset[end]

    get_data, get_time, get_label_y = _get_data_accessors(state.data)
    legend_elements = Any[]
    legend_labels = String[]
    epoch_obs = state.data isa EpochedDataState ? state.data.current_epoch : Observable(1)

    for (idx, channel) in enumerate(state.extra_channel.channels)
        channel_offset = base_offset + idx * step

        if eltype(get_data(state.data.current[], epoch_obs[], 1:1, channel)) == Bool
            highlight_data = @views _splitgroups(findall(get_data(state.data.current[], epoch_obs[], :, channel)))
            if !isempty(highlight_data[1])
                end_indices = copy(highlight_data[2])
                n_samp = length(get_time(state.data.current[], epoch_obs[], :))
                for i in eachindex(end_indices)
                    end_indices[i] == highlight_data[1][i] && (end_indices[i] = min(n_samp, end_indices[i] + 5))
                end
                state.extra_channel.data_lines[channel] = vspan!(
                    ax,
                    get_time(state.data.current[], epoch_obs[], highlight_data[1]),
                    get_time(state.data.current[], epoch_obs[], end_indices),
                    color = :black,
                    alpha = 0.5,
                    visible = true,
                )
            end
            push!(legend_elements, PolyElement(color = (:black, 0.5), strokecolor = :transparent))
            push!(legend_labels, String(channel))
        else
            time_obs = @lift(get_time($(state.data.current), $(epoch_obs), :))
            data_obs = @lift(get_data($(state.data.current), $(epoch_obs), :, $channel) .* $(state.view.amplitude_scale) .+ $channel_offset)
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
            append!(state.extra_channel.draw_observers, (time_obs, data_obs, lbl_x_obs, lbl_y_obs))

            state.extra_channel.data_lines[channel] = lines!(ax, time_obs, data_obs, color = :black, linewidth = 2)
            state.extra_channel.data_labels[channel] =
                text!(ax, lbl_x_obs, lbl_y_obs, text = String(channel), align = (:left, :center), color = :black, fontsize = 18)
        end
    end

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

_get_title(dat::EpochData) = "Epoch 1/$(n_epochs(dat))"
_get_title(::ContinuousData) = ""
_get_title(dat::ErpData) = "Epoch Average (n=$(n_epochs(dat)))"

function _plot_vertical_lines!(ax, marker, active)
    marker.line.visible = active
    marker.text.visible = active
    marker.visible = active
end

function _add_scale_indicator!(ax, state, plot_kwargs)
    scale_value = plot_kwargs[:scale_indicator_value]
    if !(plot_kwargs[:scale_indicator_position] isa Observable)
        plot_kwargs[:scale_indicator_position] = Observable{Tuple{Float64,Float64}}(Float64.(plot_kwargs[:scale_indicator_position]))
    end
    pos_obs = plot_kwargs[:scale_indicator_position]

    xlims_obs = ax.xaxis.attributes.limits
    ylims_obs = ax.yaxis.attributes.limits

    x_pos = @lift($xlims_obs[1] + ($xlims_obs[2] - $xlims_obs[1]) * $pos_obs[1])
    y_top = @lift($ylims_obs[1] + ($ylims_obs[2] - $ylims_obs[1]) * $pos_obs[2])
    y_bottom = @lift($y_top - scale_value * $(state.view.amplitude_scale))
    tick_length = @lift(($xlims_obs[2] - $xlims_obs[1]) * 0.005)

    segments = @lift([
        Point2f($x_pos, $y_bottom),
        Point2f($x_pos, $y_top),
        Point2f($x_pos - $tick_length, $y_bottom),
        Point2f($x_pos + $tick_length, $y_bottom),
        Point2f($x_pos - $tick_length, $y_top),
        Point2f($x_pos + $tick_length, $y_top),
    ])
    linesegments!(ax, segments, color = :black, linewidth = 1.5, overdraw = true)

    val_str = isinteger(scale_value) ? "$(Int(scale_value)) μV" : "$scale_value μV"
    text!(
        ax,
        @lift(Point2f($x_pos + ($xlims_obs[2] - $xlims_obs[1]) * 0.01, $y_top)),
        text = val_str,
        align = (:left, :center),
        fontsize = 14,
        color = :black,
        space = :data,
        overdraw = true,
    )

    is_dragging = Observable(false)

    on(events(ax).mousebutton, priority = 10) do event
        if event.button == Mouse.left
            if event.action == Mouse.press
                mpos = mouseposition(ax)
                x_tol = (xlims_obs[][2] - xlims_obs[][1]) * 0.05
                y_tol = (ylims_obs[][2] - ylims_obs[][1]) * 0.05
                if abs(mpos[1] - x_pos[]) < x_tol && mpos[2] >= (y_bottom[] - y_tol) && mpos[2] <= (y_top[] + y_tol)
                    is_dragging[] = true
                    return Consume(true)
                end
            elseif event.action == Mouse.release && is_dragging[]
                is_dragging[] = false
                return Consume(true)
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

# =============================================================================
# Entry Points & File Launcher
# =============================================================================

function _load_example_databrowser_data()
    raw_path = example_path("data/bdf/example1.bdf")
    lay = read_layout(example_path("layouts/biosemi/biosemi72.csv"))
    polar_to_cartesian_xy!(lay)
    dat = create_eegfun_data(read_raw_data(raw_path), lay)
    dat.file = basename(raw_path)
    return dat
end

"""
    plot_databrowser(dat::EegData, ica = nothing; screen = nothing, kwargs...)
    plot_databrowser(; launcher::Bool = false, screen = nothing, kwargs...)
    plot_databrowser(filename::String, ica = nothing; input_dir = pwd(), participant_selection = participants(), screen = nothing, kwargs...)
    plot_databrowser(data::Vector{<:EegData}, ica = nothing; screen = nothing, kwargs...)

Launch the interactive EEG data browser for continuous, epoched, or ERP datasets, with optional
dynamic ICA component subtraction, filtering, re-referencing, channel repair, and region selection.
Called with no arguments, opens the bundled `example1.bdf` dataset (`launcher = true` opens the
drag-and-drop file launcher).

$(_generate_kwargs_doc(PLOT_DATABROWSER_KWARGS))
"""
function plot_databrowser(dat::EegData, ica = nothing; screen = nothing, kwargs...)
    if string(Makie.current_backend()) == "CairoMakie"
        @minimal_warning "CairoMakie detected. For full interactivity in plot_databrowser, use GLMakie."
    end

    _set_window_title(_generate_window_title(dat))
    try
        plot_kwargs = _merge_plot_kwargs(PLOT_DATABROWSER_KWARGS, kwargs)
        fontsize_kw = isnothing(plot_kwargs[:theme_fontsize]) ? (;) : (; fontsize = plot_kwargs[:theme_fontsize])
        fig = _eegfun_figure(; size = plot_kwargs[:figure_size], figure_padding = plot_kwargs[:figure_padding], fontsize_kw...)
        ax = Axis(fig[1, 1], xlabel = plot_kwargs[:xlabel], ylabel = plot_kwargs[:ylabel], title = _get_title(dat))

        state = _create_browser_state(dat, dat.layout.data.label, ax, ica, plot_kwargs)
        state.extra_channel.fig = fig
        _setup_ui(fig, ax, state, dat, ica, plot_kwargs)
        _draw(ax, state)
        _draw_extra_channel!(ax, state)
        plot_kwargs[:show_scale_indicator] && _add_scale_indicator!(ax, state, plot_kwargs)

        screen_obj = nothing
        if plot_kwargs[:display_plot]
            if !isnothing(screen)
                display(screen, fig)
                screen_obj = screen
            else
                screen_obj = _display_figure(fig)
            end
        end

        return (fig = fig, axes = [ax], ax = ax, analysis_settings = state.analysis_settings, screen = screen_obj)
    finally
        _set_window_title("Makie")
    end
end

function plot_databrowser(; launcher::Bool = false, screen = nothing, kwargs...)
    launcher && return _plot_databrowser_launcher(; screen = screen, kwargs...)
    return plot_databrowser(_load_example_databrowser_data(); screen = screen, kwargs...)
end

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
        filepath != "" && isfile(filepath) && _load_layout(filepath)
    end

    on(btn_example.clicks) do _
        status_label.text = "Loading example1.bdf…"
        try
            dat = _load_example_databrowser_data()
            _set_window_title(_generate_window_title(dat))
            plot_databrowser(dat; screen = _create_screen(), kwargs...)
            n_loaded[] += 1
            status_label.text = "Opened example1.bdf"
        catch e
            @error "Failed to load example data" exception = (e, catch_backtrace())
            status_label.text = "Error: $(sprint(showerror, e))"
        end
    end

    on(files -> foreach(_load_file, files), events(fig).dropped_files)

    screen_obj = _display_in_screen(fig)
    _set_window_title("EegFun Data Browser")
    return (fig = fig, axes = Axis[], screen = screen_obj)
end

function _percent_decode(s::AbstractString)
    bytes = UInt8[]
    i = 1
    n = ncodeunits(s)
    while i <= n
        c = s[i]
        if c == '%' && i + 2 <= n && isxdigit(s[i+1]) && isxdigit(s[i+2])
            push!(bytes, parse(UInt8, s[(i+1):(i+2)], base = 16))
            i += 3
        else
            append!(bytes, codeunits(string(c)))
            i = nextind(s, i)
        end
    end
    return String(bytes)
end

function _sanitize_dropped_path(filepath::AbstractString)
    s = strip(filepath, ['\r', '\n', '"', '\''])
    if startswith(s, "file://")
        s = replace(s, r"^file://localhost" => "")
        s = replace(s, r"^file://" => "")
        s = _percent_decode(s)
        Sys.iswindows() && occursin(r"^/[A-Za-z]:", s) && (s = s[2:end])
    end
    return String(s)
end

function _read_any_eeg_file(filepath::String, layout::Union{Nothing,Layout} = nothing)
    ext = lowercase(splitext(filepath)[2])
    ica_data = nothing

    dat_eeg = if ext == ".jld2"
        read_data(filepath)
    elseif ext == ".set"
        res = read_eeglab(filepath)
        res isa Tuple ? ((dat_tmp, ica_data) = res; dat_tmp) : res
    else
        raw = read_raw_data(filepath)
        isnothing(layout) ? create_eegfun_data(raw) : create_eegfun_data(raw, layout)
    end

    isnothing(dat_eeg) && @minimal_error("Unsupported file format or failed to read: $filepath")
    return (dat_eeg, ica_data)
end

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
        !isnothing(ica) && (ica = read_data(ica))
        return plot_databrowser(dat_eeg, ica; screen = screen, kwargs...)
    end
    files = _require_batch_files(filename, input_dir, participant_selection)
    results = NamedTuple[]
    for file in files
        @info "Browsing: $file"
        dat_eeg = read_data(joinpath(input_dir, file))
        isnothing(dat_eeg) && continue
        push!(results, plot_databrowser(dat_eeg, nothing; screen = screen, kwargs...))
    end
    return results
end

function plot_databrowser(data::Vector{<:EegData}, ica = nothing; screen = nothing, kwargs...)
    @info "Vector of $(length(data)) datasets provided — opening a browser for each"
    results = NamedTuple[]
    display_plot = get(kwargs, :display_plot, true)
    for dat in data
        _set_window_title(_generate_window_title(dat))
        s = !isnothing(screen) ? screen : (display_plot ? _create_screen() : nothing)
        push!(results, plot_databrowser(dat, ica; screen = s, kwargs...))
    end
    return results
end
