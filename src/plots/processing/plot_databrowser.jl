const PLOT_DATABROWSER_KWARGS = Dict{Symbol,Tuple{Any,String}}(

    # Figure and layout
    :figure_padding => ((50, 0, 50, 50), "Figure padding as (left, right, bottom, top)"),

    # Axis styling
    :xlabel => ("Time (S)", "X-axis label"),
    :ylabel => ("Amplitude (μV)", "Y-axis label"),

    # UI styling
    :ui_fontsize => (18, "Font size for UI elements"),

    # Line styling
    :channel_line_width => (1, "Line width for channel lines"),
    :channel_line_alpha => (1, "Transparency for channel lines"),
    :selected_channel_color => (:black, "Color for selected channels"),
    :unselected_channel_color => (:darkgrey, "Color for unselected channels"),

    # TODO: could this be done better?
    :channel_offset_scale => (1500, "Scale factor for channel vertical offset"),
    :channel_offset_margin => (0.9, "Margin factor for channel offset range"),

    # Selection styling
    :selection_color => ((:blue, 0.1), "Color and transparency for selection rectangle"),

    # Filter parameters
    :default_hp_freq => (0.1, "Default high-pass filter frequency in Hz"),
    :default_lp_freq => (40.0, "Default low-pass filter frequency in Hz"),

    # Scale indicator
    :show_scale_indicator => (true, "Show scale indicator bar"),
    :scale_indicator_value => (100.0, "Scale indicator value in μV"),
    :scale_indicator_position => ((0.92, 0.96), "Scale indicator position as (x, y) in axis coordinates (0-1)"),
)

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


mutable struct SelectionState
    active::Observable{Bool}
    bounds::Observable{Tuple{Float64,Float64}}
    visible::Observable{Bool}
    rectangle::Makie.Poly
    selected_regions::Observable{Vector{Tuple{Float64,Float64}}}  # Store multiple regions
    region_plots::Vector{Makie.Poly}  # Store plot objects for each region
    function SelectionState(ax, plot_kwargs)
        poly_element = poly!(ax, [Point2f(0.0, 0.0)], color = plot_kwargs[:selection_color], visible = false)
        new(Observable(false), Observable((0.0, 0.0)), Observable(false), poly_element, Observable([]), [])
    end
end

mutable struct FilterState
    active::Observable{NamedTuple{(:hp, :lp),Tuple{Bool,Bool}}}
    hp_freq::Observable{Float64}
    lp_freq::Observable{Float64}
    FilterState(plot_kwargs) =
        new(Observable((hp = false, lp = false)), Observable(plot_kwargs[:default_hp_freq]), Observable(plot_kwargs[:default_lp_freq]))
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
    amplitude_scale::Observable{Float64}
    function ViewState(n_channels::Int, n_samples::Int, plot_kwargs)
        offset_scale = plot_kwargs[:channel_offset_scale]
        offset_margin = plot_kwargs[:channel_offset_margin]
        offset =
            n_channels > 1 ? LinRange(offset_scale * offset_margin, -offset_scale * offset_margin, n_channels + 2)[2:(end-1)] :
            zeros(n_channels)
        new(Observable(1:n_samples), Observable(-1500:1500), offset, Observable(0.0), Observable(false), Observable(1.0))
    end
end

mutable struct ChannelState
    labels::Vector{Symbol}
    visible::Vector{Bool}
    selected::Vector{Bool}
    individually_selected::Vector{Symbol}  # Track individually selected electrodes
    data_labels::Dict{Symbol,Makie.Text}
    data_lines::Dict{Symbol,Union{Makie.Lines,Makie.PolyElement,Any}}
    function ChannelState(channel_labels::Vector{Symbol})
        new(
            channel_labels,
            fill(true, length(channel_labels)),
            fill(false, length(channel_labels)),
            Symbol[],  # Start with empty list
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
    function ContinuousDataState(data::EegData, plot_kwargs)
        new(Observable(copy(data)), data, FilterState(plot_kwargs))
    end
end

mutable struct EpochedDataState <: AbstractDataState
    current::Observable{EegData}
    original::EegData
    filter_state::FilterState
    current_epoch::Observable{Int}
    function EpochedDataState(data::EegData, plot_kwargs)
        new(Observable(copy(data)), data, FilterState(plot_kwargs), Observable(1))
    end
end

mutable struct ExtraChannelInfo
    channel::Union{Nothing,Symbol}
    visible::Bool
    data_lines::Dict{Symbol,Union{Makie.Lines,Makie.PolyElement,Any}}
    data_labels::Dict{Symbol,Makie.Text}
    ExtraChannelInfo() = new(nothing, false, Dict(), Dict())
end


mutable struct DataBrowserState{T<:AbstractDataState}
    view::ViewState
    channels::ChannelState
    data::T
    selection::SelectionState
    markers::Vector{Marker}
    ica_original::Union{Nothing,InfoIca}
    ica_current::Union{Nothing,InfoIca}
    extra_channel::ExtraChannelInfo
    reference_state::Symbol
    channel_repair_history::Vector{Tuple{Vector{Symbol},Symbol,Matrix{Float64}}}  # (channels, method, original_data) - stack for multiple undos
    removed_ica_components::Vector{Int}  # Track removed ICA components
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
            ica_original,
            isnothing(ica_original) ? nothing : copy(ica_original),
            extra_channel,
            data.original.analysis_info.reference,
            [],  # empty repair history
            Int[],  # empty removed ICA components
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
    update_analysis_settings!(state::DataBrowserState)

Update the analysis settings observable with current state.
"""
function update_analysis_settings!(state::DataBrowserState)
    # Get current filter settings
    hp_freq = state.data.filter_state.active[].hp ? state.data.filter_state.hp_freq[] : 0.0
    lp_freq = state.data.filter_state.active[].lp ? state.data.filter_state.lp_freq[] : 0.0

    # Get repaired channels and their method
    repaired_channels = Symbol[]
    repair_method = :none
    for (channels, method, _) in state.channel_repair_history
        append!(repaired_channels, channels)
        repair_method = method  # TODO: should we check that same repair method was applied for all channel repairs?
    end

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
        state.removed_ica_components,
    )
end


"""
    _add_selected_regions!(data::EegData, selected_regions::Vector{Tuple{Float64,Float64}})

Add a boolean column :selected_region to the data indicating which time points are in selected regions.
"""
function _add_selected_regions!(dat::EegData, selected_regions::Vector{Tuple{Float64,Float64}})

    # Create boolean vector and mark where true
    selected_mask = falses(n_samples(dat))
    for (start_time, end_time) in selected_regions
        region_mask = (dat.data.time .>= start_time) .& (dat.data.time .<= end_time)
        selected_mask .|= region_mask
    end
    dat.data[!, :selected_region] = selected_mask

end

# EpochData version - apply to each epoch DataFrame
function _add_selected_regions!(dat::EpochData, selected_regions::Vector{Tuple{Float64,Float64}})
    for epoch_df in dat.data
        # Create boolean vector for this epoch
        selected_mask = falses(nrow(epoch_df))
        for (start_time, end_time) in selected_regions
            region_mask = (epoch_df.time .>= start_time) .& (epoch_df.time .<= end_time)
            selected_mask .|= region_mask
        end
        epoch_df[!, :selected_region] = selected_mask
    end
    return nothing
end




# Data browser state creation
function create_browser_state(dat::T, channel_labels, ax, ica, plot_kwargs) where {T<:EegData}
    state_type = data_state_type(T)
    initial_window = get_initial_window_size(dat)
    return DataBrowserState{state_type}(
        view = ViewState(length(channel_labels), initial_window, plot_kwargs),
        channels = ChannelState(channel_labels),
        data = state_type(dat, plot_kwargs),  # Pass kwargs to data state constructor
        selection = SelectionState(ax, plot_kwargs),
        ica_original = ica,
        plot_kwargs = plot_kwargs,
    )
end

# Helper function to get initial window size based on data type
get_initial_window_size(dat::ContinuousData) = min(5000, nrow(dat.data))  # Show reasonable window, not entire dataset
get_initial_window_size(dat::ErpData) = nrow(dat.data)  # Show whole epoch
get_initial_window_size(dat::EpochData) = nrow(dat.data[1])  # Show entire epoch

# Type mapping
data_state_type(::Type{ContinuousData}) = ContinuousDataState
data_state_type(::Type{EpochData}) = EpochedDataState
data_state_type(::Type{ErpData}) = ContinuousDataState

# Helper functions for common data access/resetting/updating
get_current_data(state::ContinuousDataState) = state.current[].data
get_current_data(state::EpochedDataState) = state.current[].data[state.current_epoch[]]
get_time_bounds(dat::ContinuousDataState) = (dat.current[].data.time[1], dat.current[].data.time[end])
get_time_bounds(dat::EpochedDataState) =
    (dat.current[].data[dat.current_epoch[]].time[1], dat.current[].data[dat.current_epoch[]].time[end])
has_column(state::ContinuousDataState, col::Symbol) = col in propertynames(state.current[].data)
has_column(state::EpochedDataState, col::Symbol) = col in propertynames(state.current[].data[state.current_epoch[]])
notify_data_update(state::AbstractDataState) = notify(state.current)
reset_to_original!(state::AbstractDataState) = state.current[] = copy(state.original)

######
# UI #
######
function setup_ui_base(fig, ax, state, dat, ica = nothing, plot_kwargs = nothing)

    deregister_interaction!(ax, :rectanglezoom) # need to turn this off!
    set_axes!(ax, state)

    # Mouse and keyboard events
    handle_mouse_events!(ax, state)
    handle_keyboard_events!(fig, ax, state)

    # Create toggles/markers/menus
    state.markers = init_markers(ax, state)
    toggles = create_toggles(fig, ax, state)
    labels_menu = create_labels_menu(fig, ax, state)
    reference_menu = create_reference_menu(fig, state, dat)

    # Create optional menus (ica/extra channels)
    ica_menu = nothing
    if !isnothing(ica) && !isnothing(state.ica_original)
        ica_menu = create_ica_menu(fig, ax, state, ica)
    end

    extra_menu = nothing
    extra_labels_result = extra_labels(state.data.original)
    if !isnothing(extra_labels_result) && !isempty(extra_labels_result)
        extra_menu = create_extra_channel_menu(fig, ax, state, dat)
    end

    return (toggles, labels_menu, reference_menu, ica_menu, extra_menu)
end

# Unified setup_ui method using multiple dispatch for the epoch menu
function setup_ui(fig, ax, state::DataBrowserState{<:AbstractDataState}, dat, ica = nothing, plot_kwargs = nothing)
    # Get common UI elements
    toggles, labels_menu, reference_menu, ica_menu, extra_menu = setup_ui_base(fig, ax, state, dat, ica, plot_kwargs)

    # Get type-specific epoch menu (or nothing)
    epoch_menu = get_epoch_menu(fig, ax, state)

    # Build the grid components
    build_grid_components!(fig, dat, state, toggles, labels_menu, reference_menu, ica_menu, extra_menu, epoch_menu)

    # Apply theme
    update_theme!(Theme(fontsize = plot_kwargs[:ui_fontsize]))
    hideydecorations!(ax, label = true)

    return state
end

function create_toggles(fig, ax, state)
    configs = [ToggleConfig("Butterfly Plot", (active) -> butterfly_plot!(ax, state))]

    # Add marker toggles based on configuration
    marker_toggle_configs = [(:triggers, "Trigger"), (:is_vEOG, "vEOG"), (:is_hEOG, "hEOG")]

    for (marker_symbol, toggle_label) in marker_toggle_configs
        if has_column(state.data, marker_symbol)
            marker_index = findfirst(m -> m.name == marker_symbol, state.markers)
            if !isnothing(marker_index)
                push!(configs, ToggleConfig(toggle_label, (active) -> plot_vertical_lines!(ax, state.markers[marker_index], active)))
            end
        end
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

    # Setup observers for toggle actions
    for toggle in toggles
        on(toggle[2].active) do active
            toggle[3](active)
        end
    end

    # Return as grid layout
    return permutedims(reduce(hcat, [[t[2], Label(fig, t[1], fontsize = 22, halign = :left)] for t in toggles]))

end

function create_menu(fig, options, default, label; kwargs...)
    menu = Menu(fig, options = options, default = default, direction = :down, fontsize = 18, width = Auto(), tellwidth = false, kwargs...)
    return hcat(menu, Label(fig, label, fontsize = 22, halign = :left, tellwidth = false))
end

function create_labels_menu(fig, ax, state)
    options = vcat(["All", "Left", "Right", "Central", "BioSemi16", "BioSemi32", "BioSemi64"], state.channels.labels)
    menu = create_menu(fig, options, "All", "Labels")

    # Flag to prevent recursive updates when programmatically setting menu selection
    updating_menu = Ref(false)

    on(menu[1].selection) do s
        # Skip if we're programmatically updating the menu
        if updating_menu[]
            return
        end
        # Group selections clear individual selections
        if s == "All"
            state.channels.individually_selected = Symbol[]
            state.channels.visible .= true
        elseif s == "Left"
            state.channels.individually_selected = Symbol[]
            state.channels.visible .= occursin.(r"\d*[13579]$", String.(state.channels.labels))
        elseif s == "Right"
            state.channels.individually_selected = Symbol[]
            state.channels.visible .= occursin.(r"\d*[24680]$", String.(state.channels.labels))
        elseif s == "Central"
            state.channels.individually_selected = Symbol[]
            state.channels.visible .= occursin.(r"z$", String.(state.channels.labels))
        elseif s == "BioSemi16"
            state.channels.individually_selected = Symbol[]
            package_layouts_dir = joinpath(@__DIR__, "..", "..", "data", "layouts")
            layout_file = find_file("biosemi16.csv", package_layouts_dir)
            if layout_file === nothing
                @minimal_error "Layout file biosemi16.csv not found in $package_layouts_dir"
                return
            end
            tmp_layout = read_layout(layout_file)
            state.channels.visible .= state.channels.labels .∈ Ref(tmp_layout.data.label)
        elseif s == "BioSemi32"
            state.channels.individually_selected = Symbol[]
            package_layouts_dir = joinpath(@__DIR__, "..", "..", "data", "layouts")
            layout_file = find_file("biosemi32.csv", package_layouts_dir)
            if layout_file === nothing
                @minimal_error "Layout file biosemi32.csv not found in $package_layouts_dir"
                return
            end
            tmp_layout = read_layout(layout_file)
            state.channels.visible .= state.channels.labels .∈ Ref(tmp_layout.data.label)
        elseif s == "BioSemi64"
            state.channels.individually_selected = Symbol[]
            package_layouts_dir = joinpath(@__DIR__, "..", "..", "data", "layouts")
            layout_file = find_file("biosemi64.csv", package_layouts_dir)
            if layout_file === nothing
                @minimal_error "Layout file biosemi64.csv not found in $package_layouts_dir"
                return
            end
            tmp_layout = read_layout(layout_file)
            state.channels.visible .= state.channels.labels .∈ Ref(tmp_layout.data.label)
        else
            # Individual electrode selection - toggle selection
            selected_sym = Symbol(s)
            if selected_sym in state.channels.individually_selected
                # Deselect if already selected
                filter!(x -> x != selected_sym, state.channels.individually_selected)
            else
                # Select if not already selected
                push!(state.channels.individually_selected, selected_sym)
            end
            # Update visibility to show only individually selected electrodes
            if isempty(state.channels.individually_selected)
                # If no individual selections, show all channels
                state.channels.visible .= true
            else
                state.channels.visible .= (state.channels.labels .∈ Ref(state.channels.individually_selected))
            end
            # Always reset menu to "All" after toggling any individual channel
            # This ensures clicking any channel (including the same one) will always trigger the callback
            updating_menu[] = true
            menu[1].selection[] = "All"
            updating_menu[] = false
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
        update_analysis_settings!(state)
    end

    return menu
end

# Update create_ica_menu to use the new struct
function create_ica_menu(fig, ax, state, ica)

    options = vcat(["None"], ica.ica_label)
    menu = create_menu(fig, options, "None", "ICA Components")

    # Create observable for removed components display
    removed_components_text = Observable("")

    # Create label to display removed components
    removed_label = Label(
        fig,
        @lift(isempty($removed_components_text) ? "" : "Removed: $($removed_components_text)"),
        fontsize = 16,
        halign = :left,
        color = :red,
        tellwidth = false,
    )

    # Function to update the removed components display
    function update_removed_display()
        if isempty(state.removed_ica_components)
            removed_components_text[] = ""
        else
            # Sort components for consistent display
            sorted_components = sort(state.removed_ica_components)
            components_str = join(string.(sorted_components), ", ")
            removed_components_text[] = components_str
        end
    end

    # Initialize display
    update_removed_display()

    on(menu[1].selection) do s

        clear_axes!(ax, [state.channels.data_lines, state.channels.data_labels])

        # Check explicitly for "None" first
        if s == "None"
            # Selected "None" - clear all removals and restore original data
            # Clear removed ICA components first
            empty!(state.removed_ica_components)

            # Reset to original ICA state
            state.ica_current = copy(state.ica_original)

            # Reset data to original (no components removed)
            reset_to_original!(state.data)
        else
            # Extract component number from selection string
            component_to_toggle_int = extract_int(String(s))

            if !isnothing(component_to_toggle_int)
                # Toggle component: if already removed, restore it; if not removed, remove it
                if component_to_toggle_int in state.removed_ica_components
                    # Restore component - remove from list
                    filter!(x -> x != component_to_toggle_int, state.removed_ica_components)
                else
                    # Remove component - add to list
                    push!(state.removed_ica_components, component_to_toggle_int)
                end

                # Always reset to original data first
                state.ica_current = copy(state.ica_original)
                reset_to_original!(state.data)

                # Then apply all currently removed components
                if !isempty(state.removed_ica_components)
                    apply_ica_removal!(state.data, state.ica_current, state.removed_ica_components)
                end

            end
        end

        # Update the removed components display
        update_removed_display()

        notify_data_update(state.data)
        draw(ax, state)
        update_analysis_settings!(state)
    end

    # Return menu and removed components label in a vertical stack
    return vcat(menu, hcat(removed_label, Label(fig, "", fontsize = 16, tellwidth = false)))

end

function create_epoch_menu(fig, ax, state)
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
            clear_axes!(ax, [state.channels.data_lines, state.channels.data_labels])
            state.data.current_epoch[] = epoch_num
            ax.title = "Epoch $(epoch_num)/$(n_epochs(state.data.original))"
            update_markers!(ax, state)
            draw(ax, state)
            draw_extra_channel!(ax, state)
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

function show_additional_menu(state, clicked_region_idx = nothing)

    # Create the menu figure
    # TODO: why does new window not always take this size?
    menu_fig = Figure(size = (300, 300))
    plot_types = ["Topoplot", "Spectrum", "Get Selected Regions"]

    menu_buttons = [Button(menu_fig[idx, 1], label = plot_type) for (idx, plot_type) in enumerate(plot_types)]

    for btn in menu_buttons
        on(btn.clicks) do _
            if btn.label[] == "Get Selected Regions"
                # Get the boolean vector of selected regions
                selected_regions_bool = get_selected_regions_bool(state)
                @info "Selected regions boolean vector:"
                @info "Length: $(length(selected_regions_bool))"
                @info "Number of selected samples: $(sum(selected_regions_bool))"
                @info "Selected regions: $(state.selection.selected_regions[])"
            else
                selected_data = subset_selected_data(state, clicked_region_idx)
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

function _channel_repair_menu(state, selected_channels, ax)
    # Get all available channels
    all_channels = state.channels.labels
    n_channels = length(all_channels)

    # TODO: this looks ok for my typical 70/72 channel setup but could be improved for other setups

    cols = 7  # 7 columns

    # Create the repair menu figure
    menu_fig = Figure(size = (700, 800))

    # Add title
    title_text = "Channel Repair Interface"
    Label(menu_fig[1, 1], title_text, fontsize = 18, halign = :center)

    # Create a scrollable area for channels
    scroll_area = menu_fig[2, 1] = GridLayout()

    # Create checkboxes and labels arrays
    channel_checkboxes = []
    channel_labels = []

    # Add all channels in 7-column layout
    for (i, ch) in enumerate(all_channels)
        row = ((i - 1) ÷ cols) + 1
        col = ((i - 1) % cols) + 1

        # Check if channel has been repaired
        is_repaired = false
        for (repaired_channels, _, _) in state.channel_repair_history
            if ch in repaired_channels
                is_repaired = true
                break
            end
        end

        # Create a horizontal layout for checkbox and label
        channel_cell = scroll_area[row, col] = GridLayout()

        # Create checkbox
        checkbox = Checkbox(channel_cell[1, 1], checked = ch in selected_channels, width = 20, height = 20)
        push!(channel_checkboxes, checkbox)

        # Create label with repair status
        label_text = is_repaired ? "$(string(ch)) ✓" : string(ch)
        label_color = is_repaired ? :green : :black
        label = Label(channel_cell[1, 2], label_text, fontsize = 12, color = label_color, halign = :left)
        push!(channel_labels, label)
    end

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

    # Method selection (radio button behavior)
    selected_method = Observable(:neighbor_interpolation)

    on(method_buttons[1].clicks) do n
        selected_method[] = :neighbor_interpolation
        method_buttons[1].buttoncolor[] = :lightblue
        method_buttons[2].buttoncolor[] = :white
    end

    on(method_buttons[2].clicks) do n
        selected_method[] = :spherical_spline
        method_buttons[1].buttoncolor[] = :white
        method_buttons[2].buttoncolor[] = :lightblue
    end

    # Initialize with neighbor interpolation selected as default
    selected_method[] = :neighbor_interpolation

    # Apply repair
    on(action_buttons[1].clicks) do n
        selected_channels = all_channels[findall(cb -> cb.checked[], channel_checkboxes)]
        if !isempty(selected_channels)
            repair_selected_channels!(state, selected_channels, selected_method[], ax)
        else
            @info "No channels selected for repair"
        end
    end

    # Undo last repair
    on(action_buttons[2].clicks) do n
        if !isempty(state.channel_repair_history)
            undo_last_repair!(state, ax)
        else
            @info "No repairs to undo"
        end
    end

    new_screen = GLMakie.Screen()
    display(new_screen, menu_fig)
end

function repair_selected_channels!(state, selected_channels, method, ax)
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
    original_data = copy(get_channel_data_matrix(state.data.current[], selected_channels))

    # Perform the repair
    if method == :neighbor_interpolation
        repair_channels!(state.data.current[], selected_channels, method = :neighbor_interpolation)
    elseif method == :spherical_spline
        repair_channels!(state.data.current[], selected_channels, method = :spherical_spline)
    end

    # Store repair in history
    push!(state.channel_repair_history, (selected_channels, method, original_data))

    # Notify that data has been updated
    notify_data_update(state.data)

    # Update analysis settings
    update_analysis_settings!(state)

    # Clear and redraw the plot
    clear_axes!(ax, [state.channels.data_lines, state.channels.data_labels])
    draw(ax, state)

    total_repairs = length(state.channel_repair_history)
    @info "Successfully repaired channels: $(join(string.(selected_channels), ", ")) using $method"
    @info "Total repairs in history: $total_repairs"
end

function undo_last_repair!(state, ax)

    if isempty(state.channel_repair_history)
        @info "No repairs to undo"
        return
    end

    # Get the last repair
    channels, method, original_data = pop!(state.channel_repair_history)

    # Restore original data
    restore_channel_data!(state.data.current[], channels, original_data)

    # Notify that data has been updated
    notify_data_update(state.data)

    # Clear and redraw the plot
    clear_axes!(ax, [state.channels.data_lines, state.channels.data_labels])
    draw(ax, state)

    remaining_repairs = length(state.channel_repair_history)
    @info "Undid repair of channels: $(join(string.(channels), ", ")) (was $method)"
    @info "Remaining repairs in history: $remaining_repairs"
end

# Helper function to get channel data matrix
function get_channel_data_matrix(data, channels)
    if hasfield(typeof(data), :data) && hasfield(typeof(data), :layout)
        channel_data = data.data[:, channels]
        return Matrix(channel_data)
    else
        throw(ArgumentError("Unsupported data type for channel repair"))
    end
end

# Helper function to restore channel data
function restore_channel_data!(data, channels, original_data)
    if hasfield(typeof(data), :data) && hasfield(typeof(data), :layout)
        data.data[:, channels] = original_data
    else
        throw(ArgumentError("Unsupported data type for channel repair"))
    end
end


# Create common sliders for both continuous and epoched data
function create_common_sliders(fig, state, dat)
    sliders = []

    # Extreme value slider
    slider_extreme = Slider(fig[1, 2], range = 0:5:100, startvalue = 0, width = 100)
    on(slider_extreme.value) do x
        state.view.crit_val[] = x
    end
    push!(sliders, hcat(slider_extreme, Label(fig, @lift("Extreme: $($(slider_extreme.value)) μV"), fontsize = 22)))

    # Define filter slider configurations
    filter_configs = [(:hp_filter, :hp_freq, 0.1:0.1:2, 0.5, "HP-Filter"), (:lp_filter, :lp_freq, 5:5:60, 20, "LP-Filter")]

    # Create filter sliders based on configuration
    for (filter_field, freq_field, range, startval, label) in filter_configs
        if getfield(dat.analysis_info, filter_field) == 0.0
            slider = Slider(fig[1, 2], range = range, startvalue = startval, width = 100)
            on(slider.value) do val
                getfield(state.data.filter_state, freq_field)[] = val
            end
            # Initialize filter state with slider default value
            getfield(state.data.filter_state, freq_field)[] = startval
            push!(sliders, hcat(slider, Label(fig, @lift("$label: $($(slider.value)) Hz"), fontsize = 22)))
        end
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
    menu = Menu(fig, options = [:none; extra_labels(dat)], default = "none", direction = :down, fontsize = 18, width = 200)

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
const KEYBOARD_ACTIONS = Dict(Keyboard.left => :left, Keyboard.right => :right, Keyboard.up => :up, Keyboard.down => :down)

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
_handle_left_navigation(ax, state, data::ContinuousDataState) = xback!(ax, state)
_handle_left_navigation(ax, state, data::EpochedDataState) = step_epoch_backward(ax, state)
_handle_right_navigation(ax, state, data::ContinuousDataState) = xforward!(ax, state)
_handle_right_navigation(ax, state, data::EpochedDataState) = step_epoch_forward(ax, state)

function xback!(ax, state::ContinuousDataBrowserState)
    state.view.xrange.val[1] - 200 < 1 && return
    state.view.xrange[] = state.view.xrange.val .- 200
    xlims!(ax, state.data.current[].data.time[state.view.xrange.val[1]], state.data.current[].data.time[state.view.xrange.val[end]])
end

function xforward!(ax, state::ContinuousDataBrowserState)
    state.view.xrange.val[1] + 200 > nrow(state.data.current[].data) && return
    state.view.xrange[] = state.view.xrange.val .+ 200
    xlims!(ax, state.data.current[].data.time[state.view.xrange.val[1]], state.data.current[].data.time[state.view.xrange.val[end]])
end

step_epoch_backward(ax, state::EpochedDataBrowserState) = step_epoch!(ax, state, -1)
step_epoch_forward(ax, state::EpochedDataBrowserState) = step_epoch!(ax, state, 1)

function step_epoch!(ax, state::EpochedDataBrowserState, direction::Int)
    clear_axes!(ax, [state.channels.data_lines, state.channels.data_labels])
    current = state.data.current_epoch[]
    total = n_epochs(state.data.original)
    state.data.current_epoch[] = clamp(current + direction, 1, total)
    ax.title = "Epoch $(state.data.current_epoch[])/$total"
    update_markers!(ax, state)
    draw(ax, state)
    draw_extra_channel!(ax, state)
end

yless!(ax, state) = yzoom!(ax, state, 1.2)
ymore!(ax, state) = yzoom!(ax, state, 0.8)

function yzoom!(ax, state, factor::Float64)
    if state.view.butterfly[]
        # In butterfly mode: adjust y-range
        y_min, y_max = state.view.yrange.val[1], state.view.yrange.val[end]
        if factor > 1.0  # Zoom in (yless)
            (y_min + 100 >= 0 || y_max - 100 <= 0) && return
            state.view.yrange[] = (y_min+100):(y_max-100)
        else  # Zoom out (ymore)
            state.view.yrange[] = (y_min-100):(y_max+100)
        end
        ylims!(ax, state.view.yrange.val[1], state.view.yrange.val[end])
    else # In non-butterfly mode: adjust amplitude scale
        state.view.amplitude_scale[] = state.view.amplitude_scale[] * factor
    end
end

function is_mouse_in_axis(ax, pos)
    bbox = ax.layoutobservables.computedbbox[]
    return bbox.origin[1] <= pos[1] <= (bbox.origin[1] + bbox.widths[1]) && bbox.origin[2] <= pos[2] <= (bbox.origin[2] + bbox.widths[2])
end

function find_clicked_region(state, mouse_x)
    # Check if click is within any of the selected regions
    regions = state.selection.selected_regions[]
    for (i, (start_time, end_time)) in enumerate(regions)
        if mouse_x >= min(start_time, end_time) && mouse_x <= max(start_time, end_time)
            return i
        end
    end
    return nothing
end

function remove_region_from_selection!(ax, state, region_idx)
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
    state.selection.rectangle.visible[] = true

    # Add this selection to the list of selected regions
    add_region_to_selection!(ax, state, state.selection.bounds[][1], mouse_x)

    # Clear the temporary selection rectangle after adding to permanent regions
    clear_x_region_selection!(state)
end

function handle_mouse_events!(ax, state)
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
        if !is_mouse_in_axis(ax, pos)
            return
        end

        mouse_x = mouseposition(ax)[1]

        if event.button == Mouse.left
            if event.action == Mouse.press
                if ctrl_pressed[]
                    mouse_y = mouseposition(ax)[2]
                    clicked_channel = find_closest_channel(ax, state, mouse_x, mouse_y)
                    if !isnothing(clicked_channel)
                        toggle_channel_visibility!(ax, state, clicked_channel)
                        return
                    end
                elseif shift_pressed[]
                    handle_left_click!(ax, state, event, mouse_x)
                end
            elseif event.action == Mouse.release
                if shift_pressed[]
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
        # Check if click is within any existing selected region
        clicked_region_idx = find_clicked_region(state, mouse_x)
        if clicked_region_idx !== nothing
            # Remove the clicked region
            remove_region_from_selection!(ax, state, clicked_region_idx)
        else
            # Start a new selection
            start_selection!(ax, state, mouse_x)
        end
    elseif event.action == Mouse.release && state.selection.active[]
        finish_selection!(ax, state, mouse_x)
    end
end

function handle_right_click!(ax, state, mouse_x)
    # Check if right-click is within any selected region
    clicked_region_idx = find_clicked_region(state, mouse_x)
    if clicked_region_idx !== nothing
        show_additional_menu(state, clicked_region_idx)
    end
    # Right-click outside regions does nothing (use 'r' key for channel repair)
end

function _show_channel_repair_menu(state::DataBrowserState{<:ContinuousDataState}, ax)
    _channel_repair_menu(state, get_selected_channels(state), ax)
end

function _show_channel_repair_menu(state::DataBrowserState{<:EpochedDataState}, ax)
    @info "Channel repair is only available for continuous data"
end

# Helper function to get selected channels
function get_selected_channels(state)
    selected_indices = findall(state.channels.selected)
    return state.channels.labels[selected_indices]
end

# Helper function to find the closest channel to a click
function find_closest_channel(ax, state, mouse_x, mouse_y)
    current_data = get_current_data(state.data)
    tolerance = 10  # 50 pixel tolerance

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

            # If this is the closest channel so far
            if distance < min_distance
                min_distance = distance
                closest_channel = idx
            end
            if distance < tolerance
                return closest_channel
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
        if event.action == Keyboard.press && event.key == Keyboard.i
            # Show help for databrowser
            show_plot_help(:databrowser)
        elseif event.action == Keyboard.press && event.key == Keyboard.r
            # Open channel repair menu
            _show_channel_repair_menu(state, ax)
        elseif event.action == Keyboard.press && event.key == Keyboard.c
            # Clear all selected regions
            clear_all_selected_regions!(ax, state)
        elseif event.action in (Keyboard.press, Keyboard.repeat) && haskey(KEYBOARD_ACTIONS, event.key)
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

function _handle_selection_movement_impl(ax, state::DataBrowserState{<:AbstractDataState}, action::Symbol)
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

function update_x_region_selection!(ax, state, x1, x2)
    ylims = ax.limits[][2]
    state.selection.rectangle[1] = Point2f[
        Point2f(Float64(x1), Float64(ylims[1])),
        Point2f(Float64(x2), Float64(ylims[1])),
        Point2f(Float64(x2), Float64(ylims[2])),
        Point2f(Float64(x1), Float64(ylims[2])),
    ]
    state.selection.rectangle.visible[] = true
end

function add_region_to_selection!(ax, state, x1, x2)
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
    region_plot = poly!(ax, region_points, color = state.plot_kwargs[:selection_color], strokecolor = :transparent)
    push!(state.selection.region_plots, region_plot)

    # Update analysis settings
    update_analysis_settings!(state)

end

function clear_x_region_selection!(state)
    # Set to a single point instead of empty vector to avoid CairoMakie issues
    state.selection.rectangle[1] = [Point2f(0.0, 0.0)]
    state.selection.bounds[] = (0.0, 0.0)
    state.selection.visible[] = false
    state.selection.rectangle.visible[] = false
end

function clear_all_selected_regions!(ax, state)
    # Clear all region plots
    for plot in state.selection.region_plots
        delete!(ax.scene, plot)
    end
    empty!(state.selection.region_plots)

    # Clear the selected regions list
    state.selection.selected_regions[] = Tuple{Float64,Float64}[]

end

"""
    get_selected_regions_bool(state::DataBrowserState) -> Vector{Bool}

Returns a boolean vector indicating which samples are within the selected regions.
The vector has the same length as the total number of samples in the data.
"""
function get_selected_regions_bool(state::DataBrowserState)
    current_data = get_current_data(state.data)
    total_samples = nrow(current_data)
    time_data = current_data.time
    bool_vector = falses(total_samples)

    for (start_time, end_time) in state.selection.selected_regions[]
        # Find the closest sample indices
        start_idx = argmin(abs.(time_data .- start_time))
        end_idx = argmin(abs.(time_data .- end_time))

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
    get_selected_regions_info(state::DataBrowserState) -> NamedTuple

Returns detailed information about the selected regions including:
- `bool_vector`: Boolean vector of selected samples
- `regions`: List of (start_time, end_time) tuples
- `n_samples`: Number of selected samples
- `n_regions`: Number of selected regions
"""
function get_selected_regions_info(state::DataBrowserState)
    bool_vector = get_selected_regions_bool(state)
    regions = state.selection.selected_regions[]

    return (bool_vector = bool_vector, regions = regions, n_samples = sum(bool_vector), n_regions = length(regions))
end

function subset_selected_data(state::ContinuousDataBrowserState, clicked_region_idx = nothing)
    # Use the clicked region if specified, otherwise use the most recent region, or fall back to bounds
    if clicked_region_idx !== nothing && 1 <= clicked_region_idx <= length(state.selection.selected_regions[])
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

function subset_selected_data(state::EpochedDataBrowserState, clicked_region_idx = nothing)
    # Get the selected region
    if clicked_region_idx !== nothing && 1 <= clicked_region_idx <= length(state.selection.selected_regions[])
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

    return convert(epoch_data, 1)
end


############
# Filtering
############
function apply_filter!(state::DataBrowserState{T}, filter_type, freq) where {T<:AbstractDataState}
    # Get the current data, apply filter, then update the observable
    current_data = state.data.current[]
    if filter_type == :hp
        highpass_filter!(current_data, freq)
    elseif filter_type == :lp
        lowpass_filter!(current_data, freq)
    end
    state.data.current[] = current_data  # Explicitly update the observable
end

function apply_filters!(state)
    # Reset to original if no filters active
    if !state.data.filter_state.active[].hp && !state.data.filter_state.active[].lp
        reset_to_original!(state.data)
        rereference!(state.data, state.reference_state)
        notify_data_update(state.data)
        return
    end

    # Always start with fresh data when applying filters to ensure clean filtering
    reset_to_original!(state.data)
    rereference!(state.data, state.reference_state)

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
    new_state = (hp = !current_state.hp, lp = current_state.lp)
    state.data.filter_state.active[] = new_state
    apply_filters!(state)
    update_analysis_settings!(state)
end

function apply_lp_filter!(state)
    current_state = state.data.filter_state.active[]
    new_state = (hp = current_state.hp, lp = !current_state.lp)
    state.data.filter_state.active[] = new_state
    apply_filters!(state)
    update_analysis_settings!(state)
end



########################
# Reference
########################
function rereference!(state::AbstractDataState, ref)
    rereference!(state.current[], ref, channels())
    notify_data_update(state)  # Notify that data has been updated
end

########################
# ICA Component Application Helpers
########################

# Apply ICA component removal based on state type
function apply_ica_removal!(state::ContinuousDataState, ica::InfoIca, components_to_remove::Vector{Int})
    remove_ica_components!(state.current[].data, ica, component_selection = components(components_to_remove))
    return nothing  # Activations are now stored in ica.removed_activations
end

function apply_ica_removal!(state::EpochedDataState, ica::InfoIca, components_to_remove::Vector{Int})
    for (i, epoch_df) in enumerate(state.current[].data)
        remove_ica_components!(epoch_df, ica, component_selection = components(components_to_remove))
    end
    return nothing  # Activations are now stored in ica.removed_activations
end


########################
# Drawing
########################
function add_marker!(markers, ax, data, col; label = nothing, trial = nothing, visible = false)
    if isnothing(trial)
        # More efficient: filter directly without findall
        mask = data[!, col] .!= 0
        marker_data = data[mask, [:time, col]]
    else
        mask = data[trial][!, col] .!= 0
        marker_data = data[trial][mask, [:time, col]]
    end

    # if no markers, return
    nrow(marker_data) == 0 && return

    label = isnothing(label) ? string.(marker_data[!, col]) : repeat([label], nrow(marker_data))

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
        y_max = state.view.yrange[][end] * 0.9
        y_min = state.view.yrange[][1] * 0.9
        step = (y_min - y_max) / (nchannels - 1)
        visible_indices = findall(state.channels.visible)
        for (i, idx) in enumerate(visible_indices)
            state.view.offset[idx] = y_max + (i - 1) * step
        end
    else
        state.view.offset[state.channels.visible] .= 0.0
    end
end

function clear_axes!(ax, datas)
    for data in datas
        for (key, value) in data
            delete!(ax, value)
        end
        empty!(data)
    end
end

function set_axes!(ax, state::DataBrowserState{<:AbstractDataState})
    @lift ylims!(ax, $(state.view.yrange)[1], $(state.view.yrange)[end])
    set_x_limits!(ax, state, state.data)
end

function set_x_limits!(ax, state, data::ContinuousDataState)
    @lift xlims!(ax, $(data.current).data.time[$(state.view.xrange)[1]], $(data.current).data.time[$(state.view.xrange)[end]])
end

function set_x_limits!(ax, state, data::EpochedDataState)
    @lift xlims!(ax, $(data.current).data[1].time[1], $(data.current).data[1].time[end])
end

# Common marker initialization
function init_markers(ax, state; marker_visible = Dict{Symbol,Bool}())
    markers = Marker[]
    data = get_current_data(state.data)

    # Define marker configurations
    marker_configs = [(:triggers, nothing), (:is_vEOG, "v"), (:is_hEOG, "h")]

    # Add markers based on configuration
    for (symbol, label) in marker_configs
        if has_column(state.data, symbol)
            add_marker!(markers, ax, data, symbol, label = label, visible = get(marker_visible, symbol, false))
        end
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

            # Channel data (compute once)
            channel_data_obs = @lift(get_data($(state.data.current), $(state.view.xrange), $col))
            channel_data_with_offset = @lift($(channel_data_obs) .* $(state.view.amplitude_scale) .+ state.view.offset[idx])

            # Check if channel is repaired
            is_repaired = false
            for (repaired_channels, _, _) in state.channel_repair_history
                if col in repaired_channels
                    is_repaired = true
                    break
                end
            end

            # Line properties (reuse channel_data_obs for efficiency)
            if is_repaired || is_selected
                # Repaired channels get black color and thicker lines
                line_color = state.plot_kwargs[:selected_channel_color]
                line_colormap = [:black]
                line_width = state.plot_kwargs[:channel_line_width] * 2  # Make repaired channels thicker
            else
                # Normal channels
                line_color = @lift(abs.($(channel_data_obs)) .>= $(state.view.crit_val))
                line_colormap = [state.plot_kwargs[:unselected_channel_color], state.plot_kwargs[:unselected_channel_color], :red]
                line_width = state.plot_kwargs[:channel_line_width]
            end

            # Update or create line
            create_line!(
                state.channels.data_lines,
                col,
                ax,
                visible_time_obs,
                channel_data_with_offset,
                line_color,
                line_colormap,
                line_width,
                state.plot_kwargs[:channel_line_alpha],
            )

            # Handle labels
            if !state.view.butterfly[]
                label_y_obs = @lift(get_label_y($(state.data.current), $col, state.view.offset[idx]))
                create_label!(state.channels.data_labels, col, ax, time_start_obs, label_y_obs, is_selected)
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

function get_data_accessors(state::EpochedDataState)
    get_data = (current, range, col) -> current.data[state.current_epoch[]][range, col]
    get_time = (current, range) -> current.data[state.current_epoch[]].time[range]
    get_label_y = (current, col, offset) -> current.data[state.current_epoch[]][!, col][1] .+ offset
    return get_data, get_time, get_label_y
end

# Helper functions for line/label management
function create_line!(data_lines, col, ax, x_obs, y_obs, color, colormap, linewidth, alpha)
    data_lines[col] = lines!(ax, x_obs, y_obs, color = color, colormap = colormap, linewidth = linewidth, alpha = alpha)
end

function create_label!(data_labels, col, ax, x_obs, y_obs, is_selected)
    data_labels[col] =
        text!(ax, x_obs, y_obs, text = String(col), align = (:left, :center), fontsize = 18, color = is_selected ? :red : :black)
end

function hide_channel_label!(data_labels, col)
    haskey(data_labels, col) && hide!(data_labels[col])
end

function hide_channel_objects!(channels, col)
    haskey(channels.data_lines, col) && hide!(channels.data_lines[col])
    haskey(channels.data_labels, col) && hide!(channels.data_labels[col])
end

# Single function with data access abstraction
function draw_extra_channel!(ax, state::DataBrowserState{<:AbstractDataState})
    clear_axes!(ax, [state.extra_channel.data_lines, state.extra_channel.data_labels])

    if state.extra_channel.visible && !isnothing(state.extra_channel.channel)
        if length(state.view.offset) > 1
            current_offset = state.view.offset[end] + (state.view.offset[end] - state.view.offset[end-1])
        else
            current_offset = state.view.offset[end] + 100.0  # Default spacing
        end
        channel = state.extra_channel.channel

        # Get data access functions based on type
        get_data, get_time, get_label_y = get_data_accessors(state.data)

        if eltype(get_data(state.data.current[], 1:1, channel)) == Bool # Boolean data - create highlights
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
        else # Regular data - create line and label
            state.extra_channel.data_lines[channel] = lines!(
                ax,
                @lift(get_time($(state.data.current), :)),
                @lift(get_data($(state.data.current), :, $channel) .* $(state.view.amplitude_scale) .+ $current_offset),
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
get_title(dat::ErpData) = "Epoch Average (n=$(n_epochs(dat)))"

function plot_vertical_lines!(ax, marker, active)
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


function plot_databrowser(dat::EegData, ica = nothing; screen = nothing, kwargs...)

    # Check if CairoMakie is being used and warn about lack of interactivity
    if string(Makie.current_backend()) == "CairoMakie"
        @minimal_warning "CairoMakie detected. For full interactivity in plot_databrowser, use GLMakie."
    end

    # Generate window title from dataset
    title_str = _generate_window_title(dat)
    set_window_title(title_str)

    # Merge user kwargs with defaults
    plot_kwargs = _merge_plot_kwargs(PLOT_DATABROWSER_KWARGS, kwargs)

    # Common fig/ax/state/ui setup
    fig = Figure(figure_padding = plot_kwargs[:figure_padding])
    ax = Axis(fig[1, 1], xlabel = plot_kwargs[:xlabel], ylabel = plot_kwargs[:ylabel], title = get_title(dat))

    state = create_browser_state(dat, dat.layout.data.label, ax, ica, plot_kwargs)
    setup_ui(fig, ax, state, dat, ica, plot_kwargs)

    draw(ax, state)
    draw_extra_channel!(ax, state)

    # Add scale indicator
    if plot_kwargs[:show_scale_indicator]
        _add_scale_indicator!(ax, state, plot_kwargs)
    end

    # Display on the provided screen if given, otherwise use default display
    if screen !== nothing
        display(screen, fig)
    else
        display(fig)
    end

    set_window_title("Makie")
    # Return the observable analysis settings
    return (fig = fig, ax = ax, analysis_settings = state.analysis_settings)
end

# can plot saved file
function plot_databrowser(filename::String, ica = nothing; screen = nothing, kwargs...)
    dat_eeg = read_data(filename)
    if !isnothing(ica)
        ica = read_data(ica)
    end
    plot_databrowser(dat_eeg, ica; screen = screen, kwargs...)
end

plot_databrowser(data::Vector{<:EegData}, ica = nothing; screen = nothing, kwargs...) =
    plot_databrowser(data[1], ica; screen = screen, kwargs...);
