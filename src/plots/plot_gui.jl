"""Shared font sizes and widget dimensions for the interactive GUI."""
const UI_STYLE = (
    label_font   = 20,
    button_font  = 18,
    textbox_font = 16,
    input_width  = 160,
    input_height = 30,
)

Base.@kwdef mutable struct GUIState
    directory::Observable{String} = Observable("")
    filename::Observable{String} = Observable("")
    participant::Observable{String} = Observable("")
    condition::Observable{String} = Observable("")
    epoch::Observable{String} = Observable("")
    plottype::Observable{String} = Observable("select")
    layout::Observable{String} = Observable("select")
    layout_file::Observable{String} = Observable("")
    layout_object::Observable{Union{Nothing,Layout}} = Observable{Union{Nothing,Layout}}(nothing)
    electrodes::Observable{Vector{String}} = Observable(String[])
    xlim::Observable{Tuple{Union{Nothing,Float64},Union{Nothing,Float64}}} =
        Observable{Tuple{Union{Nothing,Float64},Union{Nothing,Float64}}}((nothing, nothing))
    ylim::Observable{Tuple{Union{Nothing,Float64},Union{Nothing,Float64}}} =
        Observable{Tuple{Union{Nothing,Float64},Union{Nothing,Float64}}}((nothing, nothing))
    zlim::Observable{Tuple{Union{Nothing,Float64},Union{Nothing,Float64}}} =
        Observable{Tuple{Union{Nothing,Float64},Union{Nothing,Float64}}}((nothing, nothing))
    baseline_start::Observable{Union{Nothing,Float64}} = Observable{Union{Nothing,Float64}}(nothing)
    baseline_end::Observable{Union{Nothing,Float64}} = Observable{Union{Nothing,Float64}}(nothing)
    baseline_type::Observable{String} = Observable("select")
    layout_type::Observable{String} = Observable("single")
    average_channels::Observable{Bool} = Observable(false)
    invert_y::Observable{Bool} = Observable(false)
    submenu_active::Observable{Bool} = Observable(false)
    submenu_type::Observable{String} = Observable("")
    channel_menu::Union{Nothing,Menu} = nothing
    ica_components::Observable{String} = Observable("")
end

"""Create a styled button for GUI selection actions."""
function _create_select_button(parent, label)
    Button(
        parent,
        label = label,
        fontsize = UI_STYLE.button_font,
        width = UI_STYLE.input_width,
        height = UI_STYLE.input_height,
        buttoncolor = :lightgrey,
        buttoncolor_hover = :grey,
        buttoncolor_active = :green,
    )
end

"""Create a styled textbox for GUI text input."""
function _create_textbox(parent; width = nothing, placeholder = "", boxcolor_obs = nothing, halign = :center)
    bc = isnothing(boxcolor_obs) ? Observable(:white) : boxcolor_obs
    Textbox(
        parent,
        placeholder = placeholder,
        fontsize = UI_STYLE.textbox_font,
        width = isnothing(width) ? UI_STYLE.input_width : width,
        height = UI_STYLE.input_height,
        halign = halign,
        boxcolor = bc,
        cornerradius = 0,
    )
end

"""Create a styled dropdown menu for GUI option selection."""
function _create_menu(parent; options = ["Select"], width = nothing, height = nothing)
    Menu(
        parent,
        options = options,
        width = isnothing(width) ? UI_STYLE.input_width : width,
        height = isnothing(height) ? UI_STYLE.input_height : height,
        fontsize = UI_STYLE.textbox_font,
    )
end

"""Create a styled label for the GUI."""
function _create_label(parent, text; fontsize = UI_STYLE.label_font, color = :black)
    Label(parent, text, fontsize = fontsize, color = color)
end

"""Truncate a path string to a maximum length, showing the end with '...' prefix."""
function _truncate_path(path::String, max_length::Int = 20)
    if length(path) <= max_length
        return path
    end
    # Show last (max_length - 3) characters with "..." prefix (3 chars for "...")
    return "..." * path[(end-(max_length-4)):end]
end

"""
    plot_gui()

Interactive GUI for quick data plotting and visualization.

# GUI Structure

**Column 1: File & Plot Selection**
- Directory and file browser
- Plot type dropdown with hierarchical submenus:
  - Common plots: Data Browser, Epochs, ERP, ERP Image, Topography
  - Time-Frequency >: Time-Frequency, Power Spectrum
  - ICA > : Components, Activation
  - ERP Analysis > : ERP GUI
  - Diagnostic Plots >: Artifact Detection, Triggers, Channel Summary,
    Joint Probability, Correlation Heatmap, Global Field Power,
    Layout View, Filter Response

**Column 2: Data Selection & Layout**
- Participant, Condition, Epoch filters (integers, space/comma-separated)
- Layout type: Single, Single Avg, Grid, Topo
- Channel selection with multi-select support

**Column 3: Axis Settings**
- X, Y, Z axis limits (validated numeric input with range checking)
- Baseline correction window and type
- Display options (invert Y-axis)

# Supported Plot Types

See EegFun.jl documentation for complete plot type details and data requirements.

# Examples

```julia
# Launch GUI
EegFun.plot_gui()

# Navigate to data directory, select file and plot type
# Configure channels and axis limits as needed
# Click "Plot" to generate visualization
```
"""
function plot_gui()
    _set_window_title("PLOT GUI")
    # Dynamic sizing using Auto() rather than strict Fixed everywhere
    gui_fig = Figure(size = (700, 600), title = "Plot GUI", backgroundcolor = :lightgrey, figure_padding = (20, 20, 20, 20))
    main_layout = GridLayout(gui_fig[1, 1:3], rowgap = 2, colgap = 4)


    # Pre-emptively create the channel_menu so it can live in GUIState
    channel_menu = _create_menu(main_layout[10, 2], options = ["Select"])

    gui_state = GUIState(; channel_menu = channel_menu)

    _build_file_column!(main_layout, gui_state)
    _build_filter_column!(main_layout, gui_state)
    _build_axis_column!(main_layout, gui_state)

    # Set columns sizes
    colsize!(main_layout, 1, Fixed(UI_STYLE.input_width))  # pin col 1 so widgets don't drift on resize
    colsize!(main_layout, 2, Auto())
    colsize!(main_layout, 3, Auto())

    display(gui_fig)
    _set_window_title("Makie")

    return nothing
end

function _build_file_column!(main_layout, gui_state)
    # Select Directory Section
    directory_select_button = _create_select_button(main_layout[1, 1], "Select Directory")
    directory_label_text = Observable("Dir:")
    _create_label(main_layout[2, 1], directory_label_text, fontsize = UI_STYLE.textbox_font, color = :gray)

    # Select File Section
    file_select_button = _create_select_button(main_layout[3, 1], "Select File")
    file_pattern_input = _create_textbox(main_layout[4, 1], placeholder = "", halign = :left)

    # Layout Section
    layout_select_button = _create_select_button(main_layout[5, 1], "Select Layout")
    layout_label_text = Observable("File:")
    _create_label(main_layout[6, 1], layout_label_text, fontsize = UI_STYLE.textbox_font, color = :gray)

    # Plot Type Section
    _create_label(main_layout[7, 1], "Plot Type")
    plottype_options = [
        "Select",
        "Data Browser",
        "Epochs",
        "ERP",
        "Topography",
        "ERP Image",
        "Time-Frequency",
        "─────────────────",
        "ICA >",
        "ERP Analysis >",
        "Diagnostic Plots >",
    ]
    plottype_dropdown = _create_menu(main_layout[8, 1], options = plottype_options)

    # Submenu dropdown
    submenu_dropdown = _create_menu(main_layout[9, 1], options = ["---"])

    component_bc = Observable(:lightgrey)
    component_input = _create_textbox(main_layout[10, 1], placeholder = "", boxcolor_obs = component_bc)

    # Events
    on(directory_select_button.clicks) do _
        default_path = gui_state.directory[] != "" ? gui_state.directory[] : pwd()
        dir_path = pick_folder(default_path)
        if !isnothing(dir_path) && dir_path != ""
            gui_state.directory[] = dir_path
            directory_label_text[] = "." * _truncate_path(String(strip(dir_path)))
        end
    end

    on(file_select_button.clicks) do _
        default_path = gui_state.directory[] != "" ? gui_state.directory[] : pwd()
        filename = pick_file(default_path)
        if !isnothing(filename) && filename != ""
            file_pattern_input.displayed_string[] = ""
            file_pattern_input.displayed_string[] = basename(filename)
            gui_state.filename[] = filename
        end
    end

    last_validated_file = Ref{String}("__uninit__")
    function validate_file_input(value)
        val = strip(value)
        if val == last_validated_file[]
            return
        end
        last_validated_file[] = val
        if !isempty(val) && val != basename(gui_state.filename[])
            gui_state.filename[] = val
        end
    end
    on(file_pattern_input.stored_string) do value
        validate_file_input(value)
    end
    on(file_pattern_input.focused) do is_focused
        if !is_focused
            validate_file_input(file_pattern_input.displayed_string[])
        end
    end

    on(layout_select_button.clicks) do _
        default_path = gui_state.directory[] != "" ? gui_state.directory[] : pwd()
        filename = pick_file(default_path)
        if !isnothing(filename) && filename != ""
            basename_only = basename(filename)
            gui_state.layout_file[] = filename
            gui_state.layout[] = basename_only
            layout_label_text[] = basename_only
            try
                layout = read_layout(filename)
                gui_state.layout_object[] = layout
                gui_state.channel_menu.options = vcat(["Select"], string.(channel_labels(layout)))
                @info "Loaded layout with $(length(channel_labels(layout))) channels"
            catch layout_error
                @minimal_warning "Error loading layout: $layout_error"
                gui_state.layout_object[] = nothing
            end
        end
    end

    updating_submenu = Ref{Bool}(false)
    on(plottype_dropdown.selection) do selection
        if selection == "Select" || selection == "─────────────────"
            return
        end
        updating_submenu[] = true
        try
            if selection == "ICA >"
                submenu_dropdown.options = ["Select", "Components", "Activation"]
                gui_state.submenu_active[] = true
                gui_state.submenu_type[] = "ica"
                component_bc[] = :white
            elseif selection == "ERP Analysis >"
                submenu_dropdown.options = ["Select", "ERP GUI"]
                gui_state.submenu_active[] = true
                gui_state.submenu_type[] = "erp_analysis"
                component_bc[] = :lightgrey
            elseif selection == "Diagnostic Plots >"
                submenu_dropdown.options = ["Select", "Channel Summary", "Layout View"]
                gui_state.submenu_active[] = true
                gui_state.submenu_type[] = "diagnostic"
                component_bc[] = :lightgrey
            else
                gui_state.plottype[] = selection
                gui_state.submenu_active[] = false
                component_bc[] = :lightgrey
            end
        finally
            updating_submenu[] = false
        end
    end

    on(submenu_dropdown.selection) do selection
        if updating_submenu[] || isnothing(selection) || selection == "Select"
            return
        end
        gui_state.plottype[] = selection
    end

    last_validated_component = Ref{String}("__uninit__")
    function validate_component_input(value)
        val = strip(value)
        if val == last_validated_component[]
            return
        end
        last_validated_component[] = val
        gui_state.ica_components[] = val
    end
    on(component_input.stored_string) do value
        validate_component_input(value)
    end
    on(component_input.focused) do is_focused
        if !is_focused
            validate_component_input(component_input.displayed_string[])
        end
    end
end

function _build_filter_column!(main_layout, gui_state)
    participant_bc = Observable(:white)
    condition_bc = Observable(:white)
    epoch_bc = Observable(:white)

    _create_label(main_layout[1, 2], "Participant")
    participant_input = _create_textbox(main_layout[2, 2], boxcolor_obs = participant_bc)
    _create_label(main_layout[3, 2], "Condition")
    condition_input = _create_textbox(main_layout[4, 2], boxcolor_obs = condition_bc)
    _create_label(main_layout[5, 2], "Epoch")
    epoch_input = _create_textbox(main_layout[6, 2], boxcolor_obs = epoch_bc)

    _create_label(main_layout[7, 2], "Layout")
    layout_dropdown = _create_menu(main_layout[8, 2], options = ["Single", "Single Avg", "Grid", "Topo"])

    _create_label(main_layout[9, 2], "Channel(s)")
    # channel_menu is already in layout at [10, 2] by plot_gui setup

    selected_channels_text = Observable("Channels: ")
    _create_label(main_layout[11, 2], selected_channels_text, fontsize = UI_STYLE.textbox_font, color = :gray)

    # Connections
    setup_multi_int_callback(participant_input, participant_bc, "Participant", gui_state.participant)
    setup_multi_int_callback(condition_input, condition_bc, "Condition", gui_state.condition)
    setup_multi_int_callback(epoch_input, epoch_bc, "Epoch", gui_state.epoch)

    on(layout_dropdown.selection) do selection
        if !isnothing(selection)
            if selection == "Single Avg"
                gui_state.layout_type[] = "single"
                gui_state.average_channels[] = true
            else
                gui_state.layout_type[] = lowercase(selection)
                gui_state.average_channels[] = false
            end
        end
    end

    on(gui_state.channel_menu.selection) do selection
        if selection == "Select"
            gui_state.electrodes[] = String[]
        elseif selection in gui_state.electrodes[]
            gui_state.electrodes[] = filter(!=(selection), gui_state.electrodes[])
        else
            gui_state.electrodes[] = vcat(gui_state.electrodes[], selection)
        end
        if isempty(gui_state.electrodes[])
            selected_channels_text[] = "Selected: "
        else
            selected_channels_text[] = "Selected: " * join(gui_state.electrodes[], ", ")
        end
        gui_state.channel_menu.i_selected.val = 0
        gui_state.channel_menu.selection.val = nothing
    end
end

function _build_axis_column!(main_layout, gui_state)
    xmin_bc = Observable(:white)
    xmax_bc = Observable(:white)
    ymin_bc = Observable(:white)
    ymax_bc = Observable(:white)
    zmin_bc = Observable(:white)
    zmax_bc = Observable(:white)
    bl_start_bc = Observable(:white)
    bl_end_bc = Observable(:white)

    _create_label(main_layout[1, 3], "X Limits", fontsize = UI_STYLE.textbox_font)
    x_limits_layout = GridLayout(main_layout[2, 3], tellwidth = false, colgap = 8)
    xmin_input = _create_textbox(x_limits_layout[1, 1], width = 80, boxcolor_obs = xmin_bc)
    xmax_input = _create_textbox(x_limits_layout[1, 2], width = 80, boxcolor_obs = xmax_bc)

    _create_label(main_layout[3, 3], "Y Limits", fontsize = UI_STYLE.textbox_font)
    y_limits_layout = GridLayout(main_layout[4, 3], tellwidth = false, colgap = 8)
    ymin_input = _create_textbox(y_limits_layout[1, 1], width = 80, boxcolor_obs = ymin_bc)
    ymax_input = _create_textbox(y_limits_layout[1, 2], width = 80, boxcolor_obs = ymax_bc)

    _create_label(main_layout[5, 3], "Z Limits", fontsize = UI_STYLE.textbox_font)
    z_limits_layout = GridLayout(main_layout[6, 3], tellwidth = false, colgap = 8)
    zmin_input = _create_textbox(z_limits_layout[1, 1], width = 80, boxcolor_obs = zmin_bc)
    zmax_input = _create_textbox(z_limits_layout[1, 2], width = 80, boxcolor_obs = zmax_bc)

    _create_label(main_layout[7, 3], "Baseline", fontsize = UI_STYLE.textbox_font)
    baseline_layout = GridLayout(main_layout[8, 3], tellwidth = false, colgap = 8)
    baseline_start = _create_textbox(baseline_layout[1, 1], width = 80, boxcolor_obs = bl_start_bc)
    baseline_end = _create_textbox(baseline_layout[1, 2], width = 80, boxcolor_obs = bl_end_bc)

    _create_label(main_layout[9, 3], "Baseline Type TF", fontsize = UI_STYLE.textbox_font)
    baseline_type = _create_menu(
        main_layout[10, 3],
        options = ["Select", "db", "absolute", "relative", "relchange", "percent", "zscore", "normchange"],
    )

    invert_y_layout = GridLayout(main_layout[11, 3], tellwidth = false, colgap = 8)
    _create_label(invert_y_layout[1, 1], "Invert Y axis", fontsize = UI_STYLE.textbox_font)
    invert_y_checkbox = Checkbox(invert_y_layout[1, 2], checked = false)

    plot_button = Button(
        main_layout[13, 3],
        label = "PLOT",
        fontsize = UI_STYLE.button_font,
        buttoncolor = :darkgrey,
        buttoncolor_hover = :grey,
        buttoncolor_active = :green,
        width = 100,
        height = 50,
    )

    # Callbacks
    x_check = () -> check_range(gui_state.xlim, xmin_bc, xmax_bc, "X Limits")
    y_check = () -> check_range(gui_state.ylim, ymin_bc, ymax_bc, "Y Limits")
    z_check = () -> check_range(gui_state.zlim, zmin_bc, zmax_bc, "Z Limits")
    bl_check = () -> begin
        s, e = gui_state.baseline_start[], gui_state.baseline_end[]
        if !isnothing(s) && !isnothing(e) && s > e
            bl_start_bc[] = :lightcoral
            bl_end_bc[] = :lightcoral
            @minimal_warning "Invalid Baseline range: start ($s) must be less than end ($e)"
        end
    end

    setup_numeric_callback(xmin_input, xmin_bc, "X min", v -> gui_state.xlim[] = (v, gui_state.xlim[][2]); range_check = x_check)
    setup_numeric_callback(xmax_input, xmax_bc, "X max", v -> gui_state.xlim[] = (gui_state.xlim[][1], v); range_check = x_check)
    setup_numeric_callback(ymin_input, ymin_bc, "Y min", v -> gui_state.ylim[] = (v, gui_state.ylim[][2]); range_check = y_check)
    setup_numeric_callback(ymax_input, ymax_bc, "Y max", v -> gui_state.ylim[] = (gui_state.ylim[][1], v); range_check = y_check)
    setup_numeric_callback(zmin_input, zmin_bc, "Z min", v -> gui_state.zlim[] = (v, gui_state.zlim[][2]); range_check = z_check)
    setup_numeric_callback(zmax_input, zmax_bc, "Z max", v -> gui_state.zlim[] = (gui_state.zlim[][1], v); range_check = z_check)
    setup_numeric_callback(baseline_start, bl_start_bc, "Baseline start", v -> gui_state.baseline_start[] = v; range_check = bl_check)
    setup_numeric_callback(baseline_end, bl_end_bc, "Baseline end", v -> gui_state.baseline_end[] = v; range_check = bl_check)

    on(baseline_type.selection) do selection
        if !isnothing(selection) && selection != "Select"
            gui_state.baseline_type[] = selection
        end
    end

    on(invert_y_checkbox.checked) do is_checked
        gui_state.invert_y[] = is_checked
    end

    plot_dispatch = Dict{String,Function}(
        "Data Browser"    => gs -> _plot_databrowser(gs),
        "Epochs"          => gs -> _plot_epochs(gs),
        "ERP"             => gs -> _plot_erp(gs),
        "ERP Image"       => gs -> _plot_erp_image(gs),
        "Topography"      => gs -> _plot_topography(gs),
        "Time-Frequency"  => gs -> _plot_time_frequency(gs),
        "Components"      => gs -> _plot_ica(gs),
        "Activation"      => gs -> _plot_ica_activation(gs),
        "ERP GUI"         => gs -> _plot_erp_measurement_gui(gs),
        "Channel Summary" => gs -> _plot_channel_summary(gs),
        "Layout View"     => gs -> _plot_layout(gs),
    )

    on(plot_button.clicks) do _
        plot_type = gui_state.plottype[]
        if plot_type == "Select" || startswith(plot_type, "─") || plot_type == "select"
            return
        end
        handler = get(plot_dispatch, plot_type, nothing)
        if !isnothing(handler)
            handler(gui_state)
        else
            @minimal_warning "Unsupported plot type: $plot_type"
        end
    end
end




"""GUI bridge: load data and open the databrowser."""
function _plot_databrowser(gui_state)

    # Check if we have the required files
    if gui_state.filename[] == ""
        @minimal_warning "Requested plot settings incompatible: recheck!"
        return
    end

    file_ext = lowercase(splitext(gui_state.filename[])[2])
    if file_ext ∉ [".bdf", ".jld2"]
        @minimal_warning "Requested plot settings incompatible: recheck!"
        return
    end

    try
        if file_ext == ".jld2"
            # Screen must be created for single datasets to open in a separate window.
            # For Vectors, the Vector dispatch handles its own screens — creating one
            # upfront would produce a spurious empty window.
            @async begin
                dat = read_data(gui_state.filename[])
                isnothing(dat) && return
                dat = _apply_gui_filters(dat, gui_state)
                isnothing(dat) && return
                if dat isa Vector
                    plot_databrowser(dat)
                else
                    new_screen = GLMakie.Screen()
                    plot_databrowser(dat; screen = new_screen)
                end
            end
        else # Load from raw eeg data file
            layout = nothing
            if gui_state.layout_file[] != ""
                layout = read_layout(gui_state.layout_file[])
                polar_to_cartesian_xy!(layout)
            end

            dat = read_raw_data(gui_state.filename[])
            if isnothing(layout)
                dat = create_eegfun_data(dat)
            else
                dat = create_eegfun_data(dat, layout)
            end

            # Update electrode menu with actual channel labels from the loaded data
            gui_state.channel_menu.options = vcat(["Select"], string.(channel_labels(dat)))
            @async begin
                new_screen = GLMakie.Screen()
                plot_databrowser(dat; screen = new_screen)
            end
        end
    catch e
        _handle_plot_error(e, "Data Browser")
    end
end

"""GUI bridge: load epoch data and call `plot_epochs`."""
function _plot_epochs(gui_state)
    isnothing(_validate_file(gui_state)) && return

    fname = gui_state.filename[]
    is_pattern = _is_pattern(fname)

    if is_pattern
        # Pattern mode — parse participant/condition filters and forward to plot_epochs
        cond_list = _parse_gui_int_field(gui_state.condition[])
        if isnothing(cond_list)
            @minimal_warning "Requested plot settings incompatible: recheck!"
            return
        end
        cond_sel = isempty(cond_list) ? conditions() : conditions(cond_list)
        selected_channels = isempty(gui_state.electrodes[]) ? channels() : channels(Symbol.(gui_state.electrodes[]))

        res = _resolve_batch_files_helper(gui_state)
        isnothing(res) && return
        files, part_sel, input_dir_resolved = res

        @async begin
            layout_sym = Symbol(gui_state.layout_type[])
            try
                plot_epochs(
                    fname;
                    input_dir             = _gui_dir(gui_state),
                    participant_selection = part_sel,
                    condition_selection   = cond_sel,
                    channel_selection     = selected_channels,
                    layout                = layout_sym,
                    xlim                  = gui_state.xlim[],
                    ylim                  = gui_state.ylim[],
                    yreversed             = gui_state.invert_y[],
                )
            catch e
                _handle_plot_error(e, "Epochs")
            end
        end
        return
    end

    # .jld2 mode — apply condition + epoch filters from GUI
    try
        cond_list  = _parse_gui_int_field(gui_state.condition[])
        epoch_list = _parse_gui_int_field(gui_state.epoch[])
        if isnothing(cond_list) || isnothing(epoch_list)
            @minimal_warning "Requested plot settings incompatible: recheck!"
            return
        end
        cond_sel  = isempty(cond_list) ? conditions() : conditions(cond_list)
        epoch_sel = isempty(epoch_list) ? epochs() : epochs(epoch_list)

        selected_channels = isempty(gui_state.electrodes[]) ? channels() : channels(Symbol.(gui_state.electrodes[]))

        @async begin
            layout_sym = Symbol(gui_state.layout_type[])
            try
                plot_epochs(
                    fname;
                    channel_selection = selected_channels,
                    condition_selection = cond_sel,
                    epoch_selection = epoch_sel,
                    layout = layout_sym,
                    average_channels = gui_state.average_channels[],
                    xlim = gui_state.xlim[],
                    ylim = gui_state.ylim[],
                    yreversed = gui_state.invert_y[],
                )
            catch e
                _handle_plot_error(e, "Epochs")
            end
        end
    catch e
        _handle_plot_error(e, "Epochs")
    end
end

"""GUI bridge: load ERP data and call `plot_erp`."""
function _plot_erp(gui_state)
    isnothing(_validate_file(gui_state)) && return

    fname = gui_state.filename[]
    input_dir = _gui_dir(gui_state)
    is_pattern = _is_pattern(fname)

    if is_pattern
        # Pattern mode — parse all relevant GUI filters and forward to plot_erp
        cond_list = _parse_gui_int_field(gui_state.condition[])
        if isnothing(cond_list)
            @minimal_warning "Requested plot settings incompatible: recheck!"
            return
        end
        cond_sel = isempty(cond_list) ? conditions() : conditions(cond_list)
        selected_channels = isempty(gui_state.electrodes[]) ? channels() : channels(Symbol.(gui_state.electrodes[]))

        res = _resolve_batch_files_helper(gui_state)
        isnothing(res) && return
        files, part_sel, input_dir_resolved = res

        @async begin
            layout_sym = Symbol(gui_state.layout_type[])
            try
                baseline = nothing
                if !isnothing(gui_state.baseline_start[]) && !isnothing(gui_state.baseline_end[])
                    baseline = (gui_state.baseline_start[], gui_state.baseline_end[])
                end
                plot_erp(
                    fname;
                    input_dir = input_dir,
                    participant_selection = part_sel,
                    layout = layout_sym,
                    channel_selection = selected_channels,
                    condition_selection = cond_sel,
                    baseline_interval = baseline,
                    average_channels = gui_state.average_channels[],
                    xlim = gui_state.xlim[],
                    ylim = gui_state.ylim[],
                    yreversed = gui_state.invert_y[],
                )
            catch e
                _handle_plot_error(e, "ERP")
            end
        end
        return
    end

    # .jld2 mode — apply condition filter from GUI
    try
        cond_list = _parse_gui_int_field(gui_state.condition[])
        if isnothing(cond_list)
            @minimal_warning "Requested plot settings incompatible: recheck!"
            return
        end
        cond_sel = isempty(cond_list) ? conditions() : conditions(cond_list)

        selected_channels = isempty(gui_state.electrodes[]) ? channels() : channels(Symbol.(gui_state.electrodes[]))

        baseline = nothing
        if !isnothing(gui_state.baseline_start[]) && !isnothing(gui_state.baseline_end[])
            baseline = (gui_state.baseline_start[], gui_state.baseline_end[])
        end

        @async begin
            layout_sym = Symbol(gui_state.layout_type[])
            try
                plot_erp(
                    fname;
                    layout = layout_sym,
                    condition_selection = cond_sel,
                    channel_selection = selected_channels,
                    baseline_interval = baseline,
                    average_channels = gui_state.average_channels[],
                    xlim = gui_state.xlim[],
                    ylim = gui_state.ylim[],
                    yreversed = gui_state.invert_y[],
                )
            catch e
                _handle_plot_error(e, "ERP")
            end
        end
    catch e
        _handle_plot_error(e, "ERP")
    end
end

"""GUI bridge: load data and call `plot_topography`."""
function _plot_topography(gui_state)
    isnothing(_validate_file(gui_state)) && return

    fname = gui_state.filename[]
    input_dir = _gui_dir(gui_state)
    is_pattern = _is_pattern(fname)

    # Shared setup
    selected_channels = isempty(gui_state.electrodes[]) ? channels() : channels(Symbol.(gui_state.electrodes[]))

    function _do_topo_plot(data)
        data = _apply_gui_filters(data, gui_state)
        isnothing(data) && return

        # Build baseline interval from GUI
        b_lo, b_hi = gui_state.baseline_start[], gui_state.baseline_end[]
        baseline_interval = (!isnothing(b_lo) && !isnothing(b_hi)) ? (b_lo, b_hi) : nothing

        !isnothing(baseline_interval) && baseline!(data, baseline_interval)

        # xlim → interval_selection (time window to average over); ylim → color scale
        x_lo, x_hi   = gui_state.xlim[]
        interval_sel = (!isnothing(x_lo) && !isnothing(x_hi)) ? times(x_lo, x_hi) : times()
        y_lo, y_hi   = gui_state.ylim[]
        color_range  = (!isnothing(y_lo) && !isnothing(y_hi)) ? (y_lo, y_hi) : nothing
        extra_kwargs = isnothing(color_range) ? (;) : (ylim = color_range,)

        if data isa Vector{<:ErpData} || data isa ErpData
            plot_topography(data; channel_selection = selected_channels, interval_selection = interval_sel, extra_kwargs...)
        elseif data isa Vector{<:EpochData} || data isa EpochData
            epoch_num = try
                epoch_str = strip(gui_state.epoch[])
                isempty(epoch_str) ? 1 : parse(Int, epoch_str)
            catch e
                @debug "Failed to parse epoch number, defaulting to 1" exception = e
                1
            end
            plot_topography(data, epoch_num; channel_selection = selected_channels, interval_selection = interval_sel, extra_kwargs...)
        else
            @minimal_warning "Requested plot settings incompatible: recheck!"
        end
    end

    if is_pattern
        cond_list = _parse_gui_int_field(gui_state.condition[])
        if isnothing(cond_list)
            @minimal_warning "Requested plot settings incompatible: recheck!"
            return
        end

        res = _resolve_batch_files_helper(gui_state)
        isnothing(res) && return
        files, part_sel, input_dir = res

        @async begin
            try
                for file in sort(files, by = _natural_sort_key)
                    @info "Plotting: $file"
                    data = read_data(joinpath(input_dir, file))
                    isnothing(data) && continue
                    _do_topo_plot(data)
                end
            catch e
                _handle_plot_error(e, "Topography")
            end
        end
        return
    end

    # .jld2 mode
    try
        data = read_data(fname)
        if isnothing(data)
            @minimal_warning "Requested plot settings incompatible: recheck!"
            return
        end
        @async begin
            try
                _do_topo_plot(data)
            catch e
                _handle_plot_error(e, "Topography")
            end
        end
    catch e
        _handle_plot_error(e, "Topography")
    end
end
