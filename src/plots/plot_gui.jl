# Base UI styling values
const BASE_FONTS = (label = 20, button = 18, textbox = 16)
const BASE_SIZES = (input_width = 160, input_height = 30)

"""Font sizes and widget dimensions for the interactive GUI."""
struct UIStyle
    label_font::Int
    button_font::Int
    textbox_font::Int
    input_width::Int
    input_height::Int

    UIStyle() = new(BASE_FONTS.label, BASE_FONTS.button, BASE_FONTS.textbox, BASE_SIZES.input_width, BASE_SIZES.input_height)
end

"""Create a styled button for GUI selection actions."""
function _create_select_button(parent, label, style::UIStyle)
    Button(
        parent,
        label = label,
        fontsize = style.button_font,
        width = style.input_width,
        height = style.input_height,
        buttoncolor = :lightgrey,
        buttoncolor_hover = :grey,
        buttoncolor_active = :green,
    )
end

"""Create a styled textbox for GUI text input."""
function _create_textbox(parent, style::UIStyle; width = nothing, placeholder = "", boxcolor_obs = nothing, halign = :center)
    bc = isnothing(boxcolor_obs) ? Observable(:white) : boxcolor_obs
    Textbox(
        parent,
        placeholder = placeholder,
        fontsize = style.textbox_font,
        width = isnothing(width) ? style.input_width : width,
        height = style.input_height,
        halign = halign,
        boxcolor = bc,
        cornerradius = 0,
    )
end

"""Create a styled dropdown menu for GUI option selection."""
function _create_menu(parent, style::UIStyle; options = ["Select"], width = nothing, height = nothing)
    Menu(
        parent,
        options = options,
        width = isnothing(width) ? style.input_width : width,
        height = isnothing(height) ? style.input_height : height,
        fontsize = style.textbox_font,
    )
end

"""Create a styled label for the GUI."""
function _create_label(parent, text, style::UIStyle; fontsize = nothing, color = :black)
    fontsize_val = isnothing(fontsize) ? style.label_font : fontsize
    Label(parent, text, fontsize = fontsize_val, color = color)
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

    # main figure window, layout, and UI style
    _set_window_title("PLOT GUI")
    gui_fig = Figure(size = (650, 650), title = "Plot GUI", backgroundcolor = :lightgrey, figure_padding = (20, 20, 20, 20))
    main_layout = GridLayout(gui_fig[1, 1:3], rowgap = 2, colgap = 4)
    ui_style = UIStyle()
    colsize!(main_layout, 1, Fixed(ui_style.input_width))  # pin col 1 so widgets don't drift on resize

    # Select Directory Section
    directory_select_button = _create_select_button(main_layout[1, 1], "Select Directory", ui_style)
    directory_label_text = Observable("Dir:")
    _create_label(main_layout[2, 1], directory_label_text, ui_style, fontsize = ui_style.textbox_font, color = :gray)

    # Select File Section — Browse button + editable textbox (type path or pattern)
    file_select_button = _create_select_button(main_layout[3, 1], "Select File", ui_style)
    file_pattern_input = _create_textbox(main_layout[4, 1], ui_style, placeholder = "", halign = :left)

    # Layout Section
    layout_select_button = _create_select_button(main_layout[5, 1], "Select Layout", ui_style)
    layout_label_text = Observable("File:")
    _create_label(main_layout[6, 1], layout_label_text, ui_style, fontsize = ui_style.textbox_font, color = :gray)

    # Plot Type Section
    _create_label(main_layout[7, 1], "Plot Type", ui_style)
    plottype_options = [
        "Select",
        "Data Browser",
        "Epochs",
        "ERP",
        "Topography",
        "ERP Image",
        "Time-Frequency",
        "─────────────────",  # Visual separator
        "ICA >",
        "ERP Analysis >",
        "Diagnostic Plots >",
    ]
    plottype_dropdown = _create_menu(main_layout[8, 1], ui_style, options = plottype_options)

    # Submenu dropdown (starts with placeholder to appear inactive)
    submenu_dropdown = _create_menu(main_layout[9, 1], ui_style, options = ["---"])

    # Component input textbox — "hidden" (lightgrey) by default, revealed (white) when ICA is active.
    component_bc = Observable(:lightgrey)
    component_input = _create_textbox(main_layout[10, 1], ui_style, placeholder = "", boxcolor_obs = component_bc)
    # Column 2: Participant, Condition, Epoch, Channels
    # Create boxcolor Observables for validated inputs
    participant_bc = Observable(:white)
    condition_bc = Observable(:white)
    epoch_bc = Observable(:white)

    # Participant Section
    _create_label(main_layout[1, 2], "Participant", ui_style)
    participant_input = _create_textbox(main_layout[2, 2], ui_style, boxcolor_obs = participant_bc)

    # Condition Section
    _create_label(main_layout[3, 2], "Condition", ui_style)
    condition_input = _create_textbox(main_layout[4, 2], ui_style, boxcolor_obs = condition_bc)

    # Epoch Section
    _create_label(main_layout[5, 2], "Epoch", ui_style)
    epoch_input = _create_textbox(main_layout[6, 2], ui_style, boxcolor_obs = epoch_bc)

    # Layout Section (replaces Channel(s) and Average Channels)
    _create_label(main_layout[7, 2], "Layout", ui_style)
    layout_dropdown = _create_menu(main_layout[8, 2], ui_style, options = ["Single", "Single Avg", "Grid", "Topo"])

    # Channel(s) Section (moved down)
    _create_label(main_layout[9, 2], "Channel(s)", ui_style)
    channel_menu = _create_menu(main_layout[10, 2], ui_style, options = ["Select"])

    # Selected channels display
    selected_channels_text = Observable("Channels: ")
    _create_label(main_layout[11, 2], selected_channels_text, ui_style, fontsize = ui_style.textbox_font, color = :gray)

    # Column 3: Axis Settings
    # Create boxcolor Observables for validated inputs
    xmin_bc = Observable(:white)
    xmax_bc = Observable(:white)
    ymin_bc = Observable(:white)
    ymax_bc = Observable(:white)
    zmin_bc = Observable(:white)
    zmax_bc = Observable(:white)
    bl_start_bc = Observable(:white)
    bl_end_bc = Observable(:white)

    # X Limits Section
    _create_label(main_layout[1, 3], "X Limits", ui_style, fontsize = ui_style.textbox_font)
    x_limits_layout = GridLayout(main_layout[2, 3], tellwidth = false, colgap = 8)
    xmin_input = _create_textbox(x_limits_layout[1, 1], ui_style, width = 80, boxcolor_obs = xmin_bc)
    xmax_input = _create_textbox(x_limits_layout[1, 2], ui_style, width = 80, boxcolor_obs = xmax_bc)

    # Y Limits Section
    _create_label(main_layout[3, 3], "Y Limits", ui_style, fontsize = ui_style.textbox_font)
    y_limits_layout = GridLayout(main_layout[4, 3], tellwidth = false, colgap = 8)
    ymin_input = _create_textbox(y_limits_layout[1, 1], ui_style, width = 80, boxcolor_obs = ymin_bc)
    ymax_input = _create_textbox(y_limits_layout[1, 2], ui_style, width = 80, boxcolor_obs = ymax_bc)

    # Z Limits Section (currently unused, but kept for future use)
    _create_label(main_layout[5, 3], "Z Limits", ui_style, fontsize = ui_style.textbox_font)
    z_limits_layout = GridLayout(main_layout[6, 3], tellwidth = false, colgap = 8)
    zmin_input = _create_textbox(z_limits_layout[1, 1], ui_style, width = 80, boxcolor_obs = zmin_bc)
    zmax_input = _create_textbox(z_limits_layout[1, 2], ui_style, width = 80, boxcolor_obs = zmax_bc)

    # Baseline Section
    _create_label(main_layout[7, 3], "Baseline", ui_style, fontsize = ui_style.textbox_font)
    baseline_layout = GridLayout(main_layout[8, 3], tellwidth = false, colgap = 8)
    baseline_start = _create_textbox(baseline_layout[1, 1], ui_style, width = 80, boxcolor_obs = bl_start_bc)
    baseline_end = _create_textbox(baseline_layout[1, 2], ui_style, width = 80, boxcolor_obs = bl_end_bc)

    # Baseline Type
    _create_label(main_layout[9, 3], "Baseline Type TF", ui_style, fontsize = ui_style.textbox_font)
    baseline_type = _create_menu(
        main_layout[10, 3],
        ui_style,
        options = ["Select", "db", "absolute", "relative", "relchange", "percent", "zscore", "normchange"],
    )

    # Invert Y axis option
    invert_y_layout = GridLayout(main_layout[11, 3], tellwidth = false, colgap = 8)
    _create_label(invert_y_layout[1, 1], "Invert Y axis", ui_style, fontsize = ui_style.textbox_font)
    invert_y_checkbox = Checkbox(invert_y_layout[1, 2], checked = false)

    # Main plot button
    plot_button = Button(
        main_layout[13, 3],
        label = "PLOT",
        fontsize = ui_style.button_font,
        buttoncolor = :darkgrey,
        buttoncolor_hover = :grey,
        buttoncolor_active = :green,
        width = 100,
        height = 50,
    )

    # Set columns sizes
    colsize!(main_layout, 1, Relative(0.3))
    colsize!(main_layout, 2, Relative(0.3))
    colsize!(main_layout, 3, Relative(0.3))

    # Data structure to store GUI state
    gui_state = (
        directory = Observable(""),
        filename = Observable(""),
        participant = Observable(""),
        condition = Observable(""),
        epoch = Observable(""),
        plottype = Observable("select"),
        layout = Observable("select"),
        layout_file = Observable(""),
        layout_object = Observable{Any}(nothing),
        electrodes = Observable(String[]),
        xlim = Observable{Tuple{Union{Nothing,Float64},Union{Nothing,Float64}}}((nothing, nothing)),
        ylim = Observable{Tuple{Union{Nothing,Float64},Union{Nothing,Float64}}}((nothing, nothing)),
        zlim = Observable{Tuple{Union{Nothing,Float64},Union{Nothing,Float64}}}((nothing, nothing)),
        baseline_start = Observable{Union{Nothing,Float64}}(nothing),
        baseline_end = Observable{Union{Nothing,Float64}}(nothing),
        baseline_type = Observable("select"),
        layout_type = Observable("single"),
        average_channels = Observable(false),
        invert_y = Observable(false),
        submenu_active = Observable(false),
        submenu_type = Observable(""),
        channel_menu = channel_menu,
        ica_components = Observable(""),
    )

    # Plot type dispatch table
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

    function plot()
        plot_type = gui_state.plottype[]
        # Skip placeholder and separator options
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

    # Helper function for numeric input validation with visual feedback
    function setup_numeric_callback(input, bc_obs::Observable, field_name::String, update_fn::Function; range_check = nothing)
        last_validated = Ref{String}("__uninit__")

        function validate(value)
            stripped = strip(value)
            if stripped == last_validated[]
                return
            end
            last_validated[] = stripped

            if isempty(stripped)
                bc_obs[] = :white
                update_fn(nothing)
            else
                parsed = tryparse(Float64, stripped)
                if isnothing(parsed)
                    bc_obs[] = :lightcoral
                    @minimal_warning "Invalid $field_name value: \"$stripped\" (must be a number)"
                else
                    bc_obs[] = :white
                    update_fn(parsed)
                end
            end
            # Run range check after updating value
            !isnothing(range_check) && range_check()
        end

        # Validate on Enter
        on(input.stored_string) do value
            validate(value)
        end

        # Validate on defocus (clicking away)
        on(input.focused) do is_focused
            if !is_focused
                validate(input.displayed_string[])
            end
        end
    end

    # Range validation: if both values are set and min >= max, both boxes turn red
    function check_range(limit_obs, min_bc, max_bc, label)
        lo, hi = limit_obs[]
        if !isnothing(lo) && !isnothing(hi) && lo > hi
            min_bc[] = :lightcoral
            max_bc[] = :lightcoral
            @minimal_warning "Invalid $label range: min ($lo) must be less than or equal to max ($hi)"
        end
    end

    # Set up X/Y/Z limit callbacks with range checks
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

    # Set up baseline callbacks with range check
    setup_numeric_callback(baseline_start, bl_start_bc, "Baseline start", v -> gui_state.baseline_start[] = v; range_check = bl_check)
    setup_numeric_callback(baseline_end, bl_end_bc, "Baseline end", v -> gui_state.baseline_end[] = v; range_check = bl_check)

    # Open file picker when Browse button is clicked
    on(file_select_button.clicks) do _
        default_path = gui_state.directory[] != "" ? gui_state.directory[] : pwd()
        filename = pick_file(default_path)
        if !isnothing(filename) && filename != ""
            # Clear first to reset Makie's internal scroll offset, then set value
            file_pattern_input.displayed_string[] = ""
            file_pattern_input.displayed_string[] = basename(filename)
            gui_state.filename[] = filename
        end
    end

    # Textbox drives gui_state.filename[] for direct typing.
    # Guard: if the textbox value is just the basename of the already-stored full path
    # (set there by Browse), don't replace the full path with the basename.
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
                # Direct plot type (Data Browser, Epochs, ERP, etc.)
                gui_state.plottype[] = selection
                gui_state.submenu_active[] = false
                component_bc[] = :lightgrey
            end
        finally
            updating_submenu[] = false
        end
    end

    # Submenu dropdown callback — missing until now
    on(submenu_dropdown.selection) do selection
        if updating_submenu[] || isnothing(selection) || selection == "Select"
            return
        end
        gui_state.plottype[] = selection
    end


    # Layout dropdown callback
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

    # Open layout file picker when Select button is clicked
    on(layout_select_button.clicks) do _
        default_path = gui_state.directory[] != "" ? gui_state.directory[] : pwd()
        filename = pick_file(default_path)
        if !isnothing(filename) && filename != ""
            basename_only = basename(filename)
            gui_state.layout_file[] = filename
            gui_state.layout[] = basename_only
            layout_label_text[] = basename_only

            # Load layout and populate electrode menu with channels from layout
            try
                layout = read_layout(filename)
                gui_state.layout_object[] = layout
                channel_menu.options = vcat(["Select"], string.(channel_labels(layout)))
                @info "Loaded layout with $(length(channel_labels(layout))) channels"
            catch layout_error
                @minimal_warning "Error loading layout: $layout_error"
                gui_state.layout_object[] = nothing
            end
        end
    end

    # Helper function to update selected channels display
    function update_selected_channels_display()
        if isempty(gui_state.electrodes[])
            selected_channels_text[] = "Selected: "
        else
            selected_channels_text[] = "Selected: " * join(gui_state.electrodes[], ", ")
        end
    end

    on(channel_menu.selection) do selection
        if selection == "Select"
            gui_state.electrodes[] = String[]
        elseif selection in gui_state.electrodes[]
            # Already selected — remove it (deselect)
            gui_state.electrodes[] = filter(!=(selection), gui_state.electrodes[])
        else
            # Not yet selected — add it
            gui_state.electrodes[] = vcat(gui_state.electrodes[], selection)
        end
        update_selected_channels_display()
        # Reset to "no selection" (index 0) so every click — including repeated picks
        # and the "Select" clear option — always registers as a new change event.
        channel_menu.i_selected.val = 0
        channel_menu.selection.val = nothing
    end

    # Multi-number validation helper (accepts space/comma-separated integers)
    function setup_multi_int_callback(input, bc_obs::Observable, field_name::String, state_obs::Observable)
        last_validated = Ref{String}("__uninit__")

        function validate(value)
            stripped = strip(value)
            if stripped == last_validated[]
                return
            end
            last_validated[] = stripped

            if isempty(stripped)
                bc_obs[] = :white
                state_obs[] = ""
            else
                # Split on spaces/commas; each token may be an integer or a n:m range
                tokens = filter(!isempty, split(stripped, r"[\s,]+"))
                valid = all(tokens) do t
                    if occursin(':', t)
                        parts = split(t, ':')
                        length(parts) == 2 && !isnothing(tryparse(Int, parts[1])) && !isnothing(tryparse(Int, parts[2]))
                    else
                        !isnothing(tryparse(Int, t))
                    end
                end
                if !valid
                    bc_obs[] = :lightcoral
                    state_obs[] = "__INVALID__"
                    @minimal_warning "Invalid $field_name value: \"$stripped\" (must be integers or ranges, e.g. 1 2 3 or 1:50)"
                else
                    bc_obs[] = :white
                    state_obs[] = stripped
                end
            end
        end

        on(input.stored_string) do value
            validate(value)
        end

        on(input.focused) do is_focused
            if !is_focused
                validate(input.displayed_string[])
            end
        end
    end

    setup_multi_int_callback(participant_input, participant_bc, "Participant", gui_state.participant)
    setup_multi_int_callback(condition_input, condition_bc, "Condition", gui_state.condition)
    setup_multi_int_callback(epoch_input, epoch_bc, "Epoch", gui_state.epoch)

    # Component input callback (used by ICA bridges via gui_state.ica_components)
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


    # Plot option callbacks

    on(baseline_type.selection) do selection
        if !isnothing(selection) && selection != "Select"
            gui_state.baseline_type[] = selection
        end
    end

    on(invert_y_checkbox.checked) do is_checked
        gui_state.invert_y[] = is_checked
    end

    # Open directory picker when Select button is clicked
    on(directory_select_button.clicks) do _
        default_path = gui_state.directory[] != "" ? gui_state.directory[] : pwd()
        dir_path = pick_folder(default_path)
        if !isnothing(dir_path) && dir_path != ""
            gui_state.directory[] = dir_path
            directory_label_text[] = "." * _truncate_path(String(strip(dir_path)))
        end
    end

    on(plot_button.clicks) do _
        plot()
    end

    # Display the figure
    display(gui_fig)
    _set_window_title("Makie")

    return nothing
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
    is_pattern = lowercase(splitext(fname)[2]) != ".jld2"

    if is_pattern
        # Pattern mode — parse participant/condition filters and forward to plot_epochs
        part_list = _parse_gui_int_field(gui_state.participant[])
        cond_list = _parse_gui_int_field(gui_state.condition[])
        if isnothing(part_list) || isnothing(cond_list)
            @minimal_warning "Requested plot settings incompatible: recheck!"
            return
        end
        part_sel = isempty(part_list) ? participants() : participants(part_list)
        cond_sel = isempty(cond_list) ? conditions() : conditions(cond_list)

        selected_channels = isempty(gui_state.electrodes[]) ? channels() : channels(Symbol.(gui_state.electrodes[]))

        # Warn about any explicitly requested participant IDs that have no matching file
        if !isempty(part_list)
            all_files   = _find_batch_files(fname, gui_state.directory[] == "" ? pwd() : gui_state.directory[])
            avail_ids   = sort(unique(_extract_participant_id.(all_files)))
            missing_ids = filter(id -> id ∉ avail_ids, part_list)
            if !isempty(missing_ids)
                @minimal_warning "Participant ID(s) $missing_ids not found in '$fname'. Available IDs: $avail_ids"
            end
        end

        @async begin
            layout_sym = Symbol(gui_state.layout_type[])
            try
                plot_epochs(
                    fname;
                    input_dir             = gui_state.directory[] == "" ? pwd() : gui_state.directory[],
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
    input_dir = gui_state.directory[] == "" ? pwd() : gui_state.directory[]
    is_pattern = lowercase(splitext(fname)[2]) != ".jld2"

    if is_pattern
        # Pattern mode — parse all relevant GUI filters and forward to plot_erp
        part_list = _parse_gui_int_field(gui_state.participant[])
        cond_list = _parse_gui_int_field(gui_state.condition[])
        if isnothing(part_list) || isnothing(cond_list)
            @minimal_warning "Requested plot settings incompatible: recheck!"
            return
        end
        part_sel = isempty(part_list) ? participants() : participants(part_list)
        cond_sel = isempty(cond_list) ? conditions() : conditions(cond_list)

        selected_channels = isempty(gui_state.electrodes[]) ? channels() : channels(Symbol.(gui_state.electrodes[]))

        # Warn about any explicitly requested participant IDs that have no matching file
        if !isempty(part_list)
            all_files   = _find_batch_files(fname, input_dir)
            avail_ids   = sort(unique(_extract_participant_id.(all_files)))
            missing_ids = filter(id -> id ∉ avail_ids, part_list)
            if !isempty(missing_ids)
                @minimal_warning "Participant ID(s) $missing_ids not found in '$fname'. Available IDs: $avail_ids"
            end
        end

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
    input_dir = gui_state.directory[] == "" ? pwd() : gui_state.directory[]
    is_pattern = lowercase(splitext(fname)[2]) != ".jld2"

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
            catch
                1
            end
            plot_topography(data, epoch_num; channel_selection = selected_channels, interval_selection = interval_sel, extra_kwargs...)
        else
            @minimal_warning "Requested plot settings incompatible: recheck!"
        end
    end

    if is_pattern
        part_list = _parse_gui_int_field(gui_state.participant[])
        cond_list = _parse_gui_int_field(gui_state.condition[])
        if isnothing(part_list) || isnothing(cond_list)
            @minimal_warning "Requested plot settings incompatible: recheck!"
            return
        end
        part_sel = isempty(part_list) ? participants() : participants(part_list)

        # Warn about missing participant IDs
        if !isempty(part_list)
            all_files   = _find_batch_files(fname, input_dir)
            avail_ids   = sort(unique(_extract_participant_id.(all_files)))
            missing_ids = filter(id -> id ∉ avail_ids, part_list)
            if !isempty(missing_ids)
                @minimal_warning "Participant ID(s) $missing_ids not found in '$fname'. Available IDs: $avail_ids"
            end
        end

        files = _find_batch_files(fname, input_dir, part_sel)
        isempty(files) && return  # already warned above if applicable

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
