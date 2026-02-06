# Base UI styling values
const BASE_FONTS = (label = 20, button = 18, textbox = 16)
const BASE_SIZES = (input_width = 160, input_height = 22)

# UI styling parameters struct with defaults
struct UIStyle
    label_font::Int
    button_font::Int
    textbox_font::Int
    input_width::Int
    input_height::Int

    UIStyle() = new(BASE_FONTS.label, BASE_FONTS.button, BASE_FONTS.textbox, BASE_SIZES.input_width, BASE_SIZES.input_height)
end

# Helper functions to reduce code repetition
function create_select_button(parent, label, style::UIStyle)
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

function create_textbox(parent, style::UIStyle; width = nothing, placeholder = "")
    Textbox(
        parent,
        placeholder = placeholder,
        fontsize = style.textbox_font,
        width = isnothing(width) ? style.input_width : width,
        height = style.input_height,
        halign = :center,
        boxcolor = :white,
        cornerradius = 0,
    )
end

function create_menu(parent, style::UIStyle; options = ["Select"], width = nothing, height = nothing)
    Menu(
        parent,
        options = options,
        width = isnothing(width) ? style.input_width : width,
        height = isnothing(height) ? style.input_height : height,
        fontsize = style.textbox_font,
    )
end

function create_label(parent, text, style::UIStyle; fontsize = nothing, color = :black)
    fontsize_val = isnothing(fontsize) ? style.label_font : fontsize
    Label(parent, text, fontsize = fontsize_val, color = color)
end

"""Truncate a path string to a maximum length, showing the end with '...' prefix."""
function truncate_path(path::String, max_length::Int = 30)
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
  - Common plots: Data Browser, Epochs, ERP, ERP Image, Topography, GFP
  - Time-Frequency >: Time-Frequency Analysis, Power Spectrum
  - ICA >: ICA Components  
  - Diagnostic Plots >: Artifact Detection, Triggers, Correlation, Layout View, Filter Response

**Column 2: Data Selection & Layout**
- Participant, Condition, Epoch filters
- Layout type: Single, Single Avg, Grid, Topo
- Channel selection with multi-select support

**Column 3: Axis Settings**
- X, Y, Z axis limits
- Baseline correction settings
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
    set_window_title("PLOT GUI")
    gui_fig = Figure(size = (650, 550), title = "Plot GUI", backgroundcolor = :lightgrey, figure_padding = (20, 20, 20, 20))
    main_layout = GridLayout(gui_fig[1, 1:3], rowgap = 2, colgap = 4)
    ui_style = UIStyle()

    # Select Directory Section
    directory_select_button = create_select_button(main_layout[1, 1], "Select Directory", ui_style)
    directory_label_text = Observable("Dir: ")
    create_label(main_layout[2, 1], directory_label_text, ui_style, fontsize = ui_style.textbox_font, color = :gray)


    # Select File Section
    file_select_button = create_select_button(main_layout[3, 1], "Select File", ui_style)
    file_label_text = Observable("File:")
    create_label(main_layout[4, 1], file_label_text, ui_style, fontsize = ui_style.textbox_font, color = :gray)

    # Layout Section
    layout_select_button = create_select_button(main_layout[5, 1], "Select Layout", ui_style)
    layout_label_text = Observable("File: ")
    create_label(main_layout[6, 1], layout_label_text, ui_style, fontsize = ui_style.textbox_font, color = :gray)

    # File Filter Section
    create_label(main_layout[7, 1], "File Filter", ui_style)
    create_textbox(main_layout[8, 1], ui_style)

    # Plot Type Section
    create_label(main_layout[9, 1], "Plot Type", ui_style)
    plottype_options = [
        "Select",
        "Data Browser",
        "Epochs",
        "ERP",
        "ERP Image",
        "Topography",
        "Global Field Power",
        "─────────────────",  # Visual separator
        "Time-Frequency >",
        "ICA >",
        "Diagnostic Plots >",
    ]
    plottype_dropdown = create_menu(main_layout[10, 1], ui_style, options = plottype_options)

    # Submenu dropdown (starts with placeholder to appear inactive)
    submenu_dropdown = create_menu(main_layout[11, 1], ui_style, options = ["---"])

    # Column 2: Participant, Condition, Epoch, Channels
    # Participant Section
    create_label(main_layout[1, 2], "Participant", ui_style)
    participant_input = create_textbox(main_layout[2, 2], ui_style)

    # Condition Section
    create_label(main_layout[3, 2], "Condition", ui_style)
    condition_input = create_textbox(main_layout[4, 2], ui_style)

    # Epoch Section
    create_label(main_layout[5, 2], "Epoch", ui_style)
    epoch_input = create_textbox(main_layout[6, 2], ui_style)

    # Layout Section (replaces Channel(s) and Average Channels)
    create_label(main_layout[7, 2], "Layout", ui_style)
    layout_dropdown = create_menu(main_layout[8, 2], ui_style, options = ["Single", "Single Avg", "Grid", "Topo"])

    # Channel(s) Section (moved down)
    create_label(main_layout[9, 2], "Channel(s)", ui_style)
    channel_menu = create_menu(main_layout[10, 2], ui_style, options = ["Select"])

    # Selected channels display
    selected_channels_text = Observable("Channels: ")
    create_label(main_layout[11, 2], selected_channels_text, ui_style, fontsize = ui_style.textbox_font, color = :gray)

    # Column 3: Axis Settings
    # X Limits Section
    create_label(main_layout[1, 3], "X Limits", ui_style, fontsize = ui_style.textbox_font)
    x_limits_layout = GridLayout(main_layout[2, 3], tellwidth = false, colgap = 8)
    xmin_input = create_textbox(x_limits_layout[1, 1], ui_style, width = 80)
    xmax_input = create_textbox(x_limits_layout[1, 2], ui_style, width = 80)

    # Y Limits Section
    create_label(main_layout[3, 3], "Y Limits", ui_style, fontsize = ui_style.textbox_font)
    y_limits_layout = GridLayout(main_layout[4, 3], tellwidth = false, colgap = 8)
    ymin_input = create_textbox(y_limits_layout[1, 1], ui_style, width = 80)
    ymax_input = create_textbox(y_limits_layout[1, 2], ui_style, width = 80)

    # Z Limits Section (currently unused, but kept for future use)
    create_label(main_layout[5, 3], "Z Limits", ui_style, fontsize = ui_style.textbox_font)
    z_limits_layout = GridLayout(main_layout[6, 3], tellwidth = false, colgap = 8)
    zmin_input = create_textbox(z_limits_layout[1, 1], ui_style, width = 80)
    zmax_input = create_textbox(z_limits_layout[1, 2], ui_style, width = 80)

    # Baseline Section
    create_label(main_layout[7, 3], "Baseline", ui_style, fontsize = ui_style.textbox_font)
    baseline_layout = GridLayout(main_layout[8, 3], tellwidth = false, colgap = 8)
    baseline_start = create_textbox(baseline_layout[1, 1], ui_style, width = 80)
    baseline_end = create_textbox(baseline_layout[1, 2], ui_style, width = 80)

    # Baseline Type
    create_label(main_layout[9, 3], "Baseline Type TF", ui_style, fontsize = ui_style.textbox_font)
    baseline_type = create_menu(main_layout[10, 3], ui_style, options = ["Select", "absolute", "relative", "relchange", "perchange", "db"])

    # Invert Y axis option
    invert_y_layout = GridLayout(main_layout[11, 3], tellwidth = false, colgap = 8)
    create_label(invert_y_layout[1, 1], "Invert Y axis", ui_style, fontsize = ui_style.textbox_font)
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
        xlim = Observable((nothing, nothing)),
        ylim = Observable((nothing, nothing)),
        zlim = Observable((nothing, nothing)),
        baseline_start = Observable(nothing),
        baseline_end = Observable(nothing),
        baseline_type = Observable("select"),
        layout_type = Observable("single"),
        average_channels = Observable(false),
        invert_y = Observable(false),
        submenu_active = Observable(false),
        submenu_type = Observable(""),
        channel_menu = channel_menu,
    )

    function plot()
        plot_type = gui_state.plottype[]
        # Skip placeholder and separator options
        if plot_type == "Select" || startswith(plot_type, "─")
            return
        end
        if plot_type == "Data Browser"
            _plot_databrowser(gui_state)
        elseif startswith(plot_type, "Epochs")
            _plot_epochs(gui_state)
        elseif startswith(plot_type, "ERP Image")
            _plot_erp_image(gui_state)
        elseif startswith(plot_type, "ERP")
            _plot_erp(gui_state)
        elseif plot_type == "Topography"
            _plot_topography(gui_state)
        elseif plot_type == "Global Field Power"
            _plot_gfp(gui_state)
        elseif plot_type == "Time-Frequency"
            _plot_time_frequency(gui_state)
        elseif plot_type == "Power Spectrum"
            _plot_power_spectrum(gui_state)
        elseif plot_type == "ICA Components"
            _plot_ica(gui_state)
        elseif plot_type == "Filter Response"
            _plot_filter(gui_state)
        elseif plot_type == "Artifact Detection"
            _plot_artifacts(gui_state)
        elseif plot_type == "Triggers"
            _plot_triggers(gui_state)
        elseif plot_type == "Correlation Heatmap"
            _plot_correlation(gui_state)
        elseif plot_type == "Layout View"
            _plot_layout(gui_state)
        else
            println("Error: Unsupported plot type: $plot_type")
        end
    end

    # Helper function for limit input callbacks
    function setup_limit_callback(input, limit_observable, position::Symbol)
        on(input.stored_string) do value
            parsed_value = tryparse(Float64, value)
            current = limit_observable[]
            new_tuple = position == :min ? (parsed_value, current[2]) : (current[1], parsed_value)
            limit_observable[] = new_tuple
        end
    end

    # Connect callbacks
    # Open file picker when Select File button is clicked
    on(file_select_button.clicks) do _
        default_path = gui_state.directory[] !== nothing && gui_state.directory[] != "" ? gui_state.directory[] : ""
        filename = fetch(Threads.@spawn pick_file(default_path))
        if filename !== nothing && filename != ""
            basename_only = basename(filename)
            gui_state.filename[] = filename
            file_label_text[] = strip(basename_only)
        end
    end

    on(plottype_dropdown.selection) do selection
        if selection == "Select" || selection == "─────────────────"
            return
        elseif selection == "Time-Frequency >"
            # Show submenu with time-frequency plot options
            submenu_dropdown.options = ["Select", "Time-Frequency", "Power Spectrum"]
            gui_state.submenu_active[] = true
            gui_state.submenu_type[] = "timefreq"
        elseif selection == "Diagnostic Plots >"
        elseif selection == "ICA >"
            # Show submenu with ICA plot options
            submenu_dropdown.options = ["Select", "ICA Components"]
            gui_state.submenu_active[] = true
            gui_state.submenu_type[] = "ica"
            # Show submenu with diagnostic plot options
            submenu_dropdown.options = ["Select", "Artifact Detection", "Triggers", "Correlation Heatmap", "Layout View", "Filter Response"]
            gui_state.submenu_active[] = true
            gui_state.submenu_type[] = "diagnostic"
        elseif selection == "Data Browser" || selection == "Topography" || selection == "Global Field Power"
            # Direct plot types without layout variants
            gui_state.plottype[] = selection
            gui_state.submenu_active[] = false
        elseif selection == "Epochs" || selection == "ERP" || selection == "ERP Image"
            # Plot types that support layout variants
            gui_state.plottype[] = selection
            gui_state.submenu_active[] = false
            # Keep submenu showing previous selection
        else
            # Fallback for any other selections
            gui_state.plottype[] = selection
            gui_state.submenu_active[] = false
            # Keep submenu showing previous selection
        end
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
        default_path = gui_state.directory[] !== nothing && gui_state.directory[] != "" ? gui_state.directory[] : ""
        filename = fetch(Threads.@spawn pick_file(default_path))
        if filename !== nothing && filename != ""
            basename_only = basename(filename)
            gui_state.layout_file[] = filename
            gui_state.layout[] = basename_only
            layout_label_text[] = basename_only

            # Load layout and populate electrode menu with channels from layout
            try
                layout = read_layout(filename)
                gui_state.layout_object[] = layout
                channel_menu.options = vcat(["Select"], string.(channel_labels(layout)))
                println("Loaded layout with $(length(channel_labels(layout))) channels")
            catch layout_error
                println("Error loading layout: $layout_error")
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
        else
            # Add to selection if not already present
            if !(selection in gui_state.electrodes[])
                gui_state.electrodes[] = vcat(gui_state.electrodes[], selection)
            end
        end
        update_selected_channels_display()
    end

    on(participant_input.stored_string) do value
        gui_state.participant[] = value
    end

    on(condition_input.stored_string) do value
        gui_state.condition[] = value
    end

    on(epoch_input.stored_string) do value
        gui_state.epoch[] = value
    end

    # Plot option callbacks

    on(invert_y_checkbox.checked) do is_checked
        gui_state.invert_y[] = is_checked
    end

    # Open directory picker when Select button is clicked
    on(directory_select_button.clicks) do _
        dir_path = fetch(Threads.@spawn pick_folder(""))
        if dir_path !== nothing && dir_path != ""
            gui_state.directory[] = dir_path
            directory_label_text[] = "Dir: " * truncate_path(String(strip(dir_path)))
        end
    end

    # Set up limit input callbacks using helper function
    setup_limit_callback(xmin_input, gui_state.xlim, :min)
    setup_limit_callback(xmax_input, gui_state.xlim, :max)
    setup_limit_callback(ymin_input, gui_state.ylim, :min)
    setup_limit_callback(ymax_input, gui_state.ylim, :max)

    on(plot_button.clicks) do _
        plot()
    end

    # Display the figure
    display(gui_fig)
    set_window_title("Makie")

    return nothing
end



function _plot_databrowser(gui_state)

    # Check if we have the required files
    gui_state.filename[] == "" && @minimal_error "Error: No file specified!"

    # Read data file (could be BDF or other format)
    file_ext = lowercase(splitext(gui_state.filename[])[2])  # [2] is the extension (with dot)
    if file_ext ∉ [".bdf"]
        @minimal_error "Error: Unsupported file format"
    end

    try

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

        # Create a new screen/window for the plot
        @async begin
            new_screen = GLMakie.Screen()
            plot_databrowser(dat; screen = new_screen)
        end
    catch e
        println("Error creating EEG plot: $e")
        showerror(stdout, e, catch_backtrace())
    end
end

function _plot_epochs(gui_state)
    # Check if we have the required files
    gui_state.filename[] == "" && @minimal_error "Error: No file specified!"
    gui_state.layout_file[] == "" && @minimal_error "Error: No layout file selected!"

    # Read data file (should be JLD2 with EpochData)
    file_ext = lowercase(splitext(gui_state.filename[])[2])
    if file_ext != ".jld2"
        @minimal_error "Error: Epochs plot requires JLD2 file format"
    end

    try
        # Read data from JLD2 file
        data = read_data(gui_state.filename[])
        isnothing(data) && @minimal_error "Error: No data found in file"

        # Build channel selection from GUI
        selected_channels = isempty(gui_state.electrodes[]) ? channels() : channels(gui_state.electrodes[])

        # Create a new screen/window for the plot
        @async begin
            layout_sym = Symbol(gui_state.layout_type[])
            plot_epochs(
                gui_state.filename[];
                channel_selection = selected_channels,
                layout = layout_sym,
                average_channels = gui_state.average_channels[],
                xlim = gui_state.xlim[],
                ylim = gui_state.ylim[],
            )
        end
    catch e
        println("Error creating epochs plot: $e")
        showerror(stdout, e, catch_backtrace())
    end
end

function _plot_erp(gui_state)
    # Check if we have the required files
    gui_state.filename[] == "" && @minimal_error "Error: No file specified!"
    gui_state.layout_file[] == "" && @minimal_error "Error: No layout file selected!"

    # Read data file (should be JLD2 with ErpData)
    file_ext = lowercase(splitext(gui_state.filename[])[2])
    if file_ext != ".jld2"
        @minimal_error "Error: ERP plot requires JLD2 file format"
    end

    try
        # Build channel selection from GUI
        selected_channels = isempty(gui_state.electrodes[]) ? channels() : channels(gui_state.electrodes[])

        # Build baseline interval if provided
        baseline = nothing
        if gui_state.baseline_start[] !== nothing && gui_state.baseline_end[] !== nothing
            baseline = (gui_state.baseline_start[], gui_state.baseline_end[])
        end

        # Create a new screen/window for the plot
        @async begin
            layout_sym = Symbol(gui_state.layout_type[])
            plot_erp(
                gui_state.filename[];
                layout = layout_sym,
                channel_selection = selected_channels,
                baseline_interval = baseline,
                average_channels = gui_state.average_channels[],
                xlim = gui_state.xlim[],
                ylim = gui_state.ylim[],
            )
        end
    catch e
        println("Error creating ERP plot: $e")
        showerror(stdout, e, catch_backtrace())
    end
end

function _plot_topography(gui_state)
    # Check if we have the required files
    gui_state.filename[] == "" && @minimal_error "Error: No file specified!"
    gui_state.layout_file[] == "" && @minimal_error "Error: No layout file selected!"

    # Read data file (should be JLD2 with ErpData or EpochData)
    file_ext = lowercase(splitext(gui_state.filename[])[2])
    if file_ext != ".jld2"
        @minimal_error "Error: Topography plot requires JLD2 file format"
    end

    try
        # Read data from JLD2 file
        data = read_data(gui_state.filename[])
        isnothing(data) && @minimal_error "Error: No data found in file"

        # Build channel selection from GUI
        selected_channels = isempty(gui_state.electrodes[]) ? channels() : channels(gui_state.electrodes[])

        # Create a new screen/window for the plot
        @async begin
            # For topography, we need to use the layout from the loaded data or from the layout file
            if data isa Vector{<:ErpData} || data isa ErpData
                plot_topography(data; channel_selection = selected_channels)
            elseif data isa Vector{<:EpochData} || data isa EpochData
                # For EpochData, we need to specify an epoch number
                # Parse epoch from GUI input, default to 1 if empty or invalid
                epoch_num = try
                    epoch_str = strip(gui_state.epoch[])
                    isempty(epoch_str) ? 1 : parse(Int, epoch_str)
                catch
                    1  # Default to epoch 1 if parsing fails
                end
                plot_topography(data, epoch_num; channel_selection = selected_channels)
            else
                @minimal_error "Error: Topography plot requires ErpData or EpochData"
            end
        end
    catch e
        println("Error creating topography plot: $e")
        showerror(stdout, e, catch_backtrace())
    end
end
