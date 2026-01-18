# Base UI styling values
const BASE_FONTS = (label = 20, button = 18, textbox = 14)
const BASE_SIZES = (input_width = 150, input_height = 25)

# UI styling parameters struct with defaults
struct UIStyle
    label_font::Int
    button_font::Int
    textbox_font::Int
    input_width::Int
    input_height::Int

    UIStyle() =
        new(BASE_FONTS.label, BASE_FONTS.button, BASE_FONTS.textbox, BASE_SIZES.input_width, BASE_SIZES.input_height)
end

# Helper functions to reduce code repetition
function create_select_button(parent, label, style::UIStyle)
    Button(
        parent,
        label = label,
        fontsize = style.button_font,
        font = :bold,
        width = style.input_width,
        height = style.input_height,
        buttoncolor = :darkgrey,
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
    )
end

function create_label(parent, text, style::UIStyle; fontsize = nothing, bold = true, color = :black)
    fontsize_val = isnothing(fontsize) ? style.label_font : fontsize
    Label(parent, text, fontsize = fontsize_val, font = bold ? :bold : :regular, color = color)
end

function create_label_display(parent, observable, style::UIStyle)
    Label(
        parent,
        observable,
        fontsize = style.textbox_font,
        width = style.input_width,
        height = style.input_height,
        halign = :center,
    )
end


"""
    plot_gui()

Interactive GUI for quick data plotting 

"""
function plot_gui()

    # main figure window, layout, and UI style
    set_window_title("PLOT GUI")
    gui_fig = Figure(
        size = (600, 800),
        title = "Plot GUI",
        backgroundcolor = :lightgrey,
        figure_padding = (0, 0, 0, 0),  # Remove padding: (left, right, bottom, top)
    )
    main_layout = GridLayout(gui_fig[1, 1:2])
    ui_style = UIStyle()

    #########################################################
    # COLUMN 1: SETTINGS
    #########################################################
    # Select Directory Section - Place this first to create column 1
    directory_select_button = create_select_button(main_layout[1, 1], "Select Directory", ui_style)
    directory_label_text = Observable("")
    create_label_display(main_layout[2, 1], directory_label_text, ui_style)

    # Select File Section
    file_select_button = create_select_button(main_layout[3, 1], "Select File", ui_style)
    file_label_text = Observable("")
    create_label_display(main_layout[4, 1], file_label_text, ui_style)

    # Column 2: Electrode Section - Place this first to create column 2
    create_label(main_layout[1, 2], "Channel(s)", ui_style)
    channel_menu =
        Menu(main_layout[2, 2], options = ["Select"], width = ui_style.input_width, height = ui_style.input_height)

    # Selected channels display
    selected_channels_text = Observable("Selected: ")
    create_label(
        main_layout[3, 2],
        selected_channels_text,
        ui_style;
        fontsize = ui_style.textbox_font,
        bold = false,
        color = :gray,
    )

    # Settings Section
    create_label(main_layout[4, 2], "Axis Settings", ui_style)

    # Layout Section
    layout_select_button = create_select_button(main_layout[5, 1], "Select Layout", ui_style)
    layout_label_text = Observable("")
    create_label_display(main_layout[6, 1], layout_label_text, ui_style)

    # File Filter Section
    create_label(main_layout[7, 1], "File Filter", ui_style)
    create_textbox(main_layout[8, 1], ui_style)

    # Plot Configuration (moved from column 2)
    # Plot Type Section
    create_label(main_layout[9, 1], "Plot Type", ui_style)
    plottype_options = ["Select", "Data Browser"]
    plottype_dropdown = Menu(
        main_layout[10, 1],
        options = plottype_options,
        width = ui_style.input_width,
        height = ui_style.input_height,
    )

    # Participant Section
    create_label(main_layout[11, 1], "Participant", ui_style)
    participant_input = create_textbox(main_layout[12, 1], ui_style)

    # Condition Section
    create_label(main_layout[13, 1], "Condition", ui_style)
    condition_input = create_textbox(main_layout[14, 1], ui_style)

    #########################################################
    # COLUMN 2: SETTINGS
    #########################################################
    # X Limits Section
    create_label(main_layout[5, 2], "X Limits", ui_style; fontsize = ui_style.textbox_font)
    x_limits_layout = GridLayout(main_layout[6, 2], tellwidth = false, colgap = 10)
    xmin_input = create_textbox(x_limits_layout[1, 1], ui_style; width = 90)
    xmax_input = create_textbox(x_limits_layout[1, 2], ui_style; width = 90)

    # Y Limits Section
    create_label(main_layout[7, 2], "Y Limits", ui_style; fontsize = ui_style.textbox_font)
    y_limits_layout = GridLayout(main_layout[8, 2], tellwidth = false, colgap = 10)
    ymin_input = create_textbox(y_limits_layout[1, 1], ui_style; width = 90)
    ymax_input = create_textbox(y_limits_layout[1, 2], ui_style; width = 90)

    # Z Limits Section
    create_label(main_layout[9, 2], "Z Limits", ui_style; fontsize = ui_style.textbox_font)
    z_limits_layout = GridLayout(main_layout[10, 2], tellwidth = false, colgap = 10)
    create_textbox(z_limits_layout[1, 1], ui_style; width = 90)
    create_textbox(z_limits_layout[1, 2], ui_style; width = 90)

    # Baseline Section
    create_label(main_layout[11, 2], "Baseline", ui_style; fontsize = ui_style.textbox_font)
    baseline_layout = GridLayout(main_layout[12, 2], tellwidth = false, colgap = 10)
    baseline_start = create_textbox(baseline_layout[1, 1], ui_style; width = 90)
    baseline_end = create_textbox(baseline_layout[1, 2], ui_style; width = 90)

    # Baseline Type
    create_label(main_layout[13, 2], "Baseline Type TF", ui_style; fontsize = ui_style.textbox_font)
    baseline_type = Menu(
        main_layout[14, 2],
        options = ["Select", "absolute", "relative", "relchange", "perchange", "db"],
        width = ui_style.input_width,
        height = ui_style.input_height,
    )

    # Action Button Section - Outside the settings panel
    # Main plot button
    plot_button = Button(
        main_layout[15, 1:2],
        label = "PLOT",
        fontsize = ui_style.label_font,
        buttoncolor = :darkgrey,
        buttoncolor_hover = :grey,
        buttoncolor_active = :green,
        width = 100,
        height = 50,
    )

    # Set columns sizes
    colsize!(main_layout, 1, Relative(0.4))
    colsize!(main_layout, 2, Relative(0.4))

    # Data structure to store GUI state
    gui_state = (
        directory = Observable(""),
        filetype = Observable("select"),
        filename = Observable(""),
        participant = Observable(""),
        condition = Observable(""),
        additional = Observable(""),
        plottype = Observable("select"),
        layout = Observable("select"),
        layout_file = Observable(""),
        electrodes = Observable(String[]),
        xlim = Observable((nothing, nothing)),
        ylim = Observable((nothing, nothing)),
        zlim = Observable((nothing, nothing)),
        baseline_start = Observable(nothing),
        baseline_end = Observable(nothing),
        baseline_type = Observable("select"),
    )

    function plot()
        if gui_state.plottype[] == "Data Browser"
            _plot_databrowser(gui_state, channel_menu)
        else
            println("Error: Unsupported plot type. Currently only Data Browser is supported.")
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
        default_path =
            gui_state.directory[] !== nothing && gui_state.directory[] != "" ? gui_state.directory[] : ""
        filename = fetch(Threads.@spawn pick_file(default_path))
        if filename !== nothing && filename != ""
            basename_only = basename(filename)
            gui_state.filename[] = filename
            file_label_text[] = strip(basename_only)
        end
    end

    on(plottype_dropdown.selection) do selection
        gui_state.plottype[] = selection
    end

    # Open layout file picker when Select button is clicked
    on(layout_select_button.clicks) do _
        default_path =
            gui_state.directory[] !== nothing && gui_state.directory[] != "" ? gui_state.directory[] : ""
        filename = fetch(Threads.@spawn pick_file(default_path))
        if filename !== nothing && filename != ""
            basename_only = basename(filename)
            gui_state.layout_file[] = filename
            gui_state.layout[] = basename_only
            layout_label_text[] = basename_only

            # Load layout and populate electrode menu with channels from layout
            try
                layout = read_layout(filename)
                gui_state.layout_object[] = layout  # Store for later use
                channel_menu.options = vcat(["Select"], string.(channel_labels(layout)))
                println("Loaded layout with $(length(channel_labels)) channels")
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
            update_selected_channels_display()
        elseif selection != "Select"
            # Add to selection if not already present
            if !(selection in gui_state.electrodes[])
                gui_state.electrodes[] = vcat(gui_state.electrodes[], selection)
            end
            update_selected_channels_display()
        end
    end

    on(participant_input.stored_string) do value
        gui_state.participant[] = value
    end

    on(condition_input.stored_string) do value
        gui_state.condition[] = value
    end

    # Open directory picker when Select button is clicked
    on(directory_select_button.clicks) do _
        dir_path = fetch(Threads.@spawn pick_folder(""))
        if dir_path !== nothing && dir_path != ""
            gui_state.directory[] = dir_path
            directory_label_text[] = strip(dir_path)
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



function _plot_databrowser(gui_state, channel_menu)

    # Check if we have the required files
    gui_state.filename[] == "" && @minimal_error "Error: No file specified!"
    gui_state.layout_file[] == "" && @minimal_error "Error: No layout file selected!"

    # Load data file (could be BDF or other format)
    file_ext = lowercase(splitext(gui_state.filename[])[2])  # [2] is the extension (with dot)
    if file_ext ∉ [".bdf"]
        @minimal_error "Error: Unsupported file format"
    end

    try

        dat = read_bdf(gui_state.filename[])
        dat = create_eeg_dataframe(dat, gui_state.layout_file[])

        # Update electrode menu with actual channel labels from the loaded data
        channel_menu.options = vcat(["Select"], string.(channel_labels(dat)))

        # Create a new screen/window for the plot to avoid overwriting the GUI
        @async begin
            new_screen = GLMakie.Screen()
            plot_databrowser(dat; screen = new_screen)
        end
    catch e
        println("Error creating EEG plot: $e")
        showerror(stdout, e, catch_backtrace())
    end
end