# Base UI styling values
const BASE_FONTS = (title = 24, label = 20, tick = 18, slider = 16, button = 18)
const BASE_SIZES = (input_width = 150, input_height = 25)
const MIN_FONTS = (title = 16, label = 14, tick = 12, slider = 12, button = 14)
const MIN_SIZES = (input_width = 150, input_height = 25)
const SCALE_FACTORS = (input_width = 200, input_height = 30)

# UI styling parameters struct with defaults
struct UIStyle
    title_font::Observable{Int}
    label_font::Observable{Int}
    tick_font::Observable{Int}
    slider_font::Observable{Int}
    button_font::Observable{Int}
    input_width::Observable{Int}
    input_height::Observable{Int}

    # Constructor with default values
    UIStyle() = new(
        Observable(BASE_FONTS.title),
        Observable(BASE_FONTS.label),
        Observable(BASE_FONTS.tick),
        Observable(BASE_FONTS.slider),
        Observable(BASE_FONTS.button),
        Observable(BASE_SIZES.input_width),
        Observable(BASE_SIZES.input_height),
    )
end

"""
    plot_gui()

Interactive GUI for quick data plotting 

This GUI allows quick data exploration with various plotting options for EEG data.

# Example
```julia
using eegfun
eegfun.plot_my_data_gui()
```

# Returns
- `fig::Figure`: The Makie figure object containing the interactive GUI
"""
function plot_gui()

    # Create the main figure
    fig = Figure(size = (700, 800), title = "Plot GUI", backgroundcolor = :white)

    # Try and make sizing adaptive?
    function setup_adaptive_sizing(fig)
        ui_style = UIStyle()

        # Update fonts and UI elements when figure is resized
        on(fig.scene.viewport) do area
            scale_factor = area.widths[1] / 2000
            ui_style.title_font[] = max(MIN_FONTS.title, round(Int, BASE_FONTS.title * scale_factor))
            ui_style.label_font[] = max(MIN_FONTS.label, round(Int, BASE_FONTS.label * scale_factor))
            ui_style.tick_font[] = max(MIN_FONTS.tick, round(Int, BASE_FONTS.tick * scale_factor))
            ui_style.slider_font[] = max(MIN_FONTS.slider, round(Int, BASE_FONTS.slider * scale_factor))
            ui_style.button_font[] = max(MIN_FONTS.button, round(Int, BASE_FONTS.button * scale_factor))
            ui_style.input_width[] = max(MIN_SIZES.input_width, round(Int, SCALE_FACTORS.input_width * scale_factor))
            ui_style.input_height[] = max(MIN_SIZES.input_height, round(Int, SCALE_FACTORS.input_height * scale_factor))
        end

        return ui_style
    end

    ui_style = setup_adaptive_sizing(fig)

    # Create main layout with proper 3-column structure
    main_layout = GridLayout(fig[1, 1:3], colgap = 10, rowgap = 10)

    # Column 1: File & Data Selection
    # File Type Section
    Label(main_layout[1, 1], "File Type", fontsize = ui_style.label_font, font = :bold)
    filetype_options = [
        "Select",
        "*.bdf",
        "*all.mat",
        "*pre.mat",
        "*ica.mat",
        "*avg.mat",
        "GA_avg.mat",
        "*lrp.mat",
        "GA_lrp.mat",
        "*tf.mat",
        "GA_tf.mat",
        "ERP Stats",
        "TF Stats",
    ]
    filetype_dropdown = Menu(
        main_layout[2, 1],
        options = filetype_options,
        width = ui_style.input_width,
        height = ui_style.input_height,
    )

    # BDF File Name Section
    Label(main_layout[3, 1], "Raw File Name", fontsize = ui_style.label_font, font = :bold)
    bdf_filename_input = Textbox(
        main_layout[4, 1],
        placeholder = "Raw filename ...",
        fontsize = ui_style.slider_font,
        width = ui_style.input_width,
        height = ui_style.input_height,
    )

    # Participant Section
    Label(main_layout[5, 1], "Participant", fontsize = ui_style.label_font, font = :bold)
    participant_input = Textbox(
        main_layout[6, 1],
        placeholder = "Participant numbers ...",
        fontsize = ui_style.slider_font,
        width = ui_style.input_width,
        height = ui_style.input_height,
    )

    # Condition Section
    Label(main_layout[7, 1], "Condition", fontsize = ui_style.label_font, font = :bold)
    condition_input = Textbox(
        main_layout[8, 1],
        placeholder = "Condition numbers ...",
        fontsize = ui_style.slider_font,
        width = ui_style.input_width,
        height = ui_style.input_height,
    )

    # Column 2: Plot Configuration
    # Plot Type Section
    Label(main_layout[1, 2], "Plot Type", fontsize = ui_style.label_font, font = :bold)
    plottype_options = [
        "Select",
        "Data Browser",
        "Data Browser (Python)",
        "Data Summary",
        "Electrode(s)",
        "ERP Image",
        "Multiplot (Grid)",
        "Multiplot (Topo)",
        "Topoplot",
        "Boxplot",
    ]
    plottype_dropdown = Menu(
        main_layout[2, 2],
        options = plottype_options,
        width = ui_style.input_width,
        height = ui_style.input_height,
    )

    # Layout Section
    Label(main_layout[3, 2], "Layout", fontsize = ui_style.label_font, font = :bold)
    layout_options = ["Select", "BioSemi72", "BioSemi70", "BioSemi68", "BioSemi66", "BioSemi64", "Custom"]
    layout_dropdown =
        Menu(main_layout[4, 2], options = layout_options, width = ui_style.input_width, height = ui_style.input_height)

    # Electrode Selection Section
    Label(main_layout[5, 2], "Electrode", fontsize = ui_style.label_font, font = :bold)
    electrode_menu =
        Menu(main_layout[6, 2], options = ["Select"], width = ui_style.input_width, height = ui_style.input_height)

    # Multi-select electrode info
    electrode_info =
        Label(main_layout[7, 2], "Ctrl+Click for multiple selection", fontsize = ui_style.slider_font, color = :gray)

    # Column 3: Settings
    # Settings Title
    Label(main_layout[1, 3], "Settings", fontsize = ui_style.label_font, font = :bold)

    # X Limits Section
    Label(main_layout[2, 3], "X Limits", fontsize = ui_style.slider_font, font = :bold)
    x_limits_layout = GridLayout(main_layout[3, 3], tellwidth = false, colgap = 10)
    xmin_input = Textbox(
        x_limits_layout[1, 1],
        placeholder = "min",
        fontsize = ui_style.slider_font,
        width = 90,
        height = ui_style.input_height,
    )
    xmax_input = Textbox(
        x_limits_layout[1, 2],
        placeholder = "max",
        fontsize = ui_style.slider_font,
        width = 90,
        height = ui_style.input_height,
    )

    # X Topo Series
    Label(main_layout[4, 3], "X Topo Series", fontsize = ui_style.slider_font, font = :bold)
    xtopo_input = Textbox(
        main_layout[5, 3],
        placeholder = "series",
        fontsize = ui_style.slider_font,
        width = ui_style.input_width,
        height = ui_style.input_height,
    )

    # Y Limits Section
    Label(main_layout[6, 3], "Y Limits", fontsize = ui_style.slider_font, font = :bold)
    y_limits_layout = GridLayout(main_layout[7, 3], tellwidth = false, colgap = 10)
    ymin_input = Textbox(
        y_limits_layout[1, 1],
        placeholder = "min",
        fontsize = ui_style.slider_font,
        width = 90,
        height = ui_style.input_height,
    )
    ymax_input = Textbox(
        y_limits_layout[1, 2],
        placeholder = "max",
        fontsize = ui_style.slider_font,
        width = 90,
        height = ui_style.input_height,
    )

    # Z Limits Section
    Label(main_layout[8, 3], "Z Limits", fontsize = ui_style.slider_font, font = :bold)
    z_limits_layout = GridLayout(main_layout[9, 3], tellwidth = false, colgap = 10)
    zmin_input = Textbox(
        z_limits_layout[1, 1],
        placeholder = "min",
        fontsize = ui_style.slider_font,
        width = 90,
        height = ui_style.input_height,
    )
    zmax_input = Textbox(
        z_limits_layout[1, 2],
        placeholder = "max",
        fontsize = ui_style.slider_font,
        width = 90,
        height = ui_style.input_height,
    )

    # Baseline Section
    Label(main_layout[10, 3], "Baseline", fontsize = ui_style.slider_font, font = :bold)
    baseline_layout = GridLayout(main_layout[11, 3], tellwidth = false, colgap = 10)
    baseline_start = Textbox(
        baseline_layout[1, 1],
        placeholder = "start",
        fontsize = ui_style.slider_font,
        width = 90,
        height = ui_style.input_height,
    )
    baseline_end = Textbox(
        baseline_layout[1, 2],
        placeholder = "end",
        fontsize = ui_style.slider_font,
        width = 90,
        height = ui_style.input_height,
    )

    # Baseline Type
    Label(main_layout[12, 3], "Baseline Type", fontsize = ui_style.slider_font, font = :bold)
    baseline_type = Menu(
        main_layout[13, 3],
        options = ["Select", "absolute", "relative", "relchange", "perchange", "db"],
        width = ui_style.input_width,
        height = ui_style.input_height,
    )

    # Action Buttons Section - Outside the settings panel
    # Create a button layout for the third column
    button_layout = GridLayout(main_layout[12:14, 3], tellwidth = false, colgap = 15, rowgap = 10)

    # First row of buttons
    examples_button = Button(
        button_layout[1, 1],
        label = "Examples",
        fontsize = ui_style.button_font,
        width = 120,
        height = ui_style.input_height,
    )

    clear_button = Button(
        button_layout[1, 2],
        label = "Clear GUI",
        fontsize = ui_style.button_font,
        width = 120,
        height = ui_style.input_height,
    )

    # Second row of buttons
    export_button = Button(
        button_layout[2, 1],
        label = "Export",
        fontsize = ui_style.button_font,
        width = 120,
        height = ui_style.input_height,
    )

    cursor_button = Button(
        button_layout[2, 2],
        label = "Cursor",
        fontsize = ui_style.button_font,
        width = 120,
        height = ui_style.input_height,
    )

    # Third row - Main plot button
    plot_button = Button(
        button_layout[3, 1:2],
        label = "PLOT",
        fontsize = ui_style.title_font,
        buttoncolor = :lightblue,
        width = 255,
        height = 50,
    )

    # Data structure to store GUI state
    gui_state = (
        filetype = Observable("select"),
        bdf_filename = Observable(""),
        participant = Observable(""),
        condition = Observable(""),
        additional = Observable(""),
        plottype = Observable("select"),
        layout = Observable("select"),
        layout_file = Observable(""),
        electrodes = Observable(Int[]),
        xlim = Observable((nothing, nothing)),
        ylim = Observable((nothing, nothing)),
        zlim = Observable((nothing, nothing)),
        xtopo = Observable(nothing),
        baseline_start = Observable(nothing),
        baseline_end = Observable(nothing),
        baseline_type = Observable("select"),
    )

    # Callback functions
    function update_filetype(selection)
        # Handle both string and integer selection
        if isa(selection, String)
            selected_option = selection
        else
            selected_option = filetype_options[selection]
        end

        gui_state.filetype[] = selected_option
        # Update electrode labels based on file type
        update_electrode_labels(selection)

        # If BDF is selected, open file picker
        if selected_option == "*.bdf"
            @async begin
                try
                    filename = fetch(Threads.@spawn pick_file(""))
                    if filename !== nothing && filename != ""
                        # Extract just the filename (basename) from the full path
                        basename_only = basename(filename)

                        # Update the textbox with just the filename
                        bdf_filename_input.stored_string[] = basename_only
                        gui_state.bdf_filename[] = filename  # Keep full path in state for processing

                        # Try different properties that might control display
                        if hasfield(typeof(bdf_filename_input), :value)
                            bdf_filename_input.value[] = basename_only
                        end
                        if hasfield(typeof(bdf_filename_input), :displayed_string)
                            bdf_filename_input.displayed_string[] = basename_only
                        end

                        # Force update the display
                        notify(bdf_filename_input.stored_string)

                        println("Selected file: $filename")  # Debug output
                        println("Displaying basename: $basename_only")
                        println("Textbox stored_string: $(bdf_filename_input.stored_string[])")
                    end
                catch e
                    println("File picker error: $e")
                end
            end
        end
    end

    function update_plottype(selection)
        gui_state.plottype[] = plottype_options[selection]
    end

    function update_layout(selection)
        # Handle both string and integer selections
        if isa(selection, String)
            selected_layout = selection
        else
            selected_layout = layout_options[selection]
        end

        gui_state.layout[] = selected_layout

        # If "Select" is chosen, open file picker for custom layout
        if selected_layout == "Select"
            @async begin
                try
                    filename = fetch(Threads.@spawn pick_file(""))
                    if filename !== nothing && filename != ""
                        # Extract just the filename (basename) from the full path
                        basename_only = basename(filename)

                        # Store the full path in the GUI state
                        gui_state.layout_file[] = filename

                        # Update the dropdown to show the selected file
                        # Add the selected file to the options and select it
                        new_options = vcat(layout_options, basename_only)
                        layout_dropdown.options = new_options
                        layout_dropdown.selection = length(new_options)

                        println("Selected custom layout file: $filename")  # Debug output
                        println("Displaying basename: $basename_only")
                    end
                catch e
                    println("Layout file picker error: $e")
                end
            end
        else
            # For predefined layouts, set the layout file path based on the selection
            # This would map to actual layout files in your data/layouts directory
            if selected_layout == "BioSemi72"
                gui_state.layout_file[] = "data/layouts/biosemi/biosemi72.csv"
            elseif selected_layout == "BioSemi70"
                gui_state.layout_file[] = "data/layouts/biosemi/biosemi70.csv"
            elseif selected_layout == "BioSemi68"
                gui_state.layout_file[] = "data/layouts/biosemi/biosemi68.csv"
            elseif selected_layout == "BioSemi66"
                gui_state.layout_file[] = "data/layouts/biosemi/biosemi66.csv"
            elseif selected_layout == "BioSemi64"
                gui_state.layout_file[] = "data/layouts/biosemi/biosemi64.csv"
            elseif selected_layout == "Custom"
                gui_state.layout_file[] = ""  # Will need to be set via file picker
            end

            println("Selected predefined layout: $selected_layout")
            println("Layout file path: $(gui_state.layout_file[])")
        end
    end

    function update_electrode_labels(filetype_selection)
        # This would be implemented to load electrode labels based on file type
        # For now, just show a placeholder
        electrode_menu.options = ["Select", "Fp1", "Fp2", "F3", "F4", "C3", "C4", "P3", "P4", "O1", "O2"]
    end

    function execute_plot()
        println("Executing EEG plot with:")
        println("  File Type: $(gui_state.filetype[])")
        println("  BDF File: $(gui_state.bdf_filename[])")
        println("  Plot Type: $(gui_state.plottype[])")
        println("  Layout: $(gui_state.layout[])")
        println("  Layout File: $(gui_state.layout_file[])")
        println("  Electrodes: $(gui_state.electrodes[])")
        println("  X Limits: $(gui_state.xlim[])")
        println("  Y Limits: $(gui_state.ylim[])")

        # Check if we have the required files
        if gui_state.bdf_filename[] == ""
            println("Error: No BDF file selected!")
            return
        end

        if gui_state.layout_file[] == ""
            println("Error: No layout file selected!")
            return
        end

        try
            println("Loading BDF data...")
            # Load BDF data
            dat = eegfun.read_bdf(gui_state.bdf_filename[])

            println("Loading layout...")
            # Load layout
            layout = eegfun.read_layout(gui_state.layout_file[])

            println("Creating EEG dataframe...")
            # Create EEG dataframe
            dat = eegfun.create_eeg_dataframe(dat, layout)

            println("Creating plot...")
            # Create the plot based on plot type
            if gui_state.plottype[] == "Data Browser"
                eegfun.plot_databrowser(dat)
            else
                # For other plot types, you can add more cases here
                println("Plot type '$(gui_state.plottype[])' not yet implemented")
                # Fallback to databrowser for now
                eegfun.plot_databrowser(dat)
            end

            println("Plot created and displayed successfully!")

        catch e
            println("Error creating EEG plot: $e")
            println("Stacktrace:")
            for (exc, bt) in Base.catch_stack()
                showerror(stdout, exc, bt)
                println()
            end
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
    on(filetype_dropdown.selection) do selection
        update_filetype(selection)
    end

    on(plottype_dropdown.selection) do selection
        update_plottype(selection)
    end

    on(layout_dropdown.selection) do selection
        update_layout(selection)
    end

    on(electrode_menu.selection) do selection
        gui_state.electrodes[] = [selection]
    end

    on(participant_input.stored_string) do value
        gui_state.participant[] = value
    end

    on(condition_input.stored_string) do value
        gui_state.condition[] = value
    end

    # Note: bdf_filename is set directly in the file picker callback, not here
    # to avoid overriding the full path with just the basename


    # Set up limit input callbacks using helper function
    setup_limit_callback(xmin_input, gui_state.xlim, :min)
    setup_limit_callback(xmax_input, gui_state.xlim, :max)
    setup_limit_callback(ymin_input, gui_state.ylim, :min)
    setup_limit_callback(ymax_input, gui_state.ylim, :max)

    on(plot_button.clicks) do _
        execute_plot()
    end

    on(clear_button.clicks) do _
        # Reset all inputs
        filetype_dropdown.selection = 1
        plottype_dropdown.selection = 1
        layout_dropdown.selection = 1
        electrode_menu.selection = 1
        participant_input.stored_string[] = ""
        condition_input.stored_string[] = ""
        bdf_filename_input.stored_string[] = ""
        xmin_input.stored_string[] = ""
        xmax_input.stored_string[] = ""
        ymin_input.stored_string[] = ""
        ymax_input.stored_string[] = ""
        zmin_input.stored_string[] = ""
        zmax_input.stored_string[] = ""
        xtopo_input.stored_string[] = ""
        baseline_start.stored_string[] = ""
        baseline_end.stored_string[] = ""
        baseline_type.selection = 1
        additional_label.text = ""
    end

    # Display the figure
    display(fig)

    return fig
end
