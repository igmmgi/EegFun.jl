# Base UI styling values
const BASE_FONTS = (title = 24, label = 20, tick = 18, slider = 16, button = 18)
const BASE_SIZES = (input_width = 150, input_height = 25)

# UI styling parameters struct with defaults
struct UIStyle
    title_font::Int
    label_font::Int
    tick_font::Int
    slider_font::Int
    button_font::Int
    input_width::Int
    input_height::Int

    # Constructor with default values
    UIStyle() = new(
        BASE_FONTS.title,
        BASE_FONTS.label,
        BASE_FONTS.tick,
        BASE_FONTS.slider,
        BASE_FONTS.button,
        BASE_SIZES.input_width,
        BASE_SIZES.input_height,
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

    # Create the main figure with minimal padding
    gui_fig = Figure(size = (600, 800), title = "Plot GUI", backgroundcolor = :lightgrey, figure_padding = (0, 0, 0, 0))

    ui_style = UIStyle()

    # Create main layout with 2-column structure
    # We'll add spacing by using column sizing with Relative widths
    main_layout = GridLayout(gui_fig[1, 1:2])

    # Column 1: File & Data Selection
    # Select Directory Section - Place this first to create column 1
    directory_select_button = Button(
        main_layout[1, 1],
        label = "Select Directory",
        fontsize = ui_style.button_font,
        font = :bold,
        width = ui_style.input_width,
        height = ui_style.input_height,
        buttoncolor = :darkgrey,
        buttoncolor_hover = :grey,
        buttoncolor_active = :green,
    )
    directory_label_text = Observable("")
    Label(
        main_layout[2, 1],
        directory_label_text,
        fontsize = ui_style.slider_font,
        width = ui_style.input_width,
        height = ui_style.input_height,
        halign = :center,
    )

    # Select File Section
    file_select_button = Button(
        main_layout[3, 1],
        label = "Select File",
        fontsize = ui_style.button_font,
        font = :bold,
        width = ui_style.input_width,
        height = ui_style.input_height,
        buttoncolor = :darkgrey,
        buttoncolor_hover = :grey,
        buttoncolor_active = :green,
    )
    file_label_text = Observable("")
    Label(
        main_layout[4, 1],
        file_label_text,
        fontsize = ui_style.slider_font,
        width = ui_style.input_width,
        height = ui_style.input_height,
        halign = :center,
    )

    # Column 2: Electrode Section - Place this first to create column 2
    Label(main_layout[1, 2], "Channel(s)", fontsize = ui_style.label_font, font = :bold)
    channel_menu =
        Menu(main_layout[2, 2], options = ["Select"], width = ui_style.input_width, height = ui_style.input_height)

    # Selected channels display
    selected_channels_text = Observable("Selected: ")
    selected_channels_label = Label(main_layout[3, 2], selected_channels_text, fontsize = ui_style.slider_font, color = :gray)

    # Settings Section
    Label(main_layout[4, 2], "Axis Settings", fontsize = ui_style.label_font, font = :bold)

    # Now that both columns exist, set their widths to create spacing
    colsize!(main_layout, 1, Relative(0.45))  # First column takes 45% of width
    colsize!(main_layout, 2, Relative(0.45))  # Second column takes 45% of width
    # The remaining 10% creates natural spacing between columns

    # Layout Section
    layout_select_button = Button(
        main_layout[5, 1],
        label = "Select Layout",
        fontsize = ui_style.button_font,
        font = :bold,
        width = ui_style.input_width,
        height = ui_style.input_height,
        buttoncolor = :darkgrey,
        buttoncolor_hover = :grey,
        buttoncolor_active = :green,
    )
    layout_label_text = Observable("")
    Label(
        main_layout[6, 1],
        layout_label_text,
        fontsize = ui_style.slider_font,
        width = ui_style.input_width,
        height = ui_style.input_height,
        halign = :center,
    )

    # File Filter Section
    Label(main_layout[7, 1], "File Filter", fontsize = ui_style.label_font, font = :bold)
    Textbox(
        main_layout[8, 1],
        placeholder = "",
        fontsize = ui_style.slider_font,
        width = ui_style.input_width,
        height = ui_style.input_height,
        halign = :center,
        boxcolor = :white,
    )

    # Plot Configuration (moved from column 2)
    # Plot Type Section
    Label(main_layout[9, 1], "Plot Type", fontsize = ui_style.label_font, font = :bold)
    plottype_options = [
        "Select",
        "Data Browser",
    ]
    plottype_dropdown = Menu(
        main_layout[10, 1],
        options = plottype_options,
        width = ui_style.input_width,
        height = ui_style.input_height,
    )

    # Participant Section
    Label(main_layout[11, 1], "Participant", fontsize = ui_style.label_font, font = :bold)
    participant_input = Textbox(
        main_layout[12, 1],
        placeholder = "",
        fontsize = ui_style.slider_font,
        width = ui_style.input_width,
        height = ui_style.input_height,
        halign = :center,
        boxcolor = :white,
    )

    # Condition Section
    Label(main_layout[13, 1], "Condition", fontsize = ui_style.label_font, font = :bold)
    condition_input = Textbox(
        main_layout[14, 1],
        placeholder = "",
        fontsize = ui_style.slider_font,
        width = ui_style.input_width,
        height = ui_style.input_height,
        halign = :center,
        boxcolor = :white,
    )


    # Column 2: Settings (Settings title already placed above)
    # X Limits Section
    Label(main_layout[5, 2], "X Limits", fontsize = ui_style.slider_font, font = :bold)
    x_limits_layout = GridLayout(main_layout[6, 2], tellwidth = false, colgap = 10)
    xmin_input = Textbox(
        x_limits_layout[1, 1],
        placeholder = "",
        fontsize = ui_style.slider_font,
        width = 90,
        height = ui_style.input_height,
        halign = :center,
        boxcolor = :white,
    )
    xmax_input = Textbox(
        x_limits_layout[1, 2],
        placeholder = "",
        fontsize = ui_style.slider_font,
        width = 90,
        height = ui_style.input_height,
        halign = :center,
        boxcolor = :white,
    )

    # Y Limits Section
    Label(main_layout[7, 2], "Y Limits", fontsize = ui_style.slider_font, font = :bold)
    y_limits_layout = GridLayout(main_layout[8, 2], tellwidth = false, colgap = 10)
    ymin_input = Textbox(
        y_limits_layout[1, 1],
        placeholder = "",
        fontsize = ui_style.slider_font,
        width = 90,
        height = ui_style.input_height,
        halign = :center,
        boxcolor = :white,
    )
    ymax_input = Textbox(
        y_limits_layout[1, 2],
        placeholder = "",
        fontsize = ui_style.slider_font,
        width = 90,
        height = ui_style.input_height,
        halign = :center,
        boxcolor = :white,
    )

    # Z Limits Section
    Label(main_layout[9, 2], "Z Limits", fontsize = ui_style.slider_font, font = :bold)
    z_limits_layout = GridLayout(main_layout[10, 2], tellwidth = false, colgap = 10)
    Textbox(
        z_limits_layout[1, 1],
        placeholder = "",
        fontsize = ui_style.slider_font,
        width = 90,
        height = ui_style.input_height,
        halign = :center,
        boxcolor = :white,
    )
    Textbox(
        z_limits_layout[1, 2],
        placeholder = "",
        fontsize = ui_style.slider_font,
        width = 90,
        height = ui_style.input_height,
        halign = :center,
        boxcolor = :white,
    )

    # Baseline Section
    Label(main_layout[11, 2], "Baseline", fontsize = ui_style.slider_font, font = :bold)
    baseline_layout = GridLayout(main_layout[12, 2], tellwidth = false, colgap = 10)
    baseline_start = Textbox(
        baseline_layout[1, 1],
        placeholder = "",
        fontsize = ui_style.slider_font,
        width = 90,
        height = ui_style.input_height,
        halign = :center,
        boxcolor = :white,
    )
    baseline_end = Textbox(
        baseline_layout[1, 2],
        placeholder = "",
        fontsize = ui_style.slider_font,
        width = 90,
        height = ui_style.input_height,
        halign = :center,
        boxcolor = :white,
    )

    # Baseline Type
    Label(main_layout[13, 2], "Baseline Type TF", fontsize = ui_style.slider_font, font = :bold)
    baseline_type = Menu(
        main_layout[14, 2],
        options = ["Select", "absolute", "relative", "relchange", "perchange", "db"],
        width = ui_style.input_width,
        height = ui_style.input_height,
    )

    # Action Button Section - Outside the settings panel
    # Main plot button
    plot_button = Button(
        main_layout[15, 2],
        label = "PLOT",
        fontsize = ui_style.title_font,
        buttoncolor = :grey,
        width = 100,
        height = 50,
    )

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


    function execute_plot()

        # Check if we have the required files
        if gui_state.filename[] == ""
            println("Error: No file filter specified!")
            return
        end

        if gui_state.layout_file[] == ""
            println("Error: No layout file selected!")
            return
        end

        try
            # Resolve layout file path - try as-is first, then search in package layouts
            layout_file_path = gui_state.layout_file[]
            if !isfile(layout_file_path)
                # Try to find it in the package's data/layouts directory
                package_layouts_dir = joinpath(@__DIR__, "..", "..", "data", "layouts")
                if isdir(package_layouts_dir)
                    found_layout = eegfun.find_file(basename(layout_file_path), package_layouts_dir)
                    if found_layout !== nothing
                        layout_file_path = found_layout
                    end
                end
            end
            
            if !isfile(layout_file_path)
                println("Error: Layout file not found: $(gui_state.layout_file[])")
                return
            end

            println("Loading data file...")
            # Load data file (could be BDF or other format)
            file_path = gui_state.filename[]
            if endswith(lowercase(file_path), ".bdf")
                dat = eegfun.read_bdf(file_path)
            else
                println("Error: Unsupported file format. Currently only .bdf files are supported.")
                return
            end

            println("Loading layout from: $layout_file_path")
            layout = eegfun.read_layout(layout_file_path)

            println("Creating EEG dataframe...")
            dat = eegfun.create_eeg_dataframe(dat, layout)

            # Update electrode menu with actual channel labels from the loaded data
            println("Updating electrode menu with channel labels...")
            channel_labels = eegfun.channel_labels(dat)
            electrode_options = vcat(["Select"], string.(channel_labels))
            channel_menu.options = electrode_options
            println("Found $(length(channel_labels)) channels: $(join(string.(channel_labels[1:min(10, length(channel_labels))]), ", "))$(length(channel_labels) > 10 ? "..." : "")")

            # Create a new screen/window for the plot to avoid overwriting the GUI
            println("Creating plot...")
            @async begin
                try
                    # Create a new screen and pass it to plot_databrowser
                    new_screen = GLMakie.Screen()
                    eegfun.plot_databrowser(dat; screen = new_screen)
                    println("Plot created and displayed successfully!")
                catch plot_error
                    println("Error creating plot: $plot_error")
                end
            end

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
    # Open file picker when Select File button is clicked
    on(file_select_button.clicks) do _
        @async begin
            try
                # Use the selected directory if available, otherwise use empty string
                default_path = gui_state.directory[] !== nothing && gui_state.directory[] != "" ? gui_state.directory[] : ""
                filename = fetch(Threads.@spawn pick_file(default_path))
                if filename !== nothing && filename != ""
                    basename_only = basename(filename)
                    gui_state.filename[] = filename
                    file_label_text[] = strip(basename_only)
                end
            catch e
                println("File picker error: $e")
            end
        end
    end

    on(plottype_dropdown.selection) do selection
        gui_state.plottype[] = selection
    end

    # Open layout file picker when Select button is clicked
    on(layout_select_button.clicks) do _
        @async begin
            try
                # Use the selected directory if available, otherwise use empty string
                default_path = gui_state.directory[] !== nothing && gui_state.directory[] != "" ? gui_state.directory[] : ""
                filename = fetch(Threads.@spawn pick_file(default_path))
                if filename !== nothing && filename != ""
                    basename_only = basename(filename)
                    gui_state.layout_file[] = filename
                    gui_state.layout[] = basename_only
                    layout_label_text[] = basename_only
                    
                    # Load layout and populate electrode menu with channels from layout
                    try
                        layout = eegfun.read_layout(filename)
                        channel_labels = eegfun.channel_labels(layout)
                        electrode_options = vcat(["Select"], string.(channel_labels))
                        channel_menu.options = electrode_options
                        println("Loaded layout with $(length(channel_labels)) channels")
                    catch layout_error
                        println("Error loading layout: $layout_error")
                    end
                end
            catch e
                println("Layout file picker error: $e")
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
        @async begin
            try
                dir_path = fetch(Threads.@spawn pick_folder(""))
                if dir_path !== nothing && dir_path != ""
                    gui_state.directory[] = dir_path
                    directory_label_text[] = strip(dir_path)
                end
            catch e
                println("Directory picker error: $e")
            end
        end
    end

    # Set up limit input callbacks using helper function
    setup_limit_callback(xmin_input, gui_state.xlim, :min)
    setup_limit_callback(xmax_input, gui_state.xlim, :max)
    setup_limit_callback(ymin_input, gui_state.ylim, :min)
    setup_limit_callback(ymax_input, gui_state.ylim, :max)

    on(plot_button.clicks) do _
        execute_plot()
    end


    # Display the figure
    display(gui_fig)
    return nothing
end
