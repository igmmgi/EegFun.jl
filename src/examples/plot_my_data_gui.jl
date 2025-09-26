"""
    plot_my_data_gui()

Interactive GUI for EEG data plotting - Julia/Makie version of plotMyData.m

This GUI allows quick data exploration with various plotting options for EEG data.
It provides a user-friendly interface for selecting file types, plot types, 
electrodes, and various plotting parameters.

# Example
```julia
using eegfun
eegfun.plot_my_data_gui()
```

# Returns
- `fig::Figure`: The Makie figure object containing the interactive GUI
"""
function plot_my_data_gui()
    # Create the main figure
    fig = Figure(
        size = (1200, 800),
        title = "Plot My Data",
        backgroundcolor = :white
    )
    
    # Set up adaptive font and UI sizing
    function setup_adaptive_sizing(fig)
        # Create observables for different font sizes and UI elements
        title_font = Observable(24)
        label_font = Observable(20)
        tick_font = Observable(18)
        slider_font = Observable(16)
        button_font = Observable(18)
        input_width = Observable(200)
        input_height = Observable(30)
        
        # Update fonts and UI elements when figure is resized
        on(fig.scene.viewport) do area
            scale_factor = area.widths[1] / 1200  # Base on 1200px width
            title_font[] = max(16, round(Int, 24 * scale_factor))
            label_font[] = max(14, round(Int, 20 * scale_factor))
            tick_font[] = max(12, round(Int, 18 * scale_factor))
            slider_font[] = max(12, round(Int, 16 * scale_factor))
            button_font[] = max(14, round(Int, 18 * scale_factor))
            input_width[] = max(150, round(Int, 200 * scale_factor))
            input_height[] = max(25, round(Int, 30 * scale_factor))
        end
        
        return title_font, label_font, tick_font, slider_font, button_font, input_width, input_height
    end
    
    title_font, label_font, tick_font, slider_font, button_font, input_width, input_height = setup_adaptive_sizing(fig)
    
    # Create main layout with proper 3-column structure
    main_layout = GridLayout(fig[1, 1:3], 
                            colgap = 30, 
                            rowgap = 20,
                            padding = (30, 30, 30, 30))
    
    # Column 1: File & Data Selection
    # File Type Section
    Label(main_layout[1, 1], "File Type", fontsize = label_font, font = :bold)
    filetype_options = ["Select", "*.bdf", "*all.mat", "*pre.mat", "*ica.mat", "*avg.mat", "GA_avg.mat", "*lrp.mat", "GA_lrp.mat", "*tf.mat", "GA_tf.mat", "ERP Stats", "TF Stats"]
    filetype_dropdown = Menu(main_layout[2, 1], 
                            options = filetype_options,
                            width = input_width,
                            height = input_height)
    
    # BDF File Name Section
    Label(main_layout[3, 1], "BDF File Name", fontsize = label_font, font = :bold)
    bdf_filename_input = Textbox(main_layout[4, 1], 
                                placeholder = "Enter BDF filename...",
                                fontsize = slider_font,
                                width = input_width,
                                height = input_height)
    
    # Participant Section
    Label(main_layout[5, 1], "Participant", fontsize = label_font, font = :bold)
    participant_input = Textbox(main_layout[6, 1], 
                               placeholder = "Enter participant numbers...",
                               fontsize = slider_font,
                               width = input_width,
                               height = input_height)
    
    # Condition Section
    Label(main_layout[7, 1], "Condition", fontsize = label_font, font = :bold)
    condition_input = Textbox(main_layout[8, 1], 
                             placeholder = "Enter condition numbers...",
                             fontsize = slider_font,
                             width = input_width,
                             height = input_height)
    
    # Additional Inputs Button
    additional_button = Button(main_layout[9, 1], 
                              label = "Additional Inputs",
                              fontsize = button_font,
                              width = input_width,
                              height = input_height)
    
    # Additional Inputs Label
    additional_label = Label(main_layout[10, 1], "", 
                            fontsize = slider_font,
                            width = input_width,
                            height = input_height)
    
    # Column 2: Plot Configuration
    # Plot Type Section
    Label(main_layout[1, 2], "Plot Type", fontsize = label_font, font = :bold)
    plottype_options = ["Select", "Data Browser", "Data Browser (Python)", "Data Summary", "Electrode(s)", "ERP Image", "Multiplot (Grid)", "Multiplot (Topo)", "Topoplot", "Boxplot"]
    plottype_dropdown = Menu(main_layout[2, 2], 
                            options = plottype_options,
                            width = input_width,
                            height = input_height)
    
    # Layout Section
    Label(main_layout[3, 2], "Layout", fontsize = label_font, font = :bold)
    layout_options = ["Select", "BioSemi72", "BioSemi70", "BioSemi68", "BioSemi66", "BioSemi64", "Custom"]
    layout_dropdown = Menu(main_layout[4, 2], 
                          options = layout_options,
                          width = input_width,
                          height = input_height)
    
    # Electrode Selection Section
    Label(main_layout[5, 2], "Electrode", fontsize = label_font, font = :bold)
    electrode_menu = Menu(main_layout[6, 2], 
                         options = ["Select"], 
                         width = input_width,
                         height = input_height)
    
    # Multi-select electrode info
    electrode_info = Label(main_layout[7, 2], "Use Ctrl+Click for multiple selection", 
                          fontsize = slider_font,
                          color = :gray)
    
    # Column 3: Settings
    # Settings Title
    Label(main_layout[1, 3], "Settings", fontsize = label_font, font = :bold)
    
    # X Limits Section
    Label(main_layout[2, 3], "X Limits", fontsize = slider_font, font = :bold)
    x_limits_layout = GridLayout(main_layout[3, 3], 
                                tellwidth = false, 
                                colgap = 10)
    xmin_input = Textbox(x_limits_layout[1, 1], 
                        placeholder = "min", 
                        fontsize = slider_font,
                        width = 90,
                        height = input_height)
    xmax_input = Textbox(x_limits_layout[1, 2], 
                        placeholder = "max", 
                        fontsize = slider_font,
                        width = 90,
                        height = input_height)
    
    # X Topo Series
    Label(main_layout[4, 3], "X Topo Series", fontsize = slider_font, font = :bold)
    xtopo_input = Textbox(main_layout[5, 3], 
                         placeholder = "series", 
                         fontsize = slider_font,
                         width = input_width,
                         height = input_height)
    
    # Y Limits Section
    Label(main_layout[6, 3], "Y Limits", fontsize = slider_font, font = :bold)
    y_limits_layout = GridLayout(main_layout[7, 3], 
                                tellwidth = false, 
                                colgap = 10)
    ymin_input = Textbox(y_limits_layout[1, 1], 
                        placeholder = "min", 
                        fontsize = slider_font,
                        width = 90,
                        height = input_height)
    ymax_input = Textbox(y_limits_layout[1, 2], 
                        placeholder = "max", 
                        fontsize = slider_font,
                        width = 90,
                        height = input_height)
    
    # Z Limits Section
    Label(main_layout[8, 3], "Z Limits", fontsize = slider_font, font = :bold)
    z_limits_layout = GridLayout(main_layout[9, 3], 
                                tellwidth = false, 
                                colgap = 10)
    zmin_input = Textbox(z_limits_layout[1, 1], 
                        placeholder = "min", 
                        fontsize = slider_font,
                        width = 90,
                        height = input_height)
    zmax_input = Textbox(z_limits_layout[1, 2], 
                        placeholder = "max", 
                        fontsize = slider_font,
                        width = 90,
                        height = input_height)
    
    # Baseline Section
    Label(main_layout[10, 3], "Baseline", fontsize = slider_font, font = :bold)
    baseline_layout = GridLayout(main_layout[11, 3], 
                                tellwidth = false, 
                                colgap = 10)
    baseline_start = Textbox(baseline_layout[1, 1], 
                            placeholder = "start", 
                            fontsize = slider_font,
                            width = 90,
                            height = input_height)
    baseline_end = Textbox(baseline_layout[1, 2], 
                          placeholder = "end", 
                          fontsize = slider_font,
                          width = 90,
                          height = input_height)
    
    # Baseline Type
    Label(main_layout[12, 3], "Baseline Type", fontsize = slider_font, font = :bold)
    baseline_type = Menu(main_layout[13, 3], 
                        options = ["Select", "absolute", "relative", "relchange", "perchange", "db"],
                        width = input_width,
                        height = input_height)
    
    # Action Buttons Section - Outside the settings panel
    # Create a button layout for the third column
    button_layout = GridLayout(main_layout[15:17, 3], 
                              tellwidth = false, 
                              colgap = 15,
                              rowgap = 20)
    
    # First row of buttons
    examples_button = Button(button_layout[1, 1], 
                           label = "Examples",
                           fontsize = button_font,
                           width = 120,
                           height = input_height)
    
    clear_button = Button(button_layout[1, 2], 
                        label = "Clear GUI",
                        fontsize = button_font,
                        width = 120,
                        height = input_height)
    
    # Second row of buttons
    export_button = Button(button_layout[2, 1], 
                          label = "Export",
                          fontsize = button_font,
                          width = 120,
                          height = input_height)
    
    cursor_button = Button(button_layout[2, 2], 
                          label = "Cursor",
                          fontsize = button_font,
                          width = 120,
                          height = input_height)
    
    # Third row - Main plot button
    plot_button = Button(button_layout[3, 1:2], 
                        label = "PLOT",
                        fontsize = title_font,
                        buttoncolor = :lightblue,
                        width = 255,
                        height = 50)
    
    # Data structure to store GUI state
    gui_state = (
        filetype = Observable("select"),
        bdf_filename = Observable(""),
        participant = Observable(""),
        condition = Observable(""),
        additional = Observable(""),
        plottype = Observable("select"),
        layout = Observable("select"),
        electrodes = Observable(Int[]),
        xlim = Observable((nothing, nothing)),
        ylim = Observable((nothing, nothing)),
        zlim = Observable((nothing, nothing)),
        xtopo = Observable(nothing),
        baseline_start = Observable(nothing),
        baseline_end = Observable(nothing),
        baseline_type = Observable("select")
    )
    
    # Callback functions
    function update_filetype(selection)
        gui_state.filetype[] = filetype_options[selection]
        # Update electrode labels based on file type
        update_electrode_labels(selection)
    end
    
    function update_plottype(selection)
        gui_state.plottype[] = plottype_options[selection]
    end
    
    function update_layout(selection)
        gui_state.layout[] = layout_options[selection]
    end
    
    function update_electrode_labels(filetype_selection)
        # This would be implemented to load electrode labels based on file type
        # For now, just show a placeholder
        electrode_menu.options = ["Select", "Fp1", "Fp2", "F3", "F4", "C3", "C4", "P3", "P4", "O1", "O2"]
    end
    
    function execute_plot()
        println("Executing plot with:")
        println("  File Type: $(gui_state.filetype[])")
        println("  BDF File: $(gui_state.bdf_filename[])")
        println("  Plot Type: $(gui_state.plottype[])")
        println("  Layout: $(gui_state.layout[])")
        println("  Electrodes: $(gui_state.electrodes[])")
        println("  X Limits: $(gui_state.xlim[])")
        println("  Y Limits: $(gui_state.ylim[])")
        
        # Create a new figure and plot 1:10
        try
            # Create a new figure
            fig = Figure(size = (800, 600), title = "Test Plot - 1:10")
            
            # Simple plot of 1:10
            x = 1:10
            y = collect(1:10)
            
            ax = Axis(fig[1, 1], 
                     title = "Plot of 1:10",
                     xlabel = "X",
                     ylabel = "Y")
            
            lines!(ax, x, y, color = :blue, linewidth = 2)
            scatter!(ax, x, y, color = :red, markersize = 8)
            
            # Add styling
            ax.xgridvisible = true
            ax.ygridvisible = true
            ax.xgridcolor = :gray
            ax.ygridcolor = :gray
            ax.xgridwidth = 0.5
            ax.ygridwidth = 0.5
            
            # Create a new screen (window) for the plot using your method
            new_screen = getfield(Main, :GLMakie).Screen()
            display(new_screen, fig)
            
            println("Plot created and displayed in new window!")
            
        catch e
            println("Error creating plot: $e")
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
        # For now, just store single selection
        # In a full implementation, you'd handle multi-selection differently
        gui_state.electrodes[] = [selection]
    end
    
    on(participant_input.stored_string) do value
        gui_state.participant[] = value
    end
    
    on(condition_input.stored_string) do value
        gui_state.condition[] = value
    end
    
    on(bdf_filename_input.stored_string) do value
        gui_state.bdf_filename[] = value
    end
    
    on(xmin_input.stored_string) do value
        try
            gui_state.xlim[] = (parse(Float64, value), gui_state.xlim[][2])
        catch
            gui_state.xlim[] = (nothing, gui_state.xlim[][2])
        end
    end
    
    on(xmax_input.stored_string) do value
        try
            gui_state.xlim[] = (gui_state.xlim[][1], parse(Float64, value))
        catch
            gui_state.xlim[] = (gui_state.xlim[][1], nothing)
        end
    end
    
    on(ymin_input.stored_string) do value
        try
            gui_state.ylim[] = (parse(Float64, value), gui_state.ylim[][2])
        catch
            gui_state.ylim[] = (nothing, gui_state.ylim[][2])
        end
    end
    
    on(ymax_input.stored_string) do value
        try
            gui_state.ylim[] = (gui_state.ylim[][1], parse(Float64, value))
        catch
            gui_state.ylim[] = (gui_state.ylim[][1], nothing)
        end
    end
    
    on(plot_button.clicks) do _
        execute_plot()
    end
    
    on(clear_button.clicks) do _
        # Reset all inputs
        filetype_dropdown.selection = 1
        plottype_dropdown.selection = 1
        layout_dropdown.selection = 1
        electrode_menu.selection = 1
        participant_input.stored_string = ""
        condition_input.stored_string = ""
        bdf_filename_input.stored_string = ""
        xmin_input.stored_string = ""
        xmax_input.stored_string = ""
        ymin_input.stored_string = ""
        ymax_input.stored_string = ""
        zmin_input.stored_string = ""
        zmax_input.stored_string = ""
        xtopo_input.stored_string = ""
        baseline_start.stored_string = ""
        baseline_end.stored_string = ""
        baseline_type.selection = 1
        additional_label.text = ""
    end
    
    # Display the figure
    display(fig)
    
    return fig
end
