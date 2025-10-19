# Help system for interactive plots
# Provides keyboard shortcut information when 'i' key is pressed

"""
    PlotHelpInfo

Structure to hold help information for different plot types.
"""
struct PlotHelpInfo
    plot_type::String
    title::String
    keyboard_shortcuts::Vector{Tuple{String, String}}
    mouse_interactions::Vector{Tuple{String, String}}
    additional_info::Vector{String}
end

"""
    get_plot_help_info(plot_type::Symbol) -> PlotHelpInfo

Get help information for a specific plot type.
"""
function get_plot_help_info(plot_type::Symbol)::PlotHelpInfo
    help_info = Dict(
        :erp => PlotHelpInfo(
            "ERP Plot",
            "Event-Related Potential Visualization",
            [
                ("↑ (Up Arrow)", "Zoom in on Y-axis (compress Y limits)"),
                ("↓ (Down Arrow)", "Zoom out on Y-axis (expand Y limits)"),
                ("← (Left Arrow)", "Zoom in on X-axis (compress time range)"),
                ("→ (Right Arrow)", "Zoom out on X-axis (expand time range)"),
                ("i", "Show this help information"),
            ],
            [
                ("Shift + Left Click + Drag", "Select a time region (blue rectangle)"),
                ("Right Click on Selection", "Open context menu with plot options"),
                ("Ctrl + Left Click + Drag", "Select channels for analysis"),
            ],
            [
                "Use arrow keys to navigate and zoom the plot",
                "Time selections can be used for topographic analysis",
                "Channel selections persist across different time windows"
            ]
        ),
        
        :epochs => PlotHelpInfo(
            "Epochs Plot", 
            "Epoch-based Data Visualization",
            [
                ("↑ (Up Arrow)", "Zoom in on Y-axis (compress Y limits)"),
                ("↓ (Down Arrow)", "Zoom out on Y-axis (expand Y limits)"),
                ("← (Left Arrow)", "Zoom in on X-axis (compress time range)"),
                ("→ (Right Arrow)", "Zoom out on X-axis (expand time range)"),
                ("i", "Show this help information"),
            ],
            [
                ("Shift + Left Click + Drag", "Select a time region (blue rectangle)"),
                ("Right Click on Selection", "Open context menu with plot options"),
                ("Ctrl + Left Click + Drag", "Select channels for analysis"),
            ],
            [
                "Epochs are displayed with individual trial data",
                "Use selections to analyze specific time windows or channels",
                "Right-click selections for additional analysis options"
            ]
        ),
        
        :databrowser => PlotHelpInfo(
            "Data Browser",
            "Interactive Data Exploration Tool",
            [
                ("↑ (Up Arrow)", "Zoom in on Y-axis (amplitude scale)"),
                ("↓ (Down Arrow)", "Zoom out on Y-axis (amplitude scale)"),
                ("← (Left Arrow)", "Scroll backward in time"),
                ("→ (Right Arrow)", "Scroll forward in time"),
                ("i", "Show this help information"),
            ],
            [
                ("Left Click", "Clear all selections"),
                ("Ctrl + Left Click + Drag", "Select channels for analysis"),
                ("Shift + Left Click + Drag", "Select time region"),
                ("Right Click", "Context menu (planned)"),
            ],
            [
                "Use the control panel on the right to adjust filters and settings",
                "Channel selections are highlighted in the plot",
                "Time selections can be used for further analysis",
                "Toggle butterfly mode to overlay all channels"
            ]
        ),
        
        :erp_image => PlotHelpInfo(
            "ERP Image",
            "Event-Related Potential Heatmap Visualization",
            [
                ("↑ (Up Arrow)", "Zoom in on Y-axis (epochs)"),
                ("↓ (Down Arrow)", "Zoom out on Y-axis (epochs)"),
                ("← (Left Arrow)", "Zoom in on X-axis (time)"),
                ("→ (Right Arrow)", "Zoom out on X-axis (time)"),
                ("i", "Show this help information"),
            ],
            [
                ("Mouse Hover", "Context-aware controls based on mouse position"),
                ("Left Click + Drag", "Pan the view"),
            ],
            [
                "Controls are context-aware based on mouse position",
                "Hover over different axes to control different aspects",
                "ERP image shows trial-by-trial data as heatmap"
            ]
        ),
        
        :ica => PlotHelpInfo(
            "ICA Components",
            "Independent Component Analysis Visualization",
            [
                ("↑ (Up Arrow)", "Zoom in on Y-axis"),
                ("↓ (Down Arrow)", "Zoom out on Y-axis"),
                ("← (Left Arrow)", "Scroll backward in time"),
                ("→ (Right Arrow)", "Scroll forward in time"),
                ("Shift + ↑", "Previous component"),
                ("Shift + ↓", "Next component"),
                ("i", "Show this help information"),
            ],
            [
                ("Left Click", "Select/deselect components"),
                ("Ctrl + Left Click", "Multi-select components"),
            ],
            [
                "Use Shift + arrow keys to navigate between components",
                "Selected components are highlighted",
                "Time window can be adjusted with left/right arrows",
                "Topographic maps show component spatial distribution"
            ]
        ),
        
        :topography => PlotHelpInfo(
            "Topographic Plot",
            "Electrode Layout and Spatial Distribution",
            [
                ("i", "Show this help information"),
            ],
            [
                ("Left Click", "Select/deselect electrodes"),
                ("Shift + Left Click", "Multi-select electrodes"),
                ("Mouse Hover", "Highlight electrode connections"),
            ],
            [
                "Click electrodes to select them for analysis",
                "Neighbor connections are shown on hover",
                "Selected electrodes are highlighted in the plot"
            ]
        ),
        
        :power_spectrum => PlotHelpInfo(
            "Power Spectrum",
            "Frequency Domain Analysis",
            [
                ("i", "Show this help information"),
            ],
            [
                ("Checkbox Toggle", "Switch between linear and log scales"),
                ("Left Click + Drag", "Zoom into frequency range"),
            ],
            [
                "Use checkboxes to toggle between linear and log scales",
                "Interactive controls allow real-time scale changes",
                "Frequency range can be adjusted by zooming"
            ]
        ),
        
        :triggers => PlotHelpInfo(
            "Trigger Plot",
            "Event Timing and Trigger Visualization",
            [
                ("i", "Show this help information"),
            ],
            [
                ("Slider Drag", "Adjust time window position"),
                ("Slider Drag", "Adjust window size"),
                ("Left Click + Drag", "Pan the view"),
            ],
            [
                "Use sliders to control time window and position",
                "Vertical lines indicate trigger events",
                "Window size controls the visible time range"
            ]
        )
    )
    
    return get(help_info, plot_type, PlotHelpInfo(
        "Unknown Plot",
        "Interactive Plot",
        [("i", "Show this help information")],
        [],
        ["This plot type has basic interactivity"]
    ))
end

"""
    print_plot_help(plot_type::Symbol)

Print help information for a specific plot type to the console.
"""
function print_plot_help(plot_type::Symbol)
    help_info = get_plot_help_info(plot_type)
    
    println("\n" * "="^60)
    println("📊 $(help_info.title) - Interactive Controls")
    println("="^60)
    
    if !isempty(help_info.keyboard_shortcuts)
        println("\n⌨️  KEYBOARD SHORTCUTS:")
        println("-"^30)
        for (key, description) in help_info.keyboard_shortcuts
            println("  $key  →  $description")
        end
    end
    
    if !isempty(help_info.mouse_interactions)
        println("\n🖱️  MOUSE INTERACTIONS:")
        println("-"^30)
        for (action, description) in help_info.mouse_interactions
            println("  $action  →  $description")
        end
    end
    
    if !isempty(help_info.additional_info)
        println("\n💡 ADDITIONAL INFO:")
        println("-"^30)
        for info in help_info.additional_info
            println("  • $info")
        end
    end
    
    println("\n" * "="^60)
    println("Press 'i' again to hide this help")
    println("="^60 * "\n")
end

"""
    setup_help_interaction!(fig::Figure, plot_type::Symbol)

Set up help interaction for a figure. When 'i' key is pressed, shows help information.
"""
function setup_help_interaction!(fig::Figure, plot_type::Symbol)
    help_visible = Ref(false)
    
    on(events(fig).keyboardbutton) do event
        if event.action == Keyboard.press && event.key == Keyboard.i
            if help_visible[]
                # Hide help (clear console)
                print("\r" * " "^100 * "\r")  # Clear line
                help_visible[] = false
            else
                # Show help
                print_plot_help(plot_type)
                help_visible[] = true
            end
        end
    end
end

"""
    show_plot_help(plot_type::Symbol)

Convenience function to show help for a plot type without setting up interaction.
"""
function show_plot_help(plot_type::Symbol)
    print_plot_help(plot_type)
end
