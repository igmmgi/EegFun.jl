# Help system for interactive plots
# Provides keyboard shortcut information when 'i' key is pressed

"""
    PlotHelpInfo

Structure to hold help information for different plot types.
"""
struct PlotHelpInfo
    title::String
    interactions::Vector{Tuple{String, String}}
end

"""
    get_plot_help_info(plot_type::Symbol) -> PlotHelpInfo

Get help information for a specific plot type.
"""
function get_plot_help_info(plot_type::Symbol)::PlotHelpInfo
    help_info = Dict(
        :erp => PlotHelpInfo(
            "ERP Plot",
            [
                ("Up", "Zoom in on Y-axis (compress Y limits)"),
                ("Down", "Zoom out on Y-axis (expand Y limits)"),
                ("Left", "Zoom in on X-axis (compress time range)"),
                ("Right", "Zoom out on X-axis (expand time range)"),
                ("Shift + Left Click + Drag", "Select a time region (blue rectangle)"),
                ("Right Click on Selection", "Context menu"),
                ("Ctrl + Left Click + Drag", "Select channels for analysis"),
            ]
        ),
        
        :epochs => PlotHelpInfo(
            "Epochs Plot",
            [
                ("Up", "Zoom in on Y-axis (compress Y limits)"),
                ("Down", "Zoom out on Y-axis (expand Y limits)"),
                ("Left", "Zoom in on X-axis (compress time range)"),
                ("Right", "Zoom out on X-axis (expand time range)"),
                ("Shift + Left Click + Drag", "Select a time region (blue rectangle)"),
                ("Right Click on Selection", "Open context menu with plot options"),
                ("Ctrl + Left Click + Drag", "Select channels for analysis"),
            ]
        ),
        
        :databrowser => PlotHelpInfo(
            "Data Browser",
            [
                ("Up", "Zoom in on Y-axis (amplitude scale)"),
                ("Down", "Zoom out on Y-axis (amplitude scale)"),
                ("Left", "Scroll backward in time"),
                ("Right", "Scroll forward in time"),
                ("Shift + Left Click + Drag", "Select x-regions"),
                ("Shift + Left Click", "Clear selections"),
                ("Right Click within highlighted region", "Context menu"),
                ("Strg/Crtl + Left click", "Select closest channel"),
            ]
        ),
        
        :erp_image => PlotHelpInfo(
            "ERP Image",
            [
                ("Up", "Zoom in on Y-axis (epochs)"),
                ("Down", "Zoom out on Y-axis (epochs)"),
                ("Left", "Zoom in on X-axis (time)"),
                ("Right", "Zoom out on X-axis (time)"),
                ("Hover", "Context-aware controls based on mouse position"),
                ("Left Click + Drag", "Pan the view"),
            ]
        ),
        
        :ica => PlotHelpInfo(
            "ICA Components",
            [
                ("Up", "Zoom in on Y-axis"),
                ("Down", "Zoom out on Y-axis"),
                ("Left", "Scroll backward in time"),
                ("Right", "Scroll forward in time"),
                ("Shift + ↑", "Previous component"),
                ("Shift + ↓", "Next component"),
                ("Left Click", "Select/deselect components"),
                ("Ctrl + Left Click", "Multi-select components"),
            ]
        ),
        
        :topography => PlotHelpInfo(
            "Topographic Plot",
            [
                ("Up", "Increase scaling"),
                ("Down", "Decrease scaling"),
                ("Shift + Left Click + Drag", "Select region(s)"),
                ("Left Click", "Select/deselect"),
                ("Right Click", "Context menu"),
            ]
        ),
        
        :power_spectrum => PlotHelpInfo(
            "Power Spectrum",
            [
                ("Checkbox Toggle", "Switch between linear and log scales"),
                ("Left Click + Drag", "Zoom into frequency range"),
            ]
        ),
        
        :triggers => PlotHelpInfo(
            "Trigger Plot",
            [
                ("Slider Drag", "Adjust time window position"),
                ("Slider Drag", "Adjust window size"),
                ("Left Click + Drag", "Pan the view"),
            ]
        )
    )
    
    return get(help_info, plot_type, PlotHelpInfo(
        "Plot",
        []
    ))
end

"""
    print_plot_help(plot_type::Symbol)

Print help information for a specific plot type to the console.
"""
function print_plot_help(plot_type::Symbol)
    help_info = get_plot_help_info(plot_type)
    
    # If no interactions, do nothing
    if isempty(help_info.interactions)
        return
    end
    
    println("\n" * "="^40)
    println("📊 $(help_info.title)")
    println("="^40)
    
    for (action, description) in help_info.interactions
        println("  $action  →  $description")
    end
    
    println("="^40 * "\n")
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
