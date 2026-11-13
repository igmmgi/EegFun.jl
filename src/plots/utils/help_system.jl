struct PlotHelpInfo
    title::String
    interactions::Vector{Tuple{String,String}}
end

const _PLOT_HELP_INFO = Dict{Symbol,PlotHelpInfo}(
    :erp => PlotHelpInfo(
        "ERP Plot",
        [
            ("Up", "Zoom in on Y-axis (compress Y limits)"),
            ("Down", "Zoom out on Y-axis (expand Y limits)"),
            ("Left", "Zoom in on X-axis (compress time range)"),
            ("Right", "Zoom out on X-axis (expand time range)"),
            ("Left Click + Drag", "Pan the view"),
            ("Shift + Left Click + Drag", "Select a time region (blue rectangle)"),
            ("Right Click on Selection", "Context menu"),
            ("Ctrl + Left Click + Drag", "Select channels for analysis"),
            ("C", "Open control panel (baseline and condition toggles)"),
            ("L", "Toggle legend on/off"),
        ],
    ),
    :epochs => PlotHelpInfo(
        "Epochs Plot",
        [
            ("Up", "Zoom in on Y-axis (compress Y limits)"),
            ("Down", "Zoom out on Y-axis (expand Y limits)"),
            ("Left", "Zoom in on X-axis (compress time range)"),
            ("Right", "Zoom out on X-axis (expand time range)"),
            ("Left Click + Drag", "Pan the view"),
            ("Shift + Left Click + Drag", "Select a time region (blue rectangle)"),
            ("Right Click on Selection", "Open context menu with plot options"),
            ("Ctrl + Left Click + Drag", "Select channels for analysis"),
            ("C", "Open control panel (baseline and condition toggles)"),
            ("L", "Toggle legend on/off"),
        ],
    ),
    :databrowser => PlotHelpInfo(
        "Data Browser",
        [
            ("Up", "Zoom in on Y-axis (amplitude scale)"),
            ("Down", "Zoom out on Y-axis (amplitude scale)"),
            ("Left", "Scroll backward in time"),
            ("Right", "Scroll forward in time"),
            ("PageUp / PageDown", "Step backward / forward by a full window (or 5 epochs)"),
            ("Mouse Wheel", "Zoom in/out"),
            ("Right Mouse + Drag", "Pan view (when outside region)"),
            ("Shift + Left Click + Drag", "Select x-region (multiple selections allowed)"),
            ("Left Click on region", "Deselect that region"),
            ("Shift + Left Click", "Clear current selection"),
            ("C", "Clear all selected regions and channel selection"),
            ("Right Click within highlighted region", "Context menu (Plot Epochs, Plot ERP, Topoplot, Spectrum)"),
            ("Ctrl + Left click", "Select/toggle closest channel"),
            ("R", "Open channel repair menu (continuous data only)"),
        ],
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
        ],
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
        ],
    ),
    :topography => PlotHelpInfo(
        "Topographic Plot",
        [
            ("Up", "Increase scaling"),
            ("Down", "Decrease scaling"),
            ("Shift + Left Click + Drag", "Select region(s)"),
            ("Left Click", "Select/deselect"),
            ("Right Click", "Context menu"),
        ],
    ),
    :power_spectrum => PlotHelpInfo(
        "Power Spectrum",
        [("Checkbox Toggle", "Switch between linear and log scales"), ("Left Click + Drag", "Zoom into frequency range")],
    ),
    :trigger => PlotHelpInfo(
        "Trigger Plot",
        [("Slider Drag", "Adjust time interval position"), ("Slider Drag", "Adjust interval size"), ("Left Click + Drag", "Pan the view")],
    ),
)

const _EMPTY_PLOT_HELP = PlotHelpInfo("Plot", Tuple{String,String}[])

_get_plot_help_info(plot_type::Symbol)::PlotHelpInfo = get(_PLOT_HELP_INFO, plot_type, _EMPTY_PLOT_HELP)

function _print_plot_help(plot_type::Symbol)
    help_info = _get_plot_help_info(plot_type)
    isempty(help_info.interactions) && return
    println("\n" * "="^40)
    println(" $(help_info.title)")
    println("="^40)
    for (action, description) in help_info.interactions
        println("  $action  →  $description")
    end
    println("="^40 * "\n")
end

function _setup_help_interaction!(fig::Figure, plot_type::Symbol)
    on(events(fig).keyboardbutton) do event
        if event.action == Keyboard.press && event.key == Keyboard.i
            _print_plot_help(plot_type)
        end
    end
end

_show_plot_help(plot_type::Symbol) = _print_plot_help(plot_type)
