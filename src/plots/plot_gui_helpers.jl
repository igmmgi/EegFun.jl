# Utility functions for plot_gui helper functions

"""
    validate_file(gui_state, required_ext::String)

Validates that a file is specified and has the required extension.
Returns the filename if valid, otherwise shows an error dialog.
"""
function validate_file(gui_state, required_ext::String)
    gui_state.filename[] == "" && @minimal_error "Error: No file specified!"

    file_ext = lowercase(splitext(gui_state.filename[])[2])
    if file_ext != required_ext
        @minimal_error "Error: This plot requires $required_ext file format"
    end

    return gui_state.filename[]
end

"""
    get_selected_channels(gui_state)

Returns the selected channels from GUI state, or all channels if none selected.
"""
function get_selected_channels(gui_state)
    isempty(gui_state.electrodes[]) ? channels() : channels(gui_state.electrodes[])
end

"""
    handle_plot_error(e, plot_type::String)

Standardized error handling for plot functions.
"""
function handle_plot_error(e, plot_type::String)
    println("Error creating $plot_type plot: $e")
    showerror(stdout, e, catch_backtrace())
end

# Helper functions for creating plots from GUI

function _plot_erp_image(gui_state)
    validate_file(gui_state, ".jld2")

    try
        selected_channels = get_selected_channels(gui_state)

        @async begin
            layout_sym = Symbol(gui_state.layout_type[])
            plot_erp_image(
                gui_state.filename[];
                layout = layout_sym,
                channel_selection = selected_channels,
                xlim = gui_state.xlim[],
                colorbar_limits = gui_state.zlim[],
            )
        end
    catch e
        handle_plot_error(e, "ERP Image")
    end
end

function _plot_gfp(gui_state)
    validate_file(gui_state, ".jld2")

    try
        selected_channels = get_selected_channels(gui_state)

        @async begin
            plot_gfp(gui_state.filename[]; channel_selection = selected_channels, xlim = gui_state.xlim[], ylim = gui_state.ylim[])
        end
    catch e
        handle_plot_error(e, "GFP")
    end
end

function _plot_time_frequency(gui_state)
    validate_file(gui_state, ".jld2")

    try
        data = read_data(gui_state.filename[])
        isnothing(data) && @minimal_error "Error: No data found in file"

        selected_channels = get_selected_channels(gui_state)

        @async begin
            plot_time_frequency(data; channel_selection = selected_channels, xlim = gui_state.xlim[], colorbar_limits = gui_state.zlim[])
        end
    catch e
        handle_plot_error(e, "Time-Frequency")
    end
end

function _plot_power_spectrum(gui_state)
    validate_file(gui_state, ".jld2")

    try
        data = read_data(gui_state.filename[])
        isnothing(data) && @minimal_error "Error: No data found in file"

        selected_channels = get_selected_channels(gui_state)

        @async begin
            plot_frequency_spectrum(data; channel_selection = selected_channels, xlim = gui_state.xlim[], ylim = gui_state.ylim[])
        end
    catch e
        handle_plot_error(e, "Power Spectrum")
    end
end

function _plot_ica(gui_state)
    @minimal_error "Error: ICA Components plot requires ICA results.\nPlease load a file with ICA decomposition data."
end

function _plot_filter(gui_state)
    @minimal_error "Error: Filter Response plot requires filter parameters.\nThis feature is not yet implemented in the GUI."
end

function _plot_artifacts(gui_state)
    @minimal_error "Error: Artifact Detection plot requires pre-detected artifacts.\nPlease run artifact detection first and save the results."
end

function _plot_triggers(gui_state)
    gui_state.filename[] == "" && @minimal_error "Error: No file specified!"
    gui_state.layout_file[] == "" && @minimal_error "Error: No layout file selected!"

    file_ext = lowercase(splitext(gui_state.filename[])[2])
    if file_ext ∉ [".bdf"]
        @minimal_error "Error: Trigger plot currently only supports BDF format"
    end

    try
        layout = read_layout(gui_state.layout_file[])
        polar_to_cartesian_xy!(layout)

        dat = read_raw_data(gui_state.filename[])
        dat = create_eeg_dataframe(dat, layout)

        @async begin
            plot_trigger_overview(dat)
        end
    catch e
        handle_plot_error(e, "Triggers")
    end
end

function _plot_correlation(gui_state)
    validate_file(gui_state, ".jld2")

    try
        data = read_data(gui_state.filename[])
        isnothing(data) && @minimal_error "Error: No data found in file"

        @async begin
            plot_correlation_heatmap(data)
        end
    catch e
        handle_plot_error(e, "Correlation Heatmap")
    end
end

function _plot_layout(gui_state)
    gui_state.layout_file[] == "" && @minimal_error "Error: No layout file selected!"

    try
        layout = read_layout(gui_state.layout_file[])

        @async begin
            plot_layout_2d(layout)
        end
    catch e
        handle_plot_error(e, "Layout")
    end
end
