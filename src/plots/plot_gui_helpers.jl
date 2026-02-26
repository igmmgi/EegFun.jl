# Utility functions for plot_gui helper functions

"""
    _validate_file(gui_state, required_ext::String)

Validates that a file is specified and has the required extension.
Returns the filename if valid, otherwise shows an error dialog.
"""
function _validate_file(gui_state, required_ext::String)
    gui_state.filename[] == "" && @minimal_error "Error: No file specified!"

    file_ext = lowercase(splitext(gui_state.filename[])[2])
    if file_ext != required_ext
        @minimal_error "Error: This plot requires $required_ext file format"
    end

    return gui_state.filename[]
end

"""
    _gui_selected_channels(gui_state)

Returns the selected channels from GUI state, or all channels if none selected.
"""
function _gui_selected_channels(gui_state)
    isempty(gui_state.electrodes[]) ? channels() : channels(gui_state.electrodes[])
end

"""
    _handle_plot_error(e, plot_type::String)

Standardized error handling for plot functions.
"""
function _handle_plot_error(e, plot_type::String)
    println("Error creating $plot_type plot: $e")
    showerror(stdout, e, catch_backtrace())
end

# Helper functions for creating plots from GUI

function _plot_erp_image(gui_state)
    _validate_file(gui_state, ".jld2")

    try
        selected_channels = _gui_selected_channels(gui_state)

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
        _handle_plot_error(e, "ERP Image")
    end
end

function _plot_gfp(gui_state)
    _validate_file(gui_state, ".jld2")

    try
        selected_channels = _gui_selected_channels(gui_state)

        @async begin
            plot_gfp(gui_state.filename[]; channel_selection = selected_channels, xlim = gui_state.xlim[], ylim = gui_state.ylim[])
        end
    catch e
        _handle_plot_error(e, "GFP")
    end
end

function _plot_time_frequency(gui_state)
    _validate_file(gui_state, ".jld2")

    try
        data = read_data(gui_state.filename[])
        isnothing(data) && @minimal_error "Error: No data found in file"

        selected_channels = _gui_selected_channels(gui_state)

        @async begin
            plot_tf(data; channel_selection = selected_channels, xlim = gui_state.xlim[], colorbar_limits = gui_state.zlim[])
        end
    catch e
        _handle_plot_error(e, "Time-Frequency")
    end
end

function _plot_power_spectrum(gui_state)
    _validate_file(gui_state, ".jld2")

    try
        data = read_data(gui_state.filename[])
        isnothing(data) && @minimal_error "Error: No data found in file"

        selected_channels = _gui_selected_channels(gui_state)

        @async begin
            plot_frequency_spectrum(data; channel_selection = selected_channels, xlim = gui_state.xlim[], ylim = gui_state.ylim[])
        end
    catch e
        _handle_plot_error(e, "Power Spectrum")
    end
end

function _plot_ica(gui_state)
    _validate_file(gui_state, ".jld2")

    try
        data = read_data(gui_state.filename[])
        isnothing(data) && @minimal_error "Error: No ICA data found in file"

        @async begin
            plot_topography(data)
        end
    catch e
        _handle_plot_error(e, "ICA Components")
    end
end

function _plot_filter(gui_state)
    _validate_file(gui_state, ".jld2")

    try
        data = read_data(gui_state.filename[])
        isnothing(data) && @minimal_error "Error: No data found in file"

        @async begin
            plot_filter(data)
        end
    catch e
        _handle_plot_error(e, "Filter Response")
    end
end

function _plot_artifacts(gui_state)
    _validate_file(gui_state, ".jld2")

    try
        data = read_data(gui_state.filename[])
        isnothing(data) && @minimal_error "Error: No data found in file"

        @async begin
            plot_artifact_detection(data)
        end
    catch e
        _handle_plot_error(e, "Artifact Detection")
    end
end

function _plot_triggers(gui_state)
    gui_state.filename[] == "" && @minimal_error "Error: No file specified!"

    file_ext = lowercase(splitext(gui_state.filename[])[2])

    try
        if file_ext == ".jld2"
            data = read_data(gui_state.filename[])
            isnothing(data) && @minimal_error "Error: No data found in file"

            @async begin
                plot_trigger_overview(data)
            end
        elseif file_ext == ".bdf"
            gui_state.layout_file[] == "" && @minimal_error "Error: No layout file selected!"
            layout = read_layout(gui_state.layout_file[])
            polar_to_cartesian_xy!(layout)
            dat = read_raw_data(gui_state.filename[])
            dat = create_eegfun_data(dat, layout)

            @async begin
                plot_trigger_overview(dat)
            end
        else
            @minimal_error "Error: Trigger plot supports BDF and JLD2 formats"
        end
    catch e
        _handle_plot_error(e, "Triggers")
    end
end

function _plot_correlation(gui_state)
    _validate_file(gui_state, ".jld2")

    try
        data = read_data(gui_state.filename[])
        isnothing(data) && @minimal_error "Error: No data found in file"

        @async begin
            plot_correlation_heatmap(data)
        end
    catch e
        _handle_plot_error(e, "Correlation Heatmap")
    end
end

function _plot_channel_summary(gui_state)
    _validate_file(gui_state, ".jld2")

    try
        data = read_data(gui_state.filename[])
        isnothing(data) && @minimal_error "Error: No data found in file"

        @async begin
            plot_channel_summary(data)
        end
    catch e
        _handle_plot_error(e, "Channel Summary")
    end
end

function _plot_joint_probability(gui_state)
    _validate_file(gui_state, ".jld2")

    try
        data = read_data(gui_state.filename[])
        isnothing(data) && @minimal_error "Error: No data found in file"

        @async begin
            plot_joint_probability(data)
        end
    catch e
        _handle_plot_error(e, "Joint Probability")
    end
end

function _plot_erp_measurement_gui(gui_state)
    _validate_file(gui_state, ".jld2")

    try
        data = read_data(gui_state.filename[])
        isnothing(data) && @minimal_error "Error: No ERP data found in file"

        @async begin
            plot_erp_measurement_gui(data)
        end
    catch e
        _handle_plot_error(e, "ERP Measurement GUI")
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
        _handle_plot_error(e, "Layout")
    end
end
