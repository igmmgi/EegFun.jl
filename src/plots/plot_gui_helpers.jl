# Utility functions for plot_gui helper functions

"""
    _validate_file(gui_state, required_ext::String)

Validates that a file is specified and has the required extension.
Returns the filename if valid, otherwise shows an error dialog.
"""
function _validate_file(gui_state, required_ext::String = ".jld2")
    if gui_state.filename[] == ""
        @minimal_warning "Requested plot settings incompatible: recheck!"
        return nothing
    end

    file_ext = lowercase(splitext(gui_state.filename[])[2])
    # Accept: exact extension match OR no extension (treated as filename pattern)
    if file_ext != required_ext && !isempty(file_ext)
        @minimal_warning "Requested plot settings incompatible: recheck!"
        return nothing
    end

    return gui_state.filename[]
end

"""
    _gui_selected_channels(gui_state)

Returns the selected channels from GUI state, or all channels if none selected.
"""
function _gui_selected_channels(gui_state)
    isempty(gui_state.electrodes[]) ? channels() : channels(Symbol.(gui_state.electrodes[]))
end

"""
    _handle_plot_error(e, plot_type::String)

Standardized error handling for plot functions.
"""
function _handle_plot_error(e, plot_type::String)
    if e isa EegFunError
        # Already logged by @minimal_error inside the plot function — don't repeat it
        return
    end
    # Unexpected exception: one clean line, no stacktrace
    @minimal_warning "Unexpected error creating $plot_type plot: $(sprint(showerror, e))"
end

# Helper functions for creating plots from GUI

"""
    _parse_gui_int_field(str) -> Union{Vector{Int}, Nothing}

Parse a space- or comma-separated string of integers from a GUI text field.
Supports range syntax: `"1:50"` expands to integers 1…50.
Mixed input like `"1 3:5 7"` → `[1, 3, 4, 5, 7]` is also valid.
- Returns `nothing` if the field is in an invalid state (sentinel `__INVALID__`).
- Returns `Int[]` if the field is blank.
- Returns the parsed integer vector otherwise.
"""
function _parse_gui_int_field(str::String)
    str == "__INVALID__" && return nothing
    s = strip(str)
    isempty(s) && return Int[]
    result = Int[]
    for token in filter(!isempty, split(s, r"[\s,]+"))
        if occursin(':', token)
            parts = split(token, ':')
            length(parts) == 2 || return nothing   # malformed range
            lo = tryparse(Int, parts[1])
            hi = tryparse(Int, parts[2])
            (isnothing(lo) || isnothing(hi)) && return nothing
            append!(result, lo:hi)
        else
            v = tryparse(Int, token)
            isnothing(v) && return nothing
            push!(result, v)
        end
    end
    return result
end

"""
    _apply_gui_filters(dat, gui_state) -> Union{dat, Nothing}

Apply condition and epoch filters from `gui_state` to `dat`.
Returns `nothing` if any field is currently marked invalid (user has typed bad input).
- Condition filter applies only when `dat` is a `Vector` (multi-condition file).
- Epoch filter applies only when `dat` is `EpochData` or `Vector{EpochData}`.
"""
function _apply_gui_filters(dat, gui_state)
    cond_list  = _parse_gui_int_field(gui_state.condition[])
    epoch_list = _parse_gui_int_field(gui_state.epoch[])

    if isnothing(cond_list) || isnothing(epoch_list)
        @minimal_warning "Requested plot settings incompatible: recheck!"
        return nothing
    end

    if !isempty(cond_list) && dat isa Vector
        dat = subset(dat; condition_selection = conditions(cond_list))
    end

    if !isempty(epoch_list)
        if dat isa EpochData
            dat = subset(dat; epoch_selection = epochs(epoch_list))
        elseif dat isa Vector && eltype(dat) <: EpochData
            dat = subset(dat; epoch_selection = epochs(epoch_list))
        end
    end

    return dat
end

function _plot_erp_image(gui_state)
    isnothing(_validate_file(gui_state)) && return

    fname = gui_state.filename[]
    is_pattern = lowercase(splitext(fname)[2]) != ".jld2"
    selected_channels = _gui_selected_channels(gui_state)

    @async begin
        try
            layout_sym = Symbol(gui_state.layout_type[])
            z = gui_state.zlim[]
            extra_kwargs = (!isnothing(z[1]) && !isnothing(z[2])) ? (colorrange = z,) : NamedTuple()

            if is_pattern
                # Pattern mode — delegate directly; batch discovery handled by string dispatch
                plot_erp_image(
                    fname;
                    input_dir = gui_state.directory[] == "" ? pwd() : gui_state.directory[],
                    layout = layout_sym,
                    channel_selection = selected_channels,
                    xlim = gui_state.xlim[],
                    extra_kwargs...,
                )
            else
                # .jld2 mode — load data and apply GUI filters
                dat = read_data(fname)
                isnothing(dat) && return
                dat = _apply_gui_filters(dat, gui_state)
                isnothing(dat) && return
                plot_erp_image(dat; layout = layout_sym, channel_selection = selected_channels, xlim = gui_state.xlim[], extra_kwargs...)
            end
        catch e
            _handle_plot_error(e, "ERP Image")
        end
    end
end

function _plot_gfp(gui_state)
    isnothing(_validate_file(gui_state, ".jld2")) && return

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
    isnothing(_validate_file(gui_state, ".jld2")) && return

    fname = gui_state.filename[]
    input_dir = gui_state.directory[] == "" ? pwd() : gui_state.directory[]
    is_pattern = lowercase(splitext(fname)[2]) != ".jld2"

    # Shared plot kwargs
    selected_channels = _gui_selected_channels(gui_state)
    x_lo, x_hi = gui_state.xlim[]
    xlim_val = (!isnothing(x_lo) && !isnothing(x_hi)) ? (x_lo, x_hi) : nothing
    z_lo, z_hi = gui_state.zlim[]
    zlim_val = (!isnothing(z_lo) && !isnothing(z_hi)) ? (z_lo, z_hi) : nothing
    b_lo, b_hi = gui_state.baseline_start[], gui_state.baseline_end[]
    baseline_val = (!isnothing(b_lo) && !isnothing(b_hi)) ? (b_lo, b_hi) : nothing
    valid_tf_methods = ("db", "absolute", "relative", "relchange", "percent", "zscore", "normchange", "vssum")
    baseline_method_sym = let bt = gui_state.baseline_type[]
        bt ∈ valid_tf_methods ? Symbol(bt) : :db
    end
    layout_sym = let lt = gui_state.layout_type[]
        lt ∈ ("single", "grid", "topo") ? Symbol(lt) : :single
    end

    # Condition filter: keep only requested conditions (empty = all)
    cond_list = _parse_gui_int_field(gui_state.condition[])
    if isnothing(cond_list)
        @minimal_warning "Requested plot settings incompatible: recheck!"
        return
    end
    function _apply_tf_conditions(data)
        vec = data isa AbstractVector ? data : [data]
        filtered = isempty(cond_list) ? vec : filter(tf -> tf.condition ∈ cond_list, vec)
        isempty(filtered) && @minimal_warning "No conditions matching $cond_list found in file"
        return filtered
    end

    if is_pattern
        part_list = _parse_gui_int_field(gui_state.participant[])
        if isnothing(part_list)
            @minimal_warning "Requested plot settings incompatible: recheck!"
            return
        end
        part_sel = isempty(part_list) ? participants() : participants(part_list)

        # Warn about missing participant IDs
        if !isempty(part_list)
            all_files   = _find_batch_files(fname, input_dir)
            avail_ids   = sort(unique(_extract_participant_id.(all_files)))
            missing_ids = filter(id -> id ∉ avail_ids, part_list)
            if !isempty(missing_ids)
                @minimal_warning "Participant ID(s) $missing_ids not found in '$fname'. Available IDs: $avail_ids"
            end
        end

        files = _find_batch_files(fname, input_dir, part_sel)
        isempty(files) && return

        @async begin
            try
                for file in sort(files, by = _natural_sort_key)
                    @info "Plotting TF: $file"
                    data = read_data(joinpath(input_dir, file))
                    isnothing(data) && continue
                    tf_vec = _apply_tf_conditions(data)
                    isempty(tf_vec) && continue
                    tf = length(tf_vec) == 1 ? only(tf_vec) : tf_vec
                    layout_kw = tf isa TimeFreqData ? (layout = layout_sym,) : (;)
                    plot_tf(
                        tf;
                        layout_kw...,
                        channel_selection = selected_channels,
                        xlim              = xlim_val,
                        colorrange        = zlim_val,
                        baseline_interval = baseline_val,
                        baseline_method   = baseline_method_sym,
                    )
                end
            catch e
                _handle_plot_error(e, "Time-Frequency")
            end
        end
        return
    end

    # .jld2 direct mode
    try
        data = read_data(fname)
        if isnothing(data)
            @minimal_warning "Requested plot settings incompatible: recheck!"
            return
        end
        tf_vec = _apply_tf_conditions(data)
        isempty(tf_vec) && return
        tf = length(tf_vec) == 1 ? only(tf_vec) : tf_vec
        @async begin
            try
                layout_kw = tf isa TimeFreqData ? (layout = layout_sym,) : (;)
                plot_tf(
                    tf;
                    layout_kw...,
                    channel_selection = selected_channels,
                    xlim              = xlim_val,
                    colorrange        = zlim_val,
                    baseline_interval = baseline_val,
                    baseline_method   = baseline_method_sym,
                )
            catch e
                _handle_plot_error(e, "Time-Frequency")
            end
        end
    catch e
        _handle_plot_error(e, "Time-Frequency")
    end
end

function _plot_power_spectrum(gui_state)
    isnothing(_validate_file(gui_state, ".jld2")) && return

    try
        data = read_data(gui_state.filename[])
        if isnothing(data)
            @minimal_warning "Requested plot settings incompatible: recheck!"
            return
        end

        # plot_frequency_spectrum only accepts SpectrumData — unwrap vector if needed
        spectrum = data isa AbstractVector ? first(data) : data
        if !(spectrum isa SpectrumData)
            @minimal_warning "Power Spectrum requires a SpectrumData file (got $(typeof(data))). Select a file produced by compute_spectrum()."
            return
        end

        selected_channels = _gui_selected_channels(gui_state)
        _, x_hi = gui_state.xlim[]   # upper bound → max_freq

        @async begin
            try
                plot_frequency_spectrum(spectrum; channel_selection = selected_channels, max_freq = x_hi)
            catch e
                _handle_plot_error(e, "Power Spectrum")
            end
        end
    catch e
        _handle_plot_error(e, "Power Spectrum")
    end
end

function _plot_ica(gui_state)
    isnothing(_validate_file(gui_state, ".jld2")) && return

    try
        data = read_data(gui_state.filename[])
        if isnothing(data)
            @minimal_warning "Requested plot settings incompatible: recheck!"
            return
        end
        @async begin
            plot_topography(data)
        end
    catch e
        _handle_plot_error(e, "ICA Components")
    end
end

function _plot_filter(gui_state)
    isnothing(_validate_file(gui_state, ".jld2")) && return

    try
        data = read_data(gui_state.filename[])
        if isnothing(data)
            @minimal_warning "Requested plot settings incompatible: recheck!"
            return
        end
        @async begin
            plot_filter(data)
        end
    catch e
        _handle_plot_error(e, "Filter Response")
    end
end

function _plot_artifacts(gui_state)
    isnothing(_validate_file(gui_state, ".jld2")) && return

    try
        data = read_data(gui_state.filename[])
        if isnothing(data)
            @minimal_warning "Requested plot settings incompatible: recheck!"
            return
        end
        @async begin
            plot_artifact_detection(data)
        end
    catch e
        _handle_plot_error(e, "Artifact Detection")
    end
end

function _plot_triggers(gui_state)
    if gui_state.filename[] == ""
        @minimal_warning "Requested plot settings incompatible: recheck!"
        return
    end

    file_ext = lowercase(splitext(gui_state.filename[])[2])

    try
        if file_ext == ".jld2"
            data = read_data(gui_state.filename[])
            if isnothing(data)
                @minimal_warning "Requested plot settings incompatible: recheck!"
                return
            end
            @async begin
                plot_trigger_overview(data)
            end
        elseif file_ext == ".bdf"
            if gui_state.layout_file[] == ""
                @minimal_warning "Requested plot settings incompatible: recheck!"
                return
            end
            layout = read_layout(gui_state.layout_file[])
            polar_to_cartesian_xy!(layout)
            dat = read_raw_data(gui_state.filename[])
            dat = create_eegfun_data(dat, layout)
            @async begin
                plot_trigger_overview(dat)
            end
        else
            @minimal_warning "Requested plot settings incompatible: recheck!"
        end
    catch e
        _handle_plot_error(e, "Triggers")
    end
end

function _plot_correlation(gui_state)
    isnothing(_validate_file(gui_state, ".jld2")) && return

    try
        data = read_data(gui_state.filename[])
        if isnothing(data)
            @minimal_warning "Requested plot settings incompatible: recheck!"
            return
        end
        @async begin
            plot_correlation_heatmap(data)
        end
    catch e
        _handle_plot_error(e, "Correlation Heatmap")
    end
end

function _plot_channel_summary(gui_state)
    isnothing(_validate_file(gui_state, ".jld2")) && return

    try
        data = read_data(gui_state.filename[])
        if isnothing(data)
            @minimal_warning "Requested plot settings incompatible: recheck!"
            return
        end
        @async begin
            plot_channel_summary(data)
        end
    catch e
        _handle_plot_error(e, "Channel Summary")
    end
end

function _plot_joint_probability(gui_state)
    isnothing(_validate_file(gui_state, ".jld2")) && return

    try
        data = read_data(gui_state.filename[])
        if isnothing(data)
            @minimal_warning "Requested plot settings incompatible: recheck!"
            return
        end
        @async begin
            plot_joint_probability(data)
        end
    catch e
        _handle_plot_error(e, "Joint Probability")
    end
end

function _plot_erp_measurement_gui(gui_state)
    isnothing(_validate_file(gui_state, ".jld2")) && return

    try
        data = read_data(gui_state.filename[])
        if isnothing(data)
            @minimal_warning "Requested plot settings incompatible: recheck!"
            return
        end
        @async begin
            plot_erp_measurement_gui(data)
        end
    catch e
        _handle_plot_error(e, "ERP Measurement GUI")
    end
end

function _plot_layout(gui_state)
    if gui_state.layout_file[] == ""
        @minimal_warning "Requested plot settings incompatible: recheck!"
        return
    end

    try
        layout = read_layout(gui_state.layout_file[])
        @async begin
            plot_layout_2d(layout)
        end
    catch e
        _handle_plot_error(e, "Layout")
    end
end
