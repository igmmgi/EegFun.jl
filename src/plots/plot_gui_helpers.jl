# Utility functions for plot_gui helper functions

"""Check if a filename implies a pattern match (non-.jld2 extension)."""
_is_pattern(fname) = lowercase(splitext(fname)[2]) != ".jld2"

"""Get the input directory from GUI state, falling back to pwd() if empty."""
_gui_dir(gui_state) = isempty(gui_state.directory[]) ? pwd() : gui_state.directory[]
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
    if e isa EegFunError # Already logged by @minimal_error inside the plot function — don't repeat it
        return
    end
    # Unexpected exception: one clean line, no stacktrace
    @minimal_warning "Unexpected error creating $plot_type plot: $(sprint(showerror, e))"
end

"""Convert a `(lo, hi)` GUI observable pair to `Union{Nothing, Tuple{Real,Real}}`."""
_gui_lim(tup) = (!isnothing(tup[1]) && !isnothing(tup[2])) ? (tup[1], tup[2]) : nothing

"""
    _simple_jld2_plot(gui_state, plot_fn, plot_type)

Boilerplate for bridges that just load a .jld2 file and call a single-argument plot function.
"""
function _simple_jld2_plot(gui_state, plot_fn, plot_type::String)
    isnothing(_validate_file(gui_state, ".jld2")) && return
    try
        data = read_data(gui_state.filename[])
        if isnothing(data)
            @minimal_warning "Requested plot settings incompatible: recheck!"
            return
        end
        @async begin
            try
                plot_fn(data)
            catch e
                _handle_plot_error(e, plot_type)
            end
        end
    catch e
        _handle_plot_error(e, plot_type)
    end
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
        if dat isa Union{EpochData,Vector{<:EpochData}}
            dat = subset(dat; epoch_selection = epochs(epoch_list))
        end
    end

    return dat
end

"""
    _resolve_batch_files_helper(gui_state) -> Union{Nothing, Tuple{Vector{String}, Any, String}}

Resolves batch file selection, parsed participant list, and input directory.
Returns `(files, part_sel, input_dir)` or `nothing` if validation fails.
Warns if requested IDs are not found in the directory.
"""
function _resolve_batch_files_helper(gui_state)
    fname = gui_state.filename[]
    input_dir = _gui_dir(gui_state)

    part_list = _parse_gui_int_field(gui_state.participant[])
    if isnothing(part_list)
        @minimal_warning "Requested plot settings incompatible: recheck!"
        return nothing
    end

    part_sel = isempty(part_list) ? participants() : participants(part_list)

    if !isempty(part_list)
        all_files   = _find_batch_files(fname, input_dir)
        avail_ids   = sort(unique(_extract_participant_id.(all_files)))
        missing_ids = filter(id -> id ∉ avail_ids, part_list)
        if !isempty(missing_ids)
            @minimal_warning "Participant ID(s) $missing_ids not found in '$fname'. Available IDs: $avail_ids"
        end
    end

    files = _find_batch_files(fname, input_dir, part_sel)
    return isempty(files) ? nothing : (files, part_sel, input_dir)
end

# Multi-number validation helper (accepts space/comma-separated integers)
function setup_multi_int_callback(input, bc_obs::Observable, field_name::String, state_obs::Observable)
    last_validated = Ref{String}("__uninit__")

    function validate(value)
        stripped = strip(value)
        if stripped == last_validated[]
            return
        end
        last_validated[] = stripped

        if isempty(stripped)
            bc_obs[] = :white
            state_obs[] = ""
        else
            parsed = _parse_gui_int_field(stripped)
            if isnothing(parsed)
                bc_obs[] = :lightcoral
                state_obs[] = "__INVALID__"
                @minimal_warning "Invalid $field_name value: \"$stripped\" (must be integers or ranges, e.g. 1 2 3 or 1:50)"
            else
                bc_obs[] = :white
                state_obs[] = stripped
            end
        end
    end

    on(input.stored_string) do value
        validate(value)
    end

    on(input.focused) do is_focused
        if !is_focused
            validate(input.displayed_string[])
        end
    end
end

# Helper function for numeric input validation with visual feedback
function setup_numeric_callback(input, bc_obs::Observable, field_name::String, update_fn::Function; range_check = nothing)
    last_validated = Ref{String}("__uninit__")

    function validate(value)
        stripped = strip(value)
        if stripped == last_validated[]
            return
        end
        last_validated[] = stripped

        if isempty(stripped)
            bc_obs[] = :white
            update_fn(nothing)
        else
            parsed = tryparse(Float64, stripped)
            if isnothing(parsed)
                bc_obs[] = :lightcoral
                @minimal_warning "Invalid $field_name value: \"$stripped\" (must be a number)"
            else
                bc_obs[] = :white
                update_fn(parsed)
            end
        end
        # Run range check after updating value
        !isnothing(range_check) && range_check()
    end

    # Validate on Enter
    on(input.stored_string) do value
        validate(value)
    end

    # Validate on defocus (clicking away)
    on(input.focused) do is_focused
        if !is_focused
            validate(input.displayed_string[])
        end
    end
end

# Range validation: if both values are set and min >= max, both boxes turn red
function check_range(limit_obs, min_bc, max_bc, label)
    lo, hi = limit_obs[]
    if !isnothing(lo) && !isnothing(hi) && lo > hi
        min_bc[] = :lightcoral
        max_bc[] = :lightcoral
        @minimal_warning "Invalid $label range: min ($lo) must be less than or equal to max ($hi)"
    end
end

"""GUI bridge: load data and call `plot_erp_image`."""
function _plot_erp_image(gui_state)
    isnothing(_validate_file(gui_state)) && return

    fname = gui_state.filename[]
    is_pattern = _is_pattern(fname)
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
                    input_dir = _gui_dir(gui_state),
                    layout = layout_sym,
                    channel_selection = selected_channels,
                    xlim = gui_state.xlim[],
                    yreversed = gui_state.invert_y[],
                    extra_kwargs...,
                )
            else
                # .jld2 mode — load data and apply GUI filters
                dat = read_data(fname)
                isnothing(dat) && return
                dat = _apply_gui_filters(dat, gui_state)
                isnothing(dat) && return
                plot_erp_image(
                    dat;
                    layout = layout_sym,
                    channel_selection = selected_channels,
                    xlim = gui_state.xlim[],
                    yreversed = gui_state.invert_y[],
                    extra_kwargs...,
                )
            end
        catch e
            _handle_plot_error(e, "ERP Image")
        end
    end
end


"""GUI bridge: load time-frequency data and call `plot_tf`."""
function _plot_time_frequency(gui_state)
    isnothing(_validate_file(gui_state, ".jld2")) && return

    fname = gui_state.filename[]
    input_dir = _gui_dir(gui_state)
    is_pattern = _is_pattern(fname)

    # Shared plot kwargs
    selected_channels = _gui_selected_channels(gui_state)
    xlim_val = _gui_lim(gui_state.xlim[])
    zlim_val = _gui_lim(gui_state.zlim[])
    baseline_val = _gui_lim((gui_state.baseline_start[], gui_state.baseline_end[]))
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
        res = _resolve_batch_files_helper(gui_state)
        isnothing(res) && return
        files, _, input_dir = res

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
    @async begin
        try
            data = read_data(fname)
            if isnothing(data)
                @minimal_warning "Requested plot settings incompatible: recheck!"
                return
            end
            tf_vec = _apply_tf_conditions(data)
            isempty(tf_vec) && return
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
        catch e
            _handle_plot_error(e, "Time-Frequency")
        end
    end
end


"""GUI bridge: load ICA and call `plot_topography` for component maps."""
function _plot_ica(gui_state)
    fname      = gui_state.filename[]
    input_dir  = _gui_dir(gui_state)
    is_pattern = _is_pattern(fname)

    # Parse component selection from the component textbox (shared by both modes)
    comp_list = _parse_gui_int_field(gui_state.ica_components[])
    if isnothing(comp_list)
        @minimal_warning "Invalid component selection — use integers or ranges (e.g. 1:10, 1 3 5)"
        return
    end
    comp_sel = isempty(comp_list) ? components() : components(comp_list)

    if is_pattern
        # ── Batch / pattern mode ──────────────────────────────────────────────
        res = _resolve_batch_files_helper(gui_state)
        isnothing(res) && return
        files, part_sel, input_dir = res

        @async begin
            for file in sort(files, by = _natural_sort_key)
                @info "Plotting ICA: $file"
                try
                    data = read_data(joinpath(input_dir, file))
                    isnothing(data) && continue
                    if !(data isa InfoIca)
                        @minimal_warning "Skipping $file — not an InfoIca file"
                        continue
                    end
                    _set_window_title("ICA Components — " * basename(file))
                    plot_topography(data; component_selection = comp_sel)
                catch e
                    _handle_plot_error(e, "ICA Components ($file)")
                end
            end
        end

    else
        # ── Single .jld2 mode ─────────────────────────────────────────────────
        isnothing(_validate_file(gui_state, ".jld2")) && return
        try
            data = read_data(fname)
            if isnothing(data)
                @minimal_warning "Requested plot settings incompatible: recheck!"
                return
            end
            if !(data isa InfoIca)
                @minimal_warning "ICA Components requires an InfoIca file. Select a file produced by run_ica()."
                return
            end
            @async begin
                try
                    _set_window_title("ICA Components — " * basename(fname))
                    plot_topography(data; component_selection = comp_sel)
                catch e
                    _handle_plot_error(e, "ICA Components")
                end
            end
        catch e
            _handle_plot_error(e, "ICA Components")
        end
    end
end


"""
    _find_ica_source_file(ica_path, search_dirs) -> String

Locate the raw EEG source file that a given InfoIca was computed from.
Searches `search_dirs` (in order) for a file matching `ica.filename` basename
with common EEG extensions. Returns the first match, or "" if not found.
"""
function _find_ica_source_file(source_basename::String, search_dirs::Vector{String})
    extensions = [".bdf", ".jld2", ".set", ".vhdr"]
    for dir in search_dirs, ext in extensions
        candidate = joinpath(dir, source_basename * ext)
        isfile(candidate) && return candidate
    end
    return ""
end

"""
    _plot_one_ica_activation(ica_path, comp_sel, layout, gui_dir)

Load one InfoIca .jld2, locate its source EEG file, preprocess (for raw files),
and call `plot_ica_component_activation`. Used by both single and batch modes.
"""
function _plot_one_ica_activation(ica_path::String, comp_sel, layout, gui_dir::String)
    try
        ica = read_data(ica_path)
        isnothing(ica) && return
        if !(ica isa InfoIca)
            @minimal_warning "Skipping $(basename(ica_path)) — not an InfoIca file"
            return
        end

        source_basename = splitext(basename(ica.filename))[1]
        search_dirs =
            unique(filter(!isempty, [dirname(ica_path), dirname(dirname(ica_path)), dirname(dirname(dirname(ica_path))), gui_dir, pwd()]))
        source_path = _find_ica_source_file(source_basename, search_dirs)

        if isempty(source_path)
            extensions = [".bdf", ".jld2", ".set", ".vhdr"]
            searched = join(["  $d/$source_basename{$(join(extensions, ","))}" for d in search_dirs], "\n")
            @minimal_warning "Cannot locate source data for Component Activation.\n" *
                             "Searched for \"$source_basename\" with extensions $(join(extensions, ", ")) in:\n" *
                             searched
            return
        end

        src_ext = lowercase(splitext(source_path)[2])
        dat = if src_ext == ".jld2"
            read_data(source_path)
        else
            raw = read_raw_data(source_path)
            d = isnothing(layout) ? create_eegfun_data(raw) : create_eegfun_data(raw, layout)
            @info "Applying preprocessing: highpass 0.1 Hz + average rereference"
            highpass_filter!(d, 0.1)
            rereference!(d, :avg)
            d
        end

        if isnothing(dat) || !(dat isa ContinuousData)
            @minimal_warning "Source file \"$source_path\" did not load as ContinuousData."
            return
        end
        _set_window_title("Component Activation — " * basename(ica_path))
        plot_ica_component_activation(dat, ica; component_selection = comp_sel)
    catch e
        _handle_plot_error(e, "Component Activation ($(basename(ica_path)))")
    end
end

"""GUI bridge: load ICA + source data and call `plot_ica_component_activation`."""
function _plot_ica_activation(gui_state)
    fname      = gui_state.filename[]
    input_dir  = _gui_dir(gui_state)
    is_pattern = _is_pattern(fname)
    layout     = gui_state.layout_object[]

    comp_list = _parse_gui_int_field(gui_state.ica_components[])
    if isnothing(comp_list)
        @minimal_warning "Invalid component selection — use integers or ranges (e.g. 1:10, 1 3 5)"
        return
    end
    comp_sel = isempty(comp_list) ? components() : components(comp_list)

    if is_pattern
        res = _resolve_batch_files_helper(gui_state)
        isnothing(res) && return
        files, part_sel, input_dir = res

        @async begin
            for file in sort(files, by = _natural_sort_key)
                @info "Plotting Component Activation: $file"
                _plot_one_ica_activation(joinpath(input_dir, file), comp_sel, layout, input_dir)
            end
        end
    else
        isnothing(_validate_file(gui_state, ".jld2")) && return
        @async _plot_one_ica_activation(fname, comp_sel, layout, input_dir)
    end
end


"""GUI bridge: load data and call `channel_summary` + `plot_channel_summary`."""
function _plot_channel_summary(gui_state)
    fname = gui_state.filename[]
    if fname == ""
        @minimal_warning "No file selected — use 'Select File' to choose one."
        return
    end

    file_ext = lowercase(splitext(fname)[2])
    layout   = gui_state.layout_object[]

    @async begin
        try
            dat = if file_ext ∈ (".bdf", ".set", ".vhdr")
                raw = read_raw_data(fname)
                d = isnothing(layout) ? create_eegfun_data(raw) : create_eegfun_data(raw, layout)
                @info "Applying preprocessing: highpass 0.1 Hz + average rereference"
                highpass_filter!(d, 0.1)
                rereference!(d, :avg)
                d
            else
                read_data(fname)
            end
            isnothing(dat) && return
            df = channel_summary(dat)
            # channel_summary returns DataFrame or Vector{DataFrame} (epoch data)
            non_plot = (:channel, :epoch, :file, :condition)
            df_ref = df isa AbstractVector ? first(df) : df
            plot_cols = Symbol[c for c in propertynames(df_ref) if c ∉ non_plot]
            isempty(plot_cols) && return
            # If epoch column present, average across epochs and show 95% CI
            avg_kw = :epoch ∈ propertynames(df_ref) ? (average_over = :epoch,) : (;)
            plot_channel_summary(df, plot_cols; avg_kw...)
        catch e
            _handle_plot_error(e, "Channel Summary")
        end
    end
end

_plot_erp_measurement_gui(gui_state) = _simple_jld2_plot(gui_state, plot_erp_measurement_gui, "ERP Measurement GUI")

"""GUI bridge: load layout and call `plot_layout_2d`."""
function _plot_layout(gui_state)
    if gui_state.layout_file[] == ""
        @minimal_warning "Requested plot settings incompatible: recheck!"
        return
    end
    try
        layout = read_layout(gui_state.layout_file[])
        @async begin
            try
                plot_layout_2d(layout)
            catch e
                _handle_plot_error(e, "Layout")
            end
        end
    catch e
        _handle_plot_error(e, "Layout")
    end
end
