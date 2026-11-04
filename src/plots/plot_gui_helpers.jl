_is_pattern(fname) = lowercase(splitext(fname)[2]) != ".jld2"
_gui_dir(gui_state) = isempty(gui_state.directory[]) ? pwd() : gui_state.directory[]
_gui_selected_channels(gui_state) = isempty(gui_state.electrodes[]) ? channels() : channels(Symbol.(gui_state.electrodes[]))
_gui_lim(tup) = (!isnothing(tup[1]) && !isnothing(tup[2])) ? (tup[1], tup[2]) : nothing

function _validate_file(gui_state, required_ext::String = ".jld2")
    fname = gui_state.filename[]
    if fname == ""
        @minimal_warning "Requested plot settings incompatible: recheck!"
        return nothing
    end
    file_ext = lowercase(splitext(fname)[2])
    if file_ext != required_ext && !isempty(file_ext)
        @minimal_warning "Requested plot settings incompatible: recheck!"
        return nothing
    end
    return fname
end

function _handle_plot_error(e, plot_type::String)
    if e isa EegFunError
        _file_logging_active() || _log_error(e.msg)
        return
    end
    @minimal_warning "Unexpected error creating $plot_type plot: $(sprint(showerror, e))"
end

function _guarded(f, plot_type::String)
    try
        f()
    catch e
        _handle_plot_error(e, plot_type)
    end
    return nothing
end

function _guarded_async(f, plot_type::String)
    @async _guarded(f, plot_type)
    return nothing
end

function _simple_jld2_plot(gui_state, plot_fn, plot_type::String)
    isnothing(_validate_file(gui_state, ".jld2")) && return
    fname = gui_state.filename[]
    _guarded_async(plot_type) do
        data = read_data(fname)
        isnothing(data) ? @minimal_warning("Requested plot settings incompatible: recheck!") : plot_fn(data)
    end
end

function _parse_gui_int_field(str::String)
    str == "__INVALID__" && return nothing
    s = strip(str)
    isempty(s) && return Int[]
    result = Int[]
    for token in eachsplit(s, (' ', '\t', '\n', '\r', ','); keepempty = false)
        colon_idx = findfirst(':', token)
        if !isnothing(colon_idx)
            next_idx = nextind(token, colon_idx)
            findnext(':', token, next_idx) !== nothing && return nothing
            lo = tryparse(Int, SubString(token, 1, prevind(token, colon_idx)))
            hi = tryparse(Int, SubString(token, next_idx))
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

function _apply_gui_filters(dat, gui_state)
    cond_list = _parse_gui_int_field(gui_state.condition[])
    epoch_list = _parse_gui_int_field(gui_state.epoch[])
    if isnothing(cond_list) || isnothing(epoch_list)
        @minimal_warning "Requested plot settings incompatible: recheck!"
        return nothing
    end
    !isempty(cond_list) && dat isa Vector && (dat = subset(dat; condition_selection = conditions(cond_list)))
    !isempty(epoch_list) && dat isa Union{EpochData,Vector{<:EpochData}} && (dat = subset(dat; epoch_selection = epochs(epoch_list)))
    return dat
end

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
        all_files = _find_batch_files(fname, input_dir)
        avail_ids = sort(unique(_extract_participant_id.(all_files)))
        missing_ids = filter(id -> id ∉ avail_ids, part_list)
        !isempty(missing_ids) && @minimal_warning("Participant ID(s) $missing_ids not found in '$fname'. Available IDs: $avail_ids")
    end
    files = _find_batch_files(fname, input_dir, part_sel)
    return isempty(files) ? nothing : (files, part_sel, input_dir)
end

function _bind_textbox_validation!(on_valid::Function, input)
    last_validated = Ref("__uninit__")
    function run_validate(raw)
        stripped = String(strip(raw))
        stripped == last_validated[] && return
        last_validated[] = stripped
        on_valid(stripped)
    end
    on(run_validate, input.stored_string)
    on(input.focused) do is_focused
        is_focused || run_validate(input.displayed_string[])
    end
end

function setup_multi_int_callback(input, bc_obs::Observable, field_name::String, state_obs::Observable)
    _bind_textbox_validation!(input) do stripped
        if isempty(stripped)
            bc_obs[] = :white
            state_obs[] = ""
        elseif isnothing(_parse_gui_int_field(stripped))
            bc_obs[] = :lightcoral
            state_obs[] = "__INVALID__"
            @minimal_warning "Invalid $field_name value: \"$stripped\" (must be integers or ranges, e.g. 1 2 3 or 1:50)"
        else
            bc_obs[] = :white
            state_obs[] = stripped
        end
    end
end

function setup_numeric_callback(input, bc_obs::Observable, field_name::String, update_fn::Function; range_check = nothing)
    _bind_textbox_validation!(input) do stripped
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
        !isnothing(range_check) && range_check()
    end
end

function check_range(limit_obs, min_bc, max_bc, label)
    lo, hi = limit_obs[]
    if !isnothing(lo) && !isnothing(hi) && lo > hi
        min_bc[] = :lightcoral
        max_bc[] = :lightcoral
        @minimal_warning "Invalid $label range: min ($lo) must be less than or equal to max ($hi)"
    end
end

function _plot_erp_image(gui_state)
    isnothing(_validate_file(gui_state)) && return
    fname = gui_state.filename[]
    is_pattern = _is_pattern(fname)
    selected_channels = _gui_selected_channels(gui_state)

    _guarded_async("ERP Image") do
        layout_sym = Symbol(gui_state.layout_type[])
        z_lim = _gui_lim(gui_state.zlim[])
        x_lim = _gui_lim(gui_state.xlim[])
        extra_kwargs = !isnothing(z_lim) ? (colorrange = z_lim,) : (;)

        if is_pattern
            plot_erp_image(
                fname;
                input_dir = _gui_dir(gui_state),
                layout = layout_sym,
                channel_selection = selected_channels,
                xlim = x_lim,
                yreversed = gui_state.invert_y[],
                extra_kwargs...,
            )
        else
            dat = read_data(fname)
            isnothing(dat) && return
            dat = _apply_gui_filters(dat, gui_state)
            isnothing(dat) && return
            plot_erp_image(
                dat;
                layout = layout_sym,
                channel_selection = selected_channels,
                xlim = x_lim,
                yreversed = gui_state.invert_y[],
                extra_kwargs...,
            )
        end
    end
end

function _plot_time_frequency(gui_state)
    isnothing(_validate_file(gui_state, ".jld2")) && return
    fname = gui_state.filename[]
    selected_channels = _gui_selected_channels(gui_state)
    xlim_val = _gui_lim(gui_state.xlim[])
    ylim_val = _gui_lim(gui_state.ylim[])
    zlim_val = _gui_lim(gui_state.zlim[])
    baseline_val = _gui_lim((gui_state.baseline_start[], gui_state.baseline_end[]))

    bt = gui_state.baseline_type[]
    baseline_method_sym = bt ∈ ("db", "absolute", "relative", "relchange", "percent", "zscore", "normchange") ? Symbol(bt) : :db
    lt = gui_state.layout_type[]
    layout_sym = lt ∈ ("single", "grid", "topo") ? Symbol(lt) : :single

    cond_list = _parse_gui_int_field(gui_state.condition[])
    if isnothing(cond_list)
        @minimal_warning "Requested plot settings incompatible: recheck!"
        return
    end

    function _run_tf_plot(data)
        vec = data isa AbstractVector ? data : [data]
        filtered = isempty(cond_list) ? vec : filter(tf -> tf.condition ∈ cond_list, vec)
        if isempty(filtered)
            @minimal_warning "No conditions matching $cond_list found in file"
            return
        end
        tf = length(filtered) == 1 ? only(filtered) : filtered
        layout_kw = tf isa TimeFreqData ? (layout = layout_sym,) : (;)
        plot_tf(
            tf;
            layout_kw...,
            channel_selection = selected_channels,
            xlim = xlim_val,
            ylim = ylim_val,
            colorrange = zlim_val,
            baseline_interval = baseline_val,
            baseline_method = baseline_method_sym,
        )
    end

    if _is_pattern(fname)
        res = _resolve_batch_files_helper(gui_state)
        isnothing(res) && return
        files, _, input_dir = res
        _guarded_async("Time-Frequency") do
            for file in sort(files, by = _natural_sort_key)
                @info "Plotting TF: $file"
                data = read_data(joinpath(input_dir, file))
                !isnothing(data) && _run_tf_plot(data)
            end
        end
        return
    end

    _guarded_async("Time-Frequency") do
        data = read_data(fname)
        isnothing(data) ? @minimal_warning("Requested plot settings incompatible: recheck!") : _run_tf_plot(data)
    end
end

function _parse_gui_comp_sel(gui_state)
    comp_list = _parse_gui_int_field(gui_state.ica_components[])
    if isnothing(comp_list)
        @minimal_warning "Invalid component selection — use integers or ranges (e.g. 1:10, 1 3 5)"
        return nothing
    end
    return isempty(comp_list) ? components() : components(comp_list)
end

function _plot_ica(gui_state)
    fname = gui_state.filename[]
    comp_sel = _parse_gui_comp_sel(gui_state)
    isnothing(comp_sel) && return

    if _is_pattern(fname)
        res = _resolve_batch_files_helper(gui_state)
        isnothing(res) && return
        files, _, input_dir = res
        _guarded_async("ICA Components") do
            for file in sort(files, by = _natural_sort_key)
                @info "Plotting ICA: $file"
                _guarded("ICA Components ($file)") do
                    data = read_data(joinpath(input_dir, file))
                    isnothing(data) && return
                    data isa InfoIca || (@minimal_warning "Skipping $file — not an InfoIca file"; return)
                    _set_window_title("ICA Components — " * basename(file))
                    plot_topography(data; component_selection = comp_sel)
                end
            end
        end
    else
        isnothing(_validate_file(gui_state, ".jld2")) && return
        _guarded_async("ICA Components") do
            data = read_data(fname)
            isnothing(data) && (@minimal_warning "Requested plot settings incompatible: recheck!"; return)
            data isa InfoIca || (@minimal_warning "ICA Components requires an InfoIca file. Select a file produced by run_ica()."; return)
            _set_window_title("ICA Components — " * basename(fname))
            plot_topography(data; component_selection = comp_sel)
        end
    end
end

const _RAW_EEG_EXTS = (".bdf", ".jld2", ".set", ".vhdr")

function _find_ica_source_file(source_basename::String, search_dirs::Vector{String})
    for dir in search_dirs, ext in _RAW_EEG_EXTS
        candidate = joinpath(dir, source_basename * ext)
        isfile(candidate) && return candidate
    end
    return ""
end

function _load_or_preprocess_raw(path::String, layout)
    if lowercase(splitext(path)[2]) == ".jld2"
        return read_data(path)
    end
    raw = read_raw_data(path)
    d = isnothing(layout) ? create_eegfun_data(raw) : create_eegfun_data(raw, layout)
    @info "Applying preprocessing: highpass 0.1 Hz + average rereference"
    highpass_filter!(d, 0.1)
    rereference!(d, :avg)
    return d
end

function _plot_one_ica_activation(ica_path::String, comp_sel, layout, gui_dir::String)
    _guarded("Component Activation ($(basename(ica_path)))") do
        ica = read_data(ica_path)
        isnothing(ica) && return
        ica isa InfoIca || (@minimal_warning "Skipping $(basename(ica_path)) — not an InfoIca file"; return)

        source_basename = splitext(basename(ica.filename))[1]
        search_dirs =
            unique(filter(!isempty, [dirname(ica_path), dirname(dirname(ica_path)), dirname(dirname(dirname(ica_path))), gui_dir, pwd()]))
        source_path = _find_ica_source_file(source_basename, search_dirs)

        if isempty(source_path)
            ext_str = join(_RAW_EEG_EXTS, ",")
            searched = join(["  $d/$source_basename{$ext_str}" for d in search_dirs], "\n")
            @minimal_warning "Cannot locate source data for Component Activation.\nSearched for \"$source_basename\" with extensions $(join(_RAW_EEG_EXTS, ", ")) in:\n$searched"
            return
        end

        dat = _load_or_preprocess_raw(source_path, layout)
        if isnothing(dat) || !(dat isa ContinuousData)
            @minimal_warning "Source file \"$source_path\" did not load as ContinuousData."
            return
        end
        _set_window_title("Component Activation — " * basename(ica_path))
        plot_ica_component_activation(dat, ica; component_selection = comp_sel)
    end
end

function _plot_ica_activation(gui_state)
    fname = gui_state.filename[]
    input_dir = _gui_dir(gui_state)
    layout = gui_state.layout_object[]
    comp_sel = _parse_gui_comp_sel(gui_state)
    isnothing(comp_sel) && return

    if _is_pattern(fname)
        res = _resolve_batch_files_helper(gui_state)
        isnothing(res) && return
        files, _, input_dir = res
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

function _plot_channel_summary(gui_state)
    fname = gui_state.filename[]
    if fname == ""
        @minimal_warning "No file selected — use 'Select File' to choose one."
        return
    end
    file_ext = lowercase(splitext(fname)[2])
    layout = gui_state.layout_object[]

    _guarded_async("Channel Summary") do
        dat = file_ext ∈ (".bdf", ".set", ".vhdr") ? _load_or_preprocess_raw(fname, layout) : read_data(fname)
        isnothing(dat) && return
        df = channel_summary(dat)
        df_ref = df isa AbstractVector ? first(df) : df
        plot_cols = Symbol[c for c in propertynames(df_ref) if c ∉ (:channel, :epoch, :file, :condition)]
        isempty(plot_cols) && return
        avg_kw = :epoch ∈ propertynames(df_ref) ? (average_over = :epoch,) : (;)
        plot_channel_summary(df, plot_cols; avg_kw...)
    end
end

_plot_erp_measurement_gui(gui_state) = _simple_jld2_plot(gui_state, plot_erp_measurement_gui, "ERP Measurement GUI")

function _plot_layout(gui_state)
    layout_file = gui_state.layout_file[]
    if layout_file == ""
        @minimal_warning "Requested plot settings incompatible: recheck!"
        return
    end
    _guarded_async("Layout") do
        plot_layout_2d(read_layout(layout_file))
    end
end
