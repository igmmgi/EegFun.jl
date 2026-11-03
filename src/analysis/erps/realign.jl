const TriggerSpec = Union{Integer,AbstractVector{<:Integer}}

_as_triggers(t::Integer) = (t,)
_as_triggers(t::AbstractVector{<:Integer}) = t

# === CALCULATING TRIGGER INTERVALS ===

function _epoch_trigger_interval(trigs::AbstractVector, times::AbstractVector, starts, ends)::Tuple{Float64,Bool}
    start_idx = findfirst(in(starts), trigs)
    isnothing(start_idx) && return (NaN, true)
    rel_idx = findfirst(in(ends), @view trigs[start_idx:end])
    isnothing(rel_idx) && return (NaN, true)
    end_idx = start_idx + rel_idx - 1
    return (Float64(times[end_idx] - times[start_idx]), false)
end

"""
    calculate_trigger_interval!(dat::EpochData, start_triggers, end_triggers; column_name = :interval) -> EpochData
    calculate_trigger_interval!(dat_vec::AbstractVector{<:EpochData}, start_triggers, end_triggers; column_name = :interval) -> AbstractVector{<:EpochData}

Calculate the time interval (`end - start`, in seconds) between the first occurrence of `start_triggers`
and the subsequent first occurrence of `end_triggers` within each epoch, storing the result as a
constant column `column_name` (`NaN` if either trigger is missing). Accepts integer or vector trigger specs.
"""
function calculate_trigger_interval!(
    dat::EpochData,
    start_triggers::TriggerSpec,
    end_triggers::TriggerSpec;
    column_name::Union{Symbol,AbstractString} = :interval,
)
    col_sym = Symbol(column_name)
    starts = _as_triggers(start_triggers)
    ends = _as_triggers(end_triggers)
    @info "Calculating trigger intervals for EpochData (condition: $(dat.condition_name))"

    n_missing = 0
    for epoch in dat
        val, is_missing = _epoch_trigger_interval(epoch[!, :trigger], epoch[!, :time], starts, ends)
        n_missing += is_missing
        if hasproperty(epoch, col_sym) && epoch[!, col_sym] isa Vector{Float64}
            fill!(epoch[!, col_sym]::Vector{Float64}, val)
        else
            epoch[!, col_sym] = fill(val, nrow(epoch))
        end
    end

    if n_missing > 0
        @minimal_warning "Missing start or end trigger pair in $n_missing / $(n_epochs(dat)) epochs (set to NaN)."
    end
    return dat
end

function calculate_trigger_interval!(
    dat_vec::AbstractVector{<:EpochData},
    start_triggers::TriggerSpec,
    end_triggers::TriggerSpec;
    column_name::Union{Symbol,AbstractString} = :interval,
)
    for dat in dat_vec
        calculate_trigger_interval!(dat, start_triggers, end_triggers; column_name = column_name)
    end
    return dat_vec
end

"""
    calculate_trigger_interval(dat::EpochData, start_triggers, end_triggers; column_name = :interval) -> EpochData
    calculate_trigger_interval(dat_vec::AbstractVector{<:EpochData}, start_triggers, end_triggers; column_name = :interval) -> Vector{EpochData}
    calculate_trigger_interval(file_pattern::String, start_triggers, end_triggers; input_dir = pwd(), participant_selection = participants(), output_dir = nothing, column_name = :interval)

Non-mutating and batch file versions of [`calculate_trigger_interval!`](@ref).
"""
function calculate_trigger_interval(
    dat::EpochData,
    start_triggers::TriggerSpec,
    end_triggers::TriggerSpec;
    column_name::Union{Symbol,AbstractString} = :interval,
)::EpochData
    return calculate_trigger_interval!(copy(dat), start_triggers, end_triggers; column_name = column_name)
end

function calculate_trigger_interval(
    dat_vec::AbstractVector{<:EpochData},
    start_triggers::TriggerSpec,
    end_triggers::TriggerSpec;
    column_name::Union{Symbol,AbstractString} = :interval,
)::Vector{EpochData}
    return EpochData[calculate_trigger_interval(dat, start_triggers, end_triggers; column_name = column_name) for dat in dat_vec]
end

# === REALIGNMENT HELPERS ===

function _shifted_epoch_time(trigs::AbstractVector, times::AbstractVector, triggers)::Union{Vector{Float64},Nothing}
    resp_idx = findfirst(in(triggers), trigs)
    isnothing(resp_idx) && return nothing
    rt = Float64(times[resp_idx])
    if times isa Vector{Float64}
        return times .- rt
    end
    return Float64[Float64(t) - rt for t in times]
end

function _scan_realigned_epochs(dat::EpochData, triggers; warn_dropped::Bool = true)
    kept_indices = Int[]
    shifted_times = Vector{Float64}[]
    sizehint!(kept_indices, length(dat.data))
    sizehint!(shifted_times, length(dat.data))
    for (i, epoch) in enumerate(dat.data)
        st = _shifted_epoch_time(epoch[!, :trigger], epoch[!, :time], triggers)
        if !isnothing(st)
            push!(kept_indices, i)
            push!(shifted_times, st)
        elseif warn_dropped
            @minimal_warning "Epoch $i dropped: missing realignment trigger."
        end
    end
    return kept_indices, shifted_times
end

function _common_window_from_times(shifted_times::Vector{Vector{Float64}})::Tuple{Float64,Float64}
    common_start = maximum(t[1] for t in shifted_times)
    common_end = minimum(t[end] for t in shifted_times)
    if common_start >= common_end
        @minimal_error(
            "No common time interval found after realignment. " *
            "Common start: $(common_start) s, Common end: $(common_end) s. " *
            "The original epochs may not be long enough to accommodate all realignment times."
        )
    end
    return common_start, common_end
end

function _crop_start_indices(shifted_times::Vector{Vector{Float64}}, window::Tuple{Float64,Float64})
    n_ep = length(shifted_times)
    start_idxs = Vector{Int}(undef, n_ep)
    min_samples = typemax(Int)
    for i = 1:n_ep
        t = shifted_times[i]
        s_idx = find_closest_time_index(t, window[1])
        e_idx = find_closest_time_index(t, window[2])
        start_idxs[i] = s_idx
        min_samples = min(min_samples, e_idx - s_idx + 1)
    end
    return start_idxs, min_samples
end

function _build_cropped_epochs(
    orig_dfs::Vector{DataFrame},
    kept_indices::Vector{Int},
    shifted_times::Vector{Vector{Float64}},
    window::Tuple{Float64,Float64},
    sample_rate::Real,
    global_n_samples::Union{Int,Nothing} = nothing,
)::Vector{DataFrame}
    isempty(kept_indices) && return DataFrame[]
    start_idxs, local_n_samples = _crop_start_indices(shifted_times, window)
    n_samples = isnothing(global_n_samples) ? local_n_samples : global_n_samples

    fs = Float64(sample_rate)
    k_start = round(Int, shifted_times[1][start_idxs[1]] * fs)
    uniform_time = Float64[(k_start + j - 1) / fs for j = 1:n_samples]

    out = Vector{DataFrame}(undef, length(kept_indices))
    for k in eachindex(kept_indices)
        orig_df = orig_dfs[kept_indices[k]]
        s_idx = start_idxs[k]
        e_idx = s_idx + n_samples - 1
        if e_idx > nrow(orig_df)
            @minimal_error(
                "Epoch $k does not have enough samples for common interval. " *
                "Need samples up to index $e_idx, but epoch only has $(nrow(orig_df)) samples. " *
                "Window: [$(window[1]), $(window[2])]"
            )
        end
        cropped = orig_df[s_idx:e_idx, :]
        cropped[!, :time] = copy(uniform_time)
        out[k] = cropped
    end
    return out
end

function _realign_single_cropped(dat::EpochData, realignment_triggers::TriggerSpec)::Vector{DataFrame}
    triggers = _as_triggers(realignment_triggers)
    @info "Realigning epoched data to triggers: $realignment_triggers"
    kept_indices, shifted_times = _scan_realigned_epochs(dat, triggers)
    if isempty(kept_indices)
        @minimal_warning "All epochs were dropped due to missing realignment triggers!"
        return DataFrame[]
    end
    common_interval = _common_window_from_times(shifted_times)
    cropped = _build_cropped_epochs(dat.data, kept_indices, shifted_times, common_interval, dat.sample_rate)
    @info "Common time interval: $(common_interval[1]) s to $(common_interval[2]) s (strict: $(nrow(cropped[1])) samples)"
    return cropped
end

function _realign_multi_cropped(
    dat_vec::AbstractVector{<:EpochData},
    realignment_triggers::TriggerSpec;
    global_interval::Interval = nothing,
    global_n_samples::Union{Int,Nothing} = nothing,
)::Vector{Vector{DataFrame}}
    triggers = _as_triggers(realignment_triggers)
    @info "Realigning $(length(dat_vec)) conditions in Vector{EpochData}..."

    norm_interval = isnothing(global_interval) ? nothing : _normalize_interval(global_interval)
    scans = [_scan_realigned_epochs(dat, triggers) for dat in dat_vec]

    if isnothing(norm_interval)
        common_starts = Float64[]
        common_ends = Float64[]
        for (_, shifted_times) in scans
            isempty(shifted_times) && continue
            s, e = _common_window_from_times(shifted_times)
            push!(common_starts, s)
            push!(common_ends, e)
        end
        isempty(common_starts) && return [DataFrame[] for _ in dat_vec]
        norm_interval = (maximum(common_starts), minimum(common_ends))
    end

    if isnothing(global_n_samples)
        min_len = typemax(Int)
        for (_, shifted_times) in scans
            isempty(shifted_times) && continue
            _, cond_min = _crop_start_indices(shifted_times, norm_interval)
            min_len = min(min_len, cond_min)
        end
        min_len == typemax(Int) && return [DataFrame[] for _ in dat_vec]
        global_n_samples = min_len
    end

    @info "Global common time interval: $(norm_interval[1]) s to $(norm_interval[2]) s (strict: $global_n_samples samples)"
    return Vector{DataFrame}[
        _build_cropped_epochs(dat_vec[i].data, scans[i][1], scans[i][2], norm_interval, dat_vec[i].sample_rate, global_n_samples) for
        i in eachindex(dat_vec)
    ]
end

# === CORE REALIGNMENT API ===

"""
    realign!(dat::EpochData, realignment_triggers) -> EpochData
    realign!(dat_vec::AbstractVector{<:EpochData}, realignment_triggers; global_interval = nothing, global_n_samples = nothing) -> AbstractVector{<:EpochData}

Realign epoched data in-place so the first matching trigger in `realignment_triggers` becomes `t = 0`,
dropping epochs without a matching trigger and cropping all remaining epochs to a shared common time window.
Accepts an integer trigger or vector of triggers.
"""
function realign!(dat::EpochData, realignment_triggers::TriggerSpec)
    dat.data = _realign_single_cropped(dat, realignment_triggers)
    return dat
end

function realign!(
    dat_vec::AbstractVector{<:EpochData},
    realignment_triggers::TriggerSpec;
    global_interval::Interval = nothing,
    global_n_samples::Union{Int,Nothing} = nothing,
)
    cropped_sets =
        _realign_multi_cropped(dat_vec, realignment_triggers; global_interval = global_interval, global_n_samples = global_n_samples)
    for i in eachindex(dat_vec)
        dat_vec[i].data = cropped_sets[i]
    end
    return dat_vec
end

"""
    realign(dat::EpochData, realignment_triggers) -> EpochData
    realign(dat_vec::AbstractVector{<:EpochData}, realignment_triggers; global_interval = nothing, global_n_samples = nothing) -> Vector{EpochData}
    realign(file_pattern::String, realignment_triggers; input_dir = pwd(), participant_selection = participants(), output_dir = nothing)

Non-mutating and batch file versions of [`realign!`](@ref).
"""
function realign(dat::EpochData, realignment_triggers::TriggerSpec)::EpochData
    cropped = _realign_single_cropped(dat, realignment_triggers)
    return EpochData(dat.file, dat.condition, dat.condition_name, cropped, copy(dat.layout), dat.sample_rate, copy(dat.analysis_info))
end

function realign(
    dat_vec::AbstractVector{<:EpochData},
    realignment_triggers::TriggerSpec;
    global_interval::Interval = nothing,
    global_n_samples::Union{Int,Nothing} = nothing,
)::Vector{EpochData}
    cropped_sets =
        _realign_multi_cropped(dat_vec, realignment_triggers; global_interval = global_interval, global_n_samples = global_n_samples)
    return EpochData[
        EpochData(
            dat_vec[i].file,
            dat_vec[i].condition,
            dat_vec[i].condition_name,
            cropped_sets[i],
            copy(dat_vec[i].layout),
            dat_vec[i].sample_rate,
            copy(dat_vec[i].analysis_info),
        ) for i in eachindex(dat_vec)
    ]
end

# === BATCH PROCESSING ===

_default_realign_output_dir(input_dir::String, pattern::String) = joinpath(input_dir, "realigned_$(clean_pattern(pattern))")
_default_calculate_trigger_interval_output_dir(input_dir::String) = joinpath(input_dir, "epochs_calculate_trigger_interval")

function _process_realign_file(
    filepath::String,
    output_path::String,
    realignment_triggers::TriggerSpec,
    global_interval::Tuple{Float64,Float64},
    global_n_samples::Int,
)
    filename = basename(filepath)
    epochs_data = read_data(filepath)
    isnothing(epochs_data) && return BatchResult(false, filename, "No data found in file")
    (epochs_data isa EpochData || epochs_data isa Vector{<:EpochData}) || return BatchResult(false, filename, "Data is not EpochData type")

    triggers = _as_triggers(realignment_triggers)
    if epochs_data isa EpochData
        kept_indices, shifted_times = _scan_realigned_epochs(epochs_data, triggers)
        isempty(kept_indices) && return BatchResult(true, filename, "All epochs dropped (missing trigger).")
        epochs_data.data =
            _build_cropped_epochs(epochs_data.data, kept_indices, shifted_times, global_interval, epochs_data.sample_rate, global_n_samples)
        realigned_epochs = epochs_data
    else
        realigned_epochs = realign!(epochs_data, triggers; global_interval = global_interval, global_n_samples = global_n_samples)
    end

    jldsave(output_path; data = realigned_epochs, format_version = EEGFUN_FORMAT_VERSION)

    if realigned_epochs isa Vector
        n_eps = sum(n_epochs(cond) for cond in realigned_epochs)
        n_s = isempty(realigned_epochs) || isempty(realigned_epochs[1].data) ? 0 : nrow(realigned_epochs[1].data[1])
    else
        n_eps = n_epochs(realigned_epochs)
        n_s = isempty(realigned_epochs.data) ? 0 : nrow(realigned_epochs.data[1])
    end
    return BatchResult(true, filename, "Realigned $n_eps total epochs to $n_s samples each")
end

function realign(
    file_pattern::String,
    realignment_triggers::TriggerSpec;
    input_dir::String = pwd(),
    participant_selection = participants(),
    output_dir::Union{String,Nothing} = nothing,
)
    _with_batch_logging("realignment", "realign.log", input_dir, output_dir) do ctx
        @log_call "realign"
        output_dir = _prepare_output_dir!(ctx, _default_realign_output_dir(input_dir, file_pattern))

        files = _find_batch_files(file_pattern, input_dir, participant_selection)
        if isempty(files)
            @minimal_warning "No JLD2 files found matching pattern '$file_pattern' in $input_dir"
            return nothing
        end
        @info "Found $(length(files)) JLD2 files matching pattern '$file_pattern'"

        triggers = _as_triggers(realignment_triggers)
        @info "Pass 1/2: Global timeline scanning across $(length(files)) files..."
        grand_starts = Float64[]
        grand_ends = Float64[]
        all_shifted_sets = Vector{Vector{Float64}}[]

        for file in files
            epochs_data = read_data(joinpath(input_dir, file))
            isnothing(epochs_data) && continue
            vec_data = epochs_data isa Vector ? epochs_data : [epochs_data]
            for dat in vec_data
                isempty(dat.data) && continue
                _, shifted_times = _scan_realigned_epochs(dat, triggers; warn_dropped = false)
                isempty(shifted_times) && continue
                s, e = _common_window_from_times(shifted_times)
                push!(grand_starts, s)
                push!(grand_ends, e)
                push!(all_shifted_sets, shifted_times)
            end
        end

        isempty(grand_starts) && @minimal_error("No valid realigned data found across any file! (All trials likely missing the trigger)")

        global_interval = (maximum(grand_starts), minimum(grand_ends))
        global_n_samples = minimum(last(_crop_start_indices(st, global_interval)) for st in all_shifted_sets)
        @info "Grand Global Interval: $(global_interval[1]) to $(global_interval[2]) (Strict: $global_n_samples samples)"

        @info "Pass 2/2: Realigning to disk..."
        process_fn =
            (input_path, output_path) -> _process_realign_file(input_path, output_path, triggers, global_interval, global_n_samples)
        results = _run_batch_operation(process_fn, files, input_dir, output_dir; operation_name = "Realigning epochs")
        _log_batch_summary(results, output_dir)
    end
end

function _process_calculate_trigger_interval_file(
    filepath::String,
    output_path::String,
    start_triggers::TriggerSpec,
    end_triggers::TriggerSpec,
    column_name::Symbol,
)
    filename = basename(filepath)
    epochs_data = read_data(filepath)
    isnothing(epochs_data) && return BatchResult(false, filename, "No data found in file")
    (epochs_data isa EpochData || epochs_data isa AbstractVector{<:EpochData}) ||
        return BatchResult(false, filename, "Data is not EpochData type")

    calculate_trigger_interval!(epochs_data, start_triggers, end_triggers; column_name = column_name)
    jldsave(output_path; data = epochs_data, format_version = EEGFUN_FORMAT_VERSION)

    n_eps = epochs_data isa Vector ? sum(n_epochs(cond) for cond in epochs_data) : n_epochs(epochs_data)
    return BatchResult(true, filename, "Appended interval column (:$column_name) to $n_eps total epochs")
end

function calculate_trigger_interval(
    file_pattern::String,
    start_triggers::TriggerSpec,
    end_triggers::TriggerSpec;
    input_dir::String = pwd(),
    participant_selection = participants(),
    output_dir::Union{String,Nothing} = nothing,
    column_name::Union{Symbol,AbstractString} = :interval,
)
    col_sym = Symbol(column_name)
    _with_batch_logging("calculate_trigger_interval", "calculate_trigger_interval.log", input_dir, output_dir) do ctx
        @log_call "calculate_trigger_interval"
        output_dir = _prepare_output_dir!(ctx, _default_calculate_trigger_interval_output_dir(input_dir))
        process_fn =
            (input_path, output_path) ->
                _process_calculate_trigger_interval_file(input_path, output_path, start_triggers, end_triggers, col_sym)
        batch_process(process_fn, file_pattern, input_dir, output_dir, participant_selection, "Calculating trigger intervals")
    end
end
