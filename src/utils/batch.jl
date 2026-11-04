"""
    BatchResult

Result of processing a single file in a batch operation (`success`, `filename`, `message`).
"""
struct BatchResult
    success::Bool
    filename::String
    message::String
end

function _extract_participant_id(filename::String)
    name_without_ext, _ = splitext(filename)
    cu = codeunits(name_without_ext)
    i = length(cu)
    while i >= 1 && !(0x30 <= cu[i] <= 0x39)
        i -= 1
    end
    if i >= 1
        end_idx = i
        while i >= 1 && (0x30 <= cu[i] <= 0x39)
            i -= 1
        end
        return parse(Int, SubString(name_without_ext, i + 1, end_idx))
    end
    participant = Int(hash(filename) % 1_000_000_000)
    @info "  No numeric participant ID found in filename, using hash: $participant"
    return participant
end

function _find_batch_files(pattern::String, dir::String, @nospecialize(participant_selection = participants()))
    part_func = _to_participant_func(participant_selection)
    files = filter(f -> endswith(f, ".jld2") && contains(f, pattern), readdir(dir))
    ids = [_extract_participant_id(f) for f in files]
    return files[part_func(ids)]
end

_find_batch_files(pattern::String, dir::String, ids::Int) = _find_batch_files(pattern, dir, x -> x .== ids)
_find_batch_files(pattern::String, dir::String, ids::Vector{Int}) = _find_batch_files(pattern, dir, x -> [id in ids for id in x])
_find_batch_files(pattern::String, dir::String, ::Nothing) = _find_batch_files(pattern, dir, x -> trues(length(x)))

const EEGFUN_FORMAT_VERSION = 1

"""
    read_data(filepath::String)

Read `EegFunData` (or a vector of `EegFunData`) from a JLD2 file, returning `nothing` if the file
is missing, empty, or contains no `EegFunData` objects.
"""
function read_data(filepath::String)::Union{EegFunData,Vector{<:EegFunData},Nothing}
    if !isfile(filepath)
        @minimal_warning "File does not exist: $filepath"
        return nothing
    end
    jldopen(filepath, "r"; typemap = _LEGACY_JLD2_TYPEMAP) do file
        if haskey(file, "format_version")
            v = file["format_version"]
            if v > EEGFUN_FORMAT_VERSION
                @minimal_warning "File $filepath was written with EegFun format_version $v, but this version only supports up to $EEGFUN_FORMAT_VERSION."
            end
        end
        haskey(file, "data") && return _read_data(file["data"])
        keys_list = collect(keys(file))
        isempty(keys_list) && return nothing
        data = length(keys_list) == 1 ? file[keys_list[1]] : Dict(k => file[k] for k in keys_list)
        return _read_data(data)
    end
end

function _read_data(data::Vector{Any})::Union{Vector{<:EegFunData},Nothing}
    (isempty(data) || !all(x -> x isa EegFunData, data)) && return nothing
    T = mapreduce(typeof, typejoin, data)
    return T[x for x in data]
end

function _read_data(data::Dict)::Union{EegFunData,Vector{<:EegFunData},Nothing}
    values_found = EegFunData[]
    for value in values(data)
        result = _read_data(value)
        isnothing(result) && continue
        if result isa Vector{<:EegFunData}
            append!(values_found, result)
        elseif result isa EegFunData
            push!(values_found, result)
        end
    end
    return isempty(values_found) ? nothing : _single_or_vector(values_found)
end

_read_data(data::Union{EegFunData,Vector{<:EegFunData}}) = data
_read_data(::Any)::Nothing = nothing

_single_or_vector(v::Vector) = length(v) == 1 ? v[1] : v

function _condition_select(data, condition_selection::Function = conditions())
    isempty(data) && return data
    return data[condition_selection(1:length(data))]
end
_condition_select(data, condition_selection::Vector{Int}) = isempty(data) ? data : data[condition_selection]
_condition_select(data, condition_selection::Int) = _condition_select(data, [condition_selection])
_condition_select(data, ::Nothing) = _condition_select(data, x -> trues(length(x)))

function _read_all_data_core(::Type{T}, files::Vector{String}, input_dir::String) where {T}
    all_data = T[]
    sorted_files = sort(files, by = _natural_sort_key)
    n_files = length(sorted_files)
    for (i, file) in enumerate(sorted_files)
        @info "Loading: $file ($i/$n_files)"
        file_data = read_data(joinpath(input_dir, file))
        isnothing(file_data) && continue
        if file_data isa Vector{<:T}
            append!(all_data, file_data)
        elseif file_data isa T
            push!(all_data, file_data)
        else
            @minimal_warning "  Skipping $file: data is not of type $T"
        end
    end
    return all_data
end

function _split_path_pattern(path_pattern::String)
    d = dirname(path_pattern)
    return isempty(d) ? pwd() : d, basename(path_pattern)
end

"""
    read_all_data(path_pattern::String, participant_selection = participants())
    read_all_data(::Type{T}, path_pattern::String, participant_selection = participants())

Find and load EEG data from `.jld2` files in `dirname(path_pattern)` whose names contain
`basename(path_pattern)`, optionally filtered by `participant_selection`.
"""
function read_all_data(path_pattern::String, @nospecialize(participant_selection = participants()))
    input_dir, pattern = _split_path_pattern(path_pattern)
    files = _find_batch_files(pattern, input_dir, participant_selection)
    data = _read_all_data_core(EegData, files, input_dir)
    isempty(data) && return data
    T = typeof(data[1])
    return all(x -> x isa T, data) ? Vector{T}(data) : data
end

function read_all_data(::Type{T}, path_pattern::String, @nospecialize(participant_selection = participants())) where {T}
    input_dir, pattern = _split_path_pattern(path_pattern)
    files = _find_batch_files(pattern, input_dir, participant_selection)
    return _read_all_data_core(T, files, input_dir)
end

"""
    group_by_condition(items::AbstractVector{T}) -> OrderedDict{Int, Vector{T}}

Group items by their `.condition` field into a condition-sorted `OrderedDict`.
"""
function group_by_condition(items::AbstractVector{T}) where {T}
    grouped = Dict{Int,Vector{T}}()
    for item in items
        push!(get!(() -> T[], grouped, item.condition), item)
    end
    return OrderedDict(sort!(collect(grouped), by = first))
end

_validate_input_dir(dir::String) = isdir(dir) ? nothing : "Input directory does not exist: $dir"

function _validate_channel_groups(groups::Vector{Vector{Symbol}})
    isempty(groups) && return "Channel groups cannot be empty"
    for (i, group) in enumerate(groups)
        isempty(group) && return "Channel group $i is empty"
        if length(group) < 2
            @minimal_warning "Channel group $i has only $(length(group)) channel(s): $group. Consider using more channels for meaningful averaging."
        end
    end
    return nothing
end

function _validate_condition_groups(groups::AbstractVector)
    isempty(groups) && return "Condition groups cannot be empty"
    all_conditions = Set{Int}()
    for (group_idx, group) in enumerate(groups)
        unique_group = unique(group)
        if length(unique_group) != length(group)
            duplicates = [x for x in group if count(==(x), group) > 1]
            @minimal_warning "Group $group_idx contains duplicate conditions: $duplicates. Only unique conditions will be used."
            @info "  Group $group_idx: $(group) → $(unique_group) (removed duplicates)"
        end
        if groups isa Vector{Vector{Int}} || groups isa Vector{<:Vector}
            groups[group_idx] = unique_group
        end
        overlap = [c for c in unique_group if c in all_conditions]
        if !isempty(overlap)
            @minimal_warning "Condition(s) $overlap appear in multiple groups. This may not be intended."
        end
        union!(all_conditions, unique_group)
    end
    return nothing
end

function _validate_condition_pairs(pairs::AbstractVector)
    isempty(pairs) && return "Condition pairs cannot be empty"
    for (i, pair) in enumerate(pairs)
        if pair[1] == pair[2]
            @minimal_warning "Condition pair $i: ($(pair[1]), $(pair[2])) has identical conditions. Difference will be zero."
        end
    end
    return nothing
end

function _exec_batch_file(process_fn::Function, file::String, input_dir::String, output_dir::String, filename_modifier::Function)
    input_path = joinpath(input_dir, file)
    output_path = joinpath(output_dir, filename_modifier(file))
    try
        return process_fn(input_path, output_path)
    catch e
        msg = sprint(showerror, e)
        @minimal_warning "Error processing $file: $msg"
        return BatchResult(false, file, "Exception: $msg")
    end
end

function _run_batch_operation(
    process_fn::Function,
    files::Vector{String},
    input_dir::String,
    output_dir::String;
    operation_name::String = "Processing",
    filename_modifier::Function = identity,
    parallel::Bool = false,
)
    n_files = length(files)
    results = Vector{BatchResult}(undef, n_files)
    if parallel && n_files > 1 && Threads.nthreads() > 1
        Threads.@threads for i = 1:n_files
            results[i] = _exec_batch_file(process_fn, files[i], input_dir, output_dir, filename_modifier)
        end
    else
        for (i, file) in enumerate(files)
            @info "$operation_name: $file ($i/$n_files)"
            res = _exec_batch_file(process_fn, file, input_dir, output_dir, filename_modifier)
            results[i] = res
            if res.success
                @info "  ✓ $(res.message)"
            else
                @minimal_warning "  ✗ $(res.message)"
            end
        end
    end
    return results
end

"""
    batch_process(process_fn::Function, file_pattern::String, input_dir::String, output_dir::String,
                  participant_selection, operation_name::String;
                  filename_modifier::Function = identity, parallel::Bool = false)

Apply `process_fn(input_path, output_path)` across all `.jld2` files in `input_dir` matching
`file_pattern` and `participant_selection`.
"""
function batch_process(
    process_fn::Function,
    file_pattern::String,
    input_dir::String,
    output_dir::String,
    @nospecialize(participant_selection),
    operation_name::String;
    filename_modifier::Function = identity,
    parallel::Bool = false,
)
    files = _find_batch_files(file_pattern, input_dir, participant_selection)
    if isempty(files)
        @minimal_warning "No JLD2 files found matching pattern '$file_pattern' in $input_dir"
        return nothing
    end
    @info "Found $(length(files)) JLD2 files matching pattern '$file_pattern'"
    results = _run_batch_operation(
        process_fn,
        files,
        input_dir,
        output_dir;
        operation_name = operation_name,
        filename_modifier = filename_modifier,
        parallel = parallel,
    )
    return _log_batch_summary(results, output_dir)
end

function _log_batch_summary(results::Vector{BatchResult}, output_dir::String)
    n_success = count(r -> r.success, results)
    n_error = length(results) - n_success
    @info ""
    @info "Batch operation complete! Processed $n_success files successfully, $n_error errors"
    @info "Output saved to: $output_dir"
    return (success = n_success, errors = n_error)
end

function _cleanup_logging(log_file::String, output_dir::Union{String,Nothing})
    close_global_logging()
    if !isnothing(output_dir) && isfile(log_file)
        dest = joinpath(output_dir, basename(log_file))
        abspath(log_file) != abspath(dest) && mv(log_file, dest, force = true)
    end
end

mutable struct BatchContext
    output_dir::Union{String,Nothing}
end

function _with_batch_logging(f::Function, description::AbstractString, log_file::AbstractString, input_dir::AbstractString, output_dir)
    setup_global_logging(log_file)
    ctx = BatchContext(output_dir)
    try
        @info "Batch $description started at $(now())"
        error_msg = _validate_input_dir(input_dir)
        !isnothing(error_msg) && @minimal_error(error_msg)
        return f(ctx)
    finally
        _cleanup_logging(String(log_file), ctx.output_dir)
    end
end

function _prepare_output_dir!(ctx::BatchContext, default_dir::AbstractString)
    ctx.output_dir = something(ctx.output_dir, String(default_dir))
    mkpath(ctx.output_dir)
    return ctx.output_dir
end
