"""
Generic utility functions for batch processing operations.
These are reusable across different batch scripts.
"""

#=============================================================================
    DATA STRUCTURES
=============================================================================#

"""Immutable struct to hold batch processing results."""
struct BatchResult
    success::Bool
    filename::String
    message::String
end

"""Configuration for batch operations."""
struct BatchConfig
    file_pattern::String
    input_dir::String
    output_dir::String
    participants::Union{Int,Vector{Int},Nothing}
    conditions::Union{Int,Vector{Int},Nothing}
end

#=============================================================================
    PURE FUNCTIONS (no side effects)
=============================================================================#

"""Extract participant ID from filename, returns Int."""
function _extract_participant_id(filename::String)
    # Only search the filename part (without extension)
    name_without_ext, _ = splitext(filename)
   
    # extract last numeric sequence from filename (e.g., "exp1a1" -> 1, "Pract6" -> 6)
    numeric_matches = collect(eachmatch(r"\d+", name_without_ext))
    if !isempty(numeric_matches)
        return parse(Int, numeric_matches[end].match)
    end

    # use hash of filename as fallback
    participant = hash(filename) % 10000
    @info "  No numeric participant ID found in filename, using hash: $participant"
    return participant
end


"""
    find_batch_files(pattern::String, dir::String; participants=nothing)
    find_batch_files(pattern::String, dir::String; participants::Function)

Find JLD2 files matching pattern and optional participant filter.

Returns vector of filenames (not full paths).
"""
# Method for Function predicate
function _find_batch_files(pattern::String, dir::String, participant_selection::Function = participants())

    # Filter by pattern and extension
    all_files = readdir(dir)
    files = Base.filter(all_files) do f
        endswith(f, ".jld2") && contains(f, pattern)
    end

    # extract ids from filenames  and apply predicate mask
    ids = [_extract_participant_id(f) for f in files]
    mask = participant_selection(ids)
    return files[mask]
end

# create predicate like input participants()
_find_batch_files(pattern::String, dir::String, participants::Int) =
    _find_batch_files(pattern, dir, x -> x .== participants)
_find_batch_files(pattern::String, dir::String, participants::Vector{Int}) =
    _find_batch_files(pattern, dir, x -> [id in participants for id in x])
_find_batch_files(pattern::String, dir::String, participants::Nothing) =
    _find_batch_files(pattern, dir, x -> fill(true, length(x)))


"""
    load_data(filepath::String)

Load data from JLD2 file, returning the data directly, a Dict of all variables, or `nothing`.

- If file has 1 variable: returns the value directly
- If file has multiple variables: returns a Dict with all key-value pairs
- If file is empty: returns `nothing`
"""
function load_data(filepath::String)::Union{EegData, Vector{<:EegData}, InfoIca, Vector{InfoIca}, Nothing}
    jldopen(filepath, "r") do file
        keys_list = collect(keys(file))
        isempty(keys_list) && return nothing
        
        data = length(keys_list) == 1 ? file[keys_list[1]] : Dict(k => file[k] for k in keys_list)
        return _load_data(data)
    end
end

# Convert Vector{Any} to typed vector if all elements are the same type
function _load_data(data::Vector{Any})::Union{Vector{<:EegData}, Vector{InfoIca}, Nothing}
    isempty(data) && return nothing
    T = typeof(data[1])
    (T <: EegData || T <: InfoIca) && all(x -> typeof(x) == T, data) || return nothing
    return Vector{T}(data)
end

# Extract EegData or InfoIca from Dict
function _load_data(data::Dict)::Union{EegData, Vector{<:EegData}, InfoIca, Vector{InfoIca}, Nothing}
    eeg_values = EegData[]
    ica_values = InfoIca[]
    for value in values(data)
        result = _load_data(value)
        isnothing(result) && continue
        _add_to_collection!(result, eeg_values, ica_values)
    end
    # return single element if only one
    !isempty(eeg_values) && return _single_or_vector(eeg_values)
    !isempty(ica_values) && return _single_or_vector(ica_values)
    return nothing
end

# load_data is not a generic function for everything; 
# we just use it for data that is saved from eegfun
_load_data(data::Union{EegData, Vector{<:EegData}, InfoIca, Vector{InfoIca}}) = data
_load_data(::Any)::Nothing = nothing

# Helper to add loaded data to appropriate collection
_add_to_collection!(data::EegData, eeg_values, _) = push!(eeg_values, data)
_add_to_collection!(data::Vector{<:EegData}, eeg_values, _) = append!(eeg_values, data)
_add_to_collection!(data::InfoIca, _, ica_values) = push!(ica_values, data)
_add_to_collection!(data::Vector{InfoIca}, _, ica_values) = append!(ica_values, data)

# Return single element or vector
_single_or_vector(v::Vector) = length(v) == 1 ? v[1] : v

"""
    _condition_select(data, condition_selection)

Filter data by condition selection predicate.

Returns filtered data.
"""
function _condition_select(data, condition_selection::Function = conditions())
    isempty(data) && return data
    condition_indices = 1:length(data)
    mask = condition_selection(condition_indices)
    return data[mask]
end

"""
    _condition_select(data, condition_selection)

Filter data by condition selection predicate.

Returns filtered data.
"""
function _condition_select(data, condition_selection::Vector{Int})
    isempty(data) && return data
    isnothing(condition_selection) && return data
    # Ensure we always return a vector, even for single Int selection
    return data[condition_selection]
end

_condition_select(data, condition_selection::Int) = _condition_select(data, [condition_selection])
_condition_select(data, condition_selection::Nothing) = _condition_select(data, x -> fill(true, length(x)))


"""
    load_all_data(files::Vector{String}, input_dir::String)
    load_all_data(::Type{T}, files::Vector{String}, input_dir::String) where {T}
    load_all_data(file_pattern::String; input_dir, participant_selection)
    load_all_data(::Type{T}, file_pattern::String; input_dir, participant_selection) where {T}

Load all data from files into a flat vector.

# Arguments (Vector version)
- `files::Vector{String}`: Filenames to load (not full paths)
- `input_dir::String`: Directory containing the files

# Arguments (Pattern version)
- `file_pattern::String`: Pattern to match files (e.g., "erps", "epochs")
- `input_dir::String`: Directory containing the files (default: current directory)
- `participant_selection::Function`: Predicate to filter participants (default: all)

# Type Parameter
- `T` (optional): Filter to only load data of this type (e.g., `ErpData`, `EpochData`)

# Returns
- `Vector`: All data from all files combined (typed if T specified)
"""
function load_all_data(files::Vector{String}, input_dir::String)
    all_data = EegData[]
    for (i, file) in enumerate(sort(files, by = natural_sort_key))
        input_path = joinpath(input_dir, file)
        @info "Loading: $file ($i/$(length(files)))"
        file_data = load_data(input_path)
        isnothing(file_data) && continue
        if file_data isa Vector{<:EegData}
            append!(all_data, file_data)
        elseif file_data isa EegData
            push!(all_data, file_data)
        end
    end
    return all_data
end

function load_all_data(::Type{T}, files::Vector{String}, input_dir::String) where {T}
    all_data = T[]
    for (i, file) in enumerate(sort(files, by = natural_sort_key))
        input_path = joinpath(input_dir, file)
        @info "Loading: $file ($i/$(length(files)))"
        file_data = load_data(input_path)
        isnothing(file_data) && continue
        file_data isa Vector{<:T} || continue
        append!(all_data, file_data)
    end
    return all_data
end

# Pattern-based versions
function load_all_data(
    file_pattern::String;
    input_dir::String = pwd(),
    participant_selection::Function = participants(),
)
    files = _find_batch_files(file_pattern, input_dir, participant_selection)
    return load_all_data(files, input_dir)
end

function load_all_data(
    ::Type{T},
    file_pattern::String;
    input_dir::String = pwd(),
    participant_selection::Function = participants(),
) where {T}
    files = _find_batch_files(file_pattern, input_dir, participant_selection)
    return load_all_data(T, files, input_dir)
end


"""
    group_by_condition(erps::Vector{<:ErpData})

Group ERPs by their condition number.

# Arguments
- `erps::Vector{<:ErpData}`: ERPs to group

# Returns
- `OrderedDict{Int, Vector{ErpData}}`: ERPs grouped by condition number (sorted)
"""
function group_by_condition(erps::Vector{<:ErpData})
    grouped = OrderedDict{Int,Vector{ErpData}}()
    for erp in erps
        cond_num = erp.condition
        push!(get!(grouped, cond_num, ErpData[]), erp)
    end
    # Sort by condition number
    return OrderedDict(sort(collect(grouped), by = first))
end


"""
    validate_input_dir(dir::String)

Validate input directory exists, returning error message or nothing.
"""
function _validate_input_dir(dir::String)
    !isdir(dir) && return "Input directory does not exist: $dir"
    return nothing
end

"""
    validate_channel_groups(groups::Vector{Vector{Symbol}})

Validate channel groups, returning error message or nothing.
Issues warnings for groups with < 2 channels.
"""
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

"""
    validate_condition_groups(groups::Vector{Vector{Int}})

Validate condition groups, returning error message or nothing.
Modifies groups in-place to remove duplicates and warns about overlaps.
"""
function _validate_condition_groups(groups::Vector{Vector{Int}})
    isempty(groups) && return "Condition groups cannot be empty"

    all_conditions = Int[]

    for (group_idx, group) in enumerate(groups)
        # Check for duplicates within the group
        if length(group) != length(unique(group))
            duplicates = group[findall(x -> count(==(x), group) > 1, group)]
            @minimal_warning "Group $group_idx contains duplicate conditions: $duplicates. Only unique conditions will be used."
        end

        # Remove duplicates and update the group
        unique_group = unique(group)
        if length(unique_group) != length(group)
            @info "  Group $group_idx: $(group) → $(unique_group) (removed duplicates)"
        end
        groups[group_idx] = unique_group

        # Check for overlaps between groups
        overlap = intersect(all_conditions, unique_group)
        if !isempty(overlap)
            @minimal_warning "Condition(s) $overlap appear in multiple groups. This may not be intended."
        end

        append!(all_conditions, unique_group)
    end

    return nothing
end

"""
    _validate_condition_pairs(pairs)

Validate condition pairs, returning error message or nothing.
Warns about pairs with identical conditions.
Accepts both `Vector{Tuple{Int, Int}}` and `Vector{Vector{Int}}`.
"""
function _validate_condition_pairs(pairs::Union{Vector{Tuple{Int,Int}},Vector{Vector{Int}}})
    isempty(pairs) && return "Condition pairs cannot be empty"

    for (i, pair) in enumerate(pairs)
        if pair[1] == pair[2]
            @minimal_warning "Condition pair $i: ($(pair[1]), $(pair[2])) has identical conditions. Difference will be zero."
        end
    end

    return nothing
end


"""
    run_batch_operation(process_fn::Function, files::Vector{String}, 
                        input_dir::String, output_dir::String;
                        operation_name::String = "Processing")

Execute batch processing with logging and error handling.

Takes a pure processing function `process_fn(input_path, output_path)` and maps 
it over files, handling errors and logging progress.

Returns vector of `BatchResult`.
"""
function _run_batch_operation(
    process_fn::Function,
    files::Vector{String},
    input_dir::String,
    output_dir::String;
    operation_name::String = "Processing",
)
    n_files = length(files)
    results = Vector{BatchResult}(undef, n_files)

    for (i, file) in enumerate(files)
        @info "$operation_name: $file ($i/$n_files)"

        input_path = joinpath(input_dir, file)
        output_path = joinpath(output_dir, file)

        result = try
            process_fn(input_path, output_path)
        catch e
            @error "Error processing $file" exception=(e, catch_backtrace())
            BatchResult(false, file, "Exception: $(sprint(showerror, e))")
        end

        results[i] = result

        # Log individual result
        if result.success
            @info "  ✓ $(result.message)"
        else
            @minimal_warning "  ✗ $(result.message)"
        end
    end

    return results
end

"""
    log_batch_summary(results::Vector{BatchResult}, output_dir::String)

Log summary statistics from batch results.

Returns named tuple `(success=n, errors=n)`.
"""
function _log_batch_summary(results::Vector{BatchResult}, output_dir::String)
    n_success = count(r -> r.success, results)
    n_error = length(results) - n_success

    @info ""
    @info "Batch operation complete! Processed $n_success files successfully, $n_error errors"
    @info "Output saved to: $output_dir"

    return (success = n_success, errors = n_error)
end

"""
Cleanup global logging and move log file to output directory.
If output_dir is nothing, just closes logging without moving the file.
"""
function _cleanup_logging(log_file::String, output_dir::Union{String,Nothing})
    close_global_logging()

    if !isnothing(output_dir)
        log_dest = joinpath(output_dir, log_file)
        log_file != log_dest && mv(log_file, log_dest, force = true)
    end
end
