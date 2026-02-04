"""
Generic utility functions for batch processing operations.
These are reusable across different batch scripts.
"""

"""struct to hold batch processing results."""
struct BatchResult
    success::Bool
    filename::String
    message::String
end

"""configuration for batch operations."""
struct BatchConfig
    file_pattern::String
    input_dir::String
    output_dir::String
    participants::Union{Int,Vector{Int},Nothing}
    conditions::Union{Int,Vector{Int},Nothing}
end

"""Extract participant ID (number) from filename, returns Int."""
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
    _find_batch_files(pattern::String, dir::String, participant_selection::Function = participants())

Find and filter JLD2 files in a directory based on a pattern and participant IDs.

# Arguments
- `pattern::String`: Filename pattern to match.
- `dir::String`: Directory to search in.
- `participant_selection::Function`: Predicate function for filtering participant IDs (default: `participants()`).

# Returns
- `Vector{String}`: List of matching filenames.
"""
# Method for Function predicate
function _find_batch_files(pattern::String, dir::String, participant_selection::Function = participants())

    # Filter by pattern and extension
    all_files = readdir(dir)
    files = filter(all_files) do f
        endswith(f, ".jld2") && contains(f, pattern)
    end

    # extract ids from filenames  and apply predicate mask
    ids = [_extract_participant_id(f) for f in files]
    mask = participant_selection(ids)
    return files[mask]
end

# create predicate like input participants()
_find_batch_files(pattern::String, dir::String, ids::Int) = _find_batch_files(pattern, dir, x -> x .== ids)
_find_batch_files(pattern::String, dir::String, ids::Vector{Int}) = _find_batch_files(pattern, dir, x -> [id in ids for id in x])
_find_batch_files(pattern::String, dir::String, ids::Nothing) = _find_batch_files(pattern, dir, x -> fill(true, length(x)))


"""
    read_data(filepath::String)

Load data from JLD2 file, returning the data directly, a Dict of all variables, or `nothing`.

- If file has 1 variable: returns the value directly
- If file has multiple variables: returns a Dict with all key-value pairs
- If file is empty: returns `nothing`
"""
function read_data(filepath::String)::Union{EegFunData,Vector{<:EegFunData},Nothing}
    jldopen(filepath, "r") do file
        keys_list = collect(keys(file))
        isempty(keys_list) && return nothing

        data = length(keys_list) == 1 ? file[keys_list[1]] : Dict(k => file[k] for k in keys_list)
        return _read_data(data)
    end
end

# Convert Vector{Any} to typed vector if all elements are EegFun data
function _read_data(data::Vector{Any})::Union{Vector{<:EegFunData},Nothing}
    isempty(data) && return nothing
    # Ensure every element is some kind of EegFun data
    !all(x -> x isa EegFunData, data) && return nothing

    # Use a comprehension to let Julia infer the best (narrowest) type
    # e.g., Vector{ErpData} if all same, Vector{EegFunData} if mixed
    return [x for x in data]
end

# Extract EegFunData from Dict
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

# NB. load_data is not a generic function for everything; we just use it for data that is saved from EegFun
_read_data(data::Union{EegFunData,Vector{<:EegFunData}}) = data
_read_data(::Any)::Nothing = nothing

# helper to load and combine results
_add_to_collection!(data::EegFunData, values, _) = push!(values, data)
_add_to_collection!(data::Vector{<:EegFunData}, values, _) = append!(values, data)

# return single element or vector
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

function _condition_select(data, condition_selection::Vector{Int})
    isempty(data) && return data
    isnothing(condition_selection) && return data
    # Ensure we always return a vector, even for single Int selection
    return data[condition_selection]
end

_condition_select(data, condition_selection::Int) = _condition_select(data, [condition_selection])
_condition_select(data, condition_selection::Nothing) = _condition_select(data, x -> fill(true, length(x)))


# Core internal loading logic
function _read_all_data_core(::Type{T}, files::Vector{String}, input_dir::String) where {T}
    all_data = T[]
    for (i, file) in enumerate(sort(files, by = natural_sort_key))
        input_path = joinpath(input_dir, file)
        @info "Loading: $file ($i/$(length(files)))"
        file_data = read_data(input_path)
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

"""
    read_all_data(files::Vector{String}, input_dir::String)
    read_all_data(::Type{T}, files::Vector{String}, input_dir::String)

Load EEG data from a list of files into a flat vector.
"""
read_all_data(files::Vector{String}, input_dir::String) = _read_all_data_core(EegData, files, input_dir)
read_all_data(::Type{T}, files::Vector{String}, input_dir::String) where {T} = _read_all_data_core(T, files, input_dir)

"""
    read_all_data(pattern::String, input_dir::String, participant_selection = participants())
    read_all_data(::Type{T}, pattern::String, input_dir::String, participant_selection = participants())

Find and load EEG data matching a pattern and participant selection.
"""
function read_all_data(pattern::String, input_dir::String = pwd(), participant_selection::Function = participants())
    files = _find_batch_files(pattern, input_dir, participant_selection)
    return read_all_data(files, input_dir)
end

function read_all_data(::Type{T}, pattern::String, input_dir::String = pwd(), participant_selection::Function = participants()) where {T}
    files = _find_batch_files(pattern, input_dir, participant_selection)
    return read_all_data(T, files, input_dir)
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
    group_by_condition(epochs::Vector{<:EpochData})

Group Epochs by their condition number.

# Arguments
- `epochs::Vector{<:EpochData}`: Epochs to group

# Returns
- `OrderedDict{Int, Vector{EpochData}}`: Epochs grouped by condition number (sorted)
"""
function group_by_condition(epochs::Vector{<:EpochData})
    grouped = OrderedDict{Int,Vector{EpochData}}()
    for epoch in epochs
        cond_num = epoch.condition
        push!(get!(grouped, cond_num, EpochData[]), epoch)
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

    all_conditions = Set{Int}()

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
        overlap = [c for c in unique_group if c in all_conditions]
        if !isempty(overlap)
            @minimal_warning "Condition(s) $overlap appear in multiple groups. This may not be intended."
        end

        union!(all_conditions, unique_group)
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
            @error "Error processing $file" exception = (e, catch_backtrace())
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
