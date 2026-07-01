"""
    _clean_triggers(trigger_data::Vector{<:Integer})::Vector{<:Integer}

Clean trigger data by detecting only the onset (first occurrence) of each trigger value.

This function converts sustained trigger signals into single onset events. 
When a trigger is held for multiple samples, only the first sample is retained, 
with subsequent samples set to zero.

# Arguments
- `trigger_data::Vector{<:Integer}`: Raw trigger data vector

# Returns
- `Vector{<:Integer}`: Cleaned trigger data with only onset events

# Examples
```julia
# Input:  [0, 1, 1, 0, 0, 2, 2, 2, 0, 0]
# Output: [0, 1, 0, 0, 0, 2, 0, 0, 0, 0]
```
"""
function _clean_triggers(trigger_data::Vector{<:Integer})::Vector{<:Integer}
    n = length(trigger_data)
    cleaned = zeros(eltype(trigger_data), n)
    if n == 0
        return cleaned
    end

    if trigger_data[1] > 0
        cleaned[1] = trigger_data[1]
    end

    @inbounds for i = 2:n
        if trigger_data[i] > trigger_data[i-1]
            cleaned[i] = trigger_data[i]
        end
    end

    return cleaned
end


"""
    trigger_count(dat::ContinuousData)::TriggerInfo
    trigger_count(df::DataFrame)::TriggerInfo

Count occurrences of each trigger value in ContinuousData or DataFrame.

This function analyzes trigger data to provide a summary of how many times each
trigger value appears in the dataset. Zero values are excluded from the count.
Useful for validating experimental paradigms and checking trigger timing.

The returned `TriggerInfo` object displays as a formatted table and can be accessed
like a DataFrame (e.g., `info.data.trigger`, `info.data.count`).

# Arguments
- `dat::ContinuousData`: The ContinuousData object containing EEG data.
- `df::DataFrame`: A DataFrame with a `:trigger` column.

# Returns
A `TriggerInfo` object containing trigger counts. Display the object to see a formatted table.

# Examples
```julia
# Get trigger counts from ContinuousData
trigger_counts = trigger_count(dat)

# Get trigger counts from DataFrame
trigger_counts = trigger_count(df)

# Display the table
trigger_counts

# Access DataFrame properties
trigger_counts.data.trigger
trigger_counts.data.count
```
"""
function trigger_count(df::DataFrame)::TriggerInfo

    if !hasproperty(df, :trigger)
        @minimal_error "DataFrame must have a trigger column"
    end

    # Check if trigger_info column exists and pass it along
    trigger_info = hasproperty(df, :trigger_info) ? df.trigger_info : nothing
    return _trigger_count_impl([df.trigger], ["count"], trigger_info = trigger_info)
end
trigger_count(dat::ContinuousData)::TriggerInfo = trigger_count(dat.data)


"""
    _trigger_count_impl(trigger_datasets, column_names; trigger_info=nothing)

Trigger counting function with optional trigger info support.

# Arguments
- `trigger_datasets::Vector{<:Vector{<:Integer}}`: Vector of trigger datasets to analyze
- `column_names::Vector{String}`: Names for the count columns in the output DataFrame
- `trigger_info::Union{Nothing,Vector{String}}`: Optional trigger info strings

# Returns
- `TriggerInfo`: TriggerInfo object containing DataFrame with 'trigger' column and count columns (plus trigger_info if provided)
"""
function _trigger_count_impl(
    trigger_datasets::Vector{<:Vector{<:Integer}},
    column_names::Vector{String};
    trigger_info::Union{Nothing,Vector{String}} = nothing,
)
    # Get unique non-zero trigger values from all datasets
    all_triggers = vcat(trigger_datasets...)
    non_zero_triggers = sort(unique(filter(x -> x != 0, all_triggers)))

    if isempty(non_zero_triggers)
        # Return empty DataFrame with correct structure
        empty_cols = [Int[]]  # trigger column
        append!(empty_cols, [Int[] for _ in column_names])  # count columns
        if !isnothing(trigger_info)
            insert!(empty_cols, 2, String[])  # trigger_info column
        end
        column_symbols = [:trigger; !isnothing(trigger_info) ? [:trigger_info] : []; Symbol.(column_names)]
        result_df = DataFrame([col => data for (col, data) in zip(column_symbols, empty_cols)]...)
        return TriggerInfo(result_df)
    end

    # Count occurrences of each trigger value across all datasets
    result_data = Vector{Vector{Int}}(undef, length(column_names) + 1)
    result_data[1] = non_zero_triggers

    for (i, dataset) in enumerate(trigger_datasets)
        trigger_counts = Dict{Int,Int}()
        for val in dataset
            if val != 0
                trigger_counts[val] = get(trigger_counts, val, 0) + 1
            end
        end
        result_data[i+1] = [get(trigger_counts, trigger, 0) for trigger in non_zero_triggers]
    end

    # Create basic DataFrame
    column_symbols = [:trigger; Symbol.(column_names)]
    result_df = DataFrame([col => data for (col, data) in zip(column_symbols, result_data)]...)

    # Add trigger_info column if provided
    if !isnothing(trigger_info)
        trigger_info_map = Dict{Int,String}()
        for (trigger, info) in zip(trigger_datasets[1], trigger_info)
            if trigger != 0 && !haskey(trigger_info_map, trigger)
                trigger_info_map[trigger] = info
            end
        end

        # Insert trigger_info column after trigger column
        trigger_info_col = [get(trigger_info_map, trigger, "") for trigger in non_zero_triggers]
        result_df.trigger_info = trigger_info_col

        # Reorder columns: trigger, trigger_info, count columns
        count_cols = [col for col in names(result_df) if col != :trigger && col != :trigger_info]
        result_df = result_df[:, Cols(:trigger, :trigger_info, count_cols...)]
    end

    return TriggerInfo(result_df)
end


# Prettier trigger display info
function Base.show(io::IO, info::TriggerInfo)

    if isempty(info.data)
        println(io, "No non-zero triggers found in the data.")
        return
    end

    # Determine title based on column names
    if "raw_count" in names(info.data) && "cleaned_count" in names(info.data)
        title = "Trigger Count Summary (Raw vs Cleaned)"
    else
        title = "Trigger Count Summary"
    end

    # Calculate alignment based on actual DataFrame columns
    n_cols = length(names(info.data))
    alignment = if hasproperty(info.data, :trigger_info)
        # trigger (r), trigger_info (l), then rest (r)
        [:r, :l, fill(:r, n_cols - 2)...]
    else
        # All columns right-aligned
        fill(:r, n_cols)
    end

    pretty_table(io, info.data, title = title, alignment = alignment)
    println(io)
end



function trigger_count(dat::BiosemiDataFormat.BiosemiData)::TriggerInfo
    # Get cleaned trigger data (onset detection only)
    cleaned_triggers = _clean_triggers(dat.triggers.raw)
    return _trigger_count_impl([dat.triggers.raw, cleaned_triggers], ["raw_count", "cleaned_count"])
end

function trigger_count(dat::EuropeanDataFormat.EdfData)::TriggerInfo
    # EDF triggers are already annotations/onsets, no raw vs cleaned distinction needed natively
    # but we will extract them into a synthetic array just like we did for creating dataframes
    n_samples = size(dat.data, 1)
    sample_rate = dat.header.sample_rate[1]
    trigger = zeros(Int, n_samples)
    trigger_info = fill("", n_samples)

    if !isnothing(dat.triggers) && length(dat.triggers.onset) > 0
        unique_annotations = sort(unique(dat.triggers.annotation))
        value_to_trigger = Dict(ann => i for (i, ann) in enumerate(unique_annotations))
        for i = 1:length(dat.triggers.onset)
            sample_idx = round(Int, dat.triggers.onset[i] * sample_rate) + 1
            if 1 <= sample_idx <= n_samples
                trigger[sample_idx] = value_to_trigger[dat.triggers.annotation[i]]
                trigger_info[sample_idx] = dat.triggers.annotation[i]
            end
        end
    end

    return _trigger_count_impl([trigger], ["count"]; trigger_info = trigger_info)
end

function trigger_count(dat::BrainVisionDataFormat.BrainVisionData)::TriggerInfo
    if isnothing(dat.markers) || isempty(dat.markers)
        return _trigger_count_impl([Int[]], ["count"])
    end
    n_samples = size(dat.data, 1)
    trigger, trigger_info = _extract_triggers_from_markers(dat.markers, n_samples)
    return _trigger_count_impl([trigger], ["count"]; trigger_info = trigger_info)
end

function trigger_count(dat::ExtensibleDataFormat.XdfData)::TriggerInfo
    eeg_streams = [s for s in values(dat.streams) if s.header.type == "EEG"]
    if isempty(eeg_streams)
        return _trigger_count_impl([Int[]], ["count"])
    end
    eeg_stream = eeg_streams[1]
    n_samples = size(eeg_stream.time_series, 1)
    trigger = zeros(Int, n_samples)
    trigger_info = fill("", n_samples)

    marker_streams = [s for s in values(dat.streams) if s.header.type == "Markers"]
    if !isempty(marker_streams)
        marker_stream = marker_streams[1]
        for (i, t) in enumerate(marker_stream.timestamps)
            idx = find_closest_time_index(eeg_stream.timestamps, t)
            marker_val = marker_stream.time_series[i, 1]
            int_val = tryparse(Int, string(marker_val))
            trigger[idx] = int_val !== nothing ? int_val : 1
            trigger_info[idx] = string(marker_val)
        end
    end

    return _trigger_count_impl([trigger], ["count"]; trigger_info = trigger_info)
end

# =============================================================================
# TRIGGER SEQUENCE SEARCH FUNCTIONS
# =============================================================================

"""
    search_sequence(array, sequences::Vector{<:Vector}; ignore_values = [0], sort_indices = true) -> Vector{Vector{Int}}
    search_sequence(array, sequence::Vector; ignore_values = [0], sort_indices = true) -> Vector{Vector{Int}}
    search_sequence(array, value::Integer; ignore_values = [0]) -> Vector{Int}
    search_sequence(array, range::UnitRange; sort_indices = true) -> Vector{Int}
    search_sequence(array, ranges::Vector{<:UnitRange}; sort_indices = true) -> Vector{Int}

Search for trigger patterns in an array. Returns matched sample positions.

- **Multiple sequences** (`Vector{Vector}`): OR logic — matches any of the provided sequences.
- **Single sequence** (`Vector`): supports integer values, `:any` wildcard, and `UnitRange` elements.
- **Single value** (`Integer`): onset detection — returns each start position.
- **Range / Vector{UnitRange}**: returns indices of any value within the range(s).

Each match in the multi-element methods is a `Vector{Int}` of sample positions for each
trigger in the matched sequence (e.g., `[[30, 450], [800, 1220]]`).

# Examples
```julia
# Single value
search_sequence([0, 1, 1, 0, 2], 1)              # => [2]

# Exact sequence
search_sequence([1, 0, 2, 3, 1, 0, 2, 3], [1, 2]) # => [[1,3], [5,7]]

# Wildcard
search_sequence([1, 2, 3], [[1, :any, 3]])        # => [[1,2,3]]

# Range
search_sequence([1, 2, 3, 4, 5], 1:3)             # => [1, 2, 3]

# Multiple sequences (OR)
search_sequence([1, 2, 1, 3, 1], [[1, 2, 1], [1, 3, 1]])  # => [[1,2,3], [3,4,5]]
```
"""
function search_sequence(array, sequences::Vector{<:Vector}; ignore_values::Vector{Int} = [0], sort_indices::Bool = true)
    isempty(array) && return Vector{Int}[]
    isempty(sequences) && return Vector{Int}[]

    # Optimize common case of single sequence
    if length(sequences) == 1
        return search_sequence(array, sequences[1]; ignore_values = ignore_values, sort_indices = sort_indices)
    end

    all_matches = Vector{Int}[]
    for sequence in sequences
        matches = search_sequence(array, sequence; ignore_values = ignore_values, sort_indices = false)
        append!(all_matches, matches)
    end

    # Deduplicate by start position and optionally sort
    unique!(m -> m[1], all_matches)
    return sort_indices ? sort(all_matches, by = first) : all_matches
end

function search_sequence(array, sequence::Vector; ignore_values::Vector{Int} = [0], sort_indices::Bool = true)
    isempty(array) && return Vector{Int}[]
    isempty(sequence) && return Vector{Int}[]

    # Handle case where sequence is all UnitRanges (treat as range search)
    if all(x -> x isa UnitRange, sequence)
        indices = search_sequence(array, UnitRange{Int}[x for x in sequence]; sort_indices = sort_indices)
        return [[idx] for idx in indices]
    end

    # Handle single trigger case
    if length(sequence) == 1
        indices = _search_single_trigger(array, sequence[1], ignore_values)
        return [[idx] for idx in indices]
    end

    # Find starting positions for the first trigger
    idx_start_positions = _search_single_trigger(array, sequence[1], ignore_values)

    matches = Vector{Int}[]
    seq_len = length(sequence)
    max_idx = length(array) - seq_len + 1
    for idx in idx_start_positions
        if idx > max_idx
            continue
        end
        positions = _collect_sequence_positions(array, sequence, idx, ignore_values)
        if !isnothing(positions)
            push!(matches, positions)
        end
    end

    return sort_indices ? sort(matches, by = first) : matches
end

function search_sequence(array, value::Integer; ignore_values::Vector{Int} = [0])
    indices = Int[]
    for (i, val) in enumerate(array)
        if val == value && (i == 1 || array[i-1] != value) && !(val in ignore_values)
            push!(indices, i)
        end
    end
    return indices
end

function search_sequence(array, range::UnitRange; sort_indices::Bool = true)
    return search_sequence(array, [range]; sort_indices = sort_indices)
end

function search_sequence(array, ranges::Vector{<:UnitRange}; sort_indices::Bool = true)
    isempty(ranges) && return Int[]
    all_indices = Int[]
    for value in union(ranges...)
        indices = search_sequence(array, value)
        append!(all_indices, indices)
    end
    return sort_indices ? sort(all_indices) : all_indices
end

# Helper function to search for a single trigger using dispatch
_search_single_trigger(array, trigger::Integer, ignore_values::Vector{Int} = [0]) =
    search_sequence(array, trigger; ignore_values = ignore_values)
_search_single_trigger(array, trigger::UnitRange, ignore_values::Vector{Int} = [0]) = search_sequence(array, [trigger])
_search_single_trigger(array, trigger::Symbol, ignore_values::Vector{Int} = [0]) = error("Single wildcard sequences not supported")
# Set match: find onset positions of ANY value in the vector (e.g., [101, 103])
_search_single_trigger(array, trigger::Vector{<:Integer}, ignore_values::Vector{Int} = [0]) =
    sort(vcat([search_sequence(array, Int(v); ignore_values = ignore_values) for v in trigger]...))

# Helper function to collect the actual sample positions for a matched sequence
"""
    _collect_sequence_positions(array, sequence, start_idx, ignore_values) -> Union{Nothing, Vector{Int}}

Attempt to match a sequence starting at `start_idx`, collecting the actual sample
index of each matched trigger. Returns `nothing` if the sequence doesn't match,
or a `Vector{Int}` of positions if it does.
"""
function _collect_sequence_positions(array, sequence, start_idx, ignore_values)
    positions = Int[start_idx]
    current_idx = start_idx
    @inbounds for expected in sequence[2:end]
        current_idx += 1
        # Skip over ignored values
        while current_idx <= length(array) && array[current_idx] in ignore_values
            current_idx += 1
        end
        # Check if we've gone beyond the array bounds
        current_idx > length(array) && return nothing
        # Check if the current value matches the expected value
        _matches_expected(array[current_idx], expected) || return nothing
        push!(positions, current_idx)
    end
    return positions
end

# Dispatch-based pattern matching
_matches_expected(actual::Integer, expected::Integer) = actual == expected
_matches_expected(actual::Real, expected::Integer) = actual == expected
_matches_expected(actual::Integer, expected::Real) = actual == expected
_matches_expected(actual::Real, expected::Real) = actual == expected
_matches_expected(actual::Real, expected::Symbol) = expected == :any  # Wildcard matches anything
_matches_expected(actual::Real, expected::UnitRange) = actual in expected
_matches_expected(actual::Real, expected::Vector{<:Integer}) = actual in expected  # Set match: [101, 103]
_matches_expected(actual, expected) = error("Unsupported sequence type: $expected")

#=============================================================================
    EPOChed TRIGGER SEQ DISCOVERY
=============================================================================#

"""
    trigger_info(dat::Union{EpochData, Vector{EpochData}})

Analyze an epoched dataset and return the unique sequences of triggers found within each condition.

For a single `EpochData`, it returns a `Vector{NamedTuple}` with `sequence` and `t0_trigger` fields, 
representing the unique orderings of non-zero triggers and explicitly identifying which trigger
was active seamlessly at `time == 0.0`.

# Examples
```julia
info = trigger_info(epochs_good)
# Output example: (sequence = [101, 114, 239, 201], t0 = 101)
```
"""
function trigger_info(dat::EpochData)
    unique_seqs = NamedTuple{(:sequence, :t0),Tuple{Vector{Int},Int}}[]
    for epoch in dat.data
        if hasproperty(epoch, :trigger) && hasproperty(epoch, :time)
            # Find the active trigger right at time = 0
            zero_idx = find_closest_time_index(epoch.time, 0.0)
            t0_trigger = epoch.trigger[zero_idx]

            # Extract sequence, ignoring continuous zeros
            seq = filter(x -> x != 0, epoch.trigger)

            nt = (sequence = seq, t0 = t0_trigger)
            if !(nt in unique_seqs)
                push!(unique_seqs, nt)
            end
        end
    end
    return unique_seqs
end

"""
    trigger_info(dat_vec::Vector{EpochData})

Run `trigger_info` across all conditions in a `Vector{EpochData}`.
Returns a `Dict{String, ...}` mapping each `condition_name` to its list of unique
trigger sequences. See the `EpochData` method for details on the individual entries.
"""
function trigger_info(dat_vec::Vector{EpochData})
    info = Dict{String,Vector{NamedTuple{(:sequence, :t0),Tuple{Vector{Int},Int}}}}()
    for dat in dat_vec
        info[dat.condition_name] = trigger_info(dat)
    end
    return info
end

"""
    trigger_info(file_pattern::String; input_dir::String = pwd(), participant_selection = participants())

Run `trigger_info` across all `.jld2` files matching the `file_pattern`. 
This allows you to either pass a direct filename or a generic batch string like `"epochs_good"`.
Returns a merged dictionary of all unique trigger sequences across the entire experiment.
"""
function trigger_info(file_pattern::String; input_dir::String = pwd(), participant_selection::Function = participants())

    # Check if exact file
    if isfile(file_pattern) && endswith(file_pattern, ".jld2")
        dat = read_data(file_pattern)
        return isnothing(dat) ? nothing : trigger_info(dat isa Vector ? dat : [dat])
    elseif isfile(joinpath(input_dir, file_pattern)) && endswith(file_pattern, ".jld2")
        dat = read_data(joinpath(input_dir, file_pattern))
        return isnothing(dat) ? nothing : trigger_info(dat isa Vector ? dat : [dat])
    end

    # Handle as batch pattern
    files = _find_batch_files(file_pattern, input_dir, participant_selection)
    if isempty(files)
        @minimal_warning "No JLD2 files found matching pattern '$file_pattern' in $input_dir"
        return nothing
    end

    @info "Scanning $(length(files)) files for unique trigger sequences..."

    merged_info = Dict{String,Vector{NamedTuple{(:sequence, :t0),Tuple{Vector{Int},Int}}}}()

    for file in files
        input_path = joinpath(input_dir, file)
        dat = read_data(input_path)
        if !isnothing(dat)
            vec_dat = dat isa EpochData ? [dat] : dat
            info = trigger_info(vec_dat)
            for (cond, seqs) in info
                if !haskey(merged_info, cond)
                    merged_info[cond] = seqs
                else
                    # Merge uniqueness
                    for seq in seqs
                        if !(seq in merged_info[cond])
                            push!(merged_info[cond], seq)
                        end
                    end
                end
            end
        end
    end

    return merged_info
end
