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

    # find where triggers change (onset detection)
    trigger_changes = diff(vcat(0, trigger_data))
    onset_indices = trigger_changes .> 0

    # insert onsets
    cleaned = zeros(eltype(trigger_data), length(trigger_data))
    cleaned[onset_indices] = trigger_data[onset_indices]

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
        if trigger_info |> !isnothing
            insert!(empty_cols, 2, String[])  # trigger_info column
        end
        column_symbols = [:trigger; trigger_info |> !isnothing ? [:trigger_info] : []; Symbol.(column_names)]
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
    if trigger_info |> !isnothing
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




"""
    trigger_count(dat::BiosemiDataFormat.BiosemiData)::TriggerInfo

Count trigger occurrences in BioSemi data comparing raw vs cleaned counts.

This function provides a comprehensive view of trigger counts in BioSemi data by
comparing raw trigger counts with cleaned counts (onset-only). This comparison
helps verify trigger cleaning operations and identify sustained vs onset triggers.

The returned `TriggerInfo` object displays as a formatted table and can be accessed
like a DataFrame (e.g., `info.trigger`, `info.raw_count`, `info.cleaned_count`).

# Arguments
- `dat::BiosemiDataFormat.BiosemiData`: The BioSemiData object containing EEG data.

# Returns
A `TriggerInfo` object containing trigger counts for both raw and cleaned data. Display the object to see a formatted table.

# Examples
```julia
# Get trigger counts
trigger_counts = trigger_count(dat)

# Display the table
trigger_counts

# Access DataFrame properties
trigger_counts.data.trigger
trigger_counts.data.raw_count
trigger_counts.data.cleaned_count
```
"""
function trigger_count(dat::BiosemiDataFormat.BiosemiData)::TriggerInfo
    # Get cleaned trigger data (onset detection only)
    cleaned_triggers = _clean_triggers(dat.triggers.raw)
    return _trigger_count_impl([dat.triggers.raw, cleaned_triggers], ["raw_count", "cleaned_count"])
end

# =============================================================================
# TRIGGER SEQUENCE SEARCH FUNCTIONS
# =============================================================================

"""
    search_sequence(array, sequences; ignore_values::Vector{Int} = [0], sort_indices::Bool = true)

Return the matched positions of any of the specified sequences within an array.

Each match is returned as a `Vector{Int}` containing the sample indices of all
triggers in the matched sequence. For example, matching `[1, 5]` in trigger data
might return `[[30, 450], [800, 1220]]` — two matches, each with positions of
triggers 1 and 5.

# Arguments
- `array`: Array to search within
- `sequences::Vector{Vector{Union{Int,Symbol,UnitRange{Int}}}}`: Vector of sequences to find
- `ignore_values::Vector{Int}`: Values to ignore when detecting onsets (default: [0])
- `sort_indices::Bool`: Whether to sort the returned matches by start position (default: true)

# Returns
- `Vector{Vector{Int}}`: Each element is a vector of matched sample positions

# Examples
```julia
# Multiple sequences (OR logic)
matches = search_sequence([1, 2, 1, 3, 1, 3, 1], [[1, 2, 1], [1, 3, 1]])

# Wildcards
matches = search_sequence([1, 2, 3, 4, 5, 1, 4, 1], [[1, :any, 3]])

# Ranges
matches = search_sequence([1, 2, 3, 4, 5], [[1:3], [5:5]])

# Mixed sequences
matches = search_sequence([1, 2, 3, 4, 5], [[1, 2:4, 5]])
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

"""
    search_sequence(array, sequence; ignore_values::Vector{Int} = [0], sort_indices::Bool = true)

Return matched positions of a single sequence within an array, handling wildcards and ranges.

# Arguments
- `array`: Array to search within
- `sequence::Vector{Union{Int,Symbol,UnitRange{Int}}}`: Sequence pattern to find
- `ignore_values::Vector{Int}`: Values to ignore when detecting onsets (default: [0])
- `sort_indices::Bool`: Whether to sort the returned matches by start position (default: true)

# Returns
- `Vector{Vector{Int}}`: Each element is a vector of matched sample positions

# Notes
- Supports wildcards (`:any`) and ranges (`1:3`)
- First element cannot be a wildcard
- Single wildcard sequences are not supported
- Ignores specified values when detecting onsets

# Examples
```julia
matches = search_sequence([1,0,2,3,1,0,2,3], [1,2]; ignore_values=[0])  # Returns [[1,3], [5,7]]
```
"""
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

"""
    search_sequence(array, value::Integer; ignore_values::Vector{Int} = [0])

Return indices of a single integer value within an array (onset detection).

# Arguments
- `array`: Array to search within
- `value::Integer`: Integer value to find
- `ignore_values::Vector{Int}`: Values to ignore when detecting onsets (default: [0])

# Returns
- `Vector{Int}`: Indices where the value starts (onsets only)

# Examples
```julia
idx = search_sequence([0, 1, 1, 0, 2, 2], 1)  # Returns [2]
idx = search_sequence([0, 1, 0, 2, 0, 1], 1)  # Returns [2, 6]
```
"""
function search_sequence(array, value::Integer; ignore_values::Vector{Int} = [0])
    indices = Int[]
    for (i, val) in enumerate(array)
        if val == value && (i == 1 || array[i-1] != value) && !(val in ignore_values)
            push!(indices, i)
        end
    end
    return indices
end

"""
    search_sequence(array, range::UnitRange; sort_indices::Bool = true)

Return indices where array values match any value within the specified range.

# Arguments
- `array`: Array to search within
- `range::UnitRange`: Integer range to match
- `sort_indices::Bool`: Whether to sort the returned indices (default: true)

# Returns
- `Vector{Int}`: Indices where values in range match

# Examples
```julia
idx = search_sequence([1, 2, 3, 4, 5, 6], 1:3)  # Returns indices of 1,2,3
```
"""
function search_sequence(array, range::UnitRange; sort_indices::Bool = true)
    return search_sequence(array, [range]; sort_indices = sort_indices)
end

"""
    search_sequence(array, ranges::Vector{<:UnitRange}; sort_indices::Bool = true)

Return indices where array values match any of the values within specified ranges.

# Arguments
- `array`: Array to search within
- `ranges::Vector{UnitRange}`: Vector of integer ranges to match
- `sort_indices::Bool`: Whether to sort the returned indices (default: true)

# Returns
- `Vector{Int}`: Indices where values in ranges match

# Examples
```julia
idx = search_sequence([1, 2, 3, 4, 5, 6], [1:3, 5:6])  # Returns indices of 1,2,3,5,6
```
"""
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
_matches_expected(actual, expected) = error("Unsupported sequence type: $expected")

