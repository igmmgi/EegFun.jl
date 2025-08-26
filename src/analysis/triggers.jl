# =============================================================================
# TRIGGER CLEANING FUNCTIONS
# =============================================================================

"""
    _clean_triggers(trigger_data::Vector{<:Integer})::Vector{<:Integer}

Clean trigger data by detecting only the onset (first occurrence) of each trigger value.

This function converts sustained trigger signals into single onset events, which is
essential for proper EEG event marking. When a trigger is held for multiple samples,
only the first sample is retained, with subsequent samples set to zero.

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

# =============================================================================
# TRIGGER COUNTING FUNCTIONS
# =============================================================================

"""
    trigger_count(dat::ContinuousData; print_table::Bool = true)::DataFrame

Count occurrences of each trigger value in ContinuousData and optionally display results.

This function analyzes trigger data to provide a summary of how many times each
trigger value appears in the dataset. Zero values are excluded from the count.
Useful for validating experimental paradigms and checking trigger timing.

# Arguments
- `dat::ContinuousData`: The ContinuousData object containing EEG data.
- `print_table::Bool`: Whether to print the trigger count table (default: true).

# Returns
A DataFrame with columns `trigger` and `count` showing trigger values and their counts, excluding zero values.

# Examples
```julia
# Get trigger counts and print table
trigger_counts = trigger_count(dat)

# Get trigger counts without printing
trigger_counts = trigger_count(dat, print_table = false)
```
"""
function trigger_count(dat::ContinuousData; print_table::Bool = true)::DataFrame
    @assert hasproperty(dat.data, :triggers) "Data must have a triggers column"
    return _trigger_count_impl([dat.data.triggers], ["count"], print_table = print_table)
end

"""
    _count_triggers_core(trigger_datasets, column_names)

Core helper function for counting triggers across multiple datasets.

This function efficiently counts trigger occurrences across one or more datasets,
returning a unified DataFrame with counts for each trigger value.

# Arguments
- `trigger_datasets::Vector{<:Vector{<:Integer}}`: Vector of trigger datasets to analyze
- `column_names::Vector{String}`: Names for the count columns in the output DataFrame

# Returns
- `DataFrame`: DataFrame with 'trigger' column and one count column per dataset

# Notes
- Zero values are excluded from counting
- Uses single-pass counting with Dict for efficiency
- Returns empty DataFrame with correct structure when no triggers found
"""
function _count_triggers_core(trigger_datasets::Vector{<:Vector{<:Integer}}, column_names::Vector{String})
    # Get unique non-zero trigger values from all datasets
    all_triggers = vcat(trigger_datasets...)
    unique_triggers = unique(all_triggers)
    non_zero_triggers = filter(x -> x != 0, unique_triggers)

    if isempty(non_zero_triggers)
        # Return empty DataFrame with correct column structure
        empty_cols = [Int[]]  # trigger column
        append!(empty_cols, [Int[] for _ in column_names])  # count columns
        column_symbols = [:trigger; Symbol.(column_names)]
        return DataFrame([col => data for (col, data) in zip(column_symbols, empty_cols)]...)
    end

    # Count occurrences of each trigger value across all datasets
    sorted_triggers = sort(non_zero_triggers)
    result_data = Vector{Vector{Int}}(undef, length(column_names) + 1)
    result_data[1] = sorted_triggers

    for (i, dataset) in enumerate(trigger_datasets)
        # Count all triggers in one pass using a dictionary
        trigger_counts = Dict{Int, Int}()
        for val in dataset
            if val != 0
                trigger_counts[val] = get(trigger_counts, val, 0) + 1
            end
        end
        
        # Fill counts array in trigger order
        counts = [get(trigger_counts, trigger, 0) for trigger in sorted_triggers]
        result_data[i + 1] = counts
    end

    # Create DataFrame
    column_symbols = [:trigger; Symbol.(column_names)]
    return DataFrame([col => data for (col, data) in zip(column_symbols, result_data)]...)
end

"""
    _trigger_count_impl(trigger_datasets, column_names; print_table=true, title="Trigger Count Summary", 
                        headers=nothing, note=nothing)

Unified helper function for counting triggers across one or more datasets with flexible printing.

# Arguments
- `trigger_datasets::Vector{<:Vector{<:Integer}}`: Vector of trigger datasets to analyze
- `column_names::Vector{String}`: Names for the count columns in the output DataFrame
- `print_table::Bool`: Whether to print formatted table (default: true)
- `title::String`: Title for the printed table (default: "Trigger Count Summary")
- `headers::Union{Nothing,Vector{String}}`: Custom headers for table (auto-generated if nothing)
- `note::Union{Nothing,String}`: Optional note to display after table

# Returns
- `DataFrame`: DataFrame with 'trigger' column and count columns

# Examples
```julia
# Single dataset
df = _trigger_count_impl([triggers], ["count"])

# Multiple datasets (e.g., raw vs cleaned)
df = _trigger_count_impl([raw_triggers, cleaned_triggers], ["raw_count", "cleaned_count"])
```
"""
function _trigger_count_impl(
    trigger_datasets::Vector{<:Vector{<:Integer}},
    column_names::Vector{String};
    print_table::Bool = true,
    title::String = "Trigger Count Summary",
    headers::Union{Nothing,Vector{String}} = nothing,
    note::Union{Nothing,String} = nothing
)
    result_df = _count_triggers_core(trigger_datasets, column_names)

    if isempty(result_df.trigger)
        if print_table
            if length(trigger_datasets) > 1
                @minimal_warning "No non-zero triggers found in the data."
            else
                println("No non-zero triggers found in the data.")
            end
        end
        return result_df
    end

    # Print table if requested
    if print_table
        # Generate headers if not provided
        if headers === nothing
            headers = ["Trigger"; [uppercasefirst(replace(name, "_" => " ")) for name in column_names]]
        end
        
        pretty_table(
            result_df, 
            title = title, 
            header = headers, 
            alignment = [:r for _ in 1:length(headers)], 
            crop = :none
        )
        println()
        
        # Print note if provided
        if note !== nothing
            println(note)
        end
    end

    return result_df
end

# Backward compatibility method for old test signature
function _trigger_count_impl(
    trigger_data::Vector{<:Integer};
    print_table::Bool = true,
    title::String = "Trigger Count Summary"
)
    return _trigger_count_impl([trigger_data], ["count"], print_table=print_table, title=title)
end


"""
    trigger_count(dat::BiosemiDataFormat.BiosemiData; print_table::Bool = true)::DataFrame

Count trigger occurrences in BioSemi data comparing raw vs cleaned counts.

This function provides a comprehensive view of trigger counts in BioSemi data by
comparing raw trigger counts with cleaned counts (onset-only). This comparison
helps verify trigger cleaning operations and identify sustained vs onset triggers.

# Arguments
- `dat::BiosemiDataFormat.BiosemiData`: The BioSemiData object containing EEG data.
- `print_table::Bool`: Whether to print the trigger count table (default: true).

# Returns
A DataFrame with columns `trigger`, `raw_count`, and `cleaned_count` showing trigger values and their counts in both raw and cleaned data, excluding zero values.

# Examples
```julia
# Get trigger counts and print table
trigger_counts = trigger_count(dat)

# Get trigger counts without printing
trigger_counts = trigger_count(dat, print_table = false)
```
"""
function trigger_count(dat::BiosemiDataFormat.BiosemiData; print_table::Bool = true)::DataFrame
    # Get cleaned trigger data (onset detection only)
    cleaned_triggers = _clean_triggers(dat.triggers.raw)
    return _trigger_count_impl(
        [dat.triggers.raw, cleaned_triggers], 
        ["raw_count", "cleaned_count"],
        print_table = print_table,
        title = "Trigger Count Summary (Raw vs Cleaned)",
        note = "Note: Cleaned counts show only trigger onset events (sustained signals converted to single onsets)"
    )
end

# =============================================================================
# TRIGGER SEQUENCE SEARCH FUNCTIONS
# =============================================================================

"""
    search_sequences(array, sequences)

Return indices of any of the specified sequences within an array, handling wildcards and ranges.

# Arguments
- `array`: Array to search within
- `sequences::Vector{Vector{Union{Int,Symbol,UnitRange{Int}}}}`: Vector of sequences to find

# Returns
- `Vector{Int}`: Indices where any of the sequences start

# Examples
```julia
# Single sequence
idx = search_sequences([1, 2, 3, 4, 5], [[1, 2, 3]])

# Multiple sequences (OR logic)
idx = search_sequences([1, 2, 1, 3, 1, 3, 1], [[1, 2, 1], [1, 3, 1]])

# Wildcards
idx = search_sequences([1, 2, 3, 4, 5], [[1, :any, 3]])

# Ranges
idx = search_sequences([1, 2, 3, 4, 5], [[1:3], [5:5]])

# Mixed sequences
idx = search_sequences([1, 2, 3, 4, 5], [[1, 2:4, 5]])
```
"""
function search_sequences(array, sequences::Vector{Vector{Union{Int,Symbol,UnitRange{Int}}}})
    isempty(array) && return Int[]
    isempty(sequences) && return Int[]

    # Optimize common case of single sequence
    if length(sequences) == 1
        return search_single_sequence(array, sequences[1])
    end

    all_indices = Int[]
    for sequence in sequences
        indices = search_single_sequence(array, sequence)
        append!(all_indices, indices)
    end

    return unique(all_indices)
end

"""
    search_single_sequence(array, sequence)

Return indices of a single sequence within an array, handling wildcards and ranges.

# Arguments
- `array`: Array to search within
- `sequence::Vector{Union{Int,Symbol,UnitRange{Int}}}`: Sequence pattern to find

# Returns
- `Vector{Int}`: Indices where the sequence starts

# Notes
- Supports wildcards (`:any`) and ranges (`1:3`)
- First element cannot be a wildcard
- Single wildcard sequences are not supported
"""
function search_single_sequence(array, sequence::Vector{Union{Int,Symbol,UnitRange{Int}}})
    isempty(array) && return Int[]
    isempty(sequence) && @minimal_error_throw("Sequence cannot be empty")

    # Handle single trigger case
    if length(sequence) == 1
        if sequence[1] isa Int
            return search_sequence(array, sequence[1])
        elseif sequence[1] isa UnitRange
            return search_trigger_ranges(array, [sequence[1]])
        else
            @minimal_error_throw("Single wildcard sequences not supported")
        end
    end

    # Find starting positions for the first trigger
    first_trigger = sequence[1]
    if first_trigger isa Symbol
        @minimal_error_throw("First element in sequence cannot be a wildcard")
    elseif first_trigger isa UnitRange
        # For ranges, find all possible starting positions
        idx_start_positions = search_trigger_ranges(array, [first_trigger])
    else
        # Use search_sequence to find sequence starts
        idx_start_positions = search_sequence(array, first_trigger)
    end

    # sequence of values
    idx_positions = Int[]
    for idx in idx_start_positions
        good_sequence = true
        for seq = 1:(length(sequence)-1)
            if (idx + seq) > length(array)
                good_sequence = false
                break
            end

            expected_trigger = sequence[seq+1]
            actual_trigger = array[idx+seq]

            if expected_trigger isa Int
                # Exact match required
                if actual_trigger != expected_trigger
                    good_sequence = false
                    break
                end
            elseif expected_trigger == :any
                # Wildcard - matches any value
                continue
            elseif expected_trigger isa UnitRange
                # Range - check if actual trigger is in range
                if !(actual_trigger in expected_trigger)
                    good_sequence = false
                    break
                end
            else
                @minimal_error_throw("Unsupported trigger type: $expected_trigger")
            end
        end
        if good_sequence
            push!(idx_positions, idx)
        end
    end

    return idx_positions
end

"""
    search_trigger_ranges(array, trigger_ranges)

Return indices where triggers match any of the specified ranges.

# Arguments
- `array`: Array to search within
- `trigger_ranges::Vector{UnitRange{Int}}`: Vector of integer ranges to match

# Returns
- `Vector{Int}`: Indices where triggers match any of the ranges

# Examples
```julia
idx = search_trigger_ranges([1, 2, 3, 4, 5, 6], [1:3, 5:6])  # Returns indices of 1,2,3,5,6
```
"""
function search_trigger_ranges(array, trigger_ranges::Vector{UnitRange{Int}})
    isempty(array) && return Int[]
    isempty(trigger_ranges) && return Int[]

    idx_positions = Int[]
    for (i, trigger) in enumerate(array)
        for range in trigger_ranges
            if trigger in range
                push!(idx_positions, i)
                break  # Only add each position once
            end
        end
    end

    return idx_positions
end

"""
    search_sequence(array::AbstractVector, sequence::Int) -> Vector{Int}

Find starting indices of a sequence in an array (onset detection).

# Arguments
- `array::AbstractVector`: Array to search within
- `sequence::Int`: Trigger value to find (onset detection)

# Returns
- `Vector{Int}`: Indices where sequence starts (onsets only)

# Notes
- Uses onset detection (finds only first occurrence of sustained triggers)
- Returns intersection of value matches and onset positions
"""
function search_sequence(array::AbstractVector, sequence::Int)
    isempty(array) && return Int[]
    
    # Find all positions where the trigger value matches
    value_matches = findall(array .== sequence)
    
    # Find all positions where there's an onset (value increases from previous)
    onset_positions = findall(diff(vcat(0, array)) .>= 1)
    
    # Return intersection (positions that are both value matches and onsets)
    return intersect(value_matches, onset_positions)
end
