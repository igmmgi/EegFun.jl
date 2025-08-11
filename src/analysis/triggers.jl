"""
    _clean_triggers(trigger_data::Vector{<:Integer})::Vector{<:Integer}

Cleans trigger data by detecting only the onset (first occurrence) of each trigger value.
Converts sustained trigger signals into single onset events.

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
    cleaned = zeros(Int, length(trigger_data))
    cleaned[onset_indices] = trigger_data[onset_indices]

    return cleaned

end



"""
    trigger_count(dat::ContinuousData; print_table::Bool = true)::DataFrame

Counts the number of occurrences of each trigger value in the data and optionally prints a formatted table.

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
    return _trigger_count_impl(dat.data.triggers, print_table = print_table)
end

# Helper function for trigger counting logic
function _trigger_count_impl(
    trigger_data::Vector{<:Integer};
    print_table::Bool = true,
    title::String = "Trigger Count Summary",
)
    # Get unique non-zero trigger values
    unique_triggers = unique(trigger_data)
    non_zero_triggers = filter(x -> x != 0, unique_triggers)

    if isempty(non_zero_triggers)
        if print_table
            println("No non-zero triggers found in the data.")
        end
        return DataFrame()
    end

    # Count occurrences of each trigger value
    trigger_values = Int[]
    trigger_counts = Int[]

    for trigger in sort(non_zero_triggers)
        push!(trigger_values, trigger)
        push!(trigger_counts, count(x -> x == trigger, trigger_data))
    end

    # Create DataFrame
    result_df = DataFrame(trigger = trigger_values, count = trigger_counts)

    # Print table if requested
    if print_table
        pretty_table(result_df, title = title, header = ["Trigger", "Count"], alignment = [:r, :r], crop = :none)
        println()
    end

    return result_df
end

# Helper function for BioSemi data with both raw and cleaned counts
function _trigger_count_biosemi_impl(
    raw_triggers::Vector{<:Integer},
    cleaned_triggers::Vector{<:Integer};
    print_table::Bool = true,
)
    # Get unique non-zero trigger values from both raw and cleaned data
    unique_triggers = unique(vcat(raw_triggers, cleaned_triggers))
    non_zero_triggers = filter(x -> x != 0, unique_triggers)

    if isempty(non_zero_triggers)
        if print_table
            @minimal_warning "No non-zero triggers found in the data."
        end
        return DataFrame()
    end

    # Count occurrences of each trigger value in both raw and cleaned data
    trigger_values = Int[]
    raw_counts = Int[]
    cleaned_counts = Int[]

    for trigger in sort(non_zero_triggers)
        push!(trigger_values, trigger)
        push!(raw_counts, count(x -> x == trigger, raw_triggers))
        push!(cleaned_counts, count(x -> x == trigger, cleaned_triggers))
    end

    # Create DataFrame
    result_df = DataFrame(trigger = trigger_values, raw_count = raw_counts, cleaned_count = cleaned_counts)

    # Print table if requested
    if print_table
        pretty_table(
            result_df,
            title = "Trigger Count Summary (Raw vs Cleaned)",
            header = ["Trigger", "Raw Count", "Cleaned Count"],
            alignment = [:r, :r, :r],
            crop = :none,
        )
        println()
        println("Note: Cleaned counts show only trigger onset events (sustained signals converted to single onsets)")
    end

    return result_df
end


"""
    trigger_count(dat::BioSemiBDF.BioSemiData; print_table::Bool = true)::DataFrame

Counts the number of occurrences of each trigger value in BioSemi data (both raw and cleaned) and optionally prints a formatted table.

# Arguments
- `dat::BioSemiBDF.BioSemiData`: The BioSemiData object containing EEG data.
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
function trigger_count(dat::BioSemiBDF.BioSemiData; print_table::Bool = true)::DataFrame
    # Get cleaned trigger data (onset detection only)
    cleaned_triggers = _clean_triggers(dat.triggers.raw)
    return _trigger_count_biosemi_impl(dat.triggers.raw, cleaned_triggers, print_table = print_table)
end
