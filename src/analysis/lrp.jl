"""
    lrp(erp_left, erp_right; channel_selection = channels())

Calculate the lateralized readiness potential (LRP) from two ERP datasets.

The lateralized readiness potential is a measure of lateralized motor preparation
derived from event-related potentials recorded over motor cortex. This implementation
follows the approach described in:

- Coles, M. G. H. (1989). Modern mind-brain reading: Psychophysiology, physiology, 
  and cognition. Psychophysiology, 26(3), 251-269.
- de Jong, R., Wierda, M., Mulder, G., & Mulder, L. J. (1988). Use of partial 
  stimulus information in response processing. Journal of Experimental Psychology: 
  Human Perception and Performance, 14(4), 682-692.
- Oostenveld, R., Stegeman, D. F., Praamstra, P., & van Oosterom, A. (2003). 
  Brain symmetry and topographic analysis of lateralized event-related potentials. 
  Clinical Neurophysiology, 114(7), 1194-1202.

# Arguments
- `erp_left::ErpData`: ERP data for left-hand responses (ipsilateral activation)
- `erp_right::ErpData`: ERP data for right-hand responses (contralateral activation)
- `channel_selection::Function`: Channel predicate to select left/odd hemisphere channels.
   The function automatically pairs them with corresponding right/even channels.
   Default: `channels()` auto-detects all odd/even pairs (e.g., C3/C4, C1/C2, Fp1/Fp2)

# Returns
- `ErpData`: New ERP dataset containing LRP values for each channel in the pairs

# Description
For each channel pair (e.g., C3/C4), the LRP is calculated as:
- LRP_C3 = 0.5 × ((C3_right - C4_right) + (C4_left - C3_left))
- LRP_C4 = 0.5 × ((C4_right - C3_right) + (C3_left - C4_left))

Where:
- C3_right/C4_right: Activity at C3/C4 for right-hand responses (contralateral/ipsilateral)
- C3_left/C4_left: Activity at C3/C4 for left-hand responses (ipsilateral/contralateral)

This formula isolates lateralized activity by averaging the difference between 
contralateral and ipsilateral activation across both hemispheres, canceling out 
non-lateralized activity.

# Examples
```julia
# Auto-detect all lateral pairs (C3/C4, C1/C2, Fp1/Fp2, etc.)
lrp_data = lrp(erps[1], erps[2])

# Select specific left hemisphere channels (automatically pairs with right)
lrp_data = lrp(erps[1], erps[2], channel_selection = channels([:C3, :CP3]))
# This calculates LRP for C3/C4 and CP3/CP4

# Use channel predicates for pattern matching
lrp_data = lrp(erps[1], erps[2], 
               channel_selection = channels(x -> startswith.(string.(x), "C")))
# Selects all C-channels: C1/C2, C3/C4, C5/C6, etc.

# Full workflow example
using JLD2
erps = load("participant_05_erps.jld2", "erps")
lrp_result = lrp(erps[1], erps[2], channel_selection = channels([:C3]))
jldsave("participant_05_lrp.jld2"; data = lrp_result)
```

# References
See the original FieldTrip implementation: `ft_lateralizedpotential.m`
"""
function lrp(erp_left::ErpData, erp_right::ErpData; channel_selection::Function = channels())::ErpData

    @info "Calculating lateralized readiness potential (LRP)"

    # Validate inputs
    _validate_lrp_inputs(erp_left, erp_right)

    # Get selected left/odd channels and find their right/even pairs
    pairs = _get_channel_pairs_from_selection(erp_left, erp_right, channel_selection)

    if isempty(pairs)
        @minimal_error_throw(
            "No valid lateral channel pairs found. Check that selected channels have odd numbers and their even counterparts exist."
        )
    end

    @info "Using $(length(pairs)) channel pair(s): $(pairs)"

    # Calculate LRP
    lrp_data = _calculate_lrp(erp_left, erp_right, pairs)

    @info "LRP calculation complete"
    return lrp_data
end


#=============================================================================
    INTERNAL HELPER FUNCTIONS
=============================================================================#

"""
Validate that the two ERP datasets are compatible for LRP calculation.
"""
function _validate_lrp_inputs(erp_left::ErpData, erp_right::ErpData)
    # Check sample rates match
    if erp_left.sample_rate != erp_right.sample_rate
        @minimal_error_throw("Sample rates differ: left=$(erp_left.sample_rate) Hz, right=$(erp_right.sample_rate) Hz")
    end

    # Check time points match
    if nrow(erp_left.data) != nrow(erp_right.data)
        @minimal_error_throw("Number of time points differ: left=$(nrow(erp_left.data)), right=$(nrow(erp_right.data))")
    end

    # Check time vectors match
    if !all(erp_left.data.time .≈ erp_right.data.time)
        @minimal_error_throw("Time vectors differ between left and right datasets")
    end

    return nothing
end


"""
Parse channel label to extract letters and digits.
Returns (letters, digit) or (letters, nothing) if no digit found.
"""
function _parse_channel_label(label::Symbol)
    label_str = String(label)
    digits_only = Base.filter(isdigit, label_str)
    letters_only = Base.filter(isletter, label_str)

    digit = isempty(digits_only) ? nothing : parse(Int, digits_only)
    return (letters_only, digit)
end


"""
Get channel pairs from a channel selection predicate.

Takes a channel selection predicate, applies it to get left/odd channels,
validates they are odd-numbered, and pairs them with their right/even counterparts.
"""
function _get_channel_pairs_from_selection(
    erp_left::ErpData,
    erp_right::ErpData,
    channel_selection::Function,
)::Vector{Tuple{Symbol,Symbol}}

    # Get all available channels
    left_labels = channel_labels(erp_left)
    right_labels = channel_labels(erp_right)

    # If default channels() (all channels), auto-detect all lateral pairs
    if channel_selection == channels()
        return _detect_all_lateral_pairs(left_labels, right_labels)
    end

    # Apply the channel selection predicate
    selected_channels = get_selected_channels(erp_left, channel_selection, include_meta = false, include_extra = false)

    # Build pairs from selected channels
    pairs = Tuple{Symbol,Symbol}[]

    for ch_left in selected_channels
        letters, digit = _parse_channel_label(ch_left)

        # Skip if no digits
        if isnothing(digit)
            @minimal_warning "Channel $ch_left has no number, cannot pair. Skipping."
            continue
        end

        # Validate it's odd
        if iseven(digit)
            @minimal_warning "Channel $ch_left is even-numbered (expected odd). Skipping."
            continue
        end

        # Construct the paired channel name (even number = odd + 1)
        ch_right = Symbol(letters * string(digit + 1))

        # Check if the pair exists in both datasets
        if ch_right ∉ left_labels
            @minimal_warning "Right pair $ch_right for $ch_left not found in left dataset. Skipping."
            continue
        end
        if ch_right ∉ right_labels
            @minimal_warning "Right pair $ch_right for $ch_left not found in right dataset. Skipping."
            continue
        end

        push!(pairs, (ch_left, ch_right))
    end

    return pairs
end


"""
Detect all lateral channel pairs based on odd/even numbering (e.g., C3/C4, C1/C2).
Used when channel_selection = channels() (default).
"""
function _detect_all_lateral_pairs(
    left_labels::Vector{Symbol},
    right_labels::Vector{Symbol},
)::Vector{Tuple{Symbol,Symbol}}

    # Use intersection to ensure channels exist in both datasets
    common_labels = intersect(left_labels, right_labels)
    pairs = Tuple{Symbol,Symbol}[]

    for label in common_labels
        letters, digit = _parse_channel_label(label)

        # Skip if no digits or if even
        isnothing(digit) && continue
        iseven(digit) && continue

        # Construct the paired channel name (even number)
        paired_label = Symbol(letters * string(digit + 1))

        # Check if the pair exists in both datasets
        if paired_label ∈ common_labels
            push!(pairs, (label, paired_label))
        end
    end

    return pairs
end


"""
Batch LRP calculation for multiple participants from JLD2 files.
"""

#=============================================================================
    LRP-SPECIFIC VALIDATION
=============================================================================#

"""Validate LRP-specific parameters, returning error message or nothing."""
function _validate_lrp_params(condition_pairs::Vector{Tuple{Int,Int}})
    isempty(condition_pairs) && return "Condition pairs cannot be empty"

    for (left, right) in condition_pairs
        left < 1 && return "Condition indices must be positive, got left=$left"
        right < 1 && return "Condition indices must be positive, got right=$right"
    end

    return nothing
end

"""Generate default output directory name for LRP operation."""
function _default_lrp_output_dir(input_dir::String, pattern::String, condition_pairs::Vector{Tuple{Int,Int}})
    pairs_str = join(["$(l)-$(r)" for (l, r) in condition_pairs], "_")
    joinpath(input_dir, "lrp_$(pairs_str)")
end

#=============================================================================
    LRP-SPECIFIC PROCESSING
=============================================================================#

"""
Process a single ERP file through LRP calculation pipeline.
Returns BatchResult with success/failure info.
"""
function _process_lrp_file(
    filepath::String,
    output_path::String,
    condition_pairs::Vector{Tuple{Int,Int}},
    channel_selection::Function,
)
    filename = basename(filepath)

    # Load data (using load_data which finds by type)
    erps_data = load_data(filepath)
    
    if isnothing(erps_data)
        return BatchResult(false, filename, "No data variables found")
    end
    
    # Validate that data is Vector{ErpData}
    if !(erps_data isa Vector{<:ErpData})
        return BatchResult(false, filename, "Invalid data type: expected Vector{ErpData}, got $(typeof(erps_data))")
    end

    if isempty(erps_data)
        return BatchResult(false, filename, "Empty erps array")
    end

    # Calculate LRP for all condition pairs
    try
        lrp_results = lrp(erps_data, condition_pairs; channel_selection = channel_selection)

        # Save results
        jldsave(output_path; data = lrp_results)

        return BatchResult(true, filename, "Calculated LRP for $(length(lrp_results)) pair(s)")
    catch e
        return BatchResult(false, filename, "Error: $(sprint(showerror, e))")
    end
end

#=============================================================================
    MAIN API FUNCTION FOR BATCH PROCESSING
=============================================================================#

"""
    lrp(file_pattern::String, condition_pairs::Vector{Tuple{Int,Int}};
        input_dir::String = pwd(),
        channel_selection::Function = channels(),
        participants::Union{Int, Vector{Int}, Nothing} = nothing,
        output_dir::Union{String, Nothing} = nothing)

Calculate LRP from ERP data in JLD2 files and save to a new directory.

This function processes multiple participant files at once, calculating the
lateralized readiness potential for specified condition pairs.

# Arguments
- `file_pattern::String`: Pattern to match files (e.g., "erps", "erps_cleaned")
- `condition_pairs::Vector{Tuple{Int,Int}}`: Pairs of condition indices (left, right)
- `input_dir::String`: Input directory containing JLD2 files (default: current directory)
- `channel_selection::Function`: Channel predicate for selecting left/odd channels (default: all lateral pairs)
- `participants::Union{Int, Vector{Int}, Nothing}`: Participant number(s) to process (default: all)
- `output_dir::Union{String, Nothing}`: Output directory (default: creates subdirectory)

# Examples
```julia
# Calculate LRP for all participants with 16 conditions (odd=left, even=right)
pairs = [(i, i+1) for i in 1:2:15]
lrp("erps_cleaned", pairs)

# Specific participants only
lrp("erps_cleaned", [(1,2), (3,4)], participants=[1, 2, 3])

# Only C3/C4 pair
lrp("erps_cleaned", [(1,2), (3,4)], channel_selection=channels([:C3]))

# Custom output directory
lrp("erps_cleaned", [(1,2)], output_dir="/path/to/output")

# Full example workflow
pairs = [(i, i+1) for i in 1:2:15]
lrp("erps_cleaned", pairs, 
    input_dir="/data/study1",
    channel_selection=channels([:C3, :CP3]),
    participants=1:20)
```

# Output
- Creates new directory with LRP data files
- Each output file contains "lrp" variable with Vector{ErpData}
- Log file saved to output directory
"""
function lrp(
    file_pattern::String,
    condition_pairs::Vector{Tuple{Int,Int}};
    input_dir::String = pwd(),
    channel_selection::Function = channels(),
    participants::Union{Int,Vector{Int},Nothing} = nothing,
    output_dir::Union{String,Nothing} = nothing,
)

    # Setup logging
    log_file = "lrp.log"
    setup_global_logging(log_file)

    try
        @info "Batch LRP calculation started at $(now())"
        @log_call "lrp" (file_pattern, condition_pairs)

        # Validation
        if (error_msg = _validate_input_dir(input_dir)) !== nothing
            @minimal_error_throw(error_msg)
        end

        if (error_msg = _validate_lrp_params(condition_pairs)) !== nothing
            @minimal_error_throw(error_msg)
        end

        # Setup directories
        output_dir = something(output_dir, _default_lrp_output_dir(input_dir, file_pattern, condition_pairs))
        mkpath(output_dir)

        # Find files
        files = _find_batch_files(file_pattern, input_dir; participants)

        if isempty(files)
            @minimal_warning "No JLD2 files found matching pattern '$file_pattern' in $input_dir"
            return nothing
        end

        @info "Found $(length(files)) JLD2 files matching pattern '$file_pattern'"
        @info "Condition pairs: $condition_pairs"
        @info "Channel selection: $(channel_selection == channels() ? "all lateral pairs" : "custom")"

        # Create processing function with captured parameters
        process_fn =
            (input_path, output_path) -> _process_lrp_file(input_path, output_path, condition_pairs, channel_selection)

        # Execute batch operation
        results = _run_batch_operation(process_fn, files, input_dir, output_dir; operation_name = "Calculating LRP")

        _log_batch_summary(results, output_dir)

    finally
        _cleanup_logging(log_file, output_dir)
    end
end


"""
    lrp(erps::Vector{ErpData}, condition_pairs::Vector{Tuple{Int,Int}}; channel_selection = channels())

Calculate LRP for multiple condition pairs from a vector of ErpData.

This is a convenience function for batch processing multiple condition pairs.
Useful when you have many conditions organized as left/right hand response pairs.

# Arguments
- `erps::Vector{ErpData}`: Vector of ERP data for all conditions
- `condition_pairs::Vector{Tuple{Int,Int}}`: Pairs of condition indices (left, right)
- `channel_selection::Function`: Channel predicate to select left/odd channels (automatically pairs with right/even)

# Returns
- `Vector{ErpData}`: Vector of LRP results, one for each condition pair

# Examples
```julia
# For 16 conditions where odd=left, even=right
erps = load("participant_05_erps.jld2", "erps")

# Define pairs: (1,2), (3,4), (5,6), (7,8), etc.
pairs = [(i, i+1) for i in 1:2:15]

# Calculate all LRPs at once (auto-detect all lateral pairs)
lrp_results = lrp(erps, pairs)

# Only calculate LRP for C3/C4 and CP3/CP4
lrp_results = lrp(erps, pairs, channel_selection = channels([:C3, :CP3]))

# Use pattern matching for all C-channels
lrp_results = lrp(erps, pairs, 
                  channel_selection = channels(x -> startswith.(string.(x), "C")))

# Save results
jldsave("participant_05_lrp.jld2"; data = lrp_results)
```
"""
function lrp(
    erps::Vector{ErpData},
    condition_pairs::Vector{Tuple{Int,Int}};
    condition_selection::Function = conditions(),
    channel_selection::Function = channels(),
)::Vector{ErpData}

    # Apply condition_selection first
    erps_filtered = erps[get_selected_conditions(erps, condition_selection)]
    
    @info "Calculating LRP for $(length(condition_pairs)) condition pair(s) from $(length(erps_filtered)) condition(s)"

    lrp_results = ErpData[]

    for (idx, (left_cond, right_cond)) in enumerate(condition_pairs)
        # Validate condition indices (now referring to filtered array)
        if left_cond < 1 || left_cond > length(erps_filtered)
            @minimal_error_throw("Left condition index $left_cond out of range (1-$(length(erps_filtered)))")
        end
        if right_cond < 1 || right_cond > length(erps_filtered)
            @minimal_error_throw("Right condition index $right_cond out of range (1-$(length(erps_filtered)))")
        end

        @info "  Processing pair $idx: condition $left_cond (left) vs $right_cond (right)"

        # Calculate LRP for this pair
        lrp_data = lrp(erps_filtered[left_cond], erps_filtered[right_cond]; channel_selection = channel_selection)

        # Update condition information to reflect the pair number (set struct fields)
        lrp_data.condition = idx
        lrp_data.condition_name = "lrp_$(left_cond)_$(right_cond)"

        push!(lrp_results, lrp_data)
    end

    @info "LRP calculation complete for all $(length(condition_pairs)) pair(s)"
    return lrp_results
end


"""
Calculate LRP for each channel pair using the double-subtraction formula.

For a channel pair (C3, C4):
- LRP_C3 = 0.5 × ((C3_right - C4_right) + (C4_left - C3_left))
- LRP_C4 = 0.5 × ((C4_right - C3_right) + (C3_left - C4_left))
"""
function _calculate_lrp(erp_left::ErpData, erp_right::ErpData, pairs::Vector{Tuple{Symbol,Symbol}})::ErpData

    n_timepoints = nrow(erp_left.data)
    n_pairs = length(pairs)

    # Pre-allocate matrix for LRP data (2 channels per pair)
    lrp_matrix = zeros(n_timepoints, 2 * n_pairs)
    lrp_labels = Symbol[]

    # Calculate LRP for each pair
    for (idx, (ch_left, ch_right)) in enumerate(pairs)
        # Extract data for each channel from each condition
        ch_left_in_left = erp_left.data[!, ch_left]
        ch_right_in_left = erp_left.data[!, ch_right]
        ch_left_in_right = erp_right.data[!, ch_left]
        ch_right_in_right = erp_right.data[!, ch_right]

        # Calculate LRP using the double-subtraction formula
        # LRP for left channel (e.g., C3)
        lrp_matrix[:, idx] = 0.5 .* ((ch_left_in_right .- ch_right_in_right) .+ (ch_right_in_left .- ch_left_in_left))

        # LRP for right channel (e.g., C4) 
        lrp_matrix[:, idx+n_pairs] =
            0.5 .* ((ch_right_in_right .- ch_left_in_right) .+ (ch_left_in_left .- ch_right_in_left))

        # Build channel labels
        push!(lrp_labels, ch_left)
        push!(lrp_labels, ch_right)
    end

    # Remove duplicate labels (keep unique in order)
    unique_labels = Symbol[]
    unique_indices = Int[]
    for (idx, label) in enumerate(lrp_labels)
        if label ∉ unique_labels
            push!(unique_labels, label)
            push!(unique_indices, idx)
        end
    end

    # Select only unique columns
    lrp_matrix = lrp_matrix[:, unique_indices]

    # Create output DataFrame with metadata from left dataset
    meta_cols = meta_labels(erp_left)
    lrp_df = DataFrame()

    # Copy metadata columns
    for col in meta_cols
        lrp_df[!, col] = copy(erp_left.data[!, col])
    end

    # Add LRP channel data
    for (idx, label) in enumerate(unique_labels)
        lrp_df[!, label] = lrp_matrix[:, idx]
    end

    # Remove condition/condition_name columns if they exist (they're in struct now)
    cols_to_remove = [:condition, :condition_name, :n_epochs]
    for col in cols_to_remove
        if hasproperty(lrp_df, col)
            select!(lrp_df, Not(col))
        end
    end

    # Create layout with only the LRP channels
    lrp_layout = _create_lrp_layout(erp_left.layout, unique_labels)

    # Create and return LRP ErpData
    # Use minimum n_epochs as conservative estimate
    min_epochs = min(erp_left.n_epochs, erp_right.n_epochs)
    # LRP doesn't have a condition number, use 0 as placeholder
    return ErpData(erp_left.file, 0, "lrp", lrp_df, lrp_layout, erp_left.sample_rate, copy(erp_left.analysis_info), min_epochs)
end


"""
Create a layout containing only the specified LRP channels.
"""
function _create_lrp_layout(original_layout::Layout, lrp_channels::Vector{Symbol})::Layout
    # Filter layout data to include only LRP channels
    layout_df = copy(original_layout.data)
    mask = [label ∈ lrp_channels for label in layout_df.label]
    filtered_layout_df = layout_df[mask, :]

    # Filter neighbours if they exist
    filtered_neighbours = if !isnothing(original_layout.neighbours)
        filtered_dict = OrderedDict{Symbol,Neighbours}()
        for channel in lrp_channels
            if haskey(original_layout.neighbours, channel)
                filtered_dict[channel] = original_layout.neighbours[channel]
            end
        end
        isempty(filtered_dict) ? nothing : filtered_dict
    else
        nothing
    end

    return Layout(filtered_layout_df, filtered_neighbours, original_layout.criterion)
end
