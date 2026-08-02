"""
    lrp(erp_left::ErpData, erp_right::ErpData; channel_selection = channels())
    lrp(erps::Vector{ErpData}, condition_pairs::Vector{Tuple{Int,Int}}; channel_selection = channels())
    lrp(file_pattern::String, condition_pairs::Vector{Tuple{Int,Int}}; input_dir = pwd(),
    channel_selection = channels(), participant_selection = participants(), output_dir = nothing)

Calculate the lateralized readiness potential (LRP) from ERP data.

- **Two-argument form**: takes paired left/right ERP datasets, returns one `ErpData` with LRP values.
- **`Vector{ErpData}` form**: computes LRP for multiple condition pairs at once.
- **`file_pattern` form**: batch-processes JLD2 files across participants.

The LRP is a measure of lateralized motor preparation. For each channel pair (e.g., C3/C4):
- LRP_C3 = 0.5 × ((C3_right - C4_right) + (C4_left - C3_left))
- LRP_C4 = 0.5 × ((C4_right - C3_right) + (C3_left - C4_left))

This formula isolates lateralized activity by averaging the difference between
contralateral and ipsilateral activation across both hemispheres.

# References
- Coles, M. G. H. (1989). Modern mind-brain reading. *Psychophysiology*, 26(3), 251-269.
- de Jong, R. et al. (1988). Use of partial stimulus information. *JEP:HPP*, 14(4), 682-692.
- Oostenveld, R. et al. (2003). Brain symmetry and topographic analysis. *Clin. Neurophysiology*, 114(7), 1194-1202.

# Arguments (two-argument form)
- `erp_left::ErpData`: ERP data for left-hand responses (ipsilateral activation)
- `erp_right::ErpData`: ERP data for right-hand responses (contralateral activation)
- `channel_selection`: Channel predicate selecting left/odd hemisphere channels.
  Default `channels()` auto-detects all odd/even pairs (C3/C4, C1/C2, Fp1/Fp2, …)

# Examples
```julia
# Auto-detect all lateral pairs (C3/C4, C1/C2, etc.)
lrp_data = lrp(erps[1], erps[2])

# Only C3/C4 and CP3/CP4
lrp_data = lrp(erps[1], erps[2], channel_selection = channels([:C3, :CP3]))

# Multiple condition pairs from one participant file
pairs = [(i, i+1) for i in 1:2:15]
lrp_results = lrp(erps, pairs)

# Batch across all participants
lrp("erps_unrejected", [(1, 2), (3, 4)], input_dir = "/data/study1")
```
"""
function lrp(erp_left::ErpData, erp_right::ErpData; channel_selection::Function = channels())::ErpData

    @info "Calculating lateralized readiness potential (LRP)"

    # Validate inputs
    _have_same_structure(erp_left, erp_right) || @minimal_error("Left and right ERPs have inconsistent structure")

    # Get selected left/odd channels and find their right/even pairs
    pairs = _get_channel_pairs_from_selection(erp_left, erp_right, channel_selection)

    if isempty(pairs)
        @minimal_error(
            "No valid lateral channel pairs found. Check that selected channels have odd numbers and their even counterparts exist."
        )
    end

    @info "Using $(length(pairs)) channel pair(s): $(pairs)"

    # Calculate LRP
    lrp_data = _calculate_lrp(erp_left, erp_right, pairs)

    @info "LRP calculation complete"
    return lrp_data
end


# === INTERNAL HELPER FUNCTIONS ===

"""
Parse channel label to extract letters and digits.
Returns (letters, digit) or (letters, nothing) if no digit found.
"""
function _parse_channel_label(label::Symbol)
    label_str = String(label)
    digits_only = filter(isdigit, label_str)
    letters_only = filter(isletter, label_str)

    digit = isempty(digits_only) ? nothing : parse(Int, digits_only)
    return (letters_only, digit)
end


"""
Get channel pairs from a channel selection predicate.

Takes a channel selection predicate, applies it to get left/odd channels,
validates they are odd-numbered, and pairs them with their right/even counterparts.
"""
function _get_channel_pairs_from_selection(erp_left::ErpData, erp_right::ErpData, channel_selection::Function)::Vector{Tuple{Symbol,Symbol}}

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
function _detect_all_lateral_pairs(left_labels::Vector{Symbol}, right_labels::Vector{Symbol})::Vector{Tuple{Symbol,Symbol}}

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


# === LRP-SPECIFIC VALIDATION ===

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
function _default_lrp_output_dir(input_dir::String, condition_pairs::Vector{Tuple{Int,Int}})
    pairs_str = join(["$(l)-$(r)" for (l, r) in condition_pairs], "_")
    joinpath(input_dir, "lrp_$(pairs_str)")
end

# === LRP-SPECIFIC PROCESSING ===

"""
Process a single ERP file through LRP calculation pipeline.
Returns BatchResult with success/failure info.
"""
function _process_lrp_file(filepath::String, output_path::String, condition_pairs::Vector{Tuple{Int,Int}}, channel_selection::Function)
    filename = basename(filepath)

    # Read data (using read_data which finds by type)
    erps_data = read_data(filepath)

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

# === MAIN API FUNCTION FOR BATCH PROCESSING ===


function lrp(
    file_pattern::String,
    condition_pairs::Vector{Tuple{Int,Int}};
    input_dir::String = pwd(),
    channel_selection::Function = channels(),
    participant_selection::Function = participants(),
    output_dir::Union{String,Nothing} = nothing,
)

    # Setup logging
    log_file = "lrp.log"
    setup_global_logging(log_file)

    try
        @info "Batch LRP calculation started at $(now())"
        @log_call "lrp"

        # Validation
        error_msg = _validate_input_dir(input_dir)
        if !isnothing(error_msg)
            @minimal_error(error_msg)
        end

        error_msg = _validate_lrp_params(condition_pairs)
        if !isnothing(error_msg)
            @minimal_error(error_msg)
        end

        # Setup directories
        output_dir = something(output_dir, _default_lrp_output_dir(input_dir, condition_pairs))
        mkpath(output_dir)

        # Find files
        files = _find_batch_files(file_pattern, input_dir, participant_selection)

        if isempty(files)
            @minimal_warning "No JLD2 files found matching pattern '$file_pattern' in $input_dir"
            return nothing
        end

        @info "Found $(length(files)) JLD2 files matching pattern '$file_pattern'"
        @info "Condition pairs: $condition_pairs"
        @info "Channel selection: $(channel_selection == channels() ? "all lateral pairs" : "custom")"

        # Create processing function with captured parameters
        process_fn = (input_path, output_path) -> _process_lrp_file(input_path, output_path, condition_pairs, channel_selection)

        # Execute batch operation
        results = _run_batch_operation(process_fn, files, input_dir, output_dir; operation_name = "Calculating LRP")

        _log_batch_summary(results, output_dir)

    finally
        _cleanup_logging(log_file, output_dir)
    end
end



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
            @minimal_error("Left condition index $left_cond out of range (1-$(length(erps_filtered)))")
        end
        if right_cond < 1 || right_cond > length(erps_filtered)
            @minimal_error("Right condition index $right_cond out of range (1-$(length(erps_filtered)))")
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

    # Pre-allocate matrix for LRP data (2 channels per pair)
    lrp_matrix = zeros(n_timepoints, 2 * length(pairs))
    lrp_labels = Symbol[]

    # Calculate LRP for each pair
    for (idx, (ch_left, ch_right)) in enumerate(pairs)
        # Extract data for each channel from each condition
        ch_left_in_left = erp_left.data[!, ch_left]
        ch_right_in_left = erp_left.data[!, ch_right]
        ch_left_in_right = erp_right.data[!, ch_left]
        ch_right_in_right = erp_right.data[!, ch_right]

        # Calculate column indices for interleaved storage (C3, C4, CP3, CP4, ...)
        left_col = 2 * idx - 1
        right_col = 2 * idx

        # Calculate LRP using the double-subtraction formula
        # LRP for left channel (e.g., C3)
        @views left_view = lrp_matrix[:, left_col]
        @. left_view = 0.5 * ((ch_left_in_right - ch_right_in_right) + (ch_right_in_left - ch_left_in_left))

        # LRP for right channel (e.g., C4) is mathematically identical to the inverse of the left channel
        @views right_view = lrp_matrix[:, right_col]
        @. right_view = -left_view

        # Build channel labels (interleaved to match matrix structure)
        push!(lrp_labels, ch_left)
        push!(lrp_labels, ch_right)
    end

    # Create output DataFrame with metadata from left dataset
    meta_cols = meta_labels(erp_left)
    lrp_df = DataFrame()

    # Copy metadata columns
    for col in meta_cols
        lrp_df[!, col] = copy(erp_left.data[!, col])
    end

    # Add LRP channel data
    for (idx, label) in enumerate(lrp_labels)
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
    lrp_layout = _create_lrp_layout(erp_left.layout, lrp_labels)

    # Create and return LRP ErpData
    # Use minimum n_epochs as conservative estimate
    min_epochs = min(erp_left.n_epochs, erp_right.n_epochs)
    # Use condition number from left ERP
    return ErpData(
        erp_left.file,
        erp_left.condition,
        "lrp",
        lrp_df,
        lrp_layout,
        erp_left.sample_rate,
        copy(erp_left.analysis_info),
        min_epochs,
    )
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

    return Layout(filtered_layout_df, filtered_neighbours, original_layout.criterion, original_layout.criterion_type)
end
