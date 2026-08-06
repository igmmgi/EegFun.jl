"""
Batch computation of condition difference waves for ERP data.
"""

# === DIFFERENCE-SPECIFIC HELPERS ===

"""Generate default output directory name for difference operation."""
function _condition_difference_default_output_dir(
    input_dir::String,
    pattern::String,
    pairs::Union{Vector{Tuple{Int,Int}},Vector{Vector{Int}}},
)
    pairs_str = join([join(pair, "-") for pair in pairs], "_")
    joinpath(input_dir, "differences_$(clean_pattern(pattern))_$(pairs_str)")
end

# === DIFFERENCE-SPECIFIC PROCESSING ===

"""
Create a difference wave by subtracting ERP2 from ERP1.
"""
function _create_difference_wave(erp1::ErpData, erp2::ErpData, cond1::Int, cond2::Int, diff_cond::Int)
    # Validate that both ERPs have the same structure
    _have_same_structure(erp1, erp2) || @minimal_error("ERPs have inconsistent structure")

    # Get EEG channels (exclude metadata columns)
    metadata_cols = meta_labels(erp1)
    eeg_channels = setdiff(propertynames(erp1.data), metadata_cols)

    # Create a copy of erp1's data for the difference
    diff_data = copy(erp1.data)

    # Remove condition/condition_name/n_epochs columns if they exist (they're in struct now)
    cols_to_remove = [:condition, :condition_name, :n_epochs]
    for col in cols_to_remove
        if hasproperty(diff_data, col)
            select!(diff_data, Not(col))
        end
    end

    # Subtract EEG channels
    for ch in eeg_channels
        if hasproperty(erp2.data, ch)
            diff_data[!, ch] .-= erp2.data[!, ch]
        else
            @minimal_warning "Channel $ch not found in condition $cond2, keeping original values"
        end
    end

    # Update n_epochs to reflect the minimum (conservative estimate)
    min_epochs = min(erp1.n_epochs, erp2.n_epochs)
    diff_condition_name = "difference_$(cond1)_$(cond2)"

    return ErpData(
        erp1.file,
        diff_cond,
        diff_condition_name,
        diff_data,
        copy(erp1.layout),
        erp1.sample_rate,
        copy(erp1.analysis_info),
        min_epochs,
    )
end

"""
Process a single ERP file through difference wave creation.
Returns BatchResult with success/failure info.
"""
function _condition_difference_process_file(
    filepath::String,
    output_path::String,
    condition_pairs::Union{Vector{Tuple{Int,Int}},Vector{Vector{Int}}},
)
    filename = basename(filepath)

    # Read data (using read_data which finds by type)
    file_data = read_data(filepath)

    if isnothing(file_data)
        return BatchResult(false, filename, "No data variables found")
    end

    # Dispatch to the appropriate in-memory condition_difference based on data type
    if file_data isa Vector{<:ErpData} || file_data isa Vector{<:TimeFreqData}
        try
            difference_waves = condition_difference(file_data, condition_pairs)
            jldsave(output_path; data = difference_waves)
            return BatchResult(true, filename, "Created $(length(difference_waves)) difference(s)")
        catch e
            return BatchResult(false, filename, "$(sprint(showerror, e))")
        end
    else
        return BatchResult(false, filename, "Invalid data type: expected Vector{ErpData} or Vector{TimeFreqData}, got $(typeof(file_data))")
    end
end

# === IN-MEMORY API FUNCTION ===

"""
    condition_difference(data::Vector{<:ErpData}, condition_pairs) -> Vector{ErpData}
    condition_difference(data::Vector{TimeFreqData}, condition_pairs) -> Vector{TimeFreqData}
    condition_difference(file_pattern::String, condition_pairs;
    input_dir::String = pwd(), participant_selection::Function = participants(),
    output_dir = nothing)

Create condition difference waves (or time-frequency differences), or batch-process JLD2 files.

For each pair `(a, b)`, subtracts condition `b` from condition `a`. `condition_pairs` can be
`Vector{Tuple{Int,Int}}` or `Vector{Vector{Int}}`. The `file_pattern` method loads JLD2 files.

# Example
```julia
erps = average_epochs(epochs)
diff_waves = condition_difference(erps, [(1, 2)])

# Batch
condition_difference("erps_unrejected", [(1, 2)])
```
"""
function condition_difference(data::Vector{<:ErpData}, condition_pairs::Union{Vector{Tuple{Int,Int}},Vector{Vector{Int}}})::Vector{ErpData}

    @info "Creating condition difference waves..."

    difference_waves = ErpData[]

    for (pair_idx, pair) in enumerate(condition_pairs)
        cond1, cond2 = pair[1], pair[2]

        # Find ERP data for each condition
        erp1 = nothing
        erp2 = nothing
        for erp in data
            cond_num = erp.condition
            if cond_num == cond1
                erp1 = erp
            elseif cond_num == cond2
                erp2 = erp
            end
        end

        if isnothing(erp1)
            @minimal_warning "Condition $cond1 not found. Skipping pair ($cond1, $cond2)."
            continue
        end
        if isnothing(erp2)
            @minimal_warning "Condition $cond2 not found. Skipping pair ($cond1, $cond2)."
            continue
        end

        diff_wave = _create_difference_wave(erp1, erp2, cond1, cond2, pair_idx)
        push!(difference_waves, diff_wave)
    end

    if isempty(difference_waves)
        @minimal_error("No valid condition pairs found")
    end

    @info "Created $(length(difference_waves)) difference wave(s)"
    return difference_waves
end

# === BATCH API FUNCTION ===

function condition_difference(
    file_pattern::String,
    condition_pairs::Union{Vector{Tuple{Int,Int}},Vector{Vector{Int}}};
    input_dir::String = pwd(),
    participant_selection::Function = participants(),
    output_dir::Union{String,Nothing} = nothing,
)

    # Setup logging
    log_file = "condition_difference.log"
    setup_global_logging(log_file)

    result = (success = 0, errors = 0)  # Default return value
    try
        @info "Batch condition differencing started at $(now())"
        @log_call "condition_difference"

        # Validation (early return on error)
        error_msg = _validate_input_dir(input_dir)
        if !isnothing(error_msg)
            @minimal_error(error_msg)
        end



        error_msg = _validate_condition_pairs(condition_pairs)
        if !isnothing(error_msg)
            @minimal_error(error_msg)
        end

        # Setup directories
        output_dir = something(output_dir, _condition_difference_default_output_dir(input_dir, file_pattern, condition_pairs))
        mkpath(output_dir)

        # Find files
        files = _find_batch_files(file_pattern, input_dir, participant_selection)

        if isempty(files)
            @minimal_warning "No JLD2 files found matching pattern '$file_pattern' in $input_dir"
            result = (success = 0, errors = 0)
        else
            @info "Found $(length(files)) JLD2 files matching pattern '$file_pattern'"
            @info "Condition pairs: $condition_pairs"

            # Create processing function with captured parameters
            process_fn = (input_path, output_path) -> _condition_difference_process_file(input_path, output_path, condition_pairs)

            # Execute batch operation
            results = _run_batch_operation(process_fn, files, input_dir, output_dir; operation_name = "Creating differences")

            result = _log_batch_summary(results, output_dir)
        end

    finally
        _cleanup_logging(log_file, output_dir)
    end

    return result
end
