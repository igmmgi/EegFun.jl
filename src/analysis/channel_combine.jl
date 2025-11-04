"""
Batch channel combining for EEG/ERP data.
"""

#=============================================================================
    COMBINE-CHANNELS-SPECIFIC HELPERS
=============================================================================#

"""Generate default output directory name for channel combining."""
function _channel_combine_default_output_dir(input_dir::String, pattern::String)
    joinpath(input_dir, "combined_channels_$(pattern)")
end

#=============================================================================
    COMBINE-CHANNELS-SPECIFIC PROCESSING
=============================================================================#

"""
Process a single file through channel combining pipeline.
Returns BatchResult with success/failure info.
"""
function _channel_combine_process_file(
    filepath::String,
    output_path::String,
    channel_selections::Vector{<:Function},
    output_labels::Vector{Symbol},
    conditions,
    reduce::Bool,
)
    filename = basename(filepath)

    # Load data
    data = load_data(filepath)
    if isnothing(data)
        return BatchResult(false, filename, "No recognized data variable")
    end

    # Select conditions
    data = _condition_select(data, conditions)

    # Apply channel combination to each data item
    foreach(data) do item
        channel_average!(item, channel_selections = channel_selections, output_labels = output_labels, reduce = reduce)
    end

    # Determine variable name based on data type
    var_name = data isa Vector{<:ErpData} ? "erps" : "epochs"

    # Save
    save(output_path, var_name, data)

    n_groups = length(channel_selections)
    return BatchResult(true, filename, "Combined $n_groups channel group(s)")
end

#=============================================================================
    MAIN API FUNCTION
=============================================================================#

"""
    channel_combine(file_pattern::String, channel_selections::Vector{Function}; 
                    output_labels::Union{Vector{Symbol}, Nothing} = nothing,
                    input_dir::String = pwd(), 
                    participants::Union{Int, Vector{Int}, Nothing} = nothing,
                    conditions::Union{Int, Vector{Int}, Nothing} = nothing,
                    output_dir::Union{String, Nothing} = nothing,
                    reduce::Bool = false)

Batch process EEG data files to combine specified channels using predicates.

This function loads JLD2 files containing EEG data, combines specified channel groups by averaging,
and saves the resulting data with new combined channels to a new directory.

# Arguments
- `file_pattern::String`: Pattern to match JLD2 files (e.g., "erps_cleaned", "epochs_cleaned")
- `channel_selections::Vector{Function}`: Channel selection predicates (e.g., `[channels([:Fp1, :Fp2]), channels([:PO7, :PO8])]`)
- `output_labels::Union{Vector{Symbol}, Nothing}`: Labels for combined channels (default: auto-generated from selection)
- `input_dir::String`: Input directory containing JLD2 files (default: current directory)
- `participants::Union{Int, Vector{Int}, Nothing}`: Participant numbers to process (default: all)
- `conditions::Union{Int, Vector{Int}, Nothing}`: Condition numbers to process (default: all)
- `output_dir::Union{String, Nothing}`: Output directory (default: auto-generated)
- `reduce::Bool`: Whether to keep only combined channels (true) or append to existing (false, default)

# Examples
```julia
# Combine frontal and parietal channels using predicates
channel_combine("erps_cleaned", [channels([:Fp1, :Fp2]), channels([:PO7, :PO8])])

# With custom output labels
channel_combine("erps_cleaned", 
                [channels([:Fp1, :Fp2]), channels([:PO7, :PO8])],
                output_labels = [:frontal, :parietal])

# Process specific participants and conditions
channel_combine("erps_cleaned", [channels([:Fp1, :Fp2])], 
                input_dir = "/path/to/data", 
                participants = [1, 2, 3], 
                conditions = [1, 2])

# Create reduced dataset with only combined channels
channel_combine("erps_cleaned", 
                [channels([:Fp1, :Fp2]), channels([:PO7, :PO8])], 
                reduce = true)
```
"""
function channel_combine(
    file_pattern::String,
    channel_selections::Vector{<:Function};
    output_labels::Union{Vector{Symbol},Nothing} = nothing,
    input_dir::String = pwd(),
    participants::Union{Int,Vector{Int},Nothing} = nothing,
    conditions::Union{Int,Vector{Int},Nothing} = nothing,
    output_dir::Union{String,Nothing} = nothing,
    reduce::Bool = false,
)

    # Setup logging
    log_file = "channel_combine.log"
    setup_global_logging(log_file)

    try
        @info "Batch channel combining started at $(now())"
        @log_call "channel_combine" (file_pattern, channel_selections)

        # Validation (early return on error)
        if (error_msg = _validate_input_dir(input_dir)) !== nothing
            @minimal_error_throw(error_msg)
        end

        # Generate output labels if not provided
        if isnothing(output_labels)
            output_labels = [Symbol("combined_$i") for i = 1:length(channel_selections)]
        end

        # Validate output labels length matches channel selections
        if length(output_labels) != length(channel_selections)
            @minimal_error_throw(
                "Number of output_labels ($(length(output_labels))) must match number of channel_selections ($(length(channel_selections)))"
            )
        end

        # Setup directories
        output_dir = something(output_dir, _channel_combine_default_output_dir(input_dir, file_pattern))
        mkpath(output_dir)

        # Find files
        files = _find_batch_files(file_pattern, input_dir; participants)

        if isempty(files)
            @minimal_warning "No JLD2 files found matching pattern '$file_pattern' in $input_dir"
            return nothing
        end

        @info "Found $(length(files)) JLD2 files matching pattern '$file_pattern'"
        @info "Channel selections: $(length(channel_selections)) group(s)"
        @info "Output labels: $output_labels"
        @info "Reduce mode: $reduce"

        # Create processing function with captured parameters
        process_fn =
            (input_path, output_path) -> _channel_combine_process_file(
                input_path,
                output_path,
                channel_selections,
                output_labels,
                conditions,
                reduce,
            )

        # Execute batch operation
        results = _run_batch_operation(process_fn, files, input_dir, output_dir; operation_name = "Combining channels")

        _log_batch_summary(results, output_dir)

    finally
        _cleanup_logging(log_file, output_dir)
    end
end
