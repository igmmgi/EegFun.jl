"""
    _apply_rereference!(dat::DataFrame, channel_selection::Vector{Symbol}, reference_selection::Vector{Symbol})

Internal function that applies rereferencing to specified channels in a DataFrame.

# Arguments
- `dat::DataFrame`: The data to rereference
- `channel_selection::Vector{Symbol}`: Names of channels to rereference
- `reference_selection::Vector{Symbol}`: Names of channels to use for reference calculation
"""
function _apply_rereference!(dat::DataFrame, channel_selection::Vector{Symbol}, reference_selection::Vector{Symbol})
    reference = _calculate_reference(dat, reference_selection)
    @views dat[!, channel_selection] .-= reference
    return nothing
end

"""
    _apply_rereference!(dat::Vector{DataFrame}, channel_selection::Vector{Symbol}, reference_selection::Vector{Symbol})

Internal function that applies rereferencing to specified channels in a vector of DataFrames.
"""
function _apply_rereference!(dat::Vector{DataFrame}, channel_selection::Vector{Symbol}, reference_selection::Vector{Symbol})
    _apply_rereference!.(dat, Ref(channel_selection), Ref(reference_selection))
    return nothing
end

"""
    _calculate_reference(dat::DataFrame, reference_channels)

Calculate reference signal from specified channels.

# Arguments
- `dat::DataFrame`: The data containing reference channels
- `reference_channels`: Channel names or indices to use for reference calculation

# Returns
- Vector containing the average of specified reference channels
"""
function _calculate_reference(dat::DataFrame, reference_channels)
    reference = zeros(n_samples(dat))
    @inbounds for channel in reference_channels
        @views reference .+= dat[!, channel]
    end
    return reference ./ length(reference_channels)
end


"""
    rereference!(dat::Union{ContinuousData,ErpData,EpochData}, reference_channel; channel_selection::Function = channels())

Apply rereferencing to EEG data types using predicate-based channel selection.

# Arguments
- `dat::Union{ContinuousData,ErpData,EpochData}`: The EEG data to rereference
- `reference_channel`: Channels to use as reference, can be:
    - Channel names (as symbols): `[:M1, :M2]`, `[:Cz]`, etc.
    - Special symbols: `:avg` (average reference) or `:mastoid` (M1+M2)
- `channel_selection::Function`: Channel selection predicate (default: `channels()` - all EEG channels)

# Effects
- For ContinuousData/ErpData: Modifies data in-place
- For EpochData: Modifies each epoch in-place
- Applies to channels selected by the predicate
- If a channel is included in both reference and rereferenced set, it will become zero

# Notes
- The reference signal is calculated per epoch for EpochData to maintain proper signal processing
- Reference channels must exist in the EEG data layout
- Uses efficient pre-allocated vectors and @views for better performance
- For EpochData, the same reference channels are used across all epochs, but the reference signal is calculated from each epoch's data
"""
# helper function to handle special reference cases such as :avg and :mastoid
function _get_reference_channels(dat, reference_channel::Vector{Symbol})
    if reference_channel[1] == :none
        return Symbol[]  # No reference channels for :none
    elseif reference_channel[1] == :avg # all channels
        return channel_labels(dat)
    elseif reference_channel[1] == :mastoid
        return [:M1, :M2]
    end
    return reference_channel
end

function _get_reference_channels(dat::EegData, reference_channel::Symbol)
    return _get_reference_channels(dat, [reference_channel])
end

# Single method for all EEG data types
function rereference!(dat::EegData, reference_selection::Union{Symbol,Vector{Symbol}}, channel_selection::Function = channels())

    reference_channels = _get_reference_channels(dat, reference_selection)

    # If no reference channels (e.g., :none), return early without rereferencing
    if isempty(reference_channels)
        @info "No rereferencing applied (reference: $(reference_selection))"
        return
    end

    selected_channels = get_selected_channels(dat, channel_selection, include_meta = false, include_extra = false)

    # Verify reference channels exist in the data
    missing_channels = [ch for ch in reference_channels if ch ∉ channel_labels(dat)]
    if !isempty(missing_channels)
        @minimal_error_throw "Missing reference channels in data: $(missing_channels)"
    end

    # Calculate reference signal and apply rereferencing
    @info "Rereference channels: $(reference_selection) ($(print_vector(reference_channels; n_ends=3)))"
    _apply_rereference!(dat.data, selected_channels, reference_channels)

    # Store reference info
    dat.analysis_info.reference = reference_selection isa Symbol ? reference_selection : Symbol(join(reference_selection, '_'))

    return nothing

end

"""
    rereference!(dat::Vector{EpochData}, reference_selection; channel_selection=channels())

Apply rereferencing in-place to a vector of EpochData objects.

# Arguments
- `dat::Vector{EpochData}`: Vector of epoch data to rereference
- `reference_selection::Union{Symbol, Vector{Symbol}}`: Reference channels (`:avg`, `:mastoid`, or specific channels)
- `channel_selection::Function`: Channel selection predicate (default: channels() - all channels)

# Notes
- Modifies each EpochData in the vector in-place by applying rereferencing
- Reference channels can be special symbols (`:avg`, `:mastoid`) or specific channel names
- Reference signal is calculated per epoch to maintain proper signal processing
"""
function rereference!(dat::Vector{EpochData}, reference_selection::Union{Symbol,Vector{Symbol}}, channel_selection::Function = channels())
    rereference!.(dat, Ref(reference_selection), channel_selection)
    return nothing
end

"""
    rereference!(dat::Vector{ErpData}, reference_selection; channel_selection=channels())

Apply rereferencing in-place to a vector of ErpData objects.

# Arguments
- `dat::Vector{ErpData}`: Vector of ERP data to rereference
- `reference_selection::Union{Symbol, Vector{Symbol}}`: Reference channels (`:avg`, `:mastoid`, or specific channels)
- `channel_selection::Function`: Channel selection predicate (default: channels() - all channels)

# Notes
- Modifies each ErpData in the vector in-place by applying rereferencing
- Reference channels can be special symbols (`:avg`, `:mastoid`) or specific channel names
"""
function rereference!(dat::Vector{ErpData}, reference_selection::Union{Symbol,Vector{Symbol}}, channel_selection::Function = channels())
    rereference!.(dat, Ref(reference_selection), channel_selection)
    return nothing
end

# generates all non-mutating versions
@add_nonmutating rereference!


"""
Batch rereferencing for EEG/ERP data.
"""

#=============================================================================
    REREFERENCE-SPECIFIC VALIDATION
=============================================================================#

"""Generate default output directory name for rereferencing operation."""
function _default_rereference_output_dir(input_dir::String, pattern::String, ref_selection)
    ref_str = ref_selection isa Symbol ? string(ref_selection) : join(ref_selection, "_")
    joinpath(input_dir, "rereferenced_$(pattern)_$(ref_str)")
end

#=============================================================================
    REREFERENCE-SPECIFIC PROCESSING
=============================================================================#

"""
Process a single file through rereferencing pipeline.
Returns BatchResult with success/failure info.
"""
function _process_rereference_file(filepath::String, output_path::String, reference_selection, condition_selection::Function)
    filename = basename(filepath)

    # Load data
    data = read_data(filepath)
    if isnothing(data)
        return BatchResult(false, filename, "No data variables found")
    end

    # Validate that data is valid EEG data (Vector of ErpData or EpochData)
    if !(data isa Vector{<:Union{ErpData,EpochData}})
        return BatchResult(false, filename, "Invalid data type: expected Vector{ErpData} or Vector{EpochData}")
    end

    # Select conditions
    data = _condition_select(data, condition_selection)

    # Handle empty data (valid case - just save empty vector)
    if isempty(data)
        jldsave(output_path; data = data)
        ref_str = reference_selection isa Symbol ? string(reference_selection) : join(reference_selection, ", ")
        return BatchResult(true, filename, "Rereferenced to $ref_str (empty data)")
    end

    # Apply rereferencing (mutates data in-place)
    rereference!.(data, reference_selection)

    # Save (always use "data" as variable name since load_data finds by type)
    jldsave(output_path; data = data)

    ref_str = reference_selection isa Symbol ? string(reference_selection) : join(reference_selection, ", ")
    return BatchResult(true, filename, "Rereferenced to $ref_str")
end

#=============================================================================
    MAIN API FUNCTION
=============================================================================#

"""
    rereference(file_pattern::String; 
               input_dir::String = pwd(), 
               reference_selection::Union{Symbol, Vector{Symbol}} = :avg, 
               participant_selection::Function = participants(),
               condition_selection::Function = conditions(),
               output_dir::Union{String, Nothing} = nothing)

Rereference EEG/ERP data from JLD2 files and save to a new directory.

# Arguments
- `file_pattern::String`: Pattern to match files ("epochs", "erps", "cleaned", "original", or custom)
- `input_dir::String`: Input directory containing JLD2 files (default: current directory)
- `reference_selection::Union{Symbol, Vector{Symbol}}`: Reference channels, can be:
  - Special symbols: `:avg` (average reference), `:mastoid` (M1+M2)
  - Channel names: `[:Cz]`, `[:M1, :M2]`, etc.
- `participant_selection::Function`: Participant selection predicate (default: `participants()` for all)
- `condition_selection::Function`: Condition selection predicate (default: `conditions()` for all)
- `output_dir::Union{String, Nothing}`: Output directory (default: creates subdirectory based on reference settings)

# Example
```julia
# Rereference all epochs to average reference
rereference("epochs")

# Rereference specific participant to mastoid reference
rereference("epochs", reference_selection=:mastoid, participants=3)

# Rereference to Cz for specific participants and conditions
rereference("epochs", reference_selection=[:Cz], participants=[3, 4], conditions=[1, 2])
```
"""
function rereference(
    file_pattern::String;
    input_dir::String = pwd(),
    reference_selection::Union{Symbol,Vector{Symbol}} = :avg,
    participant_selection::Function = participants(),
    condition_selection::Function = conditions(),
    output_dir::Union{String,Nothing} = nothing,
)

    # Setup logging
    log_file = "rereference.log"
    setup_global_logging(log_file)

    result = (success = 0, errors = 0)  # Default return value
    try
        @info "Batch rereferencing started at $(now())"
        @log_call "rereference"

        # Validation
        if (error_msg = _validate_input_dir(input_dir)) !== nothing
            @minimal_error_throw(error_msg)
        end

        # Setup directories
        output_dir = something(output_dir, _default_rereference_output_dir(input_dir, file_pattern, reference_selection))
        mkpath(output_dir)

        # Find files
        files = _find_batch_files(file_pattern, input_dir, participant_selection)

        if isempty(files)
            @minimal_warning "No JLD2 files found matching pattern '$file_pattern' in $input_dir"
            result = (success = 0, errors = 0)
        else
            ref_str = reference_selection isa Symbol ? string(reference_selection) : join(reference_selection, ", ")
            @info "Found $(length(files)) JLD2 files matching pattern '$file_pattern'"
            @info "Reference settings: $ref_str"

            # Create processing function with captured parameters
            process_fn =
                (input_path, output_path) -> _process_rereference_file(input_path, output_path, reference_selection, condition_selection)

            # Execute batch operation
            results = _run_batch_operation(process_fn, files, input_dir, output_dir; operation_name = "Rereferencing")

            result = _log_batch_summary(results, output_dir)
        end

    finally
        _cleanup_logging(log_file, output_dir)
    end

    return result
end
