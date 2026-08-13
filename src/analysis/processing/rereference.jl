"""
    _apply_rereference!(dat::DataFrame, channel_selection::Vector{Symbol}, reference_selection::Vector{Symbol})

Internal function that applies rereferencing to specified channels in a DataFrame.

# Arguments
- `dat::DataFrame`: The data to rereference
- `channel_selection::Vector{Symbol}`: Names of channels to rereference
- `reference_selection::Vector{Symbol}`: Names of channels to use for reference calculation
"""
function _apply_rereference!(dat::DataFrame, channel_selection::Vector{Symbol}, reference_selection::Vector{Symbol})
    reference = zeros(n_samples(dat))
    _calculate_reference!(reference, dat, reference_selection)
    n = length(reference)
    for ch in channel_selection
        col = dat[!, ch]
        if col isa Vector{Float64}
            @inbounds @simd for i = 1:n
                col[i] -= reference[i]
            end
        else
            dat[!, ch] .-= reference
        end
    end
    return nothing
end

"""
    _apply_rereference!(dat::Vector{DataFrame}, channel_selection::Vector{Symbol}, reference_selection::Vector{Symbol})

Internal function that applies rereferencing to specified channels in a vector of DataFrames.
"""
function _apply_rereference!(dat::Vector{DataFrame}, channel_selection::Vector{Symbol}, reference_selection::Vector{Symbol})
    if isempty(dat)
        return nothing
    end

    n = n_samples(dat[1])

    Threads.@threads for df in dat
        reference = Vector{Float64}(undef, n)
        _calculate_reference!(reference, df, reference_selection)
        for ch in channel_selection
            col = df[!, ch]
            if col isa Vector{Float64}
                @inbounds @simd for i = 1:n
                    col[i] -= reference[i]
                end
            else
                df[!, ch] .-= reference
            end
        end
    end
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
function _calculate_reference!(reference::Vector{Float64}, dat::DataFrame, reference_channels)
    fill!(reference, 0.0)
    n = length(reference)
    @inbounds for channel in reference_channels
        col = dat[!, channel]
        if col isa Vector{Float64}
            @simd for i = 1:n
                reference[i] += col[i]
            end
        else
            reference .+= col
        end
    end

    n_refs = Float64(length(reference_channels))
    @inbounds @simd for i = 1:n
        reference[i] /= n_refs
    end
    return reference
end


"""
    rereference!(dat::Union{ContinuousData,ErpData,EpochData}, reference_selection; channel_selection = channels())
    rereference!(dat::Vector{EpochData}, reference_selection; channel_selection = channels())
    rereference!(dat::Vector{ErpData}, reference_selection; channel_selection = channels())

Rereference EEG data in-place. `reference_selection` can be `:avg`, `:mastoid`, or
a channel name / vector of channel names. The `Vector` form broadcasts across conditions.

# Examples
```julia
rereference!(dat, :avg)
rereference!(dat, :mastoid)
rereference!(epochs, [:M1, :M2])
```
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
"""
    rereference!(dat::EegData, reference_selection::Union{Symbol,Vector{Symbol}}; channel_selection::Function = channels())

Rereference EEG data to a specified channel or combination of channels.
"""
function rereference!(dat::EegData, reference_selection::Union{Symbol,Vector{Symbol}}; channel_selection::Function = channels())

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
        @minimal_error "Missing reference channels in data: $(missing_channels)"
    end

    # Calculate reference signal and apply rereferencing
    @info "Rereference channels: $(reference_selection) ($(_print_vector(reference_channels; n_ends=3)))"
    _apply_rereference!(dat.data, selected_channels, reference_channels)

    # Store reference info
    dat.analysis_info.reference = reference_selection isa Symbol ? reference_selection : Symbol(join(reference_selection, '_'))


    return nothing

end


function rereference!(dat::Vector{EpochData}, reference_selection::Union{Symbol,Vector{Symbol}}; channel_selection::Function = channels())
    rereference!.(dat, Ref(reference_selection); channel_selection = channel_selection)
    return nothing
end


function rereference!(dat::Vector{ErpData}, reference_selection::Union{Symbol,Vector{Symbol}}; channel_selection::Function = channels())
    rereference!.(dat, Ref(reference_selection); channel_selection = channel_selection)
    return nothing
end

# generates all non-mutating versions
"""
    rereference(dat, args...; kwargs...)

Non-mutating version of `rereference!`.
"""
function rereference(dat, args...; kwargs...)
    dat_copy = copy(dat)
    rereference!(dat_copy, args...; kwargs...)
    return dat_copy
end



# === REREFERENCE-SPECIFIC VALIDATION ===

"""Generate default output directory name for rereferencing operation."""
function _default_rereference_output_dir(input_dir::String, pattern::String, ref_selection)
    ref_str = ref_selection isa Symbol ? string(ref_selection) : join(ref_selection, "_")
    joinpath(input_dir, "rereferenced_$(clean_pattern(pattern))_$(ref_str)")
end

# === REREFERENCE-SPECIFIC PROCESSING ===

"""
Process a single file through rereferencing pipeline.
Returns BatchResult with success/failure info.
"""
function _process_rereference_file(filepath::String, output_path::String, reference_selection, condition_selection::Function)
    fname = basename(filepath)

    # Read data
    data = read_data(filepath)
    if isnothing(data)
        return BatchResult(false, fname, "No data variables found")
    end

    # Validate that data is valid EEG data (Vector of ErpData or EpochData)
    if !(data isa Vector{<:Union{ErpData,EpochData}})
        return BatchResult(false, fname, "Invalid data type: expected Vector{ErpData} or Vector{EpochData}")
    end

    # Select conditions
    data = _condition_select(data, condition_selection)

    # Handle empty data (valid case - just save empty vector)
    if isempty(data)
        jldsave(output_path; data = data)
        ref_str = reference_selection isa Symbol ? string(reference_selection) : join(reference_selection, ", ")
        return BatchResult(true, fname, "Rereferenced to $ref_str (empty data)")
    end

    # Apply rereferencing (mutates data in-place)
    rereference!.(data, Ref(reference_selection))

    # Save (always use "data" as variable name since read_data finds by type)
    jldsave(output_path; data = data)

    ref_str = reference_selection isa Symbol ? string(reference_selection) : join(reference_selection, ", ")
    return BatchResult(true, fname, "Rereferenced to $ref_str")
end

# === MAIN API FUNCTION ===

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
        error_msg = _validate_input_dir(input_dir)
        if !isnothing(error_msg)
            @minimal_error(error_msg)
        end

        # Setup directories
        output_dir = something(output_dir, _default_rereference_output_dir(input_dir, file_pattern, reference_selection))
        mkpath(output_dir)

        ref_str = reference_selection isa Symbol ? string(reference_selection) : join(reference_selection, ", ")
        @info "Reference settings: $ref_str"

        process_fn =
            (input_path, output_path) -> _process_rereference_file(input_path, output_path, reference_selection, condition_selection)

        result = batch_process(process_fn, file_pattern, input_dir, output_dir, participant_selection, "Rereferencing")

    finally
        _cleanup_logging(log_file, output_dir)
    end

    return result
end
