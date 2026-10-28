"""
    _apply_rereference!(dat::DataFrame, channel_selection::Vector{Symbol}, reference_selection::Vector{Symbol})

Internal function that applies rereferencing to specified channels in a DataFrame.

# Arguments
- `dat::DataFrame`: The data to rereference
- `channel_selection::Vector{Symbol}`: Names of channels to rereference
- `reference_selection::Vector{Symbol}`: Names of channels to use for reference calculation
"""
function _apply_rereference!(dat::DataFrame, channel_selection::AbstractVector{Symbol}, reference_selection::AbstractVector{Symbol})
    n = n_samples(dat)
    reference = Vector{Float64}(undef, n)
    ref_indices = [columnindex(dat, ch) for ch in reference_selection]
    ch_indices = [columnindex(dat, ch) for ch in channel_selection]
    _calculate_reference!(reference, dat, ref_indices)
    for col_idx in ch_indices
        col = dat[!, col_idx]
        if col isa Vector{Float64}
            col_f64 = col::Vector{Float64}
            @inbounds @simd for i = 1:n
                col_f64[i] -= reference[i]
            end
        else
            dat[!, col_idx] .-= reference
        end
    end
    return nothing
end

"""
    _apply_rereference!(dat::Vector{DataFrame}, channel_selection::Vector{Symbol}, reference_selection::Vector{Symbol})

Internal function that applies rereferencing to specified channels in a vector of DataFrames.
"""
function _apply_rereference!(
    dat::AbstractVector{DataFrame},
    channel_selection::AbstractVector{Symbol},
    reference_selection::AbstractVector{Symbol},
)
    if isempty(dat)
        return nothing
    end

    first_df = dat[1]
    n = n_samples(first_df)
    ref_indices = [columnindex(first_df, ch) for ch in reference_selection]
    ch_indices = [columnindex(first_df, ch) for ch in channel_selection]

    _foreach_chunk(length(dat)) do df_range
        reference = Vector{Float64}(undef, n)  # owned exclusively by this chunk/task
        for df_idx in df_range
            df = dat[df_idx]
            _calculate_reference!(reference, df, ref_indices)
            for col_idx in ch_indices
                col = df[!, col_idx]
                if col isa Vector{Float64}
                    col_f64 = col::Vector{Float64}
                    @inbounds @simd for i = 1:n
                        col_f64[i] -= reference[i]
                    end
                else
                    df[!, col_idx] .-= reference
                end
            end
        end
    end
    return nothing
end

"""
    _calculate_reference!(reference::Vector{Float64}, dat::DataFrame, reference_channels)

Calculate reference signal from specified channels.

# Arguments
- `dat::DataFrame`: The data containing reference channels
- `reference_channels`: Channel names or indices to use for reference calculation

# Returns
- Vector containing the average of specified reference channels
"""
function _calculate_reference!(reference::Vector{Float64}, dat::DataFrame, ref_indices::AbstractVector{<:Integer})
    fill!(reference, 0.0)
    n = length(reference)
    @inbounds for col_idx in ref_indices
        col = dat[!, col_idx]
        if col isa Vector{Float64}
            col_f64 = col::Vector{Float64}
            @simd for i = 1:n
                reference[i] += col_f64[i]
            end
        else
            reference .+= col
        end
    end

    inv_refs = 1.0 / length(ref_indices)
    @inbounds @simd for i = 1:n
        reference[i] *= inv_refs
    end
    return reference
end

function _calculate_reference!(reference::Vector{Float64}, dat::DataFrame, reference_channels::AbstractVector{Symbol})
    ref_indices = [columnindex(dat, ch) for ch in reference_channels]
    return _calculate_reference!(reference, dat, ref_indices)
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
function _get_reference_channels(dat::EegData, reference_channel::Symbol)
    if reference_channel == :none
        return Symbol[]  # No reference channels for :none
    elseif reference_channel == :avg # all channels
        return channel_labels(dat)
    elseif reference_channel == :mastoid
        return [:M1, :M2]
    end
    return [reference_channel]
end

function _get_reference_channels(dat, reference_channel::AbstractVector{Symbol})
    isempty(reference_channel) && return Symbol[]
    if length(reference_channel) == 1
        return _get_reference_channels(dat, reference_channel[1])
    end
    return reference_channel
end

_get_reference_channels(dat, ref::AbstractString) = _get_reference_channels(dat, Symbol(ref))
_get_reference_channels(dat, refs::AbstractVector{<:AbstractString}) = _get_reference_channels(dat, Symbol.(refs))

# Single method for all EEG data types
"""
    rereference!(dat::EegData, reference_selection; channel_selection = channels())

Rereference EEG data to a specified channel or combination of channels.
"""
function rereference!(dat::EegData, reference_selection; channel_selection = channels())

    reference_channels = _get_reference_channels(dat, reference_selection)

    # If no reference channels (e.g., :none), return early without rereferencing
    if isempty(reference_channels)
        @info "No rereferencing applied (reference: $(reference_selection))"
        return dat
    end

    selected_channels = get_selected_channels(dat, channel_selection, include_meta = false, include_extra = false)

    # Verify reference channels exist in the data
    missing_channels = [ch for ch in reference_channels if ch ∉ dat]
    if !isempty(missing_channels)
        @minimal_error "Missing reference channels in data: $(missing_channels)"
    end

    # Calculate reference signal and apply rereferencing
    @info "Rereference channels: $(reference_selection) ($(_print_vector(reference_channels; n_ends=3)))"
    _apply_rereference!(dat.data, selected_channels, reference_channels)

    # Store reference info
    dat.analysis_info.reference = reference_selection isa Symbol ? reference_selection : Symbol(join(reference_selection, '_'))

    return dat
end


function rereference!(dat::AbstractVector{<:EegData}, reference_selection; channel_selection = channels())
    for d in dat
        rereference!(d, reference_selection; channel_selection = channel_selection)
    end
    return dat
end

# generates all non-mutating versions
"""
    rereference(dat, args...; kwargs...)

Non-mutating version of `rereference!`.
"""
function rereference(dat::EegData, args...; kwargs...)
    return rereference!(copy(dat), args...; kwargs...)
end

function rereference(dat::AbstractVector{<:EegData}, args...; kwargs...)
    return [rereference(d, args...; kwargs...) for d in dat]
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
    if !(data isa AbstractVector{<:Union{ErpData,EpochData}})
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
    reference_selection::Union{Symbol,AbstractVector{Symbol}} = :avg,
    participant_selection = participants(),
    condition_selection = conditions(),
    output_dir::Union{String,Nothing} = nothing,
)

    _with_batch_logging("rereferencing", "rereference.log", input_dir, output_dir) do ctx
        @log_call "rereference"

        # Setup directories
        output_dir = _prepare_output_dir!(ctx, _default_rereference_output_dir(input_dir, file_pattern, reference_selection))

        ref_str = reference_selection isa Symbol ? string(reference_selection) : join(reference_selection, ", ")
        @info "Reference settings: $ref_str"

        process_fn =
            (input_path, output_path) -> _process_rereference_file(input_path, output_path, reference_selection, condition_selection)

        result = batch_process(process_fn, file_pattern, input_dir, output_dir, participant_selection, "Rereferencing")
    end
end
