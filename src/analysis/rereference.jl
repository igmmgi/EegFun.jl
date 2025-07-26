"""
    _apply_rereference!(dat::DataFrame, reference::Vector{<:Real}, channel_selection::Vector{Symbol})

Internal function that applies rereferencing to specified channels in a DataFrame.

# Arguments
- `dat::DataFrame`: The data to rereference
- `reference::Vector{<:Real}`: Reference signal to subtract from each channel
- `channel_selection::Vector{Symbol}`: Names of channels to rereference
"""
function _apply_rereference!(
    dat::DataFrame, 
    reference::Vector{Float64},
    channel_selection::Vector{Symbol}
)
    @inbounds for col in propertynames(dat)
        if col in channel_selection
            @views dat[!, col] .-= reference
        end
    end
end

"""
    _apply_rereference!(dat::Vector{DataFrame}, reference::Vector{<:Real}, channel_selection::Vector{Symbol})

Internal function that applies rereferencing to specified channels in a vector of DataFrames.
"""
function _apply_rereference!(
    dat::Vector{DataFrame}, 
    reference::Vector{Float64},
    channel_selection::Vector{Symbol}
)
    _apply_rereference!.(dat, Ref(reference), Ref(channel_selection))
    return nothing
end

"""
    calculate_reference(dat::DataFrame, reference_channels)

Calculate reference signal from specified channels.

# Arguments
- `dat::DataFrame`: The data containing reference channels
- `reference_channels`: Channel names or indices to use for reference calculation

# Returns
- Vector containing the average of specified reference channels
"""
function calculate_reference(dat::DataFrame, reference_channels)
    reference = zeros(n_samples(dat))
    @inbounds for channel in reference_channels
        @views reference .+= dat[!, channel]
    end
    return reference ./ length(reference_channels)
end

"""
    rereference!(dat::DataFrame, reference_channels, channel_selection::Vector{Symbol})

Apply rereferencing to specified channels in a DataFrame.

# Arguments
- `dat::DataFrame`: The data to rereference (should contain only EEG channels)
- `reference_channels`: Channels to use as reference, can be:
    - Channel names (as symbols): `[:M1, :M2]`, `[:Cz]`, etc.
    - Special symbols: `:avg` (average reference) or `:mastoid` (M1+M2)
- `channel_selection::Vector{Symbol}`: Names of channels to rereference

# Effects
- Modifies input data in-place by subtracting reference signal from specified channels
- If a channel is included in both reference and rereferenced set, it will become zero (e.g., Cz when referencing to Cz)
"""
function rereference!(
    dat::DataFrame, 
    reference_channels::Vector{Symbol},
    channel_selection::Vector{Symbol}
)
    @info "Rereferencing channels $(channel_selection) using reference: $(_print_vector(reference_channels))"
    reference = calculate_reference(dat, reference_channels)
    _apply_rereference!(dat, reference, channel_selection)
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
function _resolve_reference(dat, reference_channel::Vector{Symbol})
    if reference_channel[1] == :avg # all channels
        return channels(dat) 
    elseif reference_channel[1] == :mastoid
        return [:M1, :M2]
    end
    return reference_channel
end

function _resolve_reference(dat::EegData, reference_channel::Symbol)
    return _resolve_reference(dat, [reference_channel])
end

# Single method for all EEG data types
function rereference!(
    dat::EegData,
    reference_channel::Union{Symbol,Vector{Symbol}},
    channel_selection::Function = channels()
)
    ref_channels = _resolve_reference(dat, reference_channel)
    selected_channels = get_selected_channels(dat, channel_selection)
    
    # Verify reference channels exist in the data
    missing_channels = [ch for ch in ref_channels if ch ∉ channels(dat)]
    if !isempty(missing_channels)
        @minimal_error "Missing reference channels in data: $(missing_channels)"
    end
    
    @info "Rereferencing channels $(selected_channels) using: $(ref_channels)"
    
    # Calculate reference signal and apply rereferencing
    reference = calculate_reference(dat.data isa Vector ? dat.data[1] : dat.data, ref_channels)
    _apply_rereference!(dat.data, reference, selected_channels)
    
    # Store reference info
    dat.analysis_info.reference = reference_channel isa Symbol ? reference_channel : Symbol(join(reference_channel, '_'))
    return nothing
end

# generates all non-mutating versions
@add_nonmutating rereference!
