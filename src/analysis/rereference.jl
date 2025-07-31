"""
    _apply_rereference!(dat::DataFrame, channel_selection::Vector{Symbol}, reference_selection::Vector{Symbol})

Internal function that applies rereferencing to specified channels in a DataFrame.

# Arguments
- `dat::DataFrame`: The data to rereference
- `channel_selection::Vector{Symbol}`: Names of channels to rereference
- `reference_selection::Vector{Symbol}`: Names of channels to use for reference calculation
"""
function _apply_rereference!( dat::DataFrame, channel_selection::Vector{Symbol}, reference_selection::Vector{Symbol})
    reference = calculate_reference(dat, reference_selection)
    @views dat[!, channel_selection] .-= reference
    return nothing
end

"""
    _apply_rereference!(dat::Vector{DataFrame}, channel_selection::Vector{Symbol}, reference_selection::Vector{Symbol})

Internal function that applies rereferencing to specified channels in a vector of DataFrames.
"""
function _apply_rereference!( dat::Vector{DataFrame}, channel_selection::Vector{Symbol}, reference_selection::Vector{Symbol})
    _apply_rereference!.(dat, Ref(channel_selection), Ref(reference_selection))
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
function get_reference_channels(dat, reference_channel::Vector{Symbol})
    if reference_channel[1] == :avg # all channels
        return channel_labels(dat) 
    elseif reference_channel[1] == :mastoid
        return [:M1, :M2]
    end
    return reference_channel
end

function get_reference_channels(dat::EegData, reference_channel::Symbol)
    return get_reference_channels(dat, [reference_channel])
end

# Single method for all EEG data types
function rereference!(dat::EegData, reference_selection::Union{Symbol,Vector{Symbol}}, channel_selection::Function = channels())

    reference_channels = get_reference_channels(dat, reference_selection)
    selected_channels = get_selected_channels(dat, channel_selection, include_metadata_columns = false)
    
    # Verify reference channels exist in the data
    missing_channels = [ch for ch in reference_channels if ch ∉ channel_labels(dat)]
    if !isempty(missing_channels)
        @minimal_error "Missing reference channels in data: $(missing_channels)"
    end
    
    # Calculate reference signal and apply rereferencing
    @info "Rereferencing channels ($(print_vector(selected_channels))) using: $(reference_selection) ($(print_vector(reference_channels)))"
    _apply_rereference!(dat.data, selected_channels, reference_channels)
    
    # Store reference info
    dat.analysis_info.reference = reference_selection isa Symbol ? reference_selection : Symbol(join(reference_selection, '_'))

    return nothing

end

# generates all non-mutating versions
@add_nonmutating rereference!
