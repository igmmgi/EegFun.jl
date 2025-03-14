const ChannelLabel = Union{Symbol, String}

"""
    _apply_rereference!(dat::DataFrame, channel_labels, reference::Vector{<:Real})

Internal function that applies rereferencing to specified channels in a DataFrame.

# Arguments
- `dat::DataFrame`: The data to rereference
- `channel_labels::Vector{Symbol}`: Names of channels to rereference
- `reference::Vector{<:Real}`: Reference signal to subtract from each channel

# Effects
- Modifies the input data in-place by subtracting the reference signal from specified channels
"""
function _apply_rereference!(
    dat::DataFrame, 
    channel_labels::Vector{Symbol}, 
    reference::Vector{Float64}
)
    for col in propertynames(dat)
        if col in channel_labels
            dat[!, col] .-= reference
        end
    end
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
    return reduce(+, eachcol(dat[:, reference_channels])) ./ length(reference_channels)
end

"""
    rereference!(dat::DataFrame, channel_labels, reference_channels)

Apply rereferencing to specified channels in a DataFrame.

# Arguments
- `dat::DataFrame`: The data to rereference
- `channel_labels::Vector{Union{Symbol, String}}`: Names of channels to rereference
- `reference_channels`: Channels to use as reference, can be:
    - Channel names (as strings or symbols)
    - Special strings: "avg" (average reference) or "mastoid" (M1+M2)

# Effects
- Modifies input data in-place by subtracting reference signal from specified channels
"""
function rereference!(
    dat::DataFrame, 
    channel_labels::Vector{Symbol}, 
    reference_channels::Vector{Symbol}
)
    @info "Rereferencing using channel(s): $(print_vector_(reference_channels))"
    reference = calculate_reference(dat, reference_channels)
    _apply_rereference!(dat, channel_labels, reference)
end



"""
    rereference!(dat::Union{ContinuousData,ErpData,EpochData}, channel_labels, reference_channel)
    rereference!(dat::Union{ContinuousData,ErpData,EpochData}, reference_channel)

Apply rereferencing to EEG data types.

# Arguments
- `dat::Union{ContinuousData,ErpData,EpochData}`: The data to rereference
- `channel_labels::Vector{Union{Symbol, String}}`: Names of channels to rereference (optional)
- `reference_channel`: Channels to use as reference (see DataFrame method for options)

# Effects
- For ContinuousData/ErpData: Modifies data in-place
- For EpochData: Modifies each epoch in-place
- If channel_labels omitted, uses all channels from layout
"""

# helper function to handle special reference cases such as :avg and :mastoid
function resolve_reference(dat, reference_channel::Symbol)
    if reference_channel == :avg # all channels
        return channels(dat) 
    elseif reference_channel == :mastoid
        return [:M1, :M2]
    end
end

# Base methods for Vector{Symbol} reference
function rereference!(
    dat::Union{ContinuousData,ErpData},
    channel_labels::Vector{Symbol},
    reference_channel::Vector{Symbol},
)
    rereference!(dat.data, channel_labels, reference_channel)
    return nothing
end

function rereference!(dat::EpochData, channel_labels::Vector{Symbol}, reference_channel::Vector{Symbol})
    for epoch in dat.data
        rereference!(epoch, channel_labels, reference_channel)
    end
    return nothing
end

# Methods for Symbol reference (converts to Vector{Symbol})
function rereference!(
    dat::Union{ContinuousData,ErpData,EpochData},
    channel_labels::Vector{Symbol},
    reference_channel::Symbol,
)
    rereference!(dat, channel_labels, resolve_reference(dat, reference_channel))
    return nothing
end

# Methods for default channel labels
function rereference!(dat::Union{ContinuousData,ErpData,EpochData}, reference_channel::Union{Symbol,Vector{Symbol}})
    rereference!(dat, channels(dat), reference_channel)
    return nothing
end

# generates all non-mutating versions
@add_nonmutating rereference!





