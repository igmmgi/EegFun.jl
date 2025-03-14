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
    for col in names(dat)
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

# function rereference!(dat::Union{ContinuousData,ErpData}, 
#     channel_labels::Vector{Symbol},
#     reference_channel::Vector{Symbol})
#     rereference!(dat.data, channel_labels, reference_channel)
# end
# 
# function rereference!(
#     dat::EpochData,
#     channel_labels::Vector{Symbol},
#     reference_channel::Vector{Symbol},
# )
#     for epoch in dat.data
#         rereference!(epoch, channel_labels, reference_channel)
#     end
# end
# 
# function rereference!(dat::Union{ContinuousData,ErpData}, reference_channel::Vector{Symbol})
#     rereference!(dat, dat.layout.label, reference_channel)
# end
# 
# function rereference!(dat::EpochData, reference_channel::Vector{Symbol})
#     for epoch in dat.data
#         rereference!(epoch, dat.layout.label, reference_channel)
#     end
# end
# 
# 
# 
# function rereference!(dat::Union{ContinuousData,ErpData}, channel_labels::Vector{Symbol}, reference_channel::Symbol)
#     if reference_channel == :avg
#         reference_channel = dat.layout.label
#     elseif reference_channel == :mastoid
#         reference_channel = [:M1, :M2]
#     end
#     rereference!(dat.data, channel_labels, reference_channel)
# end
# 
# function rereference!(
#     dat::EpochData,
#     channel_labels::Vector{Symbol},
#     reference_channel::Symbol,
# )
#     if reference_channel == :avg
#         reference_channel = dat.layout.label
#     elseif reference_channel == :mastoid
#         reference_channel = [:M1, :M2]
#     end
#     for epoch in dat.data
#         rereference!(epoch, channel_labels, reference_channel)
#     end
# end
# 
# 
# function rereference!(dat::Union{ContinuousData,ErpData}, reference_channel::Symbol)
#     if reference_channel == :avg
#         reference_channel = dat.layout.label
#     elseif reference_channel == :mastoid
#         reference_channel = [:M1, :M2]
#     end
#     rereference!(dat, dat.layout.label, reference_channel)
# end
# 
# function rereference!(dat::EpochData, reference_channel::Symbol)
#     if reference_channel == :avg
#         reference_channel = dat.layout.label
#     elseif reference_channel == :mastoid
#         reference_channel = [:M1, :M2]
#     end
#     for epoch in dat.data
#         rereference!(epoch, dat.layout.label, reference_channel)
#     end
# end


macro add_nonmutating(func)
    # Get the function name without !
    nonmutating_name = Symbol(string(func)[1:end-1])
    
    return quote
        function $(esc(nonmutating_name))(args...)
            data_copy = deepcopy(first(args))
            $(esc(func))(data_copy, Base.tail(args)...)
            return data_copy
        end
    end
end


# Helper function to handle special reference cases
function resolve_reference(dat, reference_channel::Symbol)
    if reference_channel == :avg
        return dat.layout.label
    elseif reference_channel == :mastoid
        return [:M1, :M2]
    end
end
@add_nonmutating rereference!

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
    rereference!(dat, dat.layout.label, reference_channel)
    return nothing
end





"""
    rereference(dat::Union{ContinuousData,ErpData,EpochData}, channel_labels, reference_channel)
    rereference(dat::Union{ContinuousData,ErpData,EpochData}, reference_channel)

Create rereferenced copy of EEG data.

# Arguments
- `dat::Union{ContinuousData,ErpData,EpochData}`: The data to rereference
- `channel_labels::Vector{Union{Symbol, String}}`: Names of channels to rereference (optional)
- `reference_channel`: Channels to use as reference (see DataFrame method for options)

# Returns
- New rereferenced copy of input data
- If channel_labels omitted, uses all channels from layout
"""
function rereference(
    dat::Union{ContinuousData, ErpData},
    channel_labels::Vector{Symbol},
    reference_channel::Vector{Symbol},
)
    dat_out = deepcopy(dat)
    rereference!(dat_out, channel_labels, reference_channel)
    return dat_out
end

function rereference(
    dat::EpochData,
    channel_labels::Vector{Symbol},
    reference_channel::Vector{Symbol},
)
    dat_out = deepcopy(dat)
    for epoch in dat_out.data
        rereference!(epoch, channel_labels, reference_channel)
    end
    return dat_out
end

function rereference(dat::Union{ContinuousData,ErpData}, reference_channel::Vector{Symbol})
    dat_out = deepcopy(dat)
    rereference!(dat_out, dat_out.layout.label, reference_channel)
    return dat_out
end

function rereference(dat::EpochData, reference_channel::Vector{Symbol})
    dat_out = deepcopy(dat)
    for epoch in dat_out.data
        rereference!(epoch, dat_out.layout.label, reference_channel)
    end
    return dat_out
end
