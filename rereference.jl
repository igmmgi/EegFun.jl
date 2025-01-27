using Logging


"""
    _apply_rereference!(dat::DataFrame, channel_labels, reference::Vector{<:Real})

Internal function that applies rereferencing to specified channels in a DataFrame.

# Arguments
- `dat::DataFrame`: The data to rereference
- `channel_labels::Vector{<:Union{Symbol,AbstractString}}`: Names of channels to rereference
- `reference::Vector{<:Real}`: Reference signal to subtract from each channel

# Effects
- Modifies the input DataFrame in-place by subtracting the reference from each channel
"""
function _apply_rereference!(
    dat::DataFrame, 
    channel_labels::Vector{<:Union{Symbol,AbstractString}}, 
    reference::Vector{<:Real}
)
    @debug "Applying reference to channels: $(channel_labels)"
    for col in names(dat)
        if col in channel_labels
            dat[!, col] .-= reference
        end
    end
end

"""
    rereference!(dat::DataFrame, channel_labels, reference_channels::Union{Int,Vector{Int}})

Rereference data using channel numbers as reference.

# Arguments
- `dat::DataFrame`: The data to rereference
- `channel_labels::Vector{<:Union{Symbol,AbstractString}}`: Names of channels to rereference
- `reference_channels::Union{Int,Vector{Int}}`: Channel number(s) to use as reference

# Effects
- Modifies the input data in-place by subtracting the mean of reference channels
"""
function rereference!(
    dat::DataFrame, 
    channel_labels::Vector{<:Union{Symbol,AbstractString}}, 
    reference_channels::Union{Int,Vector{Int}}
)
    @info "Rereferencing using channel numbers: $(reference_channels)"
    reference_channels = channel_number_to_channel_label(channel_labels, reference_channels)
    reference = reduce(+, eachcol(dat[:, reference_channels])) ./ length(reference_channels)
    _apply_rereference!(dat, channel_labels, reference)
end

"""
    rereference!(dat::DataFrame, channel_labels, reference_channels::Union{AbstractString,Symbol})

Rereference data using a single channel name as reference.

# Arguments
- `dat::DataFrame`: The data to rereference
- `channel_labels::Vector{<:Union{Symbol,AbstractString}}`: Names of channels to rereference
- `reference_channels::Union{AbstractString,Symbol}`: Channel name to use as reference

# Effects
- Modifies the input data in-place by subtracting the mean of reference channels
"""
function rereference!(
    dat::DataFrame, 
    channel_labels::Vector{<:Union{Symbol,AbstractString}}, 
    reference_channels::Union{AbstractString,Symbol}
)
    @info "Rereferencing using single channel: $(reference_channels)"
    rereference!(dat, channel_labels, [reference_channels])
end

"""
    rereference!(dat::DataFrame, channel_labels, reference_channels::Vector{<:Union{AbstractString,Symbol}})

Rereference data using multiple channel names as reference.

# Arguments
- `dat::DataFrame`: The data to rereference
- `channel_labels::Vector{<:Union{Symbol,AbstractString}}`: Names of channels to rereference
- `reference_channels::Vector{<:Union{AbstractString,Symbol}}`: Channel names to use as reference

# Effects
- Modifies the input data in-place by subtracting the mean of reference channels
"""
function rereference!(
    dat::DataFrame, 
    channel_labels::Vector{<:Union{Symbol,AbstractString}}, 
    reference_channels::Vector{<:Union{AbstractString,Symbol}}
)
    @info "Rereferencing using multiple channels: $(reference_channels)"
    reference = reduce(+, eachcol(dat[:, reference_channels])) ./ length(reference_channels)
    _apply_rereference!(dat, channel_labels, reference)
end

# ContinuousData methods
"""
    rereference!(dat::ContinuousData, channel_labels, reference_channel)

Rereference ContinuousData using specified channels as reference.

# Arguments
- `dat::ContinuousData`: The EEG data object to rereference
- `channel_labels::Vector{<:Union{Symbol,AbstractString}}`: Names of channels to rereference
- `reference_channel`: Channel(s) to use as reference

# Effects
- Modifies the input data in-place by applying the rereferencing
"""
function rereference!(dat::ContinuousData, channel_labels, reference_channel)
    @info "Rereferencing continuous data"
    rereference!(dat.data, channel_labels, reference_channel)
end

"""
    rereference!(dat::ContinuousData, reference_channel)

Rereference ContinuousData using specified channels as reference,
applying to all channels from the layout.

# Effects
- Modifies the input data in-place by applying the rereferencing
"""
function rereference!(dat::ContinuousData, reference_channel)
    @debug "Using all channels from layout for rereferencing"
    rereference!(dat, dat.layout.label, reference_channel)
end

# ErpData methods
"""
    rereference!(dat::ErpData, channel_labels, reference_channel)

Rereference ErpData using specified channels as reference.

# Effects
- Modifies the input data in-place by applying the rereferencing
"""
function rereference!(dat::ErpData, channel_labels, reference_channel)
    @info "Rereferencing ERP data"
    rereference!(dat.data, channel_labels, reference_channel)
end

"""
    rereference!(dat::ErpData, reference_channel)

Rereference ErpData using specified channels as reference,
applying to all channels from the layout.

# Effects
- Modifies the input data in-place by applying the rereferencing
"""
function rereference!(dat::ErpData, reference_channel)
    @debug "Using all channels from layout for rereferencing"
    rereference!(dat, dat.layout.label, reference_channel)
end

# EpochData methods
"""
    rereference!(dat::EpochData, channel_labels, reference_channel)

Rereference each epoch in EpochData using specified channels as reference.

# Effects
- Modifies the input data in-place by applying the rereferencing
"""
function rereference!(dat::EpochData, channel_labels, reference_channel)
    @info "Rereferencing epoched data" n_epochs=length(dat.data)
    for (i, epoch) in enumerate(dat.data)
        @debug "Processing epoch" i
        rereference!(epoch, channel_labels, reference_channel)
    end
end

"""
    rereference!(dat::EpochData, reference_channel)

Rereference each epoch in EpochData using specified channels as reference,
applying to all channels from the layout.

# Effects
- Modifies the input data in-place by applying the rereferencing
"""
function rereference!(dat::EpochData, reference_channel)
    @debug "Using all channels from layout for rereferencing"
    rereference!(dat, dat.layout.label, reference_channel)
end

# Non-mutating versions
"""
    rereference(dat::ContinuousData, channel_labels, reference_channel)

Create a rereferenced copy of ContinuousData using specified channels as reference.
"""
function rereference(dat::ContinuousData, channel_labels, reference_channel)
    @debug "Creating rereferenced copy of continuous data"
    dat_out = deepcopy(dat)
    rereference!(dat_out, channel_labels, reference_channel)
    return dat_out
end

"""
    rereference(dat::ContinuousData, reference_channel)

Create a rereferenced copy of ContinuousData using all channels from layout.
"""
function rereference(dat::ContinuousData, reference_channel)
    @debug "Creating rereferenced copy using all channels from layout"
    dat_out = deepcopy(dat)
    rereference!(dat_out, dat_out.layout.label, reference_channel)
    return dat_out
end

"""
    rereference(dat::ErpData, channel_labels, reference_channel)

Create a rereferenced copy of ErpData using specified channels as reference.
"""
function rereference(dat::ErpData, channel_labels, reference_channel)
    @debug "Creating rereferenced copy of ERP data"
    dat_out = deepcopy(dat)
    rereference!(dat_out, channel_labels, reference_channel)
    return dat_out
end

"""
    rereference(dat::ErpData, reference_channel)

Create a rereferenced copy of ErpData using all channels from layout.
"""
function rereference(dat::ErpData, reference_channel)
    @debug "Creating rereferenced copy using all channels from layout"
    dat_out = deepcopy(dat)
    rereference!(dat_out, dat_out.layout.label, reference_channel)
    return dat_out
end

"""
    rereference(dat::EpochData, channel_labels, reference_channel)

Create a rereferenced copy of EpochData using specified channels as reference.
"""
function rereference(dat::EpochData, channel_labels, reference_channel)
    @debug "Creating rereferenced copy of epoched data"
    dat_out = deepcopy(dat)
    rereference!(dat_out, channel_labels, reference_channel)
    return dat_out
end

"""
    rereference(dat::EpochData, reference_channel)

Create a rereferenced copy of EpochData using all channels from layout.
"""
function rereference(dat::EpochData, reference_channel)
    @debug "Creating rereferenced copy using all channels from layout"
    rereference!(dat_out, dat_out.layout.label, reference_channel)
    return dat_out
end
