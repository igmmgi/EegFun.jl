"""
    _apply_rereference!(dat::DataFrame, channel_labels, reference::Vector{<:Real})

Internal function that applies rereferencing to specified channels in a DataFrame.
"""
function _apply_rereference!(
    dat::DataFrame, 
    channel_labels::Vector{<:Union{Symbol,AbstractString}}, 
    reference::Vector{<:Real}
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
"""
function calculate_reference(dat::DataFrame, reference_channels)
    reference = reduce(+, eachcol(dat[:, reference_channels])) ./ length(reference_channels)
    return reference
end



# DataFrame methods with different reference channel types
function rereference!(
    dat::DataFrame, 
    channel_labels::Vector{<:Union{Symbol,AbstractString}}, 
    reference_channels::Union{Int,Vector{Int}, UnitRange}
)
    reference_channels = channel_number_to_channel_label(channel_labels, reference_channels)
    @info "Rereferencing using channel(s): $(print_vector_(reference_channels))"
    reference = calculate_reference(dat, reference_channels)
    _apply_rereference!(dat, channel_labels, reference)
end

function rereference!(
    dat::DataFrame, 
    channel_labels::Vector{<:Union{Symbol,AbstractString}}, 
    reference_channels::Vector{<:Union{AbstractString,Symbol}}
)
    @info "Rereferencing using channel(s): $(print_vector_(reference_channels))"
    reference = calculate_reference(dat, reference_channels)
    _apply_rereference!(dat, channel_labels, reference)
end








# Methods for different data types
function rereference!(dat::Union{ContinuousData,ErpData}, channel_labels, reference_channel)
    rereference!(dat.data, channel_labels, reference_channel)
end

function rereference!(dat::EpochData, channel_labels, reference_channel)
    for epoch in dat.data
        rereference!(epoch, channel_labels, reference_channel)
    end
end

function rereference!(dat::Union{ContinuousData,ErpData,EpochData}, reference_channel)
    rereference!(dat, dat.layout.label, reference_channel)
end

# Non-mutating versions
function rereference(dat::Union{ContinuousData, ErpData, EpochData}, channel_labels, reference_channel)
    dat_out = deepcopy(dat)
    rereference!(dat_out, channel_labels, reference_channel)
    return dat_out
end

function rereference(dat::Union{ContinuousData,ErpData,EpochData}, reference_channel)
    dat_out = deepcopy(dat)
    rereference!(dat_out, dat_out.layout.label, reference_channel)
    return dat_out
end


