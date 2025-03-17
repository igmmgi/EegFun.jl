using Logging

"""
    calculate_channel_difference(dat::DataFrame, channels1::Vector{Symbol}, channels2::Vector{Symbol})

Calculate the mean difference between two sets of channels.

# Arguments
- `dat::DataFrame`: DataFrame containing EEG data
- `channels1::Vector{Symbol}`: First set of channel symbols (e.g., [:Fp1, :Fp2])
- `channels2::Vector{Symbol}`: Second set of channel symbols

# Returns
- `Vector{Float64}`: Difference between mean of channels1 and mean of channels2
"""
function calculate_channel_difference(dat::DataFrame, channels1::Vector{Symbol}, channels2::Vector{Symbol})
    @debug "Calculating channel difference: $(channels1) vs. $(channels2)"

    # Verify channels exist in data
    missing_channels = [ch for ch in vcat(channels1, channels2) if ch ∉ propertynames(dat)]
    if !isempty(missing_channels)
        @error "Missing channels in data: $(missing_channels)"
        throw(ArgumentError("Channels not found in data: $missing_channels"))
    end

    # Pre-allocate vectors for better performance
    n_samples = size(dat, 1)
    mean1 = zeros(n_samples)
    mean2 = zeros(n_samples)

    # Calculate means in-place
    @debug "Calculating mean for first channel set: $(channels1)"
    for channel in channels1
        @views mean1 .+= dat[!, channel]
    end
    mean1 ./= length(channels1)

    @debug "Calculating mean for second channel set: $(channels2)"
    for channel in channels2
        @views mean2 .+= dat[!, channel]
    end
    mean2 ./= length(channels2)

    # Calculate difference in-place
    mean1 .-= mean2
    return mean1
end

"""
    calculate_channel_difference(dat::DataFrame, channel1::Symbol, channel2::Symbol)

Calculate the difference between two individual channels.

# Arguments
- `dat::DataFrame`: DataFrame containing EEG data
- `channel1::Symbol`: First channel symbol (e.g., :Fp1)
- `channel2::Symbol`: Second channel symbol

# Returns
- `Vector{Float64}`: Difference between channel1 and channel2
"""
function calculate_channel_difference(dat::DataFrame, channel1::Symbol, channel2::Symbol)
    return calculate_channel_difference(dat, [channel1], [channel2])
end

"""
    diff_channel!(dat::DataFrame, channels1::Vector{Symbol}, 
                 channels2::Vector{Symbol}, 
                 difference_label::Symbol)

Calculate and add channel difference to DataFrame in-place.

# Arguments
- `dat::DataFrame`: DataFrame to modify
- `channels1::Vector{Symbol}`: First set of channel symbols
- `channels2::Vector{Symbol}`: Second set of channel symbols
- `difference_label::Symbol`: Column name for the difference (e.g., :vEOG)
"""
function diff_channel!(dat::DataFrame, channels1::Vector{Symbol}, channels2::Vector{Symbol}, difference_label::Symbol)
    @info "Adding new channel '$(difference_label)' as difference between $(channels1) and $(channels2)"

    # Check if difference channel already exists
    if difference_label ∈ propertynames(dat)
        @warn "Overwriting existing channel '$(difference_label)'"
    end

    dat[!, difference_label] = calculate_channel_difference(dat, channels1, channels2)
    return nothing
end

"""
    diff_channel!(dat::DataFrame, channel1::Symbol, channel2::Symbol, 
                 difference_label::Symbol)

Calculate and add difference between two channels to DataFrame in-place.

# Arguments
- `dat::DataFrame`: DataFrame to modify
- `channel1::Symbol`: First channel symbol
- `channel2::Symbol`: Second channel symbol
- `difference_label::Symbol`: Column name for the difference (e.g., :vEOG)
"""
function diff_channel!(dat::DataFrame, channel1::Symbol, channel2::Symbol, difference_label::Symbol)
    diff_channel!(dat, [channel1], [channel2], difference_label)
    return nothing
end

"""
    diff_channel!(dat::Union{ContinuousData,ErpData}, 
                 channels1::Vector{Symbol}, 
                 channels2::Vector{Symbol}, 
                 difference_label::Symbol)

Calculate and add channel difference to ContinuousData or ErpData in-place.

# Arguments
- `dat::Union{ContinuousData,ErpData}`: Data object to modify
- `channels1::Vector{Symbol}`: First set of channel symbols
- `channels2::Vector{Symbol}`: Second set of channel symbols
- `difference_label::Symbol`: Column name for the difference (e.g., :vEOG)
"""
function diff_channel!(
    dat::Union{ContinuousData,ErpData},
    channels1::Vector{Symbol},
    channels2::Vector{Symbol},
    difference_label::Symbol,
)
    @info "Computing channel difference for $(typeof(dat)) $(difference_label)"
    diff_channel!(dat.data, channels1, channels2, difference_label)
    return nothing
end

"""
    diff_channel!(dat::Union{ContinuousData,ErpData}, channel1::Symbol, 
                 channel2::Symbol, difference_label::Symbol)

Calculate and add difference between two channels to ContinuousData or ErpData in-place.

# Arguments
- `dat::Union{ContinuousData,ErpData}`: Data object to modify
- `channel1::Symbol`: First channel symbol
- `channel2::Symbol`: Second channel symbol
- `difference_label::Symbol`: Column name for the difference (e.g., :vEOG)
"""
function diff_channel!(dat::Union{ContinuousData,ErpData}, channel1::Symbol, channel2::Symbol, difference_label::Symbol)
    diff_channel!(dat, [channel1], [channel2], difference_label)
    return nothing
end

"""
    diff_channel!(dat::EpochData, channels1::Vector{Symbol}, 
                 channels2::Vector{Symbol}, difference_label::Symbol)

Calculate and add channel difference to each epoch in EpochData in-place.

# Arguments
- `dat::EpochData`: EpochData object to modify
- `channels1::Vector{Symbol}`: First set of channel symbols
- `channels2::Vector{Symbol}`: Second set of channel symbols
- `difference_label::Symbol`: Column name for the difference (e.g., :vEOG)
"""
function diff_channel!(dat::EpochData, channels1::Vector{Symbol}, channels2::Vector{Symbol}, difference_label::Symbol)
    @info "Computing channel difference for each epoch $(difference_label) n_epochs=length(dat.data)"
    for (i, epoch) in enumerate(dat.data)
        @debug "Processing epoch" i
        diff_channel!(epoch, channels1, channels2, difference_label)
    end
    @debug "Completed all epochs"
    return nothing
end

"""
    diff_channel!(dat::EpochData, channel1::Symbol, channel2::Symbol, 
                 difference_label::Symbol)

Calculate and add difference between two channels to each epoch in EpochData in-place.

# Arguments
- `dat::EpochData`: EpochData object to modify
- `channel1::Symbol`: First channel symbol
- `channel2::Symbol`: Second channel symbol
- `difference_label::Symbol`: Column name for the difference (e.g., :vEOG)
"""
function diff_channel!(dat::EpochData, channel1::Symbol, channel2::Symbol, difference_label::Symbol)
    diff_channel!(dat, [channel1], [channel2], difference_label)
    return nothing
end

# generates all non-mutating versions
@add_nonmutating diff_channel!
