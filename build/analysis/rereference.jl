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

# Examples

## Basic Usage
```julia
# Calculate average reference from all channels
avg_ref = calculate_reference(dat, :avg)

# Calculate mastoid reference (M1 + M2)
mastoid_ref = calculate_reference(dat, :mastoid)
mastoid_ref = calculate_reference(dat, [:M1, :M2])

# Calculate single channel reference
single_ref = calculate_reference(dat, [:M1])

## Channel Filtering Examples
# Reference from EEG channels only (exclude EOG, etc.)
eeg_channels = channels_not([:vEOG, :hEOG, :M1, :M2])
eeg_ref = calculate_reference(dat, eeg_channels)

# Reference from frontal channels only
frontal_channels = channels(x -> startswith.(string.(x), "F"))
frontal_ref = calculate_reference(dat, frontal_channels)

# Reference from parietal channels only
parietal_channels = channels(x -> startswith.(string.(x), "P"))
parietal_ref = calculate_reference(dat, parietal_channels)
```
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
    @info "Rereferencing using channel(s): $(_print_vector(reference_channels))"
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

# Examples

## Basic Usage
```julia
# Average reference (all channels)
rereference!(dat, :avg)

# Mastoid reference (M1 + M2)
rereference!(dat, :mastoid)

# Single channel reference
rereference!(dat, :M1)

# Multiple channel reference
rereference!(dat, [:M1, :M2])
```

## Channel-Specific Rereferencing
```julia
# Rereference only EEG channels (exclude EOG, etc.)
eeg_channels = channels_not([:vEOG, :hEOG, :M1, :M2])
rereference!(dat, eeg_channels, :avg)

# Rereference only frontal channels
frontal_channels = channels(x -> startswith.(string.(x), "F"))
rereference!(dat, frontal_channels, :mastoid)

# Rereference specific channels
rereference!(dat, [:Fp1, :Fp2, :F3, :F4], [:M1, :M2])
```

## Different Reference Types
```julia
# Average reference (most common)
rereference!(dat, :avg)

# Mastoid reference (traditional)
rereference!(dat, :mastoid)

# Single mastoid reference
rereference!(dat, :M1)

# Linked mastoids (average of M1 and M2)
rereference!(dat, [:M1, :M2])

# Custom reference channels
rereference!(dat, [:P9, :P10])  # Posterior reference
```

## Data Type Specific Examples
```julia
# Continuous data
rereference!(dat_continuous, :avg)

# ERP data
rereference!(dat_erp, :mastoid)

# Epoch data (applies to all epochs)
rereference!(dat_epochs, :avg)
```
"""

# helper function to handle special reference cases such as :avg and :mastoid
function resolve_reference(dat, reference_channel::Symbol)
    if reference_channel == :avg # all channels
        reference_channel = channels(dat) 
    elseif reference_channel == :mastoid
        reference_channel = [:M1, :M2]
    else
        reference_channel = [reference_channel]
    end
    return reference_channel
end

# helper function to handle special reference cases such as :avg and :mastoid
function resolve_reference(dat, reference_channel::Vector{Symbol})
    if reference_channel[1] == :avg # all channels
        reference_channel = channels(dat) 
    elseif reference_channel[1] == :mastoid
        reference_channel = [:M1, :M2]
    end
    dat.analysis_info.reference = reference_channel
    return reference_channel
end

# Base methods for all references
function rereference!(
    dat::SingleDataFrameEeg,
    channel_labels::Vector{Symbol},
    reference_channel::Union{Symbol,Vector{Symbol}},
)
    dat.analysis_info.reference = reference_channel
    ref_channels = resolve_reference(dat, reference_channel)
    rereference!(dat.data, channel_labels, ref_channels)
    dat.analysis_info.reference =
        reference_channel isa Symbol ? reference_channel : Symbol(join(reference_channel, '_'))
    return nothing
end

function rereference!(
    dat::MultiDataFrameEeg,
    channel_labels::Vector{Symbol},
    reference_channel::Union{Symbol,Vector{Symbol}},
)
    dat.analysis_info.reference = reference_channel
    ref_channels = resolve_reference(dat, reference_channel)
    for epoch in dat.data
        rereference!(epoch, channel_labels, ref_channels)
    end
    dat.analysis_info.reference =
        reference_channel isa Symbol ? reference_channel : Symbol(join(reference_channel, '_'))
    return nothing
end

function rereference!(dat::EegData, reference_channel::Union{Symbol,Vector{Symbol}})
    dat.analysis_info.reference = reference_channel
    rereference!(dat, channels(dat), reference_channel)
    return nothing
end

# generates all non-mutating versions
@add_nonmutating rereference!





