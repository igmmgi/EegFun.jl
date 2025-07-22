"""
    calculate_channel_difference(dat::DataFrame; channels_in1::Function = channels(), channels_in2::Function = channels())

Calculate the mean difference between two sets of channels using predicates.

# Arguments
- `dat::DataFrame`: DataFrame containing EEG data
- `channels_in1::Function`: Channel selection predicate for first set (default: channels() - all channels)
- `channels_in2::Function`: Channel selection predicate for second set (default: channels() - all channels)

# Returns
- `Vector{Float64}`: Difference between mean of channels_in1 and mean of channels_in2

# Examples
```julia
# Using symbols
calculate_channel_difference(dat, channels_in1 = channels(:Fp1), channels_in2 = channels(:Fp2))

# Using vectors of symbols
calculate_channel_difference(dat, channels_in1 = channels([:Fp1, :Fp2]), channels_in2 = channels([:IO1, :IO2]))

# Using functions
calculate_channel_difference(dat, channels_in1 = channels(x -> startswith.(string.(x), "F")), channels_in2 = channels(x -> startswith.(string.(x), "I")))

# Using default (all channels vs all channels - will be zero)
calculate_channel_difference(dat)
```
"""
function calculate_channel_difference(dat::DataFrame; channels_in1::Function = channels(), channels_in2::Function = channels())
    # Convert predicates to channel symbols using existing infrastructure
    channels_in1_selected = get_selected_channels(dat, channels_in1)
    channels_in2_selected = get_selected_channels(dat, channels_in2)
    
    @debug "Calculating channel difference: $(channels_in1_selected) vs. $(channels_in2_selected)"

    # Verify channels exist in data
    missing_channels = [ch for ch in vcat(channels_in1_selected, channels_in2_selected) if ch ∉ propertynames(dat)]
    if !isempty(missing_channels)
        @error "Missing channels in data: $(missing_channels)"
        throw(ArgumentError("Channels not found in data: $missing_channels"))
    end

    # Pre-allocate vectors for better performance
    n_samples = size(dat, 1)
    mean1 = zeros(n_samples)
    mean2 = zeros(n_samples)

    # Calculate means in-place
    @debug "Calculating mean for first channel set: $(channels_in1_selected)"
    for channel in channels_in1_selected
        @views mean1 .+= dat[!, channel]
    end
    mean1 ./= length(channels_in1_selected)

    @debug "Calculating mean for second channel set: $(channels_in2_selected)"
    for channel in channels_in2_selected
        @views mean2 .+= dat[!, channel]
    end
    mean2 ./= length(channels_in2_selected)

    # Calculate difference in-place
    mean1 .-= mean2
    return mean1
end

"""
    channel_difference!(dat::DataFrame; channels_in1::Function = channels(), channels_in2::Function = channels(), channel_out::Symbol = :diff)

Calculate and add channel difference to DataFrame in-place using predicates.

# Arguments
- `dat::DataFrame`: DataFrame to modify
- `channels_in1::Function`: Channel selection predicate for first set (default: channels() - all channels)
- `channels_in2::Function`: Channel selection predicate for second set (default: channels() - all channels)
- `channel_out::Symbol`: Column name for the difference (default: :diff)

# Examples
```julia
# Using symbols
channel_difference!(dat, channels_in1 = channels(:Fp1), channels_in2 = channels(:Fp2), channel_out = :vEOG)

# Using vectors of symbols
channel_difference!(dat, channels_in1 = channels([:Fp1, :Fp2]), channels_in2 = channels([:IO1, :IO2]), channel_out = :vEOG)

# Using functions
channel_difference!(dat, channels_in1 = channels(x -> startswith.(string.(x), "F")), channels_in2 = channels(x -> startswith.(string.(x), "I")), channel_out = :vEOG)

# Using default channel_out
channel_difference!(dat, channels_in1 = channels(:Fp1), channels_in2 = channels(:Fp2))  # Creates :diff column
```
"""
function channel_difference!(dat::DataFrame; channels_in1::Function = channels(), channels_in2::Function = channels(), channel_out::Symbol = :diff)
    @info "Adding new channel '$(channel_out)' as difference between channels"
    
    # Resolve predicates to get actual channel names for logging
    channels_in1_selected = get_selected_channels(dat, channels_in1)
    channels_in2_selected = get_selected_channels(dat, channels_in2)
    @info "Channel sets: $(channels_in1_selected) vs. $(channels_in2_selected)"

    # Check if difference channel already exists
    if channel_out ∈ propertynames(dat)
        @minimal_warning "Overwriting existing channel '$(channel_out)'"
    end

    dat[!, channel_out] = calculate_channel_difference(dat, channels_in1 = channels_in1, channels_in2 = channels_in2)
    return nothing
end

"""
    channel_difference!(dat::Union{ContinuousData,ErpData}; channels_in1::Function = channels(), channels_in2::Function = channels(), channel_out::Symbol = :diff)

Calculate and add channel difference to ContinuousData or ErpData in-place using predicates.

# Arguments
- `dat::Union{ContinuousData,ErpData}`: Data object to modify
- `channels_in1::Function`: Channel selection predicate for first set (default: channels() - all channels)
- `channels_in2::Function`: Channel selection predicate for second set (default: channels() - all channels)
- `channel_out::Symbol`: Column name for the difference (default: :diff)
"""
function channel_difference!(
    dat::Union{ContinuousData,ErpData};
    channels_in1::Function = channels(),
    channels_in2::Function = channels(),
    channel_out::Symbol = :diff,
)
    @info "Computing channel difference for $(typeof(dat)) $(channel_out)"
    channel_difference!(dat.data, channels_in1 = channels_in1, channels_in2 = channels_in2, channel_out = channel_out)
    return nothing
end

"""
    channel_difference!(dat::EpochData; channels_in1::Function = channels(), channels_in2::Function = channels(), channel_out::Symbol = :diff)

Calculate and add channel difference to each epoch in EpochData in-place using predicates.

# Arguments
- `dat::EpochData`: EpochData object to modify
- `channels_in1::Function`: Channel selection predicate for first set (default: channels() - all channels)
- `channels_in2::Function`: Channel selection predicate for second set (default: channels() - all channels)
- `channel_out::Symbol`: Column name for the difference (default: :diff)
"""
function channel_difference!(dat::EpochData; channels_in1::Function = channels(), channels_in2::Function = channels(), channel_out::Symbol = :diff)
    @info "Computing channel difference for each epoch $(channel_out) n_epochs=$(length(dat.data))"
    for (i, epoch) in enumerate(dat.data)
        @debug "Processing epoch" i
        channel_difference!(epoch, channels_in1 = channels_in1, channels_in2 = channels_in2, channel_out = channel_out)
    end
    @debug "Completed all epochs"
    return nothing
end

# generates all non-mutating versions
@add_nonmutating channel_difference!
