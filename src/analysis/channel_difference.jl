# Internal function that only accepts resolved channels and mutates the DataFrame
function _calculate_channel_difference!(dat::DataFrame, channels_in1::Vector{Symbol}, channels_in2::Vector{Symbol}, channel_out::Symbol)

    # Pre-allocate vectors for better performance
    means = zeros(n_samples(dat), 2)

    # Calculate means in-place
    for channel in channels_in1
        @views means[:, 1] .+= dat[!, channel]
    end
    means[:, 1] ./= length(channels_in1)

    for channel in channels_in2
        @views means[:, 2] .+= dat[!, channel]
    end
    means[:, 2] ./= length(channels_in2)

    # Calculate difference and assign directly to DataFrame
    dat[!, channel_out] = means[:, 1] .- means[:, 2]

    return nothing
end

# Internal function that works on a vector of DataFrames
function _calculate_channel_difference!(dat::Vector{DataFrame}, channels_in1::Vector{Symbol}, channels_in2::Vector{Symbol}, channel_out::Symbol)
    _calculate_channel_difference!.(dat, Ref(channels_in1), Ref(channels_in2), Ref(channel_out))
    return nothing
end


"""
    channel_difference!(dat::EegData; channel_selection1::Function = channels(), channel_selection2::Function = channels(), channel_out::Symbol = :diff)

Calculate and add channel difference to EEG data in-place using predicates.

# Arguments
- `dat::EegData`: EEG data object to modify (ContinuousData, ErpData, or EpochData)
- `channel_selection1::Function`: Channel selection predicate for first set (default: channels() - all channels)
- `channel_selection2::Function`: Channel selection predicate for second set (default: channels() - all channels)
- `channel_out::Symbol`: Column name for the difference (default: :diff)

# Examples
```julia
# Calculate difference between frontal and parietal channels
channel_difference!(dat, 
    channel_selection1 = channels(x -> startswith.(string.(x), "F")),  # Frontal channels
    channel_selection2 = channels(x -> startswith.(string.(x), "P"))   # Parietal channels
)

# Calculate difference between specific channels
channel_difference!(dat, 
    channel_selection1 = channels([:Fp1, :Fp2]), 
    channel_selection2 = channels([:O1, :O2]), 
    channel_out = :front_back_diff
)

# Calculate difference between all channels vs all channels (will be zero)
channel_difference!(dat)
```
"""
function channel_difference!(
    dat::EegData;
    channel_selection1::Function = channels(),
    channel_selection2::Function = channels(),
    channel_out::Symbol = :diff,
)
    # Resolve predicates to get actual channel names for logging
    channels_in1_selected = get_selected_channels(dat, channel_selection1, include_meta = false)
    channels_in2_selected = get_selected_channels(dat, channel_selection2, include_meta = false)
    @info "Channel sets: $(channels_in1_selected) vs. $(channels_in2_selected) -> $(channel_out)"

    # Verify channels exist in data
    missing_channels = [ch for ch in vcat(channels_in1_selected, channels_in2_selected) if ch ∉ all_labels(dat)]
    if !isempty(missing_channels)
        @minimal_error "Missing channels in data: $(missing_channels)"
    end

    # Check if difference channel already exists
    if channel_out ∈ all_labels(dat)
        @minimal_warning "Overwriting existing channel '$(channel_out)'"
    end

    # Calculate channel difference (dispatch handles DataFrame vs Vector{DataFrame})
    _calculate_channel_difference!(dat.data, channels_in1_selected, channels_in2_selected, channel_out)

    return nothing
end

# generates all non-mutating versions
@add_nonmutating channel_difference!
