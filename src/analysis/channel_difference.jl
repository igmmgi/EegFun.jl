# Internal function that only accepts resolved channels and mutates the DataFrame
function _calculate_channel_difference!(
    dat::DataFrame,
    channels_in1::Vector{Symbol},
    channels_in2::Vector{Symbol},
    channel_out::Symbol,
)

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
function _calculate_channel_difference!(
    dat::Vector{DataFrame},
    channels_in1::Vector{Symbol},
    channels_in2::Vector{Symbol},
    channel_out::Symbol,
)
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

"""
    calculate_eog_channels!(dat::EegData, eog_cfg::Dict)

Calculate EOG channels based on configuration dictionary. Modifies the data in place.

# Arguments
- `dat::EegData`: EegData object to modify
- `eog_cfg::Dict`: Configuration dictionary containing EOG channel settings

# Configuration Structure
The `eog_cfg` should contain:
- `"vEOG_channels"`: Vector of 3 elements: [channels1, channels2, output_name]
- `"hEOG_channels"`: Vector of 3 elements: [channels1, channels2, output_name]

# Examples
```julia
# Calculate EOG channels based on config
calculate_eog_channels!(dat, cfg["preprocess"]["eog"])

# Manual configuration
eog_cfg = Dict(
    "vEOG_channels" => [["Fp1", "Fp2"], ["IO1", "IO2"], ["vEOG"]],
    "hEOG_channels" => [["F9"], ["F10"], ["hEOG"]]
)
calculate_eog_channels!(dat, eog_cfg)
```
"""
function calculate_eog_channels!(dat::EegData, eog_cfg::EogConfig)
    # Calculate vertical EOG channels
    vEOG_cfg = eog_cfg.vEOG_channels
    channel_difference!(
        dat,
        channel_selection1 = channels(Symbol.(vEOG_cfg[1])),
        channel_selection2 = channels(Symbol.(vEOG_cfg[2])),
        channel_out = Symbol(vEOG_cfg[3][1])
    )

    # Calculate horizontal EOG channels
    hEOG_cfg = eog_cfg.hEOG_channels
    channel_difference!(
        dat,
        channel_selection1 = channels(Symbol.(hEOG_cfg[1])),
        channel_selection2 = channels(Symbol.(hEOG_cfg[2])),
        channel_out = Symbol(hEOG_cfg[3][1])
    )
end

# Convenience method for Dict input
function calculate_eog_channels!(dat::EegData, eog_cfg::Dict)
    calculate_eog_channels!(dat, EogConfig(eog_cfg))
end

"""
    detect_eog_signals!(dat::EegData, eog_cfg::EogConfig)

Detect EOG onsets for both vertical and horizontal EOG channels based on configuration.

# Arguments
- `dat::EegData`: The EEG data object
- `eog_cfg::EogConfig`: EOG configuration object containing vEOG and hEOG settings

# Example
```julia
eog_cfg = EogConfig(
    vEOG_criterion = 50.0,
    hEOG_criterion = 50.0,
    vEOG_channels = [["Fp1"], ["Fp2"], ["vEOG"]],
    hEOG_channels = [["F9"], ["F10"], ["hEOG"]]
)
detect_eog_signals!(dat, eog_cfg)
```
"""
function detect_eog_signals!(dat::EegData, eog_cfg::EogConfig)
    # Detect vertical EOG onsets
    vEOG_cfg = eog_cfg.vEOG_channels
    detect_eog_onsets!(
        dat,
        eog_cfg.vEOG_criterion,
        Symbol(vEOG_cfg[3][1]),
        Symbol("is_" * vEOG_cfg[3][1])
    )
    
    # Detect horizontal EOG onsets
    hEOG_cfg = eog_cfg.hEOG_channels
    detect_eog_onsets!(
        dat,
        eog_cfg.hEOG_criterion,
        Symbol(hEOG_cfg[3][1]),
        Symbol("is_" * hEOG_cfg[3][1])
    )
end


