"""
Channel difference for time-frequency data.

Provides TF-specific `channel_difference!` and `channel_difference` that apply
differencing to both `data_power` and `data_phase` DataFrames.
"""

"""
    channel_difference!(dat::TimeFreqData; channel_selection1::Function = channels(),
                        channel_selection2::Function = channels(),
                        channel_out::Symbol = :diff) -> Nothing

Calculate channel difference for TimeFreqData, applying to both power and phase DataFrames.

See `channel_difference!` for full documentation. This method ensures both `data_power`
and `data_phase` are updated consistently.

# Examples
```julia
# Calculate difference between two channel groups
channel_difference!(tf_data,
    channel_selection1 = channels([:Ch1, :Ch2]),
    channel_selection2 = channels([:Ch3]),
    channel_out = :laterality)
```
"""
function channel_difference!(
    dat::TimeFreqData;
    channel_selection1::Function = channels(),
    channel_selection2::Function = channels(),
    channel_out::Symbol = :diff,
)::Nothing
    # Resolve predicates
    channels_in1 = get_selected_channels(dat, channel_selection1, include_meta = false)
    channels_in2 = get_selected_channels(dat, channel_selection2, include_meta = false)
    @info "channel_difference! (TF): $(channels_in1) vs. $(channels_in2) → :$(channel_out)"

    # Verify channels exist in data
    all_selected = unique(vcat(channels_in1, channels_in2))
    missing_channels = [ch for ch in all_selected if !hasproperty(dat.data_power, ch)]
    if !isempty(missing_channels)
        @minimal_error "Missing channels in TF data: $(missing_channels)"
    end

    # Check for overwrites
    if hasproperty(dat.data_power, channel_out)
        @minimal_warning "Overwriting existing channel '$(channel_out)'"
    end

    # Apply difference to BOTH power and phase DataFrames
    _calculate_channel_difference!(dat.data_power, channels_in1, channels_in2, channel_out)
    _calculate_channel_difference!(dat.data_phase, channels_in1, channels_in2, channel_out)

    return nothing
end


"""
    channel_difference(dat::TimeFreqData; kwargs...) -> TimeFreqData

Non-mutating version of `channel_difference!` for TimeFreqData.
Creates a copy and applies the operation.

# Examples
```julia
# Calculate laterality difference (non-mutating)
result = EegFun.channel_difference(tf_data,
    channel_selection1 = EegFun.channels([:C3]),
    channel_selection2 = EegFun.channels([:C4]),
    channel_out = :laterality)
```
"""
function channel_difference(dat::TimeFreqData; kwargs...)
    data_copy = copy(dat)
    channel_difference!(data_copy; kwargs...)
    return data_copy
end

# Handle Vector{TimeFreqData}
function channel_difference(dat::Vector{TimeFreqData}; kwargs...)
    return [channel_difference(d; kwargs...) for d in dat]
end
