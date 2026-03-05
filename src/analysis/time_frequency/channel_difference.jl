"""
Channel difference for time-frequency data.

Provides TF-specific `channel_difference!` and `channel_difference` that apply
differencing to both `data_power` and `data_phase` DataFrames.
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


function channel_difference(dat::TimeFreqData; kwargs...)
    data_copy = copy(dat)
    channel_difference!(data_copy; kwargs...)
    return data_copy
end

# Handle Vector{TimeFreqData}
function channel_difference(dat::Vector{TimeFreqData}; kwargs...)
    return [channel_difference(d; kwargs...) for d in dat]
end
