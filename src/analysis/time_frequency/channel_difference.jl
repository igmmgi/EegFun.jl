"""
Channel difference for time-frequency data.

Provides TF-specific `channel_difference!` and `channel_difference` that apply
differencing to both `data_power` and `data_phase` DataFrames.
"""

function channel_difference!(
    dat::TimeFreqData;
    channel_selection1 = channels(),
    channel_selection2 = channels(),
    channel_out::Union{Symbol,AbstractString} = :diff,
)
    out_sym = Symbol(channel_out)
    # Resolve predicates
    channels_in1 = get_selected_channels(dat, channel_selection1, include_meta = false)
    channels_in2 = get_selected_channels(dat, channel_selection2, include_meta = false)
    @info "channel_difference! (TF): $(channels_in1) vs. $(channels_in2) → :$(out_sym)"

    # Verify channels exist in data
    all_selected = unique(vcat(channels_in1, channels_in2))
    missing_channels = [ch for ch in all_selected if !hasproperty(dat.data_power, ch)]
    if !isempty(missing_channels)
        @minimal_error "Missing channels in TF data: $(missing_channels)"
    end

    # Check for overwrites
    if hasproperty(dat.data_power, out_sym)
        @minimal_warning "Overwriting existing channel '$(out_sym)'"
    end

    # Apply difference to BOTH power and phase DataFrames
    _calculate_channel_difference!(dat.data_power, channels_in1, channels_in2, out_sym)
    _calculate_channel_difference!(dat.data_phase, channels_in1, channels_in2, out_sym)

    return dat
end

function channel_difference!(dat::AbstractVector{<:TimeFreqData}; kwargs...)
    for d in dat
        channel_difference!(d; kwargs...)
    end
    return dat
end

function channel_difference(dat::TimeFreqData; kwargs...)
    return channel_difference!(copy(dat); kwargs...)
end

# Handle Vector{TimeFreqData}
function channel_difference(dat::AbstractVector{<:TimeFreqData}; kwargs...)
    return [channel_difference(d; kwargs...) for d in dat]
end
