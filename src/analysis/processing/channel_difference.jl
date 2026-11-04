function _channel_diff_col(df::DataFrame, idx1::Vector{Int}, idx2::Vector{Int})::Vector{Float64}
    n = nrow(df)
    out = Vector{Float64}(undef, n)
    inv1 = 1.0 / length(idx1)
    c1 = df[!, idx1[1]]
    if c1 isa Vector{Float64}
        v1 = c1::Vector{Float64}
        @inbounds @simd for i = 1:n
            out[i] = v1[i] * inv1
        end
    else
        @inbounds @simd for i = 1:n
            out[i] = Float64(c1[i]) * inv1
        end
    end
    for k = 2:length(idx1)
        ck = df[!, idx1[k]]
        if ck isa Vector{Float64}
            vk = ck::Vector{Float64}
            @inbounds @simd for i = 1:n
                out[i] = muladd(vk[i], inv1, out[i])
            end
        else
            @inbounds @simd for i = 1:n
                out[i] = muladd(Float64(ck[i]), inv1, out[i])
            end
        end
    end

    inv2 = 1.0 / length(idx2)
    for k in eachindex(idx2)
        ck = df[!, idx2[k]]
        if ck isa Vector{Float64}
            vk = ck::Vector{Float64}
            @inbounds @simd for i = 1:n
                out[i] = muladd(-vk[i], inv2, out[i])
            end
        else
            @inbounds @simd for i = 1:n
                out[i] = muladd(-Float64(ck[i]), inv2, out[i])
            end
        end
    end
    return out
end

function _calculate_channel_difference!(dat::DataFrame, channels_in1::Vector{Symbol}, channels_in2::Vector{Symbol}, channel_out::Symbol)
    idx1 = Int[columnindex(dat, c) for c in channels_in1]
    idx2 = Int[columnindex(dat, c) for c in channels_in2]
    dat[!, channel_out] = _channel_diff_col(dat, idx1, idx2)
    return nothing
end

function _calculate_channel_difference!(
    dat::AbstractVector{<:DataFrame},
    channels_in1::Vector{Symbol},
    channels_in2::Vector{Symbol},
    channel_out::Symbol,
)
    isempty(dat) && return nothing
    idx1 = Int[columnindex(dat[1], c) for c in channels_in1]
    idx2 = Int[columnindex(dat[1], c) for c in channels_in2]
    work = length(dat) * nrow(dat[1]) * (length(idx1) + length(idx2))
    if length(dat) >= 8 && Threads.nthreads() > 1 && work >= 1_000_000
        Threads.@threads :static for e in eachindex(dat)
            dat[e][!, channel_out] = _channel_diff_col(dat[e], idx1, idx2)
        end
    else
        for df in dat
            df[!, channel_out] = _channel_diff_col(df, idx1, idx2)
        end
    end
    return nothing
end

"""
    channel_difference!(dat::EegData; channel_selection1 = channels(),
                        channel_selection2 = channels(), channel_out = :diff) -> EegData

Compute `mean(channel_selection1) - mean(channel_selection2)` and store it as `channel_out` in `dat`.
"""
function channel_difference!(
    dat::EegData;
    channel_selection1 = channels(),
    channel_selection2 = channels(),
    channel_out::Union{Symbol,AbstractString} = :diff,
)
    out_sym = Symbol(channel_out)
    chs1 = get_selected_channels(dat, channel_selection1; include_meta = false)
    chs2 = get_selected_channels(dat, channel_selection2; include_meta = false)

    isempty(chs1) && @minimal_error "No channels selected by channel_selection1"
    isempty(chs2) && @minimal_error "No channels selected by channel_selection2"

    @info "Channel sets: $(chs1) vs. $(chs2) -> $(out_sym)"

    existing = Set(all_labels(dat))
    missing_channels = [ch for ch in unique(vcat(chs1, chs2)) if ch ∉ existing]
    !isempty(missing_channels) && @minimal_error "Missing channels in data: $(missing_channels)"

    out_sym ∈ existing && @minimal_warning "Overwriting existing channel '$(out_sym)'"
    if out_sym ∈ chs1 || out_sym ∈ chs2
        @minimal_warning "Output channel '$(out_sym)' is also in the input channel sets. This may cause unexpected behavior."
    end

    for container in _data_frames(dat)
        _calculate_channel_difference!(container, chs1, chs2, out_sym)
    end
    return dat
end

"""
    channel_difference(dat, args...; kwargs...)

Non-mutating version of [`channel_difference!`](@ref).
"""
channel_difference(dat::EegData, args...; kwargs...) = channel_difference!(copy(dat), args...; kwargs...)

function channel_difference!(dat::AbstractVector{<:EegData}, args...; kwargs...)
    for d in dat
        channel_difference!(d, args...; kwargs...)
    end
    return dat
end

channel_difference(dat::AbstractVector{<:EegData}, args...; kwargs...) = [channel_difference(d, args...; kwargs...) for d in dat]

"""
    calculate_eog_channels!(dat::EegData, eog_cfg::Union{EogConfig, Dict}) -> EegData

Compute vertical and horizontal bipolar EOG channels in place according to `eog_cfg`.
"""
function calculate_eog_channels!(dat::EegData, eog_cfg::EogConfig)
    v_cfg = eog_cfg.vEOG_channels
    channel_difference!(
        dat;
        channel_selection1 = channels(Symbol.(v_cfg[1])),
        channel_selection2 = channels(Symbol.(v_cfg[2])),
        channel_out = Symbol(v_cfg[3][1]),
    )

    h_cfg = eog_cfg.hEOG_channels
    channel_difference!(
        dat;
        channel_selection1 = channels(Symbol.(h_cfg[1])),
        channel_selection2 = channels(Symbol.(h_cfg[2])),
        channel_out = Symbol(h_cfg[3][1]),
    )
    return dat
end

calculate_eog_channels!(dat::EegData, eog_cfg::Dict) = calculate_eog_channels!(dat, EogConfig(eog_cfg))

function calculate_eog_channels!(dat::AbstractVector{<:EegData}, eog_cfg)
    for d in dat
        calculate_eog_channels!(d, eog_cfg)
    end
    return dat
end

"""
    calculate_eog_channels(dat, eog_cfg)

Non-mutating version of [`calculate_eog_channels!`](@ref).
"""
calculate_eog_channels(dat::EegData, args...; kwargs...) = calculate_eog_channels!(copy(dat), args...; kwargs...)
calculate_eog_channels(dat::AbstractVector{<:EegData}, args...; kwargs...) = [calculate_eog_channels(d, args...; kwargs...) for d in dat]

"""
    detect_eog_onsets!(dat::EegData, eog_cfg::EogConfig) -> EegData

Detect vertical and horizontal EOG step onsets in place using the criteria and channel names in `eog_cfg`.
"""
function detect_eog_onsets!(dat::EegData, eog_cfg::EogConfig)
    v_out = eog_cfg.vEOG_channels[3][1]
    detect_eog_onsets!(dat, eog_cfg.vEOG_criterion, Symbol(v_out), Symbol("is_" * v_out))

    h_out = eog_cfg.hEOG_channels[3][1]
    detect_eog_onsets!(dat, eog_cfg.hEOG_criterion, Symbol(h_out), Symbol("is_" * h_out))
    return dat
end
