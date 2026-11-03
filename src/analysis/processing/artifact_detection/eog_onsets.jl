"""
    detect_eog_onsets!(dat::ContinuousData, criterion::Real, channel_in::Symbol, channel_out::Symbol; step_size::Int=20)
    detect_eog_onsets!(dat::ContinuousData, eog_cfg::Dict)

Detect electrooculogram (EOG) step onsets in `channel_in` exceeding `criterion` (μV) and record boolean flags in `channel_out`.
When passed `eog_cfg::Dict`, detects both vertical (`vEOG`) and horizontal (`hEOG`) onsets in-place.
"""
function detect_eog_onsets!(dat::ContinuousData, criterion::Real, channel_in::Symbol, channel_out::Symbol; step_size::Int = 20)
    @info "Detecting EOG onsets in channel $(channel_in) with stepsize criterion $(criterion) μV"
    channel_in ∉ propertynames(dat.data) && @minimal_error("channel $(channel_in) not found in data")

    step_samples = div(dat.sample_rate, step_size)
    eog_diff = @views diff(dat.data[1:step_samples:end, channel_in])
    eog_idx = findall(x -> abs(x) >= criterion, eog_diff)
    eog_idx = [idx for (i, idx) in enumerate(eog_idx) if i == 1 || (idx - eog_idx[i-1] > 2)] .* step_samples
    filter!(idx -> idx <= n_samples(dat), eog_idx)
    dat.data[!, channel_out] = falses(n_samples(dat))
    dat.data[eog_idx, channel_out] .= true
    return dat
end

function detect_eog_onsets!(dat::ContinuousData, eog_cfg::Dict)
    v_ch = eog_cfg["vEOG_channels"][3][1]
    detect_eog_onsets!(dat, eog_cfg["vEOG_criterion"], Symbol(v_ch), Symbol("is_", v_ch))
    h_ch = eog_cfg["hEOG_channels"][3][1]
    detect_eog_onsets!(dat, eog_cfg["hEOG_criterion"], Symbol(h_ch), Symbol("is_", h_ch))
    return dat
end

"""
    detect_eog_onsets(dat, args...; kwargs...)

Non-mutating version of [`detect_eog_onsets!`](@ref).
"""
detect_eog_onsets(dat::ContinuousData, args...; kwargs...) = detect_eog_onsets!(copy(dat), args...; kwargs...)

function detect_eog_onsets!(dat_vec::AbstractVector{<:ContinuousData}, args...; kwargs...)
    for d in dat_vec
        detect_eog_onsets!(d, args...; kwargs...)
    end
    return dat_vec
end

detect_eog_onsets(dat_vec::AbstractVector{<:ContinuousData}, args...; kwargs...) =
    [detect_eog_onsets(d, args...; kwargs...) for d in dat_vec]
