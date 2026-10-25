# === PEAK TO PEAK ARTIFACT DETECTION ===

"""
    _is_peak_to_peak(signal::AbstractVector{<:Real}, threshold::Real, window_samples::Int)

Detect moving-window peak-to-peak artifacts (max - min > threshold).
Returns a boolean mask where all samples within any offending window are true.
"""
function _is_peak_to_peak(signal::AbstractVector{<:Real}, threshold::Real, window_samples::Int)
    n = length(signal)
    n == 0 && return falses(0)
    mask = falses(n)
    w = max(1, window_samples)
    if w >= n
        return (maximum(signal) - minimum(signal) > threshold ? trues(n) : falses(n))
    end
    if w == 1
        return falses(n)
    end

    n_blocks = cld(n, w)
    pref_max = Vector{Float64}(undef, w)
    pref_min = Vector{Float64}(undef, w)
    suff_max_prev = Vector{Float64}(undef, w)
    suff_min_prev = Vector{Float64}(undef, w)
    suff_max_curr = Vector{Float64}(undef, w)
    suff_min_curr = Vector{Float64}(undef, w)

    last_marked = 0
    thresh_f = Float64(threshold)

    @inbounds for b = 1:n_blocks
        b_start = (b - 1) * w + 1
        b_end = min(b_start + w - 1, n)
        b_len = b_end - b_start + 1

        p_max = Float64(signal[b_start])
        p_min = Float64(signal[b_start])
        pref_max[1] = p_max
        pref_min[1] = p_min
        for k = 2:b_len
            s_val = Float64(signal[b_start+k-1])
            p_max = max(p_max, s_val)
            p_min = min(p_min, s_val)
            pref_max[k] = p_max
            pref_min[k] = p_min
        end

        s_max = Float64(signal[b_end])
        s_min = Float64(signal[b_end])
        suff_max_curr[b_len] = s_max
        suff_min_curr[b_len] = s_min
        for k = (b_len-1):-1:1
            s_val = Float64(signal[b_start+k-1])
            s_max = max(s_max, s_val)
            s_min = min(s_min, s_val)
            suff_max_curr[k] = s_max
            suff_min_curr[k] = s_min
        end

        if b == 1
            if b_len == w
                w_max = pref_max[w]
                w_min = pref_min[w]
                if (w_max - w_min) > thresh_f
                    mask[1:w] .= true
                    last_marked = w
                end
            end
        else
            for k = 2:w
                j_global = (b - 1) * w + k - 1
                if j_global > n
                    break
                end
                i_global = j_global - w + 1
                w_max = max(suff_max_prev[k], pref_max[k-1])
                w_min = min(suff_min_prev[k], pref_min[k-1])
                if (w_max - w_min) > thresh_f
                    s_idx = max(last_marked + 1, i_global)
                    if s_idx <= j_global
                        mask[s_idx:j_global] .= true
                        last_marked = j_global
                    end
                end
            end
            if b_len == w
                i_global = b_start
                j_global = b_end
                w_max = pref_max[w]
                w_min = pref_min[w]
                if (w_max - w_min) > thresh_f
                    s_idx = max(last_marked + 1, i_global)
                    if s_idx <= j_global
                        mask[s_idx:j_global] .= true
                        last_marked = j_global
                    end
                end
            end
        end

        suff_max_prev, suff_max_curr = suff_max_curr, suff_max_prev
        suff_min_prev, suff_min_curr = suff_min_curr, suff_min_prev
    end

    return mask
end

"""
    is_peak_to_peak!(dat::SingleDataFrameEeg, threshold::Real, window_size::Real;
                     channel_selection::Function = channels(),
                     sample_selection::Function = samples(),
                     interval_selection::Interval = times(),
                     mode::Symbol = :combined, channel_out = nothing)
    is_peak_to_peak!(dat::MultiDataFrameEeg, ...)
    is_peak_to_peak!(dat::Vector{EpochData}, ...)

Detect moving window peak-to-peak artifacts (where max-min > `threshold` inside `window_size` in seconds).
"""
function is_peak_to_peak!(
    dat::SingleDataFrameEeg,
    threshold::Real,
    window_size::Real;
    channel_selection = channels(),
    sample_selection = samples(),
    interval_selection::Interval = times(),
    mode::Symbol = :combined,
    channel_out::Union{Symbol,Nothing} = nothing,
)
    mode ∉ [:separate, :combined] && @minimal_error("mode must be :separate or :combined")
    threshold <= 0 && @minimal_error("threshold must be greater than 0")
    window_size <= 0 && @minimal_error("window_size must be greater than 0")

    selected_channels = get_selected_channels(dat, channel_selection; include_meta = false, include_extra = false)
    isempty(selected_channels) && @minimal_error("No channels selected")

    combined_sel = _combine_interval_sample(interval_selection, sample_selection)
    selected_samples = get_selected_samples(dat.data, combined_sel)

    window_samples = round(Int, window_size * dat.sample_rate)
    channel_out = something(channel_out, :is_peak_to_peak)

    if mode == :combined
        dat.data[!, channel_out] = falses(n_samples(dat))
        sample_mask = falses(n_samples(dat))
        sample_mask[selected_samples] .= true

        for ch in selected_channels
            p2p_mask = _is_peak_to_peak(dat.data[!, ch], threshold, window_samples) .& sample_mask
            dat.data[!, channel_out] .|= p2p_mask
        end
    elseif mode == :separate
        sample_mask = falses(n_samples(dat))
        sample_mask[selected_samples] .= true

        for ch in selected_channels
            p2p_mask = _is_peak_to_peak(dat.data[!, ch], threshold, window_samples) .& sample_mask
            column_name = Symbol("$(channel_out)_$(ch)")
            dat.data[!, column_name] = p2p_mask
        end
    end

    return dat
end

function is_peak_to_peak!(
    dat::MultiDataFrameEeg,
    threshold::Real,
    window_size::Real;
    channel_selection = channels(),
    sample_selection = samples(),
    interval_selection::Interval = times(),
    epoch_selection = epochs(),
    mode::Symbol = :combined,
    channel_out::Union{Symbol,Nothing} = nothing,
)
    mode ∉ [:separate, :combined] && @minimal_error("mode must be :separate or :combined")
    threshold <= 0 && @minimal_error("threshold must be greater than 0")
    window_size <= 0 && @minimal_error("window_size must be greater than 0")

    selected_channels = get_selected_channels(dat, channel_selection; include_meta = false, include_extra = false)
    isempty(selected_channels) && @minimal_error("No channels selected")
    selected_epochs = get_selected_epochs(dat, epoch_selection)
    isempty(selected_epochs) && @minimal_error("No epochs selected")

    window_samples = round(Int, window_size * dat.sample_rate)
    channel_out = something(channel_out, :is_peak_to_peak)

    # Combine interval and sample selection once (constant across all epochs)
    combined_sel = _combine_interval_sample(interval_selection, sample_selection)

    Threads.@threads for epoch_idx in selected_epochs
        epoch_df = dat.data[epoch_idx]
        selected_samples = get_selected_samples(epoch_df, combined_sel)

        epoch_df[!, channel_out] = falses(nrow(epoch_df))
        sample_mask = falses(nrow(epoch_df))
        sample_mask[selected_samples] .= true

        if mode == :combined
            for ch in selected_channels
                p2p_mask = _is_peak_to_peak(epoch_df[!, ch], threshold, window_samples) .& sample_mask
                epoch_df[!, channel_out] .|= p2p_mask
            end
        else
            for ch in selected_channels
                p2p_mask = _is_peak_to_peak(epoch_df[!, ch], threshold, window_samples) .& sample_mask
                column_name = Symbol("$(channel_out)_$(ch)")
                epoch_df[!, column_name] = p2p_mask
            end
        end
    end
    return dat
end

function is_peak_to_peak!(dat::AbstractVector{<:EpochData}, threshold::Real, window_size::Real; kwargs...)
    for d in dat
        is_peak_to_peak!(d, threshold, window_size; kwargs...)
    end
    return dat
end

"""
    is_peak_to_peak(dat::SingleDataFrameEeg, threshold::Real, window_size::Real;
                    channel_selection = channels(), sample_selection = samples(),
                    interval_selection = times(), mode::Symbol = :combined)

Detect moving window peak-to-peak artifacts without modifying data.
Returns `Vector{Bool}` (`:combined` mode) or a `DataFrame` (`:separate` mode).
"""
function is_peak_to_peak(
    dat::SingleDataFrameEeg,
    threshold::Real,
    window_size::Real;
    channel_selection = channels(),
    sample_selection = samples(),
    interval_selection::Interval = times(),
    mode::Symbol = :combined,
)
    mode ∉ [:separate, :combined] && @minimal_error("mode must be :separate or :combined")
    temp_dat = copy(dat)
    is_peak_to_peak!(temp_dat, threshold, window_size; channel_selection, sample_selection, interval_selection, mode)

    if mode == :combined
        return temp_dat.data[!, :is_peak_to_peak]
    elseif mode == :separate
        selected_channels = get_selected_channels(dat, channel_selection, include_meta = false, include_extra = false)
        cols = [Symbol("is_peak_to_peak_$(ch)") for ch in selected_channels]
        return temp_dat.data[!, cols]
    end
end

"""
    n_peak_to_peak(dat::SingleDataFrameEeg, threshold::Real, window_size::Real; ...)

Count the number of peak-to-peak artifact values across selected channels.
"""
function n_peak_to_peak(
    dat::SingleDataFrameEeg,
    threshold::Real,
    window_size::Real;
    channel_selection = channels(),
    sample_selection = samples(),
    interval_selection::Interval = times(),
    mode::Symbol = :combined,
)
    mode ∉ [:separate, :combined] && @minimal_error("mode must be :separate or :combined")
    selected_channels = get_selected_channels(dat, channel_selection; include_meta = false, include_extra = false)
    combined_sel = _combine_interval_sample(interval_selection, sample_selection)
    selected_samples = get_selected_samples(dat.data, combined_sel)
    sample_mask = falses(n_samples(dat))
    sample_mask[selected_samples] .= true
    window_samples = round(Int, window_size * dat.sample_rate)

    if mode == :combined
        combined_mask = falses(n_samples(dat))
        for ch in selected_channels
            p2p_mask = _is_peak_to_peak(dat.data[!, ch], threshold, window_samples) .& sample_mask
            combined_mask .|= p2p_mask
        end
        return sum(combined_mask)
    elseif mode == :separate
        counts = [sum(_is_peak_to_peak(dat.data[!, ch], threshold, window_samples) .& sample_mask) for ch in selected_channels]
        return DataFrame(channel = selected_channels, n_peak_to_peak = counts)
    end
end
