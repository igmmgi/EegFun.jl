# === FLATLINE (LOW VARIANCE) DETECTION ===

"""
    _is_flatline(signal::AbstractVector{<:Real}, threshold::Real, window_samples::Int)

Detect flatline artifacts (signal stuck with extremely low variance) using a sliding window.
Returns a boolean mask where all samples within any window whose std dev is below the threshold are true.
"""
function _is_flatline(signal::AbstractVector{<:Real}, threshold::Real, window_samples::Int)
    n = length(signal)
    n == 0 && return falses(0)
    mask = falses(n)
    w = max(1, window_samples)
    if w >= n
        return (std(signal) < threshold ? trues(n) : falses(n))
    end
    if w <= 1
        return falses(n)
    end

    # Initial window 1:w
    w_view = view(signal, 1:w)
    μ = mean(w_view)
    M2 = var(w_view) * (w - 1)

    inv_w = 1.0 / w
    inv_w_minus_1 = 1.0 / (w - 1)
    thresh_sq = Float64(threshold)^2

    last_marked = 0
    if (M2 * inv_w_minus_1) < thresh_sq
        mask[1:w] .= true
        last_marked = w
    end

    anchor_interval = 256

    @inbounds for i = 2:(n - w + 1)
        x_out = Float64(signal[i - 1])
        x_in  = Float64(signal[i + w - 1])

        if (i % anchor_interval) == 0
            cur_view = view(signal, i:(i+w-1))
            μ = mean(cur_view)
            M2 = var(cur_view) * (w - 1)
        else
            d = x_in - x_out
            old_μ = μ
            μ = old_μ + d * inv_w
            M2 = max(0.0, M2 + d * ((x_in - μ) + (x_out - old_μ)))
        end

        cur_var = M2 * inv_w_minus_1
        if cur_var < thresh_sq
            s_idx = max(last_marked + 1, i)
            e_idx = i + w - 1
            if s_idx <= e_idx
                mask[s_idx:e_idx] .= true
                last_marked = e_idx
            end
        end
    end

    return mask
end

"""
    is_flatline!(dat::SingleDataFrameEeg, threshold::Real, window_size::Real;
                 channel_selection::Function = channels(),
                 sample_selection::Function = samples(),
                 interval_selection::Interval = times(),
                 mode::Symbol = :combined, channel_out = nothing)
    is_flatline!(dat::MultiDataFrameEeg, ...)
    is_flatline!(dat::Vector{EpochData}, ...)

Detect flatline artifacts using a sliding window of `window_size` (in seconds).
Flags all samples in windows where standard deviation < `threshold`.
"""
function is_flatline!(
    dat::SingleDataFrameEeg,
    threshold::Real,
    window_size::Real;
    channel_selection::Function = channels(),
    sample_selection::Function = samples(),
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
    channel_out = something(channel_out, :is_flatline)

    if mode == :combined
        dat.data[!, channel_out] = falses(n_samples(dat))
        sample_mask = falses(n_samples(dat))
        sample_mask[selected_samples] .= true

        for ch in selected_channels
            flat_mask = _is_flatline(dat.data[!, ch], threshold, window_samples) .& sample_mask
            dat.data[!, channel_out] .|= flat_mask
        end
    elseif mode == :separate
        sample_mask = falses(n_samples(dat))
        sample_mask[selected_samples] .= true

        for ch in selected_channels
            flat_mask = _is_flatline(dat.data[!, ch], threshold, window_samples) .& sample_mask
            column_name = Symbol("$(channel_out)_$(ch)")
            dat.data[!, column_name] = flat_mask
        end
    end

    return dat
end

function is_flatline!(
    dat::MultiDataFrameEeg,
    threshold::Real,
    window_size::Real;
    channel_selection::Function = channels(),
    sample_selection::Function = samples(),
    interval_selection::Interval = times(),
    epoch_selection::Function = epochs(),
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
    channel_out = something(channel_out, :is_flatline)

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
                flat_mask = _is_flatline(epoch_df[!, ch], threshold, window_samples) .& sample_mask
                epoch_df[!, channel_out] .|= flat_mask
            end
        else
            for ch in selected_channels
                flat_mask = _is_flatline(epoch_df[!, ch], threshold, window_samples) .& sample_mask
                column_name = Symbol("$(channel_out)_$(ch)")
                epoch_df[!, column_name] = flat_mask
            end
        end
    end
    return dat
end

function is_flatline!(dat::Vector{EpochData}, threshold::Real, window_size::Real; kwargs...)
    for d in dat
        is_flatline!(d, threshold, window_size; kwargs...)
    end
    return dat
end

"""
    is_flatline(dat::SingleDataFrameEeg, threshold::Real, window_size::Real;
                channel_selection = channels(), sample_selection = samples(),
                interval_selection = times(), mode::Symbol = :combined)

Detect flatline artifacts without modifying data. Returns `Vector{Bool}` (`:combined` mode)
or a `DataFrame` (`:separate` mode).
"""
function is_flatline(
    dat::SingleDataFrameEeg,
    threshold::Real,
    window_size::Real;
    channel_selection::Function = channels(),
    sample_selection::Function = samples(),
    interval_selection::Interval = times(),
    mode::Symbol = :combined,
)
    mode ∉ [:separate, :combined] && @minimal_error("mode must be :separate or :combined")
    temp_dat = copy(dat)
    is_flatline!(temp_dat, threshold, window_size; channel_selection, sample_selection, interval_selection, mode)

    if mode == :combined
        return temp_dat.data[!, :is_flatline]
    elseif mode == :separate
        selected_channels = get_selected_channels(dat, channel_selection, include_meta = false, include_extra = false)
        cols = [Symbol("is_flatline_$(ch)") for ch in selected_channels]
        return temp_dat.data[!, cols]
    end
end

"""
    n_flatline(dat::SingleDataFrameEeg, threshold::Real, window_size::Real; ...)

Count the number of flatline values across selected channels.
"""
function n_flatline(
    dat::SingleDataFrameEeg,
    threshold::Real,
    window_size::Real;
    channel_selection::Function = channels(),
    sample_selection::Function = samples(),
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
            flat_mask = _is_flatline(dat.data[!, ch], threshold, window_samples) .& sample_mask
            combined_mask .|= flat_mask
        end
        return sum(combined_mask)
    elseif mode == :separate
        counts = [sum(_is_flatline(dat.data[!, ch], threshold, window_samples) .& sample_mask) for ch in selected_channels]
        return DataFrame(channel = selected_channels, n_flatline = counts)
    end
end

