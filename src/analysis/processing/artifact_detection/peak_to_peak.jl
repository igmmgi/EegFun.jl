# === PEAK TO PEAK ARTIFACT DETECTION ===

"""
    _is_peak_to_peak(signal::AbstractVector{<:Real}, threshold::Real, window_samples::Int)

Detect moving-window peak-to-peak artifacts (max - min > threshold).
Returns a boolean mask where all samples within any offending window are true.
"""
function _is_peak_to_peak(signal::AbstractVector{<:Real}, threshold::Real, window_samples::Int)
    n = length(signal)
    mask = falses(n)
    w = max(1, window_samples)

    @inbounds for i = 1:(n-w+1)
        w_view = view(signal, i:(i+w-1))
        if (maximum(w_view) - minimum(w_view)) > threshold
            mask[i:(i+w-1)] .= true
        end
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

    return nothing
end

function is_peak_to_peak!(
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
    return nothing
end

is_peak_to_peak!(dat::Vector{EpochData}, threshold::Real, window_size::Real; kwargs...) =
    is_peak_to_peak!.(dat, threshold, window_size; kwargs...)

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
    channel_selection::Function = channels(),
    sample_selection::Function = samples(),
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
            p2p_mask = _is_peak_to_peak(dat.data[!, ch], threshold, window_samples) .& sample_mask
            combined_mask .|= p2p_mask
        end
        return sum(combined_mask)
    elseif mode == :separate
        counts = [sum(_is_peak_to_peak(dat.data[!, ch], threshold, window_samples) .& sample_mask) for ch in selected_channels]
        return DataFrame(channel = selected_channels, n_peak_to_peak = counts)
    end
end

