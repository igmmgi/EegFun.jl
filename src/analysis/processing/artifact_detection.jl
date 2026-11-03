# Detection, reporting, and repair of artifacts in continuous and epoched EEG data.

function _validate_detector_args(mode::Symbol, threshold::Real, window_size::Union{Nothing,Real} = nothing)
    mode in (:separate, :combined) || @minimal_error("mode must be :separate or :combined")
    threshold > 0 || @minimal_error("threshold must be greater than 0")
    !isnothing(window_size) && window_size <= 0 && @minimal_error("window_size must be greater than 0")
    return nothing
end

function _selected_sample_mask(df::DataFrame, combined_sel)
    n = nrow(df)
    sel = get_selected_samples(df, combined_sel)
    length(sel) == n && return nothing
    mask = falses(n)
    mask[sel] .= true
    return mask
end

@inline function _apply_sample_mask!(buf::AbstractVector{Bool}, sample_mask::Union{Nothing,BitVector})
    isnothing(sample_mask) || (buf .&= sample_mask)
    return buf
end

function _flag_artifacts!(
    detect_into!::F,
    dat::SingleDataFrameEeg,
    combined_col::Symbol,
    sep_col_fn::G;
    channel_selection = channels(),
    sample_selection = samples(),
    interval_selection::Interval = times(),
    mode::Symbol = :combined,
) where {F,G}
    selected_channels = get_selected_channels(dat, channel_selection; include_meta = false, include_extra = false)
    isempty(selected_channels) && @minimal_error("No channels selected")

    n = n_samples(dat)
    sample_mask = _selected_sample_mask(dat.data, _combine_interval_sample(interval_selection, sample_selection))

    if mode == :combined
        out = falses(n)
        buf = falses(n)
        for ch in selected_channels
            detect_into!(buf, dat.data[!, ch])
            _apply_sample_mask!(buf, sample_mask)
            out .|= buf
        end
        dat.data[!, combined_col] = out
    else
        for ch in selected_channels
            col_mask = falses(n)
            detect_into!(col_mask, dat.data[!, ch])
            _apply_sample_mask!(col_mask, sample_mask)
            dat.data[!, sep_col_fn(ch)] = col_mask
        end
    end
    return dat
end

function _flag_artifacts!(
    detect_into!::F,
    dat::MultiDataFrameEeg,
    combined_col::Symbol,
    sep_col_fn::G;
    channel_selection = channels(),
    sample_selection = samples(),
    interval_selection::Interval = times(),
    epoch_selection = epochs(),
    mode::Symbol = :combined,
) where {F,G}
    selected_channels = get_selected_channels(dat, channel_selection; include_meta = false, include_extra = false)
    isempty(selected_channels) && @minimal_error("No channels selected")
    selected_epochs = get_selected_epochs(dat, epoch_selection)
    isempty(selected_epochs) && @minimal_error("No epochs selected")

    out_cols = mode == :combined ? [combined_col] : [sep_col_fn(ch) for ch in selected_channels]
    if length(selected_epochs) < length(dat.data)
        selected_set = Set(selected_epochs)
        for (idx, epoch_df) in enumerate(dat.data)
            idx in selected_set && continue
            for col in out_cols
                hasproperty(epoch_df, col) || (epoch_df[!, col] = falses(nrow(epoch_df)))
            end
        end
    end

    combined_sel = _combine_interval_sample(interval_selection, sample_selection)
    Threads.@threads for epoch_idx in selected_epochs
        epoch_df = dat.data[epoch_idx]
        n = nrow(epoch_df)
        sample_mask = _selected_sample_mask(epoch_df, combined_sel)
        if mode == :combined
            out = falses(n)
            buf = falses(n)
            for ch in selected_channels
                detect_into!(buf, epoch_df[!, ch])
                _apply_sample_mask!(buf, sample_mask)
                out .|= buf
            end
            epoch_df[!, combined_col] = out
        else
            for (ci, ch) in enumerate(selected_channels)
                col_mask = falses(n)
                detect_into!(col_mask, epoch_df[!, ch])
                _apply_sample_mask!(col_mask, sample_mask)
                epoch_df[!, out_cols[ci]] = col_mask
            end
        end
    end
    return dat
end

function _detect_artifacts(
    detect_into!::F,
    dat::SingleDataFrameEeg,
    sep_col_fn::G;
    channel_selection = channels(),
    sample_selection = samples(),
    interval_selection::Interval = times(),
    mode::Symbol = :combined,
) where {F,G}
    selected_channels = get_selected_channels(dat, channel_selection; include_meta = false, include_extra = false)
    isempty(selected_channels) && @minimal_error("No channels selected")

    n = n_samples(dat)
    sample_mask = _selected_sample_mask(dat.data, _combine_interval_sample(interval_selection, sample_selection))

    if mode == :combined
        out = zeros(Bool, n)
        buf = falses(n)
        for ch in selected_channels
            detect_into!(buf, dat.data[!, ch])
            _apply_sample_mask!(buf, sample_mask)
            out .|= buf
        end
        return out
    else
        pairs = Pair{Symbol,BitVector}[]
        sizehint!(pairs, length(selected_channels))
        for ch in selected_channels
            col_mask = falses(n)
            detect_into!(col_mask, dat.data[!, ch])
            _apply_sample_mask!(col_mask, sample_mask)
            push!(pairs, sep_col_fn(ch) => col_mask)
        end
        return DataFrame(pairs; copycols = false)
    end
end

function _detect_artifacts(
    detect_into!::F,
    dat::MultiDataFrameEeg,
    sep_col_fn::G;
    channel_selection = channels(),
    sample_selection = samples(),
    interval_selection::Interval = times(),
    epoch_selection = epochs(),
    mode::Symbol = :combined,
) where {F,G}
    selected_channels = get_selected_channels(dat, channel_selection; include_meta = false, include_extra = false)
    isempty(selected_channels) && @minimal_error("No channels selected")
    selected_epochs = get_selected_epochs(dat, epoch_selection)
    isempty(selected_epochs) && @minimal_error("No epochs selected")

    n_ep = length(dat.data)
    combined_sel = _combine_interval_sample(interval_selection, sample_selection)

    if mode == :combined
        results = [zeros(Bool, nrow(ep_df)) for ep_df in dat.data]
        Threads.@threads for epoch_idx in selected_epochs
            epoch_df = dat.data[epoch_idx]
            out = results[epoch_idx]
            buf = falses(nrow(epoch_df))
            sample_mask = _selected_sample_mask(epoch_df, combined_sel)
            for ch in selected_channels
                detect_into!(buf, epoch_df[!, ch])
                _apply_sample_mask!(buf, sample_mask)
                out .|= buf
            end
        end
        return results
    else
        sep_cols = [sep_col_fn(ch) for ch in selected_channels]
        results = Vector{DataFrame}(undef, n_ep)
        if length(selected_epochs) < n_ep
            selected_set = Set(selected_epochs)
            for idx = 1:n_ep
                idx in selected_set && continue
                n = nrow(dat.data[idx])
                results[idx] = DataFrame([col => falses(n) for col in sep_cols]; copycols = false)
            end
        end
        Threads.@threads for epoch_idx in selected_epochs
            epoch_df = dat.data[epoch_idx]
            n = nrow(epoch_df)
            sample_mask = _selected_sample_mask(epoch_df, combined_sel)
            pairs = Vector{Pair{Symbol,BitVector}}(undef, length(selected_channels))
            for (ci, ch) in enumerate(selected_channels)
                col_mask = falses(n)
                detect_into!(col_mask, epoch_df[!, ch])
                _apply_sample_mask!(col_mask, sample_mask)
                pairs[ci] = sep_cols[ci] => col_mask
            end
            results[epoch_idx] = DataFrame(pairs; copycols = false)
        end
        return results
    end
end

function _count_artifacts(
    detect_into!::F,
    dat::SingleDataFrameEeg,
    count_col::Symbol;
    channel_selection = channels(),
    sample_selection = samples(),
    interval_selection::Interval = times(),
    mode::Symbol = :combined,
) where {F}
    selected_channels = get_selected_channels(dat, channel_selection; include_meta = false, include_extra = false)
    isempty(selected_channels) && @minimal_error("No channels selected")

    n = n_samples(dat)
    sample_mask = _selected_sample_mask(dat.data, _combine_interval_sample(interval_selection, sample_selection))
    buf = falses(n)

    if mode == :combined
        out = falses(n)
        for ch in selected_channels
            detect_into!(buf, dat.data[!, ch])
            _apply_sample_mask!(buf, sample_mask)
            out .|= buf
        end
        return count(out)
    else
        counts = Vector{Int}(undef, length(selected_channels))
        for (ci, ch) in enumerate(selected_channels)
            detect_into!(buf, dat.data[!, ch])
            _apply_sample_mask!(buf, sample_mask)
            counts[ci] = count(buf)
        end
        return DataFrame(:channel => selected_channels, count_col => counts)
    end
end

function _count_artifacts(
    detect_into!::F,
    dat::MultiDataFrameEeg,
    count_col::Symbol;
    channel_selection = channels(),
    sample_selection = samples(),
    interval_selection::Interval = times(),
    epoch_selection = epochs(),
    mode::Symbol = :combined,
) where {F}
    selected_channels = get_selected_channels(dat, channel_selection; include_meta = false, include_extra = false)
    isempty(selected_channels) && @minimal_error("No channels selected")
    selected_epochs = get_selected_epochs(dat, epoch_selection)
    isempty(selected_epochs) && @minimal_error("No epochs selected")

    combined_sel = _combine_interval_sample(interval_selection, sample_selection)

    if mode == :combined
        total = 0
        for epoch_idx in selected_epochs
            epoch_df = dat.data[epoch_idx]
            n = nrow(epoch_df)
            sample_mask = _selected_sample_mask(epoch_df, combined_sel)
            out = falses(n)
            buf = falses(n)
            for ch in selected_channels
                detect_into!(buf, epoch_df[!, ch])
                _apply_sample_mask!(buf, sample_mask)
                out .|= buf
            end
            total += count(out)
        end
        return total
    else
        counts = zeros(Int, length(selected_channels))
        for epoch_idx in selected_epochs
            epoch_df = dat.data[epoch_idx]
            buf = falses(nrow(epoch_df))
            sample_mask = _selected_sample_mask(epoch_df, combined_sel)
            for (ci, ch) in enumerate(selected_channels)
                detect_into!(buf, epoch_df[!, ch])
                _apply_sample_mask!(buf, sample_mask)
                counts[ci] += count(buf)
            end
        end
        return DataFrame(:channel => selected_channels, count_col => counts)
    end
end

include("artifact_detection/eog_onsets.jl")
include("artifact_detection/extreme.jl")
include("artifact_detection/step.jl")
include("artifact_detection/flatline.jl")
include("artifact_detection/peak_to_peak.jl")
include("artifact_detection/bridged.jl")
include("artifact_detection/epoch_rejection_types.jl")
include("artifact_detection/epoch_detection.jl")
include("artifact_detection/repair.jl")
include("artifact_detection/subset_bad_data.jl")
include("artifact_detection/autoreject.jl")
