# Detection, reporting, and repair of artifacts in continuous and epoched EEG data.

######################################
######### Types & Constants ##########
######################################

const _DetectableEeg = Union{SingleDataFrameEeg,MultiDataFrameEeg}

######################################
######### Internal Functions #########
######################################

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

function _detect_df_combined!(
    out::AbstractVector{Bool},
    buf::BitVector,
    detect_into!::D,
    df::DataFrame,
    selected_channels::Vector{Symbol},
    sample_mask::Union{Nothing,BitVector},
) where {D}
    fill!(out, false)
    resize!(buf, nrow(df))
    for ch in selected_channels
        detect_into!(buf, df[!, ch])
        _apply_sample_mask!(buf, sample_mask)
        out .|= buf
    end
    return out
end

function _detect_df_separate(
    detect_into!::D,
    df::DataFrame,
    selected_channels::Vector{Symbol},
    sep_cols::Vector{Symbol},
    sample_mask::Union{Nothing,BitVector},
) where {D}
    n = nrow(df)
    pairs = Vector{Pair{Symbol,BitVector}}(undef, length(selected_channels))
    for (ci, ch) in enumerate(selected_channels)
        col_mask = falses(n)
        detect_into!(col_mask, df[!, ch])
        _apply_sample_mask!(col_mask, sample_mask)
        pairs[ci] = sep_cols[ci] => col_mask
    end
    return pairs
end

function _flag_artifacts!(
    make_detector::F,
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

    detect_into! = make_detector()
    n = n_samples(dat)
    sample_mask = _selected_sample_mask(dat.data, _combine_interval_sample(interval_selection, sample_selection))

    if mode == :combined
        dat.data[!, combined_col] = _detect_df_combined!(falses(n), falses(n), detect_into!, dat.data, selected_channels, sample_mask)
    else
        sep_cols = Symbol[sep_col_fn(ch) for ch in selected_channels]
        for (col, col_mask) in _detect_df_separate(detect_into!, dat.data, selected_channels, sep_cols, sample_mask)
            dat.data[!, col] = col_mask
        end
    end
    return dat
end

function _flag_artifacts!(
    make_detector::F,
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

    out_cols = mode == :combined ? Symbol[combined_col] : Symbol[sep_col_fn(ch) for ch in selected_channels]
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
    _foreach_chunk(length(selected_epochs)) do chunk
        detect_into! = make_detector()
        buf = BitVector()
        for k in chunk
            epoch_df = dat.data[selected_epochs[k]]
            sample_mask = _selected_sample_mask(epoch_df, combined_sel)
            if mode == :combined
                epoch_df[!, combined_col] =
                    _detect_df_combined!(falses(nrow(epoch_df)), buf, detect_into!, epoch_df, selected_channels, sample_mask)
            else
                for (col, col_mask) in _detect_df_separate(detect_into!, epoch_df, selected_channels, out_cols, sample_mask)
                    epoch_df[!, col] = col_mask
                end
            end
        end
    end
    return dat
end

function _detect_artifacts(
    make_detector::F,
    dat::SingleDataFrameEeg,
    sep_col_fn::G;
    channel_selection = channels(),
    sample_selection = samples(),
    interval_selection::Interval = times(),
    mode::Symbol = :combined,
) where {F,G}
    selected_channels = get_selected_channels(dat, channel_selection; include_meta = false, include_extra = false)
    isempty(selected_channels) && @minimal_error("No channels selected")

    detect_into! = make_detector()
    n = n_samples(dat)
    sample_mask = _selected_sample_mask(dat.data, _combine_interval_sample(interval_selection, sample_selection))

    if mode == :combined
        return _detect_df_combined!(zeros(Bool, n), falses(n), detect_into!, dat.data, selected_channels, sample_mask)
    else
        sep_cols = Symbol[sep_col_fn(ch) for ch in selected_channels]
        return DataFrame(_detect_df_separate(detect_into!, dat.data, selected_channels, sep_cols, sample_mask); copycols = false)
    end
end

function _detect_artifacts(
    make_detector::F,
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
        _foreach_chunk(length(selected_epochs)) do chunk
            detect_into! = make_detector()
            buf = BitVector()
            for k in chunk
                epoch_idx = selected_epochs[k]
                epoch_df = dat.data[epoch_idx]
                sample_mask = _selected_sample_mask(epoch_df, combined_sel)
                _detect_df_combined!(results[epoch_idx], buf, detect_into!, epoch_df, selected_channels, sample_mask)
            end
        end
        return results
    else
        sep_cols = Symbol[sep_col_fn(ch) for ch in selected_channels]
        results = Vector{DataFrame}(undef, n_ep)
        if length(selected_epochs) < n_ep
            selected_set = Set(selected_epochs)
            for idx = 1:n_ep
                idx in selected_set && continue
                n = nrow(dat.data[idx])
                results[idx] = DataFrame([col => falses(n) for col in sep_cols]; copycols = false)
            end
        end
        _foreach_chunk(length(selected_epochs)) do chunk
            detect_into! = make_detector()
            for k in chunk
                epoch_idx = selected_epochs[k]
                epoch_df = dat.data[epoch_idx]
                sample_mask = _selected_sample_mask(epoch_df, combined_sel)
                results[epoch_idx] =
                    DataFrame(_detect_df_separate(detect_into!, epoch_df, selected_channels, sep_cols, sample_mask); copycols = false)
            end
        end
        return results
    end
end

function _count_artifacts_dfs(
    detect_into!::D,
    dfs,
    selected_channels::Vector{Symbol},
    combined_sel,
    count_col::Symbol,
    mode::Symbol,
) where {D}
    buf = BitVector()
    if mode == :combined
        out = BitVector()
        total = 0
        for df in dfs
            resize!(out, nrow(df))
            sample_mask = _selected_sample_mask(df, combined_sel)
            _detect_df_combined!(out, buf, detect_into!, df, selected_channels, sample_mask)
            total += count(out)
        end
        return total
    else
        counts = zeros(Int, length(selected_channels))
        for df in dfs
            resize!(buf, nrow(df))
            sample_mask = _selected_sample_mask(df, combined_sel)
            for (ci, ch) in enumerate(selected_channels)
                detect_into!(buf, df[!, ch])
                _apply_sample_mask!(buf, sample_mask)
                counts[ci] += count(buf)
            end
        end
        return DataFrame(:channel => selected_channels, count_col => counts)
    end
end

function _count_artifacts(
    make_detector::F,
    dat::SingleDataFrameEeg,
    count_col::Symbol;
    channel_selection = channels(),
    sample_selection = samples(),
    interval_selection::Interval = times(),
    mode::Symbol = :combined,
) where {F}
    selected_channels = get_selected_channels(dat, channel_selection; include_meta = false, include_extra = false)
    isempty(selected_channels) && @minimal_error("No channels selected")
    combined_sel = _combine_interval_sample(interval_selection, sample_selection)
    return _count_artifacts_dfs(make_detector(), (dat.data,), selected_channels, combined_sel, count_col, mode)
end

function _count_artifacts(
    make_detector::F,
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
    return _count_artifacts_dfs(make_detector(), @view(dat.data[selected_epochs]), selected_channels, combined_sel, count_col, mode)
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
