# Multivariate pattern classification (MVPA / decoding) for ERP and time-frequency epochs.

function _copy_selected_col!(dest::Array{Float64,3}, ch_idx::Int, trial_idx::Int, col::Vector{Float64}, samples::UnitRange{Int})
    offset = first(samples) - 1
    @inbounds @simd for t = 1:length(samples)
        dest[ch_idx, t, trial_idx] = col[offset+t]
    end
    return nothing
end

function _copy_selected_col!(dest::Array{Float64,3}, ch_idx::Int, trial_idx::Int, col::Vector{Float64}, samples::AbstractVector{Int})
    @inbounds @simd for t in eachindex(samples)
        dest[ch_idx, t, trial_idx] = col[samples[t]]
    end
    return nothing
end

function _copy_selected_col!(dest::Array{Float64,3}, ch_idx::Int, trial_idx::Int, col::AbstractVector, samples::AbstractVector{Int})
    @inbounds for (t, s) in enumerate(samples)
        dest[ch_idx, t, trial_idx] = Float64(col[s])
    end
    return nothing
end

function _extract_selected_f64(col::Vector{Float64}, rows::AbstractVector{Int})
    out = Vector{Float64}(undef, length(rows))
    @inbounds @simd for i in eachindex(rows)
        out[i] = col[rows[i]]
    end
    return out
end

function _extract_selected_f64(col::AbstractVector, rows::AbstractVector{Int})
    return Float64[Float64(col[r]) for r in rows]
end

function _as_contiguous_range(rows::AbstractVector{Int}, n_total::Int)
    n = length(rows)
    if n > 0 && last(rows) - first(rows) == n - 1 && issorted(rows)
        return first(rows):last(rows)
    end
    n == 0 && n_total == 0 && return 1:0
    return rows
end

function _extract_epoch_arrays(
    epochs::AbstractVector{<:EpochData},
    selected_channels::Vector{Symbol},
    selected_samples::AbstractVector{Int},
)
    n_ch = length(selected_channels)
    n_times = length(selected_samples)
    n_conds = length(epochs)

    data_arrays = Vector{Array{Float64,3}}(undef, n_conds)
    n_trials_per_condition = Vector{Int}(undef, n_conds)
    sample_idx = isempty(epochs[1].data) ? selected_samples : _as_contiguous_range(selected_samples, nrow(epochs[1].data[1]))

    for (cond_idx, epoch_data) in enumerate(epochs)
        n_trials = n_epochs(epoch_data)
        n_trials_per_condition[cond_idx] = n_trials
        condition_data = Array{Float64,3}(undef, n_ch, n_times, n_trials)

        if n_trials > 0 && n_ch > 0 && n_times > 0
            first_df = epoch_data.data[1]
            ch_indices = Int[columnindex(first_df, ch) for ch in selected_channels]
            for (trial_idx, trial_df) in enumerate(epoch_data.data)
                for (ch_idx, col_idx) in enumerate(ch_indices)
                    _copy_selected_col!(condition_data, ch_idx, trial_idx, trial_df[!, col_idx], sample_idx)
                end
            end
        end

        data_arrays[cond_idx] = condition_data
    end

    return data_arrays, n_trials_per_condition
end

function _copy_tf_trial_features!(
    dest::Array{Float64,3},
    trial_idx::Int,
    cols::Vector{Vector{Float64}},
    rows::UnitRange{Int},
    n_times::Int,
    n_freqs::Int,
)
    n_ch = length(cols)
    row_offset = first(rows) - 1
    @inbounds for t_idx = 1:n_times
        base_row = row_offset + (t_idx - 1) * n_freqs
        feature_idx = 1
        for f_idx = 1:n_freqs
            row_idx = base_row + f_idx
            for ch_idx = 1:n_ch
                dest[feature_idx, t_idx, trial_idx] = cols[ch_idx][row_idx]
                feature_idx += 1
            end
        end
    end
    return nothing
end

function _copy_tf_trial_features!(
    dest::Array{Float64,3},
    trial_idx::Int,
    cols::Vector{Vector{Float64}},
    rows::AbstractVector{Int},
    n_times::Int,
    n_freqs::Int,
)
    n_ch = length(cols)
    @inbounds for t_idx = 1:n_times
        base_k = (t_idx - 1) * n_freqs
        feature_idx = 1
        for f_idx = 1:n_freqs
            row_idx = rows[base_k+f_idx]
            for ch_idx = 1:n_ch
                dest[feature_idx, t_idx, trial_idx] = cols[ch_idx][row_idx]
                feature_idx += 1
            end
        end
    end
    return nothing
end

function _extract_tf_epoch_arrays(
    epochs::AbstractVector{<:TimeFreqEpochData},
    selected_channels::Vector{Symbol},
    selected_rows::AbstractVector{Int},
)
    n_ch = length(selected_channels)
    first_df = epochs[1].data_power[1]
    unique_times = sort!(unique(_extract_selected_f64(first_df[!, :time], selected_rows)))
    unique_freqs = sort!(unique(_extract_selected_f64(first_df[!, :freq], selected_rows)))
    n_times = length(unique_times)
    n_freqs = length(unique_freqs)
    n_features = n_ch * n_freqs

    row_idx_spec = _as_contiguous_range(selected_rows, nrow(first_df))
    n_conds = length(epochs)
    data_arrays = Vector{Array{Float64,3}}(undef, n_conds)
    n_trials_per_condition = Vector{Int}(undef, n_conds)
    cols = Vector{Vector{Float64}}(undef, n_ch)

    for (cond_idx, epoch_data) in enumerate(epochs)
        n_trials = length(epoch_data.data_power)
        n_trials_per_condition[cond_idx] = n_trials
        condition_data = Array{Float64,3}(undef, n_features, n_times, n_trials)

        if n_trials > 0 && n_features > 0 && n_times > 0
            first_df_cond = epoch_data.data_power[1]
            ch_indices = Int[columnindex(first_df_cond, ch) for ch in selected_channels]
            for (trial_idx, trial_df) in enumerate(epoch_data.data_power)
                for (ch_idx, col_idx) in enumerate(ch_indices)
                    col = trial_df[!, col_idx]
                    cols[ch_idx] = col isa Vector{Float64} ? col : Vector{Float64}(col)
                end
                _copy_tf_trial_features!(condition_data, trial_idx, cols, row_idx_spec, n_times, n_freqs)
            end
        end

        data_arrays[cond_idx] = condition_data
    end

    return data_arrays, n_trials_per_condition, unique_times, n_freqs
end

function _prepare_decoding_data(epochs::AbstractVector{<:EpochData})
    rows = isempty(epochs[1].data) ? (1:0) : (1:nrow(epochs[1].data[1]))
    return _extract_epoch_arrays(epochs, channel_labels(epochs), rows)
end

function _prepare_decoding_data(epochs::AbstractVector{<:TimeFreqEpochData})
    rows = isempty(epochs[1].data_power) ? (1:0) : (1:nrow(epochs[1].data_power[1]))
    data_arrays, n_trials_per_condition, _, _ = _extract_tf_epoch_arrays(epochs, channel_labels(epochs), rows)
    return data_arrays, n_trials_per_condition
end

function _equalize_trials(
    data_arrays::Vector{Array{Float64,3}},
    n_trials_per_condition::Vector{Int},
    rng::AbstractRNG = Random.default_rng(),
)
    min_trials = minimum(n_trials_per_condition)
    for (cond_idx, data_array) in enumerate(data_arrays)
        if size(data_array, 3) > min_trials
            selected_trials = sort!(Random.randperm(rng, size(data_array, 3))[1:min_trials])
            data_arrays[cond_idx] = data_array[:, :, selected_trials]
        end
    end
    return data_arrays, fill(min_trials, length(data_arrays))
end

function _extract_timepoint_data!(
    x_all_t::Matrix{Float64},
    labels::Vector{Int},
    data_arrays::Vector{Array{Float64,3}},
    shuffled_indices::Vector{Vector{Int}},
    t::Int,
)
    row = 1
    n_features = size(x_all_t, 1)
    @inbounds for (cond_idx, cond_data) in enumerate(data_arrays)
        for trial_idx in shuffled_indices[cond_idx]
            @simd for f = 1:n_features
                x_all_t[f, row] = cond_data[f, t, trial_idx]
            end
            labels[row] = cond_idx
            row += 1
        end
    end
    return nothing
end

function _precompute_cv_splits(n_trials_per_condition::Vector{Int}, n_folds::Int)
    splits = Vector{Tuple{Vector{Int},Vector{Int}}}(undef, n_folds)
    # Integer division drops remainder trials per condition so fold sizes remain balanced across classes
    n_trials_per_fold = Int[div(n_trials, n_folds) for n_trials in n_trials_per_condition]
    n_test_total = sum(n_trials_per_fold)
    n_train_total = sum(n_trials_per_condition) - n_test_total

    for fold = 1:n_folds
        test_indices = sizehint!(Int[], n_test_total)
        train_indices = sizehint!(Int[], n_train_total)
        trial_start = 1

        for (cond_idx, n_trials) in enumerate(n_trials_per_condition)
            n_per_fold = n_trials_per_fold[cond_idx]
            if n_per_fold > 0
                fold_start = trial_start + (fold - 1) * n_per_fold
                fold_end = trial_start + fold * n_per_fold - 1
                cond_end = trial_start + n_trials - 1
                append!(test_indices, fold_start:fold_end)
                fold_start > trial_start && append!(train_indices, trial_start:(fold_start-1))
                fold_end < cond_end && append!(train_indices, (fold_end+1):cond_end)
            end
            trial_start += n_trials
        end

        splits[fold] = (train_indices, test_indices)
    end

    return splits
end

function _decode_iteration!(
    iter::Int,
    seed::UInt64,
    data_arrays::Vector{Array{Float64,3}},
    eff_trials_per_condition::Vector{Int},
    equalize_trials::Bool,
    cv_splits::Vector{Tuple{Vector{Int},Vector{Int}}},
    cost::Float64,
    all_accuracies::Array{Float64,3},
    all_predictions::Union{Nothing,Array{Float64,4}},
    tc::Array{Float64,3},
    should_stop::Threads.Atomic{Bool},
    completed::Union{Nothing,Threads.Atomic{Int}},
)
    should_stop[] && return nothing

    n_classes = length(data_arrays)
    n_features = size(data_arrays[1], 1)
    n_timepoints = size(all_accuracies, 3)
    n_folds = length(cv_splits)
    total_trials = sum(eff_trials_per_condition)
    local_max_train = maximum(length(first(s)) for s in cv_splits)
    local_max_test = maximum(length(last(s)) for s in cv_splits)

    local_rng = Random.Xoshiro(seed)
    X_all_t = Matrix{Float64}(undef, n_features, total_trials)
    labels = Vector{Int}(undef, total_trials)
    X_train_buf = Matrix{Float64}(undef, n_features, local_max_train)
    X_test_buf = Matrix{Float64}(undef, n_features, local_max_test)
    y_train_buf = Vector{Int}(undef, local_max_train)
    y_test_buf = Vector{Int}(undef, local_max_test)
    class_correct = Vector{Int}(undef, n_classes)
    class_total = Vector{Int}(undef, n_classes)
    ws = LiblinearWorkspace(n_features, local_max_train, cost, n_classes)

    shuffled_indices = Vector{Vector{Int}}(undef, n_classes)
    for c = 1:n_classes
        perm = Random.randperm(local_rng, size(data_arrays[c], 3))
        equalize_trials && resize!(perm, eff_trials_per_condition[c])
        shuffled_indices[c] = perm
    end

    save_preds = !isnothing(all_predictions)

    for t = 1:n_timepoints
        should_stop[] && break
        _extract_timepoint_data!(X_all_t, labels, data_arrays, shuffled_indices, t)

        for fold = 1:n_folds
            train_indices, test_indices = cv_splits[fold]
            n_train = length(train_indices)
            n_test = length(test_indices)

            @inbounds for (i, idx) in enumerate(train_indices)
                @simd for f = 1:n_features
                    X_train_buf[f, i] = X_all_t[f, idx]
                end
                y_train_buf[i] = labels[idx]
            end

            @inbounds for (i, idx) in enumerate(test_indices)
                @simd for f = 1:n_features
                    X_test_buf[f, i] = X_all_t[f, idx]
                end
                y_test_buf[i] = labels[idx]
            end

            # Within-fold z-score standardization strictly from training-fold statistics
            @inbounds for f = 1:n_features
                sum_f = 0.0
                @simd for i = 1:n_train
                    sum_f += X_train_buf[f, i]
                end
                μ_f = sum_f / n_train
                ss_f = 0.0
                @simd for i = 1:n_train
                    d = X_train_buf[f, i] - μ_f
                    ss_f += d * d
                end
                σ_f = n_train > 1 ? sqrt(max(ss_f / (n_train - 1), 0.0)) : 0.0
                inv_σ = (σ_f == 0.0 || !isfinite(σ_f)) ? 1.0 : (1.0 / σ_f)
                @simd for i = 1:n_train
                    X_train_buf[f, i] = (X_train_buf[f, i] - μ_f) * inv_σ
                end
                @simd for i = 1:n_test
                    X_test_buf[f, i] = (X_test_buf[f, i] - μ_f) * inv_σ
                end
            end

            w, labels_buf, nr_c = _train_prealloc_zero!(ws, @view(X_train_buf[:, 1:n_train]), @view(y_train_buf[1:n_train]))
            n_correct = 0
            if !equalize_trials
                fill!(class_correct, 0)
                fill!(class_total, 0)
            end

            l1 = labels_buf[1]
            l2 = nr_c >= 2 ? labels_buf[2] : l1
            b_bin = nr_c == 2 ? w[n_features+1] : 0.0

            @inbounds for i = 1:n_test
                pred = if nr_c == 2
                    val = b_bin
                    @simd for f = 1:n_features
                        val += w[f] * X_test_buf[f, i]
                    end
                    val > 0.0 ? l1 : l2
                else
                    best_c = 1
                    max_val = -Inf
                    for c = 1:nr_c
                        val = w[n_features*nr_c+c]
                        for f = 1:n_features
                            val += w[(f-1)*nr_c+c] * X_test_buf[f, i]
                        end
                        if val > max_val
                            max_val = val
                            best_c = c
                        end
                    end
                    Int(labels_buf[best_c])
                end

                y_true = y_test_buf[i]
                is_hit = y_true == pred
                n_correct += is_hit
                if !equalize_trials
                    class_total[y_true] += 1
                    class_correct[y_true] += is_hit
                end
                tc[t, y_true, pred] += 1.0
                if save_preds
                    (all_predictions::Array{Float64,4})[iter, fold, t, i] = Float64(pred)
                end
            end

            if equalize_trials
                all_accuracies[iter, fold, t] = n_correct / n_test
            else
                recall_sum = 0.0
                n_present = 0
                @inbounds for c = 1:n_classes
                    if class_total[c] > 0
                        recall_sum += class_correct[c] / class_total[c]
                        n_present += 1
                    end
                end
                all_accuracies[iter, fold, t] = n_present > 0 ? (recall_sum / n_present) : 0.0
            end
        end

        if !isnothing(completed) && !should_stop[]
            Threads.atomic_add!(completed, 1)
        end
        yield()
    end

    return nothing
end

function _decode_core(
    data_arrays::Vector{Array{Float64,3}},
    n_trials_per_condition::Vector{Int},
    time_points::Vector{Float64},
    selected_channels::Vector{Symbol},
    condition_names::Vector{String},
    file::String;
    n_iterations::Int,
    n_folds::Int,
    equalize_trials::Bool,
    cost::Real,
    save_raw_predictions::Bool = true,
    show_progress::Bool,
    progress_desc::String,
    rng::AbstractRNG = Random.default_rng(),
)
    n_iterations < 1 && @minimal_error("n_iterations must be >= 1, got $n_iterations")
    n_folds < 2 && @minimal_error("Need at least 2 folds for cross-validation, got $n_folds")
    cost <= 0 && @minimal_error("SVM cost parameter must be positive, got $cost")
    any(arr -> any(!isfinite, arr), data_arrays) &&
        @minimal_error("Input features contain NaN or Inf values; remove invalid samples (e.g., edge-masked TF regions) before decoding.")

    LIBSVM.LIBLINEAR.set_print(false)
    n_classes = length(data_arrays)
    n_timepoints = length(time_points)
    cost_f64 = Float64(cost)

    eff_trials_per_condition = equalize_trials ? fill(minimum(n_trials_per_condition), n_classes) : copy(n_trials_per_condition)
    min_trials_after_equalize = minimum(eff_trials_per_condition)
    if min_trials_after_equalize < n_folds
        @minimal_error(
            "Not enough trials for $n_folds-fold cross-validation. " *
            "Minimum trials per condition: $min_trials_after_equalize. " *
            "Either reduce n_folds or increase number of trials."
        )
    end

    all_accuracies = zeros(Float64, n_iterations, n_folds, n_timepoints)
    max_test_size = sum(div(n_tr, n_folds) for n_tr in eff_trials_per_condition)
    all_predictions = save_raw_predictions ? Array{Float64,4}(undef, n_iterations, n_folds, n_timepoints, max_test_size) : nothing
    iter_confusions = [zeros(Float64, n_timepoints, n_classes, n_classes) for _ = 1:n_iterations]
    cv_splits = _precompute_cv_splits(eff_trials_per_condition, n_folds)

    total_steps = n_iterations * n_timepoints
    completed = show_progress ? Threads.Atomic{Int}(0) : nothing
    if show_progress
        progress = Progress(total_steps, desc = progress_desc, showspeed = true)
        EegFun.track_progress_async(progress, completed::Threads.Atomic{Int}, total_steps)
    end

    should_stop = Threads.Atomic{Bool}(false)
    seeds = rand(rng, UInt64, n_iterations)

    try
        Threads.@threads for iter = 1:n_iterations
            _decode_iteration!(
                iter,
                seeds[iter],
                data_arrays,
                eff_trials_per_condition,
                equalize_trials,
                cv_splits,
                cost_f64,
                all_accuracies,
                all_predictions,
                iter_confusions[iter],
                should_stop,
                completed,
            )
        end
    catch e
        Threads.atomic_xchg!(should_stop, true)
        e isa InterruptException && @minimal_error("Decoding aborted by user (Interrupt).")
        rethrow(e)
    finally
        Threads.atomic_xchg!(should_stop, true)
    end

    # Average across folds and iterations; iteration SD is a Monte Carlo dispersion, not divided by sqrt(n_iterations)
    iter_scores = mean(all_accuracies, dims = 2)
    average_score = vec(mean(iter_scores, dims = 1))
    stderror = n_iterations > 1 ? vec(std(iter_scores, dims = 1)) : vec(std(all_accuracies, dims = 2) ./ sqrt(n_folds))

    confusion_matrices = zeros(Float64, n_timepoints, n_classes, n_classes)
    for ic in iter_confusions
        confusion_matrices .+= ic
    end
    @inbounds for t = 1:n_timepoints, c = 1:n_classes
        row_sum = 0.0
        for p = 1:n_classes
            row_sum += confusion_matrices[t, c, p]
        end
        if row_sum > 0.0
            inv_sum = 1.0 / row_sum
            for p = 1:n_classes
                confusion_matrices[t, c, p] *= inv_sum
            end
        end
    end

    parameters = DecodingParameters(1.0 / n_classes, n_iterations, n_folds, n_classes == 2 ? :binary : :one_vs_all, n_classes)
    return DecodedData(
        file,
        condition_names,
        time_points,
        average_score,
        selected_channels,
        parameters;
        stderror = stderror,
        confusion_matrix = confusion_matrices,
        raw_predictions = all_predictions,
    )
end

"""
    decode_libsvm(epochs::AbstractVector{EpochData}; channel_selection = channels(),
                  interval_selection = times(), n_iterations::Int = 100, n_folds::Int = 3,
                  equalize_trials::Bool = true, cost::Real = 1.0,
                  save_raw_predictions::Bool = true, show_progress::Bool = true,
                  rng::AbstractRNG = Random.default_rng())
    decode_libsvm(epochs::AbstractVector{TimeFreqEpochData}; kwargs...)
    decode_libsvm(participant_epochs::AbstractVector{<:AbstractVector{<:MultiDataFrameEeg}}; kwargs...)

Perform time-point-by-time-point cross-validated linear SVM decoding (`LIBLINEAR`) across
conditions on `EpochData` (channel features) or `TimeFreqEpochData` (channel × frequency features).

# Keyword Arguments
- `channel_selection`: Channel selection predicate (default: `channels()`)
- `interval_selection::Interval`: Time interval filter (default: `times()`)
- `n_iterations::Int`: Random fold/subsample iterations (default: `100`)
- `n_folds::Int`: Cross-validation folds (default: `3`)
- `equalize_trials::Bool`: Subsample conditions to minimum trial count per iteration (default: `true`)
- `cost::Real`: SVM regularization parameter `C > 0` (default: `1.0`)
- `save_raw_predictions::Bool`: Store 4D prediction array in `DecodedData` (default: `true`)
- `show_progress::Bool`: Show progress bar (default: `true`)
- `rng::AbstractRNG`: Random number generator for reproducible splits

# Returns
- `DecodedData` (or `Vector{DecodedData}` for multi-participant batch input)
"""
function decode_libsvm(
    epochs::AbstractVector{T};
    channel_selection = channels(),
    interval_selection::Interval = times(),
    n_iterations::Int = 100,
    n_folds::Int = 3,
    equalize_trials::Bool = true,
    cost::Real = 1.0,
    save_raw_predictions::Bool = true,
    show_progress::Bool = true,
    rng::AbstractRNG = Random.default_rng(),
) where {T<:MultiDataFrameEeg}
    isempty(epochs) && @minimal_error("Cannot decode with empty epochs vector")
    length(epochs) < 2 && @minimal_error("Need at least 2 conditions for decoding, got $(length(epochs))")
    n_iterations < 1 && @minimal_error("n_iterations must be >= 1, got $n_iterations")
    n_folds < 2 && @minimal_error("Need at least 2 folds for cross-validation, got $n_folds")
    cost <= 0 && @minimal_error("SVM cost parameter must be positive, got $cost")

    first_ep = epochs[1]
    selected_channels = get_selected_channels(first_ep, channel_selection; include_meta = false, include_extra = false)
    isempty(selected_channels) && @minimal_error("Channel selection produced no channels")

    combined_sel = _combine_interval_sample(interval_selection, samples())
    selected_rows = get_selected_samples(first_ep, combined_sel)
    isempty(selected_rows) && @minimal_error("Interval selection produced no time points")

    if T <: EpochData
        time_points = _extract_selected_f64(first_ep.data[1][!, :time], selected_rows)
        isempty(time_points) && @minimal_error("Interval selection produced no time points")
        data_arrays, n_trials_per_condition = _extract_epoch_arrays(epochs, selected_channels, selected_rows)
        @info "Starting decoding analysis" file = first_ep.file n_conditions = length(epochs) n_channels = length(selected_channels) n_timepoints =
            length(time_points) time_range = (first(time_points), last(time_points)) n_iterations = n_iterations n_folds = n_folds trials_per_condition =
            n_trials_per_condition
        desc = "Decoding (iter × time): "
    else
        data_arrays, n_trials_per_condition, time_points, n_freqs = _extract_tf_epoch_arrays(epochs, selected_channels, selected_rows)
        isempty(time_points) && @minimal_error("Interval selection produced no time points")
        n_features = size(data_arrays[1], 1)
        @info "Starting TF decoding analysis" file = first_ep.file n_conditions = length(epochs) n_channels = length(selected_channels) n_frequencies =
            n_freqs n_features = n_features n_timepoints = length(time_points) time_range = (first(time_points), last(time_points)) n_iterations =
            n_iterations n_folds = n_folds trials_per_condition = n_trials_per_condition
        desc = "TF Decoding (iter × time): "
    end

    return _decode_core(
        data_arrays,
        n_trials_per_condition,
        time_points,
        selected_channels,
        String[ep.condition_name for ep in epochs],
        first_ep.file;
        n_iterations,
        n_folds,
        equalize_trials,
        cost,
        save_raw_predictions,
        show_progress,
        progress_desc = desc,
        rng = rng,
    )
end

function decode_libsvm(participant_epochs::AbstractVector{<:AbstractVector{<:MultiDataFrameEeg}}; kwargs...)
    return [decode_libsvm(epochs; kwargs...) for epochs in participant_epochs]
end

function grand_average(dat::AbstractVector{<:DecodedData})
    isempty(dat) && @minimal_error("Cannot create grand average from empty decoded data list")
    length(dat) == 1 && return copy(dat[1])

    first_decoded = dat[1]
    first_times = first_decoded.times
    first_condition_names = first_decoded.condition_names
    first_channels = first_decoded.channels
    first_params = first_decoded.parameters

    for decoded in @view dat[2:end]
        decoded.times != first_times && @minimal_error("DecodedData objects have inconsistent time vectors")
        decoded.condition_names != first_condition_names && @minimal_error("DecodedData objects have inconsistent condition names")
        decoded.channels != first_channels && @minimal_error("DecodedData objects have inconsistent channels")
        decoded.parameters.n_classes != first_params.n_classes && @minimal_error(
            "DecodedData objects have inconsistent number of classes: $(first_params.n_classes) vs $(decoded.parameters.n_classes)"
        )
        decoded.parameters.n_iterations != first_params.n_iterations &&
            @minimal_warning "DecodedData objects have different n_iterations: $(first_params.n_iterations) vs $(decoded.parameters.n_iterations)"
        decoded.parameters.n_folds != first_params.n_folds &&
            @minimal_warning "DecodedData objects have different n_folds: $(first_params.n_folds) vs $(decoded.parameters.n_folds)"
        decoded.parameters.chance_level != first_params.chance_level &&
            @minimal_warning "DecodedData objects have different chance_level: $(first_params.chance_level) vs $(decoded.parameters.chance_level)"
    end

    n_participants = length(dat)
    n_times = length(first_times)
    inv_n = 1.0 / n_participants
    inv_sqrt_n = 1.0 / sqrt(n_participants)

    grand_avg_accuracy = zeros(Float64, n_times)
    for d in dat
        grand_avg_accuracy .+= d.average_score
    end
    grand_avg_accuracy .*= inv_n

    grand_avg_stderror = zeros(Float64, n_times)
    @inbounds for t = 1:n_times
        μ = grand_avg_accuracy[t]
        ss = 0.0
        for p = 1:n_participants
            diff = dat[p].average_score[t] - μ
            ss += diff * diff
        end
        grand_avg_stderror[t] = sqrt(ss / (n_participants - 1)) * inv_sqrt_n
    end

    grand_avg_confusion = if !isnothing(first_decoded.confusion_matrix)
        acc_cm = zeros(Float64, size(first_decoded.confusion_matrix))
        for d in dat
            acc_cm .+= d.confusion_matrix::Array{Float64,3}
        end
        acc_cm .*= inv_n
        acc_cm
    else
        nothing
    end

    return DecodedData(
        "grand_average",
        copy(first_condition_names),
        copy(first_times),
        grand_avg_accuracy,
        copy(first_channels),
        first_params;
        stderror = grand_avg_stderror,
        confusion_matrix = grand_avg_confusion,
        raw_predictions = nothing,
    )
end
