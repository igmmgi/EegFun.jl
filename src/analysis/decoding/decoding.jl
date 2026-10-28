"""
    _prepare_decoding_data(epochs::Vector{EpochData})

Prepare epoch data for decoding analysis.

Extracts data from multiple EpochData conditions.
Returns arrays ready for classification.

# Arguments
- `epochs::Vector{EpochData}`: Vector of EpochData, one per condition 

# Returns
- `data_arrays::Vector{Array{Float64, 3}}`: One array per condition [channels × time × trials]
- `n_trials_per_condition::Vector{Int}`: Number of trials per condition
"""
function _prepare_decoding_data(epochs::AbstractVector{<:EpochData})

    channels = channel_labels(epochs)
    n_channels = length(channels)
    n_times = length(time_vector(epochs))

    # Prepare data arrays for each condition
    data_arrays = Vector{Array{Float64,3}}(undef, length(epochs))
    n_trials_per_condition = Vector{Int}(undef, length(epochs))

    for (cond_idx, epoch_data) in enumerate(epochs)
        n_trials = n_epochs(epoch_data)
        n_trials_per_condition[cond_idx] = n_trials

        # [channels × time × trials]
        condition_data = Array{Float64}(undef, n_channels, n_times, n_trials)

        # Pre-resolve column indices outside trial loop
        first_epoch = epoch_data.data[1]
        ch_indices = [columnindex(first_epoch, ch) for ch in channels]

        for (trial_idx, trial_df) in enumerate(epoch_data.data)
            for (ch_idx, col_idx) in enumerate(ch_indices)
                col = trial_df[!, col_idx]
                if col isa Vector{Float64}
                    col_f64 = col::Vector{Float64}
                    @inbounds @simd for t = 1:n_times
                        condition_data[ch_idx, t, trial_idx] = col_f64[t]
                    end
                else
                    @inbounds @simd for t = 1:n_times
                        condition_data[ch_idx, t, trial_idx] = Float64(col[t])
                    end
                end
            end
        end

        data_arrays[cond_idx] = condition_data
    end

    return data_arrays, n_trials_per_condition
end

"""
    _prepare_decoding_data(epochs::Vector{TimeFreqEpochData})

Prepare time-frequency epoch data for decoding analysis.

Extracts data from multiple TimeFreqEpochData conditions.
At each time point, channel × frequency pairs are used as features for classification.

# Arguments
- `epochs::Vector{TimeFreqEpochData}`: Vector of TimeFreqEpochData, one per condition 

# Returns
- `data_arrays::Vector{Array{Float64, 3}}`: One array per condition [(channels×frequencies) × time × trials]
- `n_trials_per_condition::Vector{Int}`: Number of trials per condition
"""
function _prepare_decoding_data(epochs::AbstractVector{<:TimeFreqEpochData})

    channels = channel_labels(epochs)
    n_channels = length(channels)

    # Get unique time and frequency values from the first trial of the first condition
    first_df = epochs[1].data_power[1]
    unique_times = sort(unique(first_df.time))
    unique_freqs = sort(unique(first_df.freq))
    n_times = length(unique_times)
    n_freqs = length(unique_freqs)
    n_features = n_channels * n_freqs  # Flatten channels × frequencies

    # Prepare data arrays for each condition
    data_arrays = Vector{Array{Float64,3}}(undef, length(epochs))
    n_trials_per_condition = Vector{Int}(undef, length(epochs))

    for (cond_idx, epoch_data) in enumerate(epochs)
        n_trials = length(epoch_data.data_power)
        n_trials_per_condition[cond_idx] = n_trials

        # [(channels×frequencies) × time × trials]
        condition_data = Array{Float64}(undef, n_features, n_times, n_trials)

        first_df_cond = epoch_data.data_power[1]
        ch_indices = [columnindex(first_df_cond, ch) for ch in channels]

        # Pre-allocate references to columns
        cols = Vector{Vector{Float64}}(undef, n_channels)

        for (trial_idx, trial_df) in enumerate(epoch_data.data_power)
            for (ch_idx, col_idx) in enumerate(ch_indices)
                cols[ch_idx] = trial_df[!, col_idx]::Vector{Float64}
            end

            # Reshape: each timepoint gets all channel-frequency pairs as features.
            # TF layout: time outer, freq inner (tf_morlet/multitaper/stft)
            @inbounds for t_idx = 1:n_times
                base_row = (t_idx - 1) * n_freqs
                feature_idx = 1
                for f_idx = 1:n_freqs
                    row_idx = base_row + f_idx
                    for ch_idx = 1:n_channels
                        condition_data[feature_idx, t_idx, trial_idx] = cols[ch_idx][row_idx]
                        feature_idx += 1
                    end
                end
            end
        end

        data_arrays[cond_idx] = condition_data
    end

    return data_arrays, n_trials_per_condition
end

# === TRIAL EQUALIZATION ===


"""
    _equalize_trials(data_arrays::Vector{Array{Float64, 3}}, n_trials_per_condition::Vector{Int})

Equalize the number of trials across conditions by randomly downsampling to the minimum.
This ensures balanced classification by preventing bias toward conditions with more trials.

# Arguments
- `data_arrays::Vector{Array{Float64, 3}}`: Data arrays for each condition [channels × time × trials]
- `n_trials_per_condition::Vector{Int}`: Number of trials per condition

# Returns
- `data_arrays::Vector{Array{Float64, 3}}`: Equalized data arrays
- `n_trials_per_condition::Vector{Int}`: Updated trial counts (all equal to minimum)
"""
function _equalize_trials(
    data_arrays::Vector{Array{Float64,3}},
    n_trials_per_condition::Vector{Int},
    rng::AbstractRNG = Random.default_rng(),
)
    min_trials = minimum(n_trials_per_condition)
    for (cond_idx, data_array) in enumerate(data_arrays)
        if size(data_array, 3) > min_trials
            selected_trials = sort(Random.randperm(rng, size(data_array, 3))[1:min_trials])
            data_arrays[cond_idx] = data_array[:, :, selected_trials]
        end
    end
    n_trials_per_condition = fill(min_trials, length(data_arrays))
    return data_arrays, n_trials_per_condition
end


"""
    _shuffle_trials(data_arrays::Vector{Array{Float64, 3}}, rng::AbstractRNG)

Shuffle trials within each condition.

# Arguments
- `data_arrays::Vector{Array{Float64, 3}}`: Data arrays for each condition [channels × time × trials]
- `rng::AbstractRNG`: Random number generator

# Returns
- `shuffled_indices::Vector{Vector{Int}}`: Shuffled trial indices per condition (avoids copying data)
"""
function _shuffle_trials(data_arrays::Vector{Array{Float64,3}}, rng::AbstractRNG = Random.default_rng())
    shuffled_indices = Vector{Vector{Int}}(undef, length(data_arrays))
    for (cond_idx, data_array) in enumerate(data_arrays)
        n_trials = size(data_array, 3)
        shuffled_indices[cond_idx] = Random.randperm(rng, n_trials)
    end
    return shuffled_indices
end


"""
    _extract_timepoint_data!(x_all::Matrix{Float64}, labels::Vector{Int}, 
                              data_arrays::Vector{Array{Float64, 3}}, 
                              shuffled_indices::Vector{Vector{Int}}, t::Int)

Extract data at a specific time point into pre-allocated arrays.
Fills x_all [all_trials × channels] and labels [all_trials] in-place.

# Arguments
- `x_all::Matrix{Float64}`: Pre-allocated matrix to fill [all_trials × channels]
- `labels::Vector{Int}`: Pre-allocated vector to fill with condition labels [all_trials]
- `data_arrays::Vector{Array{Float64, 3}}`: Data arrays for each condition [channels × time × trials]
- `shuffled_indices::Vector{Vector{Int}}`: Shuffled trial indices per condition
- `t::Int`: Time point index to extract

# Returns
- Nothing (modifies `x_all` and `labels` in-place)
"""
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
        trial_indices = shuffled_indices[cond_idx]
        for trial_idx in trial_indices
            @simd for f = 1:n_features
                x_all_t[f, row] = cond_data[f, t, trial_idx]
            end
            labels[row] = cond_idx
            row += 1
        end
    end
end


"""
    _precompute_cv_splits(n_trials_per_condition::Vector{Int}, n_folds::Int)

Pre-compute all train/test indices for all cross-validation folds.

# Arguments
- `n_trials_per_condition::Vector{Int}`: Number of trials per condition
- `n_folds::Int`: Number of cross-validation folds

# Returns
- `splits::Vector{Tuple{Vector{Int}, Vector{Int}}}`: Vector of (train_indices, test_indices) tuples, one per fold
"""
function _precompute_cv_splits(n_trials_per_condition::Vector{Int}, n_folds::Int)
    splits = Vector{Tuple{Vector{Int},Vector{Int}}}(undef, n_folds)
    # Note: integer division intentionally drops remainder trials (e.g. 53 trials ÷ 5 folds = 10 per fold, 3 unused).
    # This is standard practice in EEG decoding to keep fold sizes balanced across conditions.
    n_trials_per_fold = [div(n_trials, n_folds) for n_trials in n_trials_per_condition]

    for fold = 1:n_folds
        test_indices = Vector{Int}()
        train_indices = Vector{Int}()
        sizehint!(test_indices, sum(n_trials_per_fold))
        sizehint!(train_indices, sum(n_trials_per_condition) - sum(n_trials_per_fold))

        trial_start = 1
        for (cond_idx, n_trials) in enumerate(n_trials_per_condition)
            n_per_fold = n_trials_per_fold[cond_idx]

            if n_per_fold > 0
                fold_start = trial_start + (fold - 1) * n_per_fold
                fold_end = trial_start + fold * n_per_fold - 1

                # Test set for this condition
                append!(test_indices, fold_start:fold_end)

                # Training set: all other trials from this condition
                if fold_start > trial_start
                    append!(train_indices, trial_start:(fold_start-1))
                end
                if fold_end < (trial_start + n_trials - 1)
                    append!(train_indices, (fold_end+1):(trial_start+n_trials-1))
                end
            end

            trial_start += n_trials
        end

        splits[fold] = (train_indices, test_indices)
    end

    return splits
end



# === CORE DECODING ENGINE (shared by ERP and TF decode_libsvm) ===

"""
    _decode_core(data_arrays, n_trials_per_condition, time_points, selected_channels,
                 condition_names, file; n_iterations, n_folds, equalize_trials, cost,
                 show_progress, progress_desc)

Shared decoding engine used by both ERP and TF `decode_libsvm` methods.

Performs trial equalization, cross-validated SVM classification at each time point,
and assembles the final `DecodedData` result.

# Arguments
- `data_arrays::Vector{Array{Float64,3}}`: Feature arrays per condition [features × time × trials]
- `n_trials_per_condition::Vector{Int}`: Trial counts per condition
- `time_points::Vector{Float64}`: Time vector for the decoded results
- `selected_channels::Vector{Symbol}`: Channel labels used
- `condition_names::Vector{String}`: Condition name labels
- `file::String`: Source filename for the result

# Keyword Arguments
- `n_iterations::Int`: Number of shuffle iterations
- `n_folds::Int`: Number of cross-validation folds
- `equalize_trials::Bool`: Whether to equalize trial counts
- `cost::Float64`: SVM regularization parameter
- `save_raw_predictions::Bool`: Whether to store raw predictions in DecodedData
- `show_progress::Bool`: Show progress bar
- `progress_desc::String`: Progress bar description

# Returns
- `DecodedData`: Decoding results
"""
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
    LIBSVM.LIBLINEAR.set_print(false)
    n_classes = length(data_arrays)
    n_timepoints = length(time_points)

    # Equalize trials if requested
    if equalize_trials
        data_arrays, n_trials_per_condition = _equalize_trials(data_arrays, n_trials_per_condition, rng)
    end

    # Validate that we have enough trials for cross-validation
    min_trials_after_equalize = minimum(n_trials_per_condition)
    if min_trials_after_equalize < n_folds
        @minimal_error(
            "Not enough trials for $n_folds-fold cross-validation. " *
            "Minimum trials per condition: $min_trials_after_equalize. " *
            "Either reduce n_folds or increase number of trials."
        )
    end

    # Preallocate results
    all_accuracies = zeros(Float64, n_iterations, n_folds, n_timepoints)

    # Determine max test set size for preallocation
    n_trials_per_fold = [div(n_trials, n_folds) for n_trials in n_trials_per_condition]
    max_test_size = sum(n_trials_per_fold)
    all_predictions = save_raw_predictions ? Array{Float64}(undef, n_iterations, n_folds, n_timepoints, max_test_size) : nothing
    iter_confusions = [zeros(Float64, n_timepoints, n_classes, n_classes) for _ = 1:n_iterations]

    # Pre-compute all CV splits once (reused for all iterations/timepoints)
    cv_splits = _precompute_cv_splits(n_trials_per_condition, n_folds)

    total_trials = sum(n_trials_per_condition)
    n_features = size(data_arrays[1], 1)

    # Initialize progress bar
    total_steps = n_iterations * n_timepoints
    local completed = nothing
    if show_progress
        progress = Progress(total_steps, desc = progress_desc, showspeed = true)
        completed = Base.Threads.Atomic{Int}(0)
        EegFun.track_progress_async(progress, completed, total_steps)
    end

    # Initialize atomic flag to safely abort all threads on Ctrl-C
    should_stop = Threads.Atomic{Bool}(false)

    # Pre-generate deterministic seeds on calling task
    seeds = rand(rng, UInt64, n_iterations)

    try
        # Main decoding loop
        Threads.@threads for iter = 1:n_iterations
            if should_stop[]
                continue
            end

            # Thread-local RNG initialized deterministically per iteration
            local_rng = Random.Xoshiro(seeds[iter])

            # Thread-local buffers (built as [features × trials] for cache-friendly extraction)
            X_all_t = Matrix{Float64}(undef, n_features, total_trials)
            labels = Vector{Int}(undef, total_trials)

            # Pre-allocate dense buffers to prevent LIBSVM from copying non-contiguous views
            local_max_train = maximum(length(first(split)) for split in cv_splits)
            local_max_test = maximum(length(last(split)) for split in cv_splits)
            X_train_buf = Matrix{Float64}(undef, n_features, local_max_train)
            X_test_buf = Matrix{Float64}(undef, n_features, local_max_test)
            y_train_buf = Vector{Int}(undef, local_max_train)
            y_test_buf = Vector{Int}(undef, local_max_test)

            ws = LiblinearWorkspace(n_features, local_max_train, cost, n_classes)
            tc = iter_confusions[iter]

            shuffled_indices = _shuffle_trials(data_arrays, local_rng)

            for t = 1:n_timepoints
                if should_stop[]
                    break
                end

                _extract_timepoint_data!(X_all_t, labels, data_arrays, shuffled_indices, t)

                for fold = 1:n_folds
                    train_indices, test_indices = cv_splits[fold]

                    n_train = length(train_indices)
                    @inbounds for (i, idx) in enumerate(train_indices)
                        for f = 1:n_features
                            X_train_buf[f, i] = X_all_t[f, idx]
                        end
                        y_train_buf[i] = labels[idx]
                    end

                    n_test = length(test_indices)
                    @inbounds for (i, idx) in enumerate(test_indices)
                        for f = 1:n_features
                            X_test_buf[f, i] = X_all_t[f, idx]
                        end
                        y_test_buf[i] = labels[idx]
                    end

                    w, labels_buf, nr_c = _train_prealloc_zero!(ws, @view(X_train_buf[:, 1:n_train]), @view(y_train_buf[1:n_train]))
                    n_correct = 0

                    if nr_c == 2
                        l1 = labels_buf[1]
                        l2 = labels_buf[2]
                        @inbounds for i = 1:n_test
                            val = 0.0
                            for f = 1:n_features
                                val += w[f] * X_test_buf[f, i]
                            end
                            pred = val > 0.0 ? l1 : l2
                            y_true = y_test_buf[i]
                            n_correct += (y_true == pred)
                            tc[t, y_true, pred] += 1.0
                            if save_raw_predictions
                                all_predictions[iter, fold, t, i] = Float64(pred)
                            end
                        end
                    else
                        @inbounds for i = 1:n_test
                            best_c = 1
                            max_val = -Inf
                            for c = 1:nr_c
                                val = 0.0
                                for f = 1:n_features
                                    val += w[(f-1)*nr_c+c] * X_test_buf[f, i]
                                end
                                if val > max_val
                                    max_val = val
                                    best_c = c
                                end
                            end
                            pred = labels_buf[best_c]
                            y_true = y_test_buf[i]
                            n_correct += (y_true == pred)
                            tc[t, y_true, pred] += 1.0
                            if save_raw_predictions
                                all_predictions[iter, fold, t, i] = Float64(pred)
                            end
                        end
                    end
                    all_accuracies[iter, fold, t] = n_correct / n_test
                end

                if show_progress && !should_stop[]
                    Base.Threads.atomic_add!(completed::Base.Threads.Atomic{Int}, 1)
                end

                # Yield to scheduler to allow Ctrl-C (InterruptException) to be processed
                yield()
            end
        end
    catch e
        # If the main task is interrupted, ensure we set the flag to stop background threads
        Threads.atomic_xchg!(should_stop, true)
        if e isa InterruptException
            @minimal_error("Decoding aborted by user (Interrupt).")
        else
            rethrow(e)
        end
    finally
        # Ensure flag is set so any lingering tasks die
        Threads.atomic_xchg!(should_stop, true)
    end

    # Average across iterations and folds
    average_score = vec(mean(mean(all_accuracies, dims = 2), dims = 1))
    stderror = vec(std(mean(all_accuracies, dims = 2), dims = 1) / sqrt(n_iterations))

    # Compute confusion matrices
    confusion_matrices = zeros(Float64, n_timepoints, n_classes, n_classes)
    for ic in iter_confusions
        confusion_matrices .+= ic
    end
    for t = 1:n_timepoints
        for c = 1:n_classes
            row_sum = sum(@view confusion_matrices[t, c, :])
            if row_sum > 0
                @view(confusion_matrices[t, c, :]) ./= row_sum
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

# === TIME POINTS & LOGGING DISPATCH HELPERS ===

_decoding_time_points(epochs::AbstractVector{<:EpochData}) = time_vector(epochs)
_decoding_time_points(epochs::AbstractVector{<:TimeFreqEpochData}) = sort(unique(epochs[1].data_power[1].time))

function _log_decoding_start(epochs::AbstractVector{<:EpochData}, data_arrays, time_points, n_iterations, n_folds, n_trials_per_condition)
    @info "Starting decoding analysis" file = epochs[1].file n_conditions = length(epochs) n_channels = n_channels(epochs[1]) n_timepoints =
        length(time_points) time_range = (first(time_points), last(time_points)) n_iterations = n_iterations n_folds = n_folds trials_per_condition =
        n_trials_per_condition
end

function _log_decoding_start(
    epochs::AbstractVector{<:TimeFreqEpochData},
    data_arrays,
    time_points,
    n_iterations,
    n_folds,
    n_trials_per_condition,
)
    n_chans = n_channels(epochs[1])
    n_features = size(data_arrays[1], 1)
    n_freqs = n_features ÷ n_chans
    @info "Starting TF decoding analysis" file = epochs[1].file n_conditions = length(epochs) n_channels = n_chans n_frequencies = n_freqs n_features =
        n_features n_timepoints = length(time_points) time_range = (first(time_points), last(time_points)) n_iterations = n_iterations n_folds =
        n_folds trials_per_condition = n_trials_per_condition
end


# === DECODE WITH DIRECT LIBSVM ===
"""
    decode_libsvm(epochs::Vector{EpochData}; channel_selection = channels(),
                  interval_selection = times(), n_iterations::Int = 100, n_folds::Int = 3,
                  equalize_trials::Bool = true, cost::Float64 = 1.0,
                  save_raw_predictions::Bool = true, show_progress::Bool = true)
    decode_libsvm(epochs::Vector{TimeFreqEpochData}; channel_selection = channels(),
                  interval_selection = times(), n_iterations::Int = 100, n_folds::Int = 3,
                  equalize_trials::Bool = true, cost::Float64 = 1.0,
                  save_raw_predictions::Bool = true, show_progress::Bool = true)
    decode_libsvm(participant_epochs::Vector{Vector{T}}; kwargs...) where {T<:MultiDataFrameEeg}

Perform multivariate pattern classification (decoding) using LIBSVM.

Time-point-by-time-point SVM classification using k-fold cross-validation. Works with
`EpochData` (channel features) and `TimeFreqEpochData` (channel × frequency features).
The `Vector{Vector{T}}` form batch-processes all participants at once.

# Keyword Arguments
- `channel_selection`: Channel selection predicate (default: all channels)
- `interval_selection`: Time interval selection (default: all time points)
- `n_iterations::Int`: Shuffle iterations (default: 100)
- `n_folds::Int`: Cross-validation folds (default: 3)
- `equalize_trials::Bool`: Equalize trial counts (default: true)
- `cost::Float64`: SVM regularization (default: 1.0)
- `save_raw_predictions::Bool`: Save full 4D prediction array (default: true)
- `show_progress::Bool`: Show progress bar (default: true)

# Returns
- `DecodedData` (or `Vector{DecodedData}` for batch)

# Example
```julia
decoded = decode_libsvm(epochs; channel_selection = channels(:Cz))
all_decoded = decode_libsvm(participant_epochs; n_iterations = 100)
```
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
    # Input validations
    isempty(epochs) && @minimal_error("Cannot decode with empty epochs vector")
    length(epochs) < 2 && @minimal_error("Need at least 2 conditions for decoding, got $(length(epochs))")
    n_folds < 2 && @minimal_error("Need at least 2 folds for cross-validation, got $n_folds")

    # Subset epochs by channel and interval selection
    sub_epochs = if T <: EpochData
        subset(epochs; channel_selection = channel_selection, interval_selection = interval_selection, include_extra = false)
    else
        subset(epochs; channel_selection = channel_selection, interval_selection = interval_selection)
    end
    n_channels(sub_epochs[1]) == 0 && @minimal_error("Channel selection produced no channels")

    time_points = _decoding_time_points(sub_epochs)
    isempty(time_points) && @minimal_error("Interval selection produced no time points")

    # Prepare data from subsetted epochs
    data_arrays, n_trials_per_condition = _prepare_decoding_data(sub_epochs)

    # Log decoding parameters
    _log_decoding_start(sub_epochs, data_arrays, time_points, n_iterations, n_folds, n_trials_per_condition)

    desc = T <: TimeFreqEpochData ? "TF Decoding (iter × time): " : "Decoding (iter × time): "

    return _decode_core(
        data_arrays,
        n_trials_per_condition,
        time_points,
        channel_labels(sub_epochs),
        [ep.condition_name for ep in sub_epochs],
        sub_epochs[1].file;
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



# === GRAND AVERAGE FOR DECODING RESULTS ===

function grand_average(dat::AbstractVector{<:DecodedData})

    isempty(dat) && @minimal_error("Cannot create grand average from empty decoded data list")
    length(dat) == 1 && return dat[1]

    # Validate all decoded data have same structure
    first_decoded = dat[1]
    first_times = first_decoded.times
    first_condition_names = first_decoded.condition_names
    first_channels = first_decoded.channels
    first_params = first_decoded.parameters

    for decoded in dat[2:end]
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

    # Average accuracy across participants
    all_accuracies = hcat([d.average_score for d in dat]...)
    grand_avg_accuracy = vec(mean(all_accuracies, dims = 2))

    # Compute standard error across participants
    grand_avg_stderror = vec(std(all_accuracies, dims = 2) / sqrt(length(dat)))

    # Average confusion matrices if available
    grand_avg_confusion = nothing
    if !isnothing(first_decoded.confusion_matrix)
        all_confusions = cat([d.confusion_matrix for d in dat]..., dims = 4)
        grand_avg_confusion = mean(all_confusions, dims = 4)[:, :, :, 1]
    end

    # Create grand average DecodedData
    grand_avg = DecodedData(
        "grand_average",
        first_condition_names,
        first_times,
        grand_avg_accuracy,
        first_channels,
        first_params;
        stderror = grand_avg_stderror,
        confusion_matrix = grand_avg_confusion,
        raw_predictions = nothing,  # Don't store raw predictions for grand average
    )

    return grand_avg
end
