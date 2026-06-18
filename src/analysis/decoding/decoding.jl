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
function _prepare_decoding_data(epochs::Vector{EpochData})

    channels = channel_labels(epochs)
    n_channels = length(channels)
    n_times = length(time_vector(epochs))

    # Prepare data arrays for each condition
    data_arrays = Vector{Array{Float64,3}}(undef, length(epochs))
    n_trials_per_condition = Vector{Int}(undef, length(epochs))

    for (cond_idx, epoch_data) in enumerate(epochs)
        n_trials = length(epoch_data.data)
        n_trials_per_condition[cond_idx] = n_trials

        # [channels × time × trials]
        condition_data = Array{Float64}(undef, n_channels, n_times, n_trials)

        # Pre-allocate references to columns to avoid dictionary lookups in the inner loop
        cols = Vector{Vector{Float64}}(undef, n_channels)
        
        for (trial_idx, trial_df) in enumerate(epoch_data.data)
            for (ch_idx, ch) in enumerate(channels)
                cols[ch_idx] = trial_df[!, ch]::Vector{Float64}
            end
            
            for ch_idx = 1:n_channels
                col = cols[ch_idx]
                @inbounds @simd for t = 1:n_times
                    condition_data[ch_idx, t, trial_idx] = col[t]
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
function _prepare_decoding_data(epochs::Vector{TimeFreqEpochData})

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

        # Pre-allocate references to columns
        cols = Vector{Vector{Float64}}(undef, n_channels)
        
        for (trial_idx, trial_df) in enumerate(epoch_data.data_power)
            for (ch_idx, ch) in enumerate(channels)
                cols[ch_idx] = trial_df[!, ch]::Vector{Float64}
            end

            # Reshape: each timepoint gets all channel-frequency pairs as features
            @inbounds for t_idx = 1:n_times
                feature_idx = 1
                for f_idx = 1:n_freqs
                    # Row in the DataFrame for this (time, freq) combination
                    row_idx = (f_idx - 1) * n_times + t_idx
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

# ==============================================================================
#   TRIAL EQUALIZATION
# ==============================================================================


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
function _equalize_trials(data_arrays::Vector{Array{Float64,3}}, n_trials_per_condition::Vector{Int})
    min_trials = minimum(n_trials_per_condition)
    for (cond_idx, data_array) in enumerate(data_arrays)
        if size(data_array, 3) > min_trials
            selected_trials = sort(shuffle(1:size(data_array, 3))[1:min_trials])
            data_arrays[cond_idx] = data_array[:, :, selected_trials]
        end
    end
    n_trials_per_condition = fill(min_trials, length(data_arrays))
    return data_arrays, n_trials_per_condition
end


"""
    _shuffle_trials(data_arrays::Vector{Array{Float64, 3}})

Shuffle trials within each condition.

# Arguments
- `data_arrays::Vector{Array{Float64, 3}}`: Data arrays for each condition [channels × time × trials]

# Returns
- `shuffled_indices::Vector{Vector{Int}}`: Shuffled trial indices per condition (avoids copying data)
"""
function _shuffle_trials(data_arrays::Vector{Array{Float64,3}})
    shuffled_indices = Vector{Vector{Int}}(undef, length(data_arrays))
    for (cond_idx, data_array) in enumerate(data_arrays)
        n_trials = size(data_array, 3)
        shuffled_indices[cond_idx] = shuffle(1:n_trials)
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


"""
    _compute_confusion_matrices(all_targets::Array{Float64, 4}, all_predictions::Array{Float64, 4}, 
                                 n_trials_per_fold::Vector{Int}, n_timepoints::Int, n_classes::Int,
                                 n_iterations::Int, n_folds::Int)

Compute confusion matrices averaged across iterations and folds for each time point.

# Arguments
- `all_targets::Array{Float64, 4}`: True labels [iterations × folds × timepoints × max_test_size]
- `all_predictions::Array{Float64, 4}`: Predicted labels [iterations × folds × timepoints × max_test_size]
- `n_trials_per_fold::Vector{Int}`: Number of test trials per fold per condition
- `n_timepoints::Int`: Number of time points
- `n_classes::Int`: Number of classes/conditions
- `n_iterations::Int`: Number of iterations
- `n_folds::Int`: Number of cross-validation folds

# Returns
- `confusion_matrices::Array{Float64, 3}`: Confusion matrices [timepoints × n_classes × n_classes], normalized to proportions
"""
function _compute_confusion_matrices(
    all_targets::Array{Float64,4},
    all_predictions::Array{Float64,4},
    n_trials_per_fold::Vector{Int},
    n_timepoints::Int,
    n_classes::Int,
    n_iterations::Int,
    n_folds::Int,
)
    confusion_matrices = zeros(Float64, n_timepoints, n_classes, n_classes)
    n_test = sum(n_trials_per_fold)

    for t = 1:n_timepoints
        # Directly accumulate counts into the confusion matrix
        for iter = 1:n_iterations
            for fold = 1:n_folds
                true_t = @view all_targets[iter, fold, t, 1:n_test]
                pred_t = @view all_predictions[iter, fold, t, 1:n_test]
                @inbounds @simd for i = 1:n_test
                    t_val = true_t[i]
                    p_val = pred_t[i]
                    if t_val != 0 && p_val != 0
                        confusion_matrices[t, Int(t_val), Int(p_val)] += 1.0
                    end
                end
            end
        end

        # Normalize to proportions
        for c = 1:n_classes
            row_sum = sum(@view confusion_matrices[t, c, :])
            if row_sum > 0
                @view(confusion_matrices[t, c, :]) ./= row_sum
            end
        end
    end

    return confusion_matrices
end



# ==============================================================================
#   CROSS-VALIDATION AND DECODING
# ==============================================================================

"""
    _create_confusion_matrix(y_true::Vector{Int}, y_pred::Vector{Int}, n_classes::Int)

Create confusion matrix from true and predicted labels.

# Returns
- `confusion::Matrix{Int}`: Confusion matrix [true_class × predicted_class]
"""
function _create_confusion_matrix(y_true::Vector{Int}, y_pred::Vector{Int}, n_classes::Int)
    length(y_true) == length(y_pred) || @minimal_error("y_true and y_pred must have same length")

    # Validate all labels are in valid range [1, n_classes]
    for (i, label) in enumerate(y_true)
        (label < 1 || label > n_classes) && @minimal_error("Invalid true label at index $i: $label (must be in range [1, $n_classes])")
    end
    for (i, label) in enumerate(y_pred)
        (label < 1 || label > n_classes) && @minimal_error("Invalid predicted label at index $i: $label (must be in range [1, $n_classes])")
    end

    confusion = zeros(Int, n_classes, n_classes)
    @inbounds for (true_label, pred_label) in zip(y_true, y_pred)
        confusion[true_label, pred_label] += 1
    end
    return confusion
end


"""
    libsvm_classifier(
        X_train::AbstractMatrix{Float64},
        y_train::Vector{Int},
        X_test::AbstractMatrix{Float64};
        cost::Float64 = 1.0,
    ) -> Vector{Int}

LIBSVM classifier for decoding.

This function uses LIBSVM.jl directly (no MLJ wrapper) for SVM classification.

# Arguments
- `X_train::AbstractMatrix{Float64}`: Training data [n_samples × n_features]
- `y_train::Vector{Int}`: Training labels [n_samples]
- `X_test::AbstractMatrix{Float64}`: Test data [n_samples × n_features]

# Keyword Arguments
- `cost::Float64`: Regularization parameter (default: 1.0). Larger values = less regularization.

# Returns
- `Vector{Int}`: Predicted class labels

# Examples
```julia
y_pred = libsvm_classifier(X_train, y_train, X_test, cost=1.0)
```
"""
function libsvm_classifier(
    X_train::AbstractMatrix{Float64},
    y_train::Vector{Int},
    X_test::AbstractMatrix{Float64};
    cost::Float64 = 1.0,
)::Vector{Int}

    n_train = size(X_train, 1)
    n_test = size(X_test, 1)

    n_train == 0 && return Int[]
    n_test == 0 && return Int[]

    # Get unique classes
    classes = unique(y_train)
    sort!(classes)
    n_classes = length(classes)

    if n_classes < 2
        error("Need at least 2 classes, got $(n_classes)")
    end

    # Fast path for binary classification with standard labels [1, 2]
    # This is the common case and we can skip all mapping overhead
    if classes == [1, 2]
        # X_train is [samples × features], transpose for LIBSVM
        # Use transpose view instead of materializing to save allocations
        X_train_t = transpose(X_train)
        X_test_t = transpose(X_test)

        # Train model - LIBSVM accepts dense matrices directly
        model = LIBSVM.svmtrain(X_train_t, y_train; svmtype = LIBSVM.SVC, kernel = LIBSVM.Kernel.Linear, cost = cost)

        # Predict and return directly (no remapping needed)
        y_pred, _ = LIBSVM.svmpredict(model, X_test_t)
        return y_pred
    end

    # General path for multi-class or non-standard labels
    # Create class mapping (LIBSVM expects labels starting from 1)
    max_class = maximum(classes)
    class_to_label = zeros(Int, max_class)
    @inbounds for (i, c) in enumerate(classes)
        class_to_label[c] = i
    end
    y_train_mapped = Vector{Int}(undef, length(y_train))
    @inbounds for i in eachindex(y_train)
        y_train_mapped[i] = class_to_label[y_train[i]]
    end

    # Transpose for LIBSVM (features in columns)
    X_train_t = transpose(X_train)
    X_test_t = transpose(X_test)

    # Train model
    model = LIBSVM.svmtrain(X_train_t, y_train_mapped; svmtype = LIBSVM.SVC, kernel = LIBSVM.Kernel.Linear, cost = cost)

    # Predict
    y_pred_mapped, _ = LIBSVM.svmpredict(model, X_test_t)

    # Map predictions back to original class labels
    y_pred = Vector{Int}(undef, length(y_pred_mapped))
    @inbounds for i in eachindex(y_pred_mapped)
        y_pred[i] = classes[y_pred_mapped[i]]
    end

    return y_pred
end


# ==============================================================================
#   CORE DECODING ENGINE (shared by ERP and TF decode_libsvm)
# ==============================================================================

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
    cost::Float64,
    show_progress::Bool,
    progress_desc::String,
)
    n_classes = length(data_arrays)
    n_timepoints = length(time_points)

    # Equalize trials if requested
    if equalize_trials
        data_arrays, n_trials_per_condition = _equalize_trials(data_arrays, n_trials_per_condition)
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
    all_predictions = Array{Float64}(undef, n_iterations, n_folds, n_timepoints, max_test_size)
    all_targets = similar(all_predictions)

    # Pre-compute all CV splits once (reused for all iterations/timepoints)
    cv_splits = _precompute_cv_splits(n_trials_per_condition, n_folds)

    total_trials = sum(n_trials_per_condition)
    n_features = size(data_arrays[1], 1)

    # Initialize progress bar
    total_steps = n_iterations * n_timepoints
    if show_progress
        progress = Progress(total_steps, desc = progress_desc, showspeed = true)
    end

    # Main decoding loop
    Threads.@threads for iter = 1:n_iterations
        # Thread-local buffers (built as [features × trials] for cache-friendly extraction)
        X_all_t = Matrix{Float64}(undef, n_features, total_trials)
        labels = Vector{Int}(undef, total_trials)
        
        shuffled_indices = _shuffle_trials(data_arrays)

        for t = 1:n_timepoints
            _extract_timepoint_data!(X_all_t, labels, data_arrays, shuffled_indices, t)

            for fold = 1:n_folds
                train_indices, test_indices = cv_splits[fold]

                # Double transpose cancels out in libsvm_classifier, giving contiguous columns!
                X_train_view = transpose(@view X_all_t[:, train_indices])
                y_train = labels[train_indices]
                X_test_view = transpose(@view X_all_t[:, test_indices])
                y_test = labels[test_indices]

                y_pred = libsvm_classifier(X_train_view, y_train, X_test_view; cost = cost)

                n_correct = 0
                @inbounds for i in eachindex(y_test)
                    n_correct += (y_test[i] == y_pred[i])
                end
                accuracy = n_correct / length(y_test)
                all_accuracies[iter, fold, t] = accuracy

                all_predictions[iter, fold, t, :] = y_pred
                all_targets[iter, fold, t, :] = y_test
            end

            if show_progress
                next!(progress)
            end
        end
    end

    # Average across iterations and folds
    average_score = vec(mean(mean(all_accuracies, dims = 2), dims = 1))
    stderror = vec(std(mean(all_accuracies, dims = 2), dims = 1) / sqrt(n_iterations))

    # Compute confusion matrices
    confusion_matrices =
        _compute_confusion_matrices(all_targets, all_predictions, n_trials_per_fold, n_timepoints, n_classes, n_iterations, n_folds)

    parameters = DecodingParameters(1.0 / n_classes, n_iterations, n_folds, n_classes == 2 ? :binary : :one_vs_one, n_classes)

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

# ===========================
#   DECODE WITH DIRECT LIBSVM 
# ===========================
"""
    decode_libsvm(epochs::Vector{EpochData}; channel_selection = channels(),
                  interval_selection = times(), n_iterations::Int = 100, n_folds::Int = 3,
                  equalize_trials::Bool = true, cost::Float64 = 1.0, show_progress::Bool = true)
                  decode_libsvm(epochs::Vector{TimeFreqEpochData}; channel_selection = channels(),
                  interval_selection = times(), n_iterations::Int = 100, n_folds::Int = 3,
                  equalize_trials::Bool = true, cost::Float64 = 1.0, show_progress::Bool = true)
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
    epochs::Vector{EpochData};
    channel_selection::Function = channels(),
    interval_selection::Interval = times(),
    n_iterations::Int = 100,
    n_folds::Int = 3,
    equalize_trials::Bool = true,
    cost::Float64 = 1.0,
    show_progress::Bool = true,
)
    # Input validations
    isempty(epochs) && @minimal_error("Cannot decode with empty epochs vector")
    length(epochs) < 2 && @minimal_error("Need at least 2 conditions for decoding, got $(length(epochs))")
    n_folds < 2 && @minimal_error("Need at least 2 folds for cross-validation, got $n_folds")

    # Subset epochs by channel and interval selection
    epochs = subset(epochs; channel_selection = channel_selection, interval_selection = interval_selection, include_extra = false)
    isempty(channel_labels(epochs[1])) && @minimal_error("Channel selection produced no channels")
    isempty(time_vector(epochs[1])) && @minimal_error("Interval selection produced no time points")

    # Prepare data from subsetted epochs
    data_arrays, n_trials_per_condition = _prepare_decoding_data(epochs)
    time_points = time_vector(epochs)

    # Log decoding parameters
    @info "Starting decoding analysis" file = epochs[1].file n_conditions = length(epochs) n_channels = length(channel_labels(epochs)) n_timepoints =
        length(time_points) time_range = (first(time_points), last(time_points)) n_iterations = n_iterations n_folds = n_folds trials_per_condition =
        n_trials_per_condition

    return _decode_core(
        data_arrays,
        n_trials_per_condition,
        time_points,
        channel_labels(epochs),
        [ep.condition_name for ep in epochs],
        epochs[1].file;
        n_iterations,
        n_folds,
        equalize_trials,
        cost,
        show_progress,
        progress_desc = "Decoding (iter × time): ",
    )
end

function decode_libsvm(participant_epochs::Vector{Vector{T}}; kwargs...) where {T<:MultiDataFrameEeg}
    return [decode_libsvm(epochs; kwargs...) for epochs in participant_epochs]
end


# ===========================
#   TF DECODING WITH LIBSVM
# ===========================
function decode_libsvm(
    epochs::Vector{TimeFreqEpochData};
    channel_selection::Function = channels(),
    interval_selection::Interval = times(),
    n_iterations::Int = 100,
    n_folds::Int = 3,
    equalize_trials::Bool = true,
    cost::Float64 = 1.0,
    show_progress::Bool = true,
)
    # Input validations
    isempty(epochs) && @minimal_error("Cannot decode with empty epochs vector")
    length(epochs) < 2 && @minimal_error("Need at least 2 conditions for decoding, got $(length(epochs))")
    n_folds < 2 && @minimal_error("Need at least 2 folds for cross-validation, got $n_folds")

    # Subset epochs by channel and interval selection
    epochs = subset(epochs; channel_selection = channel_selection, interval_selection = interval_selection)
    isempty(channel_labels(epochs[1])) && @minimal_error("Channel selection produced no channels")

    # Prepare data from subsetted epochs
    data_arrays, n_trials_per_condition = _prepare_decoding_data(epochs)

    # Extract unique time points from TF data (time is repeated per frequency)
    unique_time_points = sort(unique(epochs[1].data_power[1].time))

    # Log decoding parameters
    n_channels = length(channel_labels(epochs))
    n_features = size(data_arrays[1], 1)
    n_freqs = n_features ÷ n_channels
    @info "Starting TF decoding analysis" file = epochs[1].file n_conditions = length(epochs) n_channels = n_channels n_frequencies =
        n_freqs n_features = n_features n_timepoints = length(unique_time_points) time_range =
        (first(unique_time_points), last(unique_time_points)) n_iterations = n_iterations n_folds = n_folds trials_per_condition =
        n_trials_per_condition

    return _decode_core(
        data_arrays,
        n_trials_per_condition,
        unique_time_points,
        channel_labels(epochs),
        [ep.condition_name for ep in epochs],
        epochs[1].file;
        n_iterations,
        n_folds,
        equalize_trials,
        cost,
        show_progress,
        progress_desc = "TF Decoding (iter × time): ",
    )
end



# ====================================
#   GRAND AVERAGE FOR DECODING RESULTS
# ====================================

function grand_average(dat::Vector{DecodedData})

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
