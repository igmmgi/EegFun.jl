"""
    _prepare_decoding_data(epochs::Vector{EpochData})

Prepare epoch data for decoding analysis.

Extracts data from multiple EpochData conditions (already subsetted).
Returns arrays ready for classification.

# Arguments
- `epochs::Vector{EpochData}`: Vector of EpochData, one per condition (already subsetted)

# Returns
- `data_arrays::Vector{Array{Float64, 3}}`: One array per condition [channels × time × trials]
- `n_trials_per_condition::Vector{Int}`: Number of trials per condition
"""
function _prepare_decoding_data(epochs::Vector{EpochData})

    channels = channel_labels(epochs)
    times = time(epochs)
    n_channels = length(channels)
    n_times = length(times)

    # Prepare data arrays for each condition
    data_arrays = Vector{Array{Float64,3}}(undef, length(epochs))
    n_trials_per_condition = Vector{Int}(undef, length(epochs))

    for (cond_idx, epoch_data) in enumerate(epochs)
        n_trials = length(epoch_data.data)
        n_trials_per_condition[cond_idx] = n_trials

        # [channels × time × trials]
        condition_data = Array{Float64}(undef, n_channels, n_times, n_trials)

        for (trial_idx, trial_df) in enumerate(epoch_data.data)
            trial_matrix = Matrix(trial_df[!, channels])  # [time × channels]
            @inbounds for ch_idx = 1:n_channels
                condition_data[ch_idx, :, trial_idx] = @view trial_matrix[:, ch_idx]
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
    _equalize_trials(data_arrays::Vector{Array{Float64, 3}}, n_trials_per_condition::Vector{Int}, n_classes::Int, rng::AbstractRNG)

Equalize the number of trials across conditions by randomly downsampling to the minimum.

This ensures balanced classification by preventing bias toward conditions with more trials.

# Arguments
- `data_arrays::Vector{Array{Float64, 3}}`: Data arrays for each condition [channels × time × trials]
- `n_trials_per_condition::Vector{Int}`: Number of trials per condition
- `n_classes::Int`: Number of classes/conditions
- `rng::AbstractRNG`: Random number generator for reproducibility

# Returns
- `data_arrays::Vector{Array{Float64, 3}}`: Equalized data arrays
- `n_trials_per_condition::Vector{Int}`: Updated trial counts (all equal to minimum)
"""
function _equalize_trials(
    data_arrays::Vector{Array{Float64,3}},
    n_trials_per_condition::Vector{Int},
    n_classes::Int,
    rng::AbstractRNG,
)
    min_trials = minimum(n_trials_per_condition)
    for (cond_idx, data_array) in enumerate(data_arrays)
        if size(data_array, 3) > min_trials
            selected_trials = sort(shuffle(rng, 1:size(data_array, 3))[1:min_trials])
            data_arrays[cond_idx] = data_array[:, :, selected_trials]
        end
    end
    n_trials_per_condition = fill(min_trials, n_classes)
    return data_arrays, n_trials_per_condition
end

# ==============================================================================
#   HELPER FUNCTIONS FOR DECODING
# ==============================================================================

"""
    _shuffle_trials(data_arrays::Vector{Array{Float64, 3}}, rng::AbstractRNG)

Shuffle trials within each condition.

# Arguments
- `data_arrays::Vector{Array{Float64, 3}}`: Data arrays for each condition [channels × time × trials]
- `rng::AbstractRNG`: Random number generator for reproducibility

# Returns
- `shuffled_indices::Vector{Vector{Int}}`: Shuffled trial indices per condition (avoids copying data)
"""
function _shuffle_trials(data_arrays::Vector{Array{Float64,3}}, rng::AbstractRNG)
    shuffled_indices = Vector{Vector{Int}}(undef, length(data_arrays))
    for (cond_idx, data_array) in enumerate(data_arrays)
        n_trials = size(data_array, 3)
        shuffled_indices[cond_idx] = shuffle(rng, 1:n_trials)
    end
    return shuffled_indices
end

"""
    _extract_timepoint_data!(X_all::Matrix{Float64}, labels::Vector{Int}, 
                              data_arrays::Vector{Array{Float64, 3}}, 
                              shuffled_indices::Vector{Vector{Int}}, t::Int)

Extract data at a specific time point into pre-allocated arrays.
Fills X_all [all_trials × channels] and labels [all_trials] in-place.

# Arguments
- `X_all::Matrix{Float64}`: Pre-allocated matrix to fill [all_trials × channels]
- `labels::Vector{Int}`: Pre-allocated vector to fill with condition labels [all_trials]
- `data_arrays::Vector{Array{Float64, 3}}`: Data arrays for each condition [channels × time × trials]
- `shuffled_indices::Vector{Vector{Int}}`: Shuffled trial indices per condition
- `t::Int`: Time point index to extract

# Returns
- Nothing (modifies `X_all` and `labels` in-place)
"""
function _extract_timepoint_data!(
    X_all::Matrix{Float64},
    labels::Vector{Int},
    data_arrays::Vector{Array{Float64,3}},
    shuffled_indices::Vector{Vector{Int}},
    t::Int,
)
    row = 1
    @inbounds for (cond_idx, cond_data) in enumerate(data_arrays)
        trial_indices = shuffled_indices[cond_idx]
        for trial_idx in trial_indices
            copyto!(@view(X_all[row, :]), @view(cond_data[:, t, trial_idx]))
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
    total_valid = n_iterations * n_folds * n_test

    for t = 1:n_timepoints
        # Preallocate arrays for all valid predictions
        all_true_t = Vector{Int}(undef, total_valid)
        all_pred_t = Vector{Int}(undef, total_valid)
        valid_count = 0

        for iter = 1:n_iterations
            for fold = 1:n_folds
                if n_test <= size(all_targets, 4)
                    true_t = @view all_targets[iter, fold, t, 1:n_test]
                    pred_t = @view all_predictions[iter, fold, t, 1:n_test]
                    # Filter out zeros (sentinel values) and copy valid ones
                    for i = 1:n_test
                        if true_t[i] != 0 && pred_t[i] != 0
                            valid_count += 1
                            all_true_t[valid_count] = Int(true_t[i])
                            all_pred_t[valid_count] = Int(pred_t[i])
                        end
                    end
                end
            end
        end

        if valid_count > 0
            # Resize to actual valid count
            resize!(all_true_t, valid_count)
            resize!(all_pred_t, valid_count)

            confusion_matrices[t, :, :] = _create_confusion_matrix(all_true_t, all_pred_t, n_classes)
            # Normalize to proportions
            row_sums = sum(confusion_matrices[t, :, :], dims = 2)
            for c = 1:n_classes
                if row_sums[c] > 0
                    confusion_matrices[t, c, :] ./= row_sums[c]
                end
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
    length(y_true) == length(y_pred) || @minimal_error_throw("y_true and y_pred must have same length")
    
    # Validate all labels are in valid range [1, n_classes]
    for (i, label) in enumerate(y_true)
        (label < 1 || label > n_classes) && @minimal_error_throw("Invalid true label at index $i: $label (must be in range [1, $n_classes])")
    end
    for (i, label) in enumerate(y_pred)
        (label < 1 || label > n_classes) && @minimal_error_throw("Invalid predicted label at index $i: $label (must be in range [1, $n_classes])")
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


# ===========================
#   DECODE WITH DIRECT LIBSVM 
# ===========================
"""
    decode_libsvm(
        epochs::Vector{EpochData};
        channel_selection::Function = channels(),
        sample_selection::Function = samples(),
        n_iterations::Int = 100,
        n_folds::Int = 3,
        equalize_trials::Bool = true,
        cost::Float64 = 1.0,
        rng::AbstractRNG = Random.GLOBAL_RNG,
    )

Perform multivariate pattern classification (decoding) analysis using LIBSVM 

This function performs time-point-by-time-point decoding analysis on epoch data from one participant,
using cross-validation to estimate classification accuracy at each time point.

# Arguments
- `epochs::Vector{EpochData}`: Vector of EpochData from a SINGLE participant, one per condition/class.
  Example: `[epoch_cond1_participant1, epoch_cond2_participant1]` for participant 1

# Keyword Arguments
- `channel_selection::Function=channels()`: Channel selection predicate. See `channels()` for options.
  - Example: `channel_selection=channels([:Fz, :Cz, :Pz])` for specific channels
  - Example: `channel_selection=channels(:Cz)` for single channel
  - Default: all channels
- `sample_selection::Function=samples()`: Sample selection predicate. See `samples()` for options.
  - Example: `sample_selection=samples((-0.2, 0.8))` for time window from -0.2 to 0.8 seconds
  - Example: `sample_selection=samples()` for all time points (default)
- `n_iterations::Int`: Number of iterations with random shuffling (default: 100, matches erplab default)
- `n_folds::Int`: Number of cross-validation folds (default: 3, matches erplab default)
- `equalize_trials::Bool`: Whether to equalize number of trials across conditions (default: true, matches erplab)
- `cost::Float64`: SVM regularization parameter (default: 1.0). Larger values = less regularization.
- `rng::AbstractRNG`: Random number generator for reproducibility

# Returns
- `DecodedData`: Object containing decoding results for this participant

# Examples
```julia
decoded = decode_libsvm(epochs; channel_selection=channels(:Cz))
```
"""
function decode_libsvm(
    epochs::Vector{EpochData};
    channel_selection::Function = channels(),
    sample_selection::Function = samples(),
    n_iterations::Int = 100,
    n_folds::Int = 3,
    equalize_trials::Bool = true,
    cost::Float64 = 1.0,
    rng::AbstractRNG = Random.GLOBAL_RNG,
)

    # Input validations
    isempty(epochs) && @minimal_error_throw("Cannot decode with empty epochs vector")
    length(epochs) < 2 && @minimal_error_throw("Need at least 2 conditions for decoding, got $(length(epochs))")
    n_folds < 2 && @minimal_error_throw("Need at least 2 folds for cross-validation, got $n_folds")

    # Subset epochs by channel and sample selection
    epochs = subset(
        epochs;
        channel_selection = channel_selection,
        sample_selection = sample_selection,
        include_extra = false,
    )
    isempty(channel_labels(epochs[1])) && @minimal_error_throw("Channel selection produced no channels")
    isempty(epochs[1].data[1][!, :time]) && @minimal_error_throw("Sample selection produced no time points")

    # Prepare data from subsetted epochs
    data_arrays, n_trials_per_condition = _prepare_decoding_data(epochs)
    n_classes = length(epochs)
    times = time(epochs)
    selected_channels = channel_labels(epochs)
    n_timepoints = length(times)

    # Equalize trials if requested
    if equalize_trials
        data_arrays, n_trials_per_condition = _equalize_trials(data_arrays, n_trials_per_condition, n_classes, rng)
    end

    # Validate that we have enough trials for cross-validation
    min_trials_after_equalize = minimum(n_trials_per_condition)
    if min_trials_after_equalize < n_folds
        @minimal_error_throw(
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
    condition_names = [ep.condition_name for ep in epochs]

    # Pre-compute all CV splits once (n_folds different splits, reused for all iterations/timepoints)
    # cv_splits[i] contains (train_indices, test_indices) for fold i
    cv_splits = _precompute_cv_splits(n_trials_per_condition, n_folds)

    # Pre-allocate reusable arrays
    total_trials = sum(n_trials_per_condition)
    n_channels = size(data_arrays[1], 1)
    X_all = Matrix{Float64}(undef, total_trials, n_channels)
    labels = Vector{Int}(undef, total_trials)

    # Main decoding loop
    for iter = 1:n_iterations
        shuffled_indices = _shuffle_trials(data_arrays, rng)

        # For each time point
        for t = 1:n_timepoints
            # Fill pre-allocated arrays
            _extract_timepoint_data!(X_all, labels, data_arrays, shuffled_indices, t)

            # Cross-validation: use pre-computed splits for each fold
            for fold = 1:n_folds
                train_indices, test_indices = cv_splits[fold]

                # Extract training and test sets (use views to avoid copying)
                X_train_view = @view X_all[train_indices, :]
                y_train = labels[train_indices]
                X_test_view = @view X_all[test_indices, :]
                y_test = labels[test_indices]

                # Classify using direct LIBSVM (pass views directly - LIBSVM accepts AbstractMatrix)
                y_pred = libsvm_classifier(X_train_view, y_train, X_test_view; cost = cost)

                # Compute accuracy (more efficient than broadcasting)
                n_correct = 0
                @inbounds for i in eachindex(y_test)
                    n_correct += (y_test[i] == y_pred[i])
                end
                accuracy = n_correct / length(y_test)
                all_accuracies[iter, fold, t] = accuracy

                # Store predictions and targets
                n_test = length(y_test)
                if n_test <= size(all_predictions, 4)
                    all_predictions[iter, fold, t, 1:n_test] = y_pred
                    all_targets[iter, fold, t, 1:n_test] = y_test
                    if n_test < size(all_predictions, 4)
                        all_predictions[iter, fold, t, (n_test+1):end] .= 0
                        all_targets[iter, fold, t, (n_test+1):end] .= 0
                    end
                end
            end
        end
    end

    # Average across iterations and folds
    average_score = vec(mean(mean(all_accuracies, dims = 2), dims = 1))  # Average over iterations and folds
    stderror = vec(std(mean(all_accuracies, dims = 2), dims = 1) / sqrt(n_iterations))  # SE across iterations

    # Compute confusion matrices
    confusion_matrices = _compute_confusion_matrices(
        all_targets,
        all_predictions,
        n_trials_per_fold,
        n_timepoints,
        n_classes,
        n_iterations,
        n_folds,
    )

    parameters = DecodingParameters(
        1.0 / n_classes,  # chance_level
        n_iterations,
        n_folds,
        n_classes == 2 ? :binary : :one_vs_one,  # class_coding: Binary for 2 classes, one-vs-one for 3+ classes
        n_classes,
    )

    return DecodedData(
        epochs[1].file,
        condition_names,
        times,
        average_score,
        selected_channels,
        parameters;
        stderror = stderror,
        confusion_matrix = confusion_matrices,
        raw_predictions = all_predictions,
    )

end

# ====================================
#   GRAND AVERAGE FOR DECODING RESULTS
# ====================================

"""
    grand_average(decoded_list::Vector{DecodedData})

Create grand average decoding results across multiple participants.

Averages classification accuracy and standard errors across participants,
creating a single DecodedData object representing the group-level results.

# Arguments
- `decoded_list::Vector{DecodedData}`: Vector of DecodedData objects, one per participant

# Returns
- `DecodedData`: Grand average decoding results
"""
function grand_average(dat::Vector{DecodedData})

    isempty(dat) && @minimal_error_throw("Cannot create grand average from empty decoded data list")
    length(dat) == 1 && return dat[1]

    # Validate all decoded data have same structure
    first_decoded = dat[1]
    first_times = first_decoded.times
    first_condition_names = first_decoded.condition_names
    first_channels = first_decoded.channels
    first_params = first_decoded.parameters

    # TODO: should we check for this? warning vs. error?
    for decoded in dat[2:end]
        decoded.times != first_times && @minimal_error_throw("DecodedData objects have inconsistent time vectors")
        decoded.condition_names != first_condition_names &&
            @minimal_error_throw("DecodedData objects have inconsistent condition names")
        decoded.channels != first_channels && @minimal_error_throw("DecodedData objects have inconsistent channels")
        decoded.parameters.n_classes != first_params.n_classes && @minimal_error_throw(
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

