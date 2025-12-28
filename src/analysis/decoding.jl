"""
Multivariate Pattern Analysis (MVPA) and Decoding functions.

This module provides functions for performing multivariate pattern classification
on EEG/ERP data, including time-point-by-time-point decoding analysis.
"""

using LinearAlgebra
using Statistics
using StatsBase
using Random

# ==============================================================================
#   DATA PREPARATION
# ==============================================================================

"""
    _prepare_decoding_data(epochs::Vector{EpochData}, channels::Vector{Symbol}, time_range::Tuple{Float64, Float64})

Prepare epoch data for decoding analysis.

Extracts data from multiple EpochData conditions, selecting specified channels
and time range. Returns arrays ready for classification.

# Arguments
- `epochs::Vector{EpochData}`: Vector of EpochData, one per condition
- `channels::Vector{Symbol}`: Channel names to include
- `time_range::Tuple{Float64, Float64}`: (start_time, end_time) in seconds

# Returns
- `data::Vector{Array{Float64, 3}}`: One array per condition [channels × time × trials]
- `times::Vector{Float64}`: Time points in seconds
- `n_trials_per_condition::Vector{Int}`: Number of trials per condition
"""
function _prepare_decoding_data(
    epochs::Vector{EpochData},
    channels::Vector{Symbol},
    time_range::Tuple{Float64, Float64},
)
    if isempty(epochs)
        @minimal_error_throw("Cannot prepare decoding data from empty epochs vector")
    end

    # Validate all epochs have same structure
    first_epoch = epochs[1]
    sample_rate = first_epoch.sample_rate
    layout = first_epoch.layout

    # Validate channels exist in layout
    all_channels = layout.data.label  # Get channel labels from layout
    missing_channels = setdiff(channels, all_channels)
    if !isempty(missing_channels)
        @minimal_error_throw("Channels not found in layout: $(join(string.(missing_channels), ", "))")
    end

    # Get time column name (could be :time or similar)
    time_col = :time  # Assuming standard time column name
    if !hasproperty(first_epoch.data[1], time_col)
        # Try to find time column
        possible_time_cols = filter(n -> occursin("time", lowercase(string(n))), names(first_epoch.data[1]))
        if isempty(possible_time_cols)
            @minimal_error_throw("No time column found in epoch data")
        end
        time_col = possible_time_cols[1]
    end

    # Get time points from first epoch
    times_all = first_epoch.data[1][!, time_col]
    start_idx = findfirst(t -> t >= time_range[1], times_all)
    end_idx = findlast(t -> t <= time_range[2], times_all)

    if isnothing(start_idx) || isnothing(end_idx) || start_idx > end_idx
        @minimal_error_throw("Invalid time range $(time_range) for data with times $(times_all[1]) to $(times_all[end])")
    end

    times = times_all[start_idx:end_idx]

    # Prepare data arrays for each condition
    data_arrays = Vector{Array{Float64, 3}}()
    n_trials_per_condition = Int[]

    for epoch_data in epochs
        n_trials = length(epoch_data.data)
        push!(n_trials_per_condition, n_trials)

        # Preallocate: [channels × time × trials]
        condition_data = zeros(Float64, length(channels), length(times), n_trials)

        for (trial_idx, trial_df) in enumerate(epoch_data.data)
            # Extract data for this trial
            trial_times = trial_df[!, time_col]
            trial_start_idx = findfirst(t -> t >= time_range[1], trial_times)
            trial_end_idx = findlast(t -> t <= time_range[2], trial_times)

            if isnothing(trial_start_idx) || isnothing(trial_end_idx)
                @minimal_warning "Trial $trial_idx has no data in time range $(time_range), skipping"
                continue
            end

            # Extract channel data
            for (ch_idx, ch_name) in enumerate(channels)
                if hasproperty(trial_df, ch_name)
                    ch_data = trial_df[trial_start_idx:trial_end_idx, ch_name]
                    # Handle case where trial time range might be slightly different
                    if length(ch_data) == length(times)
                        condition_data[ch_idx, :, trial_idx] = ch_data
                    elseif length(ch_data) > length(times)
                        # Interpolate or truncate
                        condition_data[ch_idx, :, trial_idx] = ch_data[1:length(times)]
                    else
                        # Pad with last value or NaN
                        condition_data[ch_idx, 1:length(ch_data), trial_idx] = ch_data
                        condition_data[ch_idx, (length(ch_data)+1):end, trial_idx] .= ch_data[end]
                    end
                end
            end
        end

        push!(data_arrays, condition_data)
    end

    return data_arrays, times, n_trials_per_condition
end

# ==============================================================================
#   MLJ CLASSIFIER WRAPPER
# ==============================================================================

"""
    _get_default_mlj_model()

Get default MLJ model (LinearBinaryClassifier/LDA) for decoding.

Returns a LinearBinaryClassifier from MLJLinearModels if available,
otherwise throws an error asking user to load MLJ and specify a model.
"""
function _get_default_mlj_model()
    # MLJLinearModels is loaded at module level, so use it directly
    # Note: MLJLinearModels doesn't have LinearBinaryClassifier, use LogisticClassifier instead
    return MLJLinearModels.LogisticClassifier()
end

"""
    _mlj_classifier(model, X_train::Matrix{Float64}, y_train::Vector{Int}, X_test::Matrix{Float64}, rng::AbstractRNG)

Wrapper function to use MLJ models for classification.

This function requires MLJ.jl to be loaded. Users should run `using MLJ` before calling decode() with an MLJ model.

# Arguments
- `model`: MLJ model instance
- `X_train::Matrix{Float64}`: Training data [n_samples × n_features]
- `y_train::Vector{Int}`: Training labels
- `X_test::Matrix{Float64}`: Test data [n_samples × n_features]
- `rng::AbstractRNG`: Random number generator

# Returns
- `predictions::Vector{Int}`: Predicted class labels
"""
function _mlj_classifier(
    model,
    X_train::Matrix{Float64},
    y_train::Vector{Int},
    X_test::Matrix{Float64},
    rng::AbstractRNG,
)
    # MLJ is loaded at module level, so we can use it directly
    try
        MLJ_mod = MLJ
        
        # Convert to DataFrame for MLJ (MLJ expects DataFrames)
        feature_names = [Symbol("feature_$i") for i in 1:size(X_train, 2)]
        X_train_df = DataFrame(X_train, feature_names)
        X_test_df = DataFrame(X_test, feature_names)

        # Convert labels to categorical (MLJ expects CategoricalArray for classification)
        # CategoricalArrays is already a dependency, so use it directly
        y_train_cat = CategoricalArrays.categorical(y_train, ordered = false)

        # Create and train machine
        mach = MLJ_mod.machine(model, X_train_df, y_train_cat)
        MLJ_mod.fit!(mach, verbosity = 0)

        # Predict
        y_pred_cat = MLJ_mod.predict(mach, X_test_df)

        # Convert predictions back to integers
        # Handle different prediction types (probabilistic vs deterministic)
        y_pred = if isa(y_pred_cat[1], MLJ.UnivariateFinite)
            # For probabilistic predictions, get the mode (most likely class)
            # mode() returns CategoricalValue, use levelcode to get integer
            [CategoricalArrays.levelcode(MLJ.mode(pred)) for pred in y_pred_cat]
        else
            # For deterministic predictions, convert directly (handle CategoricalValue)
            [isa(pred, CategoricalArrays.CategoricalValue) ? CategoricalArrays.levelcode(pred) : Int(pred) for pred in y_pred_cat]
        end

        # Map back to original label indices (in case MLJ renumbered them)
        unique_train = unique(y_train)
        unique_pred = unique(y_pred)
        
        # If MLJ used different labels, map them back
        if length(unique_pred) == length(unique_train) && 
           !all(p -> p in unique_train, unique_pred)
            # Create mapping
            label_map = Dict(zip(unique_pred, unique_train))
            y_pred = [label_map[p] for p in y_pred]
        end

        return y_pred
    catch e
        if isa(e, UndefVarError) || isa(e, KeyError)
            @minimal_error_throw(
                "MLJ.jl is required for MLJ models. " *
                "Please run: using MLJ before calling decode() with an MLJ model."
            )
        else
            rethrow(e)
        end
    end
end


# ==============================================================================
#   CROSS-VALIDATION AND DECODING
# ==============================================================================

"""
    _compute_accuracy(y_true::Vector{Int}, y_pred::Vector{Int})

Compute classification accuracy.

# Returns
- `accuracy::Float64`: Proportion of correct predictions
"""
function _compute_accuracy(y_true::Vector{Int}, y_pred::Vector{Int})
    if length(y_true) != length(y_pred)
        @minimal_error_throw("True and predicted labels must have same length")
    end
    return sum(y_true .== y_pred) / length(y_true)
end

"""
    _create_confusion_matrix(y_true::Vector{Int}, y_pred::Vector{Int}, n_classes::Int)

Create confusion matrix from true and predicted labels.

# Returns
- `confusion::Matrix{Float64}`: Confusion matrix [true_class × predicted_class]
"""
function _create_confusion_matrix(y_true::Vector{Int}, y_pred::Vector{Int}, n_classes::Int)
    confusion = zeros(Int, n_classes, n_classes)
    for (true_label, pred_label) in zip(y_true, y_pred)
        if 1 <= true_label <= n_classes && 1 <= pred_label <= n_classes
            confusion[true_label, pred_label] += 1
        end
    end
    return confusion
end

"""
    decode(
        epochs::Vector{EpochData},
        channels::Vector{Symbol};
        time_range::Union{Tuple{Float64, Float64}, Nothing} = nothing,
        model = nothing,
        n_iterations::Int = 100,
        n_folds::Int = 3,
        equalize_trials::Bool = true,
        rng::AbstractRNG = Random.GLOBAL_RNG,
    )

Perform multivariate pattern classification (decoding) analysis for a SINGLE participant.

This function performs time-point-by-time-point decoding analysis on epoch data from one participant,
using cross-validation to estimate classification accuracy at each time point.

**IMPORTANT**: This function processes data from a SINGLE participant. For multiple participants:
1. Call `decode()` separately for each participant
2. Collect all `DecodedData` results
3. Use `grand_average()` to average across participants

This matches erplab's workflow where decoding is performed "within each subject independently".

Requires MLJ.jl to be loaded. All classification is performed using MLJ-compatible models.

# Arguments
- `epochs::Vector{EpochData}`: Vector of EpochData from a SINGLE participant, one per condition/class.
  Example: `[epoch_cond1_participant1, epoch_cond2_participant1]` for participant 1
- `channels::Vector{Symbol}`: Channel names to include in analysis
- `time_range::Union{Tuple{Float64, Float64}, Nothing}`: Time window for analysis (default: all available)
- `model`: MLJ model instance. If `nothing` (default), attempts to load and use LinearBinaryClassifier from MLJLinearModels.
  Examples:
  - `@load SVMClassifier pkg=LIBSVM` then pass the model (matches erplab's SVM)
  - `@load LinearBinaryClassifier pkg=MLJLinearModels` for LDA
  - `@load LogisticClassifier pkg=MLJLinearModels` for logistic regression
  - Any other MLJ-compatible classifier model
- `n_iterations::Int`: Number of iterations with random shuffling (default: 100, matches erplab default)
- `n_folds::Int`: Number of cross-validation folds (default: 3, matches erplab default)
- `equalize_trials::Bool`: Whether to equalize number of trials across conditions (default: true, matches erplab)
- `rng::AbstractRNG`: Random number generator for reproducibility

# Returns
- `DecodedData`: Object containing decoding results for this participant

# Examples
```julia
using MLJ
using LIBSVM

# Typical workflow: Decode for each participant separately
# Participant 1
epochs_p1 = [load_data("p1_cond1_epochs.jld2"), load_data("p1_cond2_epochs.jld2")]
decoded_p1 = decode(epochs_p1, [:Fz, :Cz, :Pz], time_range=(-0.2, 0.8))

# Participant 2
epochs_p2 = [load_data("p2_cond1_epochs.jld2"), load_data("p2_cond2_epochs.jld2")]
decoded_p2 = decode(epochs_p2, [:Fz, :Cz, :Pz], time_range=(-0.2, 0.8))

# ... repeat for all participants ...

# Average across participants (grand average)
all_decoded = [decoded_p1, decoded_p2, decoded_p3, ...]
grand_avg_decoded = grand_average(all_decoded)
plot_decoding(grand_avg_decoded)

# Using SVM (matches erplab exactly)
SVMClassifier = @load SVMClassifier pkg=LIBSVM
svm_model = SVMClassifier(kernel=LIBSVM.Kernel.Linear)  # Linear kernel = erplab default
decoded = decode(epochs_p1, channels, model=svm_model)
```
"""
function decode(
    epochs::Vector{EpochData},
    channels::Vector{Symbol};
    time_range::Union{Tuple{Float64, Float64}, Nothing} = nothing,
    model = nothing,
    n_iterations::Int = 100,
    n_folds::Int = 3,
    equalize_trials::Bool = true,
    rng::AbstractRNG = Random.GLOBAL_RNG,
)
    if isempty(epochs)
        @minimal_error_throw("Cannot decode with empty epochs vector")
    end

    if length(epochs) < 2
        @minimal_error_throw("Need at least 2 conditions for decoding, got $(length(epochs))")
    end

    if n_folds < 2
        @minimal_error_throw("Need at least 2 folds for cross-validation, got $n_folds")
    end

    # Get time range if not specified
    if isnothing(time_range)
        first_epoch = epochs[1]
        time_col = :time
        if hasproperty(first_epoch.data[1], time_col)
            times_all = first_epoch.data[1][!, time_col]
            time_range = (times_all[1], times_all[end])
        else
            @minimal_error_throw("Cannot determine time range automatically, please specify time_range")
        end
    end

    # Prepare data
    data_arrays, times, n_trials_per_condition = _prepare_decoding_data(epochs, channels, time_range)
    n_classes = length(epochs)
    n_timepoints = length(times)

    # Equalize trials if requested
    if equalize_trials
        min_trials = minimum(n_trials_per_condition)
        for (cond_idx, data_array) in enumerate(data_arrays)
            if size(data_array, 3) > min_trials
                # Randomly select min_trials
                selected_trials = sort(shuffle(rng, 1:size(data_array, 3))[1:min_trials])
                data_arrays[cond_idx] = data_arrays[cond_idx][:, :, selected_trials]
            end
        end
        n_trials_per_condition = fill(minimum(n_trials_per_condition), n_classes)
    end

    # Preallocate results
    all_accuracies = zeros(Float64, n_iterations, n_folds, n_timepoints)
    
    # Determine max test set size for preallocation
    n_trials_per_fold = [div(n_trials, n_folds) for n_trials in n_trials_per_condition]
    max_test_size = sum(n_trials_per_fold)
    all_predictions = Array{Float64}(undef, n_iterations, n_folds, n_timepoints, max_test_size)
    all_targets = similar(all_predictions)

    # Get metadata from first epoch
    first_epoch = epochs[1]
    condition_names = [ep.condition_name for ep in epochs]

    # Main decoding loop
    for iter in 1:n_iterations
        # Shuffle trials within each condition
        shuffled_data = []
        for (cond_idx, data_array) in enumerate(data_arrays)
            n_trials = size(data_array, 3)
            trial_order = shuffle(rng, 1:n_trials)
            push!(shuffled_data, data_array[:, :, trial_order])
        end

        # Create block assignments for cross-validation
        # Each condition contributes n_trials/n_folds trials to each fold
        n_trials_per_fold = [div(n_trials, n_folds) for n_trials in n_trials_per_condition]

        # For each time point
        for t in 1:n_timepoints
            # Extract data at this time point: [channels × trials] for each condition
            timepoint_data = []
            labels = Int[]

            for (cond_idx, cond_data) in enumerate(shuffled_data)
                # Extract [channels × trials] at time t
                tp_data = cond_data[:, t, :]  # [channels × trials]
                tp_data_transposed = transpose(tp_data)  # [trials × channels]
                push!(timepoint_data, tp_data_transposed)
                append!(labels, fill(cond_idx, size(tp_data_transposed, 1)))
            end

            # Combine all conditions: [all_trials × channels]
            X_all = vcat(timepoint_data...)

            # Cross-validation
            for fold in 1:n_folds
                # Create train/test split
                test_indices = Int[]
                train_indices = Int[]

                trial_start = 1
                for (cond_idx, n_trials) in enumerate(n_trials_per_condition)
                    n_per_fold = n_trials_per_fold[cond_idx]
                    fold_start = trial_start + (fold - 1) * n_per_fold
                    fold_end = trial_start + fold * n_per_fold - 1

                    # Test set for this condition
                    append!(test_indices, fold_start:fold_end)

                    # Training set: all other trials from this condition
                    train_cond_indices = vcat(
                        trial_start:(fold_start - 1),
                        (fold_end + 1):(trial_start + n_trials - 1),
                    )
                    append!(train_indices, train_cond_indices)

                    trial_start += n_trials
                end

                # Extract training and test sets
                X_train = X_all[train_indices, :]
                y_train = labels[train_indices]
                X_test = X_all[test_indices, :]
                y_test = labels[test_indices]

                # Train classifier using MLJ
                # If no model provided, try to load default LDA model
                model_to_use = if isnothing(model)
                    _get_default_mlj_model()
                else
                    model
                end
                
                y_pred = _mlj_classifier(model_to_use, X_train, y_train, X_test, rng)

                # Compute accuracy
                accuracy = _compute_accuracy(y_test, y_pred)
                all_accuracies[iter, fold, t] = accuracy

                # Store predictions and targets (for confusion matrix if needed)
                n_test = length(y_test)
                if n_test <= size(all_predictions, 4)
                    all_predictions[iter, fold, t, 1:n_test] = y_pred
                    all_targets[iter, fold, t, 1:n_test] = y_test
                    # Fill remaining with zeros (or could use NaN/0 as sentinel)
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

    # Compute confusion matrices (average across iterations and folds)
    confusion_matrices = zeros(Float64, n_timepoints, n_classes, n_classes)
    for t in 1:n_timepoints
        all_true_t = Int[]
        all_pred_t = Int[]
        for iter in 1:n_iterations
            for fold in 1:n_folds
                n_test = sum(n_trials_per_fold)
                if n_test <= size(all_targets, 4)
                    true_t = all_targets[iter, fold, t, 1:n_test]
                    pred_t = all_predictions[iter, fold, t, 1:n_test]
                    # Filter out zeros (sentinel values)
                    valid_mask = (true_t .!= 0) .& (pred_t .!= 0)
                    if any(valid_mask)
                        append!(all_true_t, true_t[valid_mask])
                        append!(all_pred_t, pred_t[valid_mask])
                    end
                end
            end
        end
        if !isempty(all_true_t)
            confusion_matrices[t, :, :] = _create_confusion_matrix(all_true_t, all_pred_t, n_classes)
            # Normalize to proportions
            row_sums = sum(confusion_matrices[t, :, :], dims = 2)
            for c in 1:n_classes
                if row_sums[c] > 0
                    confusion_matrices[t, c, :] ./= row_sums[c]
                end
            end
        end
    end

    # Create DecodedData object
    chance_level = 1.0 / n_classes

    # Determine method name for storage
    model_to_use = isnothing(model) ? _get_default_mlj_model() : model
    method_name = try
        model_name = string(typeof(model_to_use).name.name)
        Symbol(lowercase(model_name))
    catch
        :mlj_model
    end

    decoded = DecodedData(
        first_epoch.file,
        condition_names,
        times,
        average_score,
        channels,
        first_epoch.layout,
        first_epoch.sample_rate,
        method_name;
        stderror = stderror,
        chance_level = chance_level,
        n_iterations = n_iterations,
        n_folds = n_folds,
        class_coding = :one_vs_one,  # MLJ handles multi-class internally (one-vs-one for binary)
        confusion_matrix = confusion_matrices,
        raw_predictions = all_predictions,
        analysis_info = first_epoch.analysis_info,
    )

    return decoded
end

# ==============================================================================
#   GRAND AVERAGE FOR DECODING RESULTS
# ==============================================================================

"""
    grand_average(decoded_list::Vector{DecodedData})

Create grand average decoding results across multiple participants.

Averages classification accuracy and standard errors across participants,
creating a single DecodedData object representing the group-level results.

# Arguments
- `decoded_list::Vector{DecodedData}`: Vector of DecodedData objects, one per participant

# Returns
- `DecodedData`: Grand average decoding results

# Examples
```julia
# Load decoding results from multiple participants
decoded_participants = [load_decoded("p1_decoded.jld2"), load_decoded("p2_decoded.jld2")]
grand_avg = grand_average(decoded_participants)
plot_decoding(grand_avg)
```
"""
function grand_average(decoded_list::Vector{DecodedData})
    if isempty(decoded_list)
        @minimal_error_throw("Cannot create grand average from empty decoded data list")
    end

    if length(decoded_list) == 1
        return decoded_list[1]
    end

    # Validate all decoded data have same structure
    first_decoded = decoded_list[1]
    first_times = first_decoded.times
    first_condition_names = first_decoded.condition_names
    first_channels = first_decoded.channels
    first_method = first_decoded.method
    first_class_coding = first_decoded.class_coding
    first_n_classes = first_decoded.n_classes
    first_n_iterations = first_decoded.n_iterations
    first_n_folds = first_decoded.n_folds
    first_chance_level = first_decoded.chance_level

    for decoded in decoded_list[2:end]
        if decoded.times != first_times
            @minimal_error_throw("DecodedData objects have inconsistent time vectors")
        end
        if decoded.condition_names != first_condition_names
            @minimal_error_throw("DecodedData objects have inconsistent condition names")
        end
        if decoded.channels != first_channels
            @minimal_error_throw("DecodedData objects have inconsistent channels")
        end
        if decoded.method != first_method
            @minimal_warning "DecodedData objects have different methods: $(first_method) vs $(decoded.method)"
        end
    end

    # Average accuracy across participants
    all_accuracies = hcat([d.average_score for d in decoded_list]...)
    grand_avg_accuracy = vec(mean(all_accuracies, dims = 2))

    # Compute standard error across participants
    grand_avg_stderror = vec(std(all_accuracies, dims = 2) / sqrt(length(decoded_list)))

    # Average confusion matrices if available
    grand_avg_confusion = nothing
    if !isnothing(first_decoded.confusion_matrix)
        all_confusions = cat([d.confusion_matrix for d in decoded_list]..., dims = 4)
        grand_avg_confusion = mean(all_confusions, dims = 4)[:, :, :, 1]
    end

    # Create grand average DecodedData
    grand_avg = DecodedData(
        "grand_average",
        first_condition_names,
        first_times,
        grand_avg_accuracy,
        first_channels,
        first_decoded.layout,
        first_decoded.sample_rate,
        first_method;
        stderror = grand_avg_stderror,
        chance_level = first_chance_level,
        n_iterations = first_n_iterations,
        n_folds = first_n_folds,
        class_coding = first_class_coding,
        confusion_matrix = grand_avg_confusion,
        raw_predictions = nothing,  # Don't store raw predictions for grand average
        analysis_info = first_decoded.analysis_info,
    )

    return grand_avg
end

