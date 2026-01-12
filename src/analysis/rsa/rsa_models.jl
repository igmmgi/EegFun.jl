"""
Helper functions for creating model RDMs from various data types.

This module provides utilities to convert different types of "other" data
(behavioral, computational, theoretical) into Representational Dissimilarity
Matrices (RDMs) for comparison with neural RDMs.

# Overview

RSA compares neural RDMs (from EEG data) to model RDMs from various sources.
Different types of "other" data can be converted to RDMs:

## Static Models (single timepoint):
1. **Feature vectors** (e.g., word embeddings, image features)
   - Use `create_rdm_from_vectors()` or `create_rdm_from_matrix()`

2. **Behavioral data** (e.g., reaction times, similarity ratings)
   - Use `create_rdm_from_reaction_times()` or `create_rdm_from_similarity_ratings()`

3. **Categorical labels** (e.g., stimulus categories)
   - Use `create_rdm_from_categorical()`

## Temporal Models (multiple timepoints, e.g., eye tracking, EDA):
4. **Temporal data** (e.g., eye tracking, EDA, other physiological signals)
   - Use `create_temporal_rdm()` to compute RDMs at each time point
   - Temporal models are compared timepoint-by-timepoint to neural RDMs

5. **Pre-computed RDMs** (e.g., from other analyses)
   - Static: `Matrix{Float64}` - single RDM
   - Temporal: `Array{Float64, 3}` [time × condition × condition] - RDM at each time point

# Workflow

```julia
# Step 1: Compute neural RDM from EEG data (temporal)
epochs = [epoch_cond1, epoch_cond2, epoch_cond3]
neural_rsa = rsa(epochs, channels)  # RDMs computed at each time point

# Step 2: Create model RDMs from different data sources
# Static models (single timepoint)
static_models = Dict(
    "Reaction Times" => rts,                    # Vector{Float64} - static
    "Semantic Embeddings" => word_embeddings,   # Vector{Vector{Float64}} - static
)

# Temporal models (same temporal structure as EEG)
temporal_models = Dict(
    "Eye Tracking" => eye_data,  # Array{Float64, 3} [conditions × features × time]
    "EDA" => eda_data,           # Array{Float64, 3} [conditions × features × time]
)

# Step 3: Convert to RDMs
static_rdms, static_names = create_model_rdms(static_models)
temporal_rdms, temporal_names = create_temporal_model_rdms(temporal_models, neural_rsa.times)

# Step 4: Compare (handles both static and temporal automatically)
all_rdms = vcat(static_rdms, temporal_rdms)
all_names = vcat(static_names, temporal_names)
rsa_with_models = compare_models(neural_rsa, all_rdms, model_names=all_names)

# Step 5: Visualize results
plot_model_correlations(rsa_with_models)
```
"""

# ==============================================================================
#   CONVERT DATA TO RDMs
# ==============================================================================

"""
    create_rdm_from_vectors(
        vectors::Vector{Vector{Float64}};
        dissimilarity_measure::Symbol = :correlation,
    )

Create an RDM from a vector of feature vectors.

Each vector represents one condition/stimulus. The RDM is computed by comparing
all pairs of vectors using the specified dissimilarity measure.

# Arguments
- `vectors::Vector{Vector{Float64}}`: Vector of feature vectors, one per condition
- `dissimilarity_measure::Symbol`: Measure to use (:correlation, :spearman, :euclidean, :mahalanobis)

# Returns
- `rdm::Matrix{Float64}`: RDM matrix [condition × condition]

# Examples
```julia
# Create RDM from word embeddings
word_embeddings = [
    [0.1, 0.2, 0.3],  # word 1
    [0.2, 0.3, 0.4],  # word 2
    [0.5, 0.6, 0.7],  # word 3
]
rdm = create_rdm_from_vectors(word_embeddings, dissimilarity_measure=:euclidean)
```
"""
function create_rdm_from_vectors(
    vectors::Vector{Vector{Float64}};
    dissimilarity_measure::Symbol = :correlation,
)
    n_conditions = length(vectors)
    rdm = zeros(Float64, n_conditions, n_conditions)

    for i in 1:n_conditions
        for j in 1:n_conditions
            if i == j
                rdm[i, j] = 0.0
            else
                rdm[i, j] = _compute_dissimilarity(vectors[i], vectors[j], dissimilarity_measure)
                rdm[j, i] = rdm[i, j]  # Symmetric
            end
        end
    end

    return rdm
end

"""
    create_rdm_from_matrix(
        data_matrix::Matrix{Float64};
        dissimilarity_measure::Symbol = :correlation,
    )

Create an RDM from a data matrix where each row is a condition.

# Arguments
- `data_matrix::Matrix{Float64}`: Data matrix [conditions × features]
- `dissimilarity_measure::Symbol`: Measure to use (:correlation, :spearman, :euclidean, :mahalanobis)

# Returns
- `rdm::Matrix{Float64}`: RDM matrix [condition × condition]

# Examples
```julia
# Create RDM from image features (each row is an image, columns are features)
image_features = [
    0.1 0.2 0.3;  # image 1
    0.2 0.3 0.4;  # image 2
    0.5 0.6 0.7;  # image 3
]
rdm = create_rdm_from_matrix(image_features, dissimilarity_measure=:correlation)
```
"""
function create_rdm_from_matrix(
    data_matrix::Matrix{Float64};
    dissimilarity_measure::Symbol = :correlation,
)
    n_conditions = size(data_matrix, 1)
    rdm = zeros(Float64, n_conditions, n_conditions)

    for i in 1:n_conditions
        for j in 1:n_conditions
            if i == j
                rdm[i, j] = 0.0
            else
                pattern_i = vec(data_matrix[i, :])
                pattern_j = vec(data_matrix[j, :])
                rdm[i, j] = _compute_dissimilarity(pattern_i, pattern_j, dissimilarity_measure)
                rdm[j, i] = rdm[i, j]  # Symmetric
            end
        end
    end

    return rdm
end

"""
    create_rdm_from_distances(distances::Vector{Float64}, n_conditions::Int)

Create an RDM from a vector of pairwise distances.

The distances vector should contain all unique pairwise distances in the order:
(1,2), (1,3), ..., (1,n), (2,3), ..., (2,n), ..., (n-1,n)

# Arguments
- `distances::Vector{Float64}`: Vector of pairwise distances (upper triangular)
- `n_conditions::Int`: Number of conditions

# Returns
- `rdm::Matrix{Float64}`: RDM matrix [condition × condition]

# Examples
```julia
# Create RDM from reaction times (dissimilarity = difference in RT)
rts = [0.3, 0.5, 0.4]  # RTs for 3 conditions
distances = [abs(rts[1] - rts[2]), abs(rts[1] - rts[3]), abs(rts[2] - rts[3])]
rdm = create_rdm_from_distances(distances, 3)
```
"""
function create_rdm_from_distances(distances::Vector{Float64}, n_conditions::Int)
    expected_length = div(n_conditions * (n_conditions - 1), 2)
    if length(distances) != expected_length
        @minimal_error_throw(
            "Distances vector has length $(length(distances)), expected $expected_length for $n_conditions conditions"
        )
    end

    rdm = zeros(Float64, n_conditions, n_conditions)
    idx = 1

    for i in 1:n_conditions
        for j in (i+1):n_conditions
            rdm[i, j] = distances[idx]
            rdm[j, i] = distances[idx]
            idx += 1
        end
    end

    return rdm
end

"""
    create_rdm_from_similarity_ratings(
        similarity_matrix::Matrix{Float64};
        convert_to_dissimilarity::Bool = true,
    )

Create an RDM from a similarity matrix (e.g., from behavioral ratings).

# Arguments
- `similarity_matrix::Matrix{Float64}`: Similarity matrix [condition × condition] (higher = more similar)
- `convert_to_dissimilarity::Bool`: Whether to convert similarity to dissimilarity (default: true)

# Returns
- `rdm::Matrix{Float64}`: RDM matrix [condition × condition]

# Examples
```julia
# Create RDM from similarity ratings (1 = very similar, 7 = very different)
similarity = [
    1.0 2.0 5.0;
    2.0 1.0 4.0;
    5.0 4.0 1.0
]
rdm = create_rdm_from_similarity_ratings(similarity, convert_to_dissimilarity=true)
```
"""
function create_rdm_from_similarity_ratings(
    similarity_matrix::Matrix{Float64};
    convert_to_dissimilarity::Bool = true,
)
    if convert_to_dissimilarity
        # Normalize to [0, 1] and invert (similarity -> dissimilarity)
        max_val = maximum(similarity_matrix)
        min_val = minimum(similarity_matrix)
        if max_val == min_val
            # All values are the same
            return zeros(Float64, size(similarity_matrix))
        end
        normalized = (similarity_matrix .- min_val) ./ (max_val - min_val)
        rdm = 1.0 .- normalized
    else
        rdm = copy(similarity_matrix)
    end

    # Ensure diagonal is zero
    n = size(rdm, 1)
    for i in 1:n
        rdm[i, i] = 0.0
    end

    return rdm
end

"""
    create_rdm_from_reaction_times(rts::Vector{Float64})

Create an RDM from reaction times.

Dissimilarity is computed as the absolute difference in reaction times.

# Arguments
- `rts::Vector{Float64}`: Reaction times for each condition

# Returns
- `rdm::Matrix{Float64}`: RDM matrix [condition × condition]

# Examples
```julia
# Create RDM from reaction times
rts = [0.3, 0.5, 0.4]  # RTs for 3 conditions
rdm = create_rdm_from_reaction_times(rts)
```
"""
function create_rdm_from_reaction_times(rts::Vector{Float64})
    n_conditions = length(rts)
    rdm = zeros(Float64, n_conditions, n_conditions)

    for i in 1:n_conditions
        for j in 1:n_conditions
            if i == j
                rdm[i, j] = 0.0
            else
                rdm[i, j] = abs(rts[i] - rts[j])
                rdm[j, i] = rdm[i, j]
            end
        end
    end

    return rdm
end

"""
    create_rdm_from_categorical(categories::Vector{Int})

Create an RDM from categorical labels.

Dissimilarity is 0 if conditions are in the same category, 1 otherwise.

# Arguments
- `categories::Vector{Int}`: Category label for each condition

# Returns
- `rdm::Matrix{Float64}`: RDM matrix [condition × condition]

# Examples
```julia
# Create RDM from category labels
categories = [1, 1, 2, 2, 3]  # 5 conditions in 3 categories
rdm = create_rdm_from_categorical(categories)
```
"""
function create_rdm_from_categorical(categories::Vector{Int})
    n_conditions = length(categories)
    rdm = zeros(Float64, n_conditions, n_conditions)

    for i in 1:n_conditions
        for j in 1:n_conditions
            if i == j
                rdm[i, j] = 0.0
            else
                rdm[i, j] = categories[i] == categories[j] ? 0.0 : 1.0
                rdm[j, i] = rdm[i, j]
            end
        end
    end

    return rdm
end

# ==============================================================================
#   TEMPORAL ALIGNMENT AND RESAMPLING
# ==============================================================================

"""
    resample_temporal_data(
        temporal_data::Array{Float64, 3},
        original_times::Vector{Float64},
        target_times::Vector{Float64};
        method::Symbol = :linear,
    )

Resample temporal model data to match neural data timepoints.

Handles different sampling rates and time ranges by interpolating model data
to match the neural data's temporal structure.

# Arguments
- `temporal_data::Array{Float64, 3}`: Temporal data [conditions × features × time]
- `original_times::Vector{Float64}`: Original time points in seconds
- `target_times::Vector{Float64}`: Target time points (e.g., from neural data)
- `method::Symbol`: Interpolation method (:linear, :cubic, :nearest)

# Returns
- `resampled_data::Array{Float64, 3}`: Resampled data [conditions × features × target_time]
- `target_times::Vector{Float64}`: Target time points (same as input)

# Examples
```julia
# Model data at 100 Hz
model_data = randn(3, 2, 50)  # [conditions × features × time]
model_times = collect(-0.2:(1/100):0.3)  # 50 timepoints at 100 Hz

# Neural data at 500 Hz
neural_times = collect(-0.2:(1/500):0.8)  # 500 timepoints at 500 Hz

# Resample model data to match neural times
resampled_data, _ = resample_temporal_data(model_data, model_times, neural_times)
# Now resampled_data has shape [3 × 2 × 500]
```
"""
function resample_temporal_data(
    temporal_data::Array{Float64, 3},
    original_times::Vector{Float64},
    target_times::Vector{Float64};
    method::Symbol = :linear,
)
    n_conditions, n_features, n_original_times = size(temporal_data)
    n_target_times = length(target_times)

    if length(original_times) != n_original_times
        @minimal_error_throw(
            "Original times length ($(length(original_times))) doesn't match " *
            "temporal data time dimension ($n_original_times)"
        )
    end

    # Preallocate resampled data
    resampled_data = zeros(Float64, n_conditions, n_features, n_target_times)

    # Determine valid time range (intersection of original and target)
    valid_start = max(original_times[1], target_times[1])
    valid_end = min(original_times[end], target_times[end])

    # Find indices for valid range in target times
    target_valid_start_idx = findfirst(t -> t >= valid_start, target_times)
    target_valid_end_idx = findlast(t -> t <= valid_end, target_times)

    if isnothing(target_valid_start_idx) || isnothing(target_valid_end_idx)
        @minimal_warning(
            "No temporal overlap between model data ($(original_times[1]) to $(original_times[end]) s) " *
            "and target times ($(target_times[1]) to $(target_times[end]) s). " *
            "Returning zeros."
        )
        return resampled_data, target_times
    end

    # Interpolate each condition and feature
    for cond_idx in 1:n_conditions
        for feat_idx in 1:n_features
            # Extract time series for this condition/feature
            original_series = temporal_data[cond_idx, feat_idx, :]

            # Interpolate to target times
            if method == :linear
                # Linear interpolation
                resampled_series = zeros(Float64, n_target_times)
                for (t_idx, t_target) in enumerate(target_times)
                    if t_target < original_times[1]
                        # Extrapolate backward (use first value)
                        resampled_series[t_idx] = original_series[1]
                    elseif t_target > original_times[end]
                        # Extrapolate forward (use last value)
                        resampled_series[t_idx] = original_series[end]
                    else
                        # Interpolate
                        # Find surrounding indices
                        lower_idx = findlast(t -> t <= t_target, original_times)
                        upper_idx = findfirst(t -> t >= t_target, original_times)

                        if isnothing(lower_idx) || isnothing(upper_idx)
                            resampled_series[t_idx] = original_series[1]
                        elseif lower_idx == upper_idx
                            # Exact match
                            resampled_series[t_idx] = original_series[lower_idx]
                        else
                            # Linear interpolation
                            t_lower = original_times[lower_idx]
                            t_upper = original_times[upper_idx]
                            weight = (t_target - t_lower) / (t_upper - t_lower)
                            resampled_series[t_idx] = 
                                (1 - weight) * original_series[lower_idx] + 
                                weight * original_series[upper_idx]
                        end
                    end
                end
                resampled_data[cond_idx, feat_idx, :] = resampled_series
            elseif method == :nearest
                # Nearest neighbor interpolation
                for (t_idx, t_target) in enumerate(target_times)
                    if t_target < original_times[1]
                        resampled_data[cond_idx, feat_idx, t_idx] = original_series[1]
                    elseif t_target > original_times[end]
                        resampled_data[cond_idx, feat_idx, t_idx] = original_series[end]
                    else
                        # Find nearest timepoint
                        nearest_idx = argmin(abs.(original_times .- t_target))
                        resampled_data[cond_idx, feat_idx, t_idx] = original_series[nearest_idx]
                    end
                end
            elseif method == :cubic
                # Cubic interpolation (simplified - uses linear for now)
                @minimal_warning "Cubic interpolation not fully implemented, using linear"
                return resample_temporal_data(temporal_data, original_times, target_times; method=:linear)
            else
                @minimal_error_throw("Unknown interpolation method: $method. Use :linear, :nearest, or :cubic")
            end
        end
    end

    return resampled_data, target_times
end

"""
    create_temporal_rdm(
        temporal_data::Array{Float64, 3},
        times::Vector{Float64};
        dissimilarity_measure::Symbol = :correlation,
        align_to::Union{Vector{Float64}, Nothing} = nothing,
        interpolation_method::Symbol = :linear,
    )

Create temporal RDMs from temporal data, with optional temporal alignment.

If `align_to` is provided, the model data will be resampled to match those timepoints
before computing RDMs. This handles cases where model data has different sampling
rates than neural data.

# Arguments
- `temporal_data::Array{Float64, 3}`: Temporal data [conditions × features × time]
- `times::Vector{Float64}`: Time points in seconds for temporal_data
- `dissimilarity_measure::Symbol`: Measure to use (:correlation, :spearman, :euclidean, :mahalanobis)
- `align_to::Union{Vector{Float64}, Nothing}`: Target timepoints to align to (e.g., neural data times)
- `interpolation_method::Symbol`: Interpolation method for resampling (:linear, :nearest, :cubic)

# Returns
- `rdms::Array{Float64, 3}`: Temporal RDMs [time × condition × condition]

# Examples
```julia
# Model data at 100 Hz
model_data = randn(3, 2, 50)  # [conditions × features × time]
model_times = collect(-0.2:(1/100):0.3)  # 50 timepoints at 100 Hz

# Option 1: Use model's own timepoints
rdms = create_temporal_rdm(model_data, model_times)

# Option 2: Align to neural data timepoints (500 Hz)
neural_times = collect(-0.2:(1/500):0.8)  # 500 timepoints at 500 Hz
rdms_aligned = create_temporal_rdm(model_data, model_times, align_to=neural_times)
# Now rdms_aligned has shape [500 × 3 × 3] instead of [50 × 3 × 3]
```
"""
function create_temporal_rdm(
    temporal_data::Array{Float64, 3},
    times::Vector{Float64};
    dissimilarity_measure::Symbol = :correlation,
    align_to::Union{Vector{Float64}, Nothing} = nothing,
    interpolation_method::Symbol = :linear,
)
    # Resample if alignment is requested
    if !isnothing(align_to)
        @info "Resampling temporal model data from $(length(times)) to $(length(align_to)) timepoints"
        temporal_data, times = resample_temporal_data(
            temporal_data,
            times,
            align_to;
            method=interpolation_method,
        )
    end

    n_conditions, n_features, n_timepoints = size(temporal_data)

    if length(times) != n_timepoints
        @minimal_error_throw(
            "Time vector length ($(length(times))) doesn't match temporal data time dimension ($n_timepoints)"
        )
    end

    # Preallocate RDMs: [time × condition × condition]
    rdms = zeros(Float64, n_timepoints, n_conditions, n_conditions)

    # Compute RDM at each time point
    for t in 1:n_timepoints
        # Extract data at time t: [conditions × features]
        timepoint_data = temporal_data[:, :, t]

        # Compute RDM for this time point
        condition_patterns = Matrix{Float64}[]
        for cond_idx in 1:n_conditions
            # Extract feature vector for this condition at time t
            pattern = vec(timepoint_data[cond_idx, :])
            push!(condition_patterns, reshape(pattern, n_features, 1))
        end

        rdms[t, :, :] = _compute_rdm(condition_patterns, dissimilarity_measure)
    end

    return rdms
end


"""
    create_temporal_rdm_from_vectors(
        temporal_vectors::Vector{Vector{Vector{Float64}}},
        times::Vector{Float64};
        dissimilarity_measure::Symbol = :correlation,
    )

Create temporal RDMs from a vector of temporal feature vectors.

# Arguments
- `temporal_vectors::Vector{Vector{Vector{Float64}}}`: 
  - Outer vector: one per condition
  - Middle vector: one per time point
  - Inner vector: feature values at that time point
- `times::Vector{Float64}`: Time points in seconds
- `dissimilarity_measure::Symbol`: Measure to use

# Returns
- `rdms::Array{Float64, 3}`: Temporal RDMs [time × condition × condition]

# Examples
```julia
# Eye tracking: 3 conditions, each with 100 timepoints, 2 features (x, y) per timepoint
eye_tracking = [
    [[x1, y1], [x2, y2], ...],  # Condition 1: 100 timepoints
    [[x1, y1], [x2, y2], ...],  # Condition 2: 100 timepoints
    [[x1, y1], [x2, y2], ...],  # Condition 3: 100 timepoints
]
times = collect(-0.2:(1/250):0.2)
rdms = create_temporal_rdm_from_vectors(eye_tracking, times)
```
"""
function create_temporal_rdm_from_vectors(
    temporal_vectors::Vector{Vector{Vector{Float64}}},
    times::Vector{Float64};
    dissimilarity_measure::Symbol = :correlation,
)
    n_conditions = length(temporal_vectors)
    n_timepoints = length(times)

    # Validate structure
    for (cond_idx, cond_data) in enumerate(temporal_vectors)
        if length(cond_data) != n_timepoints
            @minimal_error_throw(
                "Condition $cond_idx has $(length(cond_data)) timepoints, expected $n_timepoints"
            )
        end
    end

    # Get number of features from first condition, first timepoint
    n_features = length(temporal_vectors[1][1])

    # Convert to array format: [conditions × features × time]
    temporal_data = zeros(Float64, n_conditions, n_features, n_timepoints)
    for cond_idx in 1:n_conditions
        for t in 1:n_timepoints
            temporal_data[cond_idx, :, t] = temporal_vectors[cond_idx][t]
        end
    end

    return create_temporal_rdm(temporal_data, times; dissimilarity_measure=dissimilarity_measure)
end

"""
    create_temporal_model_rdms(
        temporal_model_data::Dict{String, Any},
        times::Vector{Float64};
        dissimilarity_measure::Symbol = :correlation,
        align_to::Union{Vector{Float64}, Nothing} = nothing,
        interpolation_method::Symbol = :linear,
    )

Create temporal model RDMs from a dictionary of temporal data.

Automatically detects the format and converts to temporal RDMs. If `align_to` is provided,
model data will be resampled to match those timepoints (handles different sampling rates).

# Arguments
- `temporal_model_data::Dict{String, Any}`: Dictionary mapping model names to temporal data
  - `Tuple{Array{Float64, 3}, Vector{Float64}}`: (data, times) - data with its own timepoints
  - `Array{Float64, 3}`: Direct temporal data (uses `times` parameter)
  - `Vector{Vector{Vector{Float64}}}`: Temporal vectors (see `create_temporal_rdm_from_vectors`)
- `times::Vector{Float64}`: Time points in seconds (used if model data doesn't provide its own)
- `dissimilarity_measure::Symbol`: Measure to use
- `align_to::Union{Vector{Float64}, Nothing}`: Target timepoints to align to (e.g., neural data times)
- `interpolation_method::Symbol`: Interpolation method for resampling (:linear, :nearest, :cubic)

# Returns
- `model_rdms::Vector{Array{Float64, 3}}`: Vector of temporal RDMs [time × condition × condition]
- `model_names::Vector{String}`: Vector of model names

# Examples
```julia
# Option 1: Model data with same timepoints as neural data
temporal_models = Dict(
    "Eye Tracking" => eye_data,  # Array{Float64, 3} [conditions × features × time]
    "EDA" => eda_data,           # Array{Float64, 3} [conditions × features × time]
)
rdms, names = create_temporal_model_rdms(temporal_models, neural_rsa.times)

# Option 2: Model data with different sampling rates
# Eye tracking at 100 Hz, EDA at 50 Hz, neural data at 500 Hz
eye_times = collect(-0.2:(1/100):0.3)  # 100 Hz
eda_times = collect(-0.2:(1/50):0.3)   # 50 Hz
neural_times = collect(-0.2:(1/500):0.8)  # 500 Hz

temporal_models = Dict(
    "Eye Tracking" => (eye_data, eye_times),  # Tuple with data and its times
    "EDA" => (eda_data, eda_times),
)
rdms, names = create_temporal_model_rdms(
    temporal_models, 
    neural_times,  # Will be used if model doesn't provide times
    align_to=neural_times  # Resample all models to neural timepoints
)
```
"""
function create_temporal_model_rdms(
    temporal_model_data::Dict{String, Any},
    times::Vector{Float64};
    dissimilarity_measure::Symbol = :correlation,
    align_to::Union{Vector{Float64}, Nothing} = nothing,
    interpolation_method::Symbol = :linear,
)
    model_rdms = Array{Float64, 3}[]
    model_names = String[]

    for (name, data) in temporal_model_data
        rdm = nothing
        model_times = times  # Default to provided times

        if isa(data, Tuple) && length(data) == 2
            # Tuple: (data, times) - model has its own timepoints
            model_data, model_times = data
            if !isa(model_data, Array{Float64, 3})
                @minimal_error_throw(
                    "Model '$name': Tuple first element must be Array{Float64, 3}, got $(typeof(model_data))"
                )
            end
            if !isa(model_times, Vector{Float64})
                @minimal_error_throw(
                    "Model '$name': Tuple second element must be Vector{Float64}, got $(typeof(model_times))"
                )
            end
            # Use align_to if provided, otherwise use model's own times
            target_times = isnothing(align_to) ? model_times : align_to
            rdm = create_temporal_rdm(
                model_data,
                model_times;
                dissimilarity_measure=dissimilarity_measure,
                align_to=target_times,
                interpolation_method=interpolation_method,
            )
        elseif isa(data, Array{Float64, 3})
            # Direct temporal data [conditions × features × time]
            # Use align_to if provided
            target_times = isnothing(align_to) ? times : align_to
            rdm = create_temporal_rdm(
                data,
                times;
                dissimilarity_measure=dissimilarity_measure,
                align_to=target_times,
                interpolation_method=interpolation_method,
            )
        elseif isa(data, Vector{Vector{Vector{Float64}}})
            # Temporal vectors - convert to array first
            # Note: This format doesn't support different timepoints easily
            # Would need to be extended if needed
            @minimal_warning(
                "Model '$name': Vector{Vector{Vector{Float64}}} format doesn't support " *
                "different timepoints. Converting assuming times match."
            )
            # Convert to array format
            n_conditions = length(data)
            n_timepoints = length(data[1])
            n_features = length(data[1][1])
            model_data = zeros(Float64, n_conditions, n_features, n_timepoints)
            for cond_idx in 1:n_conditions
                for t in 1:n_timepoints
                    model_data[cond_idx, :, t] = data[cond_idx][t]
                end
            end
            target_times = isnothing(align_to) ? times : align_to
            rdm = create_temporal_rdm(
                model_data,
                times;
                dissimilarity_measure=dissimilarity_measure,
                align_to=target_times,
                interpolation_method=interpolation_method,
            )
        else
            @minimal_error_throw(
                "Unknown temporal model data type for '$name': $(typeof(data)). " *
                "Supported types: Array{Float64, 3} [conditions × features × time], " *
                "Tuple{Array{Float64, 3}, Vector{Float64}} (data, times), " *
                "Vector{Vector{Vector{Float64}}}"
            )
        end

        push!(model_rdms, rdm)
        push!(model_names, name)
    end

    return model_rdms, model_names
end

# ==============================================================================
#   CONVENIENCE FUNCTION FOR MULTIPLE MODEL TYPES
# ==============================================================================

"""
    create_model_rdms(
        model_data::Dict{String, Any};
        condition_order::Union{Vector{String}, Nothing} = nothing,
    )

Create multiple model RDMs from a dictionary of different model types.

This function automatically detects the type of each model and converts it to an RDM.

# Arguments
- `model_data::Dict{String, Any}`: Dictionary mapping model names to their data
  - `Vector{Vector{Float64}}`: Feature vectors → RDM via `create_rdm_from_vectors`
  - `Matrix{Float64}`: Data matrix → RDM via `create_rdm_from_matrix`
  - `Vector{Float64}`: If length matches n_conditions, treated as reaction times
  - `Matrix{Float64}` with symmetric structure: Treated as similarity matrix
- `condition_order::Union{Vector{String}, Nothing}`: Order of conditions (for validation)

# Returns
- `model_rdms::Vector{Matrix{Float64}}`: Vector of model RDMs
- `model_names::Vector{String}`: Vector of model names

# Examples
```julia
# Create RDMs from different model types
models = Dict(
    "Semantic Model" => [
        [0.1, 0.2, 0.3],  # word embeddings
        [0.2, 0.3, 0.4],
        [0.5, 0.6, 0.7],
    ],
    "Reaction Times" => [0.3, 0.5, 0.4],
    "Similarity Ratings" => [
        1.0 2.0 5.0;
        2.0 1.0 4.0;
        5.0 4.0 1.0
    ],
)
rdms, names = create_model_rdms(models)
rsa_with_models = compare_models(rsa_result, rdms, model_names=names)
```
"""
function create_model_rdms(
    model_data::Dict{String, Any};
    condition_order::Union{Vector{String}, Nothing} = nothing,
)
    model_rdms = Matrix{Float64}[]
    model_names = String[]

    for (name, data) in model_data
        rdm = nothing

        if isa(data, Vector{Vector{Float64}})
            # Feature vectors
            rdm = create_rdm_from_vectors(data)
        elseif isa(data, Matrix{Float64})
            if size(data, 1) == size(data, 2) && issymmetric(data)
                # Already an RDM or similarity matrix
                if all(diag(data) .== 0.0) && all(data .>= 0.0)
                    # Looks like an RDM already
                    rdm = copy(data)
                else
                    # Treat as similarity matrix
                    rdm = create_rdm_from_similarity_ratings(data)
                end
            else
                # Data matrix (conditions × features)
                rdm = create_rdm_from_matrix(data)
            end
        elseif isa(data, Vector{Float64})
            if !isnothing(condition_order) && length(data) == length(condition_order)
                # Reaction times
                rdm = create_rdm_from_reaction_times(data)
            else
                # Try as feature vectors (single vector per condition)
                @minimal_warning "Treating Vector{Float64} as reaction times. If this is wrong, use create_rdm_from_vectors()"
                rdm = create_rdm_from_reaction_times(data)
            end
        elseif isa(data, Vector{Int})
            # Categorical labels
            rdm = create_rdm_from_categorical(data)
        else
            @minimal_error_throw(
                "Unknown model data type for '$name': $(typeof(data)). " *
                "Supported types: Vector{Vector{Float64}}, Matrix{Float64}, Vector{Float64}, Vector{Int}"
            )
        end

        push!(model_rdms, rdm)
        push!(model_names, name)
    end

    return model_rdms, model_names
end

