"""
Helper functions for creating model RDMs from various data types.

This module provides utilities to convert different types of "other" data
(behavioral, computational, theoretical) into Representational Dissimilarity
Matrices (RDMs) for comparison with neural RDMs.

# Core Principle

**You write code to convert your data to our format, then we handle the RSA operations.**

The package defines clear input formats. You adapt your data to match these formats
using helper functions (for common cases) or custom conversion code (for your specific data).

# Input Format Specifications

## Static Models (Single Timepoint)

For data with **one value per condition** (e.g., reaction times):

**Expected Format**: `Vector{Float64}` [n_conditions]
- One value per condition
- Example: `[0.3, 0.5, 0.4]` for 3 conditions

**Helper Function**: `create_rdm_from_reaction_times(rts)`

---

For data with **multiple features per condition** (e.g., word embeddings):

**Expected Format**: `Vector{Vector{Float64}}` or `Matrix{Float64}`
- `Vector{Vector{Float64}}`: Each inner vector is one condition's features
- `Matrix{Float64}`: [n_conditions × n_features] - each row is one condition

**Helper Functions**: 
- `create_rdm_from_vectors(vectors)`
- `create_rdm_from_matrix(data_matrix)`

---

## Temporal Models (Multiple Timepoints)

For data with **many values per condition over time** (e.g., eye tracking, EDA):

**Expected Format**: `Array{Float64, 3}` [n_conditions × n_features × n_timepoints]
- Dimension 1: Conditions (e.g., Face, Object, Scene)
- Dimension 2: Features (e.g., x_position, y_position, pupil_size for eye tracking)
- Dimension 3: Time points
- **Plus**: `Vector{Float64}` [n_timepoints] - time vector in seconds

**Helper Function**: `create_temporal_rdm(data, times; align_to=neural_times)`

**Key Feature**: Automatic temporal alignment! If your model data has different
sampling rate than neural data, use `align_to` parameter to automatically resample.

---

## Pre-computed RDMs

If you already have RDMs computed:

**Static RDM**: `Matrix{Float64}` [n_conditions × n_conditions]
**Temporal RDM**: `Array{Float64, 3}` [n_timepoints × n_conditions × n_conditions]

Pass directly to `compare_models()` - no conversion needed!

# Workflow Example

```julia
# Step 1: Compute neural RDM from EEG data
neural_rsa = rsa(epochs; channel_selection=channels([:Fz, :Cz, :Pz]))

# Step 2: Convert your model data to our format

# Example A: Reaction times (1 value per condition)
rts = [0.3, 0.5, 0.4]  # Your data: 3 conditions
rt_rdm = create_rdm_from_reaction_times(rts)  # Convert to RDM

# Example B: Eye tracking (many values per condition over time)
# Your data might be in any format - you convert it:
function convert_my_eye_data(my_eye_data, condition_names)
    n_conds = length(condition_names)
    n_features = 3  # x, y, pupil
    n_times = length(my_eye_data.times)
    
    data = zeros(Float64, n_conds, n_features, n_times)
    for (i, cond) in enumerate(condition_names)
        cond_data = my_eye_data[cond]  # Your data structure
        data[i, 1, :] = cond_data.x_positions
        data[i, 2, :] = cond_data.y_positions
        data[i, 3, :] = cond_data.pupil_size
    end
    return data, my_eye_data.times
end

eye_data_array, eye_times = convert_my_eye_data(my_eye_data, condition_names)
# Package handles temporal alignment automatically!
eye_rdms = create_temporal_rdm(eye_data_array, eye_times; align_to=neural_rsa.times)

# Step 3: Compare models
rsa_with_models = compare_models(neural_rsa, [rt_rdm, eye_rdms], 
                                  model_names=["RTs", "Eye Tracking"])

# Step 4: Visualize
plot_model_correlations(rsa_with_models)
```

# Handling Different Sampling Rates

The package automatically handles temporal alignment when you provide `align_to`:

```julia
# Neural data: 500 Hz, -0.2 to 0.8s
neural_rsa = rsa(epochs)

# Eye tracking: 100 Hz, -0.2 to 0.5s (different sampling rate!)
eye_data = Array{Float64, 3}  # [conditions × features × time] at 100 Hz
eye_times = collect(-0.2:(1/100):0.5)  # 100 Hz timepoints

# Automatic resampling to match neural data
eye_rdms = create_temporal_rdm(eye_data, eye_times; align_to=neural_rsa.times)
# Result: eye_rdms now has same timepoints as neural_rsa (500 Hz, -0.2 to 0.8s)
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

Create an RDM from feature vectors (e.g., word embeddings, image features).

# Input Format

**Expected**: `Vector{Vector{Float64}}`
- Outer vector: one element per condition
- Inner vector: feature values for that condition
- All inner vectors must have the same length (same number of features)

# Arguments
- `vectors::Vector{Vector{Float64}}`: Vector of feature vectors, one per condition
- `dissimilarity_measure::Symbol`: Measure to use (:correlation, :spearman, :euclidean, :mahalanobis)

# Returns
- `rdm::Matrix{Float64}`: RDM matrix [condition × condition]

# Examples

```julia
# Word embeddings: 3 conditions, each with 100-dimensional embedding
word_embeddings = [
    [0.1, 0.2, 0.3, ..., 0.9],  # Condition 1: 100 features
    [0.2, 0.3, 0.4, ..., 1.0],  # Condition 2: 100 features
    [0.5, 0.6, 0.7, ..., 1.2],  # Condition 3: 100 features
]
rdm = create_rdm_from_vectors(word_embeddings, dissimilarity_measure=:euclidean)
```

# Custom Data Conversion

If your data is in a different format, convert it first:

```julia
# Your data might be in a DataFrame or custom structure
# Convert to Vector{Vector{Float64}} format
function convert_my_embeddings(my_data, condition_names)
    vectors = Vector{Vector{Float64}}()
    for cond_name in condition_names
        embedding = my_data[cond_name]  # Your data structure
        push!(vectors, embedding)  # Must be Vector{Float64}
    end
    return vectors
end

embeddings = convert_my_embeddings(my_data, ["Face", "Object", "Scene"])
rdm = create_rdm_from_vectors(embeddings)
```
"""
function create_rdm_from_vectors(vectors::Vector{Vector{Float64}}; dissimilarity_measure::Symbol = :correlation)
    n_conditions = length(vectors)
    rdm = zeros(Float64, n_conditions, n_conditions)

    # Only compute upper triangle for efficiency
    for i = 1:n_conditions
        rdm[i, i] = 0.0
        for j = (i+1):n_conditions
            dissim = _compute_dissimilarity(vectors[i], vectors[j], dissimilarity_measure)
            rdm[i, j] = dissim
            rdm[j, i] = dissim  # Symmetric
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

# Input Format

**Expected**: `Matrix{Float64}` [n_conditions × n_features]
- Each row: one condition
- Each column: one feature
- Example: Image features where rows are images, columns are feature dimensions

# Arguments
- `data_matrix::Matrix{Float64}`: Data matrix [conditions × features]
- `dissimilarity_measure::Symbol`: Measure to use (:correlation, :spearman, :euclidean, :mahalanobis)

# Returns
- `rdm::Matrix{Float64}`: RDM matrix [condition × condition]

# Examples

```julia
# Image features: 3 conditions (images), 4 features per image
image_features = [
    0.1  0.2  0.3  0.4;  # Image 1 (condition 1)
    0.2  0.3  0.4  0.5;  # Image 2 (condition 2)
    0.5  0.6  0.7  0.8   # Image 3 (condition 3)
]
rdm = create_rdm_from_matrix(image_features, dissimilarity_measure=:correlation)
```

# Custom Data Conversion

If your data is in a different format, convert it first:

```julia
# Your data might be in a DataFrame
# Convert to Matrix{Float64} format
function convert_my_features(my_df::DataFrame, condition_col::Symbol, feature_cols::Vector{Symbol})
    conditions = unique(my_df[!, condition_col])
    n_conds = length(conditions)
    n_features = length(feature_cols)
    
    matrix = zeros(Float64, n_conds, n_features)
    for (i, cond) in enumerate(conditions)
        cond_data = my_df[my_df[!, condition_col] .== cond, :]
        matrix[i, :] = vec(cond_data[1, feature_cols])  # First row for this condition
    end
    return matrix
end

feature_matrix = convert_my_features(my_data, :condition, [:feat1, :feat2, :feat3])
rdm = create_rdm_from_matrix(feature_matrix)
```
"""
function create_rdm_from_matrix(data_matrix::Matrix{Float64}; dissimilarity_measure::Symbol = :correlation)
    n_conditions = size(data_matrix, 1)
    rdm = zeros(Float64, n_conditions, n_conditions)

    # Only compute upper triangle for efficiency
    for i = 1:n_conditions
        rdm[i, i] = 0.0
        for j = (i+1):n_conditions
            pattern_i = vec(data_matrix[i, :])
            pattern_j = vec(data_matrix[j, :])
            dissim = _compute_dissimilarity(pattern_i, pattern_j, dissimilarity_measure)
            rdm[i, j] = dissim
            rdm[j, i] = dissim  # Symmetric
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
        @minimal_error_throw("Distances vector has length $(length(distances)), expected $expected_length for $n_conditions conditions")
    end

    rdm = zeros(Float64, n_conditions, n_conditions)
    idx = 1

    for i = 1:n_conditions
        for j = (i+1):n_conditions
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
function create_rdm_from_similarity_ratings(similarity_matrix::Matrix{Float64}; convert_to_dissimilarity::Bool = true)
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
    for i = 1:n
        rdm[i, i] = 0.0
    end

    return rdm
end

"""
    create_rdm_from_reaction_times(rts::Vector{Float64})

Create an RDM from reaction times (or any single-value-per-condition data).

# Input Format

**Expected**: `Vector{Float64}` [n_conditions]
- One value per condition
- Order must match condition order in neural data

# Arguments
- `rts::Vector{Float64}`: Reaction times (or other single values) for each condition

# Returns
- `rdm::Matrix{Float64}`: RDM matrix [condition × condition]
  - Dissimilarity = absolute difference between values
  - Diagonal = 0.0 (condition vs itself)

# Examples

```julia
# Reaction times: 3 conditions
rts = [0.3, 0.5, 0.4]  # Face=0.3s, Object=0.5s, Scene=0.4s
rdm = create_rdm_from_reaction_times(rts)

# Resulting RDM:
#            Face  Object  Scene
# Face       0.0   0.2    0.1
# Object     0.2   0.0    0.1
# Scene      0.1   0.1    0.0
```

# Custom Data Conversion

If your data is in a different format, convert it first:

```julia
# Your data might be in a DataFrame
my_data = DataFrame(condition=["Face", "Object", "Scene"], rt=[0.3, 0.5, 0.4])
rts = my_data.rt  # Extract to Vector{Float64}
rdm = create_rdm_from_reaction_times(rts)
```
"""
function create_rdm_from_reaction_times(rts::Vector{Float64})
    n_conditions = length(rts)
    rdm = zeros(Float64, n_conditions, n_conditions)

    # Only compute upper triangle for efficiency
    for i = 1:n_conditions
        rdm[i, i] = 0.0
        for j = (i+1):n_conditions
            dissim = abs(rts[i] - rts[j])
            rdm[i, j] = dissim
            rdm[j, i] = dissim  # Symmetric
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

    # Only compute upper triangle for efficiency
    for i = 1:n_conditions
        rdm[i, i] = 0.0
        for j = (i+1):n_conditions
            dissim = categories[i] == categories[j] ? 0.0 : 1.0
            rdm[i, j] = dissim
            rdm[j, i] = dissim  # Symmetric
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
    temporal_data::Array{Float64,3},
    original_times::Vector{Float64},
    target_times::Vector{Float64};
    method::Symbol = :linear,
)
    n_conditions, n_features, n_original_times = size(temporal_data)
    n_target_times = length(target_times)

    if length(original_times) != n_original_times
        @minimal_error_throw(
            "Original times length ($(length(original_times))) doesn't match " * "temporal data time dimension ($n_original_times)"
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
    for cond_idx = 1:n_conditions
        for feat_idx = 1:n_features
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
                            resampled_series[t_idx] = (1 - weight) * original_series[lower_idx] + weight * original_series[upper_idx]
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
                return resample_temporal_data(temporal_data, original_times, target_times; method = :linear)
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
        align_to::Union{Vector{Float64}, RsaData, Nothing} = nothing,
        interpolation_method::Symbol = :linear,
    )

Create temporal RDMs from temporal model data (e.g., eye tracking, EDA, other physiological signals).

**Automatic temporal alignment**: If your model data has a different sampling rate than neural data,
use `align_to` to automatically resample to match neural timepoints.

# Input Format

**Expected**: 
- `temporal_data::Array{Float64, 3}`: [n_conditions × n_features × n_timepoints]
  - Dimension 1: Conditions (e.g., Face, Object, Scene)
  - Dimension 2: Features (e.g., x_position, y_position, pupil_size for eye tracking)
  - Dimension 3: Time points
- `times::Vector{Float64}`: Time vector in seconds [n_timepoints]
  - Must match length of dimension 3 in `temporal_data`

# Arguments
- `temporal_data::Array{Float64, 3}`: Temporal data [conditions × features × time]
- `times::Vector{Float64}`: Time points in seconds for temporal_data
- `dissimilarity_measure::Symbol`: Measure to use (:correlation, :spearman, :euclidean, :mahalanobis)
- `align_to::Union{Vector{Float64}, RsaData, Nothing}`: Target timepoints to align to
  - `Vector{Float64}`: Time vector (e.g., `neural_rsa.times`)
  - `RsaData`: Automatically uses `neural_rsa.times`
  - `nothing`: Use model's own timepoints (no alignment)
- `interpolation_method::Symbol`: Interpolation method for resampling (:linear, :nearest, :cubic)

# Returns
- `rdms::Array{Float64, 3}`: Temporal RDMs [time × condition × condition]
  - If `align_to` provided: matches length of `align_to` timepoints
  - If `align_to` is `nothing`: matches length of input `times`

# Examples

## Eye Tracking Data

```julia
# Your eye tracking data: 3 conditions, 3 features (x, y, pupil), 100 timepoints at 100 Hz
n_conditions = 3
n_features = 3
n_times = 100
eye_data = zeros(Float64, n_conditions, n_features, n_times)
eye_times = collect(-0.2:(1/100):0.79)  # 100 Hz, -0.2 to 0.79s

# Fill with your data (you write this conversion code)
eye_data[1, 1, :] = face_x_positions  # Condition 1, feature 1 (x), all timepoints
eye_data[1, 2, :] = face_y_positions  # Condition 1, feature 2 (y), all timepoints
eye_data[1, 3, :] = face_pupil_size   # Condition 1, feature 3 (pupil), all timepoints
# ... repeat for other conditions

# Neural data: 500 Hz, -0.2 to 0.8s
neural_rsa = rsa(epochs)

# Automatic alignment to neural data timepoints
eye_rdms = create_temporal_rdm(eye_data, eye_times; align_to=neural_rsa)
# Result: eye_rdms has shape [500 × 3 × 3] matching neural_rsa timepoints
```

## EDA Data

```julia
# Your EDA data: 3 conditions, 1 feature (amplitude), 50 timepoints at 50 Hz
eda_data = zeros(Float64, 3, 1, 50)
eda_times = collect(-0.2:(1/50):0.78)  # 50 Hz

# Fill with your data
eda_data[1, 1, :] = face_eda_signal
eda_data[2, 1, :] = object_eda_signal
eda_data[3, 1, :] = scene_eda_signal

# Align to neural data
eda_rdms = create_temporal_rdm(eda_data, eda_times; align_to=neural_rsa.times)
```

## Custom Data Conversion

If your data is in a different format, convert it first:

```julia
# Your data might be in a custom structure
struct MyEyeData
    conditions::Dict{String, DataFrame}
    times::Vector{Float64}
end

function convert_my_eye_data(my_data::MyEyeData, condition_names::Vector{String})
    n_conds = length(condition_names)
    n_features = 3  # x, y, pupil
    n_times = length(my_data.times)
    
    data = zeros(Float64, n_conds, n_features, n_times)
    for (i, cond_name) in enumerate(condition_names)
        cond_df = my_data.conditions[cond_name]  # Your DataFrame structure
        data[i, 1, :] = cond_df.x_position
        data[i, 2, :] = cond_df.y_position
        data[i, 3, :] = cond_df.pupil_size
    end
    return data, my_data.times
end

# Use your converter
eye_array, eye_times = convert_my_eye_data(my_eye_data, ["Face", "Object", "Scene"])
eye_rdms = create_temporal_rdm(eye_array, eye_times; align_to=neural_rsa)
```

## Without Alignment (Use Model's Own Timepoints)

```julia
# If you want to keep model's original timepoints
rdms = create_temporal_rdm(model_data, model_times)  # No align_to parameter
# Result: rdms has shape [length(model_times) × n_conditions × n_conditions]
```
"""
function create_temporal_rdm(
    temporal_data::Array{Float64,3},
    times::Vector{Float64};
    dissimilarity_measure::Symbol = :correlation,
    align_to::Union{Vector{Float64},RsaData,Nothing} = nothing,
    interpolation_method::Symbol = :linear,
)
    # Handle align_to parameter - accept RsaData or Vector{Float64}
    target_times = nothing
    if !isnothing(align_to)
        if isa(align_to, RsaData)
            target_times = align_to.times
        elseif isa(align_to, Vector{Float64})
            target_times = align_to
        else
            @minimal_error_throw("align_to must be Vector{Float64} or RsaData, got $(typeof(align_to))")
        end
    end

    # Resample if alignment is requested
    if !isnothing(target_times)
        @info "Resampling temporal model data from $(length(times)) to $(length(target_times)) timepoints"
        temporal_data, times = resample_temporal_data(temporal_data, times, target_times; method = interpolation_method)
    end

    n_conditions, n_features, n_timepoints = size(temporal_data)

    if length(times) != n_timepoints
        @minimal_error_throw("Time vector length ($(length(times))) doesn't match temporal data time dimension ($n_timepoints)")
    end

    # Preallocate RDMs: [time × condition × condition]
    rdms = zeros(Float64, n_timepoints, n_conditions, n_conditions)

    # Compute RDM at each time point
    for t = 1:n_timepoints
        # Extract data at time t: [conditions × features]
        timepoint_data = temporal_data[:, :, t]

        # Compute RDM for this time point
        condition_patterns = Matrix{Float64}[]
        for cond_idx = 1:n_conditions
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
            @minimal_error_throw("Condition $cond_idx has $(length(cond_data)) timepoints, expected $n_timepoints")
        end
    end

    # Get number of features from first condition, first timepoint
    n_features = length(temporal_vectors[1][1])

    # Convert to array format: [conditions × features × time]
    temporal_data = zeros(Float64, n_conditions, n_features, n_timepoints)
    for cond_idx = 1:n_conditions
        for t = 1:n_timepoints
            temporal_data[cond_idx, :, t] = temporal_vectors[cond_idx][t]
        end
    end

    return create_temporal_rdm(temporal_data, times; dissimilarity_measure = dissimilarity_measure)
end

"""
    create_temporal_model_rdms(
        temporal_model_data::Dict{String, Any},
        times::Vector{Float64};
        dissimilarity_measure::Symbol = :correlation,
        align_to::Union{Vector{Float64}, RsaData, Nothing} = nothing,
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
- `align_to::Union{Vector{Float64}, RsaData, Nothing}`: Target timepoints to align to
  - `Vector{Float64}`: Time vector (e.g., `neural_rsa.times`)
  - `RsaData`: Automatically uses `neural_rsa.times`
  - `nothing`: Use model's own timepoints (no alignment)
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
neural_rsa = rsa(epochs)  # 500 Hz

temporal_models = Dict(
    "Eye Tracking" => (eye_data, eye_times),  # Tuple with data and its times
    "EDA" => (eda_data, eda_times),
)
rdms, names = create_temporal_model_rdms(
    temporal_models, 
    neural_rsa.times,  # Will be used if model doesn't provide times
    align_to=neural_rsa  # Resample all models to neural timepoints (can pass RsaData directly!)
)
```
"""
function create_temporal_model_rdms(
    temporal_model_data::Dict{String,Any},
    times::Vector{Float64};
    dissimilarity_measure::Symbol = :correlation,
    align_to::Union{Vector{Float64},RsaData,Nothing} = nothing,
    interpolation_method::Symbol = :linear,
)
    model_rdms = Array{Float64,3}[]
    model_names = String[]

    for (name, data) in temporal_model_data
        rdm = nothing
        model_times = times  # Default to provided times

        if isa(data, Tuple) && length(data) == 2
            # Tuple: (data, times) - model has its own timepoints
            model_data, model_times = data
            if !isa(model_data, Array{Float64,3})
                @minimal_error_throw("Model '$name': Tuple first element must be Array{Float64, 3}, got $(typeof(model_data))")
            end
            if !isa(model_times, Vector{Float64})
                @minimal_error_throw("Model '$name': Tuple second element must be Vector{Float64}, got $(typeof(model_times))")
            end
            # Use align_to if provided, otherwise use model's own times
            rdm = create_temporal_rdm(
                model_data,
                model_times;
                dissimilarity_measure = dissimilarity_measure,
                align_to = align_to,
                interpolation_method = interpolation_method,
            )
        elseif isa(data, Array{Float64,3})
            # Direct temporal data [conditions × features × time]
            # Use align_to if provided
            rdm = create_temporal_rdm(
                data,
                times;
                dissimilarity_measure = dissimilarity_measure,
                align_to = align_to,
                interpolation_method = interpolation_method,
            )
        elseif isa(data, Vector{Vector{Vector{Float64}}})
            # Temporal vectors - convert to array first
            @minimal_warning(
                "Model '$name': Vector{Vector{Vector{Float64}}} format doesn't support " *
                "different timepoints. Converting assuming times match."
            )
            # Convert to array format
            n_conditions = length(data)
            n_timepoints = length(data[1])
            n_features = length(data[1][1])
            model_data = zeros(Float64, n_conditions, n_features, n_timepoints)
            for cond_idx = 1:n_conditions
                for t = 1:n_timepoints
                    model_data[cond_idx, :, t] = data[cond_idx][t]
                end
            end
            rdm = create_temporal_rdm(
                model_data,
                times;
                dissimilarity_measure = dissimilarity_measure,
                align_to = align_to,
                interpolation_method = interpolation_method,
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
function create_model_rdms(model_data::Dict{String,Any}; condition_order::Union{Vector{String},Nothing} = nothing)
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

