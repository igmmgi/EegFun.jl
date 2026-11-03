# Helper functions for creating model RDMs from various data types.
#
# Converts different types of "other" data (behavioral, computational, 
# theoretical) into Representational Dissimilarity Matrices (RDMs) for 
# comparison with neural RDMs via a unified `create_rdm` multiple-dispatch interface.
#
# Input Format Specifications
#
# **Important**: Model RDMs expect **summarized/averaged values per condition**, not trial-by-trial data.
# You should average or summarize your behavioral/model data across trials BEFORE creating model RDMs.
#
# ## Static Models (Single Timepoint)
# 
# For data with **one summary value per condition** (e.g., mean reaction times):
# - Expected: `Vector{Float64}` [n_conditions]
# - Example: `create_rdm([0.302, 0.501, 0.403])`
# 
# For data with **categorical labels per condition**:
# - Expected: `Vector{Int}`, `Vector{Symbol}`, or `Vector{String}`
# - Example: `create_rdm([1, 1, 2, 2, 3])` or `create_rdm([:face, :face, :scene])`
# 
# For data with **multiple features per condition** (e.g., word embeddings, image features):
# - Expected: `Vector{Vector{Float64}}` or `Matrix{Float64}` [n_conditions × n_features]
# - Example: `create_rdm(feature_matrix, dissimilarity_measure=:correlation)`
# 
# For **similarity rating matrices**:
# - Expected: square symmetric `Matrix{Float64}`
# - Example: `create_rdm(similarity_matrix; is_similarity=true)`
# 
# For **condensed pairwise distance vectors**:
# - Expected: upper-triangular vector [n*(n-1)/2] + condition count
# - Example: `create_rdm(distances, 3)`
# 
# ## Temporal Models (Multiple Timepoints)
# 
# For data with **many values per condition over time** (e.g., eye tracking, EDA):
# - Expected: `Array{Float64, 3}` [n_conditions × n_features × n_timepoints] + time vector
# - Example: `create_rdm(data, times; align_to=neural_times)`
#
# ---

# === HELPER FUNCTIONS ===

"""
    _fill_rdm_from_dissimilarities!(rdm::Matrix{Float64}, compute_dissim::Function)

Fill an RDM matrix with dissimilarity values.

# Arguments
- `rdm::Matrix{Float64}`: Pre-allocated zero matrix
- `compute_dissim::Function`: Function (i, j) -> dissimilarity value
"""
function _fill_rdm_from_dissimilarities!(rdm::Matrix{Float64}, compute_dissim::Function)
    n = size(rdm, 1)
    for i = 1:n, j = (i+1):n
        dissim = compute_dissim(i, j)
        rdm[i, j] = rdm[j, i] = dissim
    end
    return rdm
end

function _categorical_rdm(categories::AbstractVector)
    n_conditions = length(categories)
    rdm = zeros(Float64, n_conditions, n_conditions)
    for i = 1:n_conditions
        rdm[i, i] = 0.0
        for j = (i+1):n_conditions
            dissim = categories[i] == categories[j] ? 0.0 : 1.0
            rdm[i, j] = dissim
            rdm[j, i] = dissim
        end
    end
    return rdm
end

# === UNIFIED RDM CREATION INTERFACE (create_rdm) ===

"""
    create_rdm(vectors::AbstractVector{<:AbstractVector{<:Real}}; dissimilarity_measure::Symbol = :correlation)

Create an RDM from feature vectors (e.g., word embeddings, image features).

# Input Format
- Outer vector: one element per condition
- Inner vector: feature values for that condition (all inner vectors must have the same length)

# Arguments
- `vectors`: Vector of feature vectors, one per condition
- `dissimilarity_measure::Symbol`: Measure to use (:correlation, :spearman, :euclidean, :mahalanobis)

# Returns
- `rdm::Matrix{Float64}`: RDM matrix [condition × condition]
"""
function create_rdm(vectors::AbstractVector{<:AbstractVector{<:Real}}; dissimilarity_measure::Symbol = :correlation)
    n_conditions = length(vectors)
    rdm = zeros(Float64, n_conditions, n_conditions)
    _fill_rdm_from_dissimilarities!(rdm, (i, j) -> _compute_dissimilarity(vectors[i], vectors[j], dissimilarity_measure))
    return rdm
end

"""
    create_rdm(
        data_matrix::AbstractMatrix{<:Real};
        dissimilarity_measure::Symbol = :correlation,
        is_similarity::Bool = false,
    )

Create an RDM from a feature matrix [conditions × features] or a similarity matrix [conditions × conditions].

# Arguments
- `data_matrix`: Matrix of condition features [conditions × features] or similarity ratings [conditions × conditions]
- `dissimilarity_measure::Symbol`: Measure to use (:correlation, :spearman, :euclidean, :mahalanobis)
- `is_similarity::Bool`: If true, converts similarity ratings to dissimilarity [0, 1] with zero diagonal

# Returns
- `rdm::Matrix{Float64}`: RDM matrix [condition × condition]
"""
function create_rdm(data_matrix::AbstractMatrix{<:Real}; dissimilarity_measure::Symbol = :correlation, is_similarity::Bool = false)
    if is_similarity
        # Normalize to [0, 1] and invert (similarity -> dissimilarity)
        max_val = maximum(data_matrix)
        min_val = minimum(data_matrix)
        if max_val == min_val
            return zeros(Float64, size(data_matrix))
        end
        normalized = (data_matrix .- min_val) ./ (max_val - min_val)
        rdm = 1.0 .- normalized
        for i = 1:size(rdm, 1)
            rdm[i, i] = 0.0
        end
        return rdm
    elseif size(data_matrix, 1) == size(data_matrix, 2) &&
           issymmetric(data_matrix) &&
           all(i -> abs(data_matrix[i, i]) <= 1e-12, axes(data_matrix, 1))
        # Pre-computed symmetric RDM with zero diagonal (including Crossnobis RDMs with negative entries)
        rdm = Matrix{Float64}(data_matrix)
        @inbounds for i in axes(rdm, 1)
            rdm[i, i] = 0.0
        end
        return rdm
    else
        n_conditions = size(data_matrix, 1)
        rdm = zeros(Float64, n_conditions, n_conditions)
        _fill_rdm_from_dissimilarities!(
            rdm,
            (i, j) -> _compute_dissimilarity(@view(data_matrix[i, :]), @view(data_matrix[j, :]), dissimilarity_measure),
        )
        return rdm
    end
end

"""
    create_rdm(distances::AbstractVector{<:Real}, n_conditions::Integer)

Create an RDM from a vector of condensed pairwise distances (upper triangular, length = n*(n-1)/2).
"""
function create_rdm(distances::AbstractVector{<:Real}, n_conditions::Integer)
    expected_length = div(n_conditions * (n_conditions - 1), 2)
    if length(distances) != expected_length
        @minimal_error("Distances vector has length $(length(distances)), expected $expected_length for $n_conditions conditions")
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
    create_rdm(values::AbstractVector{<:AbstractFloat}; method::Symbol = :difference)

Create an RDM from 1D scalar values per condition (e.g., reaction times, subjective ratings).

# Arguments
- `values`: Float vector of summary values, one per condition
- `method::Symbol`: `:difference` (default: absolute difference |v_i - v_j|) or `:categorical`
"""
function create_rdm(values::AbstractVector{<:AbstractFloat}; method::Symbol = :difference)
    if method == :difference
        n_conditions = length(values)
        rdm = zeros(Float64, n_conditions, n_conditions)
        for i = 1:n_conditions
            rdm[i, i] = 0.0
            for j = (i+1):n_conditions
                dissim = abs(values[i] - values[j])
                rdm[i, j] = dissim
                rdm[j, i] = dissim
            end
        end
        return rdm
    elseif method == :categorical
        return _categorical_rdm(values)
    else
        @minimal_error("Unknown method: :$method for float vector. Use :difference or :categorical")
    end
end

"""
    create_rdm(categories::AbstractVector{<:Union{Integer,Symbol,AbstractString}}; method::Symbol = :categorical)

Create an RDM from categorical labels (integers, symbols, or strings).

# Arguments
- `categories`: Category labels for each condition (e.g. `[1, 1, 2]`, `[:face, :face, :scene]`)
- `method::Symbol`: `:categorical` (default: 0 if same, 1 if different) or `:difference` (for integers)
"""
function create_rdm(categories::AbstractVector{<:Union{Integer,Symbol,AbstractString}}; method::Symbol = :categorical)
    if method == :categorical
        return _categorical_rdm(categories)
    elseif method == :difference && categories isa AbstractVector{<:Integer}
        n_conditions = length(categories)
        rdm = zeros(Float64, n_conditions, n_conditions)
        for i = 1:n_conditions
            rdm[i, i] = 0.0
            for j = (i+1):n_conditions
                dissim = abs(Float64(categories[i]) - Float64(categories[j]))
                rdm[i, j] = dissim
                rdm[j, i] = dissim
            end
        end
        return rdm
    else
        @minimal_error("Unsupported method: :$method for $(typeof(categories)). Use :categorical")
    end
end

# === TEMPORAL ALIGNMENT AND RESAMPLING ===
"""
    resample_temporal_data(
        temporal_data::AbstractArray{<:Real, 3},
        original_times::AbstractVector{<:Real},
        target_times::AbstractVector{<:Real};
        method::Symbol = :linear,
    )

Resample temporal model data to match neural data timepoints.

Handles different sampling rates and time ranges by interpolating model data
to match the neural data's temporal structure.

# Arguments
- `temporal_data::AbstractArray{<:Real, 3}`: Temporal data [conditions × features × time]
- `original_times::AbstractVector{<:Real}`: Original time points in seconds (sorted ascending)
- `target_times::AbstractVector{<:Real}`: Target time points (e.g., from neural data)
- `method::Symbol`: Interpolation method (`:linear` or `:nearest`)

# Returns
- `resampled_data::Array{Float64, 3}`: Resampled data [conditions × features × target_time]
- `target_times::Vector{Float64}`: Target time points

# Examples
```julia
# Model data at 100 Hz
model_data = randn(3, 2, 51)  # [conditions × features × time]
model_times = collect(-0.2:(1/100):0.3)

# Neural data at 500 Hz
neural_times = collect(-0.2:(1/500):0.8)

# Resample model data to match neural times
resampled_data, _ = EegFun.resample_temporal_data(model_data, model_times, neural_times)
```
"""
function resample_temporal_data(
    temporal_data::AbstractArray{<:Real,3},
    original_times::AbstractVector{<:Real},
    target_times::AbstractVector{<:Real};
    method::Symbol = :linear,
)
    n_conditions, n_features, n_original_times = size(temporal_data)
    n_target_times = length(target_times)
    target_times_f64 = Vector{Float64}(target_times)

    if length(original_times) != n_original_times
        @minimal_error(
            "Original times length ($(length(original_times))) doesn't match " * "temporal data time dimension ($n_original_times)"
        )
    end
    if n_original_times == 0
        @minimal_error("Cannot resample temporal data with zero original time points")
    end
    if !(method in (:linear, :nearest))
        @minimal_error("Unknown interpolation method: $method. Use :linear or :nearest")
    end

    # Preallocate resampled data
    resampled_data = zeros(Float64, n_conditions, n_features, n_target_times)
    if n_target_times == 0
        return resampled_data, target_times_f64
    end

    # Determine valid time range (intersection of original and target)
    t_orig_first = Float64(first(original_times))
    t_orig_last = Float64(last(original_times))
    valid_start = max(t_orig_first, minimum(target_times_f64))
    valid_end = min(t_orig_last, maximum(target_times_f64))

    if valid_start > valid_end || !any(t -> valid_start <= t <= valid_end, target_times_f64)
        @minimal_warning(
            "No temporal overlap between model data ($t_orig_first to $t_orig_last s) " *
            "and target times ($(first(target_times_f64)) to $(last(target_times_f64)) s). " *
            "Returning zeros."
        )
        return resampled_data, target_times_f64
    end

    # Precompute interpolation indices and weights once outside the condition × feature loops
    if method == :linear
        lower_idxs = Vector{Int}(undef, n_target_times)
        upper_idxs = Vector{Int}(undef, n_target_times)
        weights = Vector{Float64}(undef, n_target_times)

        @inbounds for t_idx = 1:n_target_times
            t_target = target_times_f64[t_idx]
            if t_target <= t_orig_first
                lower_idxs[t_idx] = 1
                upper_idxs[t_idx] = 1
                weights[t_idx] = 0.0
            elseif t_target >= t_orig_last
                lower_idxs[t_idx] = n_original_times
                upper_idxs[t_idx] = n_original_times
                weights[t_idx] = 0.0
            else
                idx = searchsortedlast(original_times, t_target)
                idx = clamp(idx, 1, n_original_times)
                t_lower = Float64(original_times[idx])
                if idx == n_original_times || t_lower == t_target
                    lower_idxs[t_idx] = idx
                    upper_idxs[t_idx] = idx
                    weights[t_idx] = 0.0
                else
                    u_idx = idx + 1
                    t_upper = Float64(original_times[u_idx])
                    denom = t_upper - t_lower
                    lower_idxs[t_idx] = idx
                    upper_idxs[t_idx] = u_idx
                    weights[t_idx] = denom > 0.0 ? (t_target - t_lower) / denom : 0.0
                end
            end
        end

        @inbounds for t_idx = 1:n_target_times
            l_idx = lower_idxs[t_idx]
            u_idx = upper_idxs[t_idx]
            w = weights[t_idx]
            one_minus_w = 1.0 - w
            for feat_idx = 1:n_features
                for cond_idx = 1:n_conditions
                    val =
                        one_minus_w * Float64(temporal_data[cond_idx, feat_idx, l_idx]) +
                        w * Float64(temporal_data[cond_idx, feat_idx, u_idx])
                    resampled_data[cond_idx, feat_idx, t_idx] = val
                end
            end
        end
    else # method == :nearest
        nearest_idxs = Vector{Int}(undef, n_target_times)
        @inbounds for t_idx = 1:n_target_times
            t_target = target_times_f64[t_idx]
            if t_target <= t_orig_first
                nearest_idxs[t_idx] = 1
            elseif t_target >= t_orig_last
                nearest_idxs[t_idx] = n_original_times
            else
                idx = clamp(searchsortedlast(original_times, t_target), 1, n_original_times)
                if idx == n_original_times
                    nearest_idxs[t_idx] = n_original_times
                else
                    d_low = abs(t_target - Float64(original_times[idx]))
                    d_high = abs(Float64(original_times[idx+1]) - t_target)
                    nearest_idxs[t_idx] = d_low <= d_high ? idx : (idx + 1)
                end
            end
        end

        @inbounds for t_idx = 1:n_target_times
            n_idx = nearest_idxs[t_idx]
            for feat_idx = 1:n_features
                for cond_idx = 1:n_conditions
                    resampled_data[cond_idx, feat_idx, t_idx] = Float64(temporal_data[cond_idx, feat_idx, n_idx])
                end
            end
        end
    end

    return resampled_data, target_times_f64
end
"""
    create_rdm(
        temporal_data::AbstractArray{<:Real, 3},
        times::AbstractVector{<:Real};
        dissimilarity_measure::Symbol = :correlation,
        align_to::Union{AbstractVector{<:Real}, RsaData, Nothing} = nothing,
        interpolation_method::Symbol = :linear,
    )

Create temporal RDMs from temporal model data (e.g., eye tracking, EDA, other physiological signals).

**Automatic temporal alignment**: If your model data has a different sampling rate than neural data,
use `align_to` to automatically resample to match neural timepoints.

# Input Format
- `temporal_data`: [conditions × features × timepoints]
- `times`: Time vector in seconds [timepoints]

# Arguments
- `temporal_data`: Temporal data [conditions × features × time]
- `times`: Time points in seconds for temporal_data
- `dissimilarity_measure::Symbol`: Measure to use (:correlation, :spearman, :euclidean, :mahalanobis)
- `align_to`: Target timepoints to align to (`Vector{Float64}` or `RsaData`, or `nothing`)
- `interpolation_method::Symbol`: Interpolation method for resampling (`:linear` or `:nearest`)

# Returns
- `rdms::Array{Float64, 3}`: Temporal RDMs [time × condition × condition]
"""
function create_rdm(
    temporal_data::AbstractArray{<:Real,3},
    times::AbstractVector{<:Real};
    dissimilarity_measure::Symbol = :correlation,
    align_to::Union{AbstractVector{<:Real},RsaData,Nothing} = nothing,
    interpolation_method::Symbol = :linear,
)
    # Handle align_to parameter - accept RsaData or AbstractVector{<:Real}
    target_times = nothing
    if !isnothing(align_to)
        if isa(align_to, RsaData)
            target_times = align_to.times
        elseif isa(align_to, AbstractVector{<:Real})
            target_times = Vector{Float64}(align_to)
        else
            @minimal_error("align_to must be Vector{<:Real} or RsaData, got $(typeof(align_to))")
        end
    end

    # Resample if alignment is requested
    t_data = Array{Float64,3}(temporal_data)
    t_times = Vector{Float64}(times)
    if !isnothing(target_times)
        @info "Resampling temporal model data from $(length(times)) to $(length(target_times)) timepoints"
        t_data, t_times = resample_temporal_data(t_data, t_times, target_times; method = interpolation_method)
    end

    n_conditions, n_features, n_timepoints = size(t_data)

    if length(t_times) != n_timepoints
        @minimal_error("Time vector length ($(length(t_times))) doesn't match temporal data time dimension ($n_timepoints)")
    end

    # Preallocate RDMs: [time × condition × condition]
    rdms = zeros(Float64, n_timepoints, n_conditions, n_conditions)

    # Compute RDM at each time point
    for t = 1:n_timepoints
        timepoint_data = @view t_data[:, :, t]
        condition_patterns = [vec(timepoint_data[cond_idx, :]) for cond_idx = 1:n_conditions]
        rdms[t, :, :] = _compute_rdm(condition_patterns, dissimilarity_measure)
    end

    return rdms
end

function create_rdm(
    temporal_vectors::AbstractVector{<:AbstractVector{<:AbstractVector{<:Real}}},
    times::AbstractVector{<:Real};
    dissimilarity_measure::Symbol = :correlation,
    align_to::Union{AbstractVector{<:Real},RsaData,Nothing} = nothing,
    interpolation_method::Symbol = :linear,
)
    n_conditions = length(temporal_vectors)
    n_timepoints = length(times)

    for (cond_idx, cond_data) in enumerate(temporal_vectors)
        if length(cond_data) != n_timepoints
            @minimal_error("Condition $cond_idx has $(length(cond_data)) timepoints, expected $n_timepoints")
        end
    end

    n_features = length(temporal_vectors[1][1])
    temporal_data = zeros(Float64, n_conditions, n_features, n_timepoints)
    for cond_idx = 1:n_conditions
        for t = 1:n_timepoints
            temporal_data[cond_idx, :, t] = temporal_vectors[cond_idx][t]
        end
    end

    return create_rdm(
        temporal_data,
        times;
        dissimilarity_measure = dissimilarity_measure,
        align_to = align_to,
        interpolation_method = interpolation_method,
    )
end

"""
    create_temporal_model_rdms(
        temporal_model_data::AbstractDict{<:AbstractString, <:Any},
        times::AbstractVector{<:Real};
        dissimilarity_measure::Symbol = :correlation,
        align_to::Union{AbstractVector{<:Real}, RsaData, Nothing} = nothing,
        interpolation_method::Symbol = :linear,
    )

Create temporal model RDMs from a dictionary of temporal data.

Automatically detects the format and converts to temporal RDMs. If `align_to` is provided,
model data will be resampled to match those timepoints (handles different sampling rates).

# Arguments
- `temporal_model_data`: Dictionary mapping model names to temporal data
  - `Tuple{Array{<:Real, 3}, Vector{<:Real}}`: (data, times) - data with its own timepoints
  - `Array{<:Real, 3}`: Direct temporal data (uses `times` parameter)
  - `Vector{Vector{Vector{Float64}}}`: Temporal vectors
- `times`: Time points in seconds (used if model data doesn't provide its own)
- `dissimilarity_measure::Symbol`: Measure to use
- `align_to`: Target timepoints to align to (`Vector{Float64}` or `RsaData`, or `nothing`)
- `interpolation_method::Symbol`: Interpolation method for resampling (`:linear` or `:nearest`)

# Returns
- `model_rdms::Vector{Array{Float64, 3}}`: Vector of temporal RDMs [time × condition × condition]
- `model_names::Vector{String}`: Vector of model names
"""
function create_temporal_model_rdms(
    temporal_model_data::AbstractDict{<:AbstractString,<:Any},
    times::AbstractVector{<:Real};
    dissimilarity_measure::Symbol = :correlation,
    align_to::Union{AbstractVector{<:Real},RsaData,Nothing} = nothing,
    interpolation_method::Symbol = :linear,
)
    model_rdms = Array{Float64,3}[]
    model_names = String[]

    entries = temporal_model_data isa Dict ? sort!(collect(pairs(temporal_model_data)); by = first) : collect(pairs(temporal_model_data))
    for (name, data) in entries
        rdm = nothing
        model_times = times  # Default to provided times

        if isa(data, Tuple) && length(data) == 2
            model_data, model_times = data
            if !isa(model_data, AbstractArray{<:Real,3})
                @minimal_error("Model '$name': Tuple first element must be Array{<:Real, 3}, got $(typeof(model_data))")
            end
            if !isa(model_times, AbstractVector{<:Real})
                @minimal_error("Model '$name': Tuple second element must be Vector{<:Real}, got $(typeof(model_times))")
            end
            rdm = create_rdm(
                model_data,
                model_times;
                dissimilarity_measure = dissimilarity_measure,
                align_to = align_to,
                interpolation_method = interpolation_method,
            )
        elseif isa(data, AbstractArray{<:Real,3}) || isa(data, AbstractVector{<:AbstractVector{<:AbstractVector{<:Real}}})
            rdm = create_rdm(
                data,
                times;
                dissimilarity_measure = dissimilarity_measure,
                align_to = align_to,
                interpolation_method = interpolation_method,
            )
        else
            @minimal_error(
                "Unknown temporal model data type for '$name': $(typeof(data)). " *
                "Supported types: Array{Float64, 3} [conditions × features × time], " *
                "Tuple{Array{Float64, 3}, Vector{Float64}} (data, times), " *
                "Vector{Vector{Vector{Float64}}}"
            )
        end

        push!(model_rdms, rdm)
        push!(model_names, String(name))
    end

    return model_rdms, model_names
end

# === CONVENIENCE FUNCTION FOR MULTIPLE MODEL TYPES ===

"""
    create_model_rdms(model_data::AbstractDict{<:AbstractString, <:Any})

Create multiple model RDMs from a dictionary of different model types.

This function automatically detects the type of each model and converts it to an RDM using `EegFun.create_rdm`.

# Arguments
- `model_data`: Dictionary mapping model names to their data
  - `Vector{<:AbstractVector}`: Feature vectors → RDM
  - `Matrix{<:Real}`: Data matrix [conditions × features] or pre-computed RDM
  - `Vector{<:Real}`: 1D values (reaction times, ratings, or category labels)
  - `Vector{<:Union{Symbol, AbstractString}}`: Categorical labels

# Returns
- `model_rdms::Vector{Matrix{Float64}}`: Vector of model RDMs
- `model_names::Vector{String}`: Vector of model names
"""
function create_model_rdms(model_data::AbstractDict{<:AbstractString,<:Any})
    model_rdms = Matrix{Float64}[]
    model_names = String[]

    entries = model_data isa Dict ? sort!(collect(pairs(model_data)); by = first) : collect(pairs(model_data))
    for (name, data) in entries
        try
            rdm = create_rdm(data)
            push!(model_rdms, rdm)
            push!(model_names, String(name))
        catch e
            e isa EegFunError && rethrow(e)
            @minimal_error(
                "Failed to create RDM for model '$name' with data type $(typeof(data)): $(sprint(showerror, e)). " *
                "Supported types: Vector{<:AbstractVector}, Matrix, Vector{Float}, Vector{Int/Symbol/String}"
            )
        end
    end

    return model_rdms, model_names
end

"""
    create_rdm(model_data::AbstractDict{<:AbstractString, <:Any})

Batch-create model RDMs from a dictionary of models. Alias for `create_model_rdms`.
"""
create_rdm(model_data::AbstractDict{<:AbstractString,<:Any}) = create_model_rdms(model_data)
