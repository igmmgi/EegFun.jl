# Model RDM construction and temporal model alignment for RSA.

function _fill_rdm_from_dissimilarities!(rdm::Matrix{Float64}, compute_dissim)
    n = size(rdm, 1)
    @inbounds for i = 1:n
        rdm[i, i] = 0.0
        for j = (i+1):n
            dissim = Float64(compute_dissim(i, j))
            rdm[i, j] = dissim
            rdm[j, i] = dissim
        end
    end
    return rdm
end

function _categorical_rdm(categories::AbstractVector)
    rdm = zeros(Float64, length(categories), length(categories))
    return _fill_rdm_from_dissimilarities!(rdm, (i, j) -> categories[i] == categories[j] ? 0.0 : 1.0)
end

function _difference_rdm(values::AbstractVector)
    rdm = zeros(Float64, length(values), length(values))
    return _fill_rdm_from_dissimilarities!(rdm, (i, j) -> abs(Float64(values[i]) - Float64(values[j])))
end

_ordered_entries(d::AbstractDict) = d isa Dict ? sort!(collect(pairs(d)); by = first) : collect(pairs(d))

"""
    create_rdm(vectors::AbstractVector{<:AbstractVector{<:Real}}; dissimilarity_measure::Symbol = :correlation)
    create_rdm(data_matrix::AbstractMatrix{<:Real}; dissimilarity_measure::Symbol = :correlation, is_similarity::Bool = false)
    create_rdm(distances::AbstractVector{<:Real}, n_conditions::Integer)
    create_rdm(values::AbstractVector{<:AbstractFloat}; method::Symbol = :difference)
    create_rdm(categories::AbstractVector{<:Union{Integer,Symbol,AbstractString}}; method::Symbol = :categorical)
    create_rdm(temporal_data::AbstractArray{<:Real,3}, times::AbstractVector{<:Real};
               dissimilarity_measure::Symbol = :correlation, align_to = nothing,
               interpolation_method::Symbol = :linear)
    create_rdm(model_data::AbstractDict{<:AbstractString,<:Any})

Construct a Representational Dissimilarity Matrix (`Matrix{Float64}`) or temporal RDM series
(`Array{Float64,3}` `[time × condition × condition]`) from behavioral, categorical, feature,
similarity, or time-resolved model data.
"""
function create_rdm(vectors::AbstractVector{<:AbstractVector{<:Real}}; dissimilarity_measure::Symbol = :correlation)
    n_conditions = length(vectors)
    rdm = zeros(Float64, n_conditions, n_conditions)
    vf = [v isa AbstractVector{Float64} ? v : Float64.(v) for v in vectors]
    return _fill_rdm_from_dissimilarities!(rdm, (i, j) -> _compute_dissimilarity(vf[i], vf[j], dissimilarity_measure))
end

function create_rdm(data_matrix::AbstractMatrix{<:Real}; dissimilarity_measure::Symbol = :correlation, is_similarity::Bool = false)
    if is_similarity
        min_val, max_val = extrema(data_matrix)
        max_val == min_val && return zeros(Float64, size(data_matrix))
        inv_span = 1.0 / (max_val - min_val)
        n = size(data_matrix, 1)
        rdm = Matrix{Float64}(undef, size(data_matrix))
        @inbounds for j in axes(data_matrix, 2), i in axes(data_matrix, 1)
            rdm[i, j] = i == j && size(data_matrix, 1) == size(data_matrix, 2) ? 0.0 : 1.0 - (data_matrix[i, j] - min_val) * inv_span
        end
        return rdm
    elseif size(data_matrix, 1) == size(data_matrix, 2) &&
           issymmetric(data_matrix) &&
           all(i -> abs(data_matrix[i, i]) <= 1e-12, axes(data_matrix, 1))
        # Preserve pre-computed symmetric RDMs (including Crossnobis RDMs with negative off-diagonal entries)
        rdm = Matrix{Float64}(data_matrix)
        @inbounds for i in axes(rdm, 1)
            rdm[i, i] = 0.0
        end
        return rdm
    else
        mat_f64 = data_matrix isa AbstractMatrix{Float64} ? data_matrix : Matrix{Float64}(data_matrix)
        n_conditions = size(mat_f64, 1)
        rdm = zeros(Float64, n_conditions, n_conditions)
        return _fill_rdm_from_dissimilarities!(
            rdm,
            (i, j) -> _compute_dissimilarity(@view(mat_f64[i, :]), @view(mat_f64[j, :]), dissimilarity_measure),
        )
    end
end

function create_rdm(distances::AbstractVector{<:Real}, n_conditions::Integer)
    expected_length = div(n_conditions * (n_conditions - 1), 2)
    length(distances) != expected_length &&
        @minimal_error("Distances vector has length $(length(distances)), expected $expected_length for $n_conditions conditions")

    rdm = zeros(Float64, n_conditions, n_conditions)
    idx = 1
    @inbounds for i = 1:n_conditions, j = (i+1):n_conditions
        v = Float64(distances[idx])
        rdm[i, j] = v
        rdm[j, i] = v
        idx += 1
    end
    return rdm
end

function create_rdm(values::AbstractVector{<:AbstractFloat}; method::Symbol = :difference)
    method == :difference && return _difference_rdm(values)
    method == :categorical && return _categorical_rdm(values)
    @minimal_error("Unknown method: :$method for float vector. Use :difference or :categorical")
end

function create_rdm(categories::AbstractVector{<:Union{Integer,Symbol,AbstractString}}; method::Symbol = :categorical)
    method == :categorical && return _categorical_rdm(categories)
    (method == :difference && categories isa AbstractVector{<:Integer}) && return _difference_rdm(categories)
    @minimal_error("Unsupported method: :$method for $(typeof(categories)). Use :categorical")
end

"""
    resample_temporal_data(
        temporal_data::AbstractArray{<:Real,3},
        original_times::AbstractVector{<:Real},
        target_times::AbstractVector{<:Real};
        method::Symbol = :linear,
    ) -> (Array{Float64,3}, Vector{Float64})

Interpolate temporal model data `[conditions × features × time]` from sorted `original_times`
onto `target_times` using `:linear` or `:nearest` interpolation (clamped outside overlap; returns
zeros with a warning when time ranges do not overlap).
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

    length(original_times) != n_original_times &&
        @minimal_error("Original times length ($(length(original_times))) doesn't match temporal data time dimension ($n_original_times)")
    n_original_times == 0 && @minimal_error("Cannot resample temporal data with zero original time points")
    method in (:linear, :nearest) || @minimal_error("Unknown interpolation method: $method. Use :linear or :nearest")

    resampled_data = zeros(Float64, n_conditions, n_features, n_target_times)
    n_target_times == 0 && return resampled_data, target_times_f64

    t_orig_first = Float64(first(original_times))
    t_orig_last = Float64(last(original_times))
    valid_start = max(t_orig_first, minimum(target_times_f64))
    valid_end = min(t_orig_last, maximum(target_times_f64))

    if valid_start > valid_end || !any(t -> valid_start <= t <= valid_end, target_times_f64)
        @minimal_warning(
            "No temporal overlap between model data ($t_orig_first to $t_orig_last s) " *
            "and target times ($(first(target_times_f64)) to $(last(target_times_f64)) s). Returning zeros."
        )
        return resampled_data, target_times_f64
    end

    if method == :linear
        @inbounds for t_idx = 1:n_target_times
            t_target = target_times_f64[t_idx]
            l_idx, u_idx, w = if t_target <= t_orig_first
                (1, 1, 0.0)
            elseif t_target >= t_orig_last
                (n_original_times, n_original_times, 0.0)
            else
                idx = clamp(searchsortedlast(original_times, t_target), 1, n_original_times)
                t_lower = Float64(original_times[idx])
                if idx == n_original_times || t_lower == t_target
                    (idx, idx, 0.0)
                else
                    denom = Float64(original_times[idx+1]) - t_lower
                    (idx, idx + 1, denom > 0.0 ? (t_target - t_lower) / denom : 0.0)
                end
            end
            one_minus_w = 1.0 - w
            for feat_idx = 1:n_features, cond_idx = 1:n_conditions
                resampled_data[cond_idx, feat_idx, t_idx] =
                    one_minus_w * Float64(temporal_data[cond_idx, feat_idx, l_idx]) + w * Float64(temporal_data[cond_idx, feat_idx, u_idx])
            end
        end
    else
        @inbounds for t_idx = 1:n_target_times
            t_target = target_times_f64[t_idx]
            n_idx = if t_target <= t_orig_first
                1
            elseif t_target >= t_orig_last
                n_original_times
            else
                idx = clamp(searchsortedlast(original_times, t_target), 1, n_original_times)
                if idx == n_original_times
                    n_original_times
                else
                    abs(t_target - Float64(original_times[idx])) <= abs(Float64(original_times[idx+1]) - t_target) ? idx : (idx + 1)
                end
            end
            for feat_idx = 1:n_features, cond_idx = 1:n_conditions
                resampled_data[cond_idx, feat_idx, t_idx] = Float64(temporal_data[cond_idx, feat_idx, n_idx])
            end
        end
    end

    return resampled_data, target_times_f64
end

function create_rdm(
    temporal_data::AbstractArray{<:Real,3},
    times::AbstractVector{<:Real};
    dissimilarity_measure::Symbol = :correlation,
    align_to::Union{AbstractVector{<:Real},RsaData,Nothing} = nothing,
    interpolation_method::Symbol = :linear,
)
    target_times = if isnothing(align_to)
        nothing
    elseif align_to isa RsaData
        align_to.times
    else
        Vector{Float64}(align_to)
    end

    t_data = temporal_data isa Array{Float64,3} ? temporal_data : Array{Float64,3}(temporal_data)
    t_times = times isa Vector{Float64} ? times : Vector{Float64}(times)
    if !isnothing(target_times)
        @info "Resampling temporal model data from $(length(times)) to $(length(target_times)) timepoints"
        t_data, t_times = resample_temporal_data(t_data, t_times, target_times; method = interpolation_method)
    end

    n_conditions, n_features, n_timepoints = size(t_data)
    length(t_times) != n_timepoints &&
        @minimal_error("Time vector length ($(length(t_times))) doesn't match temporal data time dimension ($n_timepoints)")

    rdms = zeros(Float64, n_timepoints, n_conditions, n_conditions)
    rdm_buf = zeros(Float64, n_conditions, n_conditions)
    condition_patterns = [Vector{Float64}(undef, n_features) for _ = 1:n_conditions]

    @inbounds for t = 1:n_timepoints
        for c = 1:n_conditions
            cp = condition_patterns[c]
            @simd for f = 1:n_features
                cp[f] = t_data[c, f, t]
            end
        end
        _compute_rdm!(rdm_buf, condition_patterns, dissimilarity_measure)
        for j = 1:n_conditions, i = 1:n_conditions
            rdms[t, i, j] = rdm_buf[i, j]
        end
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
        length(cond_data) != n_timepoints &&
            @minimal_error("Condition $cond_idx has $(length(cond_data)) timepoints, expected $n_timepoints")
    end

    n_features = length(temporal_vectors[1][1])
    temporal_data = zeros(Float64, n_conditions, n_features, n_timepoints)
    for cond_idx = 1:n_conditions, t = 1:n_timepoints
        temporal_data[cond_idx, :, t] = temporal_vectors[cond_idx][t]
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
        temporal_model_data::AbstractDict{<:AbstractString,<:Any},
        times::AbstractVector{<:Real};
        dissimilarity_measure::Symbol = :correlation,
        align_to = nothing,
        interpolation_method::Symbol = :linear,
    ) -> (Vector{Array{Float64,3}}, Vector{String})

Batch-create temporal model RDMs from a dictionary mapping model names to `Array{<:Real,3}`,
`(Array{<:Real,3}, times)` tuples, or nested temporal vectors.
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

    for (name, data) in _ordered_entries(temporal_model_data)
        rdm = if data isa Tuple && length(data) == 2
            model_data, model_times = data
            model_data isa AbstractArray{<:Real,3} ||
                @minimal_error("Model '$name': Tuple first element must be Array{<:Real, 3}, got $(typeof(model_data))")
            model_times isa AbstractVector{<:Real} ||
                @minimal_error("Model '$name': Tuple second element must be Vector{<:Real}, got $(typeof(model_times))")
            create_rdm(
                model_data,
                model_times;
                dissimilarity_measure = dissimilarity_measure,
                align_to = align_to,
                interpolation_method = interpolation_method,
            )
        elseif data isa AbstractArray{<:Real,3} || data isa AbstractVector{<:AbstractVector{<:AbstractVector{<:Real}}}
            create_rdm(
                data,
                times;
                dissimilarity_measure = dissimilarity_measure,
                align_to = align_to,
                interpolation_method = interpolation_method,
            )
        else
            @minimal_error(
                "Unknown temporal model data type for '$name': $(typeof(data)). " *
                "Supported types: Array{Float64, 3}, Tuple{Array{Float64, 3}, Vector{Float64}}, Vector{Vector{Vector{Float64}}}"
            )
        end

        push!(model_rdms, rdm)
        push!(model_names, String(name))
    end

    return model_rdms, model_names
end

"""
    create_model_rdms(model_data::AbstractDict{<:AbstractString,<:Any}) -> (Vector{Matrix{Float64}}, Vector{String})

Batch-create static model RDMs from a dictionary mapping model names to model inputs accepted by [`create_rdm`](@ref).
"""
function create_model_rdms(model_data::AbstractDict{<:AbstractString,<:Any})
    model_rdms = Matrix{Float64}[]
    model_names = String[]

    for (name, data) in _ordered_entries(model_data)
        try
            push!(model_rdms, create_rdm(data))
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

create_rdm(model_data::AbstractDict{<:AbstractString,<:Any}) = create_model_rdms(model_data)
