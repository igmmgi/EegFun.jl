"""
    run_ica(dat::ContinuousData; n_components = nothing, sample_selection = samples(),
    interval_selection = times(), channel_selection = channels(), include_extra = false,
    percentage_of_data = 100.0, algorithm = :infomax, params = IcaPrms())
    run_ica(epoched_data::Vector{EpochData}; n_components = nothing, sample_selection = samples(),
    channel_selection = channels(), remove_duplicates = true, percentage_of_data = 100.0,
    algorithm = :infomax, params = IcaPrms())

Run Independent Component Analysis (ICA) on EEG data.

The `ContinuousData` form runs ICA on a single recording. The `Vector{EpochData}` form
concatenates all epochs into a continuous matrix before decomposition (the standard approach
for epoched data, since ICA benefits from many data points).

Tip: **Pre-filter to ≥1 Hz** before running ICA to ensure a good decomposition.

# Key Arguments
- `n_components`: Number of components (default: n_channels − 1)
- `sample_selection`: Exclude bad samples (e.g., `samples_not(:is_extreme_value_100)`)
- `percentage_of_data`: Use a random subset of data for faster ICA (e.g., `25.0` = 25%)
- `algorithm`: `:infomax` (default) or `:infomax_extended` (handles sub- and super-Gaussian sources)
- `remove_duplicates` *(Vector form only)*: Remove duplicate samples across conditions (default: true)

# Returns
`InfoIca` with unmixing, mixing, sphere, variance, and metadata.

# Examples
```julia
# Continuous data
ica_result = run_ica(dat)
ica_result = run_ica(dat, sample_selection = samples_not(:is_extreme_value_100))
ica_result = run_ica(dat, percentage_of_data = 25.0)

# Epoched data
ica_result = run_ica(epochs_unrejected)
```
"""
function run_ica(
    dat::ContinuousData;
    n_components::Union{Nothing,Int} = nothing,
    sample_selection::Function = samples(),
    interval_selection::Interval = times(),
    channel_selection::Function = channels(),
    include_extra::Bool = false,
    percentage_of_data::Real = 100.0,
    algorithm::Symbol = :infomax,
    use_gpu::Bool = false,
    params::IcaPrms = IcaPrms(),
)
    params.use_gpu = use_gpu || params.use_gpu
    selected_channels = get_selected_channels(dat, channel_selection; include_meta = false, include_extra = include_extra)
    isempty(selected_channels) && error("No channels available after applying channel filter")

    # Combine interval and sample selection 
    combined_sel = _combine_interval_sample(interval_selection, sample_selection)

    # Get samples to use using predicate
    sample_indices = get_selected_samples(dat, combined_sel)
    isempty(sample_indices) && error("No samples available after applying sample filter")

    # Set n_components if not specified
    if isnothing(n_components)
        n_components = length(selected_channels) - 1
    elseif n_components > length(selected_channels)
        @minimal_warning "Requested $n_components components but only $(length(selected_channels)) channels available. Using $(length(selected_channels) - 1) components instead."
        n_components = length(selected_channels) - 1
    end

    @info "Running ICA (use_gpu=$(params.use_gpu)): $(length(selected_channels)) channels x $(length(sample_indices)) samples -> $(n_components) components"

    # Create subsetted layout that matches the selected channels
    ica_layout = subset_layout(dat.layout, channel_selection = channels(selected_channels))

    # Create data matrix and run ICA
    dat_for_ica = _create_ica_data_matrix(dat.data, selected_channels, sample_indices)

    # Subset data if requested
    if percentage_of_data != 100
        dat_for_ica = _select_subsample!(dat_for_ica, percentage_of_data)
    end

    # Get filename from the data
    input_filename = filename(dat)

    # Dispatch to the appropriate ICA algorithm
    ica_result =
        _run_ica_algorithm(dat_for_ica, ica_layout, input_filename; n_components = n_components, algorithm = algorithm, params = params)

    return ica_result
end



function run_ica(
    epoched_data::Vector{EpochData};
    n_components::Union{Nothing,Int} = nothing,
    sample_selection::Function = samples(),
    interval_selection::Interval = times(),
    channel_selection::Function = channels(),
    include_extra::Bool = false,
    remove_duplicates::Bool = true,
    percentage_of_data::Real = 100.0,
    algorithm::Symbol = :infomax,
    use_gpu::Bool = false,
    params::IcaPrms = IcaPrms(),
)
    params.use_gpu = use_gpu || params.use_gpu
    isempty(epoched_data) && error("Empty epoched_data vector provided")

    # Use the first EpochData object as reference for some meta-like data
    reference_epoch_data = epoched_data[1]
    for (i, epoch_data) in enumerate(epoched_data)
        if epoch_data.sample_rate != reference_epoch_data.sample_rate
            error(
                "Inconsistent sample rates: EpochData $i has $(epoch_data.sample_rate) Hz, expected $(reference_epoch_data.sample_rate) Hz",
            )
        end
    end

    # Get channel information from reference
    selected_channels = get_selected_channels(reference_epoch_data, channel_selection; include_meta = false, include_extra = include_extra)
    isempty(selected_channels) && error("No channels available after applying channel filter")

    # Concatenate all epoched data and check for duplicates
    concatenated_df = all_data(epoched_data)
    _check_epoched_data_uniqueness!(concatenated_df; remove_duplicates = remove_duplicates)

    # Combine interval and sample selection (consistent with subset() pattern)
    combined_sel = _combine_interval_sample(interval_selection, sample_selection)

    sample_indices = get_selected_samples(concatenated_df, combined_sel)
    isempty(sample_indices) && error("No samples available after applying sample filter to epoched data")

    # Create data matrix for ICA
    concatenated_matrix = _create_ica_data_matrix(concatenated_df, selected_channels, sample_indices)
    if percentage_of_data != 100
        concatenated_matrix = _select_subsample!(concatenated_matrix, percentage_of_data)
    end

    # Set n_components if not specified
    if isnothing(n_components)
        n_components = length(selected_channels) - 1
    elseif n_components > length(selected_channels)
        @minimal_warning "Requested $n_components components but only $(length(selected_channels)) channels available. Using $(length(selected_channels) - 1) components instead."
        n_components = length(selected_channels) - 1
    end

    final_samples = size(concatenated_matrix, 2)
    total_epochs = sum(length(epoch_data.data) for epoch_data in epoched_data)

    @info "Running ICA on concatenated epochs: $(length(selected_channels)) channels x $final_samples samples (from $total_epochs epochs) -> $n_components components"

    # Create subsetted layout that matches the selected channels  
    ica_layout = subset_layout(reference_epoch_data.layout, channel_selection = channels(selected_channels))

    # Get filename from the first epoch data
    input_filename = filename(reference_epoch_data)

    # Dispatch to the appropriate ICA algorithm with concatenated data
    ica_result = _run_ica_algorithm(
        concatenated_matrix,
        ica_layout,
        input_filename;
        n_components = n_components,
        algorithm = algorithm,
        params = params,
    )

    return ica_result
end

"""
    _check_epoched_data_uniqueness!(concatenated_df::DataFrame; remove_duplicates::Bool = false)

Check for duplicate samples in concatenated epoched data using the samples column.

The `samples` column contains the original sample indices from the continuous data.
When epochs have overlapping time windows, the same original samples appear multiple 
times in the concatenated data, which can bias ICA decomposition.

# Arguments
- `concatenated_df::DataFrame`: Concatenated epoch DataFrame
- `remove_duplicates::Bool`: Whether to remove duplicate samples (default: false)

# Effects
- Warns user if duplicates are found based on samples column
- Modifies DataFrame in place if remove_duplicates=true
"""
function _check_epoched_data_uniqueness!(concatenated_df::DataFrame; remove_duplicates::Bool = false)

    # Check if samples column exists
    if !hasproperty(concatenated_df, :samples)
        @debug "No samples column found - cannot check for duplicate samples"
        return nothing
    end

    # Check for duplicate sample indices
    n_original = nrow(concatenated_df)
    unique_samples = length(unique(concatenated_df.samples))
    n_duplicates = n_original - unique_samples

    if n_duplicates > 0
        duplicate_percentage = round(100 * n_duplicates / n_original, digits = 1)

        @minimal_warning """
        Found $n_duplicates duplicate samples ($duplicate_percentage%) in concatenated epoched data.
        This occurs when epochs have overlapping time windows - the same original data samples 
        appear multiple times. Duplicate samples may bias ICA decomposition.

        Total rows: $n_original
        Unique samples: $unique_samples  
        Duplicates: $n_duplicates ($duplicate_percentage%)
        """

        if remove_duplicates
            unique!(concatenated_df, :samples)
            @info "Automatically removed $n_duplicates duplicate samples (using $unique_samples unique samples)."
        end
    end

    return nothing
end

run_ica(epoched_data::EpochData; kwargs...) = run_ica([epoched_data]; kwargs...)




"""Extract selected channels and samples from a DataFrame into a (channels × samples) matrix."""
function _create_ica_data_matrix(dat::DataFrame, channels, samples)
    # Filter to only existing channels (matching original behavior)
    existing_channels = intersect(propertynames(dat), channels)
    n_channels = length(existing_channels)
    n_samples = length(samples)

    # Pre-allocate result matrix
    result = Matrix{Float64}(undef, n_channels, n_samples)

    # Use direct column access for better performance
    for (i, ch) in enumerate(existing_channels)
        col = dat[!, ch]::Vector{Float64}
        @inbounds @simd for j in eachindex(samples)
            result[i, j] = col[samples[j]]
        end
    end

    return result
end

"""
    _select_subsample!(data_matrix::Matrix{Float64}, percentage::Real)

Apply random subsampling to the data matrix for ICA speedup.
Modifies the matrix in place by returning a subsampled view.

# Arguments
- `data_matrix::Matrix{Float64}`: ICA data matrix (channels × samples)  
- `percentage::Real`: Percentage of samples to keep (0 < percentage <= 100)

# Returns
- `Matrix{Float64}`: Subsampled matrix with fewer columns
"""
function _select_subsample!(data_matrix::Matrix{Float64}, percentage::Real)
    (percentage <= 0 || percentage > 100) && error("percentage_of_data must be between 0 and 100, got $percentage")

    n_original_samples = size(data_matrix, 2)
    n_target_samples = round(Int, n_original_samples * percentage / 100)

    if n_target_samples == n_original_samples
        @info "Random sample: $n_target_samples of $n_original_samples ($(round(percentage, digits=1))%)"
        return data_matrix
    end

    # Random sample selection (columns are samples in data_matrix)
    sample_cols = randperm(n_original_samples)[1:n_target_samples]
    subsampled_matrix = data_matrix[:, sort(sample_cols)]

    @info "Random sample: $n_target_samples of $n_original_samples ($(round(percentage, digits=1))%)"

    return subsampled_matrix
end


"""Pre-allocated work arrays for the Infomax ICA iteration loop."""
mutable struct WorkArrays
    weights::Matrix{Float64}
    u::Matrix{Float64}
    y::Matrix{Float64}
    data_block::Matrix{Float64}
    oldweights::Matrix{Float64}
    startweights::Matrix{Float64}
    weights_temp::Matrix{Float64}
    wu_term::Matrix{Float64}
    delta::Matrix{Float64}
    olddelta::Matrix{Float64}
end

"""Allocate and initialise `WorkArrays` for `n_components` and a given block size."""
function create_work_arrays(n_components::Int, block_size::Int)
    weights = Matrix{Float64}(I, n_components, n_components)  # Initialize as identity matrix
    return WorkArrays(
        weights,  # weights - start as identity matrix
        zeros(n_components, block_size),  # u
        zeros(n_components, block_size),  # y
        zeros(n_components, block_size),  # data_block
        copy(weights),  # oldweights - copy of identity matrix
        copy(weights),  # startweights - copy of identity matrix
        zeros(n_components, n_components),  # weights_temp
        zeros(n_components, n_components),  # wu_term
        zeros(n_components, n_components),  # delta
        zeros(n_components, n_components),   # olddelta
    )
end


"""
    _run_ica_algorithm(dat_ica, layout; n_components, algorithm, params)

Internal dispatcher function that routes to the appropriate ICA algorithm implementation.

# Arguments
- `dat_ica::Matrix{Float64}`: Data matrix (channels × samples)
- `layout::Layout`: Layout information
- `n_components::Int`: Number of ICA components
- `algorithm::Symbol`: Algorithm to use (`:infomax`, `:infomax_extended`)
- `params::IcaPrms`: Algorithm-specific parameters

# Returns
`InfoIca` result structure
"""
function _run_ica_algorithm(
    dat_ica::Matrix{Float64},
    layout::Layout,
    filename::String;
    n_components::Int,
    algorithm::Symbol = :infomax,
    params::IcaPrms = IcaPrms(),
)
    if algorithm == :infomax
        return infomax_ica(dat_ica, layout, filename; n_components = n_components, params = params)
    elseif algorithm == :infomax_extended
        return infomax_extended_ica(dat_ica, layout, filename; n_components = n_components, params = params)
    elseif algorithm == :picard || algorithm == :picard_extended
        extended = algorithm == :picard_extended
        return picard_ica(dat_ica, layout, filename; n_components = n_components, extended = extended, params = params)
    elseif algorithm == :amica
        return amica_ica(dat_ica, layout, filename; n_components = n_components, params = params)
    else
        error("Unknown ICA algorithm: $algorithm. Supported algorithms: :infomax, :infomax_extended, :picard, :picard_extended, :amica")
    end
end


# =========================
# AMICA ICA IMPLEMENTATION
# =========================

# Define learning rate struct to avoid Amica.jl dependency
mutable struct AmicaLearningRate{T<:Real}
    lrate::T
    lrate0::T
    min_lrate::T
    max_lrate::T
    lratefact::T
    shapelrate::T
    shapelrate0::T
    shapelratefact::T
    numdecs::Int
    maxdecs::Int
    newtrate::T
    newt_ramp::Int
    minrho::T
    maxrho::T
end

function AmicaLearningRate(
    T::Type{<:Real} = Float64;
    lrate0 = 0.1,
    min_lrate = 1e-8,
    max_lrate = 1.0,
    lratefact = 0.5,
    shapelrate0 = 0.05,
    shapelratefact = 0.5,
    maxdecs = 3,
    newtrate = 1.0,
    newt_ramp = 10,
    minrho = 1.0,
    maxrho = 2.0,
)
    return AmicaLearningRate{T}(
        T(lrate0),
        T(lrate0),
        T(min_lrate),
        T(max_lrate),
        T(lratefact),
        T(shapelrate0),
        T(shapelrate0),
        T(shapelratefact),
        0,
        maxdecs,
        T(newtrate),
        newt_ramp,
        T(minrho),
        T(maxrho),
    )
end

# Native Amica Work Arrays (Zero-Allocation)

mutable struct AmicaWorkArrays{T<:Real}
    g_times_sources::Array{T,2}
    sum_z::Array{T,2}
    kp::Array{T,2}
    dmu_numer::Array{T,2}
    dmu_denom::Array{T,2}
    dbeta_denom::Array{T,2}
    dlambda_numer::Array{T,2}
    drho_numer::Array{T,2}
    newton_sigma2::Array{T,1}
    Lt_accum::Array{T,1}

    dA::Array{T,2}
    Lt::Array{T,1}
    newton_kappa::Array{T,1}
    newton_lambda::Array{T,1}

    # Preallocated block arrays
    source_signals::Array{T,2}
    y::Array{T,3}
    y_rho::Array{T,3}
    Q::Array{T,3}
    z::Array{T,3}
    fp::Array{T,3}
    scratch3::Array{T,3}
    scratch_N_n::Array{T,2}
    scratch_N_n2::Array{T,2}
end

function AmicaWorkArrays(T::Type{<:Real}, N::Int, n::Int, m::Int, block_size::Int)
    return AmicaWorkArrays{T}(
        zeros(T, n, n),
        zeros(T, n, m),
        zeros(T, n, m),
        zeros(T, n, m),
        zeros(T, n, m),
        zeros(T, n, m),
        zeros(T, n, m),
        zeros(T, n, m),
        zeros(T, n),
        zeros(T, N),
        zeros(T, n, n),
        zeros(T, N),
        zeros(T, n),
        zeros(T, n),
        zeros(T, block_size, n),
        zeros(T, block_size, n, m),
        zeros(T, block_size, n, m),
        zeros(T, block_size, n, m),
        zeros(T, block_size, n, m),
        zeros(T, block_size, n, m),
        zeros(T, block_size, n, m),
        zeros(T, block_size, n),
        zeros(T, block_size, n),
    )
end

function reset_accumulators!(wa::AmicaWorkArrays{T}) where {T}
    fill!(wa.g_times_sources, zero(T))
    fill!(wa.sum_z, zero(T))
    fill!(wa.kp, zero(T))
    fill!(wa.dmu_numer, zero(T))
    fill!(wa.dmu_denom, zero(T))
    fill!(wa.dbeta_denom, zero(T))
    fill!(wa.dlambda_numer, zero(T))
    fill!(wa.drho_numer, zero(T))
    fill!(wa.newton_sigma2, zero(T))
    fill!(wa.Lt_accum, zero(T))
end

function process_amica_blocks!(
    wa_threads::Vector{AmicaWorkArrays{T}},
    data::AbstractMatrix{T},
    W::AbstractMatrix{T},
    scale::AbstractMatrix{T},
    location::AbstractMatrix{T},
    shape::AbstractMatrix{T},
    proportions::AbstractMatrix{T},
    newton_active::Bool,
    block_size::Int,
) where {T<:Real}
    N, n = size(data)
    m = size(scale, 2)
    num_blocks = cld(N, block_size)
    n_threads = length(wa_threads)

    Threads.@threads for tid = 1:n_threads
        wa = wa_threads[tid]
        for block = tid:n_threads:num_blocks
            lower = (block - 1) * block_size + 1
            upper = min(N, block * block_size)
            n_samples = upper - lower + 1
            rng = lower:upper

            source_signals = view(wa.source_signals, 1:n_samples, :)
            y = view(wa.y, 1:n_samples, :, :)
            y_rho = view(wa.y_rho, 1:n_samples, :, :)
            Q = view(wa.Q, 1:n_samples, :, :)
            z = view(wa.z, 1:n_samples, :, :)
            fp = view(wa.fp, 1:n_samples, :, :)
            scratch3 = view(wa.scratch3, 1:n_samples, :, :)
            Qmax = view(wa.scratch_N_n, 1:n_samples, :)
            logexp = view(wa.scratch_N_n2, 1:n_samples, :)
            zsum = view(wa.scratch_N_n, 1:n_samples, :)
            g_block_sum = view(wa.scratch_N_n, 1:n_samples, :)

            mul!(source_signals, view(data, rng, :), transpose(W))

            for j = 1:m
                for i = 1:n
                    @fastmath @inbounds @simd for s = 1:n_samples
                        y[s, i, j] = scale[i, j] * (source_signals[s, i] - location[i, j])
                    end
                end
            end

            for j = 1:m
                for i = 1:n
                    @fastmath @inbounds @simd for s = 1:n_samples
                        val = y[s, i, j]
                        val = ifelse(val == zero(T), T(1e-16), val)
                        y_rho[s, i, j] = exp((shape[i, j] - T(1.0)) * log(abs(val)))
                    end
                end
            end

            for j = 1:m
                for i = 1:n
                    qconst = -log(T(2)) - SpecialFunctions.loggamma(T(1) + T(1) / shape[i, j]) + log(proportions[i, j]) + log(scale[i, j])
                    @fastmath @inbounds @simd for s = 1:n_samples
                        Q[s, i, j] = qconst - (y_rho[s, i, j] * abs(y[s, i, j]))
                    end
                end
            end

            maximum!(view(Qmax, :, :, 1:1), Q)
            @fastmath @inbounds @simd for I in eachindex(Q)
                scratch3[I] = exp(Q[I] - Qmax[I[1], I[2]])
            end
            sum!(view(logexp, :, :, 1:1), scratch3)
            @fastmath @inbounds @simd for I in eachindex(logexp)
                logexp[I] = log(logexp[I]) + Qmax[I]
            end

            @fastmath @inbounds @simd for I in eachindex(Q)
                z[I] = exp(Q[I] - logexp[I[1], I[2]]) + T(1e-15)
            end

            for i = 1:n
                @fastmath @inbounds @simd for s = 1:n_samples
                    wa.Lt_accum[lower+s-1] += logexp[s, i]
                end
            end

            sum!(view(zsum, :, :, 1:1), z)
            @fastmath @inbounds @simd for I in eachindex(z)
                z[I] = z[I] / zsum[I[1], I[2]]
            end

            for j = 1:m
                for i = 1:n
                    @fastmath @inbounds @simd for s = 1:n_samples
                        fp[s, i, j] = y_rho[s, i, j] * sign(y[s, i, j]) * shape[i, j]
                    end
                end
            end

            for j = 1:m
                for i = 1:n
                    @fastmath @inbounds @simd for s = 1:n_samples
                        sc = y_rho[s, i, j] * abs(y[s, i, j])
                        sc = ifelse(sc >= T(1.0e-16), z[s, i, j] * log(sc) * sc, T(0.0))
                        wa.drho_numer[i, j] += sc
                    end
                end
            end

            for j = 1:m
                for i = 1:n
                    @fastmath @inbounds @simd for s = 1:n_samples
                        scratch3[s, i, j] = scale[i, j] * z[s, i, j] * fp[s, i, j]
                    end
                end
            end
            sum!(view(g_block_sum, :, :, 1:1), scratch3)
            mul!(wa.g_times_sources, transpose(g_block_sum), source_signals, one(T), one(T))

            if newton_active
                for i = 1:n
                    @fastmath @inbounds @simd for s = 1:n_samples
                        wa.newton_sigma2[i] += (source_signals[s, i]^2) / N
                    end
                end
            end

            for j = 1:m
                for i = 1:n
                    @fastmath @inbounds @simd for s = 1:n_samples
                        z_val = z[s, i, j]
                        fp_val = fp[s, i, j]
                        y_val = y[s, i, j]
                        wa.sum_z[i, j] += z_val
                        wa.dmu_numer[i, j] += fp_val * z_val
                        wa.kp[i, j] += fp_val * z_val * fp_val
                        if shape[i, j] <= T(2)
                            y_safe = ifelse(y_val == zero(T), T(1e-16), y_val)
                            wa.dmu_denom[i, j] += (z_val * fp_val / y_safe) * scale[i, j]
                            wa.dbeta_denom[i, j] += fp_val * z_val * y_val
                        else
                            wa.dmu_denom[i, j] += (z_val * fp_val * fp_val) * scale[i, j]
                        end
                        wa.dlambda_numer[i, j] += z_val * (fp_val * y_val - T(1.0))^2
                    end
                end
            end
        end
    end
end

function update_amica_parameters!(
    wa_threads::Vector{AmicaWorkArrays{T}},
    data::AbstractMatrix{T},
    A::AbstractMatrix{T},
    LLdetS::T,
    scale::AbstractMatrix{T},
    location::AbstractMatrix{T},
    shape::AbstractMatrix{T},
    proportions::AbstractMatrix{T},
    newton_active::Bool,
    m::Int,
    upd_shape::Bool,
    shapelrate::T,
    minrho::T,
    maxrho::T,
    block_size::Int,
) where {T<:Real}
    N, n = size(data)
    ldet = -logabsdet(A)[1]

    for wa in wa_threads
        reset_accumulators!(wa)
    end

    W = inv(A)
    process_amica_blocks!(wa_threads, data, W, scale, location, shape, proportions, newton_active, block_size)

    # Reduce thread results into wa_threads[1]
    wa = wa_threads[1]
    for tid = 2:length(wa_threads)
        wa_t = wa_threads[tid]
        wa.Lt_accum .+= wa_t.Lt_accum
        wa.drho_numer .+= wa_t.drho_numer
        wa.g_times_sources .+= wa_t.g_times_sources
        wa.newton_sigma2 .+= wa_t.newton_sigma2
        wa.sum_z .+= wa_t.sum_z
        wa.dmu_numer .+= wa_t.dmu_numer
        wa.kp .+= wa_t.kp
        wa.dmu_denom .+= wa_t.dmu_denom
        wa.dbeta_denom .+= wa_t.dbeta_denom
        wa.dlambda_numer .+= wa_t.dlambda_numer
    end

    wa.Lt .= ldet .+ LLdetS
    wa.Lt .+= wa.Lt_accum
    wa.dA .= I(n) - wa.g_times_sources ./ N
    LL_iter = sum(wa.Lt) / (N * n)

    if m > 1
        proportions .= ifelse.(wa.sum_z .>= T(0), wa.sum_z ./ N, T(1) / N)
    end

    if newton_active
        dkap = (wa.kp ./ (proportions .* N)) .* scale .^ 2
        wa.newton_kappa .= dropdims(sum(proportions .* dkap, dims = 2), dims = 2)
        wa.newton_lambda .= dropdims(sum(proportions .* (wa.dlambda_numer ./ wa.sum_z .+ dkap .* location .^ 2), dims = 2), dims = 2)
    end

    if m > 1
        location .+= wa.dmu_numer ./ wa.dmu_denom
    end

    scale .*= sqrt.(wa.sum_z ./ wa.dbeta_denom)

    if upd_shape
        shape .=
            clamp.(
                shape .+ (shapelrate .* (1 .- (shape ./ SpecialFunctions.digamma.(1 .+ 1 ./ shape)) .* wa.drho_numer ./ wa.sum_z)),
                minrho,
                maxrho,
            )
    end

    return LL_iter
end

function update_amica_mixing!(
    A::AbstractMatrix{T},
    dA::AbstractMatrix{T},
    newton_kappa::AbstractVector{T},
    newton_lambda::AbstractVector{T},
    newton_sigma2::AbstractVector{T},
    iter::Int,
    do_newton::Bool,
    no_newton::Bool,
    newt_start_iter::Int,
    lrate::AmicaLearningRate{T},
) where {T<:Real}
    n = size(A, 1)

    if (do_newton && !no_newton && iter >= newt_start_iter)
        if iter == newt_start_iter
            lrate.numdecs = 0
        end

        B = similar(dA)
        posdef = true
        for i = 1:n
            for k = 1:n
                if i == k
                    if isfinite(newton_lambda[i]) && abs(newton_lambda[i]) > eps(T)
                        B[i, k] = dA[i, k] / newton_lambda[i]
                    else
                        B[i, k] = zero(T)
                        posdef = false
                    end
                else
                    sk1 = newton_sigma2[i] * newton_kappa[k]
                    sk2 = newton_sigma2[k] * newton_kappa[i]
                    denom = sk1 * sk2 - one(T)
                    if isfinite(sk1) && isfinite(sk2) && isfinite(denom) && denom > eps(T)
                        B[i, k] = (sk1 * dA[i, k] - dA[k, i]) / denom
                    else
                        B[i, k] = zero(T)
                        posdef = false
                    end
                end
            end
        end

        if posdef
            lrate.lrate = min(lrate.newtrate, lrate.lrate + min(T(1.0) / lrate.newt_ramp, lrate.lrate))
            A .-= lrate.lrate .* A * B
            return no_newton
        else
            no_newton = true
            lrate.lrate = min(lrate.lrate0, lrate.lrate + min(T(1.0) / lrate.newt_ramp, lrate.lrate))
            A .-= lrate.lrate .* A * dA
            return no_newton
        end
    else
        lrate.lrate = min(lrate.lrate0, lrate.lrate + min(T(1.0) / lrate.newt_ramp, lrate.lrate))
        A .-= lrate.lrate .* A * dA
        return no_newton
    end
end

function reparameterize_amica!(A::AbstractMatrix{T}, location::AbstractMatrix{T}, scale::AbstractMatrix{T}) where {T<:Real}
    tau = dropdims(sqrt.(sum(A .^ 2, dims = 1)), dims = 1)
    mask = tau .> zero(T)
    A .= ifelse.(reshape(mask, 1, size(mask)...), A ./ tau', A)
    location .= ifelse.(mask, location .* tau, location)
    scale .= ifelse.(mask, scale ./ tau, scale)
end

function calculate_amica_lrate!(
    LL::AbstractVector,
    iter::Int,
    newt_start_iter::Int,
    do_newton::Bool,
    lrate::AmicaLearningRate{T},
) where {T<:Real}
    if LL[iter] < LL[iter-1]
        if lrate.lrate <= lrate.min_lrate
            return true
        else
            lrate.lrate *= lrate.lratefact
            lrate.shapelrate *= lrate.shapelratefact
            lrate.numdecs += 1
            if lrate.numdecs >= lrate.maxdecs
                lrate.lrate0 *= lrate.lratefact
                if iter > newt_start_iter
                    lrate.shapelrate0 *= lrate.shapelratefact
                end
                if do_newton && (iter > newt_start_iter)
                    lrate.newtrate *= lrate.lratefact
                end
                lrate.numdecs = 0
            end
        end
    end
    return false
end

"""AMICA ICA implementation natively integrated for EegFun.jl"""
function amica_ica(dat_ica::Matrix{Float64}, layout::Layout, filename::String; n_components::Int, params::IcaPrms = IcaPrms())
    n_channels, n_samples = size(dat_ica)

    # PCA / Sphering identically to infomax
    original_mean = vec(mean(dat_ica, dims = 2))
    dat_ica .-= original_mean

    data_scale = sqrt(norm((dat_ica * transpose(dat_ica)) / n_samples))
    dat_ica ./= data_scale

    # Apply SVD
    F = svd(dat_ica)
    pca_components = F.U[:, 1:n_components]
    eigenvalues = (F.S[1:n_components] .^ 2) ./ max(n_samples - 1, 1)
    sphere = diagm(1.0 ./ sqrt.(eigenvalues))

    transform_matrix = sphere * transpose(pca_components)
    dat_ica_sphered = Matrix{Float64}(undef, n_components, n_samples)
    mul!(dat_ica_sphered, transform_matrix, dat_ica)
    dat_ica = dat_ica_sphered

    # AMICA native math expects data as (N_samples, n_components)
    dat_ica_t = Matrix{Float64}(dat_ica')

    LLdetS = logabsdet(sphere)[1]

    # Amica parameters
    m = 3 # number of mixtures, could be parameterized
    maxiter = params.max_iter > 0 ? params.max_iter : 2000
    lrate_val = params.l_rate > 0 ? params.l_rate : 0.1
    lrate = AmicaLearningRate(Float64; lrate0 = lrate_val)

    # Dynamic block size based on N
    cpu_block = min(nextpow(2, max(Int(floor(sqrt(n_samples / 3.0))), 1)), 512)
    block_size = min(n_samples, max(cpu_block, 32))

    rng = Xoshiro(42)

    # Initialize A to match Fortran/Amica.jl: small random ±0.005, diagonal = 1.0, then normalize columns
    Wtmp = rand(rng, Float64, n_components, n_components)
    A = 0.01 .* (0.5 .- Wtmp)
    for i = 1:n_components
        A[i, i] = 1.0
        A[:, i] = A[:, i] / norm(A[:, i])
    end

    proportions = fill(1.0 / m, n_components, m)
    shape = fill(1.5, n_components, m)

    # Initialize location to match Fortran: centered around 0 (e.g., -1, 0, 1 for m=3) + small noise ±0.05
    location = zeros(Float64, n_components, m)
    for j = 1:m
        location[:, j] .= Float64(j - 1 - (m - 1) / 2)
    end
    location .+= 0.05 .* (1.0 .- 2.0 .* rand(rng, Float64, n_components, m))

    # Initialize scale to match Fortran: values in range [0.95, 1.05]
    scale = ones(Float64, n_components, m) .+ 0.1 .* (0.5 .- rand(rng, Float64, n_components, m))

    n_threads = Threads.nthreads()
    wa_threads = [AmicaWorkArrays(Float64, n_samples, n_components, m, block_size) for _ = 1:n_threads]
    LL = fill(-Inf, maxiter)

    do_newton = true
    no_newton = false
    newt_start_iter = 50
    update_shape = true

    # Main loop
    for iter = 1:maxiter
        LL_iter = update_amica_parameters!(
            wa_threads,
            dat_ica_t,
            A,
            LLdetS,
            scale,
            location,
            shape,
            proportions,
            do_newton && iter >= newt_start_iter,
            m,
            update_shape,
            lrate.shapelrate,
            lrate.minrho,
            lrate.maxrho,
            block_size,
        )

        LL[iter] = LL_iter

        wa = wa_threads[1]
        no_newton = update_amica_mixing!(
            A,
            wa.dA,
            wa.newton_kappa,
            wa.newton_lambda,
            wa.newton_sigma2,
            iter,
            do_newton,
            no_newton,
            newt_start_iter,
            lrate,
        )

        reparameterize_amica!(A, location, scale)

        if iter > 1
            if isnan(LL[iter])
                @warn("Got NaN! Exiting ...")
                break
            end
            if calculate_amica_lrate!(LL, iter, newt_start_iter, do_newton, lrate)
                break
            end

            # Convergence check (from Amica)
            if abs(LL[iter] - LL[iter-1]) < 1e-7
                break
            end
        end

        if iter == 1 || iter % 10 == 0
            @info "AMICA iter $iter, LL = $(LL[iter])"
        end
    end

    # Final calculations
    W = inv(A)
    weights = W * sphere * pca_components'
    mixing = pinv(weights)

    # Calculate total variance explained and order
    # Since dat_ica was overwritten with the sphered data, its row variances are uniform.
    # The true variance of source i is proportional to the squared norm of the i-th row of W.
    source_vars = vec(sum(abs2, W, dims = 2))
    meanvar = vec(sum(abs2, mixing, dims = 1)) .* source_vars
    meanvar_normalized = meanvar ./ sum(meanvar)
    order = sortperm(meanvar_normalized, rev = true)

    return InfoIca(
        filename,
        weights[order, :],
        mixing[:, order],
        sphere,
        meanvar_normalized[order],
        data_scale,
        original_mean,
        [Symbol("IC$i") for i = 1:size(weights, 1)],
        Dict{Int,Matrix{Float64}}(),
        layout,
        falses(size(weights, 1)),
    )
end


# =============================================================================
# INFOMAX ICA IMPLEMENTATION
# =============================================================================

"""Standard Infomax ICA implementation (super-Gaussian sources only)."""
function infomax_ica(dat_ica::Matrix{Float64}, layout::Layout, filename::String; n_components::Int, params::IcaPrms = IcaPrms())

    # Store original mean before removing it
    original_mean = vec(mean(dat_ica, dims = 2))

    n_channels, n_samples = size(dat_ica)

    # Center and scale data
    dat_ica .-= original_mean
    scale = sqrt(norm((dat_ica * transpose(dat_ica)) / n_samples))
    dat_ica ./= scale

    # PCA reduction using high-precision SVD
    # This prevents smallest components from degrading into numerical noise
    F = svd(dat_ica)
    pca_components = F.U[:, 1:n_components]

    # Analytic sphering using exact singular values (skips a 260 MiB workspace allocation)
    eigenvalues = (F.S[1:n_components] .^ 2) ./ (n_samples - 1)
    sphere = diagm(1.0 ./ sqrt.(eigenvalues))

    # Apply transform and create final workspace matrix (n_components × n_samples)
    transform_matrix = sphere * transpose(pca_components)
    dat_ica_sphered = Matrix{Float64}(undef, n_components, n_samples)
    mul!(dat_ica_sphered, transform_matrix, dat_ica)
    dat_ica = dat_ica_sphered

    # initialize
    n_channels = size(dat_ica, 1)
    n_samples = size(dat_ica, 2)

    # Detect GPU hardware backend
    gpu_active = false
    if params.use_gpu
        if is_gpu_available()
            gpu_active = true
            @info "[GPU ACTIVATED] Running Infomax ICA on $(gpu_device_name())..."

            # Dynamically calculate GPU block size
            cpu_block = min(nextpow(2, max(Int(floor(sqrt(n_samples / 3.0))), 1)), 512)
            gpu_block = min(nextpow(2, cpu_block * 4), 1024)
            block = min(n_samples, max(gpu_block, 32))

            work = create_work_arrays(n_channels, block)
        else
            @minimal_warning "Requested GPU acceleration (use_gpu=true), but no functional GPU package (CUDA.jl, AMDGPU.jl, Metal.jl) has been loaded. Please run 'using CUDA' or 'using AMDGPU' before calling run_ica. Falling back to CPU."
            raw_block = Int(floor(sqrt(n_samples / 3.0)))
            block = min(nextpow(2, max(raw_block, 1)), 512)
            work = create_work_arrays(n_channels, block)
        end
    else
        raw_block = Int(floor(sqrt(n_samples / 3.0)))
        block = min(nextpow(2, max(raw_block, 1)), 512)
        work = create_work_arrays(n_channels, block)
    end

    if gpu_active
        _infomax_optimize_gpu!(work, dat_ica, params, block)
    else
        _infomax_optimize_cpu!(work, dat_ica, params, block)
    end

    # Final calculations
    work.weights = work.weights * sphere * pca_components'
    mixing = pinv(work.weights)

    # calculate total variance explained and order
    meanvar = vec(sum(abs2, mixing, dims = 1) .* sum(abs2, dat_ica, dims = 2)' ./ (n_components * n_samples - 1))
    meanvar_normalized = meanvar ./ sum(meanvar)
    order = sortperm(meanvar_normalized, rev = true)

    return InfoIca(
        filename,
        work.weights[order, :],
        mixing[:, order],
        sphere,
        meanvar_normalized[order],
        scale,
        original_mean,
        [Symbol("IC$i") for i = 1:size(work.weights, 1)],
        Dict{Int,Matrix{Float64}}(),
        layout,
        falses(size(work.weights, 1)),  # Regular Infomax: all super-Gaussian (all false)
    )

end

function _infomax_optimize_cpu!(work::WorkArrays, dat_ica::Matrix{Float64}, params::IcaPrms, block::Int)
    n_channels, n_samples = size(dat_ica)
    n_components = n_channels

    step = 0
    wts_blowup = false
    change = 0.0
    oldchange = 0.0
    angledelta = 0.0

    permute_indices = Vector{Int}(undef, n_samples)

    @inbounds while step < params.max_iter
        randperm!(permute_indices)

        for t = 1:block:n_samples
            block_end = min(t + block - 1, n_samples)
            block_size = block_end - t + 1

            # extract data block
            @inbounds for j = 1:block_size
                idx = permute_indices[t+j-1]
                @simd for i = 1:n_channels
                    work.data_block[i, j] = dat_ica[i, idx]
                end
            end

            # forward pass & weight update
            if block_size < block
                mul!(view(work.u, :, 1:block_size), work.weights, view(work.data_block, :, 1:block_size))
            else
                mul!(work.u, work.weights, work.data_block)
            end

            @fastmath @inbounds for i = 1:n_components
                @simd for j = 1:block_size
                    work.y[i, j] = 1.0 - 2.0 / (1.0 + exp(-work.u[i, j]))
                end
            end

            # update weights 
            if block_size < block
                mul!(work.wu_term, view(work.y, :, 1:block_size), transpose(view(work.u, :, 1:block_size)))
                mul!(work.weights_temp, work.wu_term, work.weights)
                @. work.weights += params.l_rate * (block_size * work.weights + work.weights_temp)
            else
                mul!(work.wu_term, work.y, transpose(work.u))
                mul!(work.weights_temp, work.wu_term, work.weights)
                @. work.weights += params.l_rate * (block_size * work.weights + work.weights_temp)
            end

            # boom?
            if maximum(abs, work.weights) > params.max_weight
                wts_blowup = true
                change = NaN
                break
            end
        end

        if !wts_blowup
            work.oldweights .-= work.weights
            step += 1
            work.delta .= work.oldweights
            change = dot(work.delta, work.delta)
        end

        if wts_blowup || isnan(change) || isinf(change)
            step = 0
            change = NaN
            wts_blowup = false
            params.l_rate *= params.restart_factor
            work.weights .= work.startweights
            work.oldweights .= work.startweights
            continue
        end

        if step > 2
            angledelta = acos(clamp(dot(work.delta, work.olddelta) / sqrt(change * oldchange), -1, 1))
            if params.degconst * angledelta > params.anneal_deg
                params.l_rate *= params.anneal_step
                work.olddelta .= work.delta
                oldchange = change
            end
            change < params.w_change && break
        elseif step == 1
            work.olddelta .= work.delta
            oldchange = change
        end

        work.oldweights .= work.weights
        change > params.blowup && (params.l_rate *= params.blowup_fac)

        if step == 1 || step % 10 == 0
            @info Printf.@sprintf(
                "Infomax step %d, change = %.7f, lrate = %.7f, angle = %.1f",
                step,
                change,
                params.l_rate,
                (params.degconst) * angledelta
            )
        end
    end
end

function _infomax_optimize_gpu!(work::WorkArrays, dat_ica::Matrix{Float64}, params::IcaPrms, block::Int)
    n_channels, n_samples = size(dat_ica)
    n_components = n_channels

    step = 0
    wts_blowup = false
    change = 0.0
    oldchange = 0.0
    angledelta = 0.0

    permute_indices = Vector{Int}(undef, n_samples)

    dat_ica_gpu = gpu_array(Float32.(dat_ica))
    weights_gpu = gpu_array(Float32.(work.weights))
    u_gpu = gpu_array(zeros(Float32, n_components, block))
    y_gpu = gpu_array(zeros(Float32, n_components, block))
    wu_gpu = gpu_array(zeros(Float32, n_components, n_components))
    wtemp_gpu = gpu_array(zeros(Float32, n_components, n_components))

    @inbounds while step < params.max_iter
        randperm!(permute_indices)

        perm_gpu = gpu_array(permute_indices)
        dat_ica_gpu_shuffled = dat_ica_gpu[:, perm_gpu]

        u_view_full = view(u_gpu, :, 1:block)
        y_view_full = view(y_gpu, :, 1:block)
        u_view_full_t = transpose(u_view_full)

        for t = 1:block:n_samples
            block_end = min(t + block - 1, n_samples)
            block_size = block_end - t + 1

            try
                data_block_gpu = view(dat_ica_gpu_shuffled, :, t:block_end)

                if block_size == block
                    mul!(u_view_full, weights_gpu, data_block_gpu)
                    @. y_view_full = 1.0f0 - 2.0f0 / (1.0f0 + exp(-u_view_full))
                    mul!(wu_gpu, y_view_full, u_view_full_t)
                else
                    u_view = view(u_gpu, :, 1:block_size)
                    y_view = view(y_gpu, :, 1:block_size)
                    mul!(u_view, weights_gpu, data_block_gpu)
                    @. y_view = 1.0f0 - 2.0f0 / (1.0f0 + exp(-u_view))
                    mul!(wu_gpu, y_view, transpose(u_view))
                end

                mul!(wtemp_gpu, wu_gpu, weights_gpu)
                @. weights_gpu += Float32(params.l_rate) * (Float32(block_size) * weights_gpu + wtemp_gpu)
            catch e
                @error "GPU processing failed in inner loop, throwing error." exception = (e, catch_backtrace())
                rethrow(e)
            end

            # Check boom on GPU
            if maximum(abs, weights_gpu) > params.max_weight
                wts_blowup = true
                change = NaN
                break
            end
        end

        if !wts_blowup
            work.weights .= Array(weights_gpu)
            work.oldweights .-= work.weights
            step += 1
            work.delta .= work.oldweights
            change = dot(work.delta, work.delta)
        end

        if wts_blowup || isnan(change) || isinf(change)
            step = 0
            change = NaN
            wts_blowup = false
            params.l_rate *= params.restart_factor
            work.weights .= work.startweights
            work.oldweights .= work.startweights
            weights_gpu .= gpu_array(Float32.(work.weights))
            continue
        end

        if step > 2
            angledelta = acos(clamp(dot(work.delta, work.olddelta) / sqrt(change * oldchange), -1, 1))
            if params.degconst * angledelta > params.anneal_deg
                params.l_rate *= params.anneal_step
                work.olddelta .= work.delta
                oldchange = change
            end
            change < params.w_change && break
        elseif step == 1
            work.olddelta .= work.delta
            oldchange = change
        end

        work.oldweights .= work.weights
        change > params.blowup && (params.l_rate *= params.blowup_fac)

        if step == 1 || step % 10 == 0
            @info Printf.@sprintf(
                "Infomax step %d, change = %.7f, lrate = %.7f, angle = %.1f",
                step,
                change,
                params.l_rate,
                (params.degconst) * angledelta
            )
        end
    end

    work.weights .= Array(weights_gpu)
end

# =============================================================================
# EXTENDED INFOMAX ICA IMPLEMENTATION
# =============================================================================

"""
    infomax_extended_ica(dat_ica::Matrix{Float64}, layout::Layout; n_components::Int, params::IcaPrms = IcaPrms())

Extended Infomax algorithm implementation that can separate both sub-Gaussian and super-Gaussian sources.

Extended Infomax is an enhancement of the standard Infomax algorithm that adapts to the statistical
properties of source signals. It uses kurtosis-based switching to handle both sub-Gaussian (kurtosis < 0)
and super-Gaussian (kurtosis > 0) sources, making it more versatile than standard Infomax.

# Arguments
- `dat_ica::Matrix{Float64}`: Data matrix (channels × samples), should be preprocessed
- `layout::Layout`: Layout information for channels
- `n_components::Int`: Number of ICA components to extract
- `params::IcaPrms`: ICA parameters (uses max_iter, w_change from params)

# Returns
`InfoIca` with unmixing, mixing, sphere, variance, and metadata.

# References
Lee, T. W., Girolami, M., & Sejnowski, T. J. (1999). Independent component analysis using an extended 
infomax algorithm for mixed subgaussian and supergaussian sources. Neural computation, 11(2), 417-441.
"""
function infomax_extended_ica(dat_ica::Matrix{Float64}, layout::Layout, filename::String; n_components::Int, params::IcaPrms = IcaPrms())

    # Store original mean before removing it
    original_mean = vec(mean(dat_ica, dims = 2))

    n_channels, n_samples = size(dat_ica)

    # Center and scale data
    dat_ica .-= original_mean
    scale = sqrt(norm((dat_ica * transpose(dat_ica)) / n_samples))
    dat_ica ./= scale

    # PCA reduction using high-precision SVD
    # This prevents smallest components from degrading into numerical noise
    F = svd(dat_ica)
    pca_components = F.U[:, 1:n_components]

    # Analytic sphering using exact singular values (skips a 260 MiB workspace allocation)
    eigenvalues = (F.S[1:n_components] .^ 2) ./ (n_samples - 1)
    sphere = diagm(1.0 ./ sqrt.(eigenvalues))

    # Apply transform and create final workspace matrix (n_components × n_samples)
    transform_matrix = sphere * transpose(pca_components)
    dat_ica_sphered = Matrix{Float64}(undef, n_components, n_samples)
    mul!(dat_ica_sphered, transform_matrix, dat_ica)
    dat_ica = dat_ica_sphered

    # initialize
    n_channels = size(dat_ica, 1)
    n_samples = size(dat_ica, 2)

    # Track kurtosis signs for each component (true = sub-Gaussian, false = super-Gaussian)
    is_sub_gaussian = falses(n_channels)

    # Detect GPU hardware backend
    gpu_active = false
    if params.use_gpu
        if is_gpu_available()
            gpu_active = true
            @info "[GPU ACTIVATED] Running Extended Infomax ICA on $(gpu_device_name())..."

            # Dynamically calculate GPU block size
            cpu_block = min(nextpow(2, max(Int(floor(sqrt(n_samples / 3.0))), 1)), 512)
            gpu_block = min(nextpow(2, cpu_block * 4), 1024)
            block = min(n_samples, max(gpu_block, 32))

            work = create_work_arrays(n_channels, block)
        else
            @minimal_warning "Requested GPU acceleration (use_gpu=true), but no functional GPU package (CUDA.jl, AMDGPU.jl, Metal.jl) has been loaded. Please run 'using CUDA' or 'using AMDGPU' before calling run_ica. Falling back to CPU."
            raw_block = Int(floor(sqrt(n_samples / 3.0)))
            block = min(nextpow(2, max(raw_block, 1)), 512)
            work = create_work_arrays(n_channels, block)
        end
    else
        raw_block = Int(floor(sqrt(n_samples / 3.0)))
        block = min(nextpow(2, max(raw_block, 1)), 512)
        work = create_work_arrays(n_channels, block)
    end

    if gpu_active
        _infomax_extended_optimize_gpu!(work, dat_ica, params, block, is_sub_gaussian)
    else
        _infomax_extended_optimize_cpu!(work, dat_ica, params, block, is_sub_gaussian)
    end

    # Final calculations
    work.weights = work.weights * sphere * pca_components'
    mixing = pinv(work.weights)

    # Log final kurtosis distribution
    n_sub_final = count(is_sub_gaussian)
    n_super_final = n_channels - n_sub_final
    @info "Extended Infomax completed: sup/sub-gauss ($n_super_final/$n_sub_final)"

    # calculate total variance explained and order
    meanvar = vec(sum(abs2, mixing, dims = 1) .* sum(abs2, dat_ica, dims = 2)' ./ (n_components * n_samples - 1))
    meanvar_normalized = meanvar ./ sum(meanvar)
    order = sortperm(meanvar_normalized, rev = true)

    return InfoIca(
        filename,
        work.weights[order, :],
        mixing[:, order],
        sphere,
        meanvar_normalized[order],
        scale,
        original_mean,
        [Symbol("IC$i") for i = 1:size(work.weights, 1)],
        Dict{Int,Matrix{Float64}}(),
        layout,
        is_sub_gaussian[order],
    )
end

function _infomax_extended_optimize_cpu!(
    work::WorkArrays,
    dat_ica::Matrix{Float64},
    params::IcaPrms,
    block::Int,
    is_sub_gaussian::BitVector,
)
    n_channels, n_samples = size(dat_ica)
    n_components = n_channels

    old_kurtosis = zeros(Float64, n_channels)
    extmomentum = 0.5

    step = 0
    wts_blowup = false
    change = 0.0
    oldchange = 0.0
    angledelta = 0.0

    permute_indices = Vector{Int}(undef, n_samples)

    n_super = n_channels
    n_sub = 0

    @inbounds while step < params.max_iter
        randperm!(permute_indices)

        for t = 1:block:n_samples
            block_end = min(t + block - 1, n_samples)
            block_size = block_end - t + 1

            # extract data block
            @inbounds for j = 1:block_size
                idx = permute_indices[t+j-1]
                @simd for i = 1:n_channels
                    work.data_block[i, j] = dat_ica[i, idx]
                end
            end

            # forward pass
            if block_size < block
                mul!(view(work.u, :, 1:block_size), work.weights, view(work.data_block, :, 1:block_size))
            else
                mul!(work.u, work.weights, work.data_block)
            end

            # Extended Infomax: use different nonlinearities based on kurtosis sign
            @fastmath @inbounds for i = 1:n_channels
                if !is_sub_gaussian[i] # Super-Gaussian
                    @simd for j = 1:block_size
                        work.y[i, j] = 1.0 - 2.0 / (1.0 + exp(-work.u[i, j]))
                    end
                else # Sub-Gaussian
                    @simd for j = 1:block_size
                        work.y[i, j] = -tanh(work.u[i, j])
                    end
                end
            end

            # update weights 
            if block_size < block
                mul!(work.wu_term, view(work.y, :, 1:block_size), transpose(view(work.u, :, 1:block_size)))
                mul!(work.weights_temp, work.wu_term, work.weights)
                @. work.weights += params.l_rate * (block_size * work.weights + work.weights_temp)
            else
                mul!(work.wu_term, work.y, transpose(work.u))
                mul!(work.weights_temp, work.wu_term, work.weights)
                @. work.weights += params.l_rate * (block_size * work.weights + work.weights_temp)
            end

            if maximum(abs, work.weights) > params.max_weight
                wts_blowup = true
                change = NaN
                break
            end
        end

        if !wts_blowup
            work.oldweights .-= work.weights
            step += 1
            work.delta .= work.oldweights
            change = dot(work.delta, work.delta)
        end

        if wts_blowup || isnan(change) || isinf(change)
            step = 0
            change = NaN
            wts_blowup = false
            params.l_rate *= params.restart_factor
            work.weights .= work.startweights
            work.oldweights .= work.startweights
            continue
        end

        # Extended Infomax: update kurtosis signs periodically
        if step > 10 && step % 10 == 0
            kurtsize = min(2000, n_samples)
            n_switched = 0
            kurtosis_values = Vector{Float64}(undef, n_channels)

            # CPU Kurtosis logic
            activations = Matrix{Float64}(undef, n_channels, kurtsize)
            if kurtsize < n_samples
                rp = randperm(n_samples)[1:kurtsize]
                @inbounds for j = 1:kurtsize
                    idx = rp[j]
                    for i = 1:n_channels
                        sum_val = 0.0
                        @simd for k = 1:n_channels
                            sum_val += work.weights[i, k] * dat_ica[k, idx]
                        end
                        activations[i, j] = sum_val
                    end
                end
            else
                mul!(activations, work.weights, dat_ica)
            end

            for i = 1:n_channels
                u2_sum = 0.0
                u4_sum = 0.0
                @inbounds @simd for j = 1:kurtsize
                    val = activations[i, j]
                    val2 = val * val
                    u2_sum += val2
                    u4_sum += val2 * val2
                end
                u2_mean_sq = (u2_sum / kurtsize)^2
                if u2_mean_sq > eps(Float64)
                    u4_mean = u4_sum / kurtsize
                    kurtosis_raw = (u4_mean / u2_mean_sq) - 3.0

                    if old_kurtosis[i] != 0.0
                        kurtosis = extmomentum * old_kurtosis[i] + (1.0 - extmomentum) * kurtosis_raw
                    else
                        kurtosis = kurtosis_raw
                    end
                    old_kurtosis[i] = kurtosis
                    kurtosis_values[i] = kurtosis

                    kurtosis_threshold = 0.05
                    old_is_sub = is_sub_gaussian[i]
                    is_sub_gaussian[i] = kurtosis < -kurtosis_threshold ? true : (kurtosis > kurtosis_threshold ? false : old_is_sub)

                    if old_is_sub != is_sub_gaussian[i]
                        n_switched += 1
                    end
                else
                    kurtosis_values[i] = 0.0
                end
            end

            # Update kurtosis counts
            n_sub = count(is_sub_gaussian)
            n_super = n_channels - n_sub
        end

        if step > 2
            angledelta = acos(clamp(dot(work.delta, work.olddelta) / sqrt(change * oldchange), -1, 1))
            if params.degconst * angledelta > params.anneal_deg
                params.l_rate *= params.anneal_step
                work.olddelta .= work.delta
                oldchange = change
            end
            change < params.w_change && break
        elseif step == 1
            work.olddelta .= work.delta
            oldchange = change
        end

        work.oldweights .= work.weights
        change > params.blowup && (params.l_rate *= params.blowup_fac)

        if step == 1 || step % 10 == 0
            @info Printf.@sprintf(
                "Extended-Infomax step %d, change = %.7f, lrate = %.7f, angle = %.1f, sup/sub-gauss: %d/%d",
                step,
                change,
                params.l_rate,
                (params.degconst) * angledelta,
                n_super,
                n_sub
            )
        end
    end
end

function _infomax_extended_optimize_gpu!(
    work::WorkArrays,
    dat_ica::Matrix{Float64},
    params::IcaPrms,
    block::Int,
    is_sub_gaussian::BitVector,
)
    n_channels, n_samples = size(dat_ica)
    n_components = n_channels

    old_kurtosis = zeros(Float64, n_channels)
    extmomentum = 0.5

    step = 0
    wts_blowup = false
    change = 0.0
    oldchange = 0.0
    angledelta = 0.0

    permute_indices = Vector{Int}(undef, n_samples)
    n_super = n_channels
    n_sub = 0

    dat_ica_gpu = gpu_array(Float32.(dat_ica))
    weights_gpu = gpu_array(Float32.(work.weights))
    u_gpu = gpu_array(zeros(Float32, n_components, block))
    y_gpu = gpu_array(zeros(Float32, n_components, block))
    wu_gpu = gpu_array(zeros(Float32, n_components, n_components))
    wtemp_gpu = gpu_array(zeros(Float32, n_components, n_components))
    is_sub_gaussian_gpu = gpu_array(is_sub_gaussian)

    kurtsize = min(2000, n_samples)
    activations_gpu = gpu_array(zeros(Float32, n_channels, kurtsize))

    @inbounds while step < params.max_iter
        randperm!(permute_indices)

        perm_gpu = gpu_array(permute_indices)
        dat_ica_gpu_shuffled = dat_ica_gpu[:, perm_gpu]

        u_view_full = view(u_gpu, :, 1:block)
        y_view_full = view(y_gpu, :, 1:block)
        u_view_full_t = transpose(u_view_full)

        for t = 1:block:n_samples
            block_end = min(t + block - 1, n_samples)
            block_size = block_end - t + 1

            try
                data_block_gpu = view(dat_ica_gpu_shuffled, :, t:block_end)

                if block_size == block
                    mul!(u_view_full, weights_gpu, data_block_gpu)
                    @. y_view_full = ifelse(is_sub_gaussian_gpu, -tanh(u_view_full), 1.0f0 - 2.0f0 / (1.0f0 + exp(-u_view_full)))
                    mul!(wu_gpu, y_view_full, u_view_full_t)
                else
                    u_view = view(u_gpu, :, 1:block_size)
                    y_view = view(y_gpu, :, 1:block_size)
                    mul!(u_view, weights_gpu, data_block_gpu)
                    @. y_view = ifelse(is_sub_gaussian_gpu, -tanh(u_view), 1.0f0 - 2.0f0 / (1.0f0 + exp(-u_view)))
                    mul!(wu_gpu, y_view, transpose(u_view))
                end

                mul!(wtemp_gpu, wu_gpu, weights_gpu)
                @. weights_gpu += Float32(params.l_rate) * (Float32(block_size) * weights_gpu + wtemp_gpu)
            catch e
                @error "GPU processing failed in inner loop, throwing error." exception = (e, catch_backtrace())
                rethrow(e)
            end

            # Check boom on GPU
            if maximum(abs, weights_gpu) > params.max_weight
                wts_blowup = true
                change = NaN
                break
            end
        end

        if !wts_blowup
            work.weights .= Array(weights_gpu)
            work.oldweights .-= work.weights
            step += 1
            work.delta .= work.oldweights
            change = dot(work.delta, work.delta)
        end

        if wts_blowup || isnan(change) || isinf(change)
            step = 0
            change = NaN
            wts_blowup = false
            params.l_rate *= params.restart_factor
            work.weights .= work.startweights
            work.oldweights .= work.startweights
            weights_gpu .= gpu_array(Float32.(work.weights))
            continue
        end

        # Extended Infomax: update kurtosis signs periodically
        if step > 10 && step % 10 == 0
            rp_cpu = randperm(n_samples)[1:kurtsize]
            rp_gpu = gpu_array(rp_cpu)
            dat_subset = dat_ica_gpu[:, rp_gpu]
            mul!(activations_gpu, weights_gpu, dat_subset)

            # Compute raw kurtosis entirely on the GPU using mapreduce patterns or broadcasting
            act_sq = activations_gpu .^ 2
            u2_sum = sum(act_sq, dims = 2)
            u2_mean_sq = (u2_sum ./ kurtsize) .^ 2

            act_quad = act_sq .^ 2
            u4_mean = sum(act_quad, dims = 2) ./ kurtsize

            # Fetch back to CPU for momentum smoothing
            u2_mean_sq_cpu = Array(u2_mean_sq)
            u4_mean_cpu = Array(u4_mean)

            n_switched = 0

            for i = 1:n_channels
                if u2_mean_sq_cpu[i] > eps(Float64)
                    kurtosis_raw = (u4_mean_cpu[i] / u2_mean_sq_cpu[i]) - 3.0

                    if old_kurtosis[i] != 0.0
                        kurtosis = extmomentum * old_kurtosis[i] + (1.0 - extmomentum) * kurtosis_raw
                    else
                        kurtosis = kurtosis_raw
                    end
                    old_kurtosis[i] = kurtosis

                    kurtosis_threshold = 0.05
                    old_is_sub = is_sub_gaussian[i]
                    is_sub_gaussian[i] = kurtosis < -kurtosis_threshold ? true : (kurtosis > kurtosis_threshold ? false : old_is_sub)

                    if old_is_sub != is_sub_gaussian[i]
                        n_switched += 1
                    end
                end
            end

            # Push updated boolean vector to GPU
            is_sub_gaussian_gpu .= gpu_array(is_sub_gaussian)

            # Update kurtosis counts
            n_sub = count(is_sub_gaussian)
            n_super = n_channels - n_sub
        end

        if step > 2
            angledelta = acos(clamp(dot(work.delta, work.olddelta) / sqrt(change * oldchange), -1, 1))
            if params.degconst * angledelta > params.anneal_deg
                params.l_rate *= params.anneal_step
                work.olddelta .= work.delta
                oldchange = change
            end
            change < params.w_change && break
        elseif step == 1
            work.olddelta .= work.delta
            oldchange = change
        end

        work.oldweights .= work.weights
        change > params.blowup && (params.l_rate *= params.blowup_fac)

        if step == 1 || step % 10 == 0
            @info Printf.@sprintf(
                "Extended-Infomax step %d, change = %.7f, lrate = %.7f, angle = %.1f, sup/sub-gauss: %d/%d",
                step,
                change,
                params.l_rate,
                (params.degconst) * angledelta,
                n_super,
                n_sub
            )
        end
    end

    work.weights .= Array(weights_gpu)
end



"""
    subtract_ica_components!(dat::DataFrame, ica::InfoIca; component_selection::Function = components())
    subtract_ica_components!(dat::ContinuousData, ica::InfoIca; component_selection::Function = components())
    subtract_ica_components(dat::DataFrame, ica::InfoIca; component_selection::Function = components())
    subtract_ica_components(dat::ContinuousData, ica::InfoIca; component_selection::Function = components())

Remove ICA components from data. Mutating (`!`) versions operate in-place; non-mutating versions
return a copy of the data with components removed along with an updated `InfoIca` object.

# Example
```julia
# In-place
subtract_ica_components!(dat, ica_result)

# Non-mutating
cleaned_dat, ica_updated = subtract_ica_components(dat, ica_result)

# Select specific components
subtract_ica_components!(dat, ica_result, component_selection = components([1, 3, 5]))
```
"""
function subtract_ica_components!(dat::DataFrame, ica::InfoIca; component_selection::Function = components())
    components_to_remove = get_selected_components(ica, component_selection)
    n_components = size(ica.unmixing, 1)
    if !all(1 .<= components_to_remove .<= n_components)
        throw(ArgumentError("Components must be between 1 and $n_components"))
    end

    # Get data dimensions
    n_channels = length(ica.layout.data.label)

    # Get data and scale it
    data = permutedims(Matrix(dat[!, ica.layout.data.label]))
    data .-= ica.mean
    data ./= ica.scale

    # Get removed activations before transformation and store individually
    all_removed_activations = view(ica.unmixing, components_to_remove, :) * data

    # Store each component's activations separately
    for (i, comp_idx) in enumerate(components_to_remove)
        ica.removed_activations[comp_idx] = all_removed_activations[i:i, :]
    end

    # Pre-compute the transformation matrix
    tra = Matrix(I, n_channels, n_channels) - view(ica.mixing, :, components_to_remove) * view(ica.unmixing, components_to_remove, :)

    # Apply transformation and restore scaling
    cleaned_data = tra * data
    cleaned_data .*= ica.scale
    cleaned_data .+= ica.mean

    # Update DataFrame in-place
    dat[!, ica.layout.data.label] .= permutedims(cleaned_data)

    return nothing
end

"""Contiguous 2D matrix BLAS unmixing for epoched DataFrames (Vector{DataFrame})."""
function subtract_ica_components!(dfs::Vector{DataFrame}, ica::InfoIca; component_selection::Function = components())
    isempty(dfs) && return nothing
    components_to_remove = get_selected_components(ica, component_selection)
    n_components = size(ica.unmixing, 1)
    if !all(1 .<= components_to_remove .<= n_components)
        throw(ArgumentError("Components must be between 1 and $n_components"))
    end

    labels = ica.layout.data.label
    n_channels = length(labels)
    n_epochs = length(dfs)
    n_samples_per_epoch = nrow(dfs[1])
    n_total_samples = n_samples_per_epoch * n_epochs

    # Extract all trial samples into a single contiguous matrix (n_channels × n_total_samples)
    contig_data = Matrix{Float64}(undef, n_channels, n_total_samples)
    for (ch_idx, ch) in enumerate(labels)
        for (ep_idx, df) in enumerate(dfs)
            col = df[!, ch]::Vector{Float64}
            sample_offset = (ep_idx - 1) * n_samples_per_epoch
            @inbounds @simd for i = 1:n_samples_per_epoch
                contig_data[ch_idx, sample_offset+i] = col[i]
            end
        end
    end

    # Scale data
    contig_data .-= ica.mean
    contig_data ./= ica.scale

    # Compute and store removed activations
    all_removed_activations = view(ica.unmixing, components_to_remove, :) * contig_data
    for (i, comp_idx) in enumerate(components_to_remove)
        ica.removed_activations[comp_idx] = all_removed_activations[i:i, :]
    end

    # Pre-compute transformation matrix
    tra = Matrix(I, n_channels, n_channels) - view(ica.mixing, :, components_to_remove) * view(ica.unmixing, components_to_remove, :)

    # Clean data and restore scaling
    cleaned_data = tra * contig_data
    cleaned_data .*= ica.scale
    cleaned_data .+= ica.mean

    # Copy cleaned values back into trial DataFrames
    for (ch_idx, ch) in enumerate(labels)
        for (ep_idx, df) in enumerate(dfs)
            col = df[!, ch]::Vector{Float64}
            sample_offset = (ep_idx - 1) * n_samples_per_epoch
            @inbounds @simd for i = 1:n_samples_per_epoch
                col[i] = cleaned_data[ch_idx, sample_offset+i]
            end
        end
    end

    return nothing
end

function subtract_ica_components(dat::DataFrame, ica::InfoIca; component_selection::Function = components())
    dat_out = copy(dat)
    ica_out = copy(ica)  # Use our custom copy method
    subtract_ica_components!(dat_out, ica_out, component_selection = component_selection)
    return dat_out, ica_out
end

function subtract_ica_components(dat::ContinuousData, ica::InfoIca; component_selection::Function = components())
    dat_out = copy(dat)
    ica_out = copy(ica)  # Use our custom copy method
    subtract_ica_components!(dat_out, ica_out, component_selection = component_selection)
    return dat_out, ica_out
end

function subtract_ica_components!(dat::ContinuousData, ica::InfoIca; component_selection::Function = components())
    subtract_ica_components!(dat.data, ica, component_selection = component_selection)

    components_to_subtract = get_selected_components(ica, component_selection)
    dat.analysis_info.n_ica_components_removed += length(components_to_subtract)

    return nothing
end


"""
    add_ica_components!(dat::DataFrame, ica::InfoIca; component_selection = components())
    add_ica_components!(dat::ContinuousData, ica::InfoIca; component_selection = components())
    add_ica_components(dat::DataFrame, ica::InfoIca; component_selection = components())
    add_ica_components(dat::ContinuousData, ica::InfoIca; component_selection = components())

Restore previously removed ICA components using stored activations.
Mutating (`!`) versions operate in-place; non-mutating versions return a copy and an updated `InfoIca`.
Components must have been removed first (activations stored in `ica.removed_activations`).

# Example
```julia
# In-place
add_ica_components!(dat, ica_result)

# Non-mutating
restored_dat, ica_updated = add_ica_components(dat, ica_result)

# Select specific components to restore
add_ica_components!(dat, ica_result, component_selection = components([1, 3]))
```
"""
function add_ica_components!(dat::DataFrame, ica::InfoIca; component_selection::Function = components())
    components_to_restore = get_selected_components(ica, component_selection)
    n_components = size(ica.unmixing, 1)
    if !all(1 .<= components_to_restore .<= n_components)
        throw(ArgumentError("Components must be between 1 and $n_components"))
    end

    # Check that all components to restore have stored activations
    for comp in components_to_restore
        if !haskey(ica.removed_activations, comp)
            throw(ArgumentError("Component $comp has no stored activations to restore"))
        end
    end

    # Get data and scale it
    data = permutedims(Matrix(dat[!, ica.layout.data.label]))
    data .-= ica.mean
    data ./= ica.scale

    # Get current activations
    activations = ica.unmixing * data

    # Restore each component's stored activations
    for comp in components_to_restore
        activations[comp, :] .= vec(ica.removed_activations[comp])
    end

    # Back to channel space and restore scaling
    restored_data = ica.mixing * activations
    restored_data .*= ica.scale
    restored_data .+= ica.mean

    # Update DataFrame in-place
    dat[!, ica.layout.data.label] .= permutedims(restored_data)

    # Remove restored components from the removed_activations dictionary
    for comp in components_to_restore
        delete!(ica.removed_activations, comp)
    end

    return nothing
end


function add_ica_components(dat::DataFrame, ica::InfoIca; component_selection::Function = components())
    dat_out = copy(dat)
    ica_out = copy(ica)  # Use our custom copy method
    add_ica_components!(dat_out, ica_out, component_selection = component_selection)
    return dat_out, ica_out
end


function add_ica_components(dat::ContinuousData, ica::InfoIca; component_selection::Function = components())
    dat_out = copy(dat)
    ica_out = copy(ica)  # Use our custom copy method
    add_ica_components!(dat_out, ica_out, component_selection = component_selection)
    return dat_out, ica_out
end


function add_ica_components!(dat::ContinuousData, ica::InfoIca; component_selection::Function = components())
    add_ica_components!(dat.data, ica, component_selection = component_selection)

    components_to_add = get_selected_components(ica, component_selection)
    dat.analysis_info.n_ica_components_removed -= length(components_to_add)

    return nothing
end


# Function barrier forces Julia to compile a fast inner loop for the concrete type of `src`
function _copy_col_to_matrix!(dst::Matrix{Float64}, col_idx::Int, src::AbstractVector, selected_samples::Vector{Int})
    @inbounds for (j, sample_idx) in enumerate(selected_samples)
        dst[col_idx, j] = src[sample_idx]
    end
end

"""
    _prepare_ica_data_matrix(dat::ContinuousData, ica::InfoIca, selected_samples::Vector{Int})

Prepare data matrix and calculate ICA components for selected samples.

Extracts relevant channels, permutes to channels × samples format, centers (subtracts mean), 
scales by ICA scale factor, and applies the unmixing matrix to compute components.

# Arguments
- `dat::ContinuousData`: The continuous data
- `ica::InfoIca`: The ICA result object
- `selected_samples::Vector{Int}`: Vector of sample indices to use

# Returns
- `components::Matrix{Float64}`: ICA components (n_components × n_samples)
- `n_components::Int`: Number of ICA components
"""
function _prepare_ica_data_matrix(dat::ContinuousData, ica::InfoIca, selected_samples::Vector{Int})
    relevant_cols = ica.layout.data.label
    n_samples = length(selected_samples)
    n_channels = length(relevant_cols)

    # Pre-allocate to avoid huge DataFrame indexing and permutedims overhead
    dat_matrix = Matrix{Float64}(undef, n_channels, n_samples)
    for (i, col) in enumerate(relevant_cols)
        data_col = dat.data[!, col]
        _copy_col_to_matrix!(dat_matrix, i, data_col, selected_samples)
    end

    dat_matrix .-= mean(dat_matrix, dims = 2)
    dat_matrix ./= ica.scale

    # Calculate components
    components = ica.unmixing * dat_matrix
    n_components = size(components, 1)

    return components, n_components
end



"""
    identify_eog_components(dat::ContinuousData, ica::InfoIca;
                          vEOG_channel::Symbol=:vEOG,
                          hEOG_channel::Symbol=:hEOG,
                          z_threshold::Float64=3.0,
                          two_step::Bool=true,
                          sample_selection::Function = samples())

Identify ICA components potentially related to eye movements based on z-scored correlation with EOG channels.

Uses a two-step approach by default:
1. **Step 1**: Calculate correlations between all ICA components and EOG channels, compute standard z-scores, 
   and identify primary EOG components (those exceeding both the z-score and correlation thresholds).
2. **Step 2**: Calculate spatial correlations (topographies) between remaining components and primary EOG components. 
   Identify secondary EOG components based on spatial correlation threshold. This approach validates that secondary 
   components share similar spatial patterns with primary components, which is more aligned with ICA principles.

This two-step approach helps detect secondary EOG components that share spatial topographies with primary components, 
even if their temporal correlations with EOG channels are not as strong.

# Arguments
- `dat::ContinuousData`: The continuous data containing EOG channels.
- `ica::InfoIca`: The ICA result object.

# Keyword Arguments
- `vEOG_channel::Symbol`: Name of the vertical EOG channel (default: :vEOG).
- `hEOG_channel::Symbol`: Name of the horizontal EOG channel (default: :hEOG).
- `z_threshold::Float64`: Absolute Z-score threshold for step1 identification (default: 3.0). Primary components must exceed 
  this z-score threshold to be identified in step1.
- `min_correlation::Float64`: For step1, minimum absolute correlation threshold with EOG channels (default: 0.6). Primary 
  components must exceed this correlation threshold in addition to the z-score threshold. For step2, this is the minimum 
  spatial correlation threshold with primary component topographies required to identify secondary components.
- `two_step::Bool`: If true (default), use two-step identification: step1 identifies primary components based on z-score and 
  EOG correlation, step2 identifies secondary components based on spatial correlation with primary components. If false, 
  use single-step identification with standard z-scores and correlation thresholds on all components simultaneously.
- `sample_selection::Function`: Function to select samples from `dat.data`. Defaults to `samples()`. Only 
  selected samples are used for correlation calculation.

# Returns
- `Dict{Symbol, Vector{Int}}`: Dictionary containing:
  - `:vEOG`: Vector of component indices identified for vertical eye movements.
  - `:hEOG`: Vector of component indices identified for horizontal eye movements.
- `DataFrame`: DataFrame containing detailed metrics per component:
  - `:Component`: Component index (1 to n_components)
  - `:vEOG_corr`: Absolute correlation with vertical EOG channel
  - `:vEOG_zscore` or `:vEOG_zscore_step1`: Z-score of correlation with vertical EOG channel.
    When `two_step=true`, only step1 z-scores are provided.
  - `:vEOG_spatial_corr`: (Only when `two_step=true`) Maximum spatial correlation (topography) with any primary vEOG 
    component. `NaN` for primary vEOG components identified in step1.
  - `:vEOG_spatial_corr_z`: Z-score of `:vEOG_spatial_corr`.
  - `:vEOG_temporal_corr`: (Only when `two_step=true`) Maximum lagged temporal correlation with any primary vEOG 
    component. `NaN` for primary vEOG components identified in step1.
  - `:vEOG_temporal_corr_z`: Z-score of `:vEOG_temporal_corr`.
  - `:hEOG_corr`: Absolute correlation with horizontal EOG channel
  - `:hEOG_zscore` or `:hEOG_zscore_step1`: Z-score of correlation with horizontal EOG channel.
    When `two_step=true`, only step1 z-scores are provided.
  - `:hEOG_spatial_corr`: (Only when `two_step=true`) Maximum spatial correlation (topography) with any primary hEOG 
    component. `NaN` for primary hEOG components identified in step1.
  - `:hEOG_spatial_corr_z`: Z-score of `:hEOG_spatial_corr`.
  - `:hEOG_temporal_corr`: (Only when `two_step=true`) Maximum lagged temporal correlation with any primary hEOG 
    component. `NaN` for primary hEOG components identified in step1.
  - `:hEOG_temporal_corr_z`: Z-score of `:hEOG_temporal_corr`.

# Examples
```julia
# Basic usage with default two-step approach
eog_comps, metrics = identify_eog_components(dat, ica_result)

# Use single-step identification
eog_comps, metrics = identify_eog_components(dat, ica_result, two_step=false)

# Custom threshold and sample selection
eog_comps, metrics = identify_eog_components(dat, ica_result, 
    z_threshold=2.5,
    sample_selection=samples_not(:is_extreme_value_100))
```
"""
function identify_eog_components(
    dat::ContinuousData,
    ica::InfoIca;
    vEOG_channel::Symbol = :vEOG,
    hEOG_channel::Symbol = :hEOG,
    sample_selection::Function = samples(),
    interval_selection::Interval = times(),
    z_threshold::Real = 3.0,
    min_correlation::Float64 = 0.5,
    two_step::Bool = true,
    _precomputed_components::Union{Nothing,Tuple{Matrix{Float64},Int}} = nothing,
    _precomputed_samples::Union{Nothing,Vector{Int}} = nothing,
)

    # Check basic inputs - return nothing as first value if EOG channels are missing
    # TODO: what about when only one is available?
    if !(vEOG_channel in propertynames(dat.data)) || !(hEOG_channel in propertynames(dat.data))
        return nothing, DataFrame()
    end

    # Get samples to use (skip if pre-computed)
    selected_samples = if !isnothing(_precomputed_samples)
        _precomputed_samples
    else
        combined_sel = _combine_interval_sample(interval_selection, sample_selection)
        get_selected_samples(dat, combined_sel)
    end
    if isempty(selected_samples)
        @minimal_warning "No samples remaining after applying exclude criteria. Cannot identify eye components."
        return Dict(:vEOG => Int[], :hEOG => Int[]), DataFrame()
    end

    # Function to calculate initial EOG-component correlations 
    """Compute absolute correlation of each ICA component with the given EOG signal."""
    function calculate_correlations(eog_signal)
        corrs = zeros(n_components)
        comp_buf = Vector{Float64}(undef, size(components, 2))
        for comp_idx = 1:n_components
            comp_buf .= @view components[comp_idx, :]
            corrs[comp_idx] = abs(cor(comp_buf, eog_signal))
        end
        return corrs
    end

    # Function to calculate spatial correlations between identified vEOG/hEOG components and remaining components (used in Step2)
    """Compute spatial (topography) correlations between remaining and primary EOG components."""
    function calculate_spatial_correlations(remaining_components, primary_components, mixing_matrix, min_corr_threshold)

        spatial_corrs = fill(NaN, n_components)
        secondary_components = Int[]

        if !isempty(remaining_components) && !isempty(primary_components)
            # Calculate spatial correlations (topographies) between remaining components and each primary component taking max correlation
            spatial_corrs_remaining = zeros(length(remaining_components))
            for (idx, comp_idx) in enumerate(remaining_components)
                comp_topography = mixing_matrix[:, comp_idx]
                max_corr = 0.0
                for primary_comp in primary_components
                    primary_topography = mixing_matrix[:, primary_comp]
                    corr_val = abs(cor(comp_topography, primary_topography))
                    max_corr = max(max_corr, corr_val)
                end
                spatial_corrs_remaining[idx] = max_corr
            end

            # Map spatial correlations back to original component indices
            spatial_corrs[remaining_components] = spatial_corrs_remaining

            # Find secondary components (spatial correlation > min_correlation)
            secondary_idx = findall(spatial_corrs_remaining .> min_corr_threshold)
            secondary_components = remaining_components[secondary_idx]
        end

        return spatial_corrs, secondary_components
    end

    # Function to calculate lagged correlations between components without huge crosscor allocations
    """Compute maximum absolute lagged cross-correlation between remaining and primary components."""
    function calculate_lagged_correlations(remaining_components, primary_components, comp_matrix, max_lag_samples, lag_step)
        lagged_corrs = fill(NaN, n_components)
        if isempty(remaining_components) || isempty(primary_components)
            return lagged_corrs
        end

        lags = (-max_lag_samples):lag_step:max_lag_samples
        n_s = size(comp_matrix, 2)

        # Pre-allocate contiguous buffers to eliminate cache misses from row-major stride
        comp_buf = Vector{Float64}(undef, n_s)
        eog_buf = Vector{Float64}(undef, n_s)

        for comp_idx in remaining_components
            comp_buf .= @view comp_matrix[comp_idx, :]
            m_c = mean(comp_buf)
            s_c = std(comp_buf, corrected = false)

            max_corr = 0.0
            for eog_comp_idx in primary_components
                eog_buf .= @view comp_matrix[eog_comp_idx, :]
                m_e = mean(eog_buf)
                s_e = std(eog_buf, corrected = false)
                denom = n_s * s_c * s_e

                # Manual zero-allocation cross correlation for the few specific lags
                max_corr_val = 0.0
                for lag in lags
                    sum_xy = 0.0
                    @inbounds @simd for i = max(1, 1 - lag):min(n_s, n_s - lag)
                        sum_xy += (eog_buf[i] - m_e) * (comp_buf[i+lag] - m_c)
                    end
                    c = abs(sum_xy / denom)
                    if c > max_corr_val
                        max_corr_val = c
                    end
                end

                if max_corr_val > max_corr
                    max_corr = max_corr_val
                end
            end

            lagged_corrs[comp_idx] = max_corr
        end

        return lagged_corrs
    end

    # Get EOG signals for valid samples only
    vEOG = dat.data[selected_samples, vEOG_channel]
    hEOG = dat.data[selected_samples, hEOG_channel]

    # Calculate components for valid samples (skip if pre-computed)
    components, n_components = if !isnothing(_precomputed_components)
        _precomputed_components
    else
        _prepare_ica_data_matrix(dat, ica, selected_samples)
    end

    # Step 1: Always calculate correlations and z-scores (used in both single-step and two-step modes)
    vEOG_corrs = calculate_correlations(vEOG)
    hEOG_corrs = calculate_correlations(hEOG)

    vEOG_corr_z = StatsBase.zscore(vEOG_corrs)
    hEOG_corr_z = StatsBase.zscore(hEOG_corrs)

    # Identify primary components (components meeting both z-score and correlation thresholds)
    primary_vEOG = findall((abs.(vEOG_corr_z) .> z_threshold) .& (abs.(vEOG_corrs) .> min_correlation))
    primary_hEOG = findall((abs.(hEOG_corr_z) .> z_threshold) .& (abs.(hEOG_corrs) .> min_correlation))

    # Step 2: Always calculate spatial and temporal correlations between identified components and remaining
    remaining_vEOG = setdiff(1:n_components, primary_vEOG)
    remaining_hEOG = setdiff(1:n_components, primary_hEOG)

    # Calculate vEOG/hEOG step2: spatial correlations with primary vEOG/hEOG components
    vEOG_spatial_corr, _ = calculate_spatial_correlations(remaining_vEOG, primary_vEOG, ica.mixing, min_correlation)
    hEOG_spatial_corr, _ = calculate_spatial_correlations(remaining_hEOG, primary_hEOG, ica.mixing, min_correlation)

    # Lag range: +- 100ms  
    # Convert to samples based on sample rate
    max_lag_samples = round(Int, 100.0 * dat.sample_rate / 1000.0) # 100 ms lag
    lag_step = max(1, round(Int, 5.0 * dat.sample_rate / 1000.0)) # 5 ms step

    # Calculate lagged correlations (zero-allocation, uses pre-computed component means)
    vEOG_temporal_corr = calculate_lagged_correlations(remaining_vEOG, primary_vEOG, components, max_lag_samples, lag_step)
    hEOG_temporal_corr = calculate_lagged_correlations(remaining_hEOG, primary_hEOG, components, max_lag_samples, lag_step)

    # Helper function to compute z-scores only for remaining components
    """Z-score a correlation vector using only the values at `remaining_indices`."""
    function zscore_for_remaining(corr_vector, remaining_indices)
        z = fill(NaN, length(corr_vector))
        if length(remaining_indices) >= 2
            remaining_vals = corr_vector[remaining_indices]
            if std(remaining_vals) > 0.0
                z[remaining_indices] = StatsBase.zscore(remaining_vals)
            end
        end
        return z
    end

    # Calculate z-scores for spatial and temporal correlations
    vEOG_spatial_corr_z = zscore_for_remaining(vEOG_spatial_corr, remaining_vEOG)
    vEOG_temporal_corr_z = zscore_for_remaining(vEOG_temporal_corr, remaining_vEOG)
    hEOG_spatial_corr_z = zscore_for_remaining(hEOG_spatial_corr, remaining_hEOG)
    hEOG_temporal_corr_z = zscore_for_remaining(hEOG_temporal_corr, remaining_hEOG)

    # Identify secondary components: require BOTH correlation > min_correlation AND z-score > z_threshold
    """Return indices of remaining components exceeding both correlation and z-score thresholds."""
    function identify_secondary_components(remaining_components, corr_vector, corr_z_vector, min_corr_threshold, z_threshold)
        isempty(remaining_components) && return Int[]
        corr_mask = corr_vector[remaining_components] .> min_corr_threshold
        z_mask = abs.(corr_z_vector[remaining_components]) .> z_threshold
        combined_mask = corr_mask .& z_mask
        return remaining_components[combined_mask]
    end

    # Identify secondary components for vEOG: spatial OR temporal (each requires both corr and z-score)
    secondary_vEOG_spatial =
        identify_secondary_components(remaining_vEOG, vEOG_spatial_corr, vEOG_spatial_corr_z, min_correlation, z_threshold)
    secondary_vEOG_temporal =
        identify_secondary_components(remaining_vEOG, vEOG_temporal_corr, vEOG_temporal_corr_z, min_correlation, z_threshold)

    # Identify secondary components for hEOG: spatial OR temporal (each requires both corr and z-score)
    secondary_hEOG_spatial =
        identify_secondary_components(remaining_hEOG, hEOG_spatial_corr, hEOG_spatial_corr_z, min_correlation, z_threshold)
    secondary_hEOG_temporal =
        identify_secondary_components(remaining_hEOG, hEOG_temporal_corr, hEOG_temporal_corr_z, min_correlation, z_threshold)

    # Combine spatial and temporal secondary components
    secondary_vEOG = union(secondary_vEOG_spatial, secondary_vEOG_temporal)
    secondary_hEOG = union(secondary_hEOG_spatial, secondary_hEOG_temporal)

    # Combine primary and secondary components (secondary only added if two_step is true)
    identified_vEOG = two_step ? union(primary_vEOG, secondary_vEOG) : primary_vEOG
    identified_hEOG = two_step ? union(primary_hEOG, secondary_hEOG) : primary_hEOG

    # DataFrame metrics
    metrics_df = DataFrame(
        :Component => 1:n_components,
        :vEOG_corr => vEOG_corrs,
        :vEOG_zscore => vEOG_corr_z,
        :vEOG_spatial_corr => vEOG_spatial_corr,
        :vEOG_spatial_corr_z => vEOG_spatial_corr_z,
        :vEOG_temporal_corr => vEOG_temporal_corr,
        :vEOG_temporal_corr_z => vEOG_temporal_corr_z,
        :hEOG_corr => hEOG_corrs,
        :hEOG_zscore => hEOG_corr_z,
        :hEOG_spatial_corr => hEOG_spatial_corr,
        :hEOG_spatial_corr_z => hEOG_spatial_corr_z,
        :hEOG_temporal_corr => hEOG_temporal_corr,
        :hEOG_temporal_corr_z => hEOG_temporal_corr_z,
    )

    sort!(identified_vEOG)
    sort!(identified_hEOG)

    result_dict = Dict{Symbol,Vector{Int}}(:vEOG => identified_vEOG, :hEOG => identified_hEOG)

    return result_dict, metrics_df
end

"""
    identify_ecg_components(dat::ContinuousData, ica::InfoIca;
                              min_bpm::Real=50, max_bpm::Real=110,
                              min_prominence_std::Real=4,
                              min_peaks::Int=100,
                              max_ibi_std_s::Real=0.2,
                              min_peak_ratio::Real=0.5,
                              sample_selection::Function = samples(),
    interval_selection::Interval = times(),
)

Identify ICA components potentially related to ECG artifacts based on peak detection
and interval regularity, using only samples consistent with ICA calculation.

# Arguments
- `dat::ContinuousData`: The continuous data (needed for sampling rate and sample selection columns).
- `ica::InfoIca`: The ICA result object.

# Keyword Arguments
- `min_bpm::Real`: Minimum plausible heart rate in beats per minute (default: 40).
- `max_bpm::Real`: Maximum plausible heart rate in beats per minute (default: 120).
- `min_prominence_std::Real`: Minimum peak prominence in standard deviations above mean for z-scored time series (default: 4).
- `min_peaks::Int`: Minimum number of valid inter-beat intervals required (default: 100). 
  Note: Since `num_valid_ibis` is the number of valid intervals between peaks (which is `num_peaks - 1` if all are valid),
  the check uses `num_valid_ibis >= (min_peaks - 1)` to account for this relationship.
- `max_ibi_std_s::Real`: Maximum standard deviation of the inter-beat intervals (in seconds) for component to be flagged (default: 0.2).
- `min_peak_ratio::Real`: Minimum ratio of valid inter-beat intervals to total inter-peak intervals (default: 0.5). 
  This ensures that a sufficient proportion of detected peaks fall within the plausible heart rate range.
- `sample_selection::Function`: Function to select samples from `dat.data`. Defaults to `samples()`. Only 
  selected samples are used for component calculation and peak detection.

# Returns
- `Vector{Int}`: Sorted vector of indices identified as potential ECG components.
- `DataFrame`: DataFrame containing metrics for each component (calculated on the selected samples):
  - `:Component`: Component index (1 to n).
  - `:num_peaks`: Number of detected prominent peaks.
  - `:num_valid_ibis`: Number of inter-beat intervals within the plausible BPM range.
  - `:mean_ibi_s`: Mean inter-beat interval in seconds (if num_valid_ibis > 0).
  - `:std_ibi_s`: Standard deviation of inter-beat intervals in seconds (if num_valid_ibis > 1).
  - `:peak_ratio`: Ratio of valid inter-beat intervals to total inter-peak intervals.
  - `:heart_rate_bpm`: Estimated heart rate in beats per minute (if mean_ibi_s is valid).
"""
function identify_ecg_components(
    dat::ContinuousData,
    ica::InfoIca;
    min_bpm::Real = 40,
    max_bpm::Real = 120,
    min_prominence_std::Real = 4,
    min_peaks::Int = 100,
    max_ibi_std_s::Real = 0.2,
    min_peak_ratio::Real = 0.5,
    sample_selection::Function = samples(),
    interval_selection::Interval = times(),
    _precomputed_components::Union{Nothing,Tuple{Matrix{Float64},Int}} = nothing,
    _precomputed_samples::Union{Nothing,Vector{Int}} = nothing,
)

    # Data Preparation (skip if pre-computed)
    selected_samples = if !isnothing(_precomputed_samples)
        _precomputed_samples
    else
        combined_sel = _combine_interval_sample(interval_selection, sample_selection)
        get_selected_samples(dat, combined_sel)
    end
    if isempty(selected_samples)
        @minimal_warning "No samples remaining after applying exclude criteria."
        return Int[], DataFrame()
    end

    # Calculate components for valid samples (skip if pre-computed)
    components_subset, n_components = if !isnothing(_precomputed_components)
        _precomputed_components
    else
        _prepare_ica_data_matrix(dat, ica, selected_samples)
    end

    # Early return if no components
    if n_components == 0
        return Int[],
        DataFrame(
            Component = Int[],
            num_peaks = Int[],
            num_valid_ibis = Int[],
            mean_ibi_s = Float64[],
            std_ibi_s = Float64[],
            peak_ratio = Float64[],
            heart_rate_bpm = Float64[],
        )
    end

    # Convert BPM to plausible IBI range
    min_ibi_s = 60.0 / max_bpm
    max_ibi_s = 60.0 / min_bpm

    metrics = []
    identified_ecg = Int[]

    # Pre-allocate buffer for z-scoring
    n_samples = size(components_subset, 2)
    ts_zscored = Vector{Float64}(undef, n_samples)

    for comp_idx = 1:n_components
        component_ts = @view components_subset[comp_idx, :]

        # Z-score the time series for consistent peak detection across components
        m = mean(component_ts)
        s = std(component_ts)
        ts_zscored .= (component_ts .- m) ./ s

        # Find prominent peaks 
        peak_indices = _find_peaks(ts_zscored; min_prominence_std = min_prominence_std)
        num_peaks = length(peak_indices)

        # Calculate IBI metrics
        num_valid_ibis, mean_ibi, std_ibi, peak_ratio = _calculate_ibi_metrics(peak_indices, dat.sample_rate, min_ibi_s, max_ibi_s)

        # Check if component meets ECG criteria
        has_sufficient_ibis = num_valid_ibis >= (min_peaks - 1)
        has_low_ibi_variability = std_ibi <= max_ibi_std_s
        has_good_peak_ratio = peak_ratio >= min_peak_ratio

        is_ecg = has_sufficient_ibis && has_low_ibi_variability && has_good_peak_ratio
        is_ecg && push!(identified_ecg, comp_idx)

        # Calculate heart rate if we have valid IBI
        if isnan(mean_ibi) || mean_ibi <= 0
            heart_rate_bpm = NaN
        else
            heart_rate_bpm = 60.0 / mean_ibi
        end

        # Store metrics
        push!(
            metrics,
            (
                Component = comp_idx,
                num_peaks = num_peaks,
                num_valid_ibis = num_valid_ibis,
                mean_ibi_s = mean_ibi,
                std_ibi_s = std_ibi,
                peak_ratio = peak_ratio,
                heart_rate_bpm = heart_rate_bpm,
            ),
        )
    end

    # Build metrics DataFrame
    metrics_df = DataFrame(metrics)

    sort!(identified_ecg)

    return identified_ecg, metrics_df

end

"""
    _calculate_ibi_metrics(peak_indices::Vector{Int}, sample_rate::Real, min_ibi_s::Real, max_ibi_s::Real)

Calculate inter-beat interval (IBI) metrics from peak indices.

# Arguments
- `peak_indices::Vector{Int}`: Indices of detected peaks
- `sample_rate::Real`: Sampling rate in Hz
- `min_ibi_s::Real`: Minimum valid IBI in seconds
- `max_ibi_s::Real`: Maximum valid IBI in seconds

# Returns
- `num_valid_ibis::Int`: Number of valid IBIs
- `mean_ibi::Float64`: Mean IBI in seconds (NaN if insufficient data)
- `std_ibi::Float64`: Standard deviation of IBIs in seconds (NaN if insufficient data)
- `peak_ratio::Float64`: Ratio of valid IBIs to total inter-peak intervals
"""
function _calculate_ibi_metrics(peak_indices::Vector{Int}, sample_rate::Real, min_ibi_s::Real, max_ibi_s::Real)

    num_peaks = length(peak_indices)
    num_peaks < 2 && return 0, NaN, NaN, 0.0

    # Calculate IBIs
    ibis_s = diff(peak_indices) ./ sample_rate

    # Filter valid IBIs
    valid_ibi_mask = (ibis_s .>= min_ibi_s) .& (ibis_s .<= max_ibi_s)
    valid_ibis = ibis_s[valid_ibi_mask]
    num_valid_ibis = length(valid_ibis)

    # Calculate peak ratio
    peak_ratio = num_valid_ibis / (num_peaks - 1)

    # Calculate statistics
    if num_valid_ibis > 1
        mean_ibi = mean(valid_ibis)
        std_ibi = std(valid_ibis)
    elseif num_valid_ibis == 1
        mean_ibi = valid_ibis[1]
        std_ibi = 0.0
    else
        mean_ibi = NaN
        std_ibi = NaN
    end

    return num_valid_ibis, mean_ibi, std_ibi, peak_ratio
end




"""
    identify_spatial_kurtosis_components(ica::InfoIca; z_threshold::Float64 = 3.0)

Identify ICA components with high spatial kurtosis (localized, spot-like activity).

Spatial kurtosis measures how localized the component's topography is. High spatial kurtosis
indicates that the component's activity is concentrated in a small number of channels (spot-like),
which is characteristic of channel noise or artifacts.

# Arguments
- `ica::InfoIca`: The ICA result object.

# Keyword Arguments
- `z_threshold::Float64`: Z-score threshold for identifying high spatial kurtosis components (default: 3.0).

# Returns
- `Vector{Int}`: Indices of components with high spatial kurtosis.
- `DataFrame`: DataFrame containing spatial kurtosis values and z-scores for all components.
"""
function identify_spatial_kurtosis_components(ica::InfoIca; z_threshold::Real = 3.0)

    # Calculate spatial kurtosis for each component's weights
    n_components = size(ica.mixing, 2)
    spatial_kurtosis = Float64[]
    for i = 1:n_components
        weights = ica.mixing[:, i]
        k = kurtosis(weights)
        push!(spatial_kurtosis, k)
    end

    # Calculate z-scores of spatial kurtosis values
    z_spatial_kurtosis = StatsBase.zscore(spatial_kurtosis)

    # Identify components with high spatial kurtosis (using z-scores)
    high_kurtosis_comps = findall(z_spatial_kurtosis .> z_threshold)  # Only positive deviations (localized activity)
    sort!(high_kurtosis_comps)

    # Create metrics DataFrame
    metrics_df = DataFrame(:Component => 1:n_components, :spatial_kurtosis => spatial_kurtosis, :z_spatial_kurtosis => z_spatial_kurtosis)

    return high_kurtosis_comps, metrics_df
end


"""
    identify_line_noise_components(dat::ContinuousData, ica::InfoIca;
                                 line_freq::Real=50.0,
                                 freq_bandwidth::Real=1.0,
                                 z_threshold::Float64=3.0,
                                 min_harmonic_power::Real=0.5)

Identify ICA components with strong line noise characteristics using spectral peakiness.

Spectral peakiness compares power at the line frequency to nearby flanking frequencies
(±5 Hz excluding the line band), rather than to the entire broadband spectrum. This avoids
false positives from components with generally elevated power that happens to include the
line frequency range.

# Arguments
- `dat::ContinuousData`: The continuous data.
- `ica::InfoIca`: The ICA result object.

# Keyword Arguments
- `sample_selection::Function`: Function to select samples (default: `samples()`).
- `line_freq::Real`: Line frequency in Hz (default: 50.0 for European power).
- `freq_bandwidth::Real`: Bandwidth around line frequency to consider (default: 1.0 Hz).
- `z_threshold::Float64`: Z-score threshold for identifying line noise components (default: 3.0).
- `min_harmonic_power::Real`: Minimum peakiness ratio of 2nd harmonic (100 Hz vs. its flanking frequencies) (default: 1.5).

# Returns
- `Vector{Int}`: Indices of components with strong line noise characteristics.
- `DataFrame`: DataFrame containing spectral metrics for all components:
  - `:Component`: Component index
  - `:line_power`: Mean power at the line frequency band
  - `:flanking_power`: Mean power at flanking frequencies (used as baseline for peakiness)
  - `:power_ratio`: Spectral peakiness ratio (line_power / flanking_power)
  - `:harmonic_ratio`: Spectral peakiness ratio at the 2nd harmonic
  - `:power_ratio_zscore`: Z-score of power_ratio across components
"""
function identify_line_noise_components(
    dat::ContinuousData,
    ica::InfoIca;
    sample_selection::Function = samples(),
    interval_selection::Interval = times(),
    line_freq::Real = 50.0,
    freq_bandwidth::Real = 1.0,
    z_threshold::Real = 3.0,
    min_harmonic_power::Real = 1.5,
    _precomputed_components::Union{Nothing,Tuple{Matrix{Float64},Int}} = nothing,
    _precomputed_samples::Union{Nothing,Vector{Int}} = nothing,
)

    # Get samples to use (skip if pre-computed)
    samples_to_use = if !isnothing(_precomputed_samples)
        _precomputed_samples
    else
        combined_sel = _combine_interval_sample(interval_selection, sample_selection)
        get_selected_samples(dat, combined_sel)
    end
    if isempty(samples_to_use)
        @minimal_warning "No samples remaining after applying exclude criteria. Cannot identify line noise components."
        return Int[], DataFrame()
    end

    # Calculate components for valid samples (skip if pre-computed)
    components, n_components = if !isnothing(_precomputed_components)
        _precomputed_components
    else
        _prepare_ica_data_matrix(dat, ica, samples_to_use)
    end

    # Calculate power spectrum for each component
    # Use a reasonable FFT size (power of 2, capped at 2^16 points for performance)
    n_samples = size(components, 2)
    nfft = min(nextpow(2, n_samples), 65536)  # Properly cap at 2^16 points
    freqs = FFTW.rfftfreq(nfft, dat.sample_rate)
    n_freqs = length(freqs)
    psd = zeros(n_freqs, n_components)

    # Pre-allocate FFT buffer and plan FFT
    signal_fft = zeros(Float64, nfft)
    fft_plan = FFTW.plan_rfft(signal_fft)

    for i = 1:n_components
        signal = @view components[i, :]

        # Prepare signal for FFT: truncate or zero-pad to nfft points
        if n_samples > nfft
            signal_fft .= @view signal[1:nfft]
        elseif n_samples < nfft
            signal_fft[1:n_samples] .= signal
            signal_fft[(n_samples+1):end] .= 0.0
        else
            signal_fft .= signal
        end

        # Calculate power spectral density
        psd[:, i] .= abs2.(fft_plan * signal_fft)
    end

    # Define frequency bands
    # Flanking bandwidth for spectral peakiness (Hz on each side of line freq, excluding the line band itself)
    flank_bandwidth = 5.0

    # Line frequency band (narrow band around 50 Hz)
    line_band = findall(abs.(freqs .- line_freq) .<= freq_bandwidth)

    # Flanking bands: frequencies near the line frequency but excluding the line band
    # e.g., for 50 Hz with freq_bandwidth=1 and flank_bandwidth=5: uses 45-49 Hz and 51-55 Hz
    line_flank_band =
        findall((abs.(freqs .- line_freq) .> freq_bandwidth) .& (abs.(freqs .- line_freq) .<= freq_bandwidth + flank_bandwidth))

    # 2nd harmonic bands
    harmonic2_freq = line_freq * 2
    harmonic2_band = findall(abs.(freqs .- harmonic2_freq) .<= freq_bandwidth)
    harmonic2_flank_band =
        findall((abs.(freqs .- harmonic2_freq) .> freq_bandwidth) .& (abs.(freqs .- harmonic2_freq) .<= freq_bandwidth + flank_bandwidth))

    # Calculate metrics for each component
    metrics = []
    for i = 1:n_components
        # Get power at line frequency band
        line_power = mean(psd[line_band, i])

        # Spectral peakiness: compare line frequency power to nearby flanking frequencies
        flanking_power = !isempty(line_flank_band) ? mean(psd[line_flank_band, i]) : mean(psd[:, i])

        # Power ratio: how much does 50 Hz stand out above its immediate neighbours
        power_ratio = line_power / (flanking_power + eps())

        # Harmonic peakiness: does 100 Hz also stand out above its flanking frequencies
        harmonic_ratio = if !isempty(harmonic2_band)
            harmonic2_power = mean(psd[harmonic2_band, i])
            harmonic2_flanking = !isempty(harmonic2_flank_band) ? mean(psd[harmonic2_flank_band, i]) : flanking_power
            harmonic2_power / (harmonic2_flanking + eps())
        else
            NaN  # 2nd harmonic frequency above Nyquist
        end

        # Store metrics
        push!(
            metrics,
            (
                Component = i,
                line_power = line_power,
                flanking_power = flanking_power,
                power_ratio = power_ratio,
                harmonic_ratio = harmonic_ratio,
            ),
        )
    end

    # Create metrics DataFrame
    metrics_df = DataFrame(metrics)

    # Calculate z-scores of power ratios
    power_ratio_z = StatsBase.zscore(metrics_df.power_ratio)
    metrics_df[!, :power_ratio_zscore] = power_ratio_z

    # Identify components with strong line noise characteristics
    # Step 1: High power ratio (z-score > threshold)
    high_power_ratio_mask = power_ratio_z .> z_threshold

    # Step 2: Must have significant 2nd harmonic (> min_harmonic_power)
    has_harmonic_mask = metrics_df.harmonic_ratio .> min_harmonic_power

    # Combine both criteria
    line_noise_comps = findall(high_power_ratio_mask .& has_harmonic_mask)
    sort!(line_noise_comps)

    return line_noise_comps, metrics_df

end

"""
    ArtifactComponents

A structure to hold all identified artifact components from ICA analysis.

# Fields
- `eog::Dict{Symbol, Vector{Int}}`: Dictionary with :vEOG and :hEOG keys containing identified eye movement components
- `ecg::Vector{Int}`: Vector of identified ECG/heartbeat components  
- `line_noise::Vector{Int}`: Vector of identified line noise components
- `channel_noise::Vector{Int}`: Vector of identified high spatial kurtosis (channel noise) components
"""
struct ArtifactComponents
    eog::Dict{Symbol,Vector{Int}}
    ecg::Vector{Int}
    line_noise::Vector{Int}
    channel_noise::Vector{Int}
end

"""
    combine_artifact_components(eog_comps, ecg_comps, line_noise_comps, channel_noise_comps)

Combine all identified artifact components into a single ArtifactComponents structure.

# Arguments
- `eog_comps::Dict{Symbol, Vector{Int}}`: EOG components dictionary
- `ecg_comps::Vector{Int}`: ECG components vector
- `line_noise_comps::Vector{Int}`: Line noise components vector
- `channel_noise_comps::Vector{Int}`: Channel noise components vector

# Returns
- `ArtifactComponents`: Combined structure containing all artifact components
"""
function combine_artifact_components(
    eog_comps::Dict{Symbol,Vector{Int}},
    ecg_comps::Vector{Int},
    line_noise_comps::Vector{Int},
    channel_noise_comps::Vector{Int},
)
    return ArtifactComponents(eog_comps, ecg_comps, line_noise_comps, channel_noise_comps)
end


"""
    identify_components(dat::ContinuousData, ica::InfoIca; method = :correlation,
    sample_selection::Function = samples(), interval_selection::Interval = times(), kwargs...)

Identify all types of artifact components in one unified call.

# Arguments
- `dat::ContinuousData`: The continuous data
- `ica::InfoIca`: The ICA result object
- `method::Symbol`: Component identification method:
  - `:correlation` (default): Correlation-based detection (EOG, ECG, line noise, spatial kurtosis).
    The original EegFun approach using z-scored correlations with EOG/ECG channels,
    spectral analysis for line noise, and spatial kurtosis for channel noise.
- `sample_selection::Function`: Sample selection function (default: samples())
- `kwargs...`: Additional keyword arguments passed to individual identification functions

# Returns
- `ArtifactComponents`: Combined structure containing all identified artifact components
- `Dict{Symbol, DataFrame}`: Dictionary containing metrics DataFrames for each component type:
  - `:eog_metrics`: EOG component metrics
  - `:ecg_metrics`: ECG component metrics  
  - `:line_noise_metrics`: Line noise component metrics
  - `:channel_noise_metrics`: Channel noise component metrics

# Examples
```julia
# Default correlation-based identification
artifacts, metrics = identify_components(dat, ica_result)

# With custom sample selection
artifacts, metrics = identify_components(dat, ica_result, 
    sample_selection = samples_not(:is_extreme_value_100))
```
"""
function identify_components(
    dat::ContinuousData,
    ica::InfoIca;
    method::Symbol = :correlation,
    sample_selection::Function = samples(),
    interval_selection::Interval = times(),
    kwargs...,
)
    if method == :correlation
        return _identify_components_correlation(
            dat,
            ica;
            sample_selection = sample_selection,
            interval_selection = interval_selection,
            kwargs...,
        )
    else
        error("Unknown component identification method: :$method. " * "Supported methods: :correlation")
    end
end

"""
    _identify_components_correlation(dat, ica; kwargs...) -> (ArtifactComponents, Dict)

The original correlation-based component identification method.
This is the internal implementation dispatched by `identify_components(method = :correlation)`.
"""
function _identify_components_correlation(
    dat::ContinuousData,
    ica::InfoIca;
    sample_selection::Function = samples(),
    interval_selection::Interval = times(),
    kwargs...,
)
    # Pre-compute shared data ONCE (avoids 3x redundant matrix operations)
    combined_sel = _combine_interval_sample(interval_selection, sample_selection)
    selected_samples = get_selected_samples(dat, combined_sel)
    precomputed = _prepare_ica_data_matrix(dat, ica, selected_samples)

    # Identify EOG components (pass through use_robust_zscore if provided)
    eog_comps, eog_metrics_df = identify_eog_components(
        dat,
        ica;
        sample_selection = sample_selection,
        _precomputed_components = precomputed,
        _precomputed_samples = selected_samples,
        kwargs...,
    )

    # Identify ECG components
    ecg_comps, ecg_metrics_df = identify_ecg_components(
        dat,
        ica;
        sample_selection = sample_selection,
        _precomputed_components = precomputed,
        _precomputed_samples = selected_samples,
        kwargs...,
    )

    # Identify line noise components
    line_noise_comps, line_noise_metrics_df = identify_line_noise_components(
        dat,
        ica;
        sample_selection = sample_selection,
        _precomputed_components = precomputed,
        _precomputed_samples = selected_samples,
        kwargs...,
    )

    # Identify channel noise components (spatial kurtosis - uses only ica.mixing, no data matrix needed)
    channel_noise_comps, channel_noise_metrics_df = identify_spatial_kurtosis_components(ica; kwargs...)

    # Combine all components
    artifacts = combine_artifact_components(eog_comps, ecg_comps, line_noise_comps, channel_noise_comps)

    # Combine all metrics
    metrics = Dict{Symbol,DataFrame}(
        :eog_metrics => eog_metrics_df,
        :ecg_metrics => ecg_metrics_df,
        :line_noise_metrics => line_noise_metrics_df,
        :channel_noise_metrics => channel_noise_metrics_df,
    )

    return artifacts, metrics
end

"""
    get_all_components(artifacts::ArtifactComponents)

Get all unique artifact component indices as a single vector.

# Arguments
- `artifacts::ArtifactComponents`: The artifact components structure

# Returns
- `Vector{Int}`: Sorted vector of all unique artifact component indices
"""
function get_all_ica_components(artifacts::ArtifactComponents)
    all_comps = Set{Int}()

    # Add EOG components
    for (_, comps) in artifacts.eog
        union!(all_comps, comps)
    end

    # Add other components
    union!(all_comps, artifacts.ecg)
    union!(all_comps, artifacts.line_noise)
    union!(all_comps, artifacts.channel_noise)

    return sort(collect(all_comps))
end



function Base.show(io::IO, artifacts::ArtifactComponents)
    # Get all components for and print short summary
    all_comps = get_all_ica_components(artifacts)
    println(io, "ICA Artifact Components")
    println(io, "vEOG: $(artifacts.eog[:vEOG])")
    println(io, "hEOG: $(artifacts.eog[:hEOG])")
    println(io, "ECG: $(artifacts.ecg)")
    println(io, "Line Noise: $(artifacts.line_noise)")
    println(io, "Channel Noise: $(artifacts.channel_noise)")
    println(io, "All: $all_comps")
end
