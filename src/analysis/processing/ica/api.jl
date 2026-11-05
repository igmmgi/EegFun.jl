"""
    run_ica(dat::ContinuousData; n_components = nothing, sample_selection = samples(),
            interval_selection = times(), channel_selection = channels(), include_extra = false,
            percentage_of_data = 100.0, algorithm = :infomax, use_gpu = false,
            params = IcaPrms(), rng = Random.default_rng()) -> InfoIca
    run_ica(epoched_data::AbstractVector{<:EpochData}; remove_duplicates = true, kwargs...) -> InfoIca
    run_ica(epoched_data::EpochData; kwargs...) -> InfoIca

Run Independent Component Analysis (ICA) on continuous or epoched EEG data (`:infomax`,
`:infomax_extended`, `:picard`, `:picard_extended`, or `:amica`), defaulting `n_components`
to the effective numerical rank of the data.

# Examples
```julia
ica_result = run_ica(dat)
ica_result = run_ica(dat; sample_selection = samples_not(:is_extreme_value_100), percentage_of_data = 25.0)
ica_result = run_ica(epochs_vec; algorithm = :picard_extended)
```
"""
function run_ica(
    dat::ContinuousData;
    n_components::Union{Nothing,Int} = nothing,
    sample_selection = samples(),
    interval_selection::Interval = times(),
    channel_selection = channels(),
    include_extra::Bool = false,
    percentage_of_data::Real = 100.0,
    algorithm::Symbol = :infomax,
    use_gpu::Bool = false,
    params::IcaPrms = IcaPrms(),
    rng::AbstractRNG = Random.default_rng(),
)
    prms = deepcopy(params)
    prms.use_gpu = use_gpu || prms.use_gpu

    selected_channels = get_selected_channels(dat, channel_selection; include_meta = false, include_extra = include_extra)
    isempty(selected_channels) && @minimal_error("No channels available after applying channel filter")

    combined_sel = _combine_interval_sample(interval_selection, sample_selection)
    sample_indices = get_selected_samples(dat, combined_sel)
    isempty(sample_indices) && @minimal_error("No samples available after applying sample filter")

    ica_layout = subset_layout(dat.layout, channel_selection = channels(selected_channels))
    dat_for_ica = _create_ica_data_matrix(dat.data, selected_channels, sample_indices)
    percentage_of_data != 100 && (dat_for_ica = _select_subsample!(dat_for_ica, percentage_of_data; rng = rng))

    n_components = _resolve_ica_n_components(dat_for_ica, n_components)
    @info "Running ICA (use_gpu=$(prms.use_gpu)): $(length(selected_channels)) channels x $(size(dat_for_ica, 2)) samples -> $(n_components) components"

    return _run_ica_algorithm(
        dat_for_ica,
        ica_layout,
        filename(dat);
        n_components = n_components,
        algorithm = algorithm,
        params = prms,
        rng = rng,
    )
end

function run_ica(
    epoched_data::AbstractVector{<:EpochData};
    n_components::Union{Nothing,Int} = nothing,
    sample_selection = samples(),
    interval_selection::Interval = times(),
    channel_selection = channels(),
    include_extra::Bool = false,
    remove_duplicates::Bool = true,
    percentage_of_data::Real = 100.0,
    algorithm::Symbol = :infomax,
    use_gpu::Bool = false,
    params::IcaPrms = IcaPrms(),
    rng::AbstractRNG = Random.default_rng(),
)
    prms = deepcopy(params)
    prms.use_gpu = use_gpu || prms.use_gpu
    isempty(epoched_data) && @minimal_error("Empty epoched_data vector provided")

    reference_epoch_data = epoched_data[1]
    for (i, epoch_data) in enumerate(epoched_data)
        if epoch_data.sample_rate != reference_epoch_data.sample_rate
            @minimal_error(
                "Inconsistent sample rates: EpochData $i has $(epoch_data.sample_rate) Hz, expected $(reference_epoch_data.sample_rate) Hz",
            )
        end
    end

    selected_channels = get_selected_channels(reference_epoch_data, channel_selection; include_meta = false, include_extra = include_extra)
    isempty(selected_channels) && @minimal_error("No channels available after applying channel filter")

    concatenated_df = all_data(epoched_data)
    _check_epoched_data_uniqueness!(concatenated_df; remove_duplicates = remove_duplicates)

    combined_sel = _combine_interval_sample(interval_selection, sample_selection)
    sample_indices = get_selected_samples(concatenated_df, combined_sel)
    isempty(sample_indices) && @minimal_error("No samples available after applying sample filter to epoched data")

    concatenated_matrix = _create_ica_data_matrix(concatenated_df, selected_channels, sample_indices)
    percentage_of_data != 100 && (concatenated_matrix = _select_subsample!(concatenated_matrix, percentage_of_data; rng = rng))

    n_components = _resolve_ica_n_components(concatenated_matrix, n_components)
    total_epochs = sum(n_epochs(epoch_data) for epoch_data in epoched_data)
    @info "Running ICA on concatenated epochs: $(length(selected_channels)) channels x $(size(concatenated_matrix, 2)) samples (from $total_epochs epochs) -> $n_components components"

    ica_layout = subset_layout(reference_epoch_data.layout, channel_selection = channels(selected_channels))
    return _run_ica_algorithm(
        concatenated_matrix,
        ica_layout,
        filename(reference_epoch_data);
        n_components = n_components,
        algorithm = algorithm,
        params = prms,
        rng = rng,
    )
end

run_ica(epoched_data::EpochData; kwargs...) = run_ica([epoched_data]; kwargs...)

function _check_epoched_data_uniqueness!(concatenated_df::DataFrame; remove_duplicates::Bool = false)
    if !hasproperty(concatenated_df, :sample)
        @debug "No sample column found - cannot check for duplicate samples"
        return nothing
    end

    n_original = nrow(concatenated_df)
    unique_samples = length(unique(concatenated_df.sample))
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
            unique!(concatenated_df, :sample)
            @info "Automatically removed $n_duplicates duplicate samples (using $unique_samples unique samples)."
        end
    end
    return nothing
end

function _estimate_ica_rank(dat_ica::Matrix{Float64})
    n_channels, n_samples = size(dat_ica)
    (n_channels == 0 || n_samples <= 1) && return 0
    dat = dat_ica .- mean(dat_ica, dims = 2)
    cov_mat = Matrix{Float64}(undef, n_channels, n_channels)
    mul!(cov_mat, dat, transpose(dat), 1.0 / n_samples, 0.0)
    scale = sqrt(norm(cov_mat))
    scale <= 0 && return 0

    if n_samples > 2 * n_channels
        cov_mat .*= (n_samples / (scale * scale))
        evals = eigvals!(Symmetric(cov_mat))
        isempty(evals) && return 0
        l_max = evals[end]
        l_max <= 0 && return 0
        # Use fast covariance eigenvalues unless an eigenvalue sits in the near-singular roundoff band
        ambiguous_lo = 1e-16 * l_max
        ambiguous_hi = 1e-10 * l_max
        if !any(ev -> ambiguous_lo < ev < ambiguous_hi, evals)
            s_max = sqrt(l_max)
            tol = max(1e-7 * s_max, max(n_channels, n_samples) * eps(Float64) * s_max, 1e-10)
            return count(ev -> ev > 0.0 && sqrt(ev) >= tol, evals)
        end
    end

    dat ./= scale
    S = svdvals!(dat)
    isempty(S) && return 0
    s_max = S[1]
    s_max <= 0 && return 0
    tol = max(1e-7 * s_max, max(n_channels, n_samples) * eps(Float64) * s_max, 1e-10)
    return count(>=(tol), S)
end

function _resolve_ica_n_components(dat_ica::Matrix{Float64}, n_components::Union{Nothing,Int})
    n_channels = size(dat_ica, 1)
    eff_rank = _estimate_ica_rank(dat_ica)
    eff_rank == 0 && @minimal_error("Data has effective numerical rank 0 (flat or constant channels); cannot run ICA.")

    if isnothing(n_components)
        eff_rank < n_channels && @info "Detected numerical rank $eff_rank (< $n_channels channels); defaulting n_components = $eff_rank."
        return eff_rank
    end

    n_components <= 0 && @minimal_error("n_components must be positive, got $n_components")

    if n_components > n_channels
        @minimal_warning "Requested $n_components components, which exceeds $n_channels selected channels (effective rank $eff_rank). Using $eff_rank components instead."
        return eff_rank
    elseif n_components > eff_rank
        @minimal_warning "Requested $n_components components, which exceeds effective numerical rank $eff_rank (out of $n_channels channels). Clamping n_components to $eff_rank."
        return eff_rank
    end

    return n_components
end

@inline function _extract_ica_row!(result::Matrix{Float64}, row_idx::Int, col::AbstractVector, samples)
    @inbounds @simd for j in eachindex(samples)
        result[row_idx, j] = Float64(col[samples[j]])
    end
    return nothing
end

function _create_ica_data_matrix(dat::DataFrame, channels, samples)
    existing_channels = intersect(propertynames(dat), channels)
    result = Matrix{Float64}(undef, length(existing_channels), length(samples))
    for (i, ch) in enumerate(existing_channels)
        _extract_ica_row!(result, i, dat[!, ch], samples)
    end
    return result
end

function _select_subsample!(data_matrix::Matrix{Float64}, percentage::Real; rng::AbstractRNG = Random.default_rng())
    (percentage <= 0 || percentage > 100) && @minimal_error("percentage_of_data must be between 0 and 100, got $percentage")

    n_original_samples = size(data_matrix, 2)
    n_target_samples = round(Int, n_original_samples * percentage / 100)

    if n_target_samples == n_original_samples
        @info "Random sample: $n_target_samples of $n_original_samples ($(round(percentage, digits=1))%)"
        return data_matrix
    end

    sample_cols = sort!(randperm(rng, n_original_samples)[1:n_target_samples])
    subsampled_matrix = data_matrix[:, sample_cols]
    @info "Random sample: $n_target_samples of $n_original_samples ($(round(percentage, digits=1))%)"
    return subsampled_matrix
end

struct WorkArrays
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

function create_work_arrays(n_components::Int, block_size::Int)
    weights = Matrix{Float64}(I, n_components, n_components)
    return WorkArrays(
        weights,
        zeros(n_components, block_size),
        zeros(n_components, block_size),
        zeros(n_components, block_size),
        copy(weights),
        copy(weights),
        zeros(n_components, n_components),
        zeros(n_components, n_components),
        zeros(n_components, n_components),
        zeros(n_components, n_components),
    )
end

function _run_ica_algorithm(
    dat_ica::Matrix{Float64},
    layout::Layout,
    filename::String;
    n_components::Int,
    algorithm::Symbol = :infomax,
    params::IcaPrms = IcaPrms(),
    rng::AbstractRNG = Random.default_rng(),
)
    if algorithm === :infomax
        return infomax_ica(dat_ica, layout, filename; n_components = n_components, params = params, rng = rng)
    elseif algorithm === :infomax_extended
        return infomax_extended_ica(dat_ica, layout, filename; n_components = n_components, params = params, rng = rng)
    elseif algorithm === :picard || algorithm === :picard_extended
        return picard_ica(
            dat_ica,
            layout,
            filename;
            n_components = n_components,
            extended = (algorithm === :picard_extended),
            params = params,
            rng = rng,
        )
    elseif algorithm === :amica
        return amica_ica(dat_ica, layout, filename; n_components = n_components, params = params, rng = rng)
    else
        @minimal_error(
            "Unknown ICA algorithm: $algorithm. Supported algorithms: :infomax, :infomax_extended, :picard, :picard_extended, :amica",
        )
    end
end
