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
    params.use_gpu = use_gpu || params.use_gpu
    selected_channels = get_selected_channels(dat, channel_selection; include_meta = false, include_extra = include_extra)
    isempty(selected_channels) && @minimal_error("No channels available after applying channel filter")

    # Combine interval and sample selection 
    combined_sel = _combine_interval_sample(interval_selection, sample_selection)

    # Get samples to use using predicate
    sample_indices = get_selected_samples(dat, combined_sel)
    isempty(sample_indices) && @minimal_error("No samples available after applying sample filter")

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
        dat_for_ica = _select_subsample!(dat_for_ica, percentage_of_data; rng = rng)
    end

    # Get filename from the data
    input_filename = filename(dat)

    # Dispatch to the appropriate ICA algorithm
    ica_result = _run_ica_algorithm(
        dat_for_ica,
        ica_layout,
        input_filename;
        n_components = n_components,
        algorithm = algorithm,
        params = params,
        rng = rng,
    )

    return ica_result
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
    params.use_gpu = use_gpu || params.use_gpu
    isempty(epoched_data) && @minimal_error("Empty epoched_data vector provided")

    # Use the first EpochData object as reference for some meta-like data
    reference_epoch_data = epoched_data[1]
    for (i, epoch_data) in enumerate(epoched_data)
        if epoch_data.sample_rate != reference_epoch_data.sample_rate
            @minimal_error(
                "Inconsistent sample rates: EpochData $i has $(epoch_data.sample_rate) Hz, expected $(reference_epoch_data.sample_rate) Hz",
            )
        end
    end

    # Get channel information from reference
    selected_channels = get_selected_channels(reference_epoch_data, channel_selection; include_meta = false, include_extra = include_extra)
    isempty(selected_channels) && @minimal_error("No channels available after applying channel filter")

    # Concatenate all epoched data and check for duplicates
    concatenated_df = all_data(epoched_data)
    _check_epoched_data_uniqueness!(concatenated_df; remove_duplicates = remove_duplicates)

    # Combine interval and sample selection (consistent with subset() pattern)
    combined_sel = _combine_interval_sample(interval_selection, sample_selection)

    sample_indices = get_selected_samples(concatenated_df, combined_sel)
    isempty(sample_indices) && @minimal_error("No samples available after applying sample filter to epoched data")

    # Create data matrix for ICA
    concatenated_matrix = _create_ica_data_matrix(concatenated_df, selected_channels, sample_indices)
    if percentage_of_data != 100
        concatenated_matrix = _select_subsample!(concatenated_matrix, percentage_of_data; rng = rng)
    end

    # Set n_components if not specified
    if isnothing(n_components)
        n_components = length(selected_channels) - 1
    elseif n_components > length(selected_channels)
        @minimal_warning "Requested $n_components components but only $(length(selected_channels)) channels available. Using $(length(selected_channels) - 1) components instead."
        n_components = length(selected_channels) - 1
    end

    final_samples = size(concatenated_matrix, 2)
    total_epochs = sum(n_epochs(epoch_data) for epoch_data in epoched_data)

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
        rng = rng,
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
function _select_subsample!(data_matrix::Matrix{Float64}, percentage::Real; rng::AbstractRNG = Random.default_rng())
    (percentage <= 0 || percentage > 100) && @minimal_error("percentage_of_data must be between 0 and 100, got $percentage")

    n_original_samples = size(data_matrix, 2)
    n_target_samples = round(Int, n_original_samples * percentage / 100)

    if n_target_samples == n_original_samples
        @info "Random sample: $n_target_samples of $n_original_samples ($(round(percentage, digits=1))%)"
        return data_matrix
    end

    # Random sample selection (columns are samples in data_matrix)
    sample_cols = randperm(rng, n_original_samples)[1:n_target_samples]
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
    rng::AbstractRNG = Random.default_rng(),
)
    if algorithm == :infomax
        return infomax_ica(dat_ica, layout, filename; n_components = n_components, params = params, rng = rng)
    elseif algorithm == :infomax_extended
        return infomax_extended_ica(dat_ica, layout, filename; n_components = n_components, params = params, rng = rng)
    elseif algorithm == :picard || algorithm == :picard_extended
        extended = algorithm == :picard_extended
        return picard_ica(dat_ica, layout, filename; n_components = n_components, extended = extended, params = params, rng = rng)
    elseif algorithm == :amica
        return amica_ica(dat_ica, layout, filename; n_components = n_components, params = params, rng = rng)
    else
        @minimal_error(
            "Unknown ICA algorithm: $algorithm. Supported algorithms: :infomax, :infomax_extended, :picard, :picard_extended, :amica"
        )
    end
end
