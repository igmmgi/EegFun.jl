"""
    run_ica(dat::ContinuousData;
            n_components::Union{Nothing,Int} = nothing,
            sample_selection::Function = samples(),
            channel_selection::Function = channels(),
            include_extra::Bool = false,
            percentage_of_data::Real = 100.0,
            algorithm::Symbol = :infomax,
            params::IcaPrms = IcaPrms())

Runs Independent Component Analysis (ICA) on EEG data. Preprocessing (e.g., filtering) should be applied prior to calling this function.

# Arguments
- `dat::ContinuousData`: The EEG data object.
- `n_components::Union{Nothing,Int}`: Number of ICA components (default: number of channels - 1).
- `sample_selection::Function`: Sample selector for quality filtering (default: include all samples).
- `channel_selection::Function`: Channel selector (default: layout channels).
- `include_extra::Bool`: Whether to allow channels outside the layout (e.g., EOG) in selection.
- `percentage_of_data::Real`: Percentage of good data to use for ICA (default: 100.0). Values < 100 enable faster computation by random subsampling.
- `algorithm::Symbol`: ICA algorithm to use. Options: `:infomax` (default), `:fastica`, `:sobi`, `:jade`. More algorithms coming soon.
- `params::IcaPrms`: ICA parameters (algorithm-specific).

# Returns
`InfoIca` with unmixing, mixing, sphere, variance, and metadata.

# Examples
```julia
# Basic ICA on layout channels (Infomax algorithm)
ica_result = run_ica(dat)

# Use FastICA algorithm (often faster)
ica_result = run_ica(dat, algorithm = :fastica)

# Excluding extreme samples (quality filtering)
ica_result = run_ica(dat, sample_selection = samples_not(:is_extreme_value_100))

# Speed up ICA by using only 25% of good data (4x faster)
ica_result = run_ica(dat, sample_selection = samples_not(:is_extreme_value_100), percentage_of_data = 25.0)

# Restrict to specific channels
ica_result = run_ica(dat, channel_selection = channels([:Fp1, :Fp2, :F3, :F4]))

# Include EOG channels explicitly (use with caution)
ica_result = run_ica(dat, include_extra = true, channel_selection = channels([:Fp1, :Fp2, :vEOG, :hEOG]))
```
"""
function run_ica(
    dat::ContinuousData;
    n_components::Union{Nothing,Int} = nothing,
    sample_selection::Function = samples(),
    channel_selection::Function = channels(),
    include_extra::Bool = false,
    percentage_of_data::Real = 100.0,
    algorithm::Symbol = :infomax,
    params::IcaPrms = IcaPrms(),
)
    # Create a copy of the data to avoid modifying the original
    dat_ica = copy(dat)

    selected_channels =
        get_selected_channels(dat_ica, channel_selection; include_meta = false, include_extra = include_extra)
    if isempty(selected_channels)
        error("No channels available after applying channel filter")
    end

    # Get samples to use using predicate
    sample_indices = get_selected_samples(dat_ica, sample_selection)
    if isempty(sample_indices)
        error("No samples available after applying sample filter")
    end

    # Set n_components if not specified
    if isnothing(n_components)
        n_components = length(selected_channels) - 1
    elseif n_components > length(selected_channels)
        @warn "Requested $n_components components but only $(length(selected_channels)) channels available. Using $(length(selected_channels) - 1) components instead."
        n_components = length(selected_channels) - 1
    end

    @info "Running ICA: $(length(selected_channels)) channels x $(length(sample_indices)) samples -> $(n_components) components"

    # Create subsetted layout that matches the selected channels
    ica_layout = subset_layout(dat_ica.layout, channel_selection = channels(selected_channels))

    # Create data matrix and run ICA
    dat_for_ica = create_ica_data_matrix(dat_ica.data, selected_channels, sample_indices)
    if percentage_of_data !== 100
        dat_for_ica = _select_subsample!(dat_for_ica, percentage_of_data)
    end

    # Dispatch to the appropriate ICA algorithm
    ica_result = _run_ica_algorithm(dat_for_ica, ica_layout, n_components = n_components, algorithm = algorithm, params = params)

    return ica_result
end

"""
    run_ica(epoched_data::Vector{EpochData};
            n_components::Union{Nothing,Int} = nothing,
            sample_selection::Function = samples(),
            channel_selection::Function = channels(),
            include_extra::Bool = false,
            remove_duplicates::Bool = true,
            percentage_of_data::Real = 100.0,
            algorithm::Symbol = :infomax,
            params::IcaPrms = IcaPrms())

Runs Independent Component Analysis (ICA) on concatenated epoched EEG data.

This method concatenates all epochs from all EpochData objects into a continuous 
data matrix before running ICA. This is the standard approach since ICA benefits 
from more data points and temporal structure is not critical for decomposition.

Automatically checks for duplicate samples (same original data appearing in multiple epochs)
and removes them by default, as duplicates can bias ICA decomposition.

# Arguments
- `epoched_data::Vector{EpochData}`: Vector of epoched EEG data objects to concatenate
- `n_components::Union{Nothing,Int}`: Number of ICA components (default: number of channels - 1)
- `sample_selection::Function`: Sample selector for quality filtering applied to each epoch (default: include all samples)
- `channel_selection::Function`: Channel selector (default: layout channels)
- `include_extra::Bool`: Whether to allow channels outside the layout (e.g., EOG) in selection
- `remove_duplicates::Bool`: Automatically remove duplicate samples based on samples column (default: true)
- `percentage_of_data::Real`: Percentage of good data to use for ICA (default: 100.0). Values < 100 enable faster computation by random subsampling.
- `algorithm::Symbol`: ICA algorithm to use. Options: `:infomax` (default), `:fastica`, `:sobi`, `:jade`. More algorithms coming soon.
- `params::IcaPrms`: ICA parameters (algorithm-specific)

# Returns
`InfoIca` with unmixing, mixing, sphere, variance, and metadata.

# Examples
```julia
# Basic ICA on concatenated epochs (Infomax algorithm)
epochs = extract_epochs(dat, epoch_conditions, -1, 2)
ica_result = run_ica(epochs)

# Use FastICA algorithm (often faster)
ica_result = run_ica(epochs, algorithm = :fastica)

# Excluding artifact samples from each epoch
ica_result = run_ica(epochs, sample_selection = samples_not(:is_extreme_value_100))

# Keep duplicate samples (not recommended for ICA)
ica_result = run_ica(epochs, remove_duplicates = false)
```
"""
function run_ica(
    epoched_data::Vector{EpochData};
    n_components::Union{Nothing,Int} = nothing,
    sample_selection::Function = samples(),
    channel_selection::Function = channels(),
    include_extra::Bool = false,
    remove_duplicates::Bool = true,
    percentage_of_data::Real = 100.0,
    algorithm::Symbol = :infomax,
    params::IcaPrms = IcaPrms(),
)
    if isempty(epoched_data)
        error("Empty epoched_data vector provided")
    end

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
    selected_channels = get_selected_channels(
        reference_epoch_data,
        channel_selection;
        include_meta = false,
        include_extra = include_extra,
    )
    if isempty(selected_channels)
        error("No channels available after applying channel filter")
    end

    # Concatenate all epoched data and check for duplicates
    concatenated_df = all_data(epoched_data)
    _check_epoched_data_uniqueness!(concatenated_df; remove_duplicates = remove_duplicates)

    sample_indices = get_selected_samples(concatenated_df, sample_selection)
    if isempty(sample_indices)
        error("No samples available after applying sample filter to epoched data")
    end

    # Create data matrix for ICA
    concatenated_matrix = create_ica_data_matrix(concatenated_df, selected_channels, sample_indices)
    if percentage_of_data !== 100
        concatenated_matrix = _select_subsample!(concatenated_matrix, percentage_of_data)
    end

    # Set n_components if not specified
    if isnothing(n_components)
        n_components = length(selected_channels) - 1
    elseif n_components > length(selected_channels)
        @warn "Requested $n_components components but only $(length(selected_channels)) channels available. Using $(length(selected_channels) - 1) components instead."
        n_components = length(selected_channels) - 1
    end

    final_samples = size(concatenated_matrix, 2)
    total_epochs = sum(length(epoch_data.data) for epoch_data in epoched_data)

    @info "Running ICA on concatenated epochs: $(length(selected_channels)) channels x $final_samples samples (from $total_epochs epochs) -> $n_components components"

    # Create subsetted layout that matches the selected channels  
    ica_layout = subset_layout(reference_epoch_data.layout, channel_selection = channels(selected_channels))

    # Run ICA on concatenated data
    # Dispatch to the appropriate ICA algorithm
    ica_result = _run_ica_algorithm(concatenated_matrix, ica_layout, n_components = n_components, algorithm = algorithm, params = params)

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

        @warn """
        Found $n_duplicates duplicate samples ($duplicate_percentage%) in concatenated epoched data.
        This occurs when epochs have overlapping time windows - the same original data samples 
        appear multiple times. Duplicate samples may bias ICA decomposition.

        Total rows: $n_original
        Unique samples: $unique_samples  
        Duplicates: $n_duplicates ($duplicate_percentage%)
        """

        if remove_duplicates
            unique!(concatenated_df, :samples)
            @info "Automatically removed $n_duplicates duplicate samples based on samples column. Using $unique_samples unique samples for ICA."
        end
    end

    return nothing
end

function run_ica(epoched_data::EpochData; kwargs...)
    run_ica([epoched_data]; kwargs...)
end

function IcaPrms(;
    l_rate = 0.001,
    max_iter = 512,
    w_change = 1e-6,
    anneal_deg = 60.0,
    anneal_step = 0.9,
    blowup = 1e15,
    blowup_fac = 0.8,
    max_weight = 1e8,
    restart_factor = 0.9,
    degconst = 180.0 / π,
    default_stop = 1e-6,
)
    IcaPrms(
        l_rate,
        max_iter,
        w_change,
        anneal_deg,
        anneal_step,
        blowup,
        blowup_fac,
        max_weight,
        restart_factor,
        degconst,
        default_stop,
    )
end

function create_ica_data_matrix(dat::DataFrame, channels, samples)
    # Filter to only existing channels (matching original behavior)
    existing_channels = intersect(propertynames(dat), channels)
    n_channels = length(existing_channels)
    n_samples = length(samples)

    # Pre-allocate result matrix
    result = Matrix{Float64}(undef, n_channels, n_samples)

    # Use direct column access for better performance
    for (i, ch) in enumerate(existing_channels)
        result[i, :] = dat[samples, ch]
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
    if percentage <= 0 || percentage > 100
        error("percentage_of_data must be between 0 and 100, got $percentage")
    end

    original_samples = size(data_matrix, 2)
    target_samples = round(Int, original_samples * percentage / 100)

    # Random sample selection (columns are samples in data_matrix)
    sample_cols = randperm(original_samples)[1:target_samples]
    subsampled_matrix = data_matrix[:, sample_cols]

    @info "Random subsampling: using $target_samples of $original_samples samples ($(round(percentage, digits=1))%) for $(round(100/percentage, digits=1))x speedup"

    return subsampled_matrix
end


mutable struct WorkArrays
    weights::Matrix{Float64}
    BI::Matrix{Float64}
    u::Matrix{Float64}
    y::Matrix{Float64}
    data_block::Matrix{Float64}
    oldweights::Matrix{Float64}
    startweights::Matrix{Float64}
    weights_temp::Matrix{Float64}
    bi_weights::Matrix{Float64}
    wu_term::Matrix{Float64}
    delta::Matrix{Float64}
    olddelta::Matrix{Float64}
end

function create_work_arrays(n_components::Int, block_size::Int)
    weights = Matrix{Float64}(I, n_components, n_components)  # Initialize as identity matrix
    return WorkArrays(
        weights,  # weights - start as identity matrix
        block_size * Matrix{Float64}(I, n_components, n_components),  # BI
        zeros(n_components, block_size),  # u
        zeros(n_components, block_size),  # y
        zeros(n_components, block_size),  # data_block
        copy(weights),  # oldweights - copy of identity matrix
        copy(weights),  # startweights - copy of identity matrix
        zeros(n_components, n_components),  # weights_temp
        zeros(n_components, n_components),  # bi_weights
        zeros(n_components, n_components),  # wu_term
        zeros(1, n_components^2),  # delta
        zeros(1, n_components^2),   # olddelta
    )
end

# =============================================================================
# ICA ALGORITHM DISPATCHER
# =============================================================================

"""
    _run_ica_algorithm(dat_ica, layout; n_components, algorithm, params)

Internal dispatcher function that routes to the appropriate ICA algorithm implementation.

# Arguments
- `dat_ica::Matrix{Float64}`: Data matrix (channels × samples)
- `layout::Layout`: Layout information
- `n_components::Int`: Number of ICA components
- `algorithm::Symbol`: Algorithm to use (`:infomax`, `:fastica`, `:sobi`, `:jade`)
- `params::IcaPrms`: Algorithm-specific parameters

# Returns
`InfoIca` result structure
"""
function _run_ica_algorithm(
    dat_ica::Matrix{Float64},
    layout::Layout;
    n_components::Int,
    algorithm::Symbol = :infomax,
    params::IcaPrms = IcaPrms(),
)
    if algorithm == :infomax
        return infomax_ica(dat_ica, layout, n_components = n_components, params = params)
    elseif algorithm == :fastica
        return fastica_ica(dat_ica, layout, n_components = n_components, params = params)
    elseif algorithm == :sobi
        error("SOBI not yet implemented. Coming soon!")
    elseif algorithm == :jade
        error("JADE not yet implemented. Coming soon!")
    elseif algorithm == :amica
        error("AMICA not yet implemented. Coming soon!")
    else
        error("Unknown ICA algorithm: $algorithm. Supported algorithms: :infomax, :fastica, :sobi, :jade, :amica")
    end
end

# =============================================================================
# INFOMAX ICA IMPLEMENTATION
# =============================================================================

function infomax_ica(dat_ica::Matrix{Float64}, layout::Layout; n_components::Int, params::IcaPrms = IcaPrms())

    # Store original mean before removing it
    original_mean = vec(mean(dat_ica, dims = 2))

    # Center and scale data
    dat_ica .-= original_mean
    scale = sqrt(norm((dat_ica * dat_ica') / size(dat_ica, 2)))
    dat_ica ./= scale

    # PCA reduction - optimized for speed
    n_channels, n_samples = size(dat_ica)
    F = svd(dat_ica)
    pca_components = F.U[:, 1:n_components]

    # PCA projection into workspace
    workspace = Matrix{Float64}(undef, n_components, n_samples)
    mul!(workspace, pca_components', dat_ica)

    # Sphering: reuse original dat_ica memory (resize to smaller dimensions)
    sphere = inv(sqrt(cov(workspace, dims = 2)))
    dat_ica = Matrix{Float64}(undef, n_components, n_samples)  # Resize to final dimensions
    mul!(dat_ica, sphere, workspace)

    # initialize
    n_channels = size(dat_ica, 1)
    n_samples = size(dat_ica, 2)
    # Keep original Infomax block size formula for algorithmic correctness
    block = min(Int(floor(sqrt(n_samples / 3.0))), 512)
    work = create_work_arrays(n_channels, block)

    step = 0
    wts_blowup = false
    change = 0.0
    oldchange = 0.0
    angledelta = 0.0

    # pre-allocate permutation vector
    permute_indices = Vector{Int}(undef, n_samples)

    @inbounds while step < params.max_iter
        randperm!(permute_indices)

        for t = 1:block:n_samples
            block_end = min(t + block - 1, n_samples)
            block_size = block_end - t + 1

            # extract data block
            copyto!(view(work.data_block, :, 1:block_size), view(dat_ica, :, view(permute_indices, t:block_end)))

            # forward pass
            mul!(work.u, work.weights, work.data_block)
            @. work.y = 1 - 2 / (1 + exp(-work.u))

            # update weights 
            mul!(work.wu_term, work.y, transpose(work.u))
            work.bi_weights .= work.BI .+ work.wu_term
            mul!(work.weights_temp, work.bi_weights, work.weights)
            @. work.weights += params.l_rate * work.weights_temp

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
            work.delta .= reshape(work.oldweights, 1, :)
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

        @info "Step $step, change = $change, lrate = $(params.l_rate), angle = $((params.degconst) * angledelta)"

    end

    # Final calculations
    work.weights = work.weights * sphere * pca_components'
    mixing = pinv(work.weights)

    # calculate total variance explained and order
    meanvar = vec(sum(abs2, mixing, dims = 1) .* sum(abs2, dat_ica, dims = 2)' ./ (n_components * n_samples - 1))
    meanvar_normalized = meanvar ./ sum(meanvar)
    order = sortperm(meanvar_normalized, rev = true)

    return InfoIca(
        work.weights[order, :],
        mixing[:, order],
        sphere,
        meanvar_normalized[order],
        scale,
        original_mean,
        [Symbol("IC$i") for i = 1:size(work.weights, 1)],
        OrderedDict{Int,Matrix{Float64}}(),
        layout,
    )

end

# =============================================================================
# FASTICA IMPLEMENTATION
# =============================================================================

"""
    _invsqrtm!(C::AbstractMatrix{<:Real})

Compute inv(sqrtm(C)) through symmetric eigenvalue decomposition.
In-place version that modifies C.

This is used for symmetric decorrelation in FastICA: W = W * (W'W)^{-1/2}
"""
function _invsqrtm!(C::AbstractMatrix{<:Real})
    n = size(C, 1)
    size(C, 2) == n || error("C must be a square matrix.")
    E = eigen!(Symmetric(C))
    U = E.vectors
    evs = E.values
    for i = 1:n
        @inbounds evs[i] = 1.0 / sqrt(sqrt(evs[i]))
    end
    rmul!(U, Diagonal(evs))
    return U * transpose(U)
end

"""
    fastica_ica(dat_ica::Matrix{Float64}, layout::Layout; n_components::Int, params::IcaPrms = IcaPrms())

FastICA algorithm implementation using fixed-point iteration.

FastICA is often faster than Infomax and uses a fixed-point iteration scheme
to maximize non-Gaussianity. This implementation uses the symmetric approach
(extracting all components simultaneously), which is faster and more stable than
the deflationary approach.

# Arguments
- `dat_ica::Matrix{Float64}`: Data matrix (channels × samples), should be preprocessed
- `layout::Layout`: Layout information for channels
- `n_components::Int`: Number of ICA components to extract
- `params::IcaPrms`: ICA parameters (uses max_iter, w_change from params)

# Returns
`InfoIca` with unmixing, mixing, sphere, variance, and metadata.
"""
function fastica_ica(dat_ica::Matrix{Float64}, layout::Layout; n_components::Int, params::IcaPrms = IcaPrms())
    
    # Store original mean before removing it
    original_mean = vec(mean(dat_ica, dims = 2))
    
    # Center and scale data
    dat_ica .-= original_mean
    scale = sqrt(norm((dat_ica * dat_ica') / size(dat_ica, 2)))
    dat_ica ./= scale
    
    # PCA reduction
    n_channels, n_samples = size(dat_ica)
    F = svd(dat_ica)
    pca_components = F.U[:, 1:n_components]
    
    # PCA projection into workspace
    workspace = Matrix{Float64}(undef, n_components, n_samples)
    mul!(workspace, pca_components', dat_ica)
    
    # Sphering: reuse original dat_ica memory (resize to smaller dimensions)
    # Match Infomax preprocessing exactly
    sphere = inv(sqrt(cov(workspace, dims = 2)))
    dat_ica = Matrix{Float64}(undef, n_components, n_samples)  # Resize to final dimensions
    mul!(dat_ica, sphere, workspace)
    
    # Data should already be centered after sphering, but ensure it
    # (sklearn doesn't re-center after whitening, so we match that)
    
    # Initialize unmixing matrix (orthogonal random matrix)
    n_channels = size(dat_ica, 1)
    n_samples = size(dat_ica, 2)
    
    # Initialize weights as random orthogonal matrix (n_channels × n_components)
    # Note: W is transposed compared to deflationary approach for symmetric FastICA
    W = Matrix{Float64}(undef, n_channels, n_components)
    randn!(W)
    # Normalize each column
    for j = 1:n_components
        w = view(W, :, j)
        rmul!(w, 1.0 / sqrt(sum(abs2, w)))
    end
    
    # FastICA parameters
    max_iter = params.max_iter
    # Use tolerance matching sklearn's default (1e-4) for FastICA
    # sklearn's FastICA uses tol=1e-4 by default, which is more lenient than Infomax
    # This is appropriate because symmetric FastICA converges differently than deflationary
    tol = max(params.w_change, 1e-4)  # At least 1e-4 like sklearn, but respect user's w_change if higher
    contrast_function = :tanh  # Options: :tanh, :gauss, :pow3
    
    # Vectorized contrast function computation (like MultiVariateStats)
    # For tanh: g(x) = tanh(x), g'(x) = 1 - tanh(x)^2
    function update_tanh!(U::AbstractMatrix{Float64}, E::AbstractVector{Float64})
        n, k = size(U)
        @inbounds for j in 1:k
            _s = zero(Float64)
            @fastmath for i in 1:n
                t = tanh(U[i, j])
                U[i, j] = t
                _s += 1.0 - t^2
            end
            E[j] = _s / n
        end
    end
    
    # Pre-allocated storage for symmetric FastICA (like MultiVariateStats)
    Wp = similar(W)                # previous version of W
    U  = Matrix{Float64}(undef, n_samples, n_components)  # to store w'x & g(w'x)
    Y  = Matrix{Float64}(undef, n_channels, n_components)  # to store E{x g(w'x)} for components
    E1 = Vector{Float64}(undef, n_components)              # store E{g'(w'x)} for components
    
    # Main symmetric FastICA loop (extract all components simultaneously)
    chg = Float64(NaN)
    converged = false
    
    for iteration in 1:max_iter
        copyto!(Wp, W)
        
        # Apply W of previous step: U = X' * W (all components at once)
        # U[i, j] = w_j' * x_i for sample i and component j
        mul!(U, transpose(dat_ica), W)
        
        # Compute g(w'x) --> U and E{g'(w'x)} --> E1 (vectorized for all components)
        if contrast_function == :tanh
            update_tanh!(U, E1)
        else
            # Fallback to element-wise for other contrast functions
            @inbounds for j in 1:n_components
                _s = zero(Float64)
                for i in 1:n_samples
                    u = U[i, j]
                    if contrast_function == :gauss
                        u2 = u^2
                        e = exp(-u2 / 2)
                        U[i, j] = u * e
                        _s += (1 - u2) * e
                    elseif contrast_function == :pow3
                        U[i, j] = u^3
                        _s += 3.0 * u^2
                    else  # default to tanh
                        t = tanh(u)
                        U[i, j] = t
                        _s += 1.0 - t^2
                    end
                end
                E1[j] = _s / n_samples
            end
        end
        
        # Compute E{x g(w'x)} --> Y (all components at once)
        # Y[:, j] = mean(X * U[:, j]) for component j
        rmul!(mul!(Y, dat_ica, U), 1.0 / n_samples)
        
        # Update all components: w := y - e1 * w
        # Following MultiVariateStats: w = y - e1 * w
        for j = 1:n_components
            w = view(W, :, j)
            y = view(Y, :, j)
            e1 = E1[j]
            @. w = y - e1 * w
        end
        
        # Symmetric decorrelation: W <- W * (W'W)^{-1/2}
        # This is more numerically stable than Gram-Schmidt
        # Following MultiVariateStats exactly: copyto!(W, W * _invsqrtm!(W'W))
        # Note: _invsqrtm! modifies its input, so we compute W'W inline
        copyto!(W, W * _invsqrtm!(W'W))
        
        # Compare with Wp to evaluate convergence change
        # We want to check if each component (column of W) has converged
        # W is m×k (channels × components), so W'*Wp is k×k (components × components)
        # The diagonal gives the dot product of each component with its previous version
        # Following the FastICA paper: we check max(|abs(diag(W'*Wp)) - 1|)
        # This measures how much each component has rotated (should be close to 1 if converged)
        chg = maximum(abs.(abs.(diag(W' * Wp)) .- 1))
        converged = (chg < tol)
        
        if iteration % 10 == 0 || converged
            @info "FastICA iteration $iteration, change = $chg, tolerance = $tol"
        end
        
        if converged
            break
        end
    end
    
    if !converged
        @warn "FastICA did not converge after $max_iter iterations. Final change = $chg (tolerance = $tol)"
    end
    
    # Transpose W to get n_components × n_channels (matching expected format)
    weights = transpose(W)
    
    # Final calculations - convert to original space
    # weights is the unmixing matrix in sphered space
    # Transform back: unmixing = weights * sphere * pca_components' (same as Infomax)
    unmixing = weights * sphere * pca_components'
    mixing = pinv(unmixing)
    
    # Calculate variance explained and order
    # Use sphered dat_ica for variance calculation (matching Infomax exactly)
    # This ensures consistent ordering between Infomax and FastICA
    meanvar = vec(sum(abs2, mixing, dims = 1) .* sum(abs2, dat_ica, dims = 2)' ./ (n_components * n_samples - 1))
    meanvar_normalized = meanvar ./ sum(meanvar)
    order = sortperm(meanvar_normalized, rev = true)
    
    return InfoIca(
        unmixing[order, :],
        mixing[:, order],
        sphere,
        meanvar_normalized[order],
        scale,
        original_mean,
        [Symbol("IC$i") for i = 1:size(unmixing, 1)],
        OrderedDict{Int,Matrix{Float64}}(),
        layout,
    )
end


"""
    remove_ica_components!(dat::DataFrame, ica::InfoIca, components_to_remove::Vector{Int})

Remove ICA components from data in-place and store the removed activations in the ICA result.

# Arguments
- `dat::DataFrame`: DataFrame containing the data to clean
- `ica::InfoIca`: ICA result object (will be mutated to store removed activations)
- `components_to_remove::Vector{Int}`: Vector of component indices to remove

# Returns
- `DataFrame`: The mutated data (same as input dat)

# Example
```julia
dat_cleaned = remove_ica_components!(dat.data, ica_result, [1, 3, 5])
# Removed activations are now stored in ica_result.removed_activations
```
"""
function remove_ica_components!(dat::DataFrame, ica::InfoIca; component_selection::Function = components())
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
    tra =
        Matrix(I, n_channels, n_channels) -
        view(ica.mixing, :, components_to_remove) * view(ica.unmixing, components_to_remove, :)

    # Apply transformation and restore scaling
    cleaned_data = tra * data
    cleaned_data .*= ica.scale
    cleaned_data .+= ica.mean

    # Update DataFrame in-place
    dat[!, ica.layout.data.label] .= permutedims(cleaned_data)

    return dat
end

"""
    remove_ica_components(dat::DataFrame, ica::InfoIca, components_to_remove::Vector{Int})

Remove ICA components from data and return cleaned data and updated ICA result.

# Arguments
- `dat::DataFrame`: DataFrame containing the data to clean
- `ica::InfoIca`: ICA result object
- `components_to_remove::Vector{Int}`: Vector of component indices to remove

# Returns
- `Tuple{DataFrame, InfoIca}`: Copy of input data with components removed, and copy of ICA result with stored activations

# Example
```julia
cleaned_data, ica_updated = remove_ica_components(dat.data, ica_result, [1, 3, 5])
```
"""
function remove_ica_components(dat::DataFrame, ica::InfoIca; component_selection::Function = components())
    dat_out = copy(dat)
    ica_out = copy(ica)  # Use our custom copy method
    remove_ica_components!(dat_out, ica_out, component_selection = component_selection)
    return dat_out, ica_out
end

"""
    remove_ica_components(dat::ContinuousData, ica::InfoIca, components_to_remove::Vector{Int})

Remove ICA components from ContinuousData and return cleaned data and updated ICA result.

# Arguments
- `dat::ContinuousData`: ContinuousData object containing the data to clean
- `ica::InfoIca`: ICA result object
- `components_to_remove::Vector{Int}`: Vector of component indices to remove

# Returns
- `Tuple{ContinuousData, InfoIca}`: Copy of input data with components removed, and copy of ICA result with stored activations

# Example
```julia
cleaned_dat, ica_updated = remove_ica_components(dat, ica_result, [1, 3, 5])
```
"""
function remove_ica_components(dat::ContinuousData, ica::InfoIca; component_selection::Function = components())
    dat_out = copy(dat)
    ica_out = copy(ica)  # Use our custom copy method
    remove_ica_components!(dat_out.data, ica_out, component_selection = component_selection)
    return dat_out, ica_out
end

"""
    remove_ica_components!(dat::ContinuousData, ica::InfoIca, components_to_remove::Vector{Int})

Remove ICA components from ContinuousData in-place and store the removed activations in the ICA result.

# Arguments
- `dat::ContinuousData`: ContinuousData object to clean in-place
- `ica::InfoIca`: ICA result object (will be mutated to store removed activations)
- `components_to_remove::Vector{Int}`: Vector of component indices to remove

# Returns
- `ContinuousData`: The mutated data (same as input dat)

# Example
```julia
dat_cleaned = remove_ica_components!(dat, ica_result, [1, 3, 5])
# Removed activations are now stored in ica_result.removed_activations
```
"""
function remove_ica_components!(dat::ContinuousData, ica::InfoIca; component_selection::Function = components())
    remove_ica_components!(dat.data, ica, component_selection = component_selection)
    return dat
end


"""
    restore_ica_components!(dat::DataFrame, ica::InfoIca, components_to_restore::Vector{Int})

Restore ICA components to data in-place using stored activations and update ICA result.

# Arguments
- `dat::DataFrame`: DataFrame containing the data to restore
- `ica::InfoIca`: ICA result object with stored activations (will be updated)
- `components_to_restore::Vector{Int}`: Vector of component indices to restore

# Returns
- `nothing` (data and ICA result are modified in-place)

# Example
```julia
restore_ica_components!(dat.data, ica_result, [1, 3, 5])
```
"""
function restore_ica_components!(dat::DataFrame, ica::InfoIca; component_selection::Function = components())
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

"""
    restore_ica_components(dat::DataFrame, ica::InfoIca, components_to_restore::Vector{Int})

Restore ICA components to data and return restored data and updated ICA result.

# Arguments
- `dat::DataFrame`: DataFrame containing the data to restore
- `ica::InfoIca`: ICA result object with stored activations
- `components_to_restore::Vector{Int}`: Vector of component indices to restore

# Returns
- `Tuple{DataFrame, InfoIca}`: Copy of input data with components restored, and copy of ICA result with updated removed_activations

# Example
```julia
restored_data, ica_updated = restore_ica_components(dat.data, ica_result, [1, 3, 5])
```
"""
function restore_ica_components(dat::DataFrame, ica::InfoIca; component_selection::Function = components())
    dat_out = copy(dat)
    ica_out = copy(ica)  # Use our custom copy method
    restore_ica_components!(dat_out, ica_out, component_selection = component_selection)
    return dat_out, ica_out
end

"""
    restore_ica_components(dat::ContinuousData, ica::InfoIca, components_to_restore::Vector{Int})

Restore ICA components to ContinuousData and return restored data and updated ICA result.

# Arguments
- `dat::ContinuousData`: ContinuousData object containing the data to restore
- `ica::InfoIca`: ICA result object with stored activations
- `components_to_restore::Vector{Int}`: Vector of component indices to restore

# Returns
- `Tuple{ContinuousData, InfoIca}`: Copy of input data with components restored, and copy of ICA result with updated removed_activations

# Example
```julia
restored_dat, ica_updated = restore_ica_components(dat, ica_result, [1, 3, 5])
```
"""
function restore_ica_components(dat::ContinuousData, ica::InfoIca; component_selection::Function = components())
    dat_out = copy(dat)
    ica_out = copy(ica)  # Use our custom copy method
    restore_ica_components!(dat_out.data, ica_out, component_selection = component_selection)
    return dat_out, ica_out
end

"""
    restore_ica_components!(dat::ContinuousData, ica::InfoIca, components_to_restore::Vector{Int})

Restore ICA components to ContinuousData in-place and update ICA result.

# Arguments
- `dat::ContinuousData`: ContinuousData object to restore in-place
- `ica::InfoIca`: ICA result object with stored activations (will be updated)
- `components_to_restore::Vector{Int}`: Vector of component indices to restore

# Returns
- `nothing` (data and ICA result are modified in-place)

# Example
```julia
restore_ica_components!(dat, ica_result, [1, 3, 5])
```
"""
function restore_ica_components!(dat::ContinuousData, ica::InfoIca; component_selection::Function = components())
    return restore_ica_components!(dat.data, ica, component_selection = component_selection)
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
    data_subset_df = dat.data[selected_samples, relevant_cols]
    dat_matrix = permutedims(Matrix(data_subset_df))
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
    z_threshold::Float64 = 3.0,
    min_correlation::Float64 = 0.5,
    two_step::Bool = true,
)

    # Check basic inputs - return nothing as first value if EOG channels are missing
    # TODO: what about when only one is available?
    if !(vEOG_channel in propertynames(dat.data)) || !(hEOG_channel in propertynames(dat.data))
        return nothing, DataFrame()
    end

    # Get samples to use
    selected_samples = get_selected_samples(dat, sample_selection)
    if isempty(selected_samples)
        @minimal_warning "No samples remaining after applying exclude criteria. Cannot identify eye components."
        return Dict(:vEOG => Int[], :hEOG => Int[]), DataFrame()
    end

    # Function to calculate initial EOG-component correlations 
    function calculate_correlations(eog_signal)
        corrs = zeros(n_components)
        for comp_idx = 1:n_components
            corrs[comp_idx] = abs(cor(components[comp_idx, :], eog_signal))
        end
        return corrs
    end

    # Function to calculate spatial correlations between identified vEOG/hEOG components and remaining components (used in Step2)
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

    # Function to calculate lagged correlations between idenitified vEOG/hEOG components and remaining components (used in Step2) 
    function calculate_lagged_correlations(remaining_components, primary_components, components_matrix, lp_filter, max_lag_samples, lag_step)

        lagged_corrs = fill(NaN, n_components)
        if isempty(remaining_components) || isempty(primary_components)
            return lagged_corrs
        end
        
        # Pre-filter all primary EOG components (they don't depend on remaining components)
        eog_components_filtered = Dict{Int, Vector{Float64}}()
        for eog_comp_idx in primary_components
            eog_components_filtered[eog_comp_idx] = abs.(filtfilt(lp_filter.filter_object, components_matrix[eog_comp_idx, :]))
        end
        
        for comp_idx in remaining_components
            # Apply low-pass filter to remove high-frequency noise
            comp_ts = abs.(filtfilt(lp_filter.filter_object, components_matrix[comp_idx, :]))
            
            max_corr = 0.0
            for eog_comp_idx in primary_components
                eog_comp_ts = eog_components_filtered[eog_comp_idx]
                # Use crosscor to compute correlations at all lags, then take maximum absolute value
                corrs = crosscor( eog_comp_ts, comp_ts, -max_lag_samples:lag_step:max_lag_samples)
                max_corr_val = maximum(abs.(corrs))
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

    # Calculate components for valid samples
    components, n_components = _prepare_ica_data_matrix(dat, ica, selected_samples)

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
    
    # Calculate vEOG/hEOG step2: lagged correlations with primary vEOG/hEOG components
    # Create low-pass filter for component time series (EOG artifacts are typically < 15 Hz)
    lp_filter = create_filter("lp", "iir", 10.0, dat.sample_rate; order = 3)
    
    # Lag range: +- 100ms  
    # Convert to samples based on sample rate
    max_lag_samples = round(Int, 100.0 * dat.sample_rate / 1000.0) # 100 ms lag
    lag_step = max(1, round(Int, 5.0 * dat.sample_rate / 1000.0)) # 5 ms step
    
    # Calculate lagged correlations for vEOG and hEOG
    vEOG_temporal_corr = calculate_lagged_correlations(remaining_vEOG, primary_vEOG, components, lp_filter, max_lag_samples, lag_step)
    hEOG_temporal_corr = calculate_lagged_correlations(remaining_hEOG, primary_hEOG, components, lp_filter, max_lag_samples, lag_step)
    
    # Helper function to compute z-scores only for remaining components
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
    function identify_secondary_components(remaining_components, corr_vector, corr_z_vector, min_corr_threshold, z_threshold)
        isempty(remaining_components) && return Int[]
        corr_mask = corr_vector[remaining_components] .> min_corr_threshold
        z_mask = abs.(corr_z_vector[remaining_components]) .> z_threshold
        combined_mask = corr_mask .& z_mask
        return remaining_components[combined_mask]
    end
    
    # Identify secondary components for vEOG: spatial OR temporal (each requires both corr and z-score)
    secondary_vEOG_spatial = identify_secondary_components(remaining_vEOG, vEOG_spatial_corr, vEOG_spatial_corr_z, min_correlation, z_threshold)
    secondary_vEOG_temporal = identify_secondary_components(remaining_vEOG, vEOG_temporal_corr, vEOG_temporal_corr_z, min_correlation, z_threshold)
    
    # Identify secondary components for hEOG: spatial OR temporal (each requires both corr and z-score)
    secondary_hEOG_spatial = identify_secondary_components(remaining_hEOG, hEOG_spatial_corr, hEOG_spatial_corr_z, min_correlation, z_threshold)
    secondary_hEOG_temporal = identify_secondary_components(remaining_hEOG, hEOG_temporal_corr, hEOG_temporal_corr_z, min_correlation, z_threshold)
    
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
)

    # Data Preparation 
    selected_samples = get_selected_samples(dat, sample_selection)
    if isempty(selected_samples)
        @minimal_warning "No samples remaining after applying exclude criteria."
        return Int[], DataFrame()
    end

    # Calculate components for valid samples
    components_subset, n_components = _prepare_ica_data_matrix(dat, ica, selected_samples)

    # Early return if no components
    if n_components == 0
        return Int[], DataFrame(
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

    for comp_idx = 1:n_components
        component_ts = components_subset[comp_idx, :]
        
        # Z-score the time series for consistent peak detection across components
        ts_zscored = StatsBase.zscore(component_ts)

        # Find prominent peaks 
        peak_indices = _find_peaks(ts_zscored; min_prominence_std = min_prominence_std)
        num_peaks = length(peak_indices)

        # Calculate IBI metrics
        num_valid_ibis, mean_ibi, std_ibi, peak_ratio = _calculate_ibi_metrics(
            peak_indices, dat.sample_rate, min_ibi_s, max_ibi_s
        )

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
function identify_spatial_kurtosis_components(ica::InfoIca; z_threshold::Float64 = 3.0)

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
    metrics_df = DataFrame(
        :Component => 1:n_components,
        :spatial_kurtosis => spatial_kurtosis,
        :z_spatial_kurtosis => z_spatial_kurtosis,
    )

    return high_kurtosis_comps, metrics_df
end


"""
    identify_line_noise_components(dat::ContinuousData, ica::InfoIca;
                                 exclude_samples::Union{Nothing,Vector{Symbol}} = [:is_extreme_value],
                                 line_freq::Real=50.0,
                                 freq_bandwidth::Real=1.0,
                                 z_threshold::Float64=3.0,
                                 min_harmonic_power::Real=0.5)

Identify ICA components with strong line noise characteristics.

# Arguments
- `dat::ContinuousData`: The continuous data.
- `ica::InfoIca`: The ICA result object.

# Keyword Arguments
- `exclude_samples::Union{Nothing,Vector{Symbol}}`: Optional vector of Bool columns in `dat.data` marking samples to exclude. Defaults to `[:is_extreme_value]`.
- `line_freq::Real`: Line frequency in Hz (default: 50.0 for European power).
- `freq_bandwidth::Real`: Bandwidth around line frequency to consider (default: 1.0 Hz).
- `z_threshold::Float64`: Z-score threshold for identifying line noise components (default: 3.0).
- `min_harmonic_power::Real`: Minimum power ratio of 2nd harmonic relative to fundamental (default: 0.5).

# Returns
- `Vector{Int}`: Indices of components with strong line noise characteristics.
- `DataFrame`: DataFrame containing spectral metrics for all components.
"""
function identify_line_noise_components(
    dat::ContinuousData,
    ica::InfoIca;
    sample_selection::Function = samples(),
    line_freq::Real = 50.0,
    freq_bandwidth::Real = 1.0,
    z_threshold::Float64 = 3.0,
    min_harmonic_power::Real = 0.5,
)

    # Get samples to use
    samples_to_use = get_selected_samples(dat, sample_selection)
    if isempty(samples_to_use)
        @minimal_warning "No samples remaining after applying exclude criteria. Cannot identify line noise components."
        return Int[], DataFrame()
    end

    # Calculate components for valid samples
    components, n_components = _prepare_ica_data_matrix(dat, ica, samples_to_use)

    # Calculate power spectrum for each component
    # Use a reasonable FFT size (power of 2, but not too large)
    n_samples = size(components, 2)
    nfft = nextpow(2, n_samples)  # Cap at 2^16 points
    freqs = FFTW.rfftfreq(nfft, dat.sample_rate)
    n_freqs = length(freqs)
    psd = zeros(n_freqs, n_components)

    for i = 1:n_components
        signal = components[i, :]
        
        # Prepare signal for FFT: truncate or zero-pad to nfft points
        if length(signal) > nfft
            signal_fft = signal[1:nfft]
        elseif length(signal) < nfft
            signal_fft = zeros(nfft)
            signal_fft[1:length(signal)] = signal
        else
            signal_fft = signal
        end
        
        # Calculate power spectral density
        psd[:, i] = abs2.(FFTW.rfft(signal_fft))
    end

    # Find frequency bands for line frequency and 2nd harmonic
    line_band = findall(abs.(freqs .- line_freq) .<= freq_bandwidth)
    harmonic2_freq = line_freq * 2
    harmonic2_band = findall(abs.(freqs .- harmonic2_freq) .<= freq_bandwidth)

    # Calculate metrics for each component
    metrics = []
    for i = 1:n_components
        # Get power at line frequency band
        line_power = mean(psd[line_band, i])

        # Calculate power in surrounding frequencies (excluding line frequency band)
        surrounding_bands = setdiff(1:n_freqs, line_band)
        surrounding_power = mean(psd[surrounding_bands, i])

        # Calculate power ratio (line power relative to surrounding frequencies)
        power_ratio = line_power / (surrounding_power + eps())

        # Calculate 2nd harmonic power ratio (relative to fundamental)
        harmonic_ratio = if !isempty(harmonic2_band)
            harmonic2_power = mean(psd[harmonic2_band, i])
            harmonic2_power / (line_power + eps())
        else
            NaN  # 2nd harmonic frequency above Nyquist
        end

        # Store metrics
        push!(
            metrics,
            (
                Component = i,
                line_power = line_power,
                surrounding_power = surrounding_power,
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
    identify_components(dat::ContinuousData, ica::InfoIca; sample_selection::Function = samples(), kwargs...)

Identify all types of artifact components in one unified call.

# Arguments
- `dat::ContinuousData`: The continuous data
- `ica::InfoIca`: The ICA result object
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
# Basic usage
artifacts, metrics = identify_components(dat, ica_result)

# With custom sample selection
artifacts, metrics = identify_components(dat, ica_result, 
    sample_selection = samples_not(:is_extreme_value_100))

# Access specific metrics
eog_metrics = metrics[:eog_metrics]
ecg_metrics = metrics[:ecg_metrics]
```
"""
function identify_components(dat::ContinuousData, ica::InfoIca; sample_selection::Function = samples(), kwargs...)
    # Identify EOG components (pass through use_robust_zscore if provided)
    eog_comps, eog_metrics_df = identify_eog_components(dat, ica; sample_selection = sample_selection, kwargs...)

    # Identify ECG components
    ecg_comps, ecg_metrics_df = identify_ecg_components(dat, ica; sample_selection = sample_selection, kwargs...)

    # Identify line noise components
    line_noise_comps, line_noise_metrics_df = identify_line_noise_components(dat, ica; kwargs...)

    # Identify channel noise components (spatial kurtosis)
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



"""
    Base.show(io::IO, artifacts::ArtifactComponents)

Custom printing for ArtifactComponents structure.
"""
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

