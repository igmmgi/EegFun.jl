"""
    run_ica(dat::ContinuousData;
            n_components::Union{Nothing,Int} = nothing,
            sample_selection::Function = samples(),
            channel_selection::Function = channels(),
            include_extra::Bool = false,
            percentage_of_data::Real = 100.0,
            params::IcaPrms = IcaPrms())

Runs Independent Component Analysis (ICA) on EEG data. Preprocessing (e.g., filtering) should be applied prior to calling this function.

# Arguments
- `dat::ContinuousData`: The EEG data object.
- `n_components::Union{Nothing,Int}`: Number of ICA components (default: number of channels - 1).
- `sample_selection::Function`: Sample selector for quality filtering (default: include all samples).
- `channel_selection::Function`: Channel selector (default: layout channels).
- `include_extra::Bool`: Whether to allow channels outside the layout (e.g., EOG) in selection.
- `percentage_of_data::Real`: Percentage of good data to use for ICA (default: 100.0). Values < 100 enable faster computation by random subsampling.
- `params::IcaPrms`: ICA parameters.

# Returns
`InfoIca` with unmixing, mixing, sphere, variance, and metadata.

# Examples
```julia
# Basic ICA on layout channels
ica_result = run_ica(dat)

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

    ica_result = infomax_ica(dat_for_ica, ica_layout, n_components = n_components, params = params)

    return ica_result
end

"""
    run_ica(epoched_data::Vector{EpochData};
            n_components::Union{Nothing,Int} = nothing,
            sample_selection::Function = samples(),
            channel_selection::Function = channels(),
            include_extra::Bool = false,
            remove_duplicates::Bool = true,
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
- `params::IcaPrms`: ICA parameters

# Returns
`InfoIca` with unmixing, mixing, sphere, variance, and metadata.

# Examples
```julia
# Basic ICA on concatenated epochs
epochs = extract_epochs(dat, epoch_conditions, -1, 2)
ica_result = run_ica(epochs)

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
    ica_result = infomax_ica(concatenated_matrix, ica_layout, n_components = n_components, params = params)

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
    identify_eog_components(dat::ContinuousData, ica::InfoIca;
                          vEOG_channel::Symbol=:vEOG,
                          hEOG_channel::Symbol=:hEOG,
                          z_threshold::Float64=3.0,
                          sample_selection::Function = samples())

Identify ICA components potentially related to eye movements based on z-scored correlation.

# Arguments
- `dat::ContinuousData`: The continuous data containing EOG channels.
- `ica::InfoIca`: The ICA result object.

# Keyword Arguments
- `vEOG_channel::Symbol`: Name of the vertical EOG channel (default: :vEOG).
- `hEOG_channel::Symbol`: Name of the horizontal EOG channel (default: :hEOG).
- `z_threshold::Float64`: Absolute Z-score threshold for identification (default: 3.0).
- `sample_selection::Function`: Function to select samples from `dat.data`. Defaults to `samples()`.

# Returns
- `Dict{Symbol, Vector{Int}}`: Dictionary containing:
  - `:vEOG`: Vector of indices identified for vertical eye movements.
  - `:hEOG`: Vector of indices identified for horizontal eye movements.
- `DataFrame`: DataFrame containing detailed correlation metrics per component.
"""
function identify_eog_components(
    dat::ContinuousData,
    ica::InfoIca;
    vEOG_channel::Symbol = :vEOG,
    hEOG_channel::Symbol = :hEOG,
    sample_selection::Function = samples(),
    z_threshold::Float64 = 3.0,
)

    # Check basic inputs
    if !(vEOG_channel in propertynames(dat.data))
        @minimal_error "Vertical EOG channel $vEOG_channel not found in data"
    end
    if !(hEOG_channel in propertynames(dat.data))
        @minimal_error "Horizontal EOG channel $hEOG_channel not found in data"
    end

    # Get samples to use
    selected_samples = get_selected_samples(dat, sample_selection)
    if isempty(selected_samples)
        @minimal_warning "No samples remaining after applying exclude criteria. Cannot identify eye components."
        return Dict(:vEOG => Int[], :hEOG => Int[]), DataFrame()
    end

    # Get EOG signals for valid samples only
    vEOG = dat.data[selected_samples, vEOG_channel]
    hEOG = dat.data[selected_samples, hEOG_channel]

    # Prepare data matrix for valid samples
    relevant_cols = vcat(ica.layout.data.label)
    data_subset_df = dat.data[selected_samples, relevant_cols]
    dat_matrix = permutedims(Matrix(data_subset_df))
    dat_matrix .-= mean(dat_matrix, dims = 2)
    dat_matrix ./= ica.scale

    # Calculate components for valid samples
    components = ica.unmixing * dat_matrix
    n_components = size(components, 1)

    # Function to calculate correlations for all components
    function calculate_correlations(eog_signal)
        corrs = zeros(n_components)
        for comp_idx = 1:n_components
            corrs[comp_idx] = abs(cor(components[comp_idx, :], eog_signal))
        end
        return corrs
    end

    identified_vEOG = Int[]
    identified_hEOG = Int[]
    vEOG_corr_z = Float64[]
    hEOG_corr_z = Float64[]
    vEOG_corrs = Float64[]
    hEOG_corrs = Float64[]

    # vEOG
    vEOG_corrs = calculate_correlations(vEOG)
    vEOG_corr_z = StatsBase.zscore(vEOG_corrs)
    identified_vEOG = findall(abs.(vEOG_corr_z) .> z_threshold)

    # hEOG
    hEOG_corrs = calculate_correlations(hEOG)
    hEOG_corr_z = StatsBase.zscore(hEOG_corrs)
    identified_hEOG = findall(abs.(hEOG_corr_z) .> z_threshold)

    sort!(identified_vEOG)
    sort!(identified_hEOG)

    result_dict = Dict{Symbol,Vector{Int}}(:vEOG => identified_vEOG, :hEOG => identified_hEOG)

    metrics_df = DataFrame(
        :Component => 1:n_components,
        :vEOG_corr => vEOG_corrs,
        :vEOG_zscore => vEOG_corr_z,
        :hEOG_corr => hEOG_corrs,
        :hEOG_zscore => hEOG_corr_z,
    )

    return result_dict, metrics_df
end


"""
    identify_ecg_components(dat::ContinuousData, ica::InfoIca;
                              min_bpm::Real=40, max_bpm::Real=120,
                              min_prominence_std::Real=2.5,
                              min_peaks::Int=10,
                              max_ibi_std_s::Real=0.05,
                              sample_selection::Function = samples(),
                              plot_component_index::Int = 0

Identify ICA components potentially related to EKG artifacts based on peak detection
and interval regularity, using only samples consistent with ICA calculation.

# Arguments
- `dat::ContinuousData`: The continuous data (needed for sampling rate `fs` and sample selection columns).
- `ica::InfoIca`: The ICA result object.

# Keyword Arguments
- `min_bpm::Real`: Minimum plausible heart rate in beats per minute (default: 40).
- `max_bpm::Real`: Maximum plausible heart rate in beats per minute (default: 120).
- `min_prominence_std::Real`: Minimum peak prominence in standard deviations above mean (default: 2.5).
- `min_peaks::Int`: Minimum number of prominent peaks required within plausible heart rate range (default: 10).
- `max_ibi_std_s::Real`: Maximum standard deviation of the inter-beat intervals (in seconds) for component to be flagged (default: 0.05).
- `include_samples::Union{Nothing,Vector{Symbol}}`: Optional vector of Bool columns in `dat.data` marking samples to *include*. Defaults to `nothing` (include all unless excluded).
- `exclude_samples::Union{Nothing,Vector{Symbol}}`: Optional vector of Bool columns in `dat.data` marking samples to *exclude*. Defaults to `[:is_extreme_value]`.

# Returns
- `Vector{Int}`: Sorted vector of indices identified as potential EKG components.
- `DataFrame`: DataFrame containing metrics for each component (calculated on the included samples):
  - `:Component`: Component index (1 to n).
  - `:num_peaks`: Number of detected prominent peaks.
  - `:num_valid_ibis`: Number of inter-beat intervals within the plausible BPM range.
  - `:mean_ibi_s`: Mean inter-beat interval in seconds (if num_valid_ibis > 0).
  - `:std_ibi_s`: Standard deviation of inter-beat intervals in seconds (if num_valid_ibis > 1).
  - `:is_ekg_artifact`: Boolean flag indicating if component met the criteria.
"""
function identify_ecg_components(
    dat::ContinuousData,
    ica::InfoIca;
    min_bpm::Real = 40,
    max_bpm::Real = 130,
    min_prominence_std::Real = 5,
    min_peaks::Int = 10,
    max_ibi_std_s::Real = 0.2,
    min_peak_ratio::Real = 0.7,
    sample_selection::Function = samples(),
    plot_component_index::Int = 0,
)

    # Data Preparation 
    selected_samples = get_selected_samples(dat, sample_selection)
    if isempty(selected_samples)
        @minimal_warning "No samples remaining after applying exclude criteria."
        return Int[], DataFrame()
    end

    # Process data
    relevant_cols = ica.layout.data.label
    data_subset_df = dat.data[selected_samples, relevant_cols]
    dat_matrix_subset = permutedims(Matrix(data_subset_df))
    dat_matrix_subset .-= mean(dat_matrix_subset, dims = 2)
    dat_matrix_subset ./= ica.scale

    components_subset = ica.unmixing * dat_matrix_subset
    n_components = size(ica.unmixing, 1)

    # Convert BPM to plausible IBI range
    min_ibi_s = 60.0 / max_bpm
    max_ibi_s = 60.0 / min_bpm

    # Store results
    metrics = []
    identified_ecg = Int[]

    # Loop through components
    for comp_idx = 1:n_components
        ts = components_subset[comp_idx, :]

        # Find prominent peaks 
        peak_indices = _find_peaks(ts; min_prominence_std = min_prominence_std)

        # for debugging
        if comp_idx == plot_component_index
            fig, ax = lines(ts, color = :black)
            scatter!(ax, peak_indices, ts[peak_indices], color = :red)
            display(fig)
        end

        num_peaks = length(peak_indices)
        mean_ibi = NaN
        std_ibi = NaN
        num_valid_ibis = 0
        peak_ratio = 0.0
        is_ecg = false

        if num_peaks >= 2
            # Calculate IBIs
            ibis_s = diff(peak_indices) ./ dat.sample_rate

            # Filter valid IBIs
            valid_ibi_mask = (ibis_s .>= min_ibi_s) .& (ibis_s .<= max_ibi_s)
            valid_ibis = ibis_s[valid_ibi_mask]
            num_valid_ibis = length(valid_ibis)

            # Calculate peak ratio (new metric)
            peak_ratio = num_valid_ibis / (num_peaks - 1)

            if num_valid_ibis > 1
                mean_ibi = mean(valid_ibis)
                std_ibi = std(valid_ibis)

                # Apply stricter criteria
                if num_valid_ibis >= (min_peaks - 1) && std_ibi <= max_ibi_std_s && peak_ratio >= min_peak_ratio  # Add peak ratio criterion
                    is_ecg = true
                    push!(identified_ecg, comp_idx)
                end
            elseif num_valid_ibis == 1
                mean_ibi = valid_ibis[1]
                std_ibi = 0.0
            end
        end

        # Calculate heart rate if we have valid IBI
        heart_rate_bpm = isnan(mean_ibi) || mean_ibi <= 0 ? NaN : 60.0 / mean_ibi

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
                is_ecg_artifact = is_ecg,
            ),
        )
    end

    # Finalize results
    metrics_df = DataFrame(metrics)
    if isempty(metrics)
        metrics_df = DataFrame(
            :Component => 1:n_components,
            :num_peaks => 0,
            :num_valid_ibis => 0,
            :mean_ibi_s => NaN,
            :std_ibi_s => NaN,
            :peak_ratio => 0.0,
            :heart_rate_bpm => NaN,
            :is_ecg_artifact => false,
        )
    end
    sort!(identified_ecg)

    return identified_ecg, metrics_df

end

"""
    identify_spatial_kurtosis_components(dat::ContinuousData, ica::InfoIca;
                                      exclude_samples::Union{Nothing,Vector{Symbol}} = [:is_extreme_value],
                                      z_threshold::Float64 = 3.0)

Identify ICA components with high spatial kurtosis (localized, spot-like activity).

# Arguments
- `dat::ContinuousData`: The continuous data.
- `ica::InfoIca`: The ICA result object.

# Keyword Arguments
- `exclude_samples::Union{Nothing,Vector{Symbol}}`: Optional vector of Bool columns in `dat.data` marking samples to exclude. Defaults to `[:is_extreme_value]`.
- `z_threshold::Float64`: Z-score threshold for identifying high spatial kurtosis components (default: 3.0).

# Returns
- `Vector{Int}`: Indices of components with high spatial kurtosis.
- `DataFrame`: DataFrame containing spatial kurtosis values and z-scores for all components.
"""
function identify_spatial_kurtosis_components(dat::ContinuousData, ica::InfoIca; z_threshold::Float64 = 3.0)
    # Calculate spatial kurtosis for each component's weights
    n_components = size(ica.mixing, 2)
    spatial_kurtosis = Float64[]

    for i = 1:n_components
        # Get component weights
        weights = ica.mixing[:, i]
        # Calculate kurtosis of the weights
        k = kurtosis(weights)
        push!(spatial_kurtosis, k)
    end

    # Calculate z-scores of spatial kurtosis values
    spatial_kurtosis_z = StatsBase.zscore(spatial_kurtosis)

    # Identify components with high spatial kurtosis (using z-scores)
    high_kurtosis_comps = findall(spatial_kurtosis_z .> z_threshold)  # Only positive deviations (localized activity)
    sort!(high_kurtosis_comps)

    # Create metrics DataFrame
    metrics_df = DataFrame(
        :Component => 1:n_components,
        :SpatialKurtosis => spatial_kurtosis,
        :SpatialKurtosisZScore => spatial_kurtosis_z,
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
- `min_harmonic_power::Real`: Minimum power ratio of harmonics relative to fundamental (default: 0.5).

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

    # Prepare data matrix for valid samples
    relevant_cols = vcat(ica.layout.data.label)
    data_subset_df = dat.data[samples_to_use, relevant_cols]
    dat_matrix = permutedims(Matrix(data_subset_df))
    dat_matrix .-= mean(dat_matrix, dims = 2)
    dat_matrix ./= ica.scale

    # Calculate components for valid samples
    components = ica.unmixing * dat_matrix
    n_components = size(components, 1)
    fs = dat.sample_rate

    # Calculate power spectrum for each component
    # Use a reasonable FFT size (power of 2, but not too large)
    nfft = min(nextpow(2, size(components, 2)), 2^16)  # Cap at 2^16 points
    freqs = FFTW.rfftfreq(nfft, fs)
    psd = zeros(length(freqs), n_components)

    for i = 1:n_components
        # Zero-pad or truncate to nfft points
        signal = components[i, :]
        if length(signal) > nfft
            signal = signal[1:nfft]
        elseif length(signal) < nfft
            signal = [signal; zeros(nfft - length(signal))]
        end
        psd[:, i] = abs2.(FFTW.rfft(signal))
    end

    # Find indices for line frequency and harmonics
    line_idx = findmin(abs.(freqs .- line_freq))[2]
    line_band = findall(abs.(freqs .- line_freq) .<= freq_bandwidth)

    # Calculate metrics for each component
    metrics = []
    for i = 1:n_components
        # Get power at line frequency and surrounding band
        line_power = mean(psd[line_band, i])

        # Calculate power in surrounding bands (excluding line frequency)
        surrounding_bands = setdiff(1:length(freqs), line_band)
        surrounding_power = mean(psd[surrounding_bands, i])

        # Calculate power ratio
        power_ratio = line_power / (surrounding_power + eps())

        # Check for harmonics (2x and 3x line frequency)
        harmonic_powers = Float64[]
        for h = 2:3
            harmonic_freq = line_freq * h
            harmonic_idx = findmin(abs.(freqs .- harmonic_freq))[2]
            harmonic_band = findall(abs.(freqs .- harmonic_freq) .<= freq_bandwidth)
            harmonic_power = mean(psd[harmonic_band, i])
            push!(harmonic_powers, harmonic_power / line_power)
        end

        # Store metrics
        push!(
            metrics,
            (
                Component = i,
                LinePower = line_power,
                SurroundingPower = surrounding_power,
                PowerRatio = power_ratio,
                Harmonic2Ratio = harmonic_powers[1],
                Harmonic3Ratio = harmonic_powers[2],
            ),
        )
    end

    # Create metrics DataFrame
    metrics_df = DataFrame(metrics)

    # Calculate z-scores of power ratios
    power_ratio_z = StatsBase.zscore(metrics_df.PowerRatio)
    metrics_df[!, :PowerRatioZScore] = power_ratio_z

    # Identify components with strong line noise characteristics
    line_noise_comps = findall(power_ratio_z .> z_threshold)

    # Additional check for harmonics
    if !isempty(line_noise_comps)
        harmonic_mask =
            (metrics_df.Harmonic2Ratio .> min_harmonic_power) .| (metrics_df.Harmonic3Ratio .> min_harmonic_power)
        line_noise_comps = intersect(line_noise_comps, findall(harmonic_mask))
    end

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
    # Identify EOG components
    eog_comps, eog_metrics_df = identify_eog_components(dat, ica; sample_selection = sample_selection, kwargs...)

    # Identify ECG components  
    ecg_comps, ecg_metrics_df = identify_ecg_components(dat, ica; sample_selection = sample_selection, kwargs...)

    # Identify line noise components
    line_noise_comps, line_noise_metrics_df = identify_line_noise_components(dat, ica; kwargs...)

    # Identify channel noise components (spatial kurtosis)
    channel_noise_comps, channel_noise_metrics_df = identify_spatial_kurtosis_components(dat, ica; kwargs...)

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
