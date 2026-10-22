# === Downstream ICA Components (Subtraction & Restoration) ===

"""
    subtract_ica_components!(dat::DataFrame, ica::InfoIca; component_selection::Union{Function, Nothing} = nothing)
    subtract_ica_components!(dfs::Vector{DataFrame}, ica::InfoIca; component_selection::Union{Function, Nothing} = nothing)
    subtract_ica_components!(dat::ContinuousData, ica::InfoIca; component_selection::Union{Function, Nothing} = nothing)
    subtract_ica_components!(dat::EpochData, ica::InfoIca; component_selection::Union{Function, Nothing} = nothing)
    subtract_ica_components(dat, ica::InfoIca; component_selection::Union{Function, Nothing} = nothing)

Remove ICA components from data. Mutating (`!`) versions operate in-place; non-mutating versions
return a copy of the data with components removed along with an updated `InfoIca` object.

# Arguments
- `dat`: `DataFrame`, `Vector{DataFrame}`, `ContinuousData`, or `EpochData`.
- `ica::InfoIca`: ICA model containing mixing, unmixing, and scale information.
- `component_selection::Union{Function, Nothing}`: Component selector specifying which components to remove.
  An explicit selector (e.g. `components([1, 2])` or `components(artifact_indices)`) is required to prevent
  accidental removal of all components. To intentionally subtract all components, pass `components()`.

# Returns
- Mutating versions (`!`): Returns the modified input data in-place.
- Non-mutating versions: Returns `(cleaned_data, ica_updated)`.

# Memory Note
Stored activations in `ica.removed_activations` occupy `Float64` memory proportional to
`k * n_samples` (e.g., ~14 MB per component per 30 minutes at 1 kHz). This allows exact lossless
restoration via `add_ica_components!`. If saving `ica` to disk (e.g., via JLD2 or BSON), calling
`empty!(ica.removed_activations)` reduces file size if restoration is no longer needed.

# Examples
```julia
# Remove artifact components (in-place)
subtract_ica_components!(dat, ica_result, component_selection = components([1, 3, 5]))

# Non-mutating version
cleaned_dat, ica_updated = subtract_ica_components(dat, ica_result, component_selection = components([1, 3]))
```
"""
function subtract_ica_components!(
    dat::DataFrame,
    ica::InfoIca;
    component_selection::Union{Function, Nothing} = nothing,
)
    _subtract_ica_components!(dat, ica; component_selection = component_selection)
    return dat
end

function subtract_ica_components!(
    dfs::Vector{DataFrame},
    ica::InfoIca;
    component_selection::Union{Function, Nothing} = nothing,
)
    _subtract_ica_components!(dfs, ica; component_selection = component_selection)
    return dfs
end

function subtract_ica_components!(
    dat::ContinuousData,
    ica::InfoIca;
    component_selection::Union{Function, Nothing} = nothing,
)
    actually_removed = _subtract_ica_components!(dat.data, ica; component_selection = component_selection)
    dat.analysis_info.n_ica_components_removed += length(actually_removed)
    return dat
end

function subtract_ica_components!(
    dat::EpochData,
    ica::InfoIca;
    component_selection::Union{Function, Nothing} = nothing,
)
    actually_removed = _subtract_ica_components!(dat.data, ica; component_selection = component_selection)
    dat.analysis_info.n_ica_components_removed += length(actually_removed)
    return dat
end

function subtract_ica_components(
    dat::DataFrame,
    ica::InfoIca;
    component_selection::Union{Function, Nothing} = nothing,
)
    dat_out = copy(dat)
    ica_out = copy(ica)
    subtract_ica_components!(dat_out, ica_out; component_selection = component_selection)
    return dat_out, ica_out
end

function subtract_ica_components(
    dfs::Vector{DataFrame},
    ica::InfoIca;
    component_selection::Union{Function, Nothing} = nothing,
)
    dfs_out = [copy(df) for df in dfs]
    ica_out = copy(ica)
    subtract_ica_components!(dfs_out, ica_out; component_selection = component_selection)
    return dfs_out, ica_out
end

function subtract_ica_components(
    dat::ContinuousData,
    ica::InfoIca;
    component_selection::Union{Function, Nothing} = nothing,
)
    dat_out = copy(dat)
    ica_out = copy(ica)
    subtract_ica_components!(dat_out, ica_out; component_selection = component_selection)
    return dat_out, ica_out
end

function subtract_ica_components(
    dat::EpochData,
    ica::InfoIca;
    component_selection::Union{Function, Nothing} = nothing,
)
    dat_out = copy(dat)
    ica_out = copy(ica)
    subtract_ica_components!(dat_out, ica_out; component_selection = component_selection)
    return dat_out, ica_out
end


"""
    add_ica_components!(dat::DataFrame, ica::InfoIca; component_selection::Union{Function, Nothing} = nothing)
    add_ica_components!(dfs::Vector{DataFrame}, ica::InfoIca; component_selection::Union{Function, Nothing} = nothing)
    add_ica_components!(dat::ContinuousData, ica::InfoIca; component_selection::Union{Function, Nothing} = nothing)
    add_ica_components!(dat::EpochData, ica::InfoIca; component_selection::Union{Function, Nothing} = nothing)
    add_ica_components(dat, ica::InfoIca; component_selection::Union{Function, Nothing} = nothing)

Restore previously removed ICA components using stored activations.
Mutating (`!`) versions operate in-place; non-mutating versions return a copy and an updated `InfoIca`.
Components must have been removed first (activations stored in `ica.removed_activations`).

By default (`component_selection = nothing`), all currently removed components stored in
`ica.removed_activations` are restored. You may also specify a subset (e.g. `components([1])`).

# Arguments
- `dat`: `DataFrame`, `Vector{DataFrame}`, `ContinuousData`, or `EpochData`.
- `ica::InfoIca`: The ICA object containing stored component activations.
- `component_selection::Union{Function, Nothing}`: Selector for components to restore. Defaults to all removed.

# Returns
- Mutating versions (`!`): Returns the restored input data in-place.
- Non-mutating versions: Returns `(restored_data, ica_updated)`.

# Examples
```julia
# Restore all removed components
add_ica_components!(dat, ica_result)

# Restore specific component
add_ica_components!(dat, ica_result, component_selection = components([1]))
```
"""
function add_ica_components!(
    dat::DataFrame,
    ica::InfoIca;
    component_selection::Union{Function, Nothing} = nothing,
)
    _add_ica_components!(dat, ica; component_selection = component_selection)
    return dat
end

function add_ica_components!(
    dfs::Vector{DataFrame},
    ica::InfoIca;
    component_selection::Union{Function, Nothing} = nothing,
)
    _add_ica_components!(dfs, ica; component_selection = component_selection)
    return dfs
end

function add_ica_components!(
    dat::ContinuousData,
    ica::InfoIca;
    component_selection::Union{Function, Nothing} = nothing,
)
    actually_added = _add_ica_components!(dat.data, ica; component_selection = component_selection)
    dat.analysis_info.n_ica_components_removed -= length(actually_added)
    return dat
end

function add_ica_components!(
    dat::EpochData,
    ica::InfoIca;
    component_selection::Union{Function, Nothing} = nothing,
)
    actually_added = _add_ica_components!(dat.data, ica; component_selection = component_selection)
    dat.analysis_info.n_ica_components_removed -= length(actually_added)
    return dat
end

function add_ica_components(
    dat::DataFrame,
    ica::InfoIca;
    component_selection::Union{Function, Nothing} = nothing,
)
    dat_out = copy(dat)
    ica_out = copy(ica)
    add_ica_components!(dat_out, ica_out; component_selection = component_selection)
    return dat_out, ica_out
end

function add_ica_components(
    dfs::Vector{DataFrame},
    ica::InfoIca;
    component_selection::Union{Function, Nothing} = nothing,
)
    dfs_out = [copy(df) for df in dfs]
    ica_out = copy(ica)
    add_ica_components!(dfs_out, ica_out; component_selection = component_selection)
    return dfs_out, ica_out
end

function add_ica_components(
    dat::ContinuousData,
    ica::InfoIca;
    component_selection::Union{Function, Nothing} = nothing,
)
    dat_out = copy(dat)
    ica_out = copy(ica)
    add_ica_components!(dat_out, ica_out; component_selection = component_selection)
    return dat_out, ica_out
end

function add_ica_components(
    dat::EpochData,
    ica::InfoIca;
    component_selection::Union{Function, Nothing} = nothing,
)
    dat_out = copy(dat)
    ica_out = copy(ica)
    add_ica_components!(dat_out, ica_out; component_selection = component_selection)
    return dat_out, ica_out
end


# === Internal Core Algorithms ===

"""
    _subtract_core!(X::Matrix{Float64}, ica::InfoIca, comps::Vector{Int}) -> Matrix{Float64}

In-place subtraction of selected ICA components from `X` (samples × channels).
Computes activations `acts = (Uc * (X' - mean)) / scale` and performs:
`X -= scale * acts' * Mc'` using BLAS `mul!`. Returns `acts` (k × samples).
"""
function _subtract_core!(X::Matrix{Float64}, ica::InfoIca, comps::Vector{Int})
    Uc = ica.unmixing[comps, :] # k × n_channels (contiguous copy for BLAS)
    Mc = ica.mixing[:, comps]   # n_channels × k (contiguous copy for BLAS)

    # Compute component activations: acts has shape k × n_samples
    acts = Uc * transpose(X)
    acts .-= Uc * ica.mean
    acts ./= ica.scale

    # In-place subtract artifact projection: X -= scale * acts' * Mc'
    mul!(X, transpose(acts), transpose(Mc), -ica.scale, 1.0)

    return acts
end

"""
    _add_core!(X::Matrix{Float64}, ica::InfoIca, comps::Vector{Int})

In-place restoration of selected ICA components into `X` (samples × channels)
using stored activations from `ica.removed_activations`.
Performs: `X += scale * acts' * Mc'`.
"""
function _add_core!(X::Matrix{Float64}, ica::InfoIca, comps::Vector{Int})
    n_samples = size(X, 1)
    acts = reduce(vcat, [ica.removed_activations[c] for c in comps]) # k × n_samples
    if size(acts, 2) != n_samples
        @minimal_error "Stored activations have $(size(acts, 2)) samples, but data has $n_samples samples."
    end

    Mc = ica.mixing[:, comps] # n_channels × k
    mul!(X, transpose(acts), transpose(Mc), ica.scale, 1.0)
    return nothing
end

"""Validate and filter components to be subtracted from an ICA model."""
function _validate_and_filter_subtract_components(
    ica::InfoIca,
    component_selection::Union{Function, Nothing},
)
    if isnothing(component_selection)
        @minimal_error "subtract_ica_components! requires an explicit 'component_selection' (e.g., component_selection = components([1, 2]) or components(artifacts)). To subtract all components, pass component_selection = components()."
    end

    raw_comps = get_selected_components(ica, component_selection)
    comps = sort!(unique(raw_comps))
    n_components = size(ica.unmixing, 1)
    if !all(1 .<= comps .<= n_components)
        @minimal_error "Selected components must be between 1 and $n_components"
    end

    already_removed = filter(c -> haskey(ica.removed_activations, c), comps)
    if !isempty(already_removed)
        @minimal_warning "Components $(already_removed) have already been removed from this ICA model and are in removed_activations. Skipping to avoid subtracting zeros or corrupting stored activations. Use a fresh copy(ica) or clear ica.removed_activations if applying to a new dataset."
        comps = filter(c -> !haskey(ica.removed_activations, c), comps)
    end

    return comps
end

"""Validate and filter components to be restored to an ICA model."""
function _validate_and_filter_add_components(
    ica::InfoIca,
    component_selection::Union{Function, Nothing},
)
    n_components = size(ica.unmixing, 1)
    comps = if isnothing(component_selection)
        sort(collect(keys(ica.removed_activations)))
    else
        sort!(unique(get_selected_components(ica, component_selection)))
    end

    if !all(1 .<= comps .<= n_components)
        @minimal_error "Components must be between 1 and $n_components"
    end

    # If the user selected all components (1:n_components), interpret as restoring all removed components
    if length(comps) == n_components && comps == collect(1:n_components)
        comps = sort(collect(keys(ica.removed_activations)))
    else
        for comp in comps
            if !haskey(ica.removed_activations, comp)
                @minimal_error "Component $comp has no stored activations to restore."
            end
        end
    end

    if isempty(comps)
        @minimal_warning "No stored removed activations found to restore."
        return Int[]
    end

    return comps
end

"""Internal DataFrame implementation of component subtraction."""
function _subtract_ica_components!(
    dat::DataFrame,
    ica::InfoIca;
    component_selection::Union{Function, Nothing} = nothing,
)
    comps = _validate_and_filter_subtract_components(ica, component_selection)
    isempty(comps) && return Int[]

    labels = ica.layout.data.label
    X = Matrix{Float64}(dat[!, labels])
    acts = _subtract_core!(X, ica, comps)

    for (i, comp_idx) in enumerate(comps)
        ica.removed_activations[comp_idx] = acts[i:i, :]
    end

    for (ch_idx, ch) in enumerate(labels)
        copyto!(dat[!, ch], view(X, :, ch_idx))
    end

    return comps
end

"""Internal epoched Vector{DataFrame} implementation of component subtraction."""
function _subtract_ica_components!(
    dfs::Vector{DataFrame},
    ica::InfoIca;
    component_selection::Union{Function, Nothing} = nothing,
)
    isempty(dfs) && return Int[]
    n_samples_per_epoch = nrow(dfs[1])
    if !all(nrow(df) == n_samples_per_epoch for df in dfs)
        @minimal_error "All epochs in Vector{DataFrame} must have the same number of samples ($n_samples_per_epoch)."
    end

    comps = _validate_and_filter_subtract_components(ica, component_selection)
    isempty(comps) && return Int[]

    labels = ica.layout.data.label
    n_channels = length(labels)
    n_epochs = length(dfs)
    n_total_samples = n_samples_per_epoch * n_epochs

    X = Matrix{Float64}(undef, n_total_samples, n_channels)
    for (ch_idx, ch) in enumerate(labels)
        col_out = view(X, :, ch_idx)
        for (ep_idx, df) in enumerate(dfs)
            src = df[!, ch]
            offset = (ep_idx - 1) * n_samples_per_epoch
            for s = 1:n_samples_per_epoch
                col_out[offset + s] = src[s]
            end
        end
    end

    acts = _subtract_core!(X, ica, comps)

    for (i, comp_idx) in enumerate(comps)
        ica.removed_activations[comp_idx] = acts[i:i, :]
    end

    for (ch_idx, ch) in enumerate(labels)
        col_in = view(X, :, ch_idx)
        for (ep_idx, df) in enumerate(dfs)
            dst = df[!, ch]
            offset = (ep_idx - 1) * n_samples_per_epoch
            for s = 1:n_samples_per_epoch
                dst[s] = col_in[offset + s]
            end
        end
    end

    return comps
end

"""Internal DataFrame implementation of component restoration."""
function _add_ica_components!(
    dat::DataFrame,
    ica::InfoIca;
    component_selection::Union{Function, Nothing} = nothing,
)
    comps = _validate_and_filter_add_components(ica, component_selection)
    isempty(comps) && return Int[]

    labels = ica.layout.data.label
    X = Matrix{Float64}(dat[!, labels])
    _add_core!(X, ica, comps)

    for (ch_idx, ch) in enumerate(labels)
        copyto!(dat[!, ch], view(X, :, ch_idx))
    end

    for comp in comps
        delete!(ica.removed_activations, comp)
    end

    return comps
end

"""Internal epoched Vector{DataFrame} implementation of component restoration."""
function _add_ica_components!(
    dfs::Vector{DataFrame},
    ica::InfoIca;
    component_selection::Union{Function, Nothing} = nothing,
)
    isempty(dfs) && return Int[]
    n_samples_per_epoch = nrow(dfs[1])
    if !all(nrow(df) == n_samples_per_epoch for df in dfs)
        @minimal_error "All epochs in Vector{DataFrame} must have the same number of samples ($n_samples_per_epoch)."
    end

    comps = _validate_and_filter_add_components(ica, component_selection)
    isempty(comps) && return Int[]

    labels = ica.layout.data.label
    n_channels = length(labels)
    n_epochs = length(dfs)
    n_total_samples = n_samples_per_epoch * n_epochs

    X = Matrix{Float64}(undef, n_total_samples, n_channels)
    for (ch_idx, ch) in enumerate(labels)
        col_out = view(X, :, ch_idx)
        for (ep_idx, df) in enumerate(dfs)
            src = df[!, ch]
            offset = (ep_idx - 1) * n_samples_per_epoch
            for s = 1:n_samples_per_epoch
                col_out[offset + s] = src[s]
            end
        end
    end

    _add_core!(X, ica, comps)

    for (ch_idx, ch) in enumerate(labels)
        col_in = view(X, :, ch_idx)
        for (ep_idx, df) in enumerate(dfs)
            dst = df[!, ch]
            offset = (ep_idx - 1) * n_samples_per_epoch
            for s = 1:n_samples_per_epoch
                dst[s] = col_in[offset + s]
            end
        end
    end

    for comp in comps
        delete!(ica.removed_activations, comp)
    end

    return comps
end


# Function barrier forces Julia to compile a fast inner loop for the concrete type of `src`
function _copy_col_to_matrix!(
    dst::Matrix{Float64},
    col_idx::Int,
    src::AbstractVector,
    selected_samples::Vector{Int},
)
    isempty(selected_samples) && return dst
    min_idx, max_idx = extrema(selected_samples)
    if min_idx < 1 || max_idx > length(src)
        @minimal_error "Sample indices out of bounds: range [$min_idx, $max_idx] exceeds column length $(length(src))"
    end
    @inbounds for (j, sample_idx) in enumerate(selected_samples)
        dst[col_idx, j] = src[sample_idx]
    end
    return dst
end

"""
    _prepare_ica_data_matrix(dat::ContinuousData, ica::InfoIca, selected_samples::Vector{Int})

Prepare data matrix and calculate ICA components for selected samples.

Extracts relevant channels, permutes to channels × samples format, centers (subtracts sample mean), 
scales by ICA scale factor, and applies the unmixing matrix to compute components.

# Arguments
- `dat::ContinuousData`: The continuous data
- `ica::InfoIca`: The ICA result object
- `selected_samples::Vector{Int}`: Vector of sample indices to use

# Returns
- `components::Matrix{Float64}`: ICA components (n_components × n_samples)
- `n_components::Int`: Number of ICA components

# Note
`_prepare_ica_data_matrix` centers with the mean of the selected samples (`mean(dat_matrix, dims = 2)`),
which is appropriate for sample statistics (correlations, spectra, kurtosis). This differs from
`subtract_ica_components!`, which centers with `ica.mean` (and mathematically cancels when projected back).
"""
function _prepare_ica_data_matrix(dat::ContinuousData, ica::InfoIca, selected_samples::Vector{Int})
    relevant_cols = ica.layout.data.label
    n_samples = length(selected_samples)
    n_channels = length(relevant_cols)
    n_data_samples = nrow(dat.data)

    if !isempty(selected_samples)
        min_idx, max_idx = extrema(selected_samples)
        if min_idx < 1 || max_idx > n_data_samples
            @minimal_error "Selected sample indices out of bounds: [$min_idx, $max_idx] vs data length $n_data_samples"
        end
    end

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
