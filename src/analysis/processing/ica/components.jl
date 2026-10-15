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


# === Downstream Components ===
function subtract_ica_components!(dat::DataFrame, ica::InfoIca; component_selection::Function = components())
    components_to_remove = get_selected_components(ica, component_selection)
    n_components = size(ica.unmixing, 1)
    if !all(1 .<= components_to_remove .<= n_components)
        @minimal_error("Components must be between 1 and $n_components")
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
        @minimal_error("Components must be between 1 and $n_components")
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

"""
    subtract_ica_components(dat, ica::InfoIca; component_selection::Function = components())

Non-mutating version of `subtract_ica_components!`.
"""
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
        @minimal_error("Components must be between 1 and $n_components")
    end

    # Check that all components to restore have stored activations
    for comp in components_to_restore
        if !haskey(ica.removed_activations, comp)
            @minimal_error("Component $comp has no stored activations to restore")
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
    add_ica_components(dat, ica::InfoIca; component_selection::Function = components())

Non-mutating version of `add_ica_components!`.
"""
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

