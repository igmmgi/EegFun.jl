# === UNIFIED CHANNEL REPAIR INTERFACE ===
"""
    repair_channels!(data, channels_to_repair; method=:neighbor_interpolation, kwargs...)

Repair bad EEG channels using the specified method.

# Arguments
- `data`: EEG data object (ContinuousData, EpochData, or AbstractMatrix)
- `channels_to_repair::Vector{Symbol}`: List of channel labels to repair

# Keyword Arguments
- `method::Symbol`: Repair method to use
  - `:neighbor_interpolation` (default): Weighted neighbor interpolation
  - `:spherical_spline`: Spherical spline interpolation
- `neighbours_dict::Union{OrderedDict, Nothing}`: Neighbor information (for neighbor_interpolation)
- `m::Int`: Order of Legendre polynomials (for spherical_spline, default: 4)
- `lambda::Float64`: Regularization parameter (for spherical_spline, default: 1e-5)
- `epoch_selection::Function`: Epoch selection predicate (for EpochData, default: all epochs)

# Examples
```julia
# Neighbor interpolation (default)
repair_channels!(dat, [:Fp1, :Fp2])

# Spherical spline
repair_channels!(dat, [:Fp1, :Fp2], method=:spherical_spline)

# With custom parameters
repair_channels!(dat, [:Fp1], method=:spherical_spline, m=6, lambda=1e-6)
```
"""
function repair_channels!(data, channels_to_repair; method::Symbol = :neighbor_interpolation, kwargs...)
    if method == :neighbor_interpolation
        _repair_channels_neighbor!(data, channels_to_repair; kwargs...)
    elseif method == :spherical_spline
        _repair_channels_spherical!(data, channels_to_repair; kwargs...)
    else
        throw(ArgumentError("Unknown repair method: $method. Use :neighbor_interpolation or :spherical_spline"))
    end
end

@add_nonmutating repair_channels!

# Helper function for neighbor interpolation
function _repair_channels_neighbor!(data, channels_to_repair; neighbours_dict = nothing, kwargs...)
    if isnothing(neighbours_dict)
        # For data objects, get neighbors from layout
        if hasfield(typeof(data), :layout)
            layout = data.layout
            if isnothing(layout.neighbours)
                get_layout_neighbours_xyz!(layout, 0.5)
            end
            neighbours_dict = layout.neighbours
        else
            throw(ArgumentError("neighbours_dict must be provided for AbstractMatrix data"))
        end
    end

    # Use multiple dispatch
    _repair_channels_neighbor!(data, channels_to_repair, neighbours_dict; kwargs...)
end

# === CORE IMPLEMENTATION FUNCTIONS ===

"""
    _repair_channels_neighbor!(data::AbstractMatrix, channels::Vector{Symbol}, channels_to_repair::Vector{Symbol}, neighbours_dict::OrderedDict{Symbol, Neighbours})

Repair bad channels in a data matrix using weighted neighbor interpolation.
"""
function _repair_channels_neighbor!(
    data::AbstractMatrix,
    channels::Vector{Symbol},
    channels_to_repair::Vector{Symbol},
    neighbours_dict::OrderedDict{Symbol,Neighbours},
)
    # Check that channels match data dimensions
    n_channels = size(data, 1)
    if length(channels) != n_channels
        throw(ArgumentError("Number of channels ($(length(channels))) doesn't match data dimensions ($n_channels)"))
    end

    # Create channel index lookup
    ch_indices = Dict(ch => i for (i, ch) in enumerate(channels))

    # Process each bad channel
    for bad_ch in channels_to_repair
        # Get bad channel index
        bad_idx = get(ch_indices, bad_ch, nothing)
        if isnothing(bad_idx)
            @minimal_warning "Channel $bad_ch not found in channel list, skipping"
            continue
        end

        # Get neighbors information
        neighbours = get(neighbours_dict, bad_ch, nothing)
        if isnothing(neighbours) || isempty(neighbours.channels)
            @minimal_warning "No neighbors found for channel $bad_ch, skipping"
            continue
        end
        
        # Require at least 2 neighbors for meaningful interpolation
        if length(neighbours.channels) < 2
            @minimal_warning "Channel $bad_ch has only $(length(neighbours.channels)) neighbor(s), need at least 2 for repair, skipping"
            continue
        end

        # Since check_channel_neighbors ensures ALL neighbors are good, 
        # we can use all neighbors directly (weights already normalized to sum to 1.0)
        # For each time point, calculate weighted average of neighbors
        n_timepoints = size(data, 2)
        for t = 1:n_timepoints
            weighted_sum = 0.0

            # Sum weighted neighbor values
            for (i, neighbour_ch) in enumerate(neighbours.channels)
                neigh_idx = get(ch_indices, neighbour_ch, nothing)
                if !isnothing(neigh_idx)
                    weight = neighbours.weights[i]
                    weighted_sum += data[neigh_idx, t] * weight
                end
            end

            # Update bad channel data
            data[bad_idx, t] = weighted_sum
        end

        @info "Repaired channel $bad_ch using weighted neighbor interpolation with neighbors: $(neighbours.channels)"
    end

    return nothing
end

"""
    _repair_channels_neighbor!(data::ContinuousData, channels_to_repair::Vector{Symbol}, neighbours_dict::OrderedDict{Symbol, Neighbours})

Repair bad channels in ContinuousData using weighted neighbor interpolation.
"""
function _repair_channels_neighbor!(
    data::ContinuousData,
    channels_to_repair::Vector{Symbol},
    neighbours_dict::OrderedDict{Symbol,Neighbours},
)
    # Extract channels vector from ContinuousData
    channels = data.layout.data.label

    # Extract the data matrix (transpose to get channels × time points)
    data_matrix = Matrix(data.data[:, data.layout.data.label])'

    # Call the core implementation
    _repair_channels_neighbor!(data_matrix, channels, channels_to_repair, neighbours_dict)

    # Update the data in the ContinuousData object
    data.data[:, data.layout.data.label] = data_matrix'

    return nothing
end

"""
    _repair_channels_neighbor!(data::EpochData, channels_to_repair::Vector{Symbol}, neighbours_dict::OrderedDict{Symbol, Neighbours}; epoch_selection::Function=epochs())

Repair bad channels in epoched EEG data using weighted neighbor interpolation.
"""
function _repair_channels_neighbor!(
    data::EpochData,
    channels_to_repair::Vector{Symbol},
    neighbours_dict::OrderedDict{Symbol,Neighbours};
    epoch_selection::Function = epochs(),
)
    # Get selected epochs
    selected_epochs = get_selected_epochs(data, epoch_selection)

    # Get channels vector
    channels = Symbol.(data.layout.data.label)

    # Repair each selected epoch
    for epoch_idx in selected_epochs
        epoch = data.data[epoch_idx]
        epoch_matrix = Matrix(epoch[:, channels])'

        _repair_channels_neighbor!(epoch_matrix, channels, channels_to_repair, neighbours_dict)

        # Update the epoch data
        epoch[:, channels] = epoch_matrix'
    end

    return nothing
end

"""
    _repair_channels_spherical!(data::AbstractMatrix, channels_to_repair::Vector{Symbol}, channels::Vector{Symbol}, layout::DataFrame; m::Int=4, lambda::Float64=1e-5)

Core implementation for repairing bad channels using spherical spline interpolation.
"""
function _repair_channels_spherical!(
    data::AbstractMatrix,
    channels_to_repair::Vector{Symbol},
    channels::Vector{Symbol},
    layout::DataFrame;
    m::Int = 4,
    lambda::Float64 = 1e-5,
)
    # Check that channels match data dimensions
    n_channels = size(data, 1)
    if length(channels) != n_channels
        throw(ArgumentError("Number of channels doesn't match data dimensions"))
    end

    # Extract and normalize coordinates
    coords = zeros(Float64, length(channels), 3)
    for (i, ch) in enumerate(channels)
        idx = findfirst(x -> x == ch, layout.label)
        if isnothing(idx)
            @warn "Channel $ch not found in layout, using (0,0,0)"
            continue
        end
        coords[i, :] = [layout.x3[idx], layout.y3[idx], layout.z3[idx]]
    end

    # Normalize coordinates to unit sphere
    norms = sqrt.(sum(coords .^ 2, dims = 2))
    coords = coords ./ norms

    # Process each bad channel
    for bad_ch in channels_to_repair
        # Get bad channel index
        bad_idx = findfirst(x -> x == bad_ch, channels)
        if isnothing(bad_idx)
            @minimal_warning "Channel $bad_ch not found in channel list, skipping"
            continue
        end

        # Get coordinates of bad channel
        bad_coords = coords[bad_idx, :]

        # Find all other channels (not the bad one)
        other_indices = setdiff(1:length(channels), bad_idx)
        other_coords = coords[other_indices, :]
        other_channels = channels[other_indices]

        # Calculate spherical spline interpolation
        n_timepoints = size(data, 2)
        for t = 1:n_timepoints
            # Get data from other channels at this time point
            other_data = data[other_indices, t]

            # Calculate spherical spline weights
            weights = _spherical_spline_weights(bad_coords, other_coords, m, lambda)

            # Interpolate the value
            interpolated_value = sum(weights .* other_data)

            # Update bad channel data
            data[bad_idx, t] = interpolated_value
        end

        @info "Repaired channel $bad_ch using spherical spline interpolation"
    end

    return nothing
end

"""
    _repair_channels_spherical!(data::ContinuousData, channels_to_repair::Vector{Symbol}; m::Int=4, lambda::Float64=1e-5)

Repair bad channels in ContinuousData using spherical spline interpolation.
"""
function _repair_channels_spherical!(
    data::ContinuousData,
    channels_to_repair::Vector{Symbol};
    m::Int = 4,
    lambda::Float64 = 1e-5,
)
    # Ensure 3D coordinates are available
    _ensure_coordinates_3d!(data.layout)

    channels = Symbol.(data.layout.data.label)
    data_matrix = Matrix(data.data[:, channels])'

    _repair_channels_spherical!(data_matrix, channels_to_repair, channels, data.layout.data; m = m, lambda = lambda)

    # Update the data
    data.data[:, channels] = data_matrix'

    return nothing
end

"""
    _repair_channels_spherical!(data::EpochData, channels_to_repair::Vector{Symbol}; epoch_selection::Function=epochs(), m::Int=4, lambda::Float64=1e-5)

Repair bad channels in epoched EEG data using spherical spline interpolation.
"""
function _repair_channels_spherical!(
    data::EpochData,
    channels_to_repair::Vector{Symbol};
    epoch_selection::Function = epochs(),
    m::Int = 4,
    lambda::Float64 = 1e-5,
)
    # Ensure 3D coordinates are available
    _ensure_coordinates_3d!(data.layout)

    # Get selected epochs
    selected_epochs = get_selected_epochs(data, epoch_selection)

    # Get channels vector
    channels = Symbol.(data.layout.data.label)

    # Repair each selected epoch
    for epoch_idx in selected_epochs
        epoch = data.data[epoch_idx]
        epoch_matrix = Matrix(epoch[:, channels])'

        _repair_channels_spherical!(
            epoch_matrix,
            channels_to_repair,
            channels,
            data.layout.data;
            m = m,
            lambda = lambda,
        )

        # Update the epoch data
        epoch[:, channels] = epoch_matrix'
    end

    return nothing
end

# === SPHERICAL SPLINE UTILITIES ===

"""
    _spherical_spline_weights(bad_coords, other_coords, m, lambda)

Calculate spherical spline interpolation weights.
"""
function _spherical_spline_weights(bad_coords, other_coords, m, lambda)
    n_others = size(other_coords, 1)

    # Calculate distances between bad channel and other channels
    distances = zeros(n_others)
    for i = 1:n_others
        distances[i] = acos(clamp(dot(bad_coords, other_coords[i, :]), -1.0, 1.0))
    end

    # Calculate spherical spline basis functions
    basis = zeros(n_others)
    for i = 1:n_others
        if distances[i] < 1e-10  # Very close channels
            basis[i] = 1.0
        else
            # Legendre polynomial basis
            cos_theta = cos(distances[i])
            basis[i] = _legendre_polynomial(cos_theta, m)
        end
    end

    # Add regularization
    if lambda > 0
        # Simple regularization: add small amount to diagonal
        weights = basis ./ (sum(basis) + lambda)
    else
        weights = basis ./ sum(basis)
    end

    return weights
end

"""
    _legendre_polynomial(cos_theta, m)

Calculate Legendre polynomial of order m at cos_theta.
"""
function _legendre_polynomial(cos_theta, m)
    if m == 0
        return 1.0
    elseif m == 1
        return cos_theta
    elseif m == 2
        return 0.5 * (3 * cos_theta^2 - 1)
    elseif m == 3
        return 0.5 * (5 * cos_theta^3 - 3 * cos_theta)
    elseif m == 4
        return 0.125 * (35 * cos_theta^4 - 30 * cos_theta^2 + 3)
    else
        # For higher orders, use recurrence relation
        p_prev2 = 1.0
        p_prev1 = cos_theta

        for n = 2:m
            p_current = ((2*n - 1) * cos_theta * p_prev1 - (n - 1) * p_prev2) / n
            p_prev2 = p_prev1
            p_prev1 = p_current
        end

        return p_prev1
    end
end

# =============================================================================
# PER-EPOCH CHANNEL REPAIR
# =============================================================================

"""
    repair_channels_per_epoch!(epochs::Vector{EpochData}, layout::Layout, threshold::Real, artifact_col::Symbol)

Repair channels on a per-epoch basis for epochs with artifacts.

This function identifies channels that exceed the artifact threshold in each epoch and repairs
them using neighbor interpolation if good neighbors are available. This allows trial-by-trial 
repair when channels have good neighbors.

# Arguments
- `epochs_list::Vector{EpochData}`: Vector of EpochData objects (one per condition)
- `layout::Layout`: Layout object containing neighbor information
- `threshold::Real`: Artifact detection threshold (e.g., 100.0 μV)
- `artifact_col::Symbol`: Column name for artifact flags (e.g., :is_artifact_value_100)

# Returns
- `Tuple{Vector{Symbol}, Int}`: (unique_channels_repaired, epochs_repaired_count)

# Examples
```julia
# Repair channels per epoch
channels_repaired, epochs_repaired = repair_channels_per_epoch!(
    epochs_cleaned, 
    layout, 
    100.0, 
    :is_artifact_value_100
)
```
"""
function repair_channels_per_epoch!(
    epochs_list::Vector{EpochData}, 
    layout::Layout, 
    threshold::Real, 
    artifact_col::Symbol
)::Tuple{Vector{Symbol}, Int}
    
    channels_with_artifacts = Symbol[]
    epochs_repaired = 0
    
    for (condition_idx, epoch_data) in enumerate(epochs_list)
        for epoch_idx in 1:length(epoch_data.data)
            epoch_df = epoch_data.data[epoch_idx]
            
            if hasproperty(epoch_df, artifact_col)
                if any(epoch_df[!, artifact_col])
                    @info "Condition $condition_idx, Epoch $epoch_idx: Artifacts detected"
                    
                    eeg_channels = channel_labels(epoch_data)
                    
                    # Find which specific channels have artifacts (exceed threshold) in this epoch
                    artifact_channels = Symbol[]
                    
                    for ch in eeg_channels
                        if hasproperty(epoch_df, ch)
                            ch_data = epoch_df[!, ch]
                            max_val = maximum(abs.(ch_data))
                            if max_val > threshold
                                push!(artifact_channels, ch)
                                @info "    - Channel $ch: max absolute value = $(round(max_val, sigdigits=3)), threshold = $threshold"
                            end
                        end
                    end
                    
                    if !isempty(artifact_channels)
                        # Check if any of these channels have good neighbors and can be repaired
                        repaired_in_epoch = Symbol[]
                        for ch in artifact_channels
                            if haskey(layout.neighbours, ch)
                                neighbors = layout.neighbours[ch]
                                good_neighbors = setdiff(neighbors.channels, artifact_channels)
                                if !isempty(good_neighbors)
                                    @info "    - Channel $ch: Has $(length(good_neighbors)) good neighbors: $(good_neighbors), attempting repair"
                                    push!(repaired_in_epoch, ch)
                                    # Repair this channel for this specific epoch
                                    repair_channels!(epoch_data, [ch]; epoch_selection = epochs(epoch_idx))
                                    if ch ∉ channels_with_artifacts
                                        push!(channels_with_artifacts, ch)
                                    end
                                else
                                    @info "    - Channel $ch: No good neighbors available (all neighbors also have artifacts)"
                                end
                            else
                                @info "    - Channel $ch: No neighbor information available in layout"
                            end
                        end
                        
                        if !isempty(repaired_in_epoch)
                            epochs_repaired += 1
                        end
                    end
                end
            end
        end
    end
    
    @info "Per-epoch repair summary: Repaired $(length(channels_with_artifacts)) unique channels across $(epochs_repaired) epochs"
    @info "Channels repaired: $(channels_with_artifacts)"
    
    return (channels_with_artifacts, epochs_repaired)
end

