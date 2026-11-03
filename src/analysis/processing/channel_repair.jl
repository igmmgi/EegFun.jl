# === UNIFIED CHANNEL REPAIR INTERFACE ===
_to_repair_channels(ch::Symbol) = [ch]
_to_repair_channels(ch::AbstractString) = [Symbol(ch)]
_to_repair_channels(chans::AbstractVector{Symbol}) = collect(chans)
_to_repair_channels(chans::AbstractVector{<:AbstractString}) = Symbol.(chans)
_to_repair_channels(chans) = collect(chans)

"""
    repair_channels!(data, channels_to_repair; method=:neighbor_interpolation, kwargs...)

Repair bad EEG channels using the specified method.

# Arguments
- `data`: EEG data object (`ContinuousData`, `ErpData`, `EpochData`, or `AbstractMatrix`)
- `channels_to_repair`: Channel label or list of channel labels to repair

# Keyword Arguments
- `method::Symbol`: Repair method to use
  - `:neighbor_interpolation` (default): Weighted neighbor interpolation
  - `:spherical_spline`: Spherical spline interpolation
- `neighbors_dict::Union{AbstractDict, Nothing}`: Neighbor information (for `neighbor_interpolation`)
- `m::Int`: Order of Legendre polynomials (for `spherical_spline`, default: 4)
- `lambda::Float64`: Regularization parameter (for `spherical_spline`, default: 1e-5)
- `epoch_selection::Function`: Epoch selection predicate (for `EpochData`, default: all epochs)

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
function repair_channels!(data, channels_to_repair; method::Symbol = :neighbor_interpolation, repair_info = nothing, kwargs...)
    channels_to_repair = _to_repair_channels(channels_to_repair)
    # repair_info is only used for ContinuousData/ErpData, not for EpochData (tracking is in EpochRejectionInfo)
    actually_repaired = if method == :neighbor_interpolation
        _repair_channels_neighbor!(data, channels_to_repair; repair_info = repair_info, kwargs...)
    elseif method == :spherical_spline
        if data isa EpochData
            _repair_channels_spherical!(data, channels_to_repair; kwargs...)
        else
            _repair_channels_spherical!(data, channels_to_repair; repair_info = repair_info, kwargs...)
        end
    else
        @minimal_error("Unknown repair method: $method. Use :neighbor_interpolation or :spherical_spline")
    end

    if !(data isa EpochData) && hasfield(typeof(data), :analysis_info) && !isnothing(actually_repaired)
        for ch in actually_repaired
            if ch ∉ data.analysis_info.repaired_channels
                push!(data.analysis_info.repaired_channels, ch)
            end
        end
        data.analysis_info.n_channels_repaired = length(data.analysis_info.repaired_channels)
    end

    return data
end

"""
    repair_channels(dat, args...; kwargs...)

Non-mutating version of `repair_channels!`.
"""
function repair_channels(dat, args...; kwargs...)
    return repair_channels!(copy(dat), args...; kwargs...)
end

function repair_channels!(dat::AbstractVector{<:EegData}, args...; kwargs...)
    for d in dat
        repair_channels!(d, args...; kwargs...)
    end
    return dat
end

function repair_channels(dat::AbstractVector{<:EegData}, args...; kwargs...)
    return [repair_channels(d, args...; kwargs...) for d in dat]
end

# Helper function for neighbor interpolation
function _repair_channels_neighbor!(data, channels_to_repair; neighbors_dict = nothing, repair_info = nothing, kwargs...)
    if isnothing(neighbors_dict)
        # For data objects, get neighbors from layout
        if hasfield(typeof(data), :layout)
            layout = data.layout
            if isnothing(layout.neighbors)
                get_neighbors_xyz!(layout, 0.5)
            end
            neighbors_dict = layout.neighbors
        else
            @minimal_error("neighbors_dict must be provided for AbstractMatrix data")
        end
    end

    # Use multiple dispatch - repair_info only used for ContinuousData/ErpData, not for EpochData
    if data isa EpochData
        return _repair_channels_neighbor!(data, channels_to_repair, neighbors_dict; kwargs...)
    else
        return _repair_channels_neighbor!(data, channels_to_repair, neighbors_dict; repair_info = repair_info, kwargs...)
    end
end

# === CORE IMPLEMENTATION FUNCTIONS ===

"""
    _repair_channels_neighbor!(data::AbstractMatrix, channels::AbstractVector{Symbol}, channels_to_repair::AbstractVector{Symbol}, neighbors_dict::AbstractDict{Symbol, Neighbors})

Repair bad channels in a data matrix using weighted neighbor interpolation.
Excludes all channels in `channels_to_repair` (and optional `bad_channels`) from each channel's
neighbor set, re-normalizes remaining weights, and writes interpolated traces into a temporary
matrix before copying back so repair is order-independent.
"""
function _repair_channels_neighbor!(
    data::AbstractMatrix,
    channels::AbstractVector{Symbol},
    channels_to_repair::AbstractVector{Symbol},
    neighbors_dict::AbstractDict{Symbol,Neighbors};
    repair_info = nothing,
    bad_channels::Union{Nothing,AbstractVector{Symbol}} = nothing,
)
    # Check that channels match data dimensions
    n_channels = size(data, 1)
    if length(channels) != n_channels
        @minimal_error "Number of channels ($(length(channels))) doesn't match data dimensions ($n_channels)"
    end

    # Initialize tracking if provided
    if !isnothing(repair_info)
        repair_info.repaired = Symbol[]
        repair_info.skipped = Symbol[]
    end

    unique_to_repair = unique(channels_to_repair)
    isempty(unique_to_repair) && return Symbol[]

    # Exclude all bad channels from neighbor candidates
    excluded_set = Set{Symbol}(unique_to_repair)
    if !isnothing(bad_channels)
        union!(excluded_set, bad_channels)
    end

    # Create channel index lookup
    ch_indices = Dict(ch => i for (i, ch) in enumerate(channels))

    n_timepoints = size(data, 2)
    temp_traces = Matrix{Float64}(undef, length(unique_to_repair), n_timepoints)
    repaired_slots = Tuple{Int,Int,Symbol}[]

    # Compute interpolated traces into temporary buffer without mutating `data`
    for bad_ch in unique_to_repair
        bad_idx = get(ch_indices, bad_ch, nothing)
        if isnothing(bad_idx)
            @minimal_warning "Channel $bad_ch not found in channel list, skipping"
            !isnothing(repair_info) && push!(repair_info.skipped, bad_ch)
            continue
        end

        neighbors = get(neighbors_dict, bad_ch, nothing)
        if isnothing(neighbors) || isempty(neighbors.channels)
            @minimal_warning "No neighbors found for channel $bad_ch, skipping"
            !isnothing(repair_info) && push!(repair_info.skipped, bad_ch)
            continue
        end

        clean_neigh_channels = Symbol[]
        clean_neigh_indices = Int[]
        clean_neigh_weights = Float64[]
        for (i, neighbor_ch) in enumerate(neighbors.channels)
            neighbor_ch ∈ excluded_set && continue
            neigh_idx = get(ch_indices, neighbor_ch, nothing)
            if !isnothing(neigh_idx)
                push!(clean_neigh_channels, neighbor_ch)
                push!(clean_neigh_indices, neigh_idx)
                push!(clean_neigh_weights, Float64(neighbors.weights[i]))
            end
        end

        if isempty(clean_neigh_indices)
            @minimal_warning "No clean neighbors available for channel $bad_ch, skipping"
            !isnothing(repair_info) && push!(repair_info.skipped, bad_ch)
            continue
        end

        weight_sum = sum(clean_neigh_weights)
        if weight_sum > 0.0
            clean_neigh_weights ./= weight_sum
        else
            fill!(clean_neigh_weights, 1.0 / length(clean_neigh_weights))
        end

        slot = length(repaired_slots) + 1
        @inbounds for t = 1:n_timepoints
            weighted_sum = 0.0
            for k in eachindex(clean_neigh_indices)
                weighted_sum += Float64(data[clean_neigh_indices[k], t]) * clean_neigh_weights[k]
            end
            temp_traces[slot, t] = weighted_sum
        end

        push!(repaired_slots, (slot, bad_idx, bad_ch))
        @info "Repaired channel $bad_ch using weighted neighbor interpolation with neighbors: $(clean_neigh_channels)"
    end

    # Copy interpolated traces back into `data` after all repairs have been computed
    actually_repaired = Symbol[]
    for (slot, bad_idx, bad_ch) in repaired_slots
        @inbounds for t = 1:n_timepoints
            data[bad_idx, t] = temp_traces[slot, t]
        end
        push!(actually_repaired, bad_ch)
        !isnothing(repair_info) && push!(repair_info.repaired, bad_ch)
    end

    return actually_repaired
end

"""
    _repair_channels_neighbor!(data::Union{ContinuousData,ErpData}, channels_to_repair::AbstractVector{Symbol}, neighbors_dict::AbstractDict{Symbol, Neighbors}; repair_info = nothing)

Repair bad channels in ContinuousData or ErpData using weighted neighbor interpolation.
"""
function _repair_channels_neighbor!(
    data::Union{ContinuousData,ErpData},
    channels_to_repair::AbstractVector{Symbol},
    neighbors_dict::AbstractDict{Symbol,Neighbors};
    repair_info = nothing,
    bad_channels::Union{Nothing,AbstractVector{Symbol}} = nothing,
)
    channels = Symbol.(data.layout.data.label)
    data_matrix = Matrix{Float64}(Matrix(data.data[:, channels])')

    actually_repaired = _repair_channels_neighbor!(
        data_matrix,
        channels,
        channels_to_repair,
        neighbors_dict;
        repair_info = repair_info,
        bad_channels = bad_channels,
    )

    data.data[:, channels] = data_matrix'

    return actually_repaired
end

"""
    _repair_channels_neighbor!(data::EpochData, channels_to_repair::AbstractVector{Symbol}, neighbors_dict::AbstractDict{Symbol, Neighbors}; epoch_selection::Function=epochs())

Repair bad channels in epoched EEG data using weighted neighbor interpolation.
Tracking for EpochData is handled via EpochRejectionInfo, not through a separate repair_info parameter.
"""
function _repair_channels_neighbor!(
    data::EpochData,
    channels_to_repair::AbstractVector{Symbol},
    neighbors_dict::AbstractDict{Symbol,Neighbors};
    epoch_selection = epochs(),
    bad_channels::Union{Nothing,AbstractVector{Symbol}} = nothing,
)
    selected_epochs = get_selected_epochs(data, epoch_selection)
    channels = Symbol.(data.layout.data.label)
    actually_repaired = Symbol[]

    for epoch_idx in selected_epochs
        epoch = data.data[epoch_idx]
        epoch_matrix = Matrix{Float64}(Matrix(epoch[:, channels])')

        rep = _repair_channels_neighbor!(
            epoch_matrix,
            channels,
            channels_to_repair,
            neighbors_dict;
            repair_info = nothing,
            bad_channels = bad_channels,
        )
        for ch in rep
            ch ∉ actually_repaired && push!(actually_repaired, ch)
        end

        epoch[:, channels] = epoch_matrix'
    end

    return actually_repaired
end

"""
    _repair_channels_spherical!(data::AbstractMatrix, channels_to_repair::AbstractVector{Symbol}, channels::AbstractVector{Symbol}, layout::DataFrame; m::Int=4, lambda::Float64=1e-5)

Core implementation for repairing bad channels using spherical spline interpolation.
"""
function _repair_channels_spherical!(
    data::AbstractMatrix,
    channels_to_repair::AbstractVector{Symbol},
    channels::AbstractVector{Symbol},
    layout::DataFrame;
    m::Int = 4,
    lambda::Real = 1e-5,
    n_legendre_terms::Int = 50,
    repair_info = nothing,
)
    # Check that channels match data dimensions
    n_channels = size(data, 1)
    if length(channels) != n_channels
        @minimal_error("Number of channels doesn't match data dimensions")
    end

    # Initialize tracking if provided
    if !isnothing(repair_info)
        repair_info.repaired = Symbol[]
        repair_info.skipped = Symbol[]
    end

    unique_to_repair = unique(channels_to_repair)
    isempty(unique_to_repair) && return Symbol[]

    # Extract coordinates for all channels
    coords = zeros(Float64, length(channels), 3)
    for (i, ch) in enumerate(channels)
        idx = findfirst(x -> x == ch, layout.label)
        if isnothing(idx)
            @minimal_warning "Channel $ch not found in layout, using (0,0,0)"
            continue
        end
        coords[i, :] = [layout.x3[idx], layout.y3[idx], layout.z3[idx]]
    end

    # Find indices of bad and good channels
    bad_indices = Int[]
    found_bad_channels = Symbol[]
    for bad_ch in unique_to_repair
        bad_idx = findfirst(x -> x == bad_ch, channels)
        if isnothing(bad_idx)
            @minimal_warning "Channel $bad_ch not found in channel list, skipping"
            !isnothing(repair_info) && push!(repair_info.skipped, bad_ch)
            continue
        end
        push!(bad_indices, bad_idx)
        push!(found_bad_channels, bad_ch)
    end

    if isempty(bad_indices)
        return Symbol[]
    end

    # Good channels are all channels except the bad ones
    good_indices = setdiff(1:length(channels), bad_indices)
    if isempty(good_indices)
        @minimal_warning "Cannot repair channels via spherical spline: no clean channels remain, skipping"
        if !isnothing(repair_info)
            append!(repair_info.skipped, found_bad_channels)
        end
        return Symbol[]
    end

    # Extract positions
    pos_good = coords[good_indices, :]
    pos_bad = coords[bad_indices, :]

    # Calculate interpolation matrix once for all bad channels
    interpolation_matrix = _make_interpolation_matrix(pos_good, pos_bad; m = m, lambda = lambda, n_legendre_terms = n_legendre_terms)

    # Apply interpolation into a temporary matrix before copying back
    repaired_traces = interpolation_matrix * Float64.(data[good_indices, :])
    data[bad_indices, :] = repaired_traces

    # Log and track repairs
    for bad_ch in found_bad_channels
        @info "Repaired channel $bad_ch using spherical spline interpolation"
        !isnothing(repair_info) && push!(repair_info.repaired, bad_ch)
    end

    return found_bad_channels
end

"""
    _repair_channels_spherical!(data::Union{ContinuousData,ErpData}, channels_to_repair::AbstractVector{Symbol}; m::Int=4, lambda::Float64=1e-5)

Repair bad channels in ContinuousData or ErpData using spherical spline interpolation.
"""
function _repair_channels_spherical!(
    data::Union{ContinuousData,ErpData},
    channels_to_repair::AbstractVector{Symbol};
    m::Int = 4,
    lambda::Real = 1e-5,
    repair_info = nothing,
)
    # Ensure 3D coordinates are available
    _ensure_coordinates_3d!(data.layout)

    channels = Symbol.(data.layout.data.label)
    data_matrix = Matrix{Float64}(Matrix(data.data[:, channels])')

    actually_repaired = _repair_channels_spherical!(
        data_matrix,
        channels_to_repair,
        channels,
        data.layout.data;
        m = m,
        lambda = lambda,
        repair_info = repair_info,
    )

    # Update the data
    data.data[:, channels] = data_matrix'

    return actually_repaired
end

"""
    _repair_channels_spherical!(data::EpochData, channels_to_repair::AbstractVector{Symbol}; epoch_selection::Function=epochs(), m::Int=4, lambda::Float64=1e-5)

Repair bad channels in epoched EEG data using spherical spline interpolation.
Tracking for EpochData is handled via EpochRejectionInfo, not through a separate repair_info parameter.
"""
function _repair_channels_spherical!(
    data::EpochData,
    channels_to_repair::AbstractVector{Symbol};
    epoch_selection = epochs(),
    m::Int = 4,
    lambda::Real = 1e-5,
)
    # Ensure 3D coordinates are available
    _ensure_coordinates_3d!(data.layout)

    # Get selected epochs
    selected_epochs = get_selected_epochs(data, epoch_selection)

    # Get channels vector
    channels = Symbol.(data.layout.data.label)
    actually_repaired = Symbol[]

    # Repair each selected epoch
    for epoch_idx in selected_epochs
        epoch = data.data[epoch_idx]
        epoch_matrix = Matrix{Float64}(Matrix(epoch[:, channels])')

        # Note: repair_info is ignored for EpochData - tracking is handled via EpochRejectionInfo
        rep = _repair_channels_spherical!(
            epoch_matrix,
            channels_to_repair,
            channels,
            data.layout.data;
            m = m,
            lambda = lambda,
            repair_info = nothing,
        )
        for ch in rep
            ch ∉ actually_repaired && push!(actually_repaired, ch)
        end

        # Update the epoch data
        epoch[:, channels] = epoch_matrix'
    end

    return actually_repaired
end

# === SPHERICAL SPLINE UTILITIES (Perrin et al., 1989) ===

"""Factors for Legendre polynomial series in spherical spline interpolation."""
_g_factors(m::Int, nterms::Int) = [(2n + 1) / (float(n)^m * float(n + 1)^m * 4π) for n = 1:nterms]

"""Evaluate Legendre polynomial series using forward 3-term recurrence."""
@inline function _legendre_series(x::Float64, f::Vector{Float64})
    x = clamp(x, -1.0, 1.0)
    p_prev, p = 1.0, x
    acc = f[1] * x
    @inbounds for n = 2:length(f)
        p_next = ((2n - 1) * x * p - (n - 1) * p_prev) / n
        acc += f[n] * p_next
        p_prev, p = p, p_next
    end
    return acc
end

"""In-place computation of G matrix from cosine matrix and factors."""
function _g_from_cos!(C::Matrix{Float64}, f::Vector{Float64})
    @inbounds for i in eachindex(C)
        C[i] = _legendre_series(C[i], f)
    end
    return C
end

"""
    _spline_geometry(d, ch::Vector{Int})

Compute 3D unit-sphere vectors for channels in a layout.
Prioritizes :inc and :azi spherical coordinates; falls back to 2D equidistant projection.
"""
function _spline_geometry(d, ch::Vector{Int})
    ρe = hypot.(d.x2[ch], d.y2[ch])
    if hasproperty(d, :inc) && hasproperty(d, :azi)
        rad_factor = pi / 180.0
        inc = Float64.(d.inc[ch]) .* rad_factor
        azi = Float64.(d.azi[ch]) .* rad_factor
        U = hcat(sin.(inc) .* cos.(azi), sin.(inc) .* sin.(azi), cos.(inc))
        θ = abs.(inc)
        denom = dot(ρe, ρe)
        k = denom > 0 ? dot(ρe, θ) / denom : (π / 2)
        if maximum(abs, θ .- k .* ρe) > 1e-3
            @warn "x2/y2 are not an azimuthal-equidistant image of inc/azi; spline map may be distorted." maxlog = 1
        end
    else
        # Fallback: only 2D coords available. Assume azimuthal-equidistant x2/y2 with ρ=1 at equator.
        @warn "Layout has no :inc/:azi; assuming azimuthal-equidistant x2/y2 with ρ=1 at the equator." maxlog = 1
        k = π / 2
        θ = min.(k .* ρe, π)
        sx = [ρe[i] > 0 ? sin(θ[i]) / ρe[i] : 0.0 for i = 1:length(ρe)]
        U = hcat(sx .* d.x2[ch], sx .* d.y2[ch], cos.(θ))
    end
    return U, k
end

"""
    _spline_good_to_bad(U_good::Matrix{Float64}, U_bad::Matrix{Float64}, f::Vector{Float64}, λ::Float64)
                        -> Matrix{Float64}

Construct spherical spline interpolation matrix (n_bad × n_good) from unit vectors on the sphere
using LU factorization (Perrin et al., 1989).
"""
function _spline_good_to_bad(U_good::Matrix{Float64}, U_bad::Matrix{Float64}, f::Vector{Float64}, λ::Real)
    n = size(U_good, 1)
    A = zeros(Float64, n + 1, n + 1)
    A[1:n, 1:n] .= _g_from_cos!(U_good * U_good', f)
    @inbounds for i = 1:n
        A[i, i] += λ
        A[i, n+1] = 1.0
        A[n+1, i] = 1.0
    end
    A[n+1, n+1] = 0.0
    X = lu(A) \ [Matrix{Float64}(I, n, n); zeros(Float64, 1, n)]
    Gb = _g_from_cos!(U_bad * U_good', f)
    return Gb * view(X, 1:n, :) .+ view(X, (n+1):(n+1), :)
end

"""
    _calc_g(cosang::AbstractArray, m::Int=4, n_legendre_terms::Int=50)

Calculate spherical spline G function between points on a sphere using Legendre polynomials.
"""
function _calc_g(cosang::AbstractArray, m::Int = 4, n_legendre_terms::Int = 50)
    f = _g_factors(m, n_legendre_terms)
    C = Float64.(cosang)
    return _g_from_cos!(C, f)
end

"""
    _make_interpolation_matrix(pos_from, pos_to, m=4, lambda=1e-5, n_legendre_terms=50)

Create spherical spline interpolation matrix using Perrin et al. 1989 algorithm.

# Arguments
- `pos_from`: Positions of good channels (n_good × 3)
- `pos_to`: Positions of bad channels to interpolate (n_bad × 3)
- `m`: Stiffness parameter
- `lambda`: Regularization parameter
- `n_legendre_terms`: Number of Legendre terms

# Returns
- Interpolation matrix (n_bad × n_good)
"""
function _make_interpolation_matrix(
    pos_from::Matrix{Float64},
    pos_to::Matrix{Float64};
    m::Int = 4,
    lambda::Real = 1e-5,
    n_legendre_terms::Int = 50,
)
    n_from = size(pos_from, 1)
    n_to = size(pos_to, 1)

    # Normalize positions to unit sphere
    pos_from_norm = copy(pos_from)
    pos_to_norm = copy(pos_to)

    for i = 1:n_from
        norm_v = sqrt(sum(pos_from_norm[i, :] .^ 2))
        if norm_v > 1e-12
            pos_from_norm[i, :] ./= norm_v
        else
            pos_from_norm[i, :] .= (0.0, 0.0, 1.0)
        end
    end

    for i = 1:n_to
        norm_v = sqrt(sum(pos_to_norm[i, :] .^ 2))
        if norm_v > 1e-12
            pos_to_norm[i, :] ./= norm_v
        else
            pos_to_norm[i, :] .= (0.0, 0.0, 1.0)
        end
    end

    f = _g_factors(m, n_legendre_terms)
    return _spline_good_to_bad(pos_from_norm, pos_to_norm, f, lambda)
end



# === PER-EPOCH CHANNEL REPAIR ===

"""
    repair_channels_per_epoch!(epochs::Vector{EpochData}, layout::Layout, threshold::Real, artifact_col::Symbol)

Repair channels on a per-epoch basis for epochs with artifacts.

Identify channels that exceed the artifact threshold in each epoch and repair
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
    epochs_unrejected, 
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
    artifact_col::Symbol,
)::Tuple{Vector{Symbol},Int}

    channels_with_artifacts = Symbol[]
    epochs_repaired = 0

    for (condition_idx, epoch_data) in enumerate(epochs_list)
        for epoch_idx = 1:n_epochs(epoch_data)
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
                            max_val = maximum(abs, ch_data)
                            if max_val > threshold
                                push!(artifact_channels, ch)
                                @info "    - Channel $ch: max absolute value = $(round(max_val, sigdigits=3)), threshold = $threshold"
                            end
                        end
                    end

                    if !isempty(artifact_channels)
                        # Check if any of these channels have good neighbors and can be repaired
                        repaired_in_epoch = Symbol[]
                        if !isnothing(layout.neighbors)
                            for ch in artifact_channels
                                if haskey(layout.neighbors, ch)
                                    neighbors = layout.neighbors[ch]
                                    good_neighbors = setdiff(neighbors.channels, artifact_channels)
                                    if !isempty(good_neighbors)
                                        @info "    - Channel $ch: Has $(length(good_neighbors)) good neighbors: $(good_neighbors), attempting repair"
                                        push!(repaired_in_epoch, ch)
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
                        end

                        if !isempty(repaired_in_epoch)
                            repair_channels!(
                                epoch_data,
                                repaired_in_epoch;
                                epoch_selection = epochs(epoch_idx),
                                neighbors_dict = layout.neighbors,
                                bad_channels = artifact_channels,
                            )
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
