"""
    repair_artifacts!(dat::EpochData, artifacts::EpochRejectionInfo; method::Symbol=:neighbor_interpolation, kwargs...)

Repair detected artifacts using the specified method.

# Arguments
- `dat::EpochData`: The epoch data to repair (modified in-place)
- `artifacts::EpochRejectionInfo`: Artifact information from detect_artifacts
- `method::Symbol`: Repair method to use

# Available Methods
- `:neighbor_interpolation` - Weighted neighbor interpolation (default). Uses `dat.layout.neighbours` for neighbor information.
- `:spherical_spline` - Spherical spline interpolation

# Keyword Arguments (for :spherical_spline method)
- `m::Int`: Order of Legendre polynomials (default: 4)
- `lambda::Float64`: Regularization parameter (default: 1e-5)

# Returns
Nothing (mutates dat in-place)
"""
function repair_artifacts!(dat::EpochData, artifacts::EpochRejectionInfo; method::Symbol = :neighbor_interpolation, kwargs...)

    if method ∉ [:neighbor_interpolation, :spherical_spline]
        throw(ArgumentError("Unknown repair method: $method. Available: :neighbor_interpolation, :spherical_spline"))
    end

    @info "--------------------------------"
    @info "Condition: $(dat.condition) ($(dat.condition_name)) - Repairing artifacts using method: $method"
    if method == :neighbor_interpolation
        # Determine which channels can be repaired if not already done
        if isnothing(artifacts.repaired)
            channel_repairable!(artifacts, dat.layout)
        end
        repair_artifacts_neighbor!(dat, artifacts; kwargs...)
    elseif method == :spherical_spline
        repair_artifacts_spherical_spline!(dat, artifacts; kwargs...)
    end
    return nothing
end


"""
    repair_artifacts(dat::EpochData, artifacts::EpochRejectionInfo; method::Symbol=:neighbor_interpolation, kwargs...)

Non-mutating version of repair_artifacts!. Creates a copy of the data and repairs artifacts without modifying the original.

# Arguments
- `dat::EpochData`: The epoch data to repair (NOT modified)
- `artifacts::EpochRejectionInfo`: Artifact information from detect_artifacts
- `method::Symbol`: Repair method to use (default: :neighbor_interpolation)

# Keyword Arguments
Same as repair_artifacts! for the respective methods.

# Returns
- `EpochData`: A new EpochData object with repaired artifacts

# Examples
```julia
# Basic artifact repair (creates new object)
repaired_epochs = repair_artifacts(epochs, artifacts)

# Using spherical spline method
repaired_epochs = repair_artifacts(epochs, artifacts, :spherical_spline)

# Rejecting bad epochs entirely (use reject_epochs instead)
clean_epochs = reject_epochs(epochs, artifacts)
```
"""
function repair_artifacts(dat::EpochData, artifacts::EpochRejectionInfo; method::Symbol = :neighbor_interpolation, kwargs...)
    dat_copy = copy(dat)
    repair_artifacts!(dat_copy, artifacts; method, kwargs...)
    return dat_copy
end

function repair_artifacts(dat::Vector{EpochData}, artifacts::Vector{EpochRejectionInfo}; kwargs...)
    return repair_artifacts.(dat, artifacts; kwargs...)
end

function repair_artifacts!(dat::Vector{EpochData}, artifacts::Vector{EpochRejectionInfo}; kwargs...)
    repair_artifacts!.(dat, artifacts; kwargs...)
    return nothing
end



"""
    channel_repairable!(artifacts::EpochRejectionInfo, layout::Layout)

Analyze which channels can be repaired and which cannot, based on neighbor availability.
Populates `artifacts.repaired` and `artifacts.skipped` with the analysis.

# Arguments
- `artifacts::EpochRejectionInfo`: Artifact information from detect_artifacts (mutated to add repair analysis)
- `layout::Layout`: Layout object containing neighbor information

# Returns
- `EpochRejectionInfo`: The same artifacts object (modified in-place)

# Notes
Only analyzes repairability - does not perform any repairs.
Use `repair_artifacts_neighbor!` to actually perform the repairs after this analysis.
"""
function channel_repairable!(artifacts::EpochRejectionInfo, layout::Layout)
    rejected = unique([r.epoch for r in artifacts.rejected])

    # Initialize tracking dictionaries in artifacts struct (OrderedDict to maintain sorted order)
    artifacts.repaired = OrderedDict{Int,Vector{Symbol}}()
    artifacts.skipped = OrderedDict{Int,Vector{Symbol}}()

    # Process epochs in sorted order to maintain ordering in OrderedDict
    for epoch_idx in sort(rejected)
        bad_channels = [artifact.channel for artifact in artifacts.rejected if artifact.epoch == epoch_idx]
        isempty(bad_channels) && continue

        repairable_channels = check_channel_neighbors(bad_channels, layout)

        if isempty(repairable_channels)
            if length(bad_channels) == 1
                @info "Epoch $epoch_idx: Cannot repair channel $(bad_channels[1]) (fewer than 2 neighbors)"
            else
                @info "Epoch $epoch_idx: Cannot repair channels $(bad_channels) (bad neighbors and/or fewer than 2 neighbors)"
            end
            # Track that all were skipped
            artifacts.skipped[epoch_idx] = bad_channels
        else
            skipped_channels = setdiff(bad_channels, repairable_channels)
            artifacts.repaired[epoch_idx] = repairable_channels
            if !isempty(skipped_channels)
                @info "Epoch $epoch_idx: Skipping repair of $(length(skipped_channels)) channel(s) with bad neighbors: $skipped_channels"
                artifacts.skipped[epoch_idx] = skipped_channels
            end
        end
    end

    return artifacts
end

channel_repairable!(artifacts::Vector{EpochRejectionInfo}, layout::Layout) = channel_repairable!.(artifacts, Ref(layout))


"""
    repair_artifacts_neighbor!(dat::EpochData, artifacts::EpochRejectionInfo)

Repair artifacts using weighted neighbor interpolation.
Uses `artifacts.repaired` to determine which channels to repair (should be populated by `channel_repairable!`).
Uses `dat.layout.neighbours` for neighbor information.

# Arguments
- `dat::EpochData`: The epoch data to repair (modified in-place)
- `artifacts::EpochRejectionInfo`: Artifact information with `repaired` already populated

# Returns
- `EpochData`: The repaired epoch data (same object, modified in-place)

# See also
- `channel_repairable!`: Analyze which channels can be repaired before calling this function
"""
function repair_artifacts_neighbor!(dat::EpochData, artifacts::EpochRejectionInfo)
    # Check if repaired has been populated
    if isnothing(artifacts.repaired)
        throw(ArgumentError("repaired not populated. Call channel_repairable! first."))
    end

    if isempty(artifacts.repaired)
        @info "No channels to repair (all bad channels were skipped)"
        return dat
    end

    # Process epochs in sorted order (already sorted in OrderedDict)
    for (epoch_idx, repairable_channels) in artifacts.repaired
        @info "Epoch $epoch_idx: Repairing channels $(repairable_channels) using neighbor interpolation"

        # Use unified channel repair function with epoch selection
        repair_channels!(
            dat,
            repairable_channels;
            method = :neighbor_interpolation,
            epoch_selection = epochs([epoch_idx]),
            neighbours_dict = dat.layout.neighbours,
        )
        @info "" # formatting
    end

    return dat
end

"""
    repair_artifacts_spherical_spline!(dat::EpochData, artifacts::EpochRejectionInfo; m::Int=4, lambda::Real=1e-5)

Repair artifacts using spherical spline interpolation.

# Arguments
- `dat::EpochData`: The epoch data to repair (modified in-place)
- `artifacts::EpochRejectionInfo`: Artifact information from detect_artifacts
- `m::Int`: Order of Legendre polynomials (default: 4)
- `lambda::Real`: Regularization parameter (default: 1e-5)

# Returns
- `EpochData`: The repaired epoch data (same object, modified in-place)
"""
function repair_artifacts_spherical_spline!(dat::EpochData, artifacts::EpochRejectionInfo; m::Int = 4, lambda::Real = 1e-5)

    _ensure_coordinates_3d!(dat.layout)

    # Get all rejected epochs with their bad channels
    rejected = unique([r.epoch for r in artifacts.rejected])

    for epoch_idx in rejected
        # Get bad channels for this epoch
        bad_channels = [artifact.channel for artifact in artifacts.rejected if artifact.epoch == epoch_idx]
        isempty(bad_channels) && continue

        @info "Repairing epoch $epoch_idx channels $(bad_channels) using spherical spline interpolation"

        # Use unified channel repair function with epoch selection
        repair_channels!(dat, bad_channels; method = :spherical_spline, epoch_selection = epochs([epoch_idx]), m = m, lambda = lambda)
    end

    return dat
end



