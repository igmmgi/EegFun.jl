"""
    repair_artifacts!(dat::EpochData, artifacts::EpochRejectionInfo; method::Symbol=:neighbor_interpolation, kwargs...)

Repair detected epoch artifacts in-place using `:neighbor_interpolation` (weighted neighbor interpolation from `dat.layout.neighbors`)
or `:spherical_spline` (accepts `m::Int=4`, `lambda::Real=1e-5`). Updates `artifacts.repaired` and `artifacts.skipped`.
"""
function repair_artifacts!(dat::EpochData, artifacts::EpochRejectionInfo; method::Symbol = :neighbor_interpolation, kwargs...)
    if method ∉ (:neighbor_interpolation, :spherical_spline)
        @minimal_error "Unknown repair method: $method. Available: :neighbor_interpolation, :spherical_spline"
    end

    @info "--------------------------------"
    @info "Condition: $(dat.condition) ($(dat.condition_name)) - Repairing artifacts using method: $method"
    if method == :neighbor_interpolation
        isnothing(artifacts.repaired) && channel_repairable!(artifacts, dat.layout)
        repair_artifacts_neighbor!(dat, artifacts; kwargs...)
    else
        repair_artifacts_spherical_spline!(dat, artifacts; kwargs...)
    end
    return dat
end

"""
    repair_artifacts(dat::EpochData, artifacts::EpochRejectionInfo; method::Symbol=:neighbor_interpolation, kwargs...)

Non-mutating version of [`repair_artifacts!`](@ref); leaves both `dat` and `artifacts` unmodified.
"""
function repair_artifacts(dat::EpochData, artifacts::EpochRejectionInfo; method::Symbol = :neighbor_interpolation, kwargs...)
    return repair_artifacts!(copy(dat), copy(artifacts); method = method, kwargs...)
end

repair_artifacts(dat::AbstractVector{<:EpochData}, artifacts::AbstractVector{<:EpochRejectionInfo}; kwargs...) =
    repair_artifacts.(dat, artifacts; kwargs...)

function repair_artifacts!(dat::AbstractVector{<:EpochData}, artifacts::AbstractVector{<:EpochRejectionInfo}; kwargs...)
    for (d, a) in zip(dat, artifacts)
        repair_artifacts!(d, a; kwargs...)
    end
    return dat
end

"""
    channel_repairable!(artifacts::EpochRejectionInfo, layout::Layout)

Determine which rejected channels per epoch can be repaired from `layout` neighbor information,
populating `artifacts.repaired` and `artifacts.skipped` in-place.
"""
function channel_repairable!(artifacts::EpochRejectionInfo, layout::Layout)
    artifacts.repaired = OrderedDict{Int,Vector{Symbol}}()
    artifacts.skipped = OrderedDict{Int,Vector{Symbol}}()

    for epoch_idx in sort(unique(r.epoch for r in artifacts.rejected))
        bad_channels = [r.channel for r in artifacts.rejected if r.epoch == epoch_idx]
        isempty(bad_channels) && continue

        repairable_channels = check_channel_neighbors(bad_channels, layout)
        if isempty(repairable_channels)
            if length(bad_channels) == 1
                @info "Epoch $epoch_idx: Cannot repair channel $(bad_channels[1]) (fewer than 2 neighbors)"
            else
                @info "Epoch $epoch_idx: Cannot repair channels $(bad_channels) (bad neighbors and/or fewer than 2 neighbors)"
            end
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

function repair_artifacts_neighbor!(dat::EpochData, artifacts::EpochRejectionInfo)
    isnothing(artifacts.repaired) && @minimal_error("repaired not populated. Call channel_repairable! first.")
    if isempty(artifacts.repaired)
        @info "No channels to repair (all bad channels were skipped)"
        return dat
    end

    for (epoch_idx, repairable_channels) in artifacts.repaired
        @info "Epoch $epoch_idx: Repairing channels $(repairable_channels) using neighbor interpolation"
        repair_channels!(
            dat,
            repairable_channels;
            method = :neighbor_interpolation,
            epoch_selection = epochs([epoch_idx]),
            neighbors_dict = dat.layout.neighbors,
        )
        @info ""
    end
    return dat
end

function repair_artifacts_spherical_spline!(dat::EpochData, artifacts::EpochRejectionInfo; m::Int = 4, lambda::Real = 1e-5)
    _ensure_coordinates_3d!(dat.layout)
    artifacts.repaired = OrderedDict{Int,Vector{Symbol}}()
    artifacts.skipped = OrderedDict{Int,Vector{Symbol}}()
    all_channels = channel_labels(dat)

    for epoch_idx in sort(unique(r.epoch for r in artifacts.rejected))
        bad_channels = unique(r.channel for r in artifacts.rejected if r.epoch == epoch_idx)
        isempty(bad_channels) && continue

        present_bad = [ch for ch in bad_channels if ch in all_channels]
        missing_bad = [ch for ch in bad_channels if ch ∉ all_channels]
        if isempty(setdiff(all_channels, present_bad))
            @minimal_warning "Epoch $epoch_idx: Cannot repair channels $(bad_channels) using spherical spline interpolation (all channels are bad)"
            artifacts.skipped[epoch_idx] = bad_channels
            continue
        end

        if !isempty(present_bad)
            @info "Repairing epoch $epoch_idx channels $(present_bad) using spherical spline interpolation"
            repair_channels!(dat, present_bad; method = :spherical_spline, epoch_selection = epochs([epoch_idx]), m = m, lambda = lambda)
            artifacts.repaired[epoch_idx] = present_bad
        end
        !isempty(missing_bad) && (artifacts.skipped[epoch_idx] = missing_bad)
    end
    return dat
end
