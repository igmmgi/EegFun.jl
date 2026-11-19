function _merge_prior_skipped!(repair_info::ContinuousRepairInfo, prior_skipped::Vector{Symbol})
    for ch in prior_skipped
        ch ∉ repair_info.repaired && ch ∉ repair_info.skipped && push!(repair_info.skipped, ch)
    end
    return nothing
end

function repair_channels_neighbor!(data::ContinuousData, repair_info::ContinuousRepairInfo)
    if isempty(repair_info.repaired)
        @info "No channels to repair (all bad channels were skipped)"
        return data
    end
    prior_skipped = copy(repair_info.skipped)
    @info "Repairing $(length(repair_info.repaired)) channels: $(repair_info.repaired) using neighbor interpolation"
    repair_channels!(data, repair_info.repaired; method = :neighbor_interpolation, repair_info = repair_info, bad_channels = prior_skipped)
    _merge_prior_skipped!(repair_info, prior_skipped)
    return data
end

function repair_channels_spherical!(data::ContinuousData, repair_info::ContinuousRepairInfo; m::Int = 4, lambda::Real = 1e-5)
    if isempty(repair_info.repaired)
        @info "No channels to repair (all bad channels were skipped)"
        return data
    end
    prior_skipped = copy(repair_info.skipped)
    @info "Repairing $(length(repair_info.repaired)) channels: $(repair_info.repaired) using spherical spline interpolation"
    repair_channels!(
        data,
        repair_info.repaired;
        method = :spherical_spline,
        repair_info = repair_info,
        m = m,
        lambda = lambda,
        bad_channels = prior_skipped,
    )
    _merge_prior_skipped!(repair_info, prior_skipped)
    return data
end

"""
    repair_channels!(data::ContinuousData, repair_info::ContinuousRepairInfo; method::Symbol = repair_info.method, kwargs...) -> ContinuousData

Repair channels listed in `repair_info.repaired` (populated via [`channel_repairable!`](@ref)) in place using `:neighbor_interpolation` or `:spherical_spline`.
"""
function repair_channels!(data::ContinuousData, repair_info::ContinuousRepairInfo; method::Symbol = repair_info.method, kwargs...)
    if method == :neighbor_interpolation
        repair_channels_neighbor!(data, repair_info; kwargs...)
    elseif method == :spherical_spline
        repair_channels_spherical!(data, repair_info; kwargs...)
    else
        @minimal_error "Unknown repair method: $method. Available: :neighbor_interpolation, :spherical_spline"
    end
    return data
end
