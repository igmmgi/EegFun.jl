# Pipeline bad-channel interpolation dispatchers.

######################################
########### User Functions ###########
######################################

"""
    repair_channels!(data::ContinuousData, repair_info::ContinuousRepairInfo; method::Symbol = repair_info.method, kwargs...) -> ContinuousData

Repair channels listed in `repair_info.repaired` (populated via [`channel_repairable!`](@ref)) in place using `:neighbor_interpolation` or `:spherical_spline`.
"""
function repair_channels!(data::ContinuousData, repair_info::ContinuousRepairInfo; method::Symbol = repair_info.method, kwargs...)
    method in (:neighbor_interpolation, :spherical_spline) ||
        @minimal_error "Unknown repair method: $method. Available: :neighbor_interpolation, :spherical_spline"
    if isempty(repair_info.repaired)
        @info "No channels to repair (all bad channels were skipped)"
        return data
    end
    prior_skipped = copy(repair_info.skipped)
    method_label = method == :neighbor_interpolation ? "neighbor interpolation" : "spherical spline interpolation"
    @info "Repairing $(length(repair_info.repaired)) channels: $(repair_info.repaired) using $method_label"
    repair_channels!(data, repair_info.repaired; method = method, repair_info = repair_info, bad_channels = prior_skipped, kwargs...)
    for ch in prior_skipped
        ch ∉ repair_info.repaired && ch ∉ repair_info.skipped && push!(repair_info.skipped, ch)
    end
    return data
end
