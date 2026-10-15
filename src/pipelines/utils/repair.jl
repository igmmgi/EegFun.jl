# === CONTINUOUS DATA REPAIR FUNCTIONS ===

"""
    repair_channels_neighbor!(data::ContinuousData, repair_info::ContinuousRepairInfo)

Repair channels using weighted neighbor interpolation.
Uses `repair_info.repaired` to determine which channels to repair (should be populated by `channel_repairable!`).
Updates `repair_info` during repair to track actual repairs vs skips.

# Arguments
- `data::ContinuousData`: Continuous EEG data to repair (modified in-place)
- `repair_info::ContinuousRepairInfo`: Repair tracking struct with `repaired` channels already populated

# Returns
- `Nothing`: All modifications are in-place

# See also
- `channel_repairable!`: Analyze which channels can be repaired before calling this function
"""
function repair_channels_neighbor!(data::ContinuousData, repair_info::ContinuousRepairInfo)
    if isempty(repair_info.repaired)
        @info "No channels to repair (all bad channels were skipped)"
        return nothing
    end

    @info "Repairing $(length(repair_info.repaired)) channels: $(repair_info.repaired) using neighbor interpolation"

    repair_channels!(data, repair_info.repaired; method = :neighbor_interpolation, repair_info = repair_info)

    return nothing
end

"""
    repair_channels_spherical!(data::ContinuousData, repair_info::ContinuousRepairInfo; m::Int=4, lambda::Real=1e-5)

Repair channels using spherical spline interpolation.
Uses `repair_info.repaired` to determine which channels to repair (should be populated by `channel_repairable!`).
Updates `repair_info` during repair to track actual repairs vs skips.

# Arguments
- `data::ContinuousData`: Continuous EEG data to repair (modified in-place)
- `repair_info::ContinuousRepairInfo`: Repair tracking struct with `repaired` channels already populated
- `m::Int`: Order of Legendre polynomials (default: 4)
- `lambda::Real`: Regularization parameter (default: 1e-5)

# Returns
- `Nothing`: All modifications are in-place

# See also
- `channel_repairable!`: Analyze which channels can be repaired before calling this function
"""
function repair_channels_spherical!(data::ContinuousData, repair_info::ContinuousRepairInfo; m::Int = 4, lambda::Real = 1e-5)
    if isempty(repair_info.repaired)
        @info "No channels to repair (all bad channels were skipped)"
        return nothing
    end

    @info "Repairing $(length(repair_info.repaired)) channels: $(repair_info.repaired) using spherical spline interpolation"

    repair_channels!(data, repair_info.repaired; method = :spherical_spline, repair_info = repair_info, m = m, lambda = lambda)

    return nothing
end

"""
    repair_channels!(data::ContinuousData, repair_info::ContinuousRepairInfo; method::Symbol=:neighbor_interpolation, kwargs...)

Repair channels using the information in `repair_info`.
Orchestrator function similar to `repair_artifacts!` for epoch data.

# Arguments
- `data::ContinuousData`: Continuous EEG data to repair (modified in-place)
- `repair_info::ContinuousRepairInfo`: Repair tracking struct with `repaired` channels already populated
- `method::Symbol`: Repair method to use (default: :neighbor_interpolation)

# Returns
- `Nothing`: All modifications are in-place

# Notes
- `repair_info.repaired` should be populated by `channel_repairable!` before calling this function
"""
function repair_channels!(data::ContinuousData, repair_info::ContinuousRepairInfo; method::Symbol = :neighbor_interpolation, kwargs...)
    if method == :neighbor_interpolation
        repair_channels_neighbor!(data, repair_info; kwargs...)
    elseif method == :spherical_spline
        repair_channels_spherical!(data, repair_info; kwargs...)
    else
        throw(ArgumentError("Unknown repair method: $method. Available: :neighbor_interpolation, :spherical_spline"))
    end
    return nothing
end
