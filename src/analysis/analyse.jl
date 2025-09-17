# =============================================================================
# ERP ANALYSIS
# =============================================================================

"""
    get_mean_amplitude(erp_data::ErpData, time_window::Tuple{<:Real, <:Real})

Calculates the mean amplitude for each electrode within a specified time window.

# Arguments
- `erp_data::ErpData`: ERP data structure
- `time_window::Tuple{<:Real, <:Real}`: Time window as (start_time, end_time) in seconds

# Returns
- `DataFrame`: A DataFrame with electrode labels as column names and corresponding mean amplitudes

# Examples
```julia
# Calculate mean amplitude in the 100-200ms window
mean_amp = get_mean_amplitude(erp_data, (0.1, 0.2))

# Calculate mean amplitude in the 300-500ms window
mean_amp = get_mean_amplitude(erp_data, (0.3, 0.5))
```
"""
function get_mean_amplitude(erp_data::ErpData, time_window::Tuple{<:Real,<:Real})
    # Find time indices within the window
    time_indices = findall(x -> x >= time_window[1] && x <= time_window[2], erp_data.time)

    if isempty(time_indices)
        error("No data points found within the specified time window")
    end

    # Calculate mean amplitude for each electrode
    mean_amplitudes = Dict{Symbol,Float64}()
    for electrode in erp_data.layout.label
        if haskey(erp_data.data, electrode)
            mean_amplitudes[electrode] = mean(erp_data.data[time_indices, electrode])
        end
    end

    return DataFrame(mean_amplitudes)
end
