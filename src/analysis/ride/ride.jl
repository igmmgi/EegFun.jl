"""
    extract_ride(epochs::EpochData, comps::Vector{RideComponent}; kwargs...)

Extract RIDE components from `EpochData`. 
Returns a Dictionary mapping component names to `ErpData` objects containing the decomposed components.
"""
function extract_ride(epochs::EpochData, comps::Vector{RideComponent}; kwargs...)
    isempty(epochs.data) && error("EpochData contains no trials.")

    n_trials = length(epochs.data)
    n_samples = nrow(epochs.data[1])
    channels = String.(epochs.layout.data.label)
    n_channels = length(channels)

    # 1. Convert to 3D Tensor
    data_tensor = zeros(Float64, n_samples, n_channels, n_trials)

    for t = 1:n_trials
        df = epochs.data[t]
        for (c_idx, ch) in enumerate(channels)
            data_tensor[:, c_idx, t] .= df[!, Symbol(ch)]
        end
    end

    # 2. Run the RIDE engine
    ride_res = ride_decompose(data_tensor, comps; kwargs...)

    # 3. Package back into ErpData objects
    results = Dict{Symbol,ErpData}()
    for (i, comp_name) in enumerate(ride_res.components)
        df = DataFrame(time = epochs.data[1].time)
        for (c_idx, ch) in enumerate(channels)
            df[!, Symbol(ch)] = ride_res.templates[:, c_idx, i]
        end

        results[comp_name] = ErpData(
            epochs.file,
            epochs.condition,
            epochs.condition_name * "_RIDE_$(comp_name)",
            df,
            epochs.layout,
            epochs.sample_rate,
            epochs.analysis_info,
            n_trials,
        )
    end

    return results
end
