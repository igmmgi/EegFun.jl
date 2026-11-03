"""
    extract_ride(epochs::EpochData, comps::Vector{RideComponent}; kwargs...)

Extract RIDE components from `EpochData`. 
Returns a Dictionary mapping component names to `ErpData` objects containing the decomposed components
(each with independent copies of `layout` and `analysis_info`).
"""
function extract_ride(epochs::EpochData, comps::Vector{RideComponent}; kwargs...)
    isempty(epochs.data) && @minimal_error("EpochData contains no trials.")

    n_trials = n_epochs(epochs)
    n_samples = nrow(epochs.data[1])
    channels = epochs.layout.data.label
    n_channels = length(channels)
    time_vec = epochs.data[1].time

    # 1. Convert to 3D Tensor
    data_tensor = zeros(Float64, n_samples, n_channels, n_trials)

    for t = 1:n_trials
        df = epochs.data[t]
        for (c_idx, ch) in enumerate(channels)
            data_tensor[:, c_idx, t] .= df[!, Symbol(ch)]
        end
    end

    # 2. Run the RIDE engine (passing epoch timeline & sample_rate unless overridden in kwargs)
    ride_res = ride_decompose(data_tensor, comps; times = time_vec, sample_rate = epochs.sample_rate, kwargs...)

    # 3. Package back into ErpData objects with independent layout & analysis_info copies
    results = Dict{Symbol,ErpData}()
    for (i, comp_name) in enumerate(ride_res.components)
        df = DataFrame(time = copy(time_vec))
        for (c_idx, ch) in enumerate(channels)
            df[!, Symbol(ch)] = ride_res.templates[:, c_idx, i]
        end

        results[comp_name] = ErpData(
            epochs.file,
            epochs.condition,
            epochs.condition_name * "_RIDE_$(comp_name)",
            df,
            copy(epochs.layout),
            epochs.sample_rate,
            copy(epochs.analysis_info),
            n_trials,
        )
    end

    return results
end
