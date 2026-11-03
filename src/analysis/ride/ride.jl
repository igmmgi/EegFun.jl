@inline function _copy_col_to_tensor!(tensor::Array{Float64,3}, c_idx::Int, t_idx::Int, col::AbstractVector)
    @inbounds for s in axes(tensor, 1)
        tensor[s, c_idx, t_idx] = Float64(col[s])
    end
    return nothing
end

"""
    extract_ride(epochs::EpochData, comps::Vector{RideComponent}; kwargs...) -> Dict{Symbol, ErpData}

Extract RIDE components from `EpochData`, returning a `Dict{Symbol, ErpData}` mapping each
component name to its decomposed template (with independent copies of `layout` and `analysis_info`).
"""
function extract_ride(epochs::EpochData, comps::Vector{RideComponent}; kwargs...)
    isempty(epochs.data) && @minimal_error("EpochData contains no trials.")

    n_trials = n_epochs(epochs)
    n_samples = nrow(epochs.data[1])
    ch_syms = Symbol.(epochs.layout.data.label)
    n_channels = length(ch_syms)
    time_vec = Float64.(epochs.data[1].time)

    data_tensor = Array{Float64,3}(undef, n_samples, n_channels, n_trials)
    for t = 1:n_trials
        df = epochs.data[t]
        for (c_idx, ch) in enumerate(ch_syms)
            _copy_col_to_tensor!(data_tensor, c_idx, t, df[!, ch])
        end
    end

    ride_res = ride_decompose(data_tensor, comps; times = time_vec, sample_rate = epochs.sample_rate, kwargs...)

    col_names = Vector{Symbol}(undef, n_channels + 1)
    col_names[1] = :time
    copyto!(col_names, 2, ch_syms, 1, n_channels)

    results = Dict{Symbol,ErpData}()
    for (i, comp_name) in enumerate(ride_res.components)
        cols = Vector{Any}(undef, n_channels + 1)
        cols[1] = copy(time_vec)
        for c_idx = 1:n_channels
            cols[c_idx+1] = ride_res.templates[:, c_idx, i]
        end

        results[comp_name] = ErpData(
            epochs.file,
            epochs.condition,
            epochs.condition_name * "_RIDE_$(comp_name)",
            DataFrame(cols, col_names; copycols = false),
            copy(epochs.layout),
            epochs.sample_rate,
            copy(epochs.analysis_info),
            n_trials,
        )
    end

    return results
end
