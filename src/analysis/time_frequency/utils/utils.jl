function _filter_edges!(
    eegpower::AbstractArray,
    eegconv::AbstractArray,
    num_frex::Int,
    time_indices::AbstractVector{Int},
    window_lengths_per_freq::Union{Vector{Int},Vector{Float64}},
    n_samples_per_epoch::Int,
)
    for fi = 1:num_frex

        half_nsamplefreqoi = window_lengths_per_freq[fi] / 2.0
        min_valid_threshold = half_nsamplefreqoi
        max_valid_threshold = n_samples_per_epoch - half_nsamplefreqoi

        for ti in eachindex(time_indices)
            sample_idx = time_indices[ti]
            if !(sample_idx >= min_valid_threshold && sample_idx < max_valid_threshold)
                if ndims(eegpower) == 3  # return_trials = true
                    @views eegpower[fi, ti, :] .= NaN
                    @views eegconv[fi, ti, :] .= NaN * im
                else  # return_trials = false
                    eegpower[fi, ti] = NaN
                    eegconv[fi, ti] = NaN * im
                end
            end
        end
    end
end


