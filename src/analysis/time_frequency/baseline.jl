"""
    tf_baseline!(tf_data::TimeFreqData, baseline_window; method=:db)

Apply baseline correction to TimeFreqData in-place.

# Arguments
- `tf_data::TimeFreqData`: Time-frequency data to baseline correct
- `baseline_window::Tuple{Real,Real}`: Time window for baseline (start, stop) in seconds

# Keyword Arguments
- `method::Symbol=:db`: Baseline method 
  - `:absolute`: Absolute change (power - baseline_mean) - simple subtraction, no normalization
  - `:relative`: Relative power (power / baseline_mean) - ratio, 
  - `:relchange`: Relative change ((power - baseline_mean) / baseline_mean) - fractional change 
  - `:normchange`: Normalized change ((power - baseline) / (power + baseline)) - symmetric normalization
  - `:db`: Decibel change (10 * log10(power/baseline_mean)) 
  - `:vssum`: Variance-stabilized sum ((power - baseline) / (power + baseline)) 
  - `:zscore`: Z-score normalization ((power - baseline_mean) / baseline_std) 
  - `:percent`: Percent change (100 * (power - baseline) / baseline) - convenience alias for relchange × 100

# Example
```julia
tf_baseline!(tf_data, (-0.3, 0.1); method=:db)
```
"""
function tf_baseline!(tf_data::TimeFreqData, baseline_window::Tuple{Real,Real}; method::Symbol = :db)
    # Check if baseline has already been applied
    if tf_data.baseline !== nothing
        error(
            "Baseline correction has already been applied to this data (method: $(tf_data.baseline.method), window: $(tf_data.baseline.window)). " *
            "Baseline corrections are non-linear and cannot be chained. Use the original data to apply a different baseline.",
        )
    end

    # The DataFrame is structured as: all frequencies for time 1, then all frequencies for time 2, etc.
    # unique() returns frequencies in the order they first appear, which matches freqs_out from tf_analysis
    # So we can use direct indexing: row_idx = (ti - 1) * n_freqs + fi
    times = unique(tf_data.data_power.time)
    freqs_unique = unique(tf_data.data_power.freq)  # Already in correct order (matches freqs_out)
    n_freqs = length(freqs_unique)
    n_times = length(times)

    # Find baseline time indices (MATLAB: dsearchn finds nearest time points)
    # Find the nearest time points to the baseline window boundaries
    baseline_start_idx = argmin(abs.(times .- baseline_window[1]))
    baseline_end_idx = argmin(abs.(times .- baseline_window[2]))

    # Create mask for baseline time points
    base_mask = falses(n_times)
    base_mask[baseline_start_idx:baseline_end_idx] .= true
    if !any(base_mask)
        error("Baseline window $(baseline_window) does not overlap with data times")
    end

    # Get channel columns
    ch_labels = channel_labels(tf_data)

    # Process each channel 
    for ch in ch_labels
        # Reshape to freq × time matrix for baseline calculation
        # DataFrame structure: [freq1_t1, freq2_t1, ..., freqN_t1, freq1_t2, ...]
        # freqs_unique is already in the correct order (matches DataFrame row order)
        # So row_idx = (ti - 1) * n_freqs + fi
        power_mat = zeros(n_freqs, n_times)

        for ti = 1:n_times
            for fi = 1:n_freqs
                row_idx = (ti - 1) * n_freqs + fi
                power_mat[fi, ti] = tf_data.data_power[row_idx, ch]
            end
        end

        # Compute baseline statistics per frequency
        # base_mask is a boolean vector for time points
        # power_mat is (n_freqs, n_times)
        # We want mean (and std for zscore) across time points in baseline window for each frequency
        # Skip NaN values (from edge filtering) when computing baseline statistics 
        baseline_power = zeros(n_freqs)
        for fi = 1:n_freqs
            baseline_values = power_mat[fi, base_mask]
            baseline_values_no_nan = baseline_values[.!isnan.(baseline_values)]
            if isempty(baseline_values_no_nan)
                baseline_power[fi] = NaN  # All baseline values are NaN
            else
                baseline_power[fi] = mean(baseline_values_no_nan)
            end
        end

        # Apply baseline correction 
        if method == :absolute
            # 'absolute': data - meanVals (simple subtraction)
            @info "Applying absolute baseline correction"
            power_mat .= power_mat .- reshape(baseline_power, n_freqs, 1)
        elseif method == :relative
            # 'relative': power / baseline_mean
            @info "Applying relative baseline correction"
            min_baseline = max.(baseline_power, 1e-30)
            power_mat .= power_mat ./ reshape(min_baseline, n_freqs, 1)
        elseif method == :relchange
            # 'relchange': (power - baseline_mean) / baseline_mean
            @info "Applying relchange baseline correction"
            min_baseline = max.(baseline_power, 1e-30)
            power_mat .= (power_mat .- reshape(baseline_power, n_freqs, 1)) ./ reshape(min_baseline, n_freqs, 1)
        elseif method == :normchange
            # 'normchange': (power - baseline) / (power + baseline)
            @info "Applying normchange baseline correction"
            power_mat .= (power_mat .- reshape(baseline_power, n_freqs, 1)) ./ (power_mat .+ reshape(baseline_power, n_freqs, 1))
        elseif method == :db
            # 'db': 10 * log10(power / baseline_mean)
            @info "Applying db baseline correction"
            min_power = max.(baseline_power, 1e-30)
            power_mat .= 10 .* log10.(max.(power_mat, 1e-30) ./ reshape(min_power, n_freqs, 1))
        elseif method == :vssum
            # 'vssum': (power - baseline) / (power + baseline) - same as normchange
            @info "Applying vssum baseline correction"
            power_mat .= (power_mat .- reshape(baseline_power, n_freqs, 1)) ./ (power_mat .+ reshape(baseline_power, n_freqs, 1))
        elseif method == :zscore
            # 'zscore': (power - baseline_mean) / baseline_std
            # uses population std (divide by N, not N-1): nanstd(data(:,:,baselineTimes),1, 3)
            @info "Applying z-score baseline correction"
            # Compute standard deviation for each frequency across baseline time points
            # Skip NaN values (from edge filtering) when computing std 
            baseline_std = zeros(n_freqs)
            for fi = 1:n_freqs
                baseline_values = power_mat[fi, base_mask]
                baseline_values_no_nan = baseline_values[.!isnan.(baseline_values)]
                if isempty(baseline_values_no_nan) || length(baseline_values_no_nan) < 2
                    baseline_std[fi] = NaN  # All baseline values are NaN or insufficient data
                else
                    # Use population std (divide by N, not N-1) 
                    baseline_std[fi] = std(baseline_values_no_nan, corrected = false)
                end
            end
            # Avoid division by zero - use minimum threshold for std
            min_std = max.(baseline_std, 1e-30)
            power_mat .= (power_mat .- reshape(baseline_power, n_freqs, 1)) ./ reshape(min_std, n_freqs, 1)
        elseif method == :percent
            # Convenience alias: percent change = relchange × 100
            @info "Applying percent baseline correction (convenience alias for relchange × 100)"
            min_baseline = max.(baseline_power, 1e-30)
            power_mat .= 100 .* (power_mat .- reshape(baseline_power, n_freqs, 1)) ./ reshape(min_baseline, n_freqs, 1)
        else
            error("Unknown baseline method: $method. Use :absolute, :relative, :relchange, :normchange, :db, :vssum, :zscore, or :percent")
        end

        # Write back to DataFrame (using same indexing)
        for ti = 1:n_times
            for fi = 1:n_freqs
                row_idx = (ti - 1) * n_freqs + fi
                tf_data.data_power[row_idx, ch] = power_mat[fi, ti]
            end
        end
    end

    # Store baseline information
    tf_data.baseline = BaselineInfo(method = method, window = (Float64(baseline_window[1]), Float64(baseline_window[2])))

    return nothing
end

# Non-mutating version
function tf_baseline(tf_data::TimeFreqData, baseline_window::Tuple{Real,Real}; method::Symbol = :db)
    tf_copy = copy(tf_data)  # Use custom copy method instead of deepcopy
    tf_baseline!(tf_copy, baseline_window; method = method)
    return tf_copy
end

# Vector version
function tf_baseline!(tf_data::Vector{TimeFreqData}, baseline_window::Tuple{Real,Real}; method::Symbol = :db)
    tf_baseline!.(tf_data, Ref(baseline_window); method = method)
    return nothing
end

function tf_baseline(tf_data::Vector{TimeFreqData}, baseline_window::Tuple{Real,Real}; method::Symbol = :db)
    return [tf_baseline(tf, baseline_window; method = method) for tf in tf_data]
end

"""
    tf_baseline(file_pattern::String, baseline_window;
                method=:db, input_dir=pwd(), output_dir=nothing,
                participant_selection=participants(), condition_selection=conditions())

Apply baseline correction to TimeFreqData files.

# Example
```julia
tf_baseline("tf_epochs_wavelet", (-0.3, 0.0); method=:db)
```
"""
function tf_baseline(
    file_pattern::String,
    baseline_window::Tuple{Real,Real};
    method::Symbol = :db,
    input_dir::String = pwd(),
    output_dir::Union{String,Nothing} = nothing,
    participant_selection::Function = participants(),
    condition_selection::Function = conditions(),
)
    log_file = "tf_baseline.log"
    setup_global_logging(log_file)

    try
        @info "Batch TF baseline started at $(now())"
        @log_call "tf_baseline"

        if (error_msg = _validate_input_dir(input_dir)) !== nothing
            @minimal_error_throw(error_msg)
        end

        output_dir = something(output_dir, joinpath(input_dir, "tf_baseline_$(file_pattern)"))
        mkpath(output_dir)

        files = _find_batch_files(file_pattern, input_dir, participant_selection)
        if isempty(files)
            @minimal_warning "No JLD2 files found matching pattern '$file_pattern'"
            return nothing
        end

        @info "Found $(length(files)) files, baseline: $baseline_window, method: $method"

        process_fn = (input_path, output_path) -> begin
            filename = basename(input_path)
            data = load_data(input_path)
            if isnothing(data) || !(data isa Vector{TimeFreqData})
                return BatchResult(false, filename, "Invalid data type")
            end
            data = _condition_select(data, condition_selection)
            tf_baseline!(data, baseline_window; method = method)
            jldsave(output_path; data = data)
            return BatchResult(true, filename, "Baseline corrected")
        end

        results = _run_batch_operation(process_fn, files, input_dir, output_dir; operation_name = "TF baseline")
        _log_batch_summary(results, output_dir)

    finally
        _cleanup_logging(log_file, output_dir)
    end
end
