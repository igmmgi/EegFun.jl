"""
    tf_baseline!(tf_data::TimeFreqData, baseline_interval; method=:db)

Apply baseline correction to TimeFreqData in-place.

# Arguments
- `tf_data::TimeFreqData`: Time-frequency data to baseline correct
- `baseline_interval::Tuple{Real,Real}`: Time window for baseline (start, stop) in seconds

# Keyword Arguments
- `method::Symbol=:db`: Baseline method 
  - `:absolute`: Absolute change (power - baseline_mean) - simple subtraction, no normalization
  - `:relative`: Relative power (power / baseline_mean) - ratio, 
  - `:relchange`: Relative change ((power - baseline_mean) / baseline_mean) - fractional change 
  - `:normchange`: Normalized change ((power - baseline) / (power + baseline)) - symmetric normalization
  - `:db`: Decibel change (10 * log10(power/baseline_mean)) 
  - `:vssum`: Alias for `:normchange`
  - `:zscore`: Z-score normalization ((power - baseline_mean) / baseline_std) 
  - `:percent`: Percent change (100 * (power - baseline) / baseline) - convenience alias for relchange × 100

# Example
```julia
tf_baseline!(tf_data, (-0.3, 0.1); method=:db)
```
"""
function tf_baseline!(tf_data::TimeFreqData, baseline_interval::Tuple{Real,Real}; method::Symbol = :db)
    # Check if baseline has already been applied
    if tf_data.baseline |> !isnothing
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
    # Find the nearest time points to the baseline interval boundaries
    baseline_start_idx = argmin(abs.(times .- baseline_interval[1]))
    baseline_end_idx = argmin(abs.(times .- baseline_interval[2]))

    # Create mask for baseline time points
    base_mask = falses(n_times)
    base_mask[baseline_start_idx:baseline_end_idx] .= true
    if !any(base_mask)
        error("Baseline window $(baseline_interval) does not overlap with data times")
    end

    # Get channel columns
    ch_labels = channel_labels(tf_data)
    n_channels = length(ch_labels)

    @info "Applying $(method) baseline correction ($(baseline_interval[1])s to $(baseline_interval[2])s)"

    # Build row index once
    row_index = _tf_build_row_index(tf_data.data_power)

    # Extract ALL channels into a 3D array [n_channels × n_freqs × n_times] in one pass
    power_3d = Array{Float64,3}(undef, n_channels, n_freqs, n_times)
    for (fi, f) in enumerate(freqs_unique)
        rf = round(f, digits = 6)
        for (ti, t) in enumerate(times)
            row = get(row_index, (rf, round(t, digits = 6)), nothing)
            if row |> !isnothing
                for (ci, ch) in enumerate(ch_labels)
                    power_3d[ci, fi, ti] = tf_data.data_power[row, ch]
                end
            else
                power_3d[:, fi, ti] .= NaN
            end
        end
    end

    # Compute baseline statistics per channel × frequency
    # baseline_mean: [n_channels × n_freqs]
    baseline_mean = zeros(n_channels, n_freqs)
    baseline_std_mat = nothing  # only computed for :zscore

    for ci = 1:n_channels
        for fi = 1:n_freqs
            vals = @view power_3d[ci, fi, base_mask]
            valid = vals[.!isnan.(vals)]
            baseline_mean[ci, fi] = isempty(valid) ? NaN : mean(valid)
        end
    end

    if method == :zscore
        baseline_std_mat = zeros(n_channels, n_freqs)
        for ci = 1:n_channels
            for fi = 1:n_freqs
                vals = @view power_3d[ci, fi, base_mask]
                valid = vals[.!isnan.(vals)]
                baseline_std_mat[ci, fi] = (isempty(valid) || length(valid) < 2) ? NaN : std(valid, corrected = false)
            end
        end
    end

    # Apply baseline correction to entire 3D array at once (no per-channel loop)
    # Reshape baseline_mean to [n_channels × n_freqs × 1] for broadcasting over time
    bm = reshape(baseline_mean, n_channels, n_freqs, 1)

    if method == :absolute
        power_3d .-= bm
    elseif method == :relative
        safe_bm = max.(bm, 1e-30)
        power_3d ./= safe_bm
    elseif method == :relchange
        safe_bm = max.(bm, 1e-30)
        power_3d .= (power_3d .- bm) ./ safe_bm
    elseif method == :normchange || method == :vssum
        power_3d .= (power_3d .- bm) ./ (power_3d .+ bm)
    elseif method == :db
        safe_bm = max.(bm, 1e-30)
        power_3d .= 10 .* log10.(max.(power_3d, 1e-30) ./ safe_bm)
    elseif method == :zscore
        safe_std = reshape(max.(baseline_std_mat, 1e-30), n_channels, n_freqs, 1)
        power_3d .= (power_3d .- bm) ./ safe_std
    elseif method == :percent
        safe_bm = max.(bm, 1e-30)
        power_3d .= 100 .* (power_3d .- bm) ./ safe_bm
    else
        error(
            "Unknown baseline method: $method. Use :absolute, :relative, :relchange, :normchange, :db, :vssum (alias for :normchange), :zscore, or :percent",
        )
    end

    # Write all channels back to DataFrame in one pass
    for (fi, f) in enumerate(freqs_unique)
        rf = round(f, digits = 6)
        for (ti, t) in enumerate(times)
            row = get(row_index, (rf, round(t, digits = 6)), nothing)
            if row |> !isnothing
                for (ci, ch) in enumerate(ch_labels)
                    tf_data.data_power[row, ch] = power_3d[ci, fi, ti]
                end
            end
        end
    end

    # Store baseline information
    tf_data.baseline = BaselineInfo(method = method, window = (Float64(baseline_interval[1]), Float64(baseline_interval[2])))

    return nothing
end

# Non-mutating version
function tf_baseline(tf_data::TimeFreqData, baseline_interval::Tuple{Real,Real}; method::Symbol = :db)
    tf_copy = copy(tf_data)  # Use custom copy method instead of deepcopy
    tf_baseline!(tf_copy, baseline_interval; method = method)
    return tf_copy
end

# Vector version
function tf_baseline!(tf_data::Vector{TimeFreqData}, baseline_interval::Tuple{Real,Real}; method::Symbol = :db)
    tf_baseline!.(tf_data, Ref(baseline_interval); method = method)
    return nothing
end

function tf_baseline(tf_data::Vector{TimeFreqData}, baseline_interval::Tuple{Real,Real}; method::Symbol = :db)
    return [tf_baseline(tf, baseline_interval; method = method) for tf in tf_data]
end

"""
    tf_baseline(file_pattern::String, baseline_interval;
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
    baseline_interval::Tuple{Real,Real};
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

        if (error_msg = _validate_input_dir(input_dir)) |> !isnothing
            @minimal_error(error_msg)
        end

        output_dir = something(output_dir, joinpath(input_dir, "tf_baseline_$(file_pattern)"))
        mkpath(output_dir)

        files = _find_batch_files(file_pattern, input_dir, participant_selection)
        if isempty(files)
            @minimal_warning "No JLD2 files found matching pattern '$file_pattern'"
            return nothing
        end

        @info "Found $(length(files)) files, baseline: $baseline_interval, method: $method"

        process_fn = (input_path, output_path) -> begin
            filename = basename(input_path)
            data = read_data(input_path)
            if isnothing(data) || !(data isa Vector{TimeFreqData})
                return BatchResult(false, filename, "Invalid data type")
            end
            data = _condition_select(data, condition_selection)
            tf_baseline!(data, baseline_interval; method = method)
            jldsave(output_path; data = data)
            return BatchResult(true, filename, "Baseline corrected")
        end

        results = _run_batch_operation(process_fn, files, input_dir, output_dir; operation_name = "TF baseline")
        _log_batch_summary(results, output_dir)

    finally
        _cleanup_logging(log_file, output_dir)
    end
end
