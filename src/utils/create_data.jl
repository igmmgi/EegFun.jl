"""
Functions to create some synthetic data for testing
"""

"""
    _generate_channel_labels(n_channels::Int)

Helper to generate a vector of symbols for channel labels.
"""
_generate_channel_labels(n_channels::Int) = [Symbol("Ch$i") for i = 1:n_channels]

"""
    _synthetic_signal(t; freq = 10.0, amp = 1.0, noise = 0.1)

Internal helper for generating test signals.
"""
function _synthetic_signal(t; freq = 10.0, amp = 1.0, noise = 0.1)
    return amp .* sin.(2π .* freq .* t) .+ noise .* randn(length(t))
end


"""
    create_test_continuous_data(; n::Int = 2000, fs::Int = 1000, n_channels::Int = 3, seed = nothing)

Create synthetic `ContinuousData` for testing.
"""
function create_test_continuous_data(; n = 2000, fs = 1000, n_channels = 3)
    t = collect(0:(n-1)) ./ fs
    df = DataFrame(time = t, sample = 1:n, triggers = zeros(Int, n))

    channel_labels = _generate_channel_labels(n_channels)
    for ch in channel_labels
        df[!, ch] = 0.5 .+ _synthetic_signal(t, freq = 5 + rand() * 10, amp = 1.0, noise = 0.1)
    end

    layout = create_test_layout(n_channels = n_channels)
    return EegFun.ContinuousData("test_data", copy(df, copycols = true), layout, fs, EegFun.AnalysisInfo(:none, 0.0, 0.0))
end

"""
    create_test_epoch_data(; n=2000, fs=1000, condition=1, n_epochs=10, n_channels=3, seed=nothing)

Create synthetic `EpochData` for a single condition. Returns `EpochData`.
"""
function create_test_epoch_data(; n = 2000, fs = 1000, condition = 1, n_epochs = 10, n_channels = 3)
    t = collect(0:(n-1)) ./ fs
    channel_labels = _generate_channel_labels(n_channels)
    layout = create_test_layout(n_channels = n_channels)

    dfs = DataFrame[]
    for epoch = 1:n_epochs
        df = DataFrame(
            time = t,
            sample = 1:n,
            epoch = fill(epoch, n),
            condition = fill(condition, n),
            condition_name = fill("condition_$(condition)", n),
        )

        for ch in channel_labels
            df[!, ch] = _synthetic_signal(t, freq = 10.0, amp = condition * 0.5, noise = 0.1)
        end
        push!(dfs, df)
    end

    return EegFun.EpochData("test_data", condition, "condition_$(condition)", dfs, layout, fs, EegFun.AnalysisInfo(:none, 0.0, 0.0))
end

"""
    create_test_epoch_data_vector(; conditions=1:2, n_epochs=10, n_channels=3, kwargs...)

Create a `Vector{EpochData}` across multiple conditions.
"""
function create_test_epoch_data_vector(; conditions = 1:2, n_epochs = 10, n_channels = 3, fs = 1000, n = 2000)
    return [create_test_epoch_data(condition = c, n_epochs = n_epochs, n_channels = n_channels, fs = fs, n = n) for c in conditions]
end

"""
    create_test_erp_data(; participant = 1, condition = 1, fs = 1000, n_channels = 3, seed = nothing)

Create synthetic `ErpData` for testing.
"""
function create_test_erp_data(; participant = 1, condition = 1, fs = 1000, n_channels = 3)
    time = collect(-0.5:(1/fs):2.0)
    df = DataFrame(
        time = time,
        sample = 1:length(time),
        condition = fill(condition, length(time)),
        condition_name = fill("condition_$condition", length(time)),
        participant = fill(participant, length(time)),
    )

    channel_labels = _generate_channel_labels(n_channels)
    for ch in channel_labels
        df[!, ch] = _synthetic_signal(time, freq = 0.1, amp = condition * 0.5, noise = 0.1)
    end

    layout = create_test_layout(n_channels = n_channels)
    return EegFun.ErpData("test_data", condition, "condition_$condition", df, layout, fs, EegFun.AnalysisInfo(:none, 0.0, 0.0), 10)
end

"""
    create_batch_test_erp_data(; n_conditions = 2, fs = 1000, n_channels = 3, seed = nothing)

Create a vector of ErpData objects for testing batch operations.
"""
function create_batch_test_erp_data(; n_conditions = 2, fs = 1000, n_channels = 3)
    return [create_test_erp_data(condition = c, fs = fs, n_channels = n_channels) for c = 1:n_conditions]
end


"""
    create_test_epoch_data_with_rt(; participant=1, condition=1, n_epochs=5, n_timepoints=100, n_channels=3, kwargs...)

Create synthetic `EpochData` with reaction time (RT) metadata.
"""
function create_test_epoch_data_with_rt(;
    participant = 1,
    condition = 1,
    n_epochs = 5,
    n_timepoints = 100,
    n_channels = 3,
    fs = 1000,
    epoch_start = -0.2,
    epoch_end = 0.8,
    rt_range = (0.2, 0.6),
)
    time = collect(range(epoch_start, epoch_end, length = n_timepoints))
    channel_labels = _generate_channel_labels(n_channels)

    dfs = DataFrame[]
    for epoch = 1:n_epochs
        df = DataFrame(
            time = time,
            condition = fill(condition, n_timepoints),
            condition_name = fill("condition_$condition", n_timepoints),
            participant = fill(participant, n_timepoints),
            epoch = fill(epoch, n_timepoints),
            rt = fill(rand() * (rt_range[2] - rt_range[1]) + rt_range[1], n_timepoints),
        )

        for ch in channel_labels
            df[!, ch] = _synthetic_signal(time, freq = 10.0, amp = condition * 0.5, noise = 0.1)
        end
        push!(dfs, df)
    end

    layout = create_test_layout(n_channels = n_channels)
    return EegFun.EpochData("test_data", condition, "condition_$condition", dfs, layout, fs, EegFun.AnalysisInfo(:none, 0.0, 0.0))
end


"""
    create_test_continuous_data_with_triggers(; n = 1000, fs = 1000, seed = nothing)

Create synthetic `ContinuousData` with trigger sequences.
"""
function create_test_continuous_data_with_triggers(; n = 1000, fs = 1000)
    t = collect(0:(n-1)) ./ fs
    x = _synthetic_signal(t, freq = 10.0, amp = 1.0, noise = 0.1)

    triggers = zeros(Int, n)
    triggers[300:302] .= [1, 7, 3]
    triggers[400:402] .= [1, 2, 3]
    triggers[750:752] .= [1, 2, 3]
    triggers[850] = 8
    triggers[900] = 9
    triggers[950:952] .= [1, 2, 3]

    df = DataFrame(time = t, triggers = triggers, A = x)
    layout = EegFun.Layout(DataFrame(label = [:A], inc = [0.0], azi = [0.0]), nothing, nothing)
    return EegFun.ContinuousData("test_data", copy(df, copycols = true), layout, fs, EegFun.AnalysisInfo(:none, 0.0, 0.0))
end


function create_test_continuous_data_with_artifacts(; n::Int = 1000, fs::Int = 1000)

    t = collect(0:(n-1)) ./ fs

    # Create clean signal with some artifacts
    clean_signal = sin.(2π .* 10 .* t) .* 20
    # Add some extreme values (artifacts)
    artifact_signal = copy(clean_signal)
    artifact_signal[100:110] .= 200.0  # Large positive artifact
    artifact_signal[500:505] .= -200.0  # Large negative artifact
    artifact_signal[800:802] .= 100.0  # Smaller positive artifact

    df = DataFrame(:time => t, :triggers => zeros(Int, n), :Ch1 => clean_signal, :Ch2 => artifact_signal)
    layout = EegFun.Layout(DataFrame(label = [:Ch1, :Ch2], inc = [0.0, 0.0], azi = [0.0, 0.0]), nothing, nothing)

    dat = EegFun.ContinuousData("test_data", df, layout, fs, EegFun.AnalysisInfo())

    return dat

end


"""
    create_test_continuous_data_empty_triggers(; n_samples = 1000, fs = 1000, seed = nothing)

Create `ContinuousData` with zero triggers.
"""
function create_test_continuous_data_empty_triggers(; n_samples = 1000, fs = 1000)
    t = collect(0:(n_samples-1)) ./ fs
    df = DataFrame(time = t, triggers = zeros(Int16, n_samples), channel1 = _synthetic_signal(t), channel2 = _synthetic_signal(t))

    layout = EegFun.Layout(DataFrame(label = [:channel1, :channel2], inc = [0.0, 90.0], azi = [0.0, 0.0]), nothing, nothing)
    return EegFun.ContinuousData("test_data", df, layout, fs, EegFun.AnalysisInfo(:none, 0.0, 0.0))
end

"""
    create_test_lrp_data(; participant = 1, condition = 1, n_timepoints = 100, signal_scale = 1.0, fs = 250, seed = nothing)

Create `ErpData` with lateralized channel pairs.
"""
function create_test_lrp_data(; participant = 1, condition = 1, n_timepoints = 100, signal_scale = 1.0, fs = 250)
    time = collect(range(-0.2, 0.8, length = n_timepoints))

    df = DataFrame(
        time = time,
        sample = 1:n_timepoints,
        condition = fill(condition, n_timepoints),
        condition_name = fill("condition_$condition", n_timepoints),
        participant = fill(participant, n_timepoints),
        file = fill("test_file", n_timepoints),
    )

    # Left hemisphere
    df.C3 = _synthetic_signal(time, freq = 0.1, amp = signal_scale, noise = 0.01)
    df.C1 = _synthetic_signal(time, freq = 0.1, amp = signal_scale, noise = 0.01)
    df.Fp1 = _synthetic_signal(time, freq = 0.2, amp = signal_scale, noise = 0.01)

    # Right hemisphere (shifted)
    df.C4 = _synthetic_signal(time, freq = 0.1, amp = signal_scale, noise = 0.01) .+ 0.1
    df.C2 = _synthetic_signal(time, freq = 0.1, amp = signal_scale, noise = 0.01) .+ 0.1
    df.Fp2 = _synthetic_signal(time, freq = 0.2, amp = signal_scale, noise = 0.01) .+ 0.1

    df.Fz = _synthetic_signal(time, freq = 0.15, amp = signal_scale, noise = 0.01)

    channel_labels = [:C3, :C4, :C1, :C2, :Fp1, :Fp2, :Fz]
    layout = EegFun.Layout(
        DataFrame(label = channel_labels, inc = zeros(length(channel_labels)), azi = zeros(length(channel_labels))),
        nothing,
        nothing,
    )

    return EegFun.ErpData("test_data", condition, "condition_$condition", df, layout, fs, EegFun.AnalysisInfo(:none, 0.0, 0.0), 10)
end



function create_test_layout(; n_channels::Int = 4, layout_type::Symbol = :grid)
    if layout_type == :grid # Create grid positions
        rows = ceil(Int, sqrt(n_channels))
        cols = ceil(Int, n_channels / rows)
        positions = []
        for i = 1:n_channels
            row = (i - 1) ÷ cols + 1
            col = (i - 1) % cols + 1
            x = (col - 1) / max(1, cols - 1)
            y = (row - 1) / max(1, rows - 1)
            push!(positions, (x, y))
        end
        # Create DataFrame with proper structure
        layout_data = DataFrame(
            :label => _generate_channel_labels(n_channels),
            :inc => zeros(n_channels),
            :azi => zeros(n_channels),
            :x2 => [pos[1] for pos in positions],
            :y2 => [pos[2] for pos in positions],
        )
    elseif layout_type == :topo # Create a topographic layout for 6 channels
        layout_data = DataFrame(
            label = [:Fp1, :Fp2, :Cz, :Pz, :O1, :O2],
            inc = zeros(6),
            azi = zeros(6),
            x2 = [0.0, 0.0, 0.0, 0.0, -0.5, 0.5],
            y2 = [1.0, 1.0, 0.0, -1.0, -1.0, -1.0],
            x3 = [0.0, 0.0, 0.0, 0.0, -0.5, 0.5],
            y3 = [1.0, 1.0, 0.0, -1.0, -1.0, -1.0],
            z3 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        )
    else
        error("Unknown layout_type: $layout_type")
    end

    return EegFun.Layout(layout_data, nothing, nothing)
end

function create_test_summary_data()
    return DataFrame(
        channel = [:Ch1, :Ch2, :Ch3, :Ch4, :Ch5],
        min = [-2.1, -1.8, -0.5, -0.3, -0.2],
        max = [2.3, 1.9, 0.8, 0.4, 0.3],
        std = [1.2, 1.1, 0.4, 0.2, 0.15],
        var = [1.44, 1.21, 0.16, 0.04, 0.0225],
        range = [4.4, 3.7, 1.3, 0.7, 0.5],
        zvar = [1.2, 0.8, -0.5, -1.1, -1.4],
    )
end

function create_test_summary_data_with_epochs()
    return DataFrame(
        channel = repeat([:Ch1, :Ch2, :Ch3], 3),
        epoch = repeat([1, 2, 3], 3),
        min = [-2.1, -1.8, -0.5, -2.0, -1.7, -0.4, -2.2, -1.9, -0.6],
        max = [2.3, 1.9, 0.8, 2.2, 1.8, 0.7, 2.4, 2.0, 0.9],
        std = [1.2, 1.1, 0.4, 1.1, 1.0, 0.3, 1.3, 1.2, 0.5],
        var = [1.44, 1.21, 0.16, 1.21, 1.0, 0.09, 1.69, 1.44, 0.25],
        range = [4.4, 3.7, 1.3, 4.2, 3.5, 1.1, 4.6, 3.9, 1.5],
        zvar = [1.2, 0.8, -0.5, 0.9, 0.6, -0.8, 1.5, 1.1, -0.2],
    )
end

"""
    create_test_epoch_data_with_artifacts(; participant=1, condition=1, n_epochs=10, n_timepoints=100, n_channels=3, n_bad_epochs=1, seed=nothing)

Create `EpochData` with large amplitude noise in some epochs for testing artifact detection.
"""
function create_test_epoch_data_with_artifacts(;
    participant = 1,
    condition = 1,
    n_epochs = 10,
    n_timepoints = 100,
    n_channels = 3,
    fs = 1000,
    n_bad_epochs = 1,
)
    dfs = DataFrame[]
    bad_indices = collect(1:n_bad_epochs)

    time = collect(0:(n_timepoints-1)) ./ fs

    for epoch = 1:n_epochs
        df = DataFrame(
            time = time,
            epoch = fill(epoch, n_timepoints),
            participant = fill(participant, n_timepoints),
            condition = fill(condition, n_timepoints),
            condition_name = fill("condition_$condition", n_timepoints),
        )

        is_bad = epoch <= n_bad_epochs

        for ch_idx = 1:n_channels
            ch_name = Symbol("Ch$ch_idx")
            signal = _synthetic_signal(time, freq = 10.0)

            if is_bad
                signal .+= 50.0 .* randn(n_timepoints)
                signal[1:min(20, n_timepoints)] .+= 200.0
                if n_timepoints > 70
                    signal[50:70] .+= 300.0
                end
            end
            df[!, ch_name] = signal
        end
        push!(dfs, df)
    end

    layout = create_test_layout(n_channels = n_channels)
    epoch_data = EegFun.EpochData("test_data", condition, "condition_$condition", dfs, layout, fs, EegFun.AnalysisInfo(:none, 0.0, 0.0))

    return epoch_data, bad_indices
end

"""
    generate_signal(n_trials, time_window, sample_rate, signal_freqs, signal_amps, signal_times, noise)

Generate synthetic signals with known frequency components for testing TF analysis.

# Arguments
- `n_trials`: Number of trials to generate
- `time_window`: [start_time, end_time] in seconds
- `sample_rate`: Sampling rate in Hz
- `signal_freqs`: Vector of frequencies for each component
- `signal_amps`: Vector of amplitudes for each component
- `signal_times`: Vector of [start_time, end_time] for each component
- `noise`: Amplitude of random noise

# Returns
- `times`: Time vector
- `signal`: Matrix (samples × trials)
"""
function generate_signal(
    n_trials::Int,
    time_window::Vector{Float64},
    sample_rate::Real,
    signal_freqs::Vector{<:Real},
    signal_amps::Vector{<:Real},
    signal_times::Vector{Vector{Float64}},
    noise::Real,
)
    # Pre-compute time vector
    dt = 1 / sample_rate
    x_time = range(time_window[1], time_window[2] - dt, step = dt)
    n_samples = length(x_time)

    # Pre-allocate output matrix with noise
    signal_trials = zeros(n_samples, n_trials)
    if noise > 0
        @views for trial = 1:n_trials
            signal_trials[:, trial] .= randn(n_samples) .* noise
        end
    end

    # Pre-allocate temporary arrays
    temp_signal = zeros(n_samples)

    # Find time 0 in the time array (same for all trials and frequencies)
    # signal_times are relative to time 0 (event/stimulus time), not time_window[1]
    time_zero_idx = argmin(abs.(x_time))
    time_zero = x_time[time_zero_idx]

    # Add frequency components to all trials
    @inbounds for trial = 1:n_trials
        # Reset temporary signal array
        fill!(temp_signal, 0.0)

        # Add each frequency component
        for (freq, amp, times) in zip(signal_freqs, signal_amps, signal_times)
            # Convert signal_times from relative to time 0 to absolute times
            abs_times = [time_zero + times[1], time_zero + times[2]]

            # Find time window indices (now using absolute times)
            time_mask = (x_time .>= abs_times[1]) .& (x_time .< abs_times[2])
            window_times = @view x_time[time_mask]

            # Generate signal for this component
            @. temp_signal[time_mask] += amp * sin(2π * freq * window_times)
        end

        # Add to output
        signal_trials[:, trial] .+= temp_signal
    end

    return x_time, signal_trials
end

"""
    signal_to_data(times, signal, channel_name::Symbol, sample_rate::Real;
                     file::String="synthetic", condition::Int=1, condition_name::String="test") -> Union{EpochData, ErpData}

Internal function: Convert signal data from `generate_signal` to `EpochData` or `ErpData` format.
Returns `ErpData` when there is only one epoch, otherwise returns `EpochData`.

"""
function signal_to_data(
    times,
    signal,
    channel_name::Symbol,
    sample_rate::Real;
    file::String = "synthetic",
    condition::Int = 1,
    condition_name::String = "test",
)
    # Ensure signal is 2D (samples × trials)
    if ndims(signal) == 1
        signal = reshape(signal, :, 1)
    end

    n_samples, n_trials = size(signal)

    # Convert times to vector if needed
    times_vec = collect(Float64, times)
    @assert length(times_vec) == n_samples "Time vector length ($(length(times_vec))) must match signal samples ($n_samples)"

    # Create minimal layout for the channel
    layout_df = DataFrame(:label => [channel_name])
    layout = Layout(layout_df, nothing, nothing)

    # Create AnalysisInfo
    analysis_info = AnalysisInfo()

    # If only one epoch, return ErpData
    if n_trials == 1
        erp_df = DataFrame(:time => times_vec, channel_name => signal[:, 1])
        return ErpData(file, condition, condition_name, erp_df, layout, Int(sample_rate), analysis_info, 1)
    else
        trial_dfs = Vector{DataFrame}(undef, n_trials)
        for trial_idx = 1:n_trials
            trial_dfs[trial_idx] = DataFrame(:time => times_vec, channel_name => signal[:, trial_idx], :epoch => trial_idx)
        end
        return EpochData(file, condition, condition_name, trial_dfs, layout, Int(sample_rate), analysis_info)
    end
end

