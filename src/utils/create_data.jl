_generate_channel_labels(n_channels::Int) = [Symbol("Ch$i") for i = 1:n_channels]

function _with_seed(f, seed::Union{Integer,Nothing})
    isnothing(seed) && return f()
    saved = copy(Random.default_rng())
    Random.seed!(seed)
    try
        return f()
    finally
        copy!(Random.default_rng(), saved)
    end
end

_synthetic_signal(t; freq = 10.0, amp = 1.0, noise = 0.1) = amp .* sin.(2π .* freq .* t) .+ noise .* randn(length(t))

function create_test_continuous_data(; n = 2000, fs = 1000, n_channels = 3, seed::Union{Int,Nothing} = nothing)
    return _with_seed(seed) do
        t = collect(0:(n-1)) ./ fs
        df = DataFrame(time = t, sample = 1:n, trigger = zeros(Int, n))
        for ch in _generate_channel_labels(n_channels)
            df[!, ch] = 0.5 .+ _synthetic_signal(t, freq = 5 + rand() * 10, amp = 1.0, noise = 0.1)
        end
        return ContinuousData("test_data", df, create_test_layout(n_channels = n_channels), fs, AnalysisInfo(sample_rate = fs))
    end
end

function create_test_epoch_data(; n = 2000, fs = 1000, condition = 1, n_epochs = 10, n_channels = 3, seed::Union{Int,Nothing} = nothing)
    return _with_seed(seed) do
        t = collect(0:(n-1)) ./ fs
        labels = _generate_channel_labels(n_channels)
        layout = create_test_layout(n_channels = n_channels)
        dfs = Vector{DataFrame}(undef, n_epochs)
        for epoch = 1:n_epochs
            df = DataFrame(time = copy(t), sample = 1:n, epoch = fill(epoch, n))
            for ch in labels
                df[!, ch] = _synthetic_signal(t, freq = 10.0, amp = condition * 0.5, noise = 0.1)
            end
            dfs[epoch] = df
        end
        return EpochData("test_data", condition, "condition_$(condition)", dfs, layout, fs, AnalysisInfo(sample_rate = fs))
    end
end

function create_test_epoch_data_vector(;
    conditions = 1:2,
    n_epochs = 10,
    n_channels = 3,
    fs = 1000,
    n = 2000,
    seed::Union{Int,Nothing} = nothing,
)
    return _with_seed(seed) do
        [create_test_epoch_data(condition = c, n_epochs = n_epochs, n_channels = n_channels, fs = fs, n = n) for c in conditions]
    end
end

function create_test_erp_data(; participant = 1, condition = 1, fs = 1000, n_channels = 3, seed::Union{Int,Nothing} = nothing)
    return _with_seed(seed) do
        time = collect(-0.5:(1/fs):2.0)
        df = DataFrame(time = time, sample = 1:length(time), participant = fill(participant, length(time)))
        for ch in _generate_channel_labels(n_channels)
            df[!, ch] = _synthetic_signal(time, freq = 0.1, amp = condition * 0.5, noise = 0.1)
        end
        return ErpData(
            "test_data",
            condition,
            "condition_$condition",
            df,
            create_test_layout(n_channels = n_channels),
            fs,
            AnalysisInfo(sample_rate = fs),
            10,
        )
    end
end

function create_test_batch_erp_data(; n_conditions = 2, fs = 1000, n_channels = 3, seed::Union{Int,Nothing} = nothing)
    return _with_seed(seed) do
        [create_test_erp_data(condition = c, fs = fs, n_channels = n_channels) for c = 1:n_conditions]
    end
end

function create_test_continuous_data_with_triggers(; n = 1000, fs = 1000, seed::Union{Int,Nothing} = nothing)
    return _with_seed(seed) do
        t = collect(0:(n-1)) ./ fs
        x = _synthetic_signal(t, freq = 10.0, amp = 1.0, noise = 0.1)
        trigger = zeros(Int, n)
        trigger[300:302] .= (1, 7, 3)
        trigger[400:402] .= (1, 2, 3)
        trigger[750:752] .= (1, 2, 3)
        trigger[850] = 8
        trigger[900] = 9
        trigger[950:952] .= (1, 2, 3)
        df = DataFrame(time = t, trigger = trigger, A = x)
        layout = Layout(DataFrame(label = [:A], inc = [0.0], azi = [0.0]), nothing, nothing, nothing)
        return ContinuousData("test_data", df, layout, fs, AnalysisInfo(sample_rate = fs))
    end
end

function create_test_continuous_data_with_artifacts(; n::Int = 1000, fs::Int = 1000)
    t = collect(0:(n-1)) ./ fs
    clean_signal = sin.(2π .* 10 .* t) .* 20
    artifact_signal = copy(clean_signal)
    artifact_signal[100:110] .= 200.0
    artifact_signal[500:505] .= -200.0
    artifact_signal[800:802] .= 100.0
    df = DataFrame(:time => t, :trigger => zeros(Int, n), :Ch1 => clean_signal, :Ch2 => artifact_signal)
    layout = Layout(DataFrame(label = [:Ch1, :Ch2], inc = [0.0, 0.0], azi = [0.0, 0.0]), nothing, nothing, nothing)
    return ContinuousData("test_data", df, layout, fs, AnalysisInfo(sample_rate = fs))
end

function create_test_continuous_data_empty_triggers(; n_samples = 1000, fs = 1000, seed::Union{Int,Nothing} = nothing)
    return _with_seed(seed) do
        t = collect(0:(n_samples-1)) ./ fs
        df = DataFrame(time = t, trigger = zeros(Int16, n_samples), channel1 = _synthetic_signal(t), channel2 = _synthetic_signal(t))
        layout = Layout(DataFrame(label = [:channel1, :channel2], inc = [0.0, 90.0], azi = [0.0, 0.0]), nothing, nothing, nothing)
        return ContinuousData("test_data", df, layout, fs, AnalysisInfo(sample_rate = fs))
    end
end

function create_test_lrp_data(;
    participant = 1,
    condition = 1,
    n_timepoints = 100,
    signal_scale = 1.0,
    fs = 250,
    seed::Union{Int,Nothing} = nothing,
)
    return _with_seed(seed) do
        time = collect(range(-0.2, 0.8, length = n_timepoints))
        df = DataFrame(
            time = time,
            sample = 1:n_timepoints,
            participant = fill(participant, n_timepoints),
            file = fill("test_file", n_timepoints),
        )
        df.C3 = _synthetic_signal(time, freq = 0.1, amp = signal_scale, noise = 0.01)
        df.C1 = _synthetic_signal(time, freq = 0.1, amp = signal_scale, noise = 0.01)
        df.Fp1 = _synthetic_signal(time, freq = 0.2, amp = signal_scale, noise = 0.01)
        df.C4 = _synthetic_signal(time, freq = 0.1, amp = signal_scale, noise = 0.01) .+ 0.1
        df.C2 = _synthetic_signal(time, freq = 0.1, amp = signal_scale, noise = 0.01) .+ 0.1
        df.Fp2 = _synthetic_signal(time, freq = 0.2, amp = signal_scale, noise = 0.01) .+ 0.1
        df.Fz = _synthetic_signal(time, freq = 0.15, amp = signal_scale, noise = 0.01)
        labels = [:C3, :C4, :C1, :C2, :Fp1, :Fp2, :Fz]
        layout = Layout(DataFrame(label = labels, inc = zeros(length(labels)), azi = zeros(length(labels))), nothing, nothing, nothing)
        return ErpData("test_data", condition, "condition_$condition", df, layout, fs, AnalysisInfo(sample_rate = fs), 10)
    end
end

function create_test_layout(; n_channels::Int = 4, layout_type::Symbol = :grid)
    if layout_type == :grid
        rows = ceil(Int, sqrt(max(n_channels, 1)))
        cols = ceil(Int, max(n_channels, 1) / rows)
        x2 = Vector{Float64}(undef, n_channels)
        y2 = Vector{Float64}(undef, n_channels)
        denom_x = max(1, cols - 1)
        denom_y = max(1, rows - 1)
        for i = 1:n_channels
            row = (i - 1) ÷ cols + 1
            col = (i - 1) % cols + 1
            x2[i] = (col - 1) / denom_x
            y2[i] = (row - 1) / denom_y
        end
        layout_data = DataFrame(
            :label => _generate_channel_labels(n_channels),
            :inc => zeros(n_channels),
            :azi => zeros(n_channels),
            :x2 => x2,
            :y2 => y2,
        )
        return Layout(layout_data, nothing, nothing, nothing)
    elseif layout_type == :topo
        layout_data = DataFrame(
            label = [:Fp1, :Fp2, :Cz, :Pz, :O1, :O2],
            inc = zeros(6),
            azi = zeros(6),
            x2 = [-0.3, 0.3, 0.0, 0.0, -0.3, 0.3],
            y2 = [0.9, 0.9, 0.0, -0.6, -0.9, -0.9],
            x3 = [-0.3, 0.3, 0.0, 0.0, -0.3, 0.3],
            y3 = [0.9, 0.9, 0.0, -0.6, -0.9, -0.9],
            z3 = [0.2, 0.2, 1.0, 0.8, 0.2, 0.2],
        )
        return Layout(layout_data, nothing, nothing, nothing)
    else
        @minimal_error("Unknown layout_type: $layout_type")
    end
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

function create_test_epoch_data_with_artifacts(;
    participant = 1,
    condition = 1,
    n_epochs = 10,
    n_timepoints = 100,
    n_channels = 3,
    fs = 1000,
    n_bad_epochs = 1,
    seed::Union{Int,Nothing} = nothing,
)
    return _with_seed(seed) do
        dfs = Vector{DataFrame}(undef, n_epochs)
        bad_indices = collect(1:n_bad_epochs)
        time = collect(0:(n_timepoints-1)) ./ fs
        labels = _generate_channel_labels(n_channels)
        for epoch = 1:n_epochs
            df = DataFrame(
                time = copy(time),
                sample = 1:n_timepoints,
                epoch = fill(epoch, n_timepoints),
                participant = fill(participant, n_timepoints),
            )
            is_bad = epoch <= n_bad_epochs
            for ch_name in labels
                signal = _synthetic_signal(time, freq = 10.0)
                if is_bad
                    signal .+= 50.0 .* randn(n_timepoints)
                    signal[1:min(20, n_timepoints)] .+= 200.0
                    n_timepoints > 70 && (signal[50:70] .+= 300.0)
                end
                df[!, ch_name] = signal
            end
            dfs[epoch] = df
        end
        epoch_data = EpochData(
            "test_data",
            condition,
            "condition_$condition",
            dfs,
            create_test_layout(n_channels = n_channels),
            fs,
            AnalysisInfo(sample_rate = fs),
        )
        return epoch_data, bad_indices
    end
end

"""
    generate_signal(n_trials, time_window, sample_rate, signal_freqs, signal_amps, signal_times, noise)

Generate synthetic trials containing windowed sinusoids plus Gaussian noise.
Returns `(times, signal)` where `signal` is `(n_samples × n_trials)`.
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
    dt = 1 / sample_rate
    x_time = range(time_window[1], time_window[2] - dt, step = dt)
    n_samples = length(x_time)
    time_zero = x_time[find_closest_time_index(x_time, 0.0)]

    base_signal = zeros(Float64, n_samples)
    for (freq, amp, times) in zip(signal_freqs, signal_amps, signal_times)
        i1 = searchsortedfirst(x_time, time_zero + times[1])
        i2 = searchsortedfirst(x_time, time_zero + times[2]) - 1
        two_pi_f = 2π * Float64(freq)
        amp_f = Float64(amp)
        @inbounds @simd for i = i1:i2
            base_signal[i] += amp_f * sin(two_pi_f * x_time[i])
        end
    end

    signal_trials = Matrix{Float64}(undef, n_samples, n_trials)
    if noise > 0
        Random.randn!(signal_trials)
        noise_f = Float64(noise)
        @inbounds for trial = 1:n_trials
            @simd for i = 1:n_samples
                signal_trials[i, trial] = muladd(signal_trials[i, trial], noise_f, base_signal[i])
            end
        end
    else
        @inbounds for trial = 1:n_trials
            @simd for i = 1:n_samples
                signal_trials[i, trial] = base_signal[i]
            end
        end
    end

    return x_time, signal_trials
end

function signal_to_data(
    times,
    signal,
    channel_name::Symbol,
    sample_rate::Real;
    file::String = "synthetic",
    condition::Int = 1,
    condition_name::String = "test",
)
    ndims(signal) == 1 && (signal = reshape(signal, :, 1))
    n_samples, n_trials = size(signal)
    times_vec = collect(Float64, times)
    length(times_vec) == n_samples || @minimal_error("Time vector length ($(length(times_vec))) must match signal samples ($n_samples)")

    layout = Layout(DataFrame(:label => [channel_name], :inc => [0.0], :azi => [0.0]), nothing, nothing, nothing)
    fs_int = round(Int, sample_rate)
    analysis_info = AnalysisInfo(sample_rate = Float64(sample_rate))

    if n_trials == 1
        erp_df = DataFrame(:time => times_vec, channel_name => signal[:, 1])
        return ErpData(file, condition, condition_name, erp_df, layout, fs_int, analysis_info, 1)
    end
    trial_dfs = [DataFrame(:time => copy(times_vec), channel_name => signal[:, t], :epoch => fill(t, n_samples)) for t = 1:n_trials]
    return EpochData(file, condition, condition_name, trial_dfs, layout, fs_int, analysis_info)
end

function create_test_tf_data(;
    participant::Int = 1,
    file::String = "participant$participant",
    condition::Int = 1,
    condition_name::String = "condition_$condition",
    n_channels::Int = 3,
    frequencies::Vector{Float64} = [4.0, 8.0, 12.0, 16.0],
    time_points::Vector{Float64} = [0.0, 0.1, 0.2, 0.3, 0.4],
    power_offset::Float64 = 0.0,
    noise::Float64 = 0.1,
    fs::Int = 256,
)
    labels = _generate_channel_labels(n_channels)
    n_freqs = length(frequencies)
    n_time = length(time_points)
    n_rows = n_freqs * n_time
    time_col = repeat(time_points, inner = n_freqs)
    freq_col = repeat(frequencies, outer = n_time)

    power_df = DataFrame(time = time_col, freq = freq_col)
    phase_df = DataFrame(time = copy(time_col), freq = copy(freq_col))
    for ch in labels
        power_df[!, ch] = power_offset .+ noise .* randn(n_rows)
        phase_df[!, ch] = zeros(n_rows)
    end

    return TimeFreqData(
        file,
        condition,
        condition_name,
        power_df,
        phase_df,
        create_test_layout(n_channels = n_channels),
        fs,
        :wavelet,
        nothing,
        AnalysisInfo(sample_rate = fs),
    )
end

function create_test_tf_epoch_data(;
    participant::Int = 1,
    file::String = "participant$participant",
    condition::Int = 1,
    condition_name::String = "condition_$condition",
    n_channels::Int = 3,
    n_epochs::Int = 10,
    frequencies::Vector{Float64} = [4.0, 8.0, 12.0, 16.0],
    time_points::Vector{Float64} = [0.0, 0.1, 0.2, 0.3, 0.4],
    power_offset::Float64 = 0.0,
    noise::Float64 = 0.1,
    fs::Int = 256,
)
    labels = _generate_channel_labels(n_channels)
    n_freqs = length(frequencies)
    n_time = length(time_points)
    n_rows = n_freqs * n_time
    time_col = repeat(time_points, inner = n_freqs)
    freq_col = repeat(frequencies, outer = n_time)
    layout = create_test_layout(n_channels = n_channels)

    power_dfs = Vector{DataFrame}(undef, n_epochs)
    phase_dfs = Vector{DataFrame}(undef, n_epochs)
    for epoch = 1:n_epochs
        power_df = DataFrame(time = copy(time_col), freq = copy(freq_col))
        phase_df = DataFrame(time = copy(time_col), freq = copy(freq_col))
        for ch in labels
            power_df[!, ch] = power_offset .+ noise .* randn(n_rows)
            phase_df[!, ch] = zeros(n_rows)
        end
        power_dfs[epoch] = power_df
        phase_dfs[epoch] = phase_df
    end

    return TimeFreqEpochData(
        file,
        condition,
        condition_name,
        power_dfs,
        phase_dfs,
        layout,
        fs,
        :wavelet,
        nothing,
        AnalysisInfo(sample_rate = fs),
    )
end

function create_synthetic_epochs(
    participant_id,
    condition_id,
    condition_name,
    n_epochs;
    n_timepoints = 500,
    channels = [:Fz, :Cz, :Pz],
    signal_strength = 1.0,
    noise_level = 0.3,
    sample_rate::Union{Int,Nothing} = nothing,
    seed::Union{Int,Nothing} = nothing,
)
    return _with_seed(seed) do
        time =
            isnothing(sample_rate) ? range(-0.2, 0.8, length = n_timepoints) : range(-0.2, step = 1.0 / sample_rate, length = n_timepoints)
        fs = isnothing(sample_rate) ? (n_timepoints > 1 ? round(Int, 1 / step(time)) : 500) : sample_rate
        n_ch = length(channels)
        ch_syms = collect(Symbol, channels)
        inc_vals = n_ch == 0 ? Float64[] : (n_ch == 1 ? [0.0] : collect(range(90.0, -90.0, length = n_ch)))
        layout = Layout(DataFrame(label = ch_syms, inc = inc_vals, azi = zeros(Float64, n_ch)), nothing, nothing, nothing)

        time_vec = collect(Float64, time)
        sample_vec = collect(1:n_timepoints)
        amp = condition_id == 1 ? Float64(signal_strength) : -Float64(signal_strength)
        inv_2sigma2 = 1.0 / (2 * 0.05^2)
        base_sig = [amp * exp(-((t - 0.2)^2) * inv_2sigma2) for t in time_vec]
        noise_f = Float64(noise_level)

        col_names = Vector{Symbol}(undef, 3 + n_ch)
        col_names[1] = :time
        col_names[2] = :sample
        col_names[3] = :epoch
        @inbounds for j = 1:n_ch
            col_names[3+j] = ch_syms[j]
        end

        epochs = Vector{DataFrame}(undef, n_epochs)
        for epoch = 1:n_epochs
            cols = Vector{Any}(undef, 3 + n_ch)
            cols[1] = copy(time_vec)
            cols[2] = copy(sample_vec)
            cols[3] = fill(epoch, n_timepoints)
            for j = 1:n_ch
                col = Vector{Float64}(undef, n_timepoints)
                Random.randn!(col)
                @inbounds @simd for i = 1:n_timepoints
                    col[i] = muladd(noise_f, col[i], base_sig[i])
                end
                cols[3+j] = col
            end
            epochs[epoch] = DataFrame(cols, col_names; copycols = false)
        end

        return EpochData(
            "participant_$(participant_id)",
            condition_id,
            condition_name,
            epochs,
            layout,
            fs,
            AnalysisInfo(sample_rate = Float64(fs)),
        )
    end
end
