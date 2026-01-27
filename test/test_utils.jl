"""
Shared test utilities for EegFun.jl test suite

    Functions to create some synthetic data for testing
"""

using DataFrames
using EegFun
using Random

"""
    _get_rng(seed)

Helper to get a random number generator.
"""
function _get_rng(seed)
    if seed isa Integer
        return MersenneTwister(seed)
    elseif seed isa AbstractRNG
        return seed
    else
        return Random.GLOBAL_RNG
    end
end

"""
    _generate_channel_labels(n_channels::Int)

Helper to generate a vector of symbols for channel labels.
"""
_generate_channel_labels(n_channels::Int) = [Symbol("Ch$i") for i = 1:n_channels]

"""
    _synthesize_signal(t; freq = 10.0, amp = 1.0, noise = 0.1, rng = Random.GLOBAL_RNG)

Internal helper for generating test signals.
"""
function _synthesize_signal(t; freq = 10.0, amp = 1.0, noise = 0.1, rng = Random.GLOBAL_RNG)
    return amp .* sin.(2π .* freq .* t) .+ noise .* randn(rng, length(t))
end

"""
    _create_default_analysis_info()

Helper to create a default AnalysisInfo object.
"""
function _create_default_analysis_info()
    return EegFun.AnalysisInfo(:none, 0.0, 0.0)
end

"""
    create_test_data(; n::Int = 2000, fs::Int = 1000, n_channels::Int = 3, seed = nothing)

Create synthetic `ContinuousData` for testing.
"""
function create_test_data(; n = 2000, fs = 1000, n_channels = 3, seed = nothing)
    rng = _get_rng(seed)
    t = collect(0:(n-1)) ./ fs
    df = DataFrame(time = t, sample = 1:n, triggers = zeros(Int, n))

    channel_labels = _generate_channel_labels(n_channels)
    for ch in channel_labels
        df[!, ch] = 0.5 .+ _synthesize_signal(t, freq = 5 + rand(rng) * 10, amp = 1.0, noise = 0.1, rng = rng)
    end

    layout = create_test_layout(n_channels = n_channels)
    return EegFun.ContinuousData("test_data", copy(df, copycols = true), layout, fs, _create_default_analysis_info())
end

"""
    create_test_epoch_data(; n=2000, fs=1000, condition=1, n_epochs=10, n_channels=3, seed=nothing)

Create synthetic `EpochData` for a single condition. Returns `EpochData`.
"""
function create_test_epoch_data(; n = 2000, fs = 1000, condition = 1, n_epochs = 10, n_channels = 3, seed = nothing)
    rng = _get_rng(seed)
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
            df[!, ch] = _synthesize_signal(t, freq = 10.0, amp = condition * 0.5, noise = 0.1, rng = rng)
        end
        push!(dfs, df)
    end

    return EegFun.EpochData("test_data", condition, "condition_$(condition)", dfs, layout, fs, _create_default_analysis_info())
end

"""
    create_test_epoch_data_vector(; conditions=1:2, n_epochs=10, n_channels=3, kwargs...)

Create a `Vector{EpochData}` across multiple conditions.
"""
function create_test_epoch_data_vector(; conditions = 1:2, n_epochs = 10, n_channels = 3, fs = 1000, n = 2000, seed = nothing)
    return [
        create_test_epoch_data(condition = c, n_epochs = n_epochs, n_channels = n_channels, fs = fs, n = n, seed = seed) for c in conditions
    ]
end

"""
    create_test_erp_data(; participant = 1, condition = 1, fs = 1000, n_channels = 3, seed = nothing)

Create synthetic `ErpData` for testing.
"""
function create_test_erp_data(; participant = 1, condition = 1, fs = 1000, n_channels = 3, seed = nothing)
    rng = _get_rng(seed)
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
        df[!, ch] = _synthesize_signal(time, freq = 0.1, amp = condition * 0.5, noise = 0.1, rng = rng)
    end

    layout = create_test_layout(n_channels = n_channels)
    return EegFun.ErpData("test_data", condition, "condition_$condition", df, layout, fs, _create_default_analysis_info(), 10)
end

"""
    create_batch_test_erp_data(; n_conditions = 2, fs = 1000, n_channels = 3, seed = nothing)

Create a vector of ErpData objects for testing batch operations.
"""
function create_batch_test_erp_data(; n_conditions = 2, fs = 1000, n_channels = 3, seed = nothing)
    return [create_test_erp_data(condition = c, fs = fs, n_channels = n_channels, seed = seed) for c = 1:n_conditions]
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
    seed = nothing,
)
    rng = _get_rng(seed)
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
            rt = fill(rand(rng) * (rt_range[2] - rt_range[1]) + rt_range[1], n_timepoints),
        )

        for ch in channel_labels
            df[!, ch] = _synthesize_signal(time, freq = 10.0, amp = condition * 0.5, noise = 0.1, rng = rng)
        end
        push!(dfs, df)
    end

    layout = create_test_layout(n_channels = n_channels)
    return EegFun.EpochData("test_data", condition, "condition_$condition", dfs, layout, fs, _create_default_analysis_info())
end


"""
    create_test_continuous_data_with_triggers(; n = 1000, fs = 1000, seed = nothing)

Create synthetic `ContinuousData` with trigger sequences.
"""
function create_test_continuous_data_with_triggers(; n = 1000, fs = 1000, seed = nothing)
    rng = _get_rng(seed)
    t = collect(0:(n-1)) ./ fs
    x = _synthesize_signal(t, freq = 10.0, amp = 1.0, noise = 0.1, rng = rng)

    triggers = zeros(Int, n)
    triggers[300:302] .= [1, 7, 3]
    triggers[400:402] .= [1, 2, 3]
    triggers[750:752] .= [1, 2, 3]
    triggers[850] = 8
    triggers[900] = 9
    triggers[950:952] .= [1, 2, 3]

    df = DataFrame(time = t, triggers = triggers, A = x)
    layout = EegFun.Layout(DataFrame(label = [:A], inc = [0.0], azi = [0.0]), nothing, nothing)
    return EegFun.ContinuousData("test_data", copy(df, copycols = true), layout, fs, _create_default_analysis_info())
end


"""
    create_empty_trigger_data(; n_samples = 1000, fs = 1000, seed = nothing)

Create `ContinuousData` with zero triggers.
"""
function create_empty_trigger_data(; n_samples = 1000, fs = 1000, seed = nothing)
    rng = _get_rng(seed)
    t = collect(0:(n_samples-1)) ./ fs
    df = DataFrame(
        time = t,
        triggers = zeros(Int16, n_samples),
        channel1 = _synthesize_signal(t, rng = rng),
        channel2 = _synthesize_signal(t, rng = rng),
    )

    layout = EegFun.Layout(DataFrame(label = [:channel1, :channel2], inc = [0.0, 90.0], azi = [0.0, 0.0]), nothing, nothing)
    return EegFun.ContinuousData("test_data", df, layout, fs, _create_default_analysis_info())
end

"""
    create_test_lrp_erp(; participant = 1, condition = 1, n_timepoints = 100, signal_scale = 1.0, fs = 250, seed = nothing)

Create `ErpData` with lateralized channel pairs.
"""
function create_test_lrp_erp(; participant = 1, condition = 1, n_timepoints = 100, signal_scale = 1.0, fs = 250, seed = nothing)
    rng = _get_rng(seed)
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
    df.C3 = _synthesize_signal(time, freq = 0.1, amp = signal_scale, noise = 0.01, rng = rng)
    df.C1 = _synthesize_signal(time, freq = 0.1, amp = signal_scale, noise = 0.01, rng = rng)
    df.Fp1 = _synthesize_signal(time, freq = 0.2, amp = signal_scale, noise = 0.01, rng = rng)

    # Right hemisphere (shifted)
    df.C4 = _synthesize_signal(time, freq = 0.1, amp = signal_scale, noise = 0.01, rng = rng) .+ 0.1
    df.C2 = _synthesize_signal(time, freq = 0.1, amp = signal_scale, noise = 0.01, rng = rng) .+ 0.1
    df.Fp2 = _synthesize_signal(time, freq = 0.2, amp = signal_scale, noise = 0.01, rng = rng) .+ 0.1

    df.Fz = _synthesize_signal(time, freq = 0.15, amp = signal_scale, noise = 0.01, rng = rng)

    channel_labels = [:C3, :C4, :C1, :C2, :Fp1, :Fp2, :Fz]
    layout = EegFun.Layout(
        DataFrame(label = channel_labels, inc = zeros(length(channel_labels)), azi = zeros(length(channel_labels))),
        nothing,
        nothing,
    )

    return EegFun.ErpData("test_data", condition, "condition_$condition", df, layout, fs, _create_default_analysis_info(), 10)
end



function create_test_layout(; n_channels::Int = 4, layout_type::Symbol = :grid)
    if layout_type == :grid
        # Create a simple test layout with n_channels
        channels = [Symbol("Ch$i") for i = 1:n_channels]
        # Create grid positions
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
            :label => channels,
            :inc => zeros(n_channels),
            :azi => zeros(n_channels),
            :x2 => [pos[1] for pos in positions],
            :y2 => [pos[2] for pos in positions],
        )
    elseif layout_type == :topo
        # Create a topographic layout for 6 channels
        channels = [:Fp1, :Fp2, :Cz, :Pz, :O1, :O2]
        layout_data = DataFrame(
            label = channels,
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
    df = DataFrame(
        channel = [:Ch1, :Ch2, :Ch3, :Ch4, :Ch5],
        min = [-2.1, -1.8, -0.5, -0.3, -0.2],
        max = [2.3, 1.9, 0.8, 0.4, 0.3],
        std = [1.2, 1.1, 0.4, 0.2, 0.15],
        var = [1.44, 1.21, 0.16, 0.04, 0.0225],
        range = [4.4, 3.7, 1.3, 0.7, 0.5],
        zvar = [1.2, 0.8, -0.5, -1.1, -1.4],
    )
    return df
end

function create_test_summary_data_with_epochs()
    df = DataFrame(
        channel = repeat([:Ch1, :Ch2, :Ch3], 3),
        epoch = repeat([1, 2, 3], 3),
        min = [-2.1, -1.8, -0.5, -2.0, -1.7, -0.4, -2.2, -1.9, -0.6],
        max = [2.3, 1.9, 0.8, 2.2, 1.8, 0.7, 2.4, 2.0, 0.9],
        std = [1.2, 1.1, 0.4, 1.1, 1.0, 0.3, 1.3, 1.2, 0.5],
        var = [1.44, 1.21, 0.16, 1.21, 1.0, 0.09, 1.69, 1.44, 0.25],
        range = [4.4, 3.7, 1.3, 4.2, 3.5, 1.1, 4.6, 3.9, 1.5],
        zvar = [1.2, 0.8, -0.5, 0.9, 0.6, -0.8, 1.5, 1.1, -0.2],
    )
    return df
end


"""
    create_epoch_test_data(; fs = 1000, seed = nothing)

Create synthetic `EpochData` with 3 channels and 3 epochs.
"""
function create_epoch_test_data(; fs = 1000, seed = nothing)
    rng = _get_rng(seed)
    time = collect(range(-0.2, 0.8, length = 100))
    n = length(time)
    channel_labels = [:A, :B, :C]

    dfs = DataFrame[]
    for epoch = 1:3
        df = DataFrame(time = time, epoch = fill(epoch, n))
        df.A = _synthesize_signal(time, freq = 10.0, rng = rng)
        df.B = _synthesize_signal(time, freq = 10.0, rng = rng)
        df.C = _synthesize_signal(time, freq = 20.0, rng = rng)
        push!(dfs, df)
    end

    layout = EegFun.Layout(DataFrame(label = channel_labels, inc = [0.0, 0.0, 0.0], azi = [0.0, 0.0, 0.0]), nothing, nothing)
    return EegFun.EpochData("test_data", 1, "condition_1", dfs, layout, fs, _create_default_analysis_info())
end


"""
    create_test_epochs_with_artifacts(; participant=1, condition=1, n_epochs=10, n_timepoints=100, n_channels=3, n_bad_epochs=1, seed=nothing)

Create `EpochData` with large amplitude noise in some epochs for testing artifact detection.
"""
function create_test_epochs_with_artifacts(;
    participant = 1,
    condition = 1,
    n_epochs = 10,
    n_timepoints = 100,
    n_channels = 3,
    fs = 1000,
    n_bad_epochs = 1,
    seed = nothing,
)
    rng = _get_rng(seed)
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
            signal = _synthesize_signal(time, freq = 10.0, rng = rng)

            if is_bad
                signal .+= 50.0 .* randn(rng, n_timepoints)
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
    epoch_data = EegFun.EpochData("test_data", condition, "condition_$condition", dfs, layout, fs, _create_default_analysis_info())

    return epoch_data, bad_indices
end
