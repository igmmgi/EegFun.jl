"""
Shared test utilities for eegfun.jl test suite
"""

using DataFrames
using eegfun
using Random

function create_test_data(; n::Int = 2000, fs::Int = 1000, n_channels::Int = 3)

    t = collect(0:(n-1)) ./ fs
    # Create DataFrame with time, triggers, and sample
    df = DataFrame(time = t, sample = 1:n, triggers = zeros(Int, n))

    # Add random channel data
    channel_labels = [Symbol("Ch$i") for i = 1:n_channels]
    for ch in channel_labels
        df[!, ch] = 0.5 .+ sin.(2π .* (5 + rand() * 10) .* t) .+ 0.1 .* randn(length(t))
    end

    # just fake layout
    layout = eegfun.Layout(
        DataFrame(label = channel_labels, inc = zeros(length(channel_labels)), azi = zeros(length(channel_labels))),
        nothing,
        nothing,
    )

    dat = eegfun.ContinuousData("test_data", copy(df, copycols = true), layout, fs, eegfun.AnalysisInfo())
    return dat
end

# Helper function to create test epoch data
function create_test_epoch_data(;
    n::Int = 2000,
    fs::Int = 1000,
    conditions::Int = 1,
    n_epochs::Int = 10,
    n_channels::Int = 3,
)
    # Create conditions from 1 to condition
    conditions = 1:conditions

    # Create time vector
    t = collect(0:(n-1)) ./ fs

    # Generate channel labels based on n_channels
    channel_labels = [Symbol("Ch$i") for i = 1:n_channels]

    # Create layout (same for all conditions)
    layout = eegfun.Layout(
        DataFrame(label = channel_labels, inc = zeros(length(channel_labels)), azi = zeros(length(channel_labels))),
        nothing,
        nothing,
    )

    # Create analysis info
    analysis_info = eegfun.AnalysisInfo(:none, 0.0, 0.0)

    # Create EpochData for each condition
    condition_data = eegfun.EpochData[]
    for condition in conditions
        dfs = DataFrame[]
        for epoch = 1:n_epochs
            df = DataFrame(time = t, sample = 1:n, epoch = fill(epoch, n))

            # Add EEG channels with some signal structure
            for ch in channel_labels
                # Create some signal with condition-specific amplitude and epoch-specific noise
                signal = sin.(2π * 0.1 * t) .* (condition * 0.5) .+ randn(n) * 0.1
                df[!, ch] = signal
            end

            push!(dfs, df)
        end

        push!(
            condition_data,
            eegfun.EpochData("test_data", condition, "condition_$(condition)", dfs, layout, fs, analysis_info),
        )
    end

    # Return single EpochData if only one condition, otherwise return Vector{EpochData}
    return length(condition_data) == 1 ? condition_data[1] : condition_data
end






# Helper function to create test ERP data
function create_test_erp_data(participant::Int, condition::Int; fs::Int = 1000, n_channels::Int = 3)
    # Create DataFrame with columns in correct order
    time = collect(-0.5:(1/fs):2.0)
    df = DataFrame(
        time = time,
        sample = 1:length(time),
        condition = fill(condition, length(time)),
        condition_name = fill("condition_$condition", length(time)),
        participant = fill(participant, length(time)),
    )

    # chanel labels
    channel_labels = [Symbol("Ch$i") for i = 1:n_channels]
    for ch in channel_labels
        df[!, ch] = sin.(2π * 0.1 * time) .* (condition * 0.5) .+ randn(length(time)) * 0.1
    end

    # Create layout
    layout = eegfun.Layout(
        DataFrame(label = channel_labels, inc = zeros(length(channel_labels)), azi = zeros(length(channel_labels))),
        nothing,
        nothing,
    )

    # Create analysis info
    analysis_info = eegfun.AnalysisInfo(:none, 0.0, 0.0)

    return eegfun.ErpData("test_data", condition, "condition_$condition", df, layout, fs, analysis_info, 10)
end


# Helper function to create test ERP data for batch processing
function create_batch_test_erp_data(n_conditions::Int = 2; fs::Int = 1000, n_channels::Int = 3)
    erps = eegfun.ErpData[]
    for cond = 1:n_conditions
        erp_data = create_test_erp_data(1, cond, fs = fs, n_channels = n_channels)
        push!(erps, erp_data)
    end
    return erps
end



#####################################




# Helper function to create test epoch data with reaction times
function create_test_epoch_data_with_rt(
    participant::Int,
    condition::Int,
    n_epochs::Int = 5,
    n_timepoints::Int = 100,
    n_channels::Int = 3;
    epoch_start::Float64 = -0.2,
    epoch_end::Float64 = 0.8,
    rt_range::Tuple{Float64,Float64} = (0.2, 0.6),
)
    # Create time vector using the provided epoch_start and epoch_end
    time = collect(range(epoch_start, epoch_end, length = n_timepoints))

    # Create channel data with some condition-specific patterns
    channel_data = Dict{Symbol,Any}()
    channel_data[:time] = time

    # Add metadata columns
    channel_data[:condition] = fill(condition, n_timepoints)
    channel_data[:condition_name] = fill("condition_$condition", n_timepoints)
    channel_data[:participant] = fill(participant, n_timepoints)

    # Generate channel labels based on n_channels
    channel_labels = [Symbol("Ch$i") for i = 1:n_channels]

    # Add EEG channels with condition-specific patterns
    for (i, ch) in enumerate(channel_labels)
        # Create some signal with condition-specific amplitude
        signal = sin.(2π * 0.1 * time) .* (condition * 0.5) .+ randn(n_timepoints) * 0.1
        channel_data[ch] = signal
    end

    # Create multiple epochs with varying reaction times
    dfs = DataFrame[]
    for epoch = 1:n_epochs
        df = DataFrame()
        df.time = channel_data[:time]
        df.condition = channel_data[:condition]
        df.condition_name = channel_data[:condition_name]
        df.participant = channel_data[:participant]
        df.epoch = fill(epoch, n_timepoints)

        # Generate random reaction time within the specified range
        rt = rand() * (rt_range[2] - rt_range[1]) + rt_range[1]
        df.rt = fill(rt, n_timepoints)

        for ch in channel_labels
            df[!, ch] = channel_data[ch] .+ randn(n_timepoints) * 0.05  # Add some epoch-specific noise
        end

        push!(dfs, df)
    end

    # Create layout
    layout = eegfun.Layout(
        DataFrame(label = channel_labels, inc = zeros(length(channel_labels)), azi = zeros(length(channel_labels))),
        nothing,
        nothing,
    )

    # Create analysis info
    analysis_info = eegfun.AnalysisInfo(:none, 0.0, 0.0)

    return eegfun.EpochData("test_data", 1, "condition_1", dfs, layout, 1000, analysis_info)
end

# Helper function to create continuous data with triggers
function create_test_continuous_data_with_triggers(; n::Int = 1000, fs::Int = 1000)
    t = collect(0:(n-1)) ./ fs
    # Create some test signal
    x = sin.(2π .* 10 .* t) .+ 0.1 .* randn(length(t))

    # Create triggers at specific points - create exactly 2 sequences [1,2,3]
    triggers = zeros(Int, n)
    # First sequence [1,2,3] starting at sample 400 - consecutive
    triggers[400:402] .= [1, 2, 3]  # Consecutive sequence [1,2,3]
    # Second sequence [1,2,3] starting at sample 750 - consecutive
    triggers[750:752] .= [1, 2, 3]  # Consecutive sequence [1,2,3]
    # Some individual triggers for other tests
    # triggers[700] = 1  # Individual trigger 1 (removed to avoid sequence conflicts)
    # triggers[800] = 2  # Individual trigger 2 (removed to avoid sequence conflicts)
    triggers[850] = 8  # Individual trigger 8 (for position constraints test)
    triggers[900] = 9  # Individual trigger 9 (for position constraints test)
    # Add some wildcard sequences [1,7,3] for testing - consecutive
    triggers[300:302] .= [1, 7, 3]  # Consecutive sequence [1,7,3]
    # Add another sequence [1,2,3] after trigger 9 for position constraints test
    triggers[950:952] .= [1, 2, 3]  # Consecutive sequence [1,2,3] after trigger 9

    df = DataFrame(time = t, triggers = triggers, A = x)
    layout = eegfun.Layout(DataFrame(label = [:A], inc = [0.0], azi = [0.0]), nothing, nothing)
    dat = eegfun.ContinuousData("test_data", copy(df, copycols = true), layout, fs, eegfun.AnalysisInfo())
    return dat
end


function create_empty_trigger_data(; n_samples::Int = 1000, fs::Int = 1000)
    """Create data with no triggers for edge case testing"""
    time = collect(0:(n_samples-1)) ./ fs
    triggers = zeros(Int16, n_samples)

    df = DataFrame(time = time, triggers = triggers, channel1 = randn(n_samples), channel2 = randn(n_samples))

    layout =
        eegfun.Layout(DataFrame(label = [:channel1, :channel2], inc = [0.0, 90.0], azi = [0.0, 0.0]), nothing, nothing)

    dat = eegfun.ContinuousData("test_data", df, layout, fs, eegfun.AnalysisInfo())
    return dat
end

# Helper function for LRP testing with lateral channel pairs
function create_test_lrp_erp(participant::Int, condition::Int, n_timepoints::Int = 100, signal_scale::Float64 = 1.0)
    # Create time vector
    time = collect(range(-0.2, 0.8, length = n_timepoints))

    # Create channel data with lateral pairs (C3/C4, C1/C2, Fp1/Fp2)
    df = DataFrame()
    df.time = time
    df.sample = 1:n_timepoints
    df.condition = fill(condition, n_timepoints)
    df.condition_name = fill("condition_$condition", n_timepoints)
    df.participant = fill(participant, n_timepoints)
    df.file = fill("test_file", n_timepoints)

    # Add lateral channel pairs with known patterns
    # Left hemisphere (odd numbers)
    df.C3 = signal_scale .* sin.(2π * 0.1 * time) .+ 0.01 .* randn(n_timepoints)
    df.C1 = signal_scale .* cos.(2π * 0.1 * time) .+ 0.01 .* randn(n_timepoints)
    df.Fp1 = signal_scale .* sin.(2π * 0.2 * time) .+ 0.01 .* randn(n_timepoints)

    # Right hemisphere (even numbers)
    df.C4 = signal_scale .* sin.(2π * 0.1 * time .+ π / 4) .+ 0.01 .* randn(n_timepoints)
    df.C2 = signal_scale .* cos.(2π * 0.1 * time .+ π / 4) .+ 0.01 .* randn(n_timepoints)
    df.Fp2 = signal_scale .* sin.(2π * 0.2 * time .+ π / 4) .+ 0.01 .* randn(n_timepoints)

    # Add midline channel (should not be detected as a pair)
    df.Fz = signal_scale .* sin.(2π * 0.15 * time) .+ 0.01 .* randn(n_timepoints)

    # Create Layout with lateral channels
    channel_labels = [:C3, :C4, :C1, :C2, :Fp1, :Fp2, :Fz]
    layout_df =
        DataFrame(label = channel_labels, inc = zeros(length(channel_labels)), azi = zeros(length(channel_labels)))
    layout = eegfun.Layout(layout_df, nothing, nothing)

    # Create AnalysisInfo
    analysis_info = eegfun.AnalysisInfo()

    # Create ErpData
    return eegfun.ErpData("test_data", condition, "condition_$condition", df, layout, 250, analysis_info, 10)
end



# Helper function for creating test layouts
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
            :x => [pos[1] for pos in positions],
            :y => [pos[2] for pos in positions],
            :x2 => [pos[1] for pos in positions],  # For compatibility
            :y2 => [pos[2] for pos in positions],   # For compatibility
        )
    elseif layout_type == :topo
        # Create a topographic layout for 6 channels
        channels = [:Fp1, :Fp2, :Cz, :Pz, :O1, :O2]
        layout_data = DataFrame(
            label = channels,
            x2 = [0.0, 0.0, 0.0, 0.0, -0.5, 0.5],
            y2 = [1.0, 1.0, 0.0, -1.0, -1.0, -1.0],
            x3 = [0.0, 0.0, 0.0, 0.0, -0.5, 0.5],
            y3 = [1.0, 1.0, 0.0, -1.0, -1.0, -1.0],
            z3 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        )
    else
        error("Unknown layout_type: $layout_type")
    end

    return eegfun.Layout(layout_data, nothing, nothing)
end

# Helper functions for creating test summary data
function create_test_summary_data()
    # Create a DataFrame with channel summary statistics
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
    # Create a DataFrame with epoch information for testing averaging
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

# Helper function to create test epoch data with multiple conditions

# Helper function to create test epoch data for channel_summary tests
function create_epoch_test_data()
    """Create EpochData with 3 channels (A, B, C) and multiple epochs for channel_summary tests"""
    # Create time vector
    time = collect(range(-0.2, 0.8, length = 100))

    # Create multiple epochs
    dfs = DataFrame[]
    for epoch = 1:3
        df = DataFrame()
        df.time = time
        df.epoch = fill(epoch, 100)

        # Create 3 channels with different patterns
        df.A = sin.(2π * 0.1 * time) .+ 0.1 .* randn(100)
        df.B = cos.(2π * 0.1 * time) .+ 0.1 .* randn(100)
        df.C = sin.(2π * 0.2 * time) .+ 0.1 .* randn(100)

        push!(dfs, df)
    end

    # Create layout
    layout =
        eegfun.Layout(DataFrame(label = [:A, :B, :C], inc = [0.0, 0.0, 0.0], azi = [0.0, 0.0, 0.0]), nothing, nothing)

    # Create analysis info
    analysis_info = eegfun.AnalysisInfo(:none, 0.0, 0.0)

    # Create EpochData
    return eegfun.EpochData("test_data", 1, "condition_1", dfs, layout, 1000, analysis_info)
end

"""
    create_test_epochs_with_artifacts(participant::Int, condition::Int, n_epochs::Int, n_timepoints::Int, n_channels::Int; n_bad_epochs::Int = 1)

Create test epoch data with artifacts for testing artifact detection.
"""
function create_test_epochs_with_artifacts(
    participant::Int,
    condition::Int,
    n_epochs::Int,
    n_timepoints::Int,
    n_channels::Int;
    n_bad_epochs::Int = 1,
)
    dfs = DataFrame[]
    bad_indices = Int[]

    # Create time vector
    time = collect(0:(n_timepoints-1)) ./ 1000  # Convert to seconds

    for epoch = 1:n_epochs
        df = DataFrame()
        df.time = time
        df.epoch = fill(epoch, n_timepoints)
        df.participant = fill(participant, n_timepoints)
        df.condition = fill(condition, n_timepoints)
        df.condition_name = fill("condition_$condition", n_timepoints)

        # Add artifacts to some epochs (only the first n_bad_epochs)
        is_bad_epoch = epoch <= n_bad_epochs
        if is_bad_epoch
            push!(bad_indices, epoch)
        end

        # Create channels with different patterns
        for ch = 1:n_channels
            channel_name = Symbol("Ch$ch")

            # Create base signal
            base_signal = sin.(2π * 0.1 * time) .+ 0.1 .* randn(n_timepoints)

            # Add artifacts to bad epochs
            if is_bad_epoch
                # Create extremely strong artifacts that will definitely be detected
                base_signal .+= 50.0 .* randn(n_timepoints)  # Very high amplitude noise
                base_signal[1:min(20, n_timepoints)] .+= 200.0  # Extremely high amplitude at start
                if n_timepoints > 19
                    base_signal[(end-19):end] .+= 200.0  # Extremely high amplitude at end
                end
                if n_timepoints >= 70
                    base_signal[50:70] .+= 300.0  # Extremely high amplitude in middle
                elseif n_timepoints >= 50
                    base_signal[50:end] .+= 300.0  # Extremely high amplitude in middle
                end
            end

            df[!, channel_name] = base_signal
        end

        push!(dfs, df)
    end

    # Create layout
    layout = eegfun.Layout(
        DataFrame(
            label = [Symbol("Ch$i") for i = 1:n_channels],
            inc = fill(0.0, n_channels),
            azi = fill(0.0, n_channels),
        ),
        nothing,
        nothing,
    )

    # Create analysis info
    analysis_info = eegfun.AnalysisInfo(:none, 0.0, 0.0)

    # Create EpochData
    epoch_data = eegfun.EpochData("test_data", condition, "condition_$condition", dfs, layout, 1000, analysis_info)

    return epoch_data, bad_indices
end
