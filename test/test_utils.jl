"""
Shared test utilities for eegfun.jl test suite
"""

using DataFrames
using eegfun

# Helper function to create test ERP data
function create_test_erp_data(participant::Int, condition::Int, n_timepoints::Int = 100, n_channels::Int = 3)
    # Create time vector
    time = collect(range(-0.2, 0.8, length = n_timepoints))

    # Create channel data with some condition-specific patterns
    channel_data = Dict{Symbol,Any}()
    channel_data[:time] = time

    # Add metadata columns
    channel_data[:condition] = fill(condition, n_timepoints)
    channel_data[:condition_name] = fill("condition_$condition", n_timepoints)
    channel_data[:participant] = fill(participant, n_timepoints)

    # Add EEG channels with condition-specific patterns
    for (i, ch) in enumerate([:Fz, :Cz, :Pz][1:min(n_channels, 3)])
        # Create some signal with condition-specific amplitude
        signal = sin.(2π * 0.1 * time) .* (condition * 0.5) .+ randn(n_timepoints) * 0.1
        channel_data[ch] = signal
    end

    # Create DataFrame with columns in correct order
    df = DataFrame()
    df.time = channel_data[:time]
    df.condition = channel_data[:condition]
    df.condition_name = channel_data[:condition_name]
    df.participant = channel_data[:participant]
    for ch in [:Fz, :Cz, :Pz][1:min(n_channels, 3)]
        df[!, ch] = channel_data[ch]
    end

    # Create Layout with required columns
    channel_labels = [:Fz, :Cz, :Pz][1:min(n_channels, 3)]
    layout_df = DataFrame(
        label = channel_labels,
        inc = fill(0.0, length(channel_labels)),  # incidence angles
        azi = fill(0.0, length(channel_labels)),   # azimuth angles
    )
    layout = eegfun.Layout(layout_df, nothing, nothing)

    # Create AnalysisInfo
    analysis_info = eegfun.AnalysisInfo()

    # Create ErpData
    return eegfun.ErpData(df, layout, 250.0, analysis_info, 10)
end

# Helper function to create test epoch data
function create_test_epoch_data(
    participant::Int,
    condition::Int,
    n_epochs::Int = 5,
    n_timepoints::Int = 100,
    n_channels::Int = 3,
)
    epochs = []

    for epoch = 1:n_epochs
        # Create time vector
        time = collect(range(-0.2, 0.8, length = n_timepoints))

        # Create channel data
        channel_data = Dict{Symbol,Any}()
        channel_data[:time] = time
        channel_data[:condition] = fill(condition, n_timepoints)
        channel_data[:condition_name] = fill("condition_$condition", n_timepoints)
        channel_data[:participant] = fill(participant, n_timepoints)
        channel_data[:epoch] = fill(epoch, n_timepoints)

        # Add EEG channels
        for (i, ch) in enumerate([:Fz, :Cz, :Pz][1:min(n_channels, 3)])
            signal = sin.(2π * 0.1 * time) .* (condition * 0.5) .+ randn(n_timepoints) * 0.1
            channel_data[ch] = signal
        end

        # Create DataFrame with columns in correct order
        df = DataFrame()
        df.time = channel_data[:time]
        df.condition = channel_data[:condition]
        df.condition_name = channel_data[:condition_name]
        df.participant = channel_data[:participant]
        df.epoch = channel_data[:epoch]
        for ch in [:Fz, :Cz, :Pz][1:min(n_channels, 3)]
            df[!, ch] = channel_data[ch]
        end
        push!(epochs, df)
    end

    # Create EpochData
    layout = eegfun.Layout(
        DataFrame(label = [:Fz, :Cz, :Pz][1:min(n_channels, 3)], inc = [0.0, 0.0, 0.0], azi = [0.0, 0.0, 0.0]),
        nothing,
        nothing,
    )
    return eegfun.EpochData(epochs, layout, 250.0, eegfun.AnalysisInfo())
end
