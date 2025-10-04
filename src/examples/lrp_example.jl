"""
Example: Calculating Lateralized Readiness Potential (LRP)

This example demonstrates how to calculate the LRP from two ERP datasets
representing left-hand and right-hand responses.

The LRP is a measure of lateralized motor preparation derived from ERPs
recorded over motor cortex (typically C3/C4 electrode pair).
"""

using eegfun
using DataFrames

# Create example ERP data for left-hand responses
function create_example_left_erp()
    n_timepoints = 256
    sample_rate = 256.0
    time = collect(range(-0.5, 0.5, length = n_timepoints))

    # Simulate motor preparation: contralateral (right hemisphere) activation
    # is stronger than ipsilateral (left hemisphere) for left-hand responses
    c3_activity = 2.0 .+ 1.0 .* exp.(-((time .- 0.1) .^ 2) ./ 0.01)  # Ipsilateral: weaker
    c4_activity = 2.0 .+ 3.0 .* exp.(-((time .- 0.1) .^ 2) ./ 0.01)  # Contralateral: stronger

    # Add some noise
    c3_activity .+= 0.1 .* randn(n_timepoints)
    c4_activity .+= 0.1 .* randn(n_timepoints)

    # Create DataFrame
    df = DataFrame(
        time = time,
        sample = 1:n_timepoints,
        condition = fill(1, n_timepoints),
        condition_name = fill("left_hand", n_timepoints),
        participant = fill(1, n_timepoints),
        file = fill("example", n_timepoints),
        C3 = c3_activity,
        C4 = c4_activity,
    )

    # Create Layout
    layout = Layout(DataFrame(
        label = [:C3, :C4],
        inc = [0.0, 0.0],
        azi = [-45.0, 45.0],  # Left and right positions
    ), nothing, nothing)

    return ErpData(df, layout, sample_rate, AnalysisInfo(), 50)
end

# Create example ERP data for right-hand responses
function create_example_right_erp()
    n_timepoints = 256
    sample_rate = 256.0
    time = collect(range(-0.5, 0.5, length = n_timepoints))

    # Simulate motor preparation: contralateral (left hemisphere) activation
    # is stronger than ipsilateral (right hemisphere) for right-hand responses
    c3_activity = 2.0 .+ 3.0 .* exp.(-((time .- 0.1) .^ 2) ./ 0.01)  # Contralateral: stronger
    c4_activity = 2.0 .+ 1.0 .* exp.(-((time .- 0.1) .^ 2) ./ 0.01)  # Ipsilateral: weaker

    # Add some noise
    c3_activity .+= 0.1 .* randn(n_timepoints)
    c4_activity .+= 0.1 .* randn(n_timepoints)

    # Create DataFrame
    df = DataFrame(
        time = time,
        sample = 1:n_timepoints,
        condition = fill(2, n_timepoints),
        condition_name = fill("right_hand", n_timepoints),
        participant = fill(1, n_timepoints),
        file = fill("example", n_timepoints),
        C3 = c3_activity,
        C4 = c4_activity,
    )

    # Create Layout
    layout = Layout(DataFrame(
        label = [:C3, :C4],
        inc = [0.0, 0.0],
        azi = [-45.0, 45.0],  # Left and right positions
    ), nothing, nothing)

    return ErpData(df, layout, sample_rate, AnalysisInfo(), 50)
end

# Main example
function run_lrp_example()
    println("=" ^ 60)
    println("Lateralized Readiness Potential (LRP) Example")
    println("=" ^ 60)
    println()

    # Create example datasets
    println("Creating example ERP datasets...")
    erp_left = create_example_left_erp()
    erp_right = create_example_right_erp()

    println("  Left-hand response ERP:")
    println("    - Channels: ", channel_labels(erp_left))
    println("    - Time points: ", nrow(erp_left.data))
    println("    - Duration: ", duration(erp_left), " seconds")
    println()

    println("  Right-hand response ERP:")
    println("    - Channels: ", channel_labels(erp_right))
    println("    - Time points: ", nrow(erp_right.data))
    println("    - Duration: ", duration(erp_right), " seconds")
    println()

    # Calculate LRP
    println("Calculating LRP...")
    lrp_data = lrp(erp_left, erp_right)
    println()

    println("LRP Results:")
    println("  - Channels: ", channel_labels(lrp_data))
    println("  - Time points: ", nrow(lrp_data.data))
    println("  - Sample rate: ", lrp_data.sample_rate, " Hz")
    println("  - Number of epochs: ", lrp_data.n_epochs)
    println()

    # Show some LRP values
    println("Sample LRP values at different time points:")
    println("=" ^ 60)

    # Find peak LRP around expected motor preparation time (100 ms)
    time_idx_100ms = findfirst(x -> x >= 0.1, lrp_data.data.time)
    time_idx_0ms = findfirst(x -> x >= 0.0, lrp_data.data.time)
    time_idx_200ms = findfirst(x -> x >= 0.2, lrp_data.data.time)

    if !isnothing(time_idx_0ms)
        println("At t = 0 ms (stimulus onset):")
        println("  C3 LRP: ", round(lrp_data.data.C3[time_idx_0ms], digits = 3))
        println("  C4 LRP: ", round(lrp_data.data.C4[time_idx_0ms], digits = 3))
        println()
    end

    if !isnothing(time_idx_100ms)
        println("At t = 100 ms (expected motor preparation peak):")
        println("  C3 LRP: ", round(lrp_data.data.C3[time_idx_100ms], digits = 3))
        println("  C4 LRP: ", round(lrp_data.data.C4[time_idx_100ms], digits = 3))
        println()
    end

    if !isnothing(time_idx_200ms)
        println("At t = 200 ms:")
        println("  C3 LRP: ", round(lrp_data.data.C3[time_idx_200ms], digits = 3))
        println("  C4 LRP: ", round(lrp_data.data.C4[time_idx_200ms], digits = 3))
        println()
    end

    println("=" ^ 60)
    println("Interpretation:")
    println("  - Positive LRP values indicate greater contralateral activation")
    println("  - The LRP should peak around motor preparation time")
    println("  - Both C3 and C4 LRP should show similar patterns in symmetric data")
    println("=" ^ 60)
    println()

    return lrp_data
end

# Run the example if this file is executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    run_lrp_example()
end
