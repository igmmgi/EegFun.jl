#!/usr/bin/env julia
"""
Manual test for MVPA/Decoding functionality

This script tests the decoding workflow following erplab's approach:
1. Create test epoch data for multiple conditions
2. Decode for each participant separately
3. Average across participants (grand average)
4. Plot results

Tests both binary and multi-class classification with different models.
"""

using eegfun
using DataFrames
using Random

# Set seed for reproducibility
Random.seed!(42)

println("\n" * "="^70)
println("MVPA/DECODING MANUAL TEST")
println("="^70)

# ============================================================================
# SETUP: Create Test Data
# ============================================================================

println("\n[1/6] Creating test epoch data...")

# Create test epoch data for multiple participants and conditions
# Each participant will have epochs for 2 conditions (binary classification)
n_participants = 3
n_epochs_per_condition = 30  # Enough for cross-validation
n_timepoints = 200  # -0.2 to 0.8s at 1000Hz
sample_rate = 1000
channels = [:Fz, :Cz, :Pz, :POz]  # 4 channels for MVPA

# Create layout (must match channel names exactly)
layout = eegfun.Layout(
    DataFrame(
        label = channels,  # [:Fz, :Cz, :Pz, :POz]
        inc = [90.0, 0.0, -90.0, -90.0],  # Rough positions
        azi = [0.0, 0.0, 0.0, 0.0],
    ),
    nothing,
    nothing,
)

analysis_info = eegfun.AnalysisInfo(:none, 0.0, 0.0)

# Create time vector
time = collect(range(-0.2, 0.8, length = n_timepoints))

# Function to create epochs for one participant and one condition
function create_participant_condition_epochs(
    participant_id::Int,
    condition_id::Int,
    condition_name::String,
    n_epochs::Int,
    signal_strength::Float64,
)
    epochs = DataFrame[]
    for epoch = 1:n_epochs
        df = DataFrame(time = time, epoch = fill(epoch, n_timepoints))

        # Create condition-specific signal patterns
        # Condition 1: Positive deflection around 200ms
        # Condition 2: Negative deflection around 200ms
        for (ch_idx, ch) in enumerate(channels)
            # Base signal with condition-specific pattern
            if condition_id == 1
                # Condition 1: Positive ERP component
                signal = signal_strength * exp.(-((time .- 0.2).^2) / (2 * 0.05^2)) .+
                        0.3 * signal_strength * exp.(-((time .- 0.4).^2) / (2 * 0.08^2))
            else
                # Condition 2: Negative ERP component
                signal = -signal_strength * exp.(-((time .- 0.2).^2) / (2 * 0.05^2)) .+
                        0.2 * signal_strength * exp.(-((time .- 0.5).^2) / (2 * 0.1^2))
            end

            # Add channel-specific variation and noise
            channel_factor = 1.0 + (ch_idx - 1) * 0.1
            noise = 0.3 * randn(n_timepoints)  # SNR ~ 3:1
            df[!, ch] = signal .* channel_factor .+ noise
        end

        push!(epochs, df)
    end

    return eegfun.EpochData(
        "participant_$(participant_id)",
        condition_id,
        condition_name,
        epochs,
        layout,
        sample_rate,
        analysis_info,
    )
end

# Create epochs for all participants
all_participant_epochs = []
for p in 1:n_participants
    # Participant p, Condition 1
    cond1 = create_participant_condition_epochs(p, 1, "Face", n_epochs_per_condition, 0.05)
    # Participant p, Condition 2
    cond2 = create_participant_condition_epochs(p, 2, "Object", n_epochs_per_condition, 0.05)
    push!(all_participant_epochs, [cond1, cond2])
end

println("  ✓ Created epoch data for $n_participants participants")
println("  ✓ Each participant has $(n_epochs_per_condition) epochs per condition")
println("  ✓ Conditions: Face vs Object")

# ==========================================
# TEST 1: Binary Classification with Default 
# ==========================================

# model_method = :logistic
# model_method = :svm
model_method = :lda

# Decode for Participant 1 - default model is used automatically!
epochs_p1 = all_participant_epochs[1]
decoded_p1 = eegfun.decode(
    epochs_p1,
    channels;
    model = model_method,
    time_range = (-0.2, 0.8),
    n_iterations = 10,  # Reduced for testing
    n_folds = 3,
    equalize_trials = true,
)

# Decode for Participant 2 - default model used automatically
epochs_p2 = all_participant_epochs[2]
decoded_p2 = eegfun.decode(
    epochs_p2,
    channels;
    model = model_method,
    time_range = (-0.2, 0.8),
    n_iterations = 10,
    n_folds = 3,
    equalize_trials = true,
)

# Decode for Participant 3 - default model used automatically
epochs_p3 = all_participant_epochs[3]
decoded_p3 = eegfun.decode(
    epochs_p3,
    channels;
    model = model_method,
    time_range = (-0.2, 0.8),
    n_iterations = 10,
    n_folds = 3,
    equalize_trials = true,
)

# Grand average
all_decoded = [decoded_p1, decoded_p2, decoded_p3]
grand_avg = eegfun.grand_average(all_decoded)

# Plot grand average
fig1 = eegfun.plot_decoding(grand_avg, title = "Grand Average: Face vs Object ($(model_method))")



# Get some basic data with initial preprocessing steps (high-pass filter, epoch)
data_file = joinpath(@__DIR__, "..", "..", "..", "AttentionExp", "recoded", "Flank_C_11.bdf")
layout_file = eegfun.read_layout("./data/layouts/biosemi/biosemi72.csv");
eegfun.polar_to_cartesian_xy!(layout_file)
dat = eegfun.read_bdf(data_file);
dat = eegfun.create_eeg_dataframe(dat, layout_file);
eegfun.rereference!(dat, :avg)
eegfun.filter_data!(dat, "hp", 1)
eegfun.is_extreme_value!(dat, 500);
eegfun.mark_epoch_windows!(dat, [1, 2, 3, 4], [-0.5, 3.0]) # simple epoch marking with trigger 1 and 3

# EPOCHS
epoch_cfg = [
    eegfun.EpochCondition(name = "ExampleEpoch1", trigger_sequences = [[1], [3]]),
    eegfun.EpochCondition(name = "ExampleEpoch2", trigger_sequences = [[2], [4]]),
]
epochs = eegfun.extract_epochs(dat, epoch_cfg, -0.5, 3)
# Decode for Participant 1 - default model is used automatically!
epochs_p1 = epochs[1:2]

decoded_p1 = eegfun.decode(
    epochs_p1,
    channels;
    model = model_method,
    time_range = (-0.5, 3.0),
    n_iterations = 100,  # Reduced for testing
    n_folds = 3,
    equalize_trials = true,
)
fig1 = eegfun.plot_decoding(decoded_p1, title = "Participant 1: Face vs Object ($(model_method))")



using eegfun
using Glob  # For finding files (or use readdir with filter)

# Configuration
data_dir = joinpath(@__DIR__, "..", "..", "..", "AttentionExp", "recoded")
layout_file = eegfun.read_layout("./data/layouts/biosemi/biosemi72.csv")
eegfun.polar_to_cartesian_xy!(layout_file)

# Epoch configuration
epoch_cfg = [
    eegfun.EpochCondition(name = "ExampleEpoch1", trigger_sequences = [[1], [3]]),
    eegfun.EpochCondition(name = "ExampleEpoch2", trigger_sequences = [[2], [4]]),
]

# Decoding parameters
model_method = :logistic  # or :svm, :lda
time_range = (-0.5, 3.0)
n_iterations = 100
n_folds = 3

# Find all .bdf files in directory
bdf_files = glob("*.bdf", data_dir)
# Alternative if you don't have Glob: bdf_files = filter(f -> endswith(f, ".bdf"), readdir(data_dir, join=true))

println("Found $(length(bdf_files)) files to process")

# Store decoded results for each participant
all_decoded = eegfun.DecodedData[]

# Process each file
for (file_idx, data_file) in enumerate(bdf_files)
    participant_id = basename(data_file)
    println("\n[$(file_idx)/$(length(bdf_files))] Processing: $participant_id")
    
    try
        # Preprocessing (same for all files)
        dat = eegfun.read_bdf(data_file)
        dat = eegfun.create_eeg_dataframe(dat, layout_file)
        eegfun.rereference!(dat, :avg)
        eegfun.filter_data!(dat, "hp", 1)
        eegfun.is_extreme_value!(dat, 500)
        eegfun.mark_epoch_windows!(dat, [1, 2, 3, 4], [-0.5, 3.0])
        
        # Extract epochs
        epochs = eegfun.extract_epochs(dat, epoch_cfg, -0.5, 3)
        
        # Check we have both conditions
        if length(epochs) < 2
            println("  ⚠ Skipping: Only $(length(epochs)) condition(s) found")
            continue
        end
        
        # Decode
        epochs_for_decoding = epochs[1:2]  # Use first two conditions
        decoded = eegfun.decode(
            epochs_for_decoding,
            channels;
            model = model_method,
            time_range = time_range,
            n_iterations = n_iterations,
            n_folds = n_folds,
            equalize_trials = true,
        )
        
        push!(all_decoded, decoded)
        println("  ✓ Decoding complete: Max accuracy = $(round(maximum(decoded.average_score), digits=3))")
        
    catch e
        println("  ✗ Error processing $participant_id: $e")
        # Continue with next file
        continue
    end
end

# Create grand average
if isempty(all_decoded)
    error("No participants successfully decoded!")
end

# Convert to proper type (in case it's Vector{Any})
all_decoded = convert(Vector{eegfun.DecodedData}, all_decoded)

println("\n" * "="^70)
println("Creating grand average across $(length(all_decoded)) participants...")
grand_avg = eegfun.grand_average(all_decoded)

println("  ✓ Grand average created")
println("    Max accuracy: $(round(maximum(grand_avg.average_score), digits=3))")
println("    Time range: $(round(grand_avg.times[1], digits=2)) to $(round(grand_avg.times[end], digits=2)) s")

# Plot grand average
fig = eegfun.plot_decoding(
    grand_avg, 
    title = "Grand Average: $(epoch_cfg[1].name) vs $(epoch_cfg[2].name) ($(model_method))"
)
println("  ✓ Plot created")
