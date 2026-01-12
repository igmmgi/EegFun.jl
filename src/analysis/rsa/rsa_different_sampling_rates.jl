#!/usr/bin/env julia
"""
Example: RSA with Different Sampling Rates

This example demonstrates how to handle temporal model data with different
sampling rates than neural data (e.g., EEG at 500 Hz vs eye tracking at 100 Hz).
"""

using eegfun
using DataFrames
using Random

Random.seed!(42)

println("="^70)
println("RSA WITH DIFFERENT SAMPLING RATES - EXAMPLE")
println("="^70)

# ============================================================================
# SETUP: Create example neural data (EEG at 500 Hz)
# ============================================================================

println("\n[1/4] Creating example neural data (EEG at 500 Hz)...")

channels = [:Fz, :Cz, :Pz]
layout = eegfun.Layout(
    DataFrame(label = channels, inc = [90.0, 0.0, -90.0], azi = [0.0, 0.0, 0.0]),
    nothing,
    nothing,
)

# Create 3 conditions: Face, Object, Scene
n_epochs = 20
neural_sample_rate = 500  # 500 Hz
neural_times = collect(-0.2:(1/neural_sample_rate):0.8)
n_neural_timepoints = length(neural_times)

epochs = eegfun.EpochData[]
for (cond_idx, cond_name) in enumerate(["Face", "Object", "Scene"])
    epoch_dfs = DataFrame[]
    for epoch_idx in 1:n_epochs
        df = DataFrame(time = neural_times)
        for ch in channels
            signal = sin.(2π * (0.5 + cond_idx * 0.1) * neural_times) .+ 0.1 * randn(length(neural_times))
            df[!, ch] = signal
        end
        push!(epoch_dfs, df)
    end
    push!(epochs, eegfun.EpochData(
        "test_data",
        cond_idx,
        cond_name,
        epoch_dfs,
        layout,
        neural_sample_rate,
        eegfun.AnalysisInfo(),
    ))
end

println("  ✓ Created EEG epoch data at $(neural_sample_rate) Hz")
println("    Time range: $(neural_times[1]) to $(neural_times[end]) s")
println("    Number of timepoints: $n_neural_timepoints")

# ============================================================================
# STEP 1: Compute Neural RDM
# ============================================================================

println("\n[2/4] Computing neural RDM from EEG data...")

neural_rsa = eegfun.rsa(
    epochs,
    channels;
    time_range = (-0.2, 0.8),
    dissimilarity_measure = :correlation,
)

println("  ✓ Neural RDM computed")
println("    Shape: [$(size(neural_rsa.rdm, 1)) timepoints × $(size(neural_rsa.rdm, 2)) conditions × $(size(neural_rsa.rdm, 3)) conditions]")

# ============================================================================
# STEP 2: Create Temporal Models with Different Sampling Rates
# ============================================================================

println("\n[3/4] Creating temporal model data with different sampling rates...")

# Example 1: Eye tracking at 100 Hz (different from EEG's 500 Hz)
println("\n  Temporal Model 1: Eye Tracking at 100 Hz")
eye_sample_rate = 100
eye_times = collect(-0.2:(1/eye_sample_rate):0.3)  # Different time range too
n_eye_timepoints = length(eye_times)

# Simulate eye tracking: 3 conditions × 2 features (x, y) × n_eye_timepoints
eye_tracking_data = zeros(Float64, 3, 2, n_eye_timepoints)
for cond_idx in 1:3
    for (t_idx, t) in enumerate(eye_times)
        eye_tracking_data[cond_idx, 1, t_idx] = sin(2π * 0.1 * t) + cond_idx * 0.1 + 0.1 * randn()
        eye_tracking_data[cond_idx, 2, t_idx] = cos(2π * 0.1 * t) + cond_idx * 0.1 + 0.1 * randn()
    end
end

println("    ✓ Created eye tracking data")
println("      Sampling rate: $(eye_sample_rate) Hz")
println("      Timepoints: $n_eye_timepoints")
println("      Shape: [$(size(eye_tracking_data, 1)) conditions × $(size(eye_tracking_data, 2)) features × $n_eye_timepoints timepoints]")

# Example 2: EDA at 50 Hz (different from both EEG and eye tracking)
println("\n  Temporal Model 2: EDA at 50 Hz")
eda_sample_rate = 50
eda_times = collect(-0.2:(1/eda_sample_rate):0.4)  # Different time range
n_eda_timepoints = length(eda_times)

# Simulate EDA: 3 conditions × 1 feature (conductance) × n_eda_timepoints
eda_data = zeros(Float64, 3, 1, n_eda_timepoints)
for cond_idx in 1:3
    for (t_idx, t) in enumerate(eda_times)
        eda_data[cond_idx, 1, t_idx] = exp(-t / 0.5) * cond_idx * 0.1 + 0.05 * randn()
    end
end

println("    ✓ Created EDA data")
println("      Sampling rate: $(eda_sample_rate) Hz")
println("      Timepoints: $n_eda_timepoints")
println("      Shape: [$(size(eda_data, 1)) conditions × $(size(eda_data, 2)) features × $n_eda_timepoints timepoints]")

# ============================================================================
# STEP 3: Create Temporal RDMs with Alignment
# ============================================================================

println("\n[4/4] Creating temporal RDMs with automatic resampling...")

# Option 1: Manual resampling
println("\n  Option 1: Manual resampling")
eye_rdms_manual = eegfun.create_temporal_rdm(
    eye_tracking_data,
    eye_times,
    align_to=neural_rsa.times,  # Resample to match neural data
    dissimilarity_measure=:euclidean,
    interpolation_method=:linear,
)
println("    ✓ Eye tracking RDM created (resampled from $n_eye_timepoints to $(size(eye_rdms_manual, 1)) timepoints)")

# Option 2: Using convenience function with tuples
println("\n  Option 2: Using convenience function with different sampling rates")
temporal_models = Dict(
    "Eye Tracking (100 Hz)" => (eye_tracking_data, eye_times),  # Tuple: (data, times)
    "EDA (50 Hz)" => (eda_data, eda_times),
)

temporal_rdms, temporal_names = eegfun.create_temporal_model_rdms(
    temporal_models,
    neural_rsa.times,  # Base times (used if model doesn't provide times)
    align_to=neural_rsa.times,  # Resample all to neural timepoints
    dissimilarity_measure=:euclidean,
    interpolation_method=:linear,
)

println("    ✓ Created $(length(temporal_names)) temporal RDMs")
for (idx, name) in enumerate(temporal_names)
    println("      - $name: [$(size(temporal_rdms[idx], 1)) × $(size(temporal_rdms[idx], 2)) × $(size(temporal_rdms[idx], 3))]")
end

# ============================================================================
# STEP 4: Compare All Models
# ============================================================================

println("\n[5/5] Comparing neural RDM to all models...")

# Combine all models (they're all aligned to neural timepoints now)
all_model_rdms = [eye_rdms_manual, temporal_rdms...]
all_model_names = ["Eye Tracking (Manual)", temporal_names...]

rsa_with_models = eegfun.compare_models(
    neural_rsa,
    all_model_rdms,
    model_names=all_model_names,
    correlation_type=:spearman,
    n_permutations=100,
)

println("  ✓ Comparison complete")

# ============================================================================
# SUMMARY
# ============================================================================

println("\n" * "="^70)
println("SUMMARY")
println("="^70)

println("\nNeural Data (EEG):")
println("  - Sampling rate: $(neural_sample_rate) Hz")
println("  - Timepoints: $n_neural_timepoints")
println("  - Time range: $(neural_times[1]) to $(neural_times[end]) s")

println("\nModel Data:")
println("  1. Eye Tracking")
println("     - Original sampling rate: $(eye_sample_rate) Hz")
println("     - Original timepoints: $n_eye_timepoints")
println("     - Original time range: $(eye_times[1]) to $(eye_times[end]) s")
println("     - Resampled to: $(size(eye_rdms_manual, 1)) timepoints (neural rate)")
println("  2. EDA")
println("     - Original sampling rate: $(eda_sample_rate) Hz")
println("     - Original timepoints: $n_eda_timepoints")
println("     - Original time range: $(eda_times[1]) to $(eda_times[end]) s")
println("     - Resampled to: $(size(temporal_rdms[2], 1)) timepoints (neural rate)")

println("\nKey Points:")
println("  ✓ Models with different sampling rates are automatically resampled")
println("  ✓ Models with different time ranges are handled (extrapolation)")
println("  ✓ Linear interpolation used by default (can use :nearest or :cubic)")
println("  ✓ All models aligned to neural data timepoints before comparison")
println("  ✓ RDMs computed after resampling, ensuring proper temporal alignment")

println("\n" * "="^70)

