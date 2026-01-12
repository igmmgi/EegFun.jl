#!/usr/bin/env julia
"""
Example: RSA with Static vs Temporal Models

This example demonstrates how to handle both:
- Static models (single timepoint): e.g., reaction times, semantic embeddings
- Temporal models (multiple timepoints): e.g., eye tracking, EDA

Both can be compared to temporal neural RDMs (from EEG).
"""

using eegfun
using DataFrames
using Random

Random.seed!(42)

println("="^70)
println("RSA WITH STATIC VS TEMPORAL MODELS - EXAMPLE")
println("="^70)

# ============================================================================
# SETUP: Create example neural data (temporal)
# ============================================================================

println("\n[1/5] Creating example neural data (EEG - temporal)...")

channels = [:Fz, :Cz, :Pz]
layout = eegfun.Layout(
    DataFrame(label = channels, inc = [90.0, 0.0, -90.0], azi = [0.0, 0.0, 0.0]),
    nothing,
    nothing,
)

# Create 3 conditions: Face, Object, Scene
n_epochs = 20
sample_rate = 250
times = collect(-0.2:(1/sample_rate):0.8)
n_timepoints = length(times)

epochs = eegfun.EpochData[]
for (cond_idx, cond_name) in enumerate(["Face", "Object", "Scene"])
    epoch_dfs = DataFrame[]
    for epoch_idx in 1:n_epochs
        df = DataFrame(time = times)
        for ch in channels
            signal = sin.(2π * (0.5 + cond_idx * 0.1) * times) .+ 0.1 * randn(length(times))
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
        sample_rate,
        eegfun.AnalysisInfo(),
    ))
end

println("  ✓ Created EEG epoch data for 3 conditions")
println("    Time range: $(times[1]) to $(times[end]) s ($n_timepoints timepoints)")

# ============================================================================
# STEP 1: Compute Neural RDM (temporal)
# ============================================================================

println("\n[2/5] Computing neural RDM from EEG data (temporal)...")

neural_rsa = eegfun.rsa(
    epochs,
    channels;
    time_range = (-0.2, 0.8),
    dissimilarity_measure = :correlation,
)

println("  ✓ Neural RDM computed")
println("    Shape: [$(size(neural_rsa.rdm, 1)) timepoints × $(size(neural_rsa.rdm, 2)) conditions × $(size(neural_rsa.rdm, 3)) conditions]")

# ============================================================================
# STEP 2: Create Static Models (single timepoint)
# ============================================================================

println("\n[3/5] Creating static model RDMs (single timepoint)...")

# Example 1: Reaction times (static - single value per condition)
println("\n  Static Model 1: Reaction Times")
rts = [0.35, 0.42, 0.38]  # RTs for Face, Object, Scene
rt_rdm = eegfun.create_rdm_from_reaction_times(rts)
println("    ✓ Created static RDM from reaction times")
println("      Shape: [$(size(rt_rdm, 1)) × $(size(rt_rdm, 2))] (static)")

# Example 2: Semantic embeddings (static - single vector per condition)
println("\n  Static Model 2: Semantic Embeddings")
word_embeddings = [
    [0.1, 0.2, 0.3, 0.4],  # Face
    [0.2, 0.3, 0.4, 0.5],  # Object
    [0.5, 0.6, 0.7, 0.8],  # Scene
]
semantic_rdm = eegfun.create_rdm_from_vectors(word_embeddings, dissimilarity_measure=:euclidean)
println("    ✓ Created static RDM from semantic embeddings")
println("      Shape: [$(size(semantic_rdm, 1)) × $(size(semantic_rdm, 2))] (static)")

# ============================================================================
# STEP 3: Create Temporal Models (multiple timepoints)
# ============================================================================

println("\n[4/5] Creating temporal model RDMs (multiple timepoints)...")

# Example 1: Eye tracking data (temporal - same timepoints as EEG)
println("\n  Temporal Model 1: Eye Tracking")
# Simulate eye tracking: 3 conditions × 2 features (x, y) × n_timepoints
eye_tracking_data = zeros(Float64, 3, 2, n_timepoints)
for cond_idx in 1:3
    for t in 1:n_timepoints
        # Simulate eye position over time
        eye_tracking_data[cond_idx, 1, t] = sin(2π * 0.1 * times[t]) + cond_idx * 0.1 + 0.1 * randn()
        eye_tracking_data[cond_idx, 2, t] = cos(2π * 0.1 * times[t]) + cond_idx * 0.1 + 0.1 * randn()
    end
end
eye_rdms = eegfun.create_temporal_rdm(eye_tracking_data, times, dissimilarity_measure=:euclidean)
println("    ✓ Created temporal RDM from eye tracking")
println("      Shape: [$(size(eye_rdms, 1)) timepoints × $(size(eye_rdms, 2)) conditions × $(size(eye_rdms, 3)) conditions]")

# Example 2: EDA (Electrodermal Activity) data (temporal)
println("\n  Temporal Model 2: EDA (Electrodermal Activity)")
# Simulate EDA: 3 conditions × 1 feature (conductance) × n_timepoints
eda_data = zeros(Float64, 3, 1, n_timepoints)
for cond_idx in 1:3
    for t in 1:n_timepoints
        # Simulate EDA response over time
        eda_data[cond_idx, 1, t] = exp(-times[t] / 0.5) * cond_idx * 0.1 + 0.05 * randn()
    end
end
eda_rdms = eegfun.create_temporal_rdm(eda_data, times, dissimilarity_measure=:euclidean)
println("    ✓ Created temporal RDM from EDA")
println("      Shape: [$(size(eda_rdms, 1)) timepoints × $(size(eda_rdms, 2)) conditions × $(size(eda_rdms, 3)) conditions]")

# ============================================================================
# STEP 4: Compare Neural RDM to All Models (Static + Temporal)
# ============================================================================

println("\n[5/5] Comparing neural RDM to all models (static + temporal)...")

# Combine static and temporal models
# Note: compare_models automatically handles both types
all_model_rdms = [
    rt_rdm,          # Static: Matrix{Float64}
    semantic_rdm,    # Static: Matrix{Float64}
    eye_rdms,        # Temporal: Array{Float64, 3}
    eda_rdms,        # Temporal: Array{Float64, 3}
]

all_model_names = [
    "Reaction Times (Static)",
    "Semantic Embeddings (Static)",
    "Eye Tracking (Temporal)",
    "EDA (Temporal)",
]

println("\n  Comparing to $(length(all_model_names)) models:")
println("    - 2 static models (same RDM at all timepoints)")
println("    - 2 temporal models (RDM computed at each timepoint)")

rsa_with_all_models = eegfun.compare_models(
    neural_rsa,
    all_model_rdms,
    model_names=all_model_names,
    correlation_type=:spearman,
    n_permutations=100,  # Reduced for example
)

println("  ✓ Comparison complete")

# ============================================================================
# STEP 5: Visualize Results
# ============================================================================

println("\n[6/6] Visualizing results...")

# Plot model correlations
fig = eegfun.plot_model_correlations(
    rsa_with_all_models,
    title="Neural RDM vs Static & Temporal Models",
    colors=[:blue, :red, :green, :orange],
)
println("  ✓ Model correlation plot created")

# ============================================================================
# SUMMARY
# ============================================================================

println("\n" * "="^70)
println("SUMMARY")
println("="^70)

println("\nNeural Data (EEG):")
println("  - Temporal: RDMs computed at each timepoint")
println("  - Shape: [$(size(neural_rsa.rdm, 1)) timepoints × $(size(neural_rsa.rdm, 2)) × $(size(neural_rsa.rdm, 3))]")

println("\nStatic Models (single timepoint):")
println("  1. Reaction Times")
println("     - Single value per condition")
println("     - RDM shape: [$(size(rt_rdm, 1)) × $(size(rt_rdm, 2))]")
println("     - Same RDM used at all timepoints")
println("  2. Semantic Embeddings")
println("     - Single vector per condition")
println("     - RDM shape: [$(size(semantic_rdm, 1)) × $(size(semantic_rdm, 2))]")
println("     - Same RDM used at all timepoints")

println("\nTemporal Models (multiple timepoints):")
println("  1. Eye Tracking")
println("     - Multiple features (x, y) measured over time")
println("     - RDM shape: [$(size(eye_rdms, 1)) timepoints × $(size(eye_rdms, 2)) × $(size(eye_rdms, 3))]")
println("     - Different RDM at each timepoint")
println("  2. EDA (Electrodermal Activity)")
println("     - Single feature (conductance) measured over time")
println("     - RDM shape: [$(size(eda_rdms, 1)) timepoints × $(size(eda_rdms, 2)) × $(size(eda_rdms, 3))]")
println("     - Different RDM at each timepoint")

println("\nComparison:")
println("  - Static models: Compared to neural RDM at each timepoint")
println("    (same model RDM, different neural RDM)")
println("  - Temporal models: Compared timepoint-by-timepoint")
println("    (both model and neural RDMs change over time)")

println("\nKey Points:")
println("  ✓ Static models are replicated across all timepoints")
println("  ✓ Temporal models must have same number of timepoints as neural data")
println("  ✓ Both types can be compared in a single call to compare_models()")
println("  ✓ Results show how well each model matches neural patterns over time")

println("\n" * "="^70)

