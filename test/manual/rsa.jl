#!/usr/bin/env julia
"""
Manual test for RSA (Representational Similarity Analysis) functionality

This script tests the RSA workflow:
1. Create test epoch data for multiple conditions
2. Compute RDMs for each participant separately
3. Compare to model RDMs
4. Average across participants (grand average)
5. Plot results

Tests both basic RSA and model comparison.
"""

using eegfun
using DataFrames
using Random

# Set seed for reproducibility
Random.seed!(42)

println("\n" * "="^70)
println("RSA (REPRESENTATIONAL SIMILARITY ANALYSIS) MANUAL TEST")
println("="^70)

# ============================================================================
# SETUP: Create Test Data
# ============================================================================

println("\n[1/5] Creating test epoch data...")

# Create test epoch data for multiple participants and conditions
# Each participant will have epochs for 3 conditions
n_participants = 3
n_epochs_per_condition = 30
channels = [:Fz, :Cz, :Pz, :POz]

# Helper function to create epoch data for a participant and condition
function create_participant_condition_epochs(
    participant_id::Int,
    condition_id::Int,
    condition_name::String,
    n_epochs::Int,
    base_freq::Float64,
)
    # Create layout
    layout = eegfun.Layout(
        DataFrame(
            label = channels,
            inc = [90.0, 0.0, -90.0, -90.0],
            azi = [0.0, 0.0, 0.0, 0.0],
        ),
        nothing,
        nothing,
    )

    # Create epochs with condition-specific patterns
    epochs = Vector{DataFrame}()
    sample_rate = 250
    time_window = (-0.2, 0.8)
    times = collect(time_window[1]:(1/sample_rate):time_window[2])

    for epoch_idx in 1:n_epochs
        epoch_data = DataFrame(time = times)

        # Add condition-specific signal patterns
        for (ch_idx, ch_name) in enumerate(channels)
            # Base signal with condition-specific frequency
            signal = sin.(2π * base_freq * times)
            
            # Add condition-specific phase shift
            phase_shift = condition_id * π / 4
            signal = sin.(2π * base_freq * times .+ phase_shift)
            
            # Add channel-specific amplitude
            amplitude = 1.0 + 0.3 * ch_idx
            
            # Add noise
            noise = 0.2 * randn(length(times))
            
            # Add participant-specific offset
            participant_offset = 0.1 * participant_id
            
            epoch_data[!, ch_name] = amplitude .* signal .+ noise .+ participant_offset
        end

        push!(epochs, epoch_data)
    end

    # Create EpochData object
    epoch_data_obj = eegfun.EpochData(
        "test_p$(participant_id)_$(condition_name).jld2",
        condition_id,
        condition_name,
        epochs,
        layout,
        sample_rate,
        eegfun.AnalysisInfo(),
    )

    return epoch_data_obj
end

# Create epochs for all participants and conditions
all_participant_epochs = Vector{Vector{eegfun.EpochData}}()

for p in 1:n_participants
    participant_epochs = [
        create_participant_condition_epochs(p, 1, "Face", n_epochs_per_condition, 2.0),
        create_participant_condition_epochs(p, 2, "Object", n_epochs_per_condition, 1.8),
        create_participant_condition_epochs(p, 3, "Scene", n_epochs_per_condition, 1.6),
    ]
    push!(all_participant_epochs, participant_epochs)
end

println("  ✓ Created epoch data for $n_participants participants")
println("  ✓ Each participant has $n_epochs_per_condition epochs per condition")
println("  ✓ Conditions: Face, Object, Scene")

# ============================================================================
# TEST 1: Basic RSA Analysis
# ============================================================================

println("\n[2/5] Test 1: Basic RSA analysis with correlation-based dissimilarity...")

# Perform RSA for each participant
all_rsa_results = Vector{eegfun.RsaData}()

for (p_idx, participant_epochs) in enumerate(all_participant_epochs)
    rsa_result = eegfun.rsa(
        participant_epochs,
        channels;
        time_range = (-0.2, 0.8),
        dissimilarity_measure = :correlation,
        average_trials = true,
    )
    push!(all_rsa_results, rsa_result)
    println("  ✓ Participant $p_idx RSA complete")
end

# Compute grand average
println("\n[3/5] Computing grand average RSA...")
grand_avg_rsa = eegfun.grand_average(all_rsa_results)
println("  ✓ Grand average created")
println("    Conditions: $(join(grand_avg_rsa.condition_names, ", "))")
println("    Time range: $(round(grand_avg_rsa.times[1], digits=2)) to $(round(grand_avg_rsa.times[end], digits=2)) s")

# Plot grand average RDM at a specific time point
println("\n[4/5] Plotting RSA results...")
fig1 = eegfun.plot_rdm(grand_avg_rsa, time_point = 0.3, title = "Grand Average RDM at 0.3s")
println("  ✓ RDM plot created")

# Plot average RDM
fig2 = eegfun.plot_rdm(grand_avg_rsa, title = "Grand Average RDM (averaged across time)")
println("  ✓ Average RDM plot created")

# ============================================================================
# TEST 2: Model Comparison
# ============================================================================

println("\n[5/5] Test 2: Model comparison...")

if isnothing(grand_avg_rsa)
    println("  ⚠ Skipping model comparison - grand average RSA not available")
else
    try
        # Create example model RDMs
        # Model 1: Face and Object are similar, Scene is different
        model1_rdm = [
            0.0 0.3 0.8;
            0.3 0.0 0.7;
            0.8 0.7 0.0
        ]

        # Model 2: All conditions equally dissimilar
        model2_rdm = [
            0.0 0.5 0.5;
            0.5 0.0 0.5;
            0.5 0.5 0.0
        ]

        # Compare grand average RSA to models
        grand_avg_with_models = eegfun.compare_models(
            grand_avg_rsa,
            [model1_rdm, model2_rdm];
            model_names = ["Semantic Model", "Equal Distance Model"],
            correlation_type = :spearman,
            n_permutations = 100,  # Reduced for testing
        )

        println("  ✓ Model comparison complete")
        println("    Models compared: $(join(grand_avg_with_models.model_names, ", "))")

        # Plot model correlations
        fig3 = eegfun.plot_model_correlations(
            grand_avg_with_models,
            title = "Model Correlations Over Time",
            colors = [:red, :blue],
        )
        println("  ✓ Model correlation plot created")

    catch e
        println("  ✗ Error: $e")
        rethrow(e)
    end
end

println("\n" * "="^70)
println("RSA TESTS COMPLETE")
println("="^70)

println("\nTests performed:")
println("  ✓ Basic RSA with correlation-based dissimilarity")
println("  ✓ Grand averaging across participants")
println("  ✓ Model comparison with permutation testing")
println("  ✓ RDM visualization")
println("  ✓ Model correlation visualization")

println("\nAll plots should be displayed above.")
println("Check that:")
println("  - RDMs show condition-specific patterns")
println("  - Model correlations vary over time")
println("  - Significant correlations are marked (if p-values computed)")

