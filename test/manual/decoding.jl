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
using BenchmarkTools

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
    cond1 = create_participant_condition_epochs(p, 1, "Face", n_epochs_per_condition, 0.1)
    # Participant p, Condition 2
    cond2 = create_participant_condition_epochs(p, 2, "Object", n_epochs_per_condition, 0.1)
    push!(all_participant_epochs, [cond1, cond2])
end

println("  ✓ Created epoch data for $n_participants participants")
println("  ✓ Each participant has $(n_epochs_per_condition) epochs per condition")
println("  ✓ Conditions: Face vs Object")

# ==========================================
# TEST 1: Binary Classification with Default 
# ==========================================

# model_method = :logistic
model_method = :svm
# model_method = :lda
# Decode for Participant 1 - default model is used automatically!
epochs_p1 = all_participant_epochs[1]

@btime decoded_p1 = eegfun.decode(
    epochs_p1;
    model = model_method,
    n_iterations = 100,  # Reduced for testing
    n_folds = 3,
    equalize_trials = true,
)

# Decode for Participant 2 - default model used automatically
epochs_p2 = all_participant_epochs[2]
decoded_p2 = eegfun.decode(
    epochs_p2,
    model = model_method,
    n_iterations = 100,
    n_folds = 3,
    equalize_trials = true,
)
# Decode for Participant 3 - default model used automatically
epochs_p3 = all_participant_epochs[3]
decoded_p3 = eegfun.decode(
    epochs_p3,
    model = model_method,
    n_iterations = 100,
    n_folds = 3,
    equalize_trials = true,
)
# Grand average
all_decoded = [decoded_p1, decoded_p2, decoded_p3]
grand_avg_mlj = eegfun.grand_average(all_decoded)
# Plot grand average
fig1 = eegfun.plot_decoding(grand_avg_mlj, title = "Grand Average: Face vs Object ($(model_method))")
# fig1 = eegfun.plot_decoding(all_decoded, title = "Grand Average: Face vs Object ($(model_method))")





decoded_p1 = eegfun.decode_pegasos(
    epochs_p1;
    n_iterations = 100,  # Reduced for testing
    n_folds = 3,
    equalize_trials = true,
)
# Decode for Participant 2 - default model used automatically
decoded_p2 = eegfun.decode_pegasos(
    epochs_p2,
    n_iterations = 100,
    n_folds = 3,
    equalize_trials = true,
)
# Decode for Participant 3 - default model used automatically
decoded_p3 = eegfun.decode_pegasos(
    epochs_p3,
    n_iterations = 100,
    n_folds = 3,
    equalize_trials = true,
)
# Grand average
all_decoded = [decoded_p1, decoded_p2, decoded_p3]
grand_avg_pegasos = eegfun.grand_average(all_decoded)
# Plot grand average
fig2 = eegfun.plot_decoding(grand_avg_pegasos, title = "Grand Average: Face vs Object (Pegasos)")
# fig1 = eegfun.plot_decoding(all_decoded, title = "Grand Average: Face vs Object ($(model_method))")




@btime decoded_p1 = eegfun.decode_libsvm(
    epochs_p1;
    n_iterations = 10,  # Reduced for testing
    n_folds = 3,
    equalize_trials = true,
)
# Decode for Participant 2 - default model used automatically
decoded_p2 = eegfun.decode_libsvm(
    epochs_p2,
    n_iterations = 100,
    n_folds = 3,
    equalize_trials = true,
)
# Decode for Participant 3 - default model used automatically
decoded_p3 = eegfun.decode_libsvm(
    epochs_p3,
    n_iterations = 100,
    n_folds = 3,
    equalize_trials = true,
)
# Grand average
all_decoded = [decoded_p1, decoded_p2, decoded_p3]
grand_avg_pegasos = eegfun.grand_average(all_decoded)
# Plot grand average
fig2 = eegfun.plot_decoding(grand_avg_pegasos, title = "Grand Average: Face vs Object (Pegasos)")
# fig1 = eegfun.plot_decoding(all_decoded, title = "Grand Average: Face vs Object ($(model_method))")













# Quick statistical test on synthetic data
println("\n[Quick Test] Statistical testing on synthetic data...")
stats_synthetic = eegfun.test_against_chance(
    all_decoded,
    alpha = 0.05,
    tail = :right,
    correction_method = :bonferroni,  # No correction for quick test
)
println("  ✓ Found $(sum(stats_synthetic.significant_mask)) significant time points")

# Plot with significance
fig_stats_synthetic = eegfun.plot_decoding(
    grand_avg,
    stats_synthetic,
    title = "Synthetic Data with Significance",
    show_significance = true,
)


println("\n[Quick Test] Statistical testing on synthetic data...")
stats_synthetic = eegfun.test_against_chance_cluster(
    all_decoded,
    alpha = 0.2,
    tail = :right,
)
println("  ✓ Found $(sum(stats_synthetic.significant_mask)) significant time points")

# Plot with significance
fig_stats_synthetic = eegfun.plot_decoding(
    grand_avg,
    stats_synthetic,
    title = "Synthetic Data with Significance",
    show_significance = true,
)



# ============================================================================
# STATISTICAL TESTING
# ============================================================================

println("\n" * "="^70)
println("Statistical Testing")
println("="^70)

# Test 1: One-sample t-test against chance (with Bonferroni correction)
println("\n[1/2] One-sample t-test against chance (Bonferroni correction)...")
stats_bonferroni = eegfun.test_against_chance(
    all_decoded,
    alpha = 0.05,
    tail = :right,  # Test if accuracy > chance
    correction_method = :bonferroni,
)
println("  ✓ T-test complete")
println("    Significant time points: $(sum(stats_bonferroni.significant_mask)) / $(length(stats_bonferroni.significant_mask))")
println("    Min p-value: $(round(minimum(stats_bonferroni.p_values), digits=4))")
println("    Max t-value: $(round(maximum(stats_bonferroni.t_statistics), digits=3))")

# Plot with significance markers (Bonferroni)
fig_stats_bonf = eegfun.plot_decoding(
    grand_avg,
    stats_bonferroni,
    title = "Grand Average with Significance (Bonferroni)",
    show_significance = true,
    sig_color = :yellow,
    sig_alpha = 0.3,
)
println("  ✓ Plot with significance markers created")

# Test 2: Cluster-based permutation test
println("\n[2/2] Cluster-based permutation test...")
stats_cluster = eegfun.test_against_chance_cluster(
    all_decoded,
    alpha = 0.05,
    tail = :right,
    n_permutations = 1000,  # Use more for publication (e.g., 10000)
    cluster_statistic = :sum,  # or :max
    random_seed = 42,  # For reproducibility
    show_progress = true,
)
println("  ✓ Cluster permutation test complete")
if !isnothing(stats_cluster.clusters) && !isempty(stats_cluster.clusters)
    n_sig_clusters = sum(c.is_significant for c in stats_cluster.clusters)
    println("    Significant clusters: $n_sig_clusters / $(length(stats_cluster.clusters))")
    for (i, cluster) in enumerate(stats_cluster.clusters)
        if cluster.is_significant
            println("      Cluster $i: $(round(cluster.time_range[1], digits=3))-$(round(cluster.time_range[2], digits=3)) s, p=$(round(cluster.p_value, digits=4))")
        end
    end
else
    println("    No clusters found")
end
println("    Significant time points: $(sum(stats_cluster.significant_mask)) / $(length(stats_cluster.significant_mask))")

# Plot with cluster-based significance markers
fig_stats_cluster = eegfun.plot_decoding(
    grand_avg,
    stats_cluster,
    title = "Grand Average with Cluster-Based Significance",
    show_significance = true,
    sig_color = :green,
    sig_alpha = 0.3,
)
println("  ✓ Plot with cluster significance markers created")

println("\n" * "="^70)
println("All tests complete!")
println("="^70)















































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
n_iterations = 10
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
        # eegfun.is_extreme_value!(dat, 500)
        # eegfun.mark_epoch_windows!(dat, [1, 2, 3, 4], [-0.5, 3.0])
        
        # Extract epochs
        epochs = eegfun.extract_epochs(dat, epoch_cfg, -0.2, 2.5)
        
        # Check we have both conditions
        if length(epochs) < 2
            println("  ⚠ Skipping: Only $(length(epochs)) condition(s) found")
            continue
        end
        
        # Decode
        epochs_for_decoding = epochs[1:2]  # Use first two conditions
        decoded = eegfun.decode(
            epochs_for_decoding,
            model = model_method,
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


# ============================================================================
# STATISTICAL TESTING
# ============================================================================

println("\n" * "="^70)
println("Statistical Testing")
println("="^70)

# Test 1: One-sample t-test against chance (with Bonferroni correction)
println("\n[1/2] One-sample t-test against chance (Bonferroni correction)...")
stats_bonferroni = eegfun.test_against_chance(
    all_decoded,
    alpha = 0.05,
    tail = :right,  # Test if accuracy > chance
    correction_method = :bonferroni,
)
println("  ✓ T-test complete")
println("    Significant time points: $(sum(stats_bonferroni.significant_mask)) / $(length(stats_bonferroni.significant_mask))")
println("    Min p-value: $(round(minimum(stats_bonferroni.p_values), digits=4))")
println("    Max t-value: $(round(maximum(stats_bonferroni.t_statistics), digits=3))")

# Plot with significance markers (Bonferroni)
fig_stats_bonf = eegfun.plot_decoding(
    grand_avg,
    stats_bonferroni,
    title = "Grand Average with Significance (Bonferroni)",
    show_significance = true,
    sig_color = :yellow,
    sig_alpha = 0.3,
)
println("  ✓ Plot with significance markers created")

# Test 2: Cluster-based permutation test
println("\n[2/2] Cluster-based permutation test...")
stats_cluster = eegfun.test_against_chance_cluster(
    all_decoded,
    alpha = 0.05,
    tail = :right,
    n_permutations = 1000,  # Use more for publication (e.g., 10000)
    cluster_statistic = :sum,  # or :max
    random_seed = 42,  # For reproducibility
    show_progress = true,
)
println("  ✓ Cluster permutation test complete")
if !isnothing(stats_cluster.clusters) && !isempty(stats_cluster.clusters)
    n_sig_clusters = sum(c.is_significant for c in stats_cluster.clusters)
    println("    Significant clusters: $n_sig_clusters / $(length(stats_cluster.clusters))")
    for (i, cluster) in enumerate(stats_cluster.clusters)
        if cluster.is_significant
            println("      Cluster $i: $(round(cluster.time_range[1], digits=3))-$(round(cluster.time_range[2], digits=3)) s, p=$(round(cluster.p_value, digits=4))")
        end
    end
else
    println("    No clusters found")
end
println("    Significant time points: $(sum(stats_cluster.significant_mask)) / $(length(stats_cluster.significant_mask))")

# Plot with cluster-based significance markers
fig_stats_cluster = eegfun.plot_decoding(
    grand_avg,
    stats_cluster,
    title = "Grand Average with Cluster-Based Significance",
    show_significance = true,
    sig_color = :green,
    sig_alpha = 0.3,
)
println("  ✓ Plot with cluster significance markers created")

println("\n" * "="^70)
println("All tests complete!")
println("="^70)

