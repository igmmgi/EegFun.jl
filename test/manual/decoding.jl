using EegFun
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
layout = EegFun.Layout(DataFrame(
    label = channels,  # [:Fz, :Cz, :Pz, :POz]
    inc = [90.0, 0.0, -90.0, -90.0],  # Rough positions
    azi = [0.0, 0.0, 0.0, 0.0],
), nothing, nothing)

analysis_info = EegFun.AnalysisInfo(:none, 0.0, 0.0)

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
                signal =
                    signal_strength * exp.(-((time .- 0.2) .^ 2) / (2 * 0.05^2)) .+
                    0.3 * signal_strength * exp.(-((time .- 0.4) .^ 2) / (2 * 0.08^2))
            else
                # Condition 2: Negative ERP component
                signal =
                    -signal_strength * exp.(-((time .- 0.2) .^ 2) / (2 * 0.05^2)) .+
                    0.2 * signal_strength * exp.(-((time .- 0.5) .^ 2) / (2 * 0.1^2))
            end

            # Add channel-specific variation and noise
            channel_factor = 1.0 + (ch_idx - 1) * 0.1
            noise = 0.3 * randn(n_timepoints)  # SNR ~ 3:1
            df[!, ch] = signal .* channel_factor .+ noise
        end

        push!(epochs, df)
    end

    return EegFun.EpochData("participant_$(participant_id)", condition_id, condition_name, epochs, layout, sample_rate, analysis_info)
end

# Create epochs for all participants
all_participant_epochs = []
for p = 1:n_participants
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
epochs_p1 = all_participant_epochs[1]
epochs_p2 = all_participant_epochs[2]
epochs_p3 = all_participant_epochs[3]


decoded_p1 = EegFun.decode_libsvm(
    epochs_p1;
    n_iterations = 10,  # Reduced for testing
    n_folds = 3,
    equalize_trials = true,
)
# Decode for Participant 2 - default model used automatically
decoded_p2 = EegFun.decode_libsvm(epochs_p2, n_iterations = 10, n_folds = 3, equalize_trials = true)
# Decode for Participant 3 - default model used automatically
decoded_p3 = EegFun.decode_libsvm(epochs_p3, n_iterations = 10, n_folds = 3, equalize_trials = true)
# Grand average
all_decoded = [decoded_p1, decoded_p2, decoded_p3]
grand_avg_decoded = EegFun.grand_average(all_decoded)

# fig = EegFun.plot_decoding(all_decoded, title = "Grand Average: Face vs Object ($(model_method))")
fig = EegFun.plot_decoding(grand_avg_decoded, title = "Test Decoding")


# ============================================================================
# STATISTICAL TESTING
# ============================================================================
# Quick statistical test on synthetic data
println("\n[Quick Test] Statistical testing on synthetic data...")
stats_synthetic = EegFun.test_against_chance(
    all_decoded,
    alpha = 0.05,
    tail = :right,
    # correction_method = :bonferroni,  # No correction for quick test
    correction_method = :none,  # No correction for quick test
)
println("  ✓ Found $(sum(stats_synthetic.significant_mask)) significant time points")

# Plot with significance
fig_stats_synthetic =
    EegFun.plot_decoding(grand_avg_decoded, stats_synthetic, title = "Synthetic Data with Significance", show_significance = true)


println("\n[Quick Test] Statistical testing on synthetic data...")
stats_synthetic = EegFun.test_against_chance_cluster(all_decoded, alpha = 0.2)
println("  ✓ Found $(sum(stats_synthetic.significant_mask)) significant time points")

# Plot with significance
fig_stats_synthetic =
    EegFun.plot_decoding(grand_avg_decoded, stats_synthetic, title = "Synthetic Data with Significance", show_significance = true)



using EegFun
using Glob  # For finding files (or use readdir with filter)

# Configuration
data_dir = joinpath(@__DIR__, "..", "..", "..", "AttentionExp", "recoded")
layout_file = EegFun.read_layout("./data/layouts/biosemi/biosemi72.csv")
EegFun.polar_to_cartesian_xy!(layout_file)

# Epoch configuration
epoch_cfg = [
    EegFun.EpochCondition(name = "ExampleEpoch1", trigger_sequences = [[1], [3]]),
    EegFun.EpochCondition(name = "ExampleEpoch2", trigger_sequences = [[2], [4]]),
]

# Find all .bdf files in directory
bdf_files = glob("*.bdf", data_dir)
println("Found $(length(bdf_files)) files to process")

# Process each file
all_decoded = EegFun.DecodedData[]
for (file_idx, data_file) in enumerate(bdf_files)
    participant_id = basename(data_file)
    println("\n[$(file_idx)/$(length(bdf_files))] Processing: $participant_id")

    try
        # Preprocessing (same for all files)
        dat = EegFun.read_raw_data(data_file)
        dat = EegFun.create_eeg_dataframe(dat, layout_file)
        EegFun.rereference!(dat, :avg)
        EegFun.highpass_filter!(dat, 1)

        # Extract epochs
        epochs = EegFun.extract_epochs(dat, epoch_cfg, -0.2, 2.5)

        # Decode
        decoded = EegFun.decode_libsvm(
            epochs,
            channel_selection = EegFun.channels([:PO7, :PO8]),
            n_iterations = 10,
            n_folds = 3,
            equalize_trials = true,
        )

        push!(all_decoded, decoded)
        println("  ✓ Decoding complete: Max accuracy = $(round(maximum(decoded.average_score), digits=3))")

    catch e
        println("  ✗ Error processing $participant_id: $e")
        continue
    end
end

grand_avg_decoded = EegFun.grand_average(all_decoded)
println("  ✓ Grand average created")
println("    Max accuracy: $(round(maximum(grand_avg_decoded.average_score), digits=3))")
println("    Time range: $(round(grand_avg_decoded.times[1], digits=2)) to $(round(grand_avg_decoded.times[end], digits=2)) s")

# Plot grand average
fig = EegFun.plot_decoding(grand_avg_decoded, title = "Grand Average: Test Decoding")
println("  ✓ Plot created")


# ============================================================================
# STATISTICAL TESTING
# ============================================================================

println("\n" * "="^70)
println("Statistical Testing")
println("="^70)

# Test 1: One-sample t-test against chance (with Bonferroni correction)
println("\n[1/2] One-sample t-test against chance (Bonferroni correction)...")
stats_bonferroni = EegFun.test_against_chance(
    all_decoded,
    alpha = 0.05,
    tail = :right,  # Test if accuracy > chance
    correction_method = :none,
    # correction_method = :bonferroni,
)
println("  ✓ T-test complete")
println("    Significant time points: $(sum(stats_bonferroni.significant_mask)) / $(length(stats_bonferroni.significant_mask))")
println("    Min p-value: $(round(minimum(stats_bonferroni.p_values), digits=4))")
println("    Max t-value: $(round(maximum(stats_bonferroni.t_statistics), digits=3))")

# Plot with significance markers (Bonferroni)
fig_stats_bonf = EegFun.plot_decoding(
    grand_avg_decoded,
    stats_bonferroni,
    title = "Grand Average with Significance (Bonferroni)",
    show_significance = true,
    sig_color = :yellow,
    sig_alpha = 0.3,
)
println("  ✓ Plot with significance markers created")

# Test 2: Cluster-based permutation test
println("\n[2/2] Cluster-based permutation test...")
stats_cluster = EegFun.test_against_chance_cluster(
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
            println(
                "      Cluster $i: $(round(cluster.time_range[1], digits=3))-$(round(cluster.time_range[2], digits=3)) s, p=$(round(cluster.p_value, digits=4))",
            )
        end
    end
else
    println("    No clusters found")
end
println("    Significant time points: $(sum(stats_cluster.significant_mask)) / $(length(stats_cluster.significant_mask))")

# Plot with cluster-based significance markers
fig_stats_cluster = EegFun.plot_decoding(
    grand_avg_decoded,
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

