# Demo: Multivariate Pattern Analysis (MVPA) / Decoding
# Shows time-resolved classification using LIBSVM, with synthetic data creation,
# grand averaging, and statistical testing (uncorrected, Bonferroni, cluster-based).

using EegFun
using DataFrames
using Random

@info EegFun.section("MVPA/DECODING MANUAL TEST")

# Create synthetic data for 10 participants
# Adjustable difficulty via signal_strength and noise_level
difficulty = "hard"

if difficulty == "easy" # maybe a bit extreme! :-)
    signal_strength, noise_level = 2.0, 0.1
elseif difficulty == "medium"
    signal_strength, noise_level = 1.0, 0.3
else  # hard
    signal_strength, noise_level = 0.25, 0.75
end

println("  Using difficulty: $difficulty (signal=$signal_strength, noise=$noise_level)")
all_synthetic = [
    [
        EegFun.create_synthetic_epochs(p, 1, "Cond1", 100; signal_strength = signal_strength, noise_level = noise_level),
        EegFun.create_synthetic_epochs(p, 2, "Cond2", 100; signal_strength = signal_strength, noise_level = noise_level),
    ] for p = 1:10
]
EegFun.plot_epochs(all_synthetic[1]) # VP1
EegFun.plot_epochs(all_synthetic[2]) # VP2
EegFun.plot_epochs(all_synthetic[3]) # VP3
# and so on

# Decode synthetic data (batch method)
decoded_synthetic = EegFun.decode_libsvm(all_synthetic; n_iterations = 20, n_folds = 3)
grand_avg_synthetic = EegFun.grand_average(decoded_synthetic)

EegFun.plot_decoding(decoded_synthetic)    # every "VP"
EegFun.plot_decoding(grand_avg_synthetic)  # grand average

# Test and plot with different methods
stats_none = EegFun.test_against_chance(decoded_synthetic, alpha = 0.05, correction_method = :none)
EegFun.plot_decoding(grand_avg_synthetic, stats_none, title = "Synthetic Data: No Correction")

stats_bonf = EegFun.test_against_chance(decoded_synthetic, alpha = 0.05, correction_method = :bonferroni)
EegFun.plot_decoding(grand_avg_synthetic, stats_bonf, title = "Synthetic Data: Bonferroni")

stats_cluster = EegFun.test_against_chance_cluster(decoded_synthetic, alpha = 0.05)
EegFun.plot_decoding(grand_avg_synthetic, stats_cluster, title = "Synthetic Data: Cluster-based")

# ============================================================================
# OPTION 2: REAL DATA TEST
# ============================================================================
@info EegFun.section("Real Data Test")

# Prepare decoding data using prepare_decoding (like prepare_stats for statistics)
println("Preparing data...")
participant_epochs = EegFun.prepare_decoding(
    "epochs_good",
    input_dir = "./resources/data/julia/erps",
    participant_selection = EegFun.participants(),
    condition_selection = EegFun.conditions([1, 2]),  # Compare conditions 1 and 2
    channel_selection = EegFun.channels(),            # All channels
    interval_selection = EegFun.times((-0.2, 1.5)),   # Time interval
)

# Decode all participants (batch method)
all_decoded = EegFun.decode_libsvm(participant_epochs; n_iterations = 5, n_folds = 3, equalize_trials = true)

# Grand average
grand_avg_decoded = EegFun.grand_average(all_decoded)

# Test and plot with different methods
stats_none = EegFun.test_against_chance(all_decoded, alpha = 0.05, correction_method = :none)
EegFun.plot_decoding(grand_avg_decoded, stats_none, title = "Real Data: No Correction")

stats_bonf = EegFun.test_against_chance(all_decoded, alpha = 0.05, correction_method = :bonferroni)
EegFun.plot_decoding(grand_avg_decoded, stats_bonf, title = "Real Data: Bonferroni")

stats_cluster = EegFun.test_against_chance_cluster(all_decoded, alpha = 0.05)
EegFun.plot_decoding(grand_avg_decoded, stats_cluster, title = "Real Data: Cluster-based")
