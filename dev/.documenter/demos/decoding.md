
# MVPA Decoding {#MVPA-Decoding}

Multivariate pattern analysis (MVPA) decoding workflow with SVM classification and statistical testing.

## Overview {#Overview}

Demonstrates Multivariate pattern analysis (MVPA) decoding workflow with SVM classification and statistical testing.

## Source Code {#Source-Code}

```julia
using EegFun
using DataFrames
using Random

@info EegFun.section("MVPA/DECODING MANUAL TEST")

# Simplified synthetic epoch creation for fast testing
# Adjustable difficulty via signal_strength and noise_level
function create_synthetic_epochs(
    participant_id,
    condition_id,
    condition_name,
    n_epochs;
    n_timepoints = 500,
    channels = [:Fz, :Cz, :Pz],
    signal_strength = 1.0,
    noise_level = 0.3,
)
    time = range(-0.2, 0.8, length = n_timepoints)
    layout = EegFun.Layout(DataFrame(label = channels, inc = [90.0, 0.0, -90.0], azi = [0.0, 0.0, 0.0]), nothing, nothing)

    epochs = DataFrame[]
    for epoch = 1:n_epochs
        df = DataFrame(time = collect(time), epoch = fill(epoch, n_timepoints))
        for ch in channels
            # Condition-specific signal: positive for cond1, negative for cond2
            signal = (condition_id == 1 ? signal_strength : -signal_strength) * exp.(-((time .- 0.2) .^ 2) / (2 * 0.05^2))
            df[!, ch] = signal .+ noise_level * randn(n_timepoints)
        end
        push!(epochs, df)
    end

    return EegFun.EpochData(
        "participant_$(participant_id)",
        condition_id,
        condition_name,
        epochs,
        layout,
        200,
        EegFun.AnalysisInfo(:none, 0.0, 0.0),
    )
end

# Create synthetic data for 3 participants
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
        create_synthetic_epochs(p, 1, "Cond1", 100; signal_strength = signal_strength, noise_level = noise_level),
        create_synthetic_epochs(p, 2, "Cond2", 100; signal_strength = signal_strength, noise_level = noise_level),
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
    input_dir = "/home/ian/Documents/Julia/output_data",
    participant_selection = EegFun.participants(),
    condition_selection = EegFun.conditions([1, 2]),  # Compare conditions 1 and 2
    channel_selection = EegFun.channels(),            # All channels
    sample_selection = EegFun.samples((-0.2, 1.5)),   # Time window
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
```


## See Also {#See-Also}
- [API Reference](../reference/index.md)
  
