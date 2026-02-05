# Decoding

This demo demonstrates time-resolved multivariate pattern analysis (MVPA) for decoding experimental conditions from EEG data.

This demo demonstrates time-resolved multivariate pattern analysis (MVPA) for decoding experimental conditions from EEG data.

### What is MVPA Decoding?

MVPA uses machine learning classifiers to decode experimental conditions from spatial patterns of brain activity:

- **Multivariate**: Uses all channels simultaneously (not one channel at a time)
- **Pattern analysis**: Finds distributed spatial patterns that discriminate conditions
- **Time-resolved**: Decodes at each time point separately to track when information is available

This reveals **when** and **how well** neural patterns can distinguish between experimental conditions.

### Why Use Decoding?

**Information content**:

Decoding tells you **when information is present** in neural patterns, even if traditional ERPs don't show clear differences.

**Spatial patterns**:

Uses distributed activity across channels, potentially more sensitive than univariate approaches.

**Temporal dynamics**:

Track when discriminative information emerges, peaks, and decays across the trial.

### Workflow

**1. Prepare data**:

```julia
participant_epochs = prepare_decoding(
    "epochs_good",
    condition_selection = conditions([1, 2]),
    sample_selection = samples((-0.2, 1.5))
)
```

**2. Run decoding**:

```julia
all_decoded = decode_libsvm(
    participant_epochs,
    n_iterations = 5,
    n_folds = 3
)
```

Uses cross-validated support vector machine (SVM) classification.

**3. Grand average**:

```julia
grand_avg = EegFun.grand_average(all_decoded)
```

Average decoding accuracy across participants.

**4. Statistical testing**:

```julia
stats = test_against_chance(all_decoded, alpha = 0.05)
stats_cluster = test_against_chance_cluster(all_decoded, alpha = 0.05)
```

**5. Visualization**:

```julia
plot_decoding(grand_avg, stats)
```

### Cross-Validation

Decoding uses **k-fold cross-validation**:

1. Split data into k folds
2. Train on k-1 folds
3. Test on held-out fold
4. Repeat for all folds
5. Average accuracy across folds

This prevents overfitting and gives unbiased accuracy estimates.

### Statistical Testing

**Multiple Comparison Correction**:

| Method | Description |
|--------|-------------|
| **:none** | No correction (liberal) |
| **:bonferroni** | Divide alpha by number of time points (conservative) |
| **:cluster** | Cluster-based permutation testing (recommended) |

**Cluster-based testing**:

Identifies contiguous time windows where decoding is above chance while controlling family-wise error rate.

### Interpreting Results

**Decoding accuracy**:

- **50%** = Chance level (for 2-class problems)
- **60-70%** = Moderate decoding (information present)
- **>80%** = Strong decoding (highly discriminative patterns)

**Temporal profile**:

- **Early peaks** (< 200 ms): Sensory processing
- **Mid-latency** (200-400 ms): Perceptual/cognitive processing
- **Late sustained** (> 400 ms): Decision-making, motor preparation

**Significance**:

Only interpret time points that survive statistical testing with appropriate correction.

### Demo Structure

**

Synthetic data**:

Creates artificial data with controllable signal-to-noise ratio to validate the pipeline.

**Real data**:

Applies decoding to actual experimental data, comparing two conditions across participants.

**Multiple corrections**:

Demonstrates different statistical correction methods for comparison.

### Best Practices

**Data requirements**:

- **Balanced classes**: Equal number of trials per condition (use `equalize_trials = true`)
- **Sufficient trials**: At least 30-50 trials per condition
- **Clean data**: Artifact rejection before decoding

**Cross-validation settings**:

- **n_folds**: 3-10 folds (fewer for small trial counts)
- **n_iterations**: 5-20 iterations (more = more stable, but slower)

**Statistical testing**:

- Use cluster-based correction as default
- Bonferroni is very conservative for time-series data
- Report corrected p-values and time windows


## Code Examples

::: details Show Code

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

:::

## See Also

- [API Reference](../reference/index.md)
