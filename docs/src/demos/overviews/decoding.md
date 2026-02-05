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
