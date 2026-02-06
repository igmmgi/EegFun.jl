# Plot Joint Probability

This demo demonstrates using joint probability analysis to detect bad channels based on multi-dimensional statistical outliers.

This demo demonstrates using joint probability analysis to detect bad channels based on multi-dimensional statistical outliers.

### What is Joint Probability?

Joint probability analysis detects outlier channels by examining how each channel's distribution relates to all other channels simultaneously:

- **Multi-dimensional comparison**: Unlike single metrics (variance, kurtosis), joint probability considers the full distribution
- **Statistical outlier detection**: Channels that don't fit the overall pattern are flagged
- **More sensitive**: Catches unusual channels that may pass single-metric tests

### How It Works

1. **Calculate channel statistics**: Extract features from each channel (amplitude, variance, etc.)
2. **Build joint distribution**: Model the typical pattern across all channels
3. **Compute probabilities**: For each channel, calculate how likely its pattern is
4. **Identify outliers**: Channels with very low probability are suspicious

### Interpreting Results

**Low probability values** (close to 0):

- Channel doesn't fit typical pattern
- Potential bad channel
- May have unusual artifacts or noise

**Consistent probability values** (similar across channels):

- Data quality is uniform
- No obvious outliers

**Spatial clusters** of low probability:

- Systematic recording issues
- Poor contact in a region
- Regional artifacts

### Use Cases

**Pre-processing quality control**:

- Automated bad channel detection
- Objective identification criteria
- Reduces manual inspection time

**Complement other metrics**:

- Use with variance and kurtosis
- Catches different types of artifacts
- More comprehensive quality assessment

**Research validation**:

- Document channel exclusion decisions
- Consistent criteria across datasets
- Reproducible preprocessing

### Threshold Selection

The default threshold identifies statistical outliers, but you may need to adjust based on:

- **Data quality**: Noisy data may need more lenient thresholds
- **Channel count**: More channels → more statistical power
- **Analysis goals**: Conservative vs aggressive cleaning

### Advantages Over Single Metrics

**Variance alone**: Misses channels with unusual distributions but normal variance

**Kurtosis alone**: Misses channels with shifted means or different scales

**Joint probability**: Combines multiple dimensions to catch subtle issues

## Workflow Summary

This demo shows joint probability analysis:

### 1. Load and Preprocess

- Read continuous data
- Apply basic preprocessing (reference, filter)

### 2. Compute Joint Probability

- Calculate multi-dimensional statistics
- Build probability distribution
- Returns metrics for each channel

### 3. Visualize Results

- Plot probability values per channel
- Identify low-probability outliers
- Assess spatial patterns

### 4. Use Results for QC

- Flag suspicious channels
- Combine with other metrics
- Make informed exclusion decisions


## Code Examples

::: details Show Code

```julia
using EegFun

# read raw data
dat = EegFun.read_raw_data("./resources/data/bdf/example1.bdf");

# read and preprate layout file
layout_file = EegFun.read_layout("./resources/layouts/biosemi/biosemi72.csv");

# create EegFun data structure (EegFun.ContinuousData)
dat = EegFun.create_eegfun_data(dat, layout_file);

# Some minimal preprocessing (average reference and highpass filter)
EegFun.rereference!(dat, :avg)
EegFun.highpass_filter!(dat, 1)

# test plot_channel_summary
jp = EegFun.channel_joint_probability(dat)
EegFun.plot_joint_probability(jp)
```

:::

## See Also

- [API Reference](../reference/index.md)
