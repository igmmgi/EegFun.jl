# Plot Channel Summary

This demo shows how to compute and visualize channel summary statistics for quality control and data characterization.

This demo shows how to compute and visualize channel summary statistics for quality control and data characterization.

### What is Channel Summary?

Channel summary computes aggregate statistics across time for each channel:

- **Descriptive statistics**: Mean, median, min, max, standard deviation
- **Distribution metrics**: Variance, kurtosis, skewness
- **Range metrics**: Absolute range, interquartile range

These statistics help identify problematic channels, assess data quality, and characterize recording conditions.

### Available Metrics

The `channel_summary()` function computes:

| Metric | Description | 
|--------|-------------|
| **mean** | Average amplitude |
| **median** | Median amplitude |
| **std** | Standard deviation |
| **var** | Variance |
| **zvar** | Z-scored variance |
| **min** | Minimum value |
| **max** | Maximum value |
| **range** | Max - Min |
| **kurtosis** | Distribution tail weight |

### Visualization Options

**Single Metric**:
```julia
plot_channel_summary(cs, :range)
```
Shows one statistic across all channels as a bar plot

**Multiple Metrics**:
```julia
plot_channel_summary(cs, [:min, :max, :std, :range])
```
Creates multiple subplots for comparison

**Epoched Data**:
```julia
plot_channel_summary(cs, :range, average_over = :epoch)
```
Averages statistics across epochs before plotting

### Common Use Cases

**1. Quality Control**:

- **High variance/range**: Noisy or artifactual channels
- **Extreme min/max**: Electrode saturation or poor contact
- **High z-variance**: Outlier channels requiring interpolation

**2. Preprocessing Validation**:

- **After filtering**: Verify reduced variance in target frequencies
- **After baseline correction**: Check mean values near zero
- **After rereferencing**: Assess reference choice effectiveness

**3. Hardware Issues**:

- **Systematic patterns**: Poor grounding, bridging between neighboring channels
- **Single channel anomalies**: Loose connection, bad electrode

### Interpretation

**Healthy channels** typically show:

- Moderate, similar variance across scalp
- No extreme outliers in min/max
- Range appropriate for recording gain settings

**Problematic patterns**:

**Very high variance** (>2-3× median):

- Likely artifact (muscle, movement, poor contact)
- Consider interpolation or exclusion

**Very low variance** (<0.5× median):

- May indicate flat/dead channel
- Check hardware connections

**Extreme values** (saturated at amplifier limits):

- Indicates clipping/saturation
- Cannot be recovered, must exclude

### Workflow Summary

This demo demonstrates:

1. **Load and preprocess** continuous data
2. **Compute summaries** with flexible channel/sample selection
3. **Plot single or multiple metrics** for visualization
4. **Apply to epoched data** with epoch averaging

The channel selection examples show how to:

- Select specific channels by name
- Use pattern matching (e.g., midline channels ending in "z")
- Apply custom sample selections


## Code Examples

::: details Show Code

```julia
using EegFun

# read raw data
dat = EegFun.read_raw_data("./resources/data/bdf/example1.bdf");

# read and preprate layout file
layout_file = EegFun.read_layout("./resources/layouts/biosemi/biosemi72.csv");

# create EegFun data structure (EegFun.ContinuousData)
dat = EegFun.create_eeg_dataframe(dat, layout_file);

# Some minimal preprocessing (average reference and highpass filter)
EegFun.rereference!(dat, :avg)
EegFun.highpass_filter!(dat, 1)

# basic channel summary statistics
cs = EegFun.channel_summary(dat)
cs = EegFun.channel_summary(dat, channel_selection = EegFun.channels([:Fp1, :Fp2]))
cs = EegFun.channel_summary(dat, channel_selection = EegFun.channels([:Fp1, :Fp2]), sample_selection = x -> x.sample .< 2000)
cs = EegFun.channel_summary(dat, channel_selection = x -> endswith.(string.(x), "z")) # all midline channels 
cs = EegFun.channel_summary(dat, channel_selection = x -> .!(endswith.(string.(x), "z"))) # all non-midline channels 

# Plotting Channel Summaries
EegFun.plot_channel_summary(cs, :range)
EegFun.plot_channel_summary(cs, :min)
EegFun.plot_channel_summary(cs, :min, bar_color = :red)
EegFun.plot_channel_summary(cs, [:min, :max, :std, :range, :var, :zvar])

#################################
# Epoched DataFrameEeg
#################################
# some epoched data
epoch_cfg = [EegFun.EpochCondition(name = "ExampleEpoch1", trigger_sequences = [[1]])]
epochs = EegFun.extract_epochs(dat, epoch_cfg, (-2, 4))

cs = EegFun.channel_summary(epochs[1])

EegFun.plot_channel_summary(cs, :range, average_over = :epoch)
EegFun.plot_channel_summary(cs, [:min, :max, :std, :range, :var, :zvar], average_over = :epoch)
```

:::

## See Also

- [API Reference](../reference/index.md)
