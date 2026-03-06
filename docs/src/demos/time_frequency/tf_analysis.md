# TF Analysis Workflow

This demo shows how to work with time-frequency data after decomposition — including condition operations, channel operations, subsetting, and batch processing.

## Prerequisites

Time-frequency data (`TimeFreqData`) is produced by one of the decomposition methods:

- [`tf_morlet`](tf_morlet.md) — Morlet wavelet decomposition
- [`tf_multitaper`](tf_multitaper.md) — Multitaper spectral estimation
- [`tf_stft`](tf_stft.md) — Short-Time Fourier Transform

Each `TimeFreqData` object contains two DataFrames: `data_power` (spectral power) and `data_phase` (phase angles). All operations below act on **both** DataFrames consistently.

## Key Functions

| Function | Purpose | Mutating? |
| --- | --- | --- |
| `condition_average` | Average power/phase across conditions | No |
| `condition_difference` | Subtract one condition's TF from another | No |
| `channel_average!` / `channel_average` | Average across channels | Yes / No |
| `channel_difference!` / `channel_difference` | Compute channel differences | Yes / No |
| `subset` | Select channels, time points, or conditions | No |

## Condition Operations

Condition averaging and differencing work identically to their ERP counterparts — they auto-detect `TimeFreqData` and dispatch appropriately:

- **Averaging**: Combine conditions into a single mean TF representation
- **Differencing**: Compute condition contrasts (e.g., target minus standard)

## Channel Operations

- **`channel_average`**: Create region-of-interest (ROI) averages (e.g., mean frontal power)
- **`channel_difference`**: Compute laterality indices (e.g., left − right)
- Both update `data_power` and `data_phase` consistently

## Batch Processing

The batch versions of `condition_average` and `condition_difference` automatically detect whether JLD2 files contain `ErpData` or `TimeFreqData` and dispatch accordingly. Use the same API for both:

```julia
# Works for both ERP and TF files
EegFun.condition_average("tf_morlet", [[1, 2], [3, 4]], input_dir = "./output")
EegFun.condition_difference("tf_morlet", [(1, 2), (3, 4)], input_dir = "./output")
```

## Subsetting

Use `subset` to select specific channels, time intervals, or conditions. For `TimeFreqData`, both `data_power` and `data_phase` are subset consistently:

```julia
sub = EegFun.subset(tf, channel_selection = EegFun.channels([:Fz, :Cz]))
sub = EegFun.subset(tfs, condition_selection = EegFun.conditions(1, 2))
```

## Code Examples

::: details Show Code

```julia
# Demo: TF Analysis Workflow
# Shows how to perform condition and channel operations on time-frequency data.

using EegFun

#######################################################################
@info EegFun.section("1. Generate TF data from epochs")
#######################################################################

# Load and preprocess data
dat = EegFun.read_raw_data("./resources/data/bdf/example1.bdf")
layout = EegFun.read_layout("./resources/layouts/biosemi/biosemi72.csv")
EegFun.polar_to_cartesian_xy!(layout)
dat = EegFun.create_test_eegfun_data(dat, layout)
EegFun.rereference!(dat, :avg)
EegFun.highpass_filter!(dat, 0.1)

# Extract epochs for four conditions
epoch_cfg = [
    EegFun.EpochCondition(name = "Cond1", trigger_sequences = [[1]]),
    EegFun.EpochCondition(name = "Cond2", trigger_sequences = [[2]]),
    EegFun.EpochCondition(name = "Cond3", trigger_sequences = [[3]]),
    EegFun.EpochCondition(name = "Cond4", trigger_sequences = [[5]]),
]
epochs = EegFun.extract_epochs(dat, epoch_cfg, (-0.5, 1.5))

# Morlet decomposition per condition
tfs = [EegFun.tf_morlet(epochs[i], frequencies = logrange(2, 40, length = 30), cycles = (3, 7)) for i in eachindex(epochs)]

# Plot single condition
EegFun.plot_tf(
    tfs[1];
    baseline_interval = (-0.4, -0.1),
    baseline_method = :db,
    xlim = (-0.2, 1.0),
    colorrange = (-3, 3),
)

#######################################################################
@info EegFun.section("2. Condition Average")
#######################################################################

# Average conditions 1+2 and 3+4
avg_tfs = EegFun.condition_average(tfs, [[1, 2], [3, 4]])
@info "Created $(length(avg_tfs)) average TF objects"

EegFun.plot_tf(
    avg_tfs[1];
    baseline_interval = (-0.4, -0.1),
    baseline_method = :db,
    xlim = (-0.2, 1.0),
    colorrange = (-3, 3),
)

#######################################################################
@info EegFun.section("3. Condition Difference")
#######################################################################

# Difference: condition 1 minus condition 2
diff_tfs = EegFun.condition_difference(tfs, [(1, 2), (3, 4)])
@info "Created $(length(diff_tfs)) difference TF objects"

EegFun.plot_tf(
    diff_tfs[1];
    baseline_interval = (-0.4, -0.1),
    baseline_method = :db,
    xlim = (-0.2, 1.0),
    colorrange = (-3, 3),
)

#######################################################################
@info EegFun.section("4. Channel Operations")
#######################################################################

# Average across frontal channels
tf_frontal = EegFun.channel_average(
    tfs[1],
    channel_groups = [[:Fz, :F3, :F4]],
    output_labels = [:Frontal],
    reduce = true,
)

EegFun.plot_tf(
    tf_frontal;
    baseline_interval = (-0.4, -0.1),
    baseline_method = :db,
    xlim = (-0.2, 1.0),
)

# Channel difference (laterality)
tf_lat = EegFun.channel_difference(
    tfs[1],
    channels_in1 = [:C3],
    channels_in2 = [:C4],
    channel_out = :Laterality,
)

#######################################################################
@info EegFun.section("5. Subsetting")
#######################################################################

# Subset by channels
tf_sub = EegFun.subset(tfs[1], channel_selection = EegFun.channels([:Fz, :Cz, :Pz]))
@info "Subset has $(length(EegFun.channel_labels(tf_sub))) channels"

# Subset vector by condition
tf_sub_vec = EegFun.subset(tfs, condition_selection = EegFun.conditions(1, 2))
@info "Subset has $(length(tf_sub_vec)) conditions"
```

:::

## See Also

- [TF Morlet](tf_morlet.md) — Morlet wavelet decomposition
- [TF Multitaper](tf_multitaper.md) — Multitaper spectral estimation
- [TF STFT](tf_stft.md) — Short-Time Fourier Transform
- [Condition Operations](../erp/condition_operations.md) — ERP condition operations
- [API Reference](../../reference/index.md)
