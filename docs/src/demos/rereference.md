# Rereference

## Overview

## Overview

This demo demonstrates different re-referencing schemes for EEG data.

### Why Re-reference?

The reference electrode determines the voltage baseline:
- **Recording reference**: Often active during recording (Cz, mastoids, etc.)
- **Analysis reference**: Should match research question and conventions

### Common Reference Schemes

**Average Reference:**
- Subtract mean of all electrodes
- Assumes zero average scalp potential
- Most common for high-density EEG

**Mastoid Reference:**
- Reference to average of left/right mastoids
- Common in ERP research
- Reduces activity from reference sites

**Linked Mastoids:**
- Average of earlobes/mastoids
- Traditional ERP standard

**Cz Reference:**
- Single central electrode
- Useful for motor/sensory studies

### Best Practices

- Apply same reference across conditions/participants
- Re-reference before interpolating bad channels
- Consider effect on topographies
- Document reference choice in publications


## Code Examples

::: details Show Code

```julia
using EegFun

# read raw data
dat = EegFun.read_raw_data("./resources/data/example1.bdf");

# read and preprate layout file
layout_file = EegFun.read_layout("./resources/layouts/biosemi/biosemi72.csv");
EegFun.polar_to_cartesian_xy!(layout_file)

# create EegFun data structure (EegFun.ContinuousData)
dat = EegFun.create_eeg_dataframe(dat, layout_file);

# Some minimal preprocessing (average reference and highpass filter)
EegFun.highpass_filter!(dat, 1)


EegFun.rereference!(dat, :Fp1)
dat.data

EegFun.rereference!(dat, :AF7)
dat.data

EegFun.rereference!(dat, :Fp1)
dat.data

EegFun.plot_databrowser(dat)


# Create some epoched data
epoch_cfg = [
    EegFun.EpochCondition(name = "ExampleEpoch1", trigger_sequences = [[1]]),
    EegFun.EpochCondition(name = "ExampleEpoch2", trigger_sequences = [[2]]),
]
epochs = EegFun.extract_epochs(dat, epoch_cfg, (-0.2, 1.0))  # -200 to 1000 ms

EegFun.rereference!(epochs, :Fp1)
EegFun.all_data(epochs[1])
EegFun.all_data(epochs[2])

EegFun.rereference!(epochs, :F1)
EegFun.all_data(epochs[1])
EegFun.all_data(epochs[2])


# ERPs
erps = EegFun.average_epochs(epochs)

EegFun.rereference!(erps, :Fp1)
EegFun.all_data(erps[1])
EegFun.all_data(erps[2])

EegFun.rereference!(erps, :F1)
EegFun.all_data(erps[1])
EegFun.all_data(erps[2])```

:::

## See Also

- [API Reference](../reference/index.md)
