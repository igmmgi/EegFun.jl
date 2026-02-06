# Rereference

This demo demonstrates applying different re-referencing schemes to EEG data.

This demo demonstrates applying different re-referencing schemes to EEG data.

### What is Re-referencing?

Re-referencing changes the voltage baseline by subtracting a reference signal from all electrodes. Since EEG measures potential differences, the choice of reference affects:

- Signal amplitudes and topographies
- Component interpretations
- Statistical comparisons

### Why Re-reference?

**Recording reference**: The reference used during data acquisition (e.g., Cz, linked mastoids, or manufacturer default)

**Analysis reference**: The reference appropriate for your research question and field conventions

Re-referencing allows you to change from the recording reference to a more appropriate analysis reference.

### Common Reference Schemes

**Average Reference (`:avg`)**:

- Subtracts mean of all electrodes from each channel
- Assumes zero average scalp potential
- Most common for high-density EEG (64+ channels)
- Reference-free topographies

**Mastoid Reference**:

- References to average of left/right mastoid electrodes (e.g., `:M1`, `:M2`)
- Common in ERP research
- Reduces contamination from reference site activity

**Single Electrode Reference**:

- References to a specific electrode (e.g., `:Cz`, `:Fz`)
- Useful for specific research questions
- Activity at reference site becomes invisible (single reference electrode = 0)

### Best Practices

- Apply the same reference across all conditions and participants
- Re-reference **before** interpolating bad channels (bad channels affect average reference)
- Document your reference choice in publications
- Use average reference for high-density recordings

### Re-referencing is Reversible

You can change references multiple times. EegFun tracks the current reference and updates accordingly.

## Workflow Summary

This demo shows re-referencing workflows:

### 1. Continuous Data

- Load and preprocess data
- Apply different reference schemes (`:Fp1`, `:AF7`)
- Visualize in databrowser to see effects

### 2. Epoched Data

- Create epochs from continuous data
- Apply references to segmented trials
- Data structure automatically updated

### 3. ERP Data

- Average epochs into ERPs
- Apply references to averaged data
- Reference changes preserved across data types


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
EegFun.highpass_filter!(dat, 1)

# Rereference to :Fp1
EegFun.rereference!(dat, :Fp1)
dat.data

# Rereference to :AF7
EegFun.rereference!(dat, :AF7)
dat.data

# We can also see the influence of reference using the databrowser
EegFun.plot_databrowser(dat)


################################
# Epoching
################################

# Create some epoched data
epoch_cfg = [
    EegFun.EpochCondition(name = "ExampleEpoch1", trigger_sequences = [[1]]),
    EegFun.EpochCondition(name = "ExampleEpoch2", trigger_sequences = [[2]]),
]
epochs = EegFun.extract_epochs(dat, epoch_cfg, (-0.2, 1.0))  # -200 to 1000 ms

EegFun.rereference!(epochs, :Fp1)
EegFun.all_data(epochs[1]) # first condition
EegFun.all_data(epochs[2]) # second condition


################################
# ERPs
################################
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
