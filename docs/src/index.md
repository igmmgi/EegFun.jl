````@raw html
---
layout: home
hero:
  name: EegFun.jl
  text: EEG/ERP analysis in Julia
  image:
    src: /EegFunLogo.png
  actions:
    - theme: alt
      text: Get Started
      link: /tutorials/getting-started
    - theme: alt
      text: View on GitHub
      link: https://github.com/igmmgi/EegFun.jl

features:
  - icon: 
      src: /eeg-cap.svg
    title: Raw data readers
    details: Biosemi, BrainVision Analyser, ...
  - icon: 
      src: /epochs.svg
    title: Data Processing
    details: filtering, referencing, artifact detection/correction, ICA, epoching ...
  - icon: 
      src: /eeg-wave.svg
    title: Analysis
    details: epochs, time-frequency analysis, ERP measurements, statistical testing
  - icon:
      src: /plots-grid.svg
    title: Visualization
    details: Interactive plots, topographic maps, and publication-quality figures with Makie.jl
---
````
## Quick Start

Install EegFun.jl from the Julia REPL:

```julia
using Pkg
Pkg.add("EegFun")
```

Load and preprocess EEG data:

```julia
using EegFun
using GLMakie  # For plotting

# Load EEG data
dat = EegFun.read_raw_data("your_data.bdf")
layout = EegFun.read_layout("biosemi64.csv")
dat = EegFun.create_eeg_dataframe(dat, layout)

# Basic preprocessing
EegFun.highpass_filter!(dat, 1)      # High-pass filter at 1 Hz
EegFun.rereference!(dat, :avg)       # Average reference
EegFun.is_extreme_value!(dat, 100)   # Mark extreme values

# Create epochs and compute ERPs
epoch_cfg = [EegFun.EpochCondition(name = "Target", trigger_sequences = [[1]])]
epochs = EegFun.extract_epochs(dat, 1, epoch_cfg[1], -0.2, 0.8)
erps = EegFun.average_epochs(epochs)
fig, ax = EegFun.plot_erp(erps)
```

## Documentation

:::tip Learn EegFun.jl
Start with our [Getting Started Tutorial](tutorials/getting-started.md) to learn the basics
:::

| Section | Description |
|---------|-------------|
| [Tutorials](tutorials/getting-started.md) | Step-by-step guides to learn EegFun.jl from scratch |
| [How-To Guides](how-to/filter-data.md) | Task-focused solutions to specific problems |
| [Explanations](explanations/data-structures.md) | Conceptual deep-dives into EEG analysis |
| [API Reference](reference/index.md) | Complete function and type documentation |

## Getting Help
- Report bugs on [GitHub Issues](https://github.com/igmmgi/EegFun.jl/issues)
