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

Read and preprocess EEG data:

```julia
using EegFun

# Load EEG data and layout
dat = EegFun.read_raw_data("your_data.bdf")
layout = EegFun.read_layout("biosemi64.csv")

# Create Eegfun datatype
dat = EegFun.create_eeg_dataframe(dat, layout)

# Basic preprocessing
EegFun.highpass_filter!(dat, 1)      # High-pass filter at 1 Hz
EegFun.rereference!(dat, :avg)       # Average reference
EegFun.is_extreme_value!(dat, 100)   # Mark extreme values

# Continuous data Browser
Eegfun.plot_databrowser(dat)

# Create epochs and compute ERPs
epoch_cfg = [
  EegFun.EpochCondition(name = "Cond1", trigger_sequences = [[1]]),
  EegFun.EpochCondition(name = "Cond2", trigger_sequences = [[2]])
  ]
epochs = EegFun.extract_epochs(dat, epoch_cfg, (-0.2, 0.8))

# Epoch data browser
EegFun.plot_databrowser(epochs[1]) # Cond1

# Plot epochs
EegFun.plot_epochs(epochs, layout = :grid)
EegFun.plot_epochs(epochs, layout = :topo)
EegFun.plot_epochs(epochs, layout = :single, channel_selection = EegFun.channels(:PO7))
EegFun.plot_epochs(epochs, layout = :single, channel_selection = EegFun.channels([:PO7, :PO8]))

# ERPs
erps = EegFun.average_epochs(epochs)

EegFun.plot_erp(erps, layout = :grid)
EegFun.plot_erp(erps, layout = :topo)
EegFun.plot_topoplot(erps, interval_selection = (0.1, 0.2)) # between 100 and 200 msA
```

## Documentation

:::tip Learn EegFun.jl
[Getting Started Tutorial](tutorials/getting-started.md) 
:::

| Section | Description |
|---------|-------------|
| [Tutorials](tutorials/getting-started.md) | Step-by-step guides |
| [API Reference](reference/index.md) | Complete function and type documentation |

## Getting Help
- Report bugs on [GitHub Issues](https://github.com/igmmgi/EegFun.jl/issues)
