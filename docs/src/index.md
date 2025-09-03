# eegfun.jl

Welcome to eegfun.jl, a comprehensive Julia package for EEG data analysis and processing.

## Overview

eegfun.jl provides a complete toolkit for analyzing electroencephalogram (EEG) data, including:

- **Data Loading**: Support for Biosemi BDF files and other EEG formats
- **Preprocessing**: Filtering, referencing, artifact detection and removal
- **Analysis**: Time-frequency analysis, ICA, channel statistics
- **Visualization**: Interactive plots, topographic maps, and ERP visualizations
- **Event-Related Potentials**: ERP analysis and visualization

## Installation

```julia
using Pkg
Pkg.add("eegfun")
```

## Quick Start

```julia
using eegfun
using GLMakie  # For plotting

# Load EEG data
dat = eegfun.read_bdf("your_data.bdf")
layout = eegfun.read_layout("biosemi64.csv")
dat = eegfun.create_eeg_dataframe(dat, layout)

# Basic preprocessing
eegfun.filter_data!(dat, "hp", 1)      # High-pass filter at 1 Hz
eegfun.rereference!(dat, :avg)         # Average reference
eegfun.is_extreme_value!(dat, 100)     # Mark extreme values

# Create epochs around events
epoch_cfg = [eegfun.EpochCondition(name = "Target", trigger_sequences = [[1]])]
epochs = eegfun.extract_epochs(dat, 1, epoch_cfg[1], -0.2, 0.8)

# Compute and plot ERPs
erps = eegfun.average_epochs(epochs)
fig, ax = eegfun.plot_erp(erps)
```

## Documentation Structure

- **[API Reference](api.md)**: Complete function and type documentation

## Index

```@index
```
