# EegFun.jl

Welcome to EegFun.jl, a Julia package for EEG data analysis and processing.

## Overview

EegFun.jl provides a toolkit for analyzing electroencephalogram (EEG) data, including:

- **Data Loading**: Support for Biosemi BDF files and other EEG formats
- **Preprocessing**: Filtering, referencing, artifact detection and removal
- **Analysis**: Time-frequency analysis, ICA, channel statistics
- **Visualization**: Interactive plots, topographic maps, and ERP visualizations
- **Event-Related Potentials**: ERP analysis and visualization

## Installation

```julia
using Pkg
Pkg.add("EegFun")
```

## Quick Start

```julia
using EegFun
using GLMakie  # For plotting

# Load EEG data
dat = EegFun.read_raw_data("your_data.bdf")
layout = EegFun.read_layout("biosemi64.csv")
dat = EegFun.create_eeg_dataframe(dat, layout)

# Basic preprocessing
EegFun.highpass_filter!(dat, 1)      # High-pass filter at 1 Hz
EegFun.rereference!(dat, :avg)         # Average reference
EegFun.is_extreme_value!(dat, 100)     # Mark extreme values

# Create epochs around events
epoch_cfg = [EegFun.EpochCondition(name = "Target", trigger_sequences = [[1]])]
epochs = EegFun.extract_epochs(dat, 1, epoch_cfg[1], -0.2, 0.8)

# Compute and plot ERPs
erps = EegFun.average_epochs(epochs)
fig, ax = EegFun.plot_erp(erps)
```

## Documentation Structure

- **[API Reference](api.md)**: Complete function and type documentation

## Index

```@index
```
