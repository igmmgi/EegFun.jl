# EeggFun.jl Documentation

Welcome to EegFun.jl, a comprehensive Julia package for EEG data analysis and processing.

## Overview

EegFun.jl provides a complete toolkit for analyzing electroencephalogram (EEG) data, including:

- **Data Loading**: Support for various EEG file formats
- **Preprocessing**: Filtering, referencing, artifact removal
- **Analysis**: Time-frequency analysis, connectivity, statistics
- **Visualization**: Interactive plots and topographic maps
- **Event-Related Potentials**: ERP analysis and visualization

## Example Plots

Here's a simple example of what you can create with EegFun.jl:

![Simple Plot Example](assets/images/simple_plot.png)

## Getting Started

### Installation

```julia
using Pkg
Pkg.add("EegFun")
```

### Quick Start

```julia
using EegFun

# Load EEG data
data = load_eeg("your_file.bdf")

# Basic preprocessing
data = filter_data(data, 0.1, 40.0)  # Bandpass filter
data = reference_data(data, :average)  # Average reference
data = remove_bad_channels(data)       # Remove bad channels

# Create epochs around events
epochs = epoch_data(data, events, (-0.2, 0.8))

# Compute and plot ERPs
erp = compute_erp(epochs)
plot_erp(erp, ["Fz", "Cz", "Pz"])
```

## Documentation Structure

- **[API Reference](api/types.md)**: Complete function and type documentation
- **[Examples](examples.md)**: Practical usage examples and tutorials

For more detailed information, see the API reference and examples sections. 