# eegfun Documentation

Welcome to the eegfun documentation! This package provides comprehensive tools for EEG data analysis in Julia.

## Quick Start

```julia
using eegfun

# Load your EEG data
dat = load_eeg_data("your_data.bdf", "biosemi64.csv")

# Basic analysis
summary = channel_summary(dat)
corr_matrix = correlation_matrix(dat)

# Quality control
is_extreme_value!(dat, 100)
detect_eog_onsets!(dat, 50.0, :vEOG, :is_vEOG)

# ICA
ica_result = run_ica(dat, samples = samples_not(:is_extreme_value_100))
```

## Documentation Sections

- **[Data Loading & Types](@ref)** - How to load and work with EEG data
- **[Preprocessing](@ref)** - Filtering, baseline correction, and artifact detection
- **[Analysis](@ref)** - Channel summaries, correlation matrices, and statistical analysis
- **[ICA](@ref)** - Independent Component Analysis
- **[Plotting](@ref)** - Visualization tools for EEG data
- **[Utilities](@ref)** - Helper functions and utilities

## Key Features

- **Flexible filtering**: Use predicate functions for channel and sample selection
- **Quality control**: Built-in artifact detection and rejection
- **Modern plotting**: Makie-based visualization
- **Comprehensive analysis**: From basic summaries to advanced ICA 