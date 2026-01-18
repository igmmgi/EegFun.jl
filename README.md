# eegfun

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://igmmgi.github.io/eegfun.jl/)
[![Build Status](https://github.com/igmmgi/eegfun/workflows/Documentation/badge.svg)](https://github.com/igmmgi/eegfun/actions)
[![CI](https://github.com/igmmgi/eegfun/workflows/Tests/badge.svg)](https://github.com/igmmgi/eegfun/actions)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

Alpha level 0.1 :-)

## ✅ Complete

### Data Loading & Preprocessing
- ✓ Basic EEG data reading (Biosemi BDF, BrainVision)
- ✓ Layout system for different electrode configurations
- ✓ Filtering (high-pass, low-pass, bandpass, bandstop)
- ✓ Baseline correction
- ✓ Rereferencing (average, single channel, custom)
- ✓ Resampling
- ✓ Channel repair (neighbor interpolation, spherical spline)
- ✓ Artifact detection and rejection
- ✓ ICA component calculation/removal

### Analysis
- ✓ Epoching and ERP analysis
- ✓ Time-frequency analysis (wavelet transforms)
- ✓ Cluster-based permutation tests (spatial, temporal, spatiotemporal)
- ✓ Analytic t-tests with multiple comparison correction
- ✓ Decoding/MVPA (multivariate pattern analysis) with MLJ integration
- ✓ RSA (Representational Similarity Analysis) with static and temporal models
- ✓ GFP (Global Field Power)
- ✓ Grand averaging
- ✓ Jackknife averaging
- ✓ ERP measurements (peak detection, latency, amplitude)

### Visualization
- ✓ Interactive databrowser with analysis options (Makie)
- ✓ Multiple plot types (ERP, topography, spectrum, time-frequency, etc.)
- ✓ Interactive plots with zooming, panning, channel selection
- ✓ Decoding accuracy plots
- ✓ RSA visualization (RDMs, model correlations)

### Utilities
- ✓ Batch processing pipelines
- ✓ Configuration system (TOML-based)
- ✓ Logging system
- ✓ Data export/import (JLD2)

## 📋 TODO / Known Issues

- ☐ Batch functions need to deal with bad input options (e.g., conditions missing)
- ☐ Add tests for decoding/MVPA functionality
- ☐ Add tests for RSA functionality
- ☐ Add more file formats (EDF, SET, FIF)
- ☐ Improve documentation and examples
- ☐ Consolidate codebase (see CONSOLIDATION_PLAN.md)

