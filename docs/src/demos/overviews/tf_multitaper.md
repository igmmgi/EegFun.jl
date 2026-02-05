## Overview

This demo demonstrates time-frequency analysis using the multitaper method.

### Multitaper Method

Reduce variance in spectral estimates using multiple tapers:
- **Tapers**: Data windows with optimal properties (Slepian sequences)
- **Multiple estimates**: Average across tapers for reduced variance
- **Time-bandwidth product**: Controls frequency smoothing

### Key Parameters

**Time-bandwidth product (NW):**
- Controls spectral concentration
- Typical values: 2-4
- Higher = more frequency smoothing, less variance

**Number of tapers:**
- Usually 2×NW - 1
- More tapers = smoother spectra

### Advantages

- Superior variance reduction vs. single-taper methods
- Minimal spectral leakage
- Optimal for narrow-band analyses
- Gold standard for power spectral estimation

### Applications

- High-quality power spectra
- Coherence analysis
- Stationary oscillatory activity
- Resting-state analyses
