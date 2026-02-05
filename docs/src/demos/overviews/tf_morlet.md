## Overview

This demo demonstrates time-frequency analysis using Morlet wavelets.

### Morlet Wavelet Transform

Decompose signals into time-varying frequency content:
- **Wavelets**: Small waves that are time-limited and frequency-specific
- **Variable resolution**: Trade-off between time and frequency precision
- **Number of cycles**: Control time-frequency resolution

### Cycle Configurations

**Fixed cycles** (e.g., 5 cycles):
- Constant spectral bandwidth
- Broader temporal windows at low frequencies
- Sharper temporal windows at high frequencies

**Variable cycles** (e.g., 3-10 cycles):
- Adaptive resolution across frequencies
- Better low-frequency temporal precision
- Better high-frequency spectral precision

### Applications

- Oscillatory activity (alpha, theta, gamma)
- Induced vs. evoked responses
- Event-related synchronization/desynchronization (ERS/ERD)
- Phase-amplitude coupling

### Advantages

- Excellent time-frequency resolution trade-off
- Widely used standard in EEG research
- Flexible parameterization
