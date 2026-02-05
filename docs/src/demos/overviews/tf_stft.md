## Overview

This demo demonstrates time-frequency analysis using the Short-Time Fourier Transform (STFT).

### STFT Method

Classic approach to time-frequency decomposition:
1. **Window the signal**: Apply time-limited window (e.g., Hann, Hamming)
2. **Compute FFT**: Calculate spectrum for windowed segment
3. **Slide window**: Move window and repeat
4. **Build spectrogram**: Stack spectra across time

### Key Parameters

**Window size:**
- Longer windows: Better frequency resolution, worse time resolution
- Shorter windows: Better time resolution, worse frequency resolution
- Typical: 200-500 ms windows for EEG

**Overlap:**
- Amount of window overlap (e.g., 50%)
- More overlap: Smoother time course
- Less overlap: Faster computation

**Window type:**
- Hann, Hamming, Tukey windows
- Controls spectral leakage

### Applications

- Simple, fast time-frequency analysis
- Broad-band power changes
- Exploratory analyses
- When constant resolution across frequencies is acceptable

### Trade-offs

- Fixed time-frequency resolution (cf. wavelets which adapt)
- Simpler than wavelets, but less flexible
- Widely understood and implemented
