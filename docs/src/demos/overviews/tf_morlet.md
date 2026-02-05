This demo demonstrates time-frequency analysis using Morlet wavelets to decompose EEG signals into time-varying frequency content.

### What are Morlet Wavelets?

Morlet wavelets are oscillating waves with Gaussian envelopes that provide localized time-frequency analysis:

- **Wavelet**: Small wave that is both time-limited and frequency-specific
- **Gaussian envelope**: Smooth windowing that minimizes spectral leakage
- **Time-frequency trade-off**: Cannot have perfect precision in both simultaneously

### Key Parameter: Number of Cycles

The `cycles` parameter controls the time-frequency resolution trade-off:

**Fixed cycles** (e.g., `cycles = 5`):
- Same spectral bandwidth at all frequencies
- **Low frequencies**: Broad temporal windows (poor time resolution)
- **High frequencies**: Narrow temporal windows (good time resolution)

**Variable cycles** (e.g., `cycles = (3, 10)`):
- Adaptive resolution across frequencies
- **Low frequencies**: Fewer cycles = better temporal precision
- **High frequencies**: More cycles = better spectral precision
- Balances resolution optimally across the spectrum

### Frequency Spacing

**Linear spacing** (`frequencies = 1:1:40`):
- Equal Hz steps
- Good for narrow frequency ranges
- Easier interpretation

**Logarithmic spacing** (`frequencies = logrange(1, 40, length = 30)`):
- Proportional frequency steps
- Better for wide ranges (e.g., 2-80 Hz)
- Use with `ylogscale = true` in plots

### Baseline Correction

Apply baseline correction to isolate event-related changes:

```julia
baseline_window = (-0.5, -0.2)  # Pre-stimulus baseline
baseline_method = :db           # Decibel conversion (10*log10(activity/baseline))
```

**Baseline methods**:
- **:db**: Decibel scale (relative power change)
- **:percent**: Percent change from baseline
- **:zscore**: Z-score normalization

### Edge Effects

The `filter_edges` and `pad` parameters control edge artifacts:

```julia
tf_morlet(epochs, pad = :both)        # Mirror padding
tf_morlet(epochs, filter_edges = true) # Remove edge samples
```

### Workflow Summary

This demo shows:

1. **Synthetic validation**: Known frequencies (5, 25, 35 Hz) for sanity checking
2. **Fixed cycles**: Constant resolution across frequencies
3. **Variable cycles**: Adaptive resolution with `cycles = (3, 10)`
4. **Linear vs log spacing**: Different frequency sampling strategies
5. **Baseline correction**: Isolating event-related power changes

### Advantages

- **Flexible resolution**: Adapt time-frequency trade-off via cycles
- **Standard method**: Widely used in EEG research
- **Interpretable**: Direct relationship between cycles and resolution

### Common Applications

- Event-related synchronization/desynchronization (ERS/ERD)
- Oscillatory dynamics (alpha, theta, gamma)
- Phase-amplitude coupling
- Induced vs evoked responses
