This demo demonstrates time-frequency analysis using the multitaper method for high-quality spectral estimation.

### What is the Multitaper Method?

The multitaper method improves spectral estimates by averaging across multiple orthogonal tapers (windows):

- **Slepian tapers**: Mathematically optimal windows that minimize spectral leakage
- **Multiple estimates**: Each taper provides an independent spectral estimate
- **Averaging**: Reduces variance while preserving spectral features
- **Time-bandwidth product (NW)**: Controls frequency smoothing

### Key Parameters

**Cycles**:

Similar to Morlet wavelets, controls the time window length:

```julia
tf_multitaper(epochs, cycles = 5)
```


### Frequency Spacing

**Linear spacing** (`frequencies = 1:1:40`):

- Equal Hz steps across range
- Good for narrow frequency bands

**Logarithmic spacing** (`frequencies = logrange(1, 40, length = 30)`):

- Proportional steps (constant ratios)
- Better for wide ranges (2-80 Hz)
- Use with `ylogscale = true` for visualization

### Baseline Correction

```julia
plot_time_frequency(
    tf_data,
    baseline_window = (-0.5, -0.2),
    baseline_method = :db
)
```

Isolates event-related changes relative to pre-stimulus baseline.

### Workflow Summary

This demo shows:

1. **Synthetic validation**: Known frequencies (5, 25, 35 Hz)
2. **Fixed cycles**: Constant resolution analysis
3. **Linear vs log spacing**: Different frequency sampling strategies
4. **Baseline correction**: Event-related power changes
5. **Edge filtering**: Removing edge artifacts

### Practical Tips

**Choose cycles based on needs**:

- **3-5 cycles**: Better temporal resolution
- **7-10 cycles**: Better frequency resolution
- **5 cycles**: Good balance for most EEG work

**Frequency sampling**:

- Dense sampling (every 0.5 Hz) for detailed spectra
- Coarse sampling (every 2-5 Hz) for faster computation

**Time steps**:

- Control temporal smoothness of the output
- Smaller steps = smoother time course but slower computation

### Further Reading

Cohen, M. X. (2014). *Analyzing Neural Time Series Data: Theory and Practice*. Chapter 16: Multitapers
