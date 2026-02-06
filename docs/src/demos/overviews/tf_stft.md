This demo demonstrates time-frequency analysis using the Short-Time Fourier Transform (STFT), a classic approach to spectrograms.

### What is STFT?

STFT applies the Fast Fourier Transform (FFT) to windowed segments of the signal:

1. **Window the signal**: Apply time-limited window (e.g., Hann, Tukey)
2. **Compute FFT**: Calculate frequency spectrum for that segment
3. **Slide window**: Move forward in time and repeat
4. **Build spectrogram**: Stack spectra to create time-frequency representation

### Key Parameters

**Window length** (`window_length`):

Controls time-frequency resolution:

```julia
tf_stft(epochs, window_length = 0.5)  # 500 ms windows
```

- **Longer windows** (e.g., 500 ms): Better frequency resolution, worse time resolution
- **Shorter windows** (e.g., 200 ms): Better time resolution, worse frequency resolution

**Cycles** (alternative):

Instead of specifying window length, specify cycles at each frequency:

```julia
tf_stft(epochs, cycles = 7)
```

Creates frequency-adaptive windows similar to Morlet wavelets.

**Time steps**:

Controls temporal resolution of output:

```julia
tf_stft(epochs, time_steps = 0.005)  # 5 ms steps
```

Smaller = smoother time course, larger = faster computation.

### Frequency Spacing

**Linear** (`frequencies = 2:1:80`):

- Equal Hz spacing
- Natural for narrow bands

**Logarithmic** (`frequencies = logrange(2, 80, length = 100)`):

- Proportional spacing
- Better for wide ranges
- Use `ylogscale = true` in plots

### Baseline Correction

```julia
plot_time_frequency(
    tf_data,
    baseline_window = (-0.5, -0.2),
    baseline_method = :db,
    colorrange = (-3, 3)
)
```

Shows relative power changes from pre-stimulus baseline.

### Edge Effects

**Filter edges** (`filter_edges = true`):

- Removes edge samples where windowing causes artifacts
- Recommended for cleaner results

### Workflow Summary

This demo shows:

1. **Synthetic validation**: Known frequencies (5, 25, 35 Hz)
2. **Fixed window length**: Constant resolution analysis
3. **Cycles mode**: Adaptive resolution
4. **Linear vs log spacing**: Different frequency sampling
5. **Baseline correction**: Event-related power changes
6. **Edge filtering**: Removing edge artifacts

### Further Reading

Cohen, M. X. (2014). *Analyzing Neural Time Series Data: Theory and Practice*. Chapter 15: Short-Time FFT
