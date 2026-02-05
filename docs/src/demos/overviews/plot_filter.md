This demo demonstrates visualizing filter frequency and phase responses to verify filter characteristics before applying them to data.

### What is Filter Visualization?

Filter response plots show how a filter affects different frequencies:

- **Magnitude response**: Attenuation (in dB) at each frequency
- **Phase response**: Time delays introduced across frequencies
- **Cutoff characteristics**: Transition band steepness and rolloff

### Why Visualize Filters?

**Verify design parameters**:

- Confirm cutoff frequencies are correct
- Check passband and stopband behavior
- Ensure appropriate attenuation levels

**Identify potential issues**:

- Excessive ripple in passband
- Slow transition bands
- Phase distortion effects

**Documentation**:

- Include in methods sections
- Show filter characteristics clearly
- Support reproducibility

### Filter Types Supported

**High-pass filters**:

- Remove slow drifts and DC offset
- Typical: 0.1-1 Hz cutoff
- Preserve task-related activity

**Low-pass filters**:

- Remove high-frequency noise
- Typical: 30-40 Hz cutoff
- Anti-aliasing before downsampling

**Band-pass filters**:

- Isolate specific frequency ranges
- Combine high-pass and low-pass
- Focus on frequency bands of interest

### Filter Methods

**IIR (Infinite Impulse Response)**:

- Butterworth filters
- Efficient computation
- Steeper rolloff with fewer coefficients
- Can introduce phase distortion

**FIR (Finite Impulse Response)**:

- Linear phase (no distortion)
- Requires more coefficients
- Computationally more expensive
- Symmetric impulse response

### Visualization Features

**Reference lines**:

- Common attenuation levels (-3dB, -6dB, -12dB)
- Highlight filter characteristics
- Customizable positions and styling

**Customization**:

- Line colors and widths
- Title and labels
- Frequency resolution (n_points)
- Reference line positions

## Workflow Summary

This demo shows filter response visualization:

### 1. Create Filters

- Lowpass IIR filter (30 Hz)
- Highpass IIR filter (1 Hz)
- Lowpass FIR filter (40 Hz)
- Various cutoff frequencies

### 2. Plot Responses

- Basic magnitude response plots
- Compare IIR vs FIR characteristics
- Visualize different cutoffs

### 3. Customize Appearance

- Custom colors and line widths
- Reference lines at specific dB levels
- Titles and styling options

### 4. Verify Characteristics

- Check cutoff frequency accuracy
- Assess transition band steepness
- Evaluate filter suitability
