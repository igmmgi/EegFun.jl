This demo shows how to apply temporal filters to EEG data to remove unwanted frequency components.

### Why Filter?

- **High-pass filter** removes slow drifts and DC offsets that obscure the signal of interest
- **Low-pass filter** removes high-frequency noise (muscle activity, line noise harmonics)
- Together they define the frequency band you want to analyse

### Key Functions

| Function | Purpose | Typical Use |
| --- | --- | --- |
| `highpass_filter!` | Remove frequencies below cutoff | 0.1 Hz for ERPs, 1 Hz for ICA |
| `lowpass_filter!` | Remove frequencies above cutoff | 30–40 Hz for ERPs |
| `create_highpass_filter` | Create a filter object for inspection | Visualize filter response |
| `create_lowpass_filter` | Create a filter object for inspection | Compare filter settings |
| `plot_filter` | Plot filter frequency response | Verify filter characteristics |

### IIR vs FIR

- **IIR** (default): Butterworth filter — fast, minimal parameters, good for most uses
- **FIR**: Windowed-sinc filter — linear phase, better when precise timing matters

## Workflow Summary

### High-Pass Filtering

- Standard 0.1 Hz for ERP analysis
- Stronger 1 Hz for ICA preprocessing

### Low-Pass Filtering

- 30 Hz for typical ERP analysis
- Anti-aliasing before resampling

### Channel-Specific Filtering

- Apply filters to selected channels only

### Filter Inspection

- Create filter objects and visualise frequency response
- Compare IIR vs FIR characteristics
