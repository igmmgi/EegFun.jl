## Overview

This demo shows how to resample EEG data to different sampling rates.

### Why Resample?

Change the sampling rate for various reasons:
- **Downsampling**: Reduce file size, speed up processing
- **Upsampling**: Match sampling rates across datasets  
- **Standard rates**: Convert to conventional rates (e.g., 250, 500, 1000 Hz)

### Considerations

**Downsampling:**
- Apply low-pass filter first (anti-aliasing)
- Common target: 250-500 Hz for ERP studies
- Avoid aliasing high-frequency content

**Upsampling:**
- Cannot add information
- Useful for synchronization
- Interpolates between samples

### Best Practices

- Downsample after filtering
- Keep rate above 2× highest frequency of interest (Nyquist)
- Document original and resampled rates
- Typical ERP range: 250-500 Hz adequate
