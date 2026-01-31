# How to Filter EEG Data

This guide shows you how to apply different types of filters to your EEG data.

## High-Pass Filtering

Remove slow drifts and DC offsets:

```julia
using EegFun

# Apply 1 Hz high-pass filter (removes frequencies below 1 Hz)
EegFun.highpass_filter!(dat, 1.0)

# More aggressive filtering
EegFun.highpass_filter!(dat, 0.1)  # Very slow drifts only
```

**When to use**: Always apply high-pass filtering (0.1-1 Hz) to remove slow drifts.

## Low-Pass Filtering

Remove high-frequency noise:

```julia
# Apply 40 Hz low-pass filter (removes frequencies above 40 Hz)
EegFun.lowpass_filter!(dat, 40.0)

# For very clean data or sleep studies
EegFun.lowpass_filter!(dat, 30.0)
```

**When to use**: Typically 30-50 Hz to remove muscle artifacts and line noise.

## Band-Pass Filtering

Isolate a specific frequency range:

```julia
# Isolate alpha band (8-12 Hz)
# TODO: Add bandpass_filter! function example when available
# EegFun.bandpass_filter!(dat, 8.0, 12.0)
```

## Notch Filtering

Remove line noise (50 or 60 Hz):

```julia
# TODO: Add notch filter example
# EegFun.notch_filter!(dat, 50.0)  # European line noise
# EegFun.notch_filter!(dat, 60.0)  # North American line noise
```

## Choosing Filter Parameters

> **TODO**: Add guidance on:
> - Filter order selection
> - Cutoff frequency selection
> - Phase distortion considerations
> - Filter type comparison (Butterworth, etc.)

## Verifying Filter Response

> **TODO**: Add example of plotting filter frequency response

## See Also

- [Filtering concepts](../explanations/filtering.md) - Understand filter theory
- [Basic preprocessing tutorial](../tutorials/basic-preprocessing.md) - Complete workflow
