# How to Filter EEG Data

This guide shows you how to apply different types of filters to your EEG data.

## High-Pass Filtering

Remove slow drifts and DC offsets:

```julia
using EegFun

# Apply 1 Hz high-pass filter (removes frequencies below 1 Hz)
EegFun.highpass_filter!(dat, 1.0)

# Standard 0.1 Hz filter (removes very slow drifts only)
EegFun.highpass_filter!(dat, 0.1)
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
# Apply high-pass + low-pass to isolate a band
EegFun.highpass_filter!(dat, 8.0)
EegFun.lowpass_filter!(dat, 12.0)
```

> [!NOTE]
> There is no dedicated `bandpass_filter!` — combine a high-pass and low-pass filter instead.

## Choosing Filter Parameters

**Cutoff frequency** — For ERP work, 0.1 Hz high-pass is standard. Use 1 Hz only
for ICA preprocessing. Low-pass at 30-40 Hz removes most muscle noise.

**Filter order** — EegFun defaults to order 1 for high-pass and order 3 for low-pass
(Butterworth IIR). Higher orders give sharper rolloff but increase ringing and
phase distortion. The defaults are safe for most EEG work.

**IIR vs FIR** — IIR (default) is fast and suitable for most use cases. FIR filters
have linear phase but require many more taps; use `filter_method="fir"` if you
need guaranteed zero phase distortion beyond what `filtfilt` provides.

**Zero-phase filtering** — `filtfilt` (default) applies the filter forwards and
backwards, eliminating phase distortion but effectively doubling the filter order.

```julia
# Customise filter parameters
EegFun.highpass_filter!(dat, 0.1; order=2, filter_method="iir")
EegFun.lowpass_filter!(dat, 30.0; order=4, filter_method="fir")
```

## Verifying Filter Response

```julia
# Create a filter and inspect it
filter_info = EegFun.create_highpass_filter(0.1, 256.0)
EegFun.print_filter_characteristics(filter_info)
EegFun.plot_filter_response(filter_info)
```

## See Also

- [API Reference](../reference/index.md) - Complete function documentation
