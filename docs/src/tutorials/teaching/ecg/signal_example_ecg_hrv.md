# ECG — Heart Rate Variability (HRV) & Respiratory Sinus Arrhythmia Demo

Interactive Heart Rate Variability (HRV) simulator demonstrating how respiration modulates beat-to-beat R–R intervals under relaxed (high HRV) versus stressed (low HRV) autonomic states.

## What it shows

- **Respiratory Sinus Arrhythmia (RSA)**: Heart rate accelerates during inhalation and decelerates during exhalation.
- **Autonomic flexibility**: Comparing large, flexible R–R oscillations at rest against rigid, low-variability rhythms under sympathetic stress.

## Code

```julia
using EegFun
using GLMakie

EegFun.signal_example_ecg_hrv()
```
