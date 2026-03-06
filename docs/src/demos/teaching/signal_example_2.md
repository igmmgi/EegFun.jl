# Signal Example 2 — Signal Composition & Filtering

Interactive multi-signal composer demonstrating how complex waveforms are built from simple sine waves, how noise affects a signal, and how filters shape the frequency content.

## What it shows

| Row | Plot | Description |
|-----|------|-------------|
| 1–3 | Signal 1–3 | Individual sine waves (y-axes linked for direct amplitude comparison) |
| 4 | Noise | Additive Gaussian noise |
| 5 | Combined Signal | Sum of all signals plus noise, with optional filtered overlay |
| 6 | Frequency Domain | Live power spectrum of the combined signal |

The x-axes of all five time-domain plots are linked — panning or zooming one panel updates all of them simultaneously.

## Teaching tips

- Start with a single signal (Amp > 0 for Signal 1 only) and observe its clean spectral peak.
- Add a second frequency to show how signals superpose in the time domain and produce two peaks in the spectrum.
- Increase noise to show how it raises the spectral floor.
- Enable the LP filter and sweep its cutoff to demonstrate how filtering removes high-frequency content while preserving low-frequency structure.

## Controls

| Control | Range | Description |
|---------|-------|-------------|
| Freq (×3) | 0–80 Hz | Frequency of each sine wave |
| Amp (×3) | 0–10 | Amplitude of each sine wave |
| Phase (×3) | −π to π | Phase offset of each sine wave |
| Noise | 0–2 | Standard deviation of additive Gaussian noise |
| ☐ LP Filter | 0–100 Hz | Low-pass filter cutoff |
| ☐ HP Filter | 0–2 Hz | High-pass filter cutoff |

## See Also

- [Signal Example 1](signal_example_1.md) — Nyquist theorem and signal reconstruction

## Code

```julia
using EegFun
EegFun.signal_example_2()
```
