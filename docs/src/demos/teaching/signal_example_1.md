# Signal Example 1 — Sampling & Reconstruction

Interactive demonstration of the **Nyquist–Shannon sampling theorem**.

Use this demo to show students how sampling rate affects the faithful reconstruction of a sine wave, and how the two main reconstruction methods compare.

## What it shows

| Feature | Description |
|---------|-------------|
| Original signal | A continuous sine wave at the chosen frequency |
| Sampled points | The discrete samples taken at the chosen sampling rate |
| Linear reconstruction | Straight-line interpolation between samples |
| Sinc reconstruction | Whittaker–Shannon (ideal) interpolation |

**Linear reconstruction** connects adjacent samples with straight lines. It is simple but introduces distortion — particularly at low sampling rates relative to the signal frequency.

**Sinc reconstruction** uses the Whittaker–Shannon formula to reconstruct the signal exactly (within the Nyquist limit). This is the method used in real digital signal processing systems.

## Teaching tips

- Start with a sampling rate well above Nyquist (e.g. 10× the signal frequency) and observe that both methods agree with the original signal.
- Gradually reduce the sampling rate toward and then below the Nyquist limit (~2× the signal frequency) to show reconstruction breakdown.
- Enable both reconstructions simultaneously to highlight the superiority of sinc interpolation near the Nyquist limit.

## Controls

| Control | Range | Description |
|---------|-------|-------------|
| Signal Frequency | 1–50 Hz | Frequency of the underlying sine wave |
| Sampling Rate | 2–200 Hz | Number of samples per second |
| ☐ Show Original | — | Toggle the continuous reference signal |
| ☐ Linear | — | Toggle linear interpolation overlay |
| ☐ Sinc | — | Toggle sinc reconstruction overlay |

## See Also

- [Sinc Interpolation for Signal Reconstruction (Wolfram Demonstrations)](https://demonstrations.wolfram.com/SincInterpolationForSignalReconstruction/)
- [Signal Example 2](signal_example_2.md) — multi-signal composition and filtering

## Code

```julia
using EegFun
EegFun.signal_example_1()
```
