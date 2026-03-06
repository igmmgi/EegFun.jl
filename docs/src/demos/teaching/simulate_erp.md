# Simulate ERP — Signal Averaging Demo

Interactive ERP simulator for teaching how trial averaging extracts signals from noise.

## What it shows

| Element | Description |
|---------|-------------|
| Grey lines | Individual simulated EEG trials with realistic 1/f background noise |
| Blue line | Trial-average ERP — clarifies as more trials are added |
| Components | Up to 5 independent cosine-shaped ERP components |

**The core insight:** The signal-to-noise ratio of the ERP average scales with √N (where N is the number of trials). Doubling the number of trials improves SNR by ~41%.

## Teaching tips

- Start with 1 trial and a single active component — the ERP equals the single trial.
- Add noise and increase trials to watch the ERP "emerge" from the background activity.
- Activate a second component at a different latency to show how components superpose.
- Introduce **latency jitter** to show how trial-to-trial variability smears and attenuates the average — a key confound in real ERP research.
- Introduce **amplitude jitter** to show how variability in peak amplitude scales the average downward.

## Controls

| Control | Range | Description |
|---------|-------|-------------|
| Number of Trials | 1–500 | Trials to average |
| Noise Amplitude | 0–20 | Background 1/f EEG noise level |
| Active toggle (×5) | on/off | Enable or disable each component |
| Freq (×5) | 0.1–5.0 Hz | Component shape (cosine frequency) |
| Amp (×5) | −10 to 10 μV | Peak amplitude |
| Latency (×5) | 0–1000 ms | Peak latency |
| Amp Jitter (×5) | 0–20 | Trial-to-trial amplitude variability (σ) |
| Lat Jitter (×5) | 0–50 ms | Trial-to-trial latency variability (σ) |

## See Also

- [Signal Example 1 — Nyquist Theorem](signal_example_1.md)
- [Signal Example 2 — Signal Composition](signal_example_2.md)

## Code

```julia
using EegFun
EegFun.simulate_erp()
```
