# ECG — Einthoven's Triangle & Vectorcardiogram Demo

Interactive Electrocardiography (ECG) simulator illustrating cardiac dipole propagation (P, Q, R, S, T waves), Einthoven's triangle, and bipolar limb leads (Lead I, Lead II, Lead III).

## What it shows

- **Cardiac dipole projection**: How a 2D/3D depolarization vector projects onto Lead I (0°), Lead II (+60°), and Lead III (+120°).
- **Mean electrical axis rotation**: Real-time waveform morphology changes across leads as the heart axis shifts from −45° to +120°.
- **Vectorcardiogram loop**: Complete P, QRS, and T loops synchronized with the ECG strip.

## Code

```julia
using EegFun
using GLMakie

EegFun.signal_example_ecg_einthoven()
```
