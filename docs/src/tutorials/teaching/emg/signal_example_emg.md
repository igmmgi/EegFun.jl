# Surface EMG — Muscle Activation & Processing Cascade Demo

Interactive Surface Electromyography (sEMG) simulator illustrating motor-unit superposition and the canonical 4-stage processing pipeline: **Raw → Rectify → Filter (RMS) → Integrate**.

## What it shows

1. **Raw Surface EMG [μV]**: Zero-mean stochastic interference pattern formed by asynchronous motor-unit action potential trains.
2. **Full-Wave Rectification [μV]**: Absolute value $|V(t)|$ preserving signal energy.
3. **Linear / RMS Envelope [μV]**: Sliding-window smoothing demonstrating the trade-off between temporal resolution and variance.
4. **Integrated EMG [μV·s]**: Cumulative muscle activation across the contraction burst.

## Code

```julia
using EegFun
using GLMakie

EegFun.signal_example_emg()
```
