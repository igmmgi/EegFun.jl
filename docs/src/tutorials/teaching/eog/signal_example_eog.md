# EOG — Corneo-Retinal Dipole & Bipolar Derivation Demo

Interactive Electrooculography (hEOG) simulator illustrating how horizontal eye rotations re-orient the corneo-retinal standing dipole and generate deflections in a bipolar difference amplifier.

## What it shows

- **Corneo-retinal standing potential**: The cornea is positive (+) relative to the retina (−), forming a ~0.4–1.0 mV dipole in each eye.
- **Bipolar derivation**: Right outer canthus ($V_{(+)}$) minus left outer canthus ($V_{(-)}$) cancels common-mode background EEG while doubling opposite-polarity eye-movement potentials.
- **Saccade playback**: Real-time saccade animation (0° → target gaze angle → 0°) with synchronized waveform playhead and optional dipole electric field lines.

## Controls

| Control | Description |
|---------|-------------|
| Gaze Angle | Rotation angle from −45° (left) to +45° (right) |
| Duration | Saccade fixation hold duration (100–1200 ms) |
| Noise | Physiological baseline noise (0–100 μV) |
| Play Saccade | Animates a real-time saccade sweep |
| Field Lines | Toggles bilateral dipole electric field contours |
| Monopolar (CMR) | Overlays unipolar right/left temple traces to illustrate Common-Mode Rejection |

## Code

```julia
using EegFun
using GLMakie

EegFun.signal_example_heog()
```
