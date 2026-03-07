# Signal Example — Dot Product & Frequency Detection

Interactive demo showing the single core idea behind the DFT: multiplying a signal
by a test sinusoid and summing the result.

![Signal Example Dot Product](/tutorials/signal_example_dotproduct.png)

## What it shows

| Row | Plot | Description |
|-----|------|-------------|
| 1 | Signal + Test sinusoid | The signal (blue) and the unit-amplitude test sinusoid (red) overlaid. |
| 2 | Element-wise product | Signal × test sinusoid, with the dot product value and Freq/Phase match status shown. |

## Key Concept

The DFT answers "how much of frequency *f* is in this signal?" using:

```math
\text{dot product} = \sum_t \, x[t] \cdot \sin(2\pi f t)
```

**Freq match, Phase match** → product is all-positive → large sum → large dot product.

**Freq mismatch** → product alternates +/- → cancels to near zero.

**Phase mismatch (90°)** → product alternates even at matching frequency → dot product ≈ 0.
This is why the real DFT uses a *complex* test sinusoid (sine + cosine) to be phase-independent.

## Controls

| Control | Range | Description |
|---------|-------|-------------|
| Signal Freq | 1–40 Hz | Frequency of the signal |
| Amplitude | 0.1–2 | Signal amplitude — dot product scales proportionally |
| Phase (°) | 0–355° | Signal phase offset — try 90° at matching frequency |
| Test Freq | 0.5–60 Hz | The DFT "probe" frequency |
| Noise | 0–1 | Additive Gaussian noise |

## Things to Try

1. **Frequency match:** set Signal Freq = Test Freq = 10 Hz → product all-positive, dot product ≈ amplitude.
2. **Frequency mismatch:** move Test Freq to 15 Hz → product alternates, dot product ≈ 0.
3. **Phase sensitivity:** set both frequencies to 10 Hz, then drag Phase to 90° → dot product drops to ≈ 0 even though frequencies match. This motivates the complex DFT.
4. **Amplitude scaling:** raise Amplitude to 2 → dot product doubles. The test signal is always ±1; the result reflects only the signal amplitude.
5. **Noise:** add noise → individual products become noisy but the sum (dot product) is still large at the matching frequency (averaging suppresses noise).

## See Also

- [Signal Example (Composition)](signal_example_composition.md) — building signals from sine waves
- [Signal Example (Spectrum)](signal_example_spectrum.md) — the full FFT runs this for every frequency simultaneously
- Cohen, M. X. (2014). *Analyzing Neural Time Series Data*. MIT Press. — Chapter 11

## Code

```julia
using EegFun
EegFun.signal_example_dotproduct()
```
