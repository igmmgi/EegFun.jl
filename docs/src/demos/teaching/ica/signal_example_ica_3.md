# Signal Example — ICA 3: Blind Source Separation

Interactive demonstration of **Independent Component Analysis (ICA)** — 3 sources, rotation geometry, and scatter plots.

> [!TIP]
> This is **Part 3** of the ICA series. Start with [Part 1: What Is a Mixture?](signal_example_ica_1.md) or [Part 2: Mixing & Unmixing](signal_example_ica_2.md) first.

## What it shows

| Feature | Description |
|---------|-------------|
| True Sources (Row 1) | Three independent signals: S1 (oscillatory bursts, Gabor-shaped), S2 (blink artifacts), S3 (muscle burst) |
| EEG Recordings (Row 2) | The mixed signals as they would appear at electrodes — each is a weighted combination of all sources |
| Recovered Components (Row 3) | The unmixed signals, either via manual rotation sliders or the Auto-ICA (Infomax) algorithm |
| Scatter Plots | Three pairwise joint distributions (S1 vs S2, S1 vs S3, S2 vs S3) — see below |
| Excess Kurtosis | Displayed for sources and recovered components |

### The Scatter Plots

The scatter plots are the most informative panel:

- **Independent, non-Gaussian signals** → distinctive **cross or T-shape** (the signals are statistically independent)
- **Mixed signals** → correlated **blob** (everything jumbled together)
- **Correctly unmixed** → cross shape restored

## Things to Try

1. Set all mixing angles > 0° and watch the scatter plots collapse from crosses into blobs.
2. Drag the unmixing sliders manually to match the mixing angles — the sources reappear.
3. Reset unmixing angles to 0°, then click **Infomax** — ICA automatically recovers the sources.
4. Increase the **Noise** slider — Infomax still works because the blink and muscle sources are highly non-Gaussian.
5. Compare manual rotation with Infomax — Infomax finds a general unmixing matrix, not just a rotation.
6. Try **Extended** Infomax, which handles both sub-Gaussian and super-Gaussian sources.

## Controls

| Control | Range | Description |
|---------|-------|-------------|
| Mix α (S1-S2) | 0–90° | Rotation angle blending S1 and S2 into the recordings |
| Mix β (S1-S3) | 0–90° | Rotation angle blending S1 and S3 into the recordings |
| Mix γ (S2-S3) | 0–90° | Rotation angle blending S2 and S3 into the recordings |
| Unmix φ (S1-S2) | 0–90° | Inverse rotation to recover S1/S2 (manual mode) |
| Unmix ψ (S1-S3) | 0–90° | Inverse rotation to recover S1/S3 (manual mode) |
| Unmix χ (S2-S3) | 0–90° | Inverse rotation to recover S2/S3 (manual mode) |
| Burst Freq | 1–20 Hz | Frequency of the oscillatory source (S1) |
| Noise | 0–1 | Additive Gaussian noise on all sources |
| Infomax | — | Run standard Infomax ICA (Bell & Sejnowski, 1995) |
| Extended | — | Run Extended Infomax (Lee et al., 1999) |

::: details Under the Hood — The Mathematics (click to expand)

### The Mixing Model

With N sources and N sensors, mixing is an N×N matrix multiplication. For 3 sources, we compose three pairwise 2D rotations (α, β, γ):

**Mixing**: R₂₃(γ) · R₁₃(β) · R₁₂(α) · S
**Unmixing**: R₁₂ᵀ(φ) · R₁₃ᵀ(ψ) · R₂₃ᵀ(χ) · M

When φ=α, ψ=β, χ=γ → perfect recovery.

### Why Non-Gaussianity?

ICA exploits the **Central Limit Theorem in reverse**: mixtures of independent signals are *more Gaussian* than the original signals. So finding the unmixing that maximises non-Gaussianity (measured by kurtosis or negentropy) recovers the independent sources.

- **Excess kurtosis = 0**: Gaussian distribution
- **Excess kurtosis > 0**: Super-Gaussian (sharp peaks, heavy tails — like blink spikes)
- **Excess kurtosis < 0**: Sub-Gaussian (flat-topped distributions)

The Infomax algorithm (Bell & Sejnowski, 1995) used in EegFun maximises the total non-Gaussianity of the output components.

:::

## See Also

- [Part 0: Matrix Math Basics](signal_example_ica_0.md) — What is an Unmixing Matrix?
- [Part 1: What Is a Mixture?](signal_example_ica_1.md) — Simplest introduction to the problem
- [Part 2: Mixing & Unmixing](signal_example_ica_2.md) — 2-source introduction to ICA
- [ICA: Central Limit Theorem](signal_example_ica_clt.md) — Deep Dive: Why mixing creates Gaussian blobs
- [Part 4: Sphering (Whitening)](signal_example_ica_4.md) — Why ICA optimization is just finding an angle
- [Part 5: Inside the Black Box](signal_example_ica_5.md) — The optimization landscape and gradient ascent
- [ICA Demo](../preprocessing/ica.md) — Full ICA workflow on real EEG data
- [EEGLAB: Running ICA](https://eeglab.org/tutorials/06_RejectArtifacts/RunICA.html) — Practical guide to ICA in EEGLAB
- [ICA for Dummies](https://arnauddelorme.com/ica_for_dummies/) — Arnaud Delorme's accessible introduction to ICA
- Bell, A. J., & Sejnowski, T. J. (1995). An information-maximization approach to blind separation. *Neural Computation*, *7*(6), 1129–1159.
- Hyvärinen, A., & Oja, E. (2000). Independent component analysis: algorithms and applications. *Neural Networks*, *13*(4-5), 411–430.

## Code

```julia
using EegFun
EegFun.signal_example_ica_3()
```
