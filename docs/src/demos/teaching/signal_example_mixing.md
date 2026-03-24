# Signal Example — Mixing & Unmixing (Beginner ICA)

Interactive demonstration of **signal mixing and ICA unmixing**.

## The Cocktail Party Problem

Imagine you're at a party with two people talking. You have two microphones, but each one picks up *both* voices — just at different volumes depending on how close each person is to each mic.

That's exactly what happens with EEG:

- Your **brain** generates an electrical signal
- Your **eyes** produce large signal changes every time you blink
- Each **electrode** on your scalp picks up a *mixture* of both

The result? Every electrode's recording is a jumbled combination of brain activity and blink artifacts. You never see the clean sources directly.

**ICA** (Independent Component Analysis) is the algorithm that untangles the mixture back into the original sources — like being able to isolate each person's voice from the party recording.

## What This Demo Shows

| Column | What you see |
| :--- | :--- |
| **Left — "What's really happening"** | The two true source signals: a brain oscillation and eye blink artifacts |
| **Middle — "What electrodes record"** | The mixed signals, with the mixing equation shown in the title (e.g. `E1 = 1.00 × Brain + 0.50 × Blink`) |
| **Right — "What ICA recovers"** | The unmixed signals recovered by ICA |

Below the plots, three panels show the **Matrix Equations**:

1. **Mixing M**: How sources `S` become electrodes `E` (`M × S = E`).
2. **Exact Inverse M⁻¹**: The mathematically perfect "reverse recipe".
3. **Your Weights / ICA Found**: The current unmixing recipe being applied.

## Things to Try

### Step 1: See the Problem

1. Set **Mix Amount** to 0.5. Look at the middle column — the blink spikes are now mixed into the brain electrode (E1), and brain oscillations contaminate the blink electrode (E2).
2. Look at the **Mixing M** panel. It's a colored 2×2 grid showing exactly how much each source bleeds into each electrode.

### Step 2: Try to Fix It Yourself

1. Use the **"Subtract E2 from E1"** slider. You're trying to remove the blink contamination from electrode 1 by subtracting a scaled version of electrode 2.
2. Watch the **Shape Match** — can you get to 100%? It's surprisingly hard!
3. Notice that even if you get 100% Shape, the **Amplitude Match** might be low. This is because a simple subtraction doesn't correctly invert the whole mixing matrix.

### Step 3: Let ICA Do It

1. Click **Unmix!** — ICA finds the optimal 2×2 matrix automatically.
2. Watch the manual sliders *jump* to the values ICA found — now you can see exactly what weights it used.
3. Compare the **ICA Found Matrix** with the **Exact Inverse M⁻¹** — ICA essentially "reverse engineered" the mixing recipe without ever seeing the original sources!

## Controls

| Control | Range | Description |
| :--- | :--- | :--- |
| Mix Amount | 0–1 | How much the sources bleed across electrodes |
| Brain Freq | 4–20 Hz | Frequency of the oscillatory brain signal |
| Blink Size | 0.5–5 | Amplitude of eye blink artifacts |
| Noise | 0–1 | Additive sensor noise |
| Subtract E2→E1 | 0–1.5 | Manual unmixing: remove E2 contribution from E1 |
| Subtract E1→E2 | 0–1.5 | Manual unmixing: remove E1 contribution from E2 |
| Unmix! | — | Run ICA to find optimal unmixing weights |
| Reset | — | Clear all unmixing |

## Why Does ICA Work?

ICA relies on one key insight: **brain signals and eye blinks are statistically independent** — knowing the state of one tells you nothing about the state of the other.

The Central Limit Theorem says that mixing independent, non-Gaussian sources always produces something *more* Gaussian (more random-looking) than either source alone. ICA exploits this in reverse: it searches for the unmixing matrix that makes the outputs *as non-Gaussian as possible* — which turns out to be equivalent to making them statistically independent.

> [!TIP]
> This demo shows that unmixing is just **Matrix Inversion** — see also the [ICA (Blind Source Separation) demo](signal_example_ica.md).

## See Also

- [ICA (Blind Source Separation)](signal_example_ica.md) — 3-source version with rotation geometry and scatter plots
- [ICA Demo](../preprocessing/ica.md) — Full ICA workflow on real EEG data
- [Signal Example — Composition](../teaching/signal_example_composition.md) — Building complex waveforms from sine waves
- [EEGLAB: Running ICA](https://eeglab.org/tutorials/06_RejectArtifacts/RunICA.html) — Practical guide to ICA in EEGLAB
- [ICA for Dummies](https://arnauddelorme.com/ica_for_dummies/) — Arnaud Delorme's accessible introduction to ICA

## Code

```julia
using EegFun
EegFun.signal_example_mixing()
```
