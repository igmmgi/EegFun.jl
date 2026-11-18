# Signal Example — ICA 5: Sphering (Whitening)

This interactive visualization illustrates a key preprocessing step in ICA: **sphering (whitening)**.

> [!TIP]
> This is Part 5 of the ICA series. It connects the correlated mixtures of Part 4 to the rotational optimization space of Part 6.

## The Problem
Raw EEG mixtures form correlated, unequally scaled, off-center distributions. Estimating an unconstrained unmixing matrix directly requires solving for centering, scaling, decorrelation, and rotation all at once.

## The Solution: Sphering
Instead of estimating everything simultaneously, ICA algorithms first apply **sphering** (whitening).

This demo steps the data through the three stages of sphering:

1. **Centering**: Subtract the mean of each channel so the distribution is centered at `(0, 0)`.
2. **PCA Rotation**: Rotate the principal axes of covariance to align with the coordinate axes, removing pairwise correlation.
3. **Variance Scaling (Sphering)**: Rescale each axis to unit variance so the covariance matrix becomes the identity matrix.

## Why Do We Do This?
Once the data are sphered (uncorrelated with unit variance in every direction), any remaining linear unmixing transformation is an **orthogonal rotation**.

Dragging the **Morph** slider to `3.0` (Sphered) produces this isotropic search space, which leads directly into [Part 6: Inside the Black Box](signal_example_ica_optimization.md), where rotating a projection vector through the sphered data locates the independent components.

## Code

```julia
using EegFun
EegFun.signal_example_ica_sphering()
```

## See Also
- [Part 4: 3 Sources & Geometry](signal_example_ica_geometry.md) — Scatter plots and correlation blobs
- [ICA: Central Limit Theorem](signal_example_ica_clt.md) — Why mixed blobs become Gaussian
- [Part 6: Inside the Black Box](signal_example_ica_optimization.md) — Gradient ascent rotation
