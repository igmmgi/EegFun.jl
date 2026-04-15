# Signal Example — ICA 4: Sphering (Whitening)

This interactive visualization reveals the hidden mathematical shortcut that makes ICA possible: **Sphering**. 

> [!TIP]
> This is **Part 4** of the ICA series. It bridges the gap between the correlated mixtures of Part 3 and the optimization search space of Part 5.

## The Problem
If you look at raw EEG mixtures, they form skewed, stretched, off-center blobs. If an algorithm blindly tried to untangle that, it would have to calculate huge mathematical matrices involving stretching, skewing, shifting, and rotating all simultaneously. 

That is too difficult for a fast optimization algorithm.

## The Solution: Sphering
Instead of doing everything at once, all modern ICA algorithms apply a standard "Preprocessing" pipeline first called Sphering (or Whitening). 

This demo lets you physically drag the data through these three exact mathematical steps:

1. **Centering**: Subtract the mean value so the data block is physically anchored exactly precisely on `(0,0)`.
2. **PCA Rotation**: Run Principal Component Analysis to find the longest arms of the blob, and physically rotate the blob so it aligns securely with the X and Y axes. This removes any cross-correlation.
3. **Variance Scaling (Sphering)**: Shrink or stretch the X and Y axes depending on their variance until the entire blob is a perfect mathematical sphere.

## Why Do We Do This?
Because once the data is transformed into a perfect sphere, the algorithm's job is over functionally: the only property left to search is **Rotation**. 

By dragging the "Morph" slider all the way to 3.0 (Sphered), you are setting up the perfect circular search space. This directly leads into [Part 5: Inside the Black Box](signal_example_ica_5.md), where we will simply rotate a line through this sphered data to locate the independent components!

## Code

```julia
using EegFun
EegFun.signal_example_ica_4()
```

## See Also
- [Part 3: 3 Sources & Geometry](signal_example_ica_3.md) — Scatter plots and correlation blobs
- [ICA: Central Limit Theorem](signal_example_ica_clt.md) — Why mixed blobs become Gaussian
- [Part 5: Inside the Black Box](signal_example_ica_5.md) — Gradient ascent rotation
