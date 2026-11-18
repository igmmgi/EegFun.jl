# Signal Example — ICA 1: Central Limit Theorem (CLT)

Before jumping into the heavy geometry and optimization algorithms of ICA, we need to answer a fundamental question: **How does an algorithm mathematically know when it has successfully unmixed a signal?**

> [!TIP]
> This demo bridges the gap between the concept of mixtures and the objective function the ICA algorithm uses to untangle them. 

## The Central Limit Theorem in EEG
The Central Limit Theorem (CLT) states that if you take multiple independent variables and add them together, the distribution of their *sum* will always be **more Gaussian** (more like a normal bell-curve) than any of the original variables.

Brain sources (like a blink artifact, a heart beat, or a rigid alpha burst) are usually highly structured, repeating events. Statistically, they are very "pointy" (Super-Gaussian) or very "flat" (Sub-Gaussian). They do not look like perfect bell curves.

However, because a scalp electrode records a volume-conducted mixture of many independent brain, muscle, and eye sources, it sums them together. By the CLT, the resulting scalp EEG distribution becomes much closer to Gaussian.

## The Core Logic of ICA
Because of the Central Limit Theorem, ICA applies this relationship in reverse:
**Mixtures tend toward Gaussianity; underlying sources are typically non-Gaussian.**

To isolate an independent source, the algorithm searches for linear projections whose distributions depart maximally from a Gaussian bell curve. This is why ICA is often described as maximizing non-Gaussianity (or kurtosis/negentropy).

## The Interactive Sandbox
This demonstration illustrates the Central Limit Theorem directly:

- It starts from a single uniform distribution (sub-Gaussian, excess kurtosis = -1.20).
- The slider adds further independent uniform distributions together.
- Summing 2 uniform variables produces a triangular distribution.
- Summing 5 or more rapidly approaches a smooth Gaussian bell curve (excess kurtosis approaching 0.0).

This illustrates the statistical property that makes blind source separation possible.

## Code

```julia
using EegFun
EegFun.signal_example_ica_clt()
```
