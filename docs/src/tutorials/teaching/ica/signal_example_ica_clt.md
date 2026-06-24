# Signal Example — Central Limit Theorem (CLT)

Before jumping into the heavy geometry and optimization algorithms of ICA, we need to answer a fundamental question: **How does an algorithm mathematically know when it has successfully unmixed a signal?**

> [!TIP]
> This demo bridges the gap between the concept of mixtures and the objective function the ICA algorithm uses to untangle them. 

## The Central Limit Theorem in EEG
The Central Limit Theorem (CLT) states that if you take multiple independent variables and add them together, the distribution of their *sum* will always be **more Gaussian** (more like a normal bell-curve) than any of the original variables.

Brain sources (like a blink artifact, a heart beat, or a rigid alpha burst) are usually highly structured, repeating events. Statistically, they are very "pointy" (Super-Gaussian) or very "flat" (Sub-Gaussian). They do not look like perfect bell curves.

However, because a scalp electrode records a huge volume-conducted mixture of hundreds of these independent brain, muscle, and eye sources, it effectively sums them all up. Because of the CLT, the resulting EEG signal on your screen looks incredibly Gaussian! 

## The Core Logic of ICA
Because of the Central Limit Theorem, ICA algorithms operate on one beautifully simple reverse rule:
**Mixtures are Gaussian. True sources are not.**

If an algorithm wants to find a pure brain source, all it has to do is computationally rotate the data backwards until the resulting distribution stops looking like a bell curve! This is why you will hear the phrase "ICA maximizes Non-Gaussianity" or "ICA maximizes Kurtosis". 

## The Interactive Sandbox
This demonstration is the purest possible physical proof of the Central Limit Theorem.

- It initializes a single flat distribution (Sub-Gaussian, Kurtosis = -1.20).
- The slider lets you physically add more independent flat distributions together.
- Adding just 2 shapes forms a triangle.
- Adding 5+ shapes instantly forces the histogram perfectly into a mathematically smooth Gaussian Bell Curve (kurtosis dropping to 0.0).

This conclusively proves the statistical behavior that makes Blind Source Separation possible. 

## Code

```julia
using EegFun
EegFun.signal_example_ica_clt()
```
