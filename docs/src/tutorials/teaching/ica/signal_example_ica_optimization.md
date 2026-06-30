# Signal Example — ICA 6: Inside the Black Box

This interactive visualization demystifies the ICA algorithm by showing how it physically searches for independent sources by treating it as an optimization ("hill climbing") problem.

> [!TIP]
> This is Part 6 of the ICA series. Before attempting this, ensure you understand why the search space is purely rotational by looking at **Part 5: Sphering**.

## What is the Black Box doing?

Earlier demos showed that ICA "magically" finds an unmixing matrix to recover the sources. Mathematically, it does this by making the projected data as **non-Gaussian** as possible.

The "Central Limit Theorem in reverse" tells us:
- Mixed signals are more Gaussian (like a standard bell curve).
- The original independent sources are non-Gaussian (often "peaky" or super-Gaussian, having high Kurtosis).

So to find the sources, ICA just rotates its unmixing weights, checks the "peakiness" (Kurtosis or Entropy), and keeps turning in the direction where it gets more peaky.

## The Search Space (The Landscape)

This demo removes time entirely and just looks at the distribution of the mixed values:
1. **The Search Space (Sphered Data)**: Instead of plotting $x_1(t)$ against $t$, we plot $x_1(t)$ against $x_2(t)$. You'll notice the data is "Sphered" (or whitened). Sphering is a preprocessing step that forces the mixture blob to be perfectly round (variance of 1 in all directions). By perfectly rounding out the data first, ICA's job is massively simplified: it only has to **rotate** the line to find the sources, rather than stretching and squishing it as well. The red line $\mathbf{w}$ is our rotatable unmixing vector.
2. **The Resulting 1D Data (Top Right)**: If we collapse the 2D blob down onto the red line, what is the shape of the histogram? When the line cuts through an "arm", the histogram becomes very peaky.
3. **The Objective Landscape (Bottom)**: We can calculate the "Peakiness" (Kurtosis) for *every possible angle* of $\mathbf{w}$ from 0° to 180°. This creates a hill-landscape curve. The red dot represents our current angle.

## Things to Try

1. **Manual Search**: Drag the **Angle slider** manually. Watch the red dot move across the landscape. Notice how the histogram peaks strongly precisely when the red dot hits the top of a hill.
2. **Single Step**: Set the angle manually to 15°. Now click **Step Algorithm**. Watch the red line visibly rotate and "climb" up the slope of the hill. This is directly executing the gradient ascent learning rule!
3. **Convergence**: Click **Run to Convergence** to watch the algorithm automatically race to the peak (the isolated source) and stop when the slope is zero.

## Code

```julia
using EegFun
EegFun.signal_example_ica_optimization()
```
