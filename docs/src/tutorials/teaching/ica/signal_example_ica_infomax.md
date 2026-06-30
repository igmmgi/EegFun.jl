# Signal Example — ICA 7: Infomax (Information Theory)

The standard algorithm in EEGLAB is **Infomax ICA**, developed by Bell & Sejnowski in 1995.

If you read the original algorithms, you won't immediately see references to "Kurtosis" or the "Central Limit Theorem" (which we used in Parts 4 and 5). Instead, Infomax is anchored purely in **Information Theory**, which mathematically defines Independence.

This demo bridges the gap, proving that hunting for Non-Gaussianity and minimizing Mutual Information are technically hunting the exact same mathematical shapes!

## The Math of Infomax
Information Theory defines statistical independence mathematically: Two variables are completely independent if and only if their **Mutual Information is zero**.

If you look at the slide from your EEGLAB lectures, the equation for Joint Entropy is:
$H(y_1, y_2) = H(y_1) + H(y_2) - I(y_1, y_2)$

Because the "Sphering" process in Part 5 locks our data into a geometrically round matrix, the total Joint Entropy $H(y_1, y_2)$ technically becomes a constant. 
Because of this, trying to drive Mutual Information $I$ to zero is mathematically equivalent to trying to shrink $H(y_1)$ and $H(y_2)$ down to their absolute minimum possible Entropy! 

And in nature, what shape has the maximum possible entropy? A Gaussian Bell curve. 
Therefore, if you look for the minimum possible entropy, you are finding the least-Gaussian shapes! 

## The Sandbox
This demo explicitly calculates the equations from the EEGLAB slide in real-time.
1. **The Outputs**: As you rotate the slider, the data is pushed out into Component 1 and Component 2. We graph their live distributions and calculate their Marginal Entropies.
2. **The Joint Map**: The center plot physically graphs $y_1$ against $y_2$. 
3. **The Proof**: The bottom slide maps out the **Mutual Information** metric $I(y_1, y_2)$. 

As you drag the rotation slider, you will watch the center map untangle itself into a perfectly independent cross-shape. Notice that at that exact moment, the Mutual Information curve mathematically crashes to its absolute floor (Zero bits). 

## Code

```julia
using EegFun
using StatsBase

EegFun.signal_example_ica_infomax()
```
