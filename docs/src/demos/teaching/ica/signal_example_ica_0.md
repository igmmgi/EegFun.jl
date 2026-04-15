# Signal Example — ICA 0: Matrix Math (U = WX)

Before introducing the geometric complexities of the Cocktail Party problem or trying to optimize statistical functions, we must define exactly what an "Unmixing Matrix" technically is.

> [!TIP]
> This is **Part 0** of the ICA series. It provides the most fundamental mathematical prerequisite for understanding component analysis.

## The Mechanism of $U = WX$

In the Independent Component Analysis equation:
* $\mathbf{X}$ is your block of EEG data, where each row is an electrode.
* $\mathbf{W}$ is the Unmixing Matrix ICA wants to find.
* $\mathbf{U}$ is the resulting independent components.

Many people assume this is an abstract equation, but it is actually a physical machine driven by **Matrix Multiplication**. 

Matrix multiplication works column-by-column across time. At a single timepoint $t$, we can take a "snapshot" of the data as a single column vector. We multiply the **Rows** of the unmixing matrix by that **Column** of data (taking the dot product) to calculate the resulting component outputs for that exact instant in time.

## The Interactive Sandbox
This demonstration physicalizes this exact "Row-by-Column" mechanism using a larger, $3 \times 3$ data architecture (3 EEG channels, 3 extracted components). 

Instead of watching a time-series play out automatically, we have **frozen time**.
Using the **Time Cursor** slider, you manually select a single timepoint $t$.
At the top of the screen, you will see a massive, live Matrix Math execution board:

```text
[ W_11  W_12  W_13 ]       [ X_1(t) ]       [ (W_11*X_1) + (W_12*X_2) + (W_13*X_3) ]       [ U_1(t) ]
                   ×                  =                                            =  
[ W_21  W_22  W_23 ]       [ X_2(t) ]       [ (W_21*X_1) + (W_22*X_2) + (W_23*X_3) ]       [ U_2(t) ]
                   ×                  =                                            =
[ W_31  W_32  W_33 ]       [ X_3(t) ]       [ (W_31*X_1) + (W_32*X_2) + (W_33*X_3) ]       [ U_3(t) ]
```

1. **Watch the inputs**: As you drag the time cursor across the raw EEG waveforms, the 3 values in the $X$ column vector instantly update.
2. **Watch the math**: The center panel actively calculates the cross-multiplication across all 3 variables.
3. **Watch the output**: The resulting $U$ component numbers pop out on the right.
4. **Draw the wave**: Those exact $U$ numbers are structurally plotted as the newest dots on the right-hand graphs. 

By dragging the Time Cursor left-to-right, you are literally executing matrix multiplication frame-by-frame, effectively "drawing" the resulting Component wave $U(t)$ yourself.

## Experiment
You can click and drag the numbers inside the $3 \times 3$ $\mathbf{W}$ matrix to change their weight. Change $W_{11}$ to see how it forces the math board to multiply differently, instantly altering the final drawn component wave.

## Code

```julia
using EegFun
EegFun.signal_example_ica_0()
```
