# EDA — Tonic/Phasic Decomposition, Habituation & Overlap Demo

Interactive multi-trial Electrodermal Activity (EDA) simulator illustrating tonic Skin Conductance Level (SCL) and drift, phasic Event-Related Skin Conductance Responses (ER-SCRs), spontaneous fluctuations (NS-SCRs / SFs), habituation across repeated stimuli, and temporal response superposition.

## What it shows

- **Tonic vs. Phasic Decomposition**: Separating slow baseline SCL and drift from phasic stimulus-locked ER-SCRs and spontaneous fluctuations (SFs).
- **Orienting & Habituation (Barry & Sokolov, 1993)**: Exponential decline of $\sqrt{\text{SCR}}$ across repeated stimuli ($S_1\text{–}S_6$) and dishabituation to a novel stimulus.
- **Temporal Superposition (Overlap)**: How short inter-stimulus intervals (ISIs) cause new SCRs to pile onto the recovery tail of preceding responses, distorting raw trough-to-peak amplitude measurements.

## Code

```julia
using EegFun
using GLMakie

EegFun.signal_example_eda_habituation()
```

