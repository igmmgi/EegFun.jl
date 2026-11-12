# EDA — Eccrine Sweat Gland Biophysics & SCR Measures Demo

Interactive Electrodermal Activity (EDA) simulator illustrating eccrine sweat gland innervation, Edelberg's (1993) poral valve model, parallel skin conductance ($G_{\text{total}} = \sum G_i = 1 / R_{\text{total}}$), and standard event-related Skin Conductance Response (ER-SCR) measures.

## What it shows

- **Edelberg's Poral Valve Model**: How sympathetic cholinergic bursts fill eccrine sweat ducts, raising intraductal pressure until pores open and skin conductance surges.
- **Conductance ($\mu\text{S}$) vs. Resistance ($\text{k}\Omega$)**: Why parallel sweat glands sum linearly in conductance ($G = 1/R$) while equal resistance drops correspond to vastly different sympathetic responses depending on baseline SCL.
- **Standard ER-SCR Morphology (Dawson et al., 2007)**: Live measurement markers for **Latency** ($1\text{–}3\text{ s}$), **Rise Time** ($1\text{–}3\text{ s}$), **SCR Amplitude** ($\mu\text{S}$), and **50% Half-Recovery Time ($\text{Rec } t/2$)** ($2\text{–}10\text{ s}$).

## Code

```julia
using EegFun
using GLMakie

EegFun.signal_example_eda_scr()
```

