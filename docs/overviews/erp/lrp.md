This demo shows how to calculate the Lateralised Readiness Potential (LRP) from ERP data.

### What is the LRP?

The LRP is an ERP measure of lateralised motor preparation. It is computed from the difference in activity between hemispheres contralateral and ipsilateral to the responding hand, and is widely used in cognitive and motor neuroscience to study response preparation timing.

### How is it Calculated?

For each channel pair (e.g., C3/C4), the LRP is:

```text
LRP_C3 = 0.5 × ((C3_right − C4_right) + (C4_left − C3_left))
```

Where "left" and "right" refer to the hand used for the response.

### Key Functions

| Function | Purpose |
| --- | --- |
| `lrp(erp_left, erp_right)` | Single-participant LRP |
| `lrp(file_pattern, pairs)` | Batch LRP across participants |

### Channel Pairing

- By default, all odd/even channel pairs are detected automatically (e.g., C3/C4, CP1/CP2)
- Use `channel_selection` to restrict to specific pairs

## Workflow Summary

### 1. Single-Participant LRP

- Calculate LRP from left and right response conditions

### 2. Batch Processing

- Process all participants with specified condition pairs

### 3. Visualisation

- Plot LRP waveforms with `plot_erp`
