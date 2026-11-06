This demo shows how to use the interactive filter GUI to explore the effect of different filter settings on ERP waveforms.

### Why an Interactive Filter Tool?

- **Teaching** — show students how filters affect ERP signals in real time
- **Parameter selection** — experiment with cutoff frequency and filter order before committing
- **Method comparison** — toggle between IIR (Butterworth) and FIR filters, and compare two-pass `filtfilt`, MNE/EEGLAB-style one-pass zero-phase FIR (`filt_zerophase`), and causal one-pass `filt`

### Key Features

| Feature | Description |
| --- | --- |
| Lowpass / Highpass toggles | Enable each filter independently |
| Cutoff sliders | Adjust frequency in real time |
| Order / FIR tap readout | Adjust IIR filter steepness or inspect live FIR tap count and group delay |
| Method & Function menus | Compare Butterworth IIR vs FIR Hamming across `filtfilt`, `filt_zerophase`, and `filt` |
| Mirror toggle | Enable mirror padding to reduce edge artifacts |
| Multi-condition | Compare filter effects across conditions side-by-side |

### What You'll Learn

1. Launching the GUI for a single ERP
2. Focusing on a specific channel
3. Comparing filter effects across multiple conditions
