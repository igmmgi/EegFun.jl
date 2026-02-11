This demo shows how to create jackknife (leave-one-out) averages for robust statistical testing of ERP latency measures.

### What is Jackknife Averaging?

The jackknife technique creates N grand averages from N participants, where each average excludes one participant. This produces smoother waveforms with clearly defined peaks, enabling reliable measurement of onset or peak latencies that would be noisy in individual-participant data.

### Why Use It?

Standard ERP peak latency measures are unreliable in noisy single-participant waveforms. The jackknife method:

- Produces cleaner waveforms by averaging many participants
- Enables latency measurement on smooth grand-average-like waveforms
- Requires a simple correction to the t-statistic: t_corrected = t / (n − 1)

### Key Functions

| Function | Purpose |
| --- | --- |
| `jackknife_average(erps)` | Single-condition jackknife from a vector of ERPs |
| `jackknife_average(file_pattern)` | Batch jackknife across participant files |

### Reference

Miller, Patterson, & Ulrich (1998). Jackknife-based method for measuring LRP onset latency differences. *Psychophysiology*, 35, 99–115.

## Workflow Summary

### 1. Single-Condition Jackknife

- Create leave-one-out averages from a vector of participant ERPs

### 2. Batch Processing

- Process all participant files and save jackknife averages

### 3. Typical Pipeline

- Calculate LRP → Jackknife average → Measure onset latency → Apply correction to t-test
