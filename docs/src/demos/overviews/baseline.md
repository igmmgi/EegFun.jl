This demo demonstrates baseline correction methods for ERP data, including different baseline windows and correction approaches for continuous, epoch, and ERP data.

### What is Baseline Correction?

Baseline correction removes the pre-stimulus mean from each trial, ensuring that activity is measured relative to a neutral reference period. This is essential for ERP analysis because:

- **Removes slow drifts**: Eliminates DC offsets that vary across trials
- **Standardizes pre-stimulus activity**: Ensures comparable starting points across conditions
- **Enables amplitude interpretation**: Makes measurements interpretable relative to baseline activity

### Common Baseline Windows

**Pre-stimulus period** (most common):

- Standard: -200 to 0 ms before stimulus onset
- Captures typical pre-stimulus activity
- Assumes stable activity before stimulus

### Baseline Methods

**Mean correction** (default):

- Subtract the mean voltage during baseline period
- Standard approach in ERP research
- Appropriate when baseline is stable

## Best Practices

**Baseline window selection**:

- Choose based on experimental design and paradigm timing
- Standard: -200 to 0 ms for most ERP paradigms
- Document chosen window in methods section
- Keep consistent across conditions and participants

**Timing considerations**:

- Ensure baseline period is free from stimulus overlap
- Avoid including activity from previous trials
- Verify baseline period is stable (visually inspect)

**When to baseline**:

- Apply after filtering and before averaging
- Can apply to epoched or average ERPs

## IMPORTANT

Baseline correction assumes the pre-stimulus period represents neutral brain activity. This assumption may be violated in designs with:

- Very short inter-trial intervals
- Anticipatory activity (e.g., motor preparation)
- Sustained activity from previous trials

## Workflow Summary

This demo demonstrates baseline correction for different data types:

### 1. Baseline Continuous Data 

- Load raw BioSemi data
- Apply baseline to entire continuous recording
- Visualize DC offset removal in databrowser
- Baseline to specific timepoint (t=0)

### 2. Baseline Epoch Data

- Extract epochs around experimental events
- Apply baseline at different timepoints (t=0, t=0.5)
- Visualize effect on individual trials
- Compare baseline choices

### 3. Baseline ERP Data

- Average epochs into ERPs
- Baseline to single timepoint (t=0, t=0.5)
- Baseline to time window (-200 to 0 ms)
- Compare effects across baseline choices
