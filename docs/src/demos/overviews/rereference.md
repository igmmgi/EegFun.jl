This demo demonstrates applying different re-referencing schemes to EEG data.

### What is Re-referencing?

Re-referencing changes the voltage baseline by subtracting a reference signal from all electrodes. Since EEG measures potential differences, the choice of reference affects:

- Signal amplitudes and topographies
- Component interpretations
- Statistical comparisons

### Why Re-reference?

**Recording reference**: The reference used during data acquisition (e.g., Cz, linked mastoids, or manufacturer default)

**Analysis reference**: The reference appropriate for your research question and field conventions

Re-referencing allows you to change from the recording reference to a more appropriate analysis reference.

### Common Reference Schemes

**Average Reference (`:avg`)**:

- Subtracts mean of all electrodes from each channel
- Assumes zero average scalp potential
- Most common for high-density EEG (64+ channels)
- Reference-free topographies

**Mastoid Reference**:

- References to average of left/right mastoid electrodes (e.g., `:M1`, `:M2`)
- Common in ERP research
- Reduces contamination from reference site activity

**Single Electrode Reference**:

- References to a specific electrode (e.g., `:Cz`, `:Fz`)
- Useful for specific research questions
- Activity at reference site becomes invisible (single reference electrode = 0)

### Best Practices

- Apply the same reference across all conditions and participants
- Re-reference **before** interpolating bad channels (bad channels affect average reference)
- Document your reference choice in publications
- Use average reference for high-density recordings

### Re-referencing is Reversible

You can change references multiple times. EegFun tracks the current reference and updates accordingly.

## Workflow Summary

This demo shows re-referencing workflows:

### 1. Continuous Data

- Load and preprocess data
- Apply different reference schemes (`:Fp1`, `:AF7`)
- Visualize in databrowser to see effects

### 2. Epoched Data

- Create epochs from continuous data
- Apply references to segmented trials
- Data structure automatically updated

### 3. ERP Data

- Average epochs into ERPs
- Apply references to averaged data
- Reference changes preserved across data types
