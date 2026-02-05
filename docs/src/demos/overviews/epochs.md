This demo demonstrates extracting epochs (segmented trials) from continuous EEG data around experimental events.

### What is Epoching?

Epoching segments continuous EEG data into time-locked trials around specific events (e.g., stimulus presentations, button presses). This is essential for:

- Analyzing event-related brain activity
- Computing ERPs (averaged epochs)
- Trial-based statistical analyses
- Artifact rejection on individual trials

### Epoch Windows

The epoch window defines the time range extracted around each event:

**Pre-stimulus period**:

- Typically -200 to 0 ms before event
- Used for baseline correction
- Captures pre-stimulus activity

**Post-stimulus period**:

- Duration depends on expected response timing
- Example: 0 to 1000 ms for typical ERP analysis
- Longer windows for slow components (e.g., P300, LPP)
- Longer windows for time-frequency analysis (e.g., alpha, theta, gamma)

### Epoch Conditions

Define which triggers correspond to which experimental conditions:

- **Simple conditions**: Single trigger value
- **Sequence conditions**: Multiple triggers in sequence
- **Multiple conditions**: Different trial types (e.g., congruent/incongruent)

### Baseline Correction

Apply baseline correction to epochs:

- **Pre-stimulus baseline**: Most common (-200 to 0 ms)
- **Entire epoch**: For induced responses
- **Post-stimulus**: Specialized analyses

## Workflow Summary

This demo shows epoch extraction and processing:

### 1. Prepare Continuous Data

- Load and preprocess raw data
- Apply filtering and re-referencing
- Apply baseline correction to continuous data

### 2. Define and Extract Epochs

- Define epoch conditions with trigger specifications
- Extract epochs with specified time window
- Separate trials by condition

### 3. Visualize Epochs

- Plot individual epoch timecourses
- Select specific channels for visualization
- Compare conditions

### 4. Baseline and Average

- Apply baseline correction to epochs
- Average epochs into ERPs
- Compare conditions on same plot
