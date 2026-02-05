This demo demonstrates visualizing event markers and triggers in continuous EEG data to verify timing and event sequences.

### What are Triggers?

Triggers (also called event markers or stimulus codes) are time-stamped codes that mark when experimental events occurred during recording:

- Stimulus presentations
- Participant responses
- Experimental conditions
- Trial boundaries
- Hardware events
- etc.

### Trigger Cleaning

EegFun automatically cleans triggers by default removing consecutive duplicates. For example, the raw sequence `0 0 1 1 0 0 2 2 2 0` becomes `0 0 1 0 0 0 2 0 0 0`. This ensures each trigger represents a single event rather than a sustained hardware signal.

### Trigger Visualization Functions

**trigger_count**:

- Summary statistics of all trigger codes
- Counts of each trigger type
- Identifies missing or unexpected triggers

**plot_trigger_overview**:

- Visual representation of trigger occurrences
- Color-coded by trigger type
- Shows distribution across recording

**plot_trigger_timing**:

- Inter-trigger intervals (ITIs)
- Timing precision verification

### Use Cases

**Quality control**:

- Verify triggers were recorded correctly
- Confirm expected trigger counts
- Identify missing or duplicate triggers

**Timing analysis**:

- Check inter-stimulus intervals
- Verify experimental timing

**Troubleshooting**:

- Identify spurious triggers
- Find timing drift or jitter

### Filtering Triggers

Use `ignore_triggers` to exclude specific codes:

- Filter out hardware markers
- Remove boundary codes
- Focus on experimental events only

## Workflow Summary

This demo shows trigger visualization workflows:

### 1. Count Triggers

- Load raw data
- Count triggers before processing
- Verify expected trigger codes exist

### 2. Create Data Structure

- Load layout and create EegFun structure
- Count triggers again to verify preservation

### 3. Visualize Overview

- Plot trigger distribution
- Optionally ignore certain trigger codes
- Assess trigger patterns

### 4. Analyze Timing

- Plot inter-trigger intervals
- Verify timing consistency
- Identify timing issues
