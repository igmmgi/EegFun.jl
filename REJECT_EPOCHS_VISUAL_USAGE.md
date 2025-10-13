

# Visual Epoch Rejection - Usage Guide

## Overview

The `detect_bad_epochs_interactive()` function provides an interactive graphical interface for manually rejecting epochs based on visual inspection. It displays epochs in a grid layout with checkboxes for marking each epoch as good or bad, along with navigation buttons for scrolling through pages of data.

## Why Visual Rejection?

While automatic rejection is fast and objective, visual inspection allows you to:

1. **Catch subtle artifacts**: Some artifacts are obvious to the eye but hard to detect automatically
2. **Make nuanced decisions**: Keep borderline epochs that might be rejected automatically
3. **Quality control**: Verify that automatic rejection worked correctly
4. **Inspect specific patterns**: Look for systematic issues in your data
5. **Learn about your data**: Understand what types of artifacts are present

## Basic Usage

### Quick Start

```julia
using eegfun, JLD2

# Load epoched data
epochs = load("participant_1_epochs.jld2", "epochs")

# Launch interactive rejection interface
state = detect_bad_epochs_interactive(epochs)

# GUI window opens - review epochs:
# ✓ CHECK boxes to mark epochs for REJECTION
# ✗ UNCHECK boxes to KEEP epochs
# Use Next/Previous buttons to navigate

# After reviewing all epochs, extract results:
rejected_indices = get_rejected_epochs(state)
clean_epochs = get_clean_epochs(state)

# Save cleaned data
save("participant_1_epochs_cleaned.jld2", "epochs", clean_epochs)

# Optionally save rejection report
save_rejection_decisions(state, "participant_1_rejection_report.txt")
```

## Interface Layout

The interactive window displays:

### Top Section
- **Title**: "Interactive Epoch Rejection"
- **Info Bar**: Shows current page, rejection count, and reminder that checked = reject

### Main Grid
- **Epoch Plots**: Grid of epoch waveforms (default: 3 rows × 4 columns = 12 epochs per page)
- **Each plot shows**: All selected channels overlaid with different colors
- **Vertical dashed lines**: Mark time = 0 (stimulus/event onset)
- **Horizontal dashed lines**: Mark amplitude = 0

### Checkboxes
- **Below each plot**: Checkbox with label showing epoch number and status
- **Unchecked** (default): "Epoch N: Keep ✓" in green = KEEP this epoch
- **Checked**: "Epoch N: REJECT ✗" in red = REJECT this epoch

### Navigation Buttons (Bottom)
- **|◀ First**: Jump to first page
- **◀ Previous**: Go to previous page
- **Next ▶**: Go to next page
- **Last ▶|**: Jump to last page

## Step-by-Step Workflow

### 1. Launch Interface

```julia
epochs = load("participant_1_epochs.jld2", "epochs")
state = detect_bad_epochs_interactive(epochs)
```

### 2. Review Epochs

For each page:
1. Visually inspect each epoch plot
2. Look for:
   - Large amplitude spikes or jumps
   - Drift or slow waves
   - High frequency noise
   - Flatlines or dropouts
   - Systematic artifacts (blinks, movements)

3. **Check the box** for epochs that look bad
4. **Leave unchecked** for good epochs

### 3. Navigate Through Data

- Use **Next** button to proceed to next page
- Use **Previous** to go back and review
- Use **First**/**Last** to jump around
- The info bar shows progress: "Page X/Y"

### 4. Extract Results

```julia
# After reviewing all epochs:
rejected_indices = get_rejected_epochs(state)
clean_epochs = get_clean_epochs(state)

# Check results
println("Rejected $(length(rejected_indices)) of $(state.n_total_epochs) epochs")
println("Rejection rate: $(round(100 * length(rejected_indices) / state.n_total_epochs, digits=1))%")
```

### 5. Save Results

```julia
# Save cleaned data
save("participant_1_epochs_cleaned.jld2", "epochs", clean_epochs)

# Save rejection report
save_rejection_decisions(state, "participant_1_rejection_report.txt")
```

## Customization Options

### Display Specific Channels

```julia
# Show only frontal channels
state = reject_epochs_interactive(epochs,
    channel_selection = channels(x -> startswith.(string.(x), "F")))

# Show specific channels
state = reject_epochs_interactive(epochs,
    channel_selection = channels([:Fz, :Cz, :Pz]))

# Show all channels except EOG
eeg_only = setdiff(channel_labels(epochs), [:HEOG, :VEOG])
state = reject_epochs_interactive(epochs,
    channel_selection = channels(eeg_only))
```

### Change Grid Layout

```julia
# Show more epochs per page (6 × 3 = 18 epochs)
state = reject_epochs_interactive(epochs,
    epochs_per_page = 18,
    grid_size = (6, 3))

# Show fewer epochs per page with larger plots (2 × 3 = 6 epochs)
state = reject_epochs_interactive(epochs,
    epochs_per_page = 6,
    grid_size = (2, 3))

# Single column for detailed review (6 × 1 = 6 epochs)
state = reject_epochs_interactive(epochs,
    epochs_per_page = 6,
    grid_size = (6, 1))
```

## Complete Example Workflows

### Basic Workflow

```julia
using eegfun, JLD2

# Load epochs
epochs = load("participant_05_epochs.jld2", "epochs")

println("Starting with $(length(epochs.data)) epochs")

# Launch GUI
state = detect_bad_epochs_interactive(epochs)

# ... review epochs visually ...
# ... mark bad ones ...
# ... navigate through all pages ...

# When done, extract results
rejected = get_rejected_epochs(state)
println("Rejected epochs: $rejected")

# Get cleaned data
clean_epochs = get_clean_epochs(state)
println("Remaining: $(length(clean_epochs.data)) epochs")

# Save
save("participant_05_epochs_cleaned.jld2", "epochs", clean_epochs)
save_rejection_decisions(state, "participant_05_manual_rejection.txt")
```

### Combined Automatic + Visual

```julia
using eegfun, JLD2

# Load raw epochs
epochs = load("participant_1_epochs.jld2", "epochs")

# Step 1: Automatic rejection (removes obvious artifacts)
epochs_auto, auto_info = reject_epochs_automatic(epochs, 2.0)
println("After automatic: $(length(epochs_auto.data)) epochs")

# Step 2: Visual review of remaining epochs
state = reject_epochs_interactive(epochs_auto)

# ... review remaining epochs ...

# Step 3: Final cleaned data
final_epochs = get_clean_epochs(state)

# Combined rejection info
total_rejected = auto_info.n_original - length(final_epochs.data)
println("Total rejected: $total_rejected of $(auto_info.n_original)")
println("  Automatic: $(length(auto_info.rejected_epochs))")
println("  Manual: $(sum(state.rejected))")

# Save
save("participant_1_epochs_final.jld2", "epochs", final_epochs)
```

### Review Only Specific Conditions

```julia
# If you have multiple conditions, review each separately
epochs = load("participant_1_epochs.jld2", "epochs")

# Get condition 1 epochs
cond1_indices = [i for (i, epoch) in enumerate(epochs.data) if epoch.condition[1] == 1]
cond1_epochs = EpochData(
    epochs.data[cond1_indices],
    epochs.layout,
    epochs.sample_rate,
    epochs.analysis_info
)

# Review condition 1
state1 = reject_epochs_interactive(cond1_epochs)

# ... review ...

clean_cond1 = get_clean_epochs(state1)

# Repeat for other conditions...
```

## Tips for Effective Review

### 1. Start with Overview

Before checking any boxes:
- Quickly flip through all pages
- Get a sense of overall data quality
- Note any systematic issues

### 2. Look for Patterns

**Good epochs typically show:**
- Smooth waveforms without sharp discontinuities
- Reasonable amplitude range (depends on your data)
- Consistent baseline
- Clear ERP components if expected

**Bad epochs typically show:**
- Sharp spikes or jumps (electrode pops, movements)
- Very large amplitudes (muscle artifacts)
- Slow drifts (skin potentials, sweat)
- High frequency noise throughout
- Flat lines (signal dropout)
- Obvious eye blinks (if not removed by ICA)

### 3. Be Consistent

- Decide on rejection criteria before starting
- Apply same standards throughout
- When in doubt, err on the side of rejection (if you have enough trials)
- Or err on keeping (if trials are scarce)

### 4. Use Multiple Passes

**Pass 1**: Quick review, reject obvious bad epochs
**Pass 2**: Go back and review borderline cases
**Pass 3**: Final check of rejected epochs

### 5. Take Breaks

- Visual inspection can be tiring
- Take breaks every 15-20 minutes
- Come back with fresh eyes if unsure

## Saving Your Progress

The rejection state is stored in memory, so you can:

```julia
# During review, check current status
println(state)
# Shows: current page, rejection count

# Save intermediate progress
using JLD2
save("rejection_state_backup.jld2", "state", state)

# Later, load and continue (would need to recreate GUI)
# This is mainly useful for the rejection vector
```

## Quality Control

### After Review

```julia
# Get summary statistics
rejected = get_rejected_epochs(state)
rejection_rate = 100 * length(rejected) / state.n_total_epochs

println("Rejection Summary:")
println("  Total epochs: $(state.n_total_epochs)")
println("  Rejected: $(length(rejected)) ($(round(rejection_rate, digits=1))%)")
println("  Kept: $(state.n_total_epochs - length(rejected))")

# Check if rejection rate is reasonable
if rejection_rate > 40
    @warn "High rejection rate (>40%). Check data quality or rejection criteria."
elseif rejection_rate < 5
    @warn "Low rejection rate (<5%). Consider if review was thorough enough."
end
```

### Visualize Rejected vs Kept

```julia
# Plot average of rejected vs kept epochs
using CairoMakie

rejected_data = epochs.data[state.rejected]
kept_data = epochs.data[.!state.rejected]

if !isempty(rejected_data) && !isempty(kept_data)
    rejected_epochs_obj = EpochData(rejected_data, epochs.layout, epochs.sample_rate, epochs.analysis_info)
    kept_epochs_obj = EpochData(kept_data, epochs.layout, epochs.sample_rate, epochs.analysis_info)
    
    rejected_erp = average_epochs(rejected_epochs_obj)
    kept_erp = average_epochs(kept_epochs_obj)
    
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel = "Time (s)", ylabel = "μV", title = "Rejected vs Kept Epochs")
    lines!(ax, kept_erp.data.time, kept_erp.data.Cz, label = "Kept", color = :green)
    lines!(ax, rejected_erp.data.time, rejected_erp.data.Cz, label = "Rejected", color = :red)
    axislegend(ax)
    display(fig)
end
```

## Troubleshooting

### Issue: Window not responsive

**Cause**: Makie rendering or event loop issue

**Solution**: Close and restart Julia session, try again

### Issue: Can't see all channels clearly

**Solution**: Adjust channel selection
```julia
# Show fewer channels
state = reject_epochs_interactive(epochs,
    channel_selection = channels([:Fz, :Cz, :Pz]))
```

### Issue: Plots too small

**Solution**: Show fewer epochs per page
```julia
state = reject_epochs_interactive(epochs,
    epochs_per_page = 6,
    grid_size = (2, 3))
```

### Issue: Too many epochs to review

**Solution 1**: Use automatic rejection first
```julia
epochs_auto, _ = reject_epochs_automatic(epochs, 2.0)
state = reject_epochs_interactive(epochs_auto)
```

**Solution 2**: Review subset
```julia
# Review every Nth epoch
subset_indices = collect(1:5:length(epochs.data))
subset_epochs = EpochData(
    epochs.data[subset_indices],
    epochs.layout,
    epochs.sample_rate,
    epochs.analysis_info
)
state = reject_epochs_interactive(subset_epochs)
```

### Issue: Made mistake, want to undo

**Solution**: Checkboxes can be toggled
- Just uncheck a box to un-reject an epoch
- Navigate back to previous pages
- Make changes at any time before extracting final results

## Keyboard Shortcuts

Currently the interface uses mouse/button clicks. Future versions may add:
- Arrow keys for navigation
- Space bar to toggle current epoch
- 'A' to accept all on page
- 'R' to reject all on page

## Best Practices

### 1. Establish Criteria

Before starting, decide what makes an epoch "bad":
- Amplitude threshold (e.g., >100 μV)?
- Specific artifact types?
- Based on specific channels?

### 2. Document Decisions

```julia
# Keep notes on rejection criteria
notes = """
Rejection criteria for Participant 1:
- Epochs with amplitude > 100 μV on any channel
- Obvious eye blinks not caught by ICA
- Baseline drift > 50 μV
- High frequency noise throughout epoch
- Rejection rate: 15%
"""

# Include in report
save_rejection_decisions(state, "participant_1_rejection.txt")
# Add notes to the file manually
```

### 3. Be Systematic

- Review in order (don't skip around randomly)
- Complete each page before moving on
- Do a final review of rejected epochs
- Check summary statistics

### 4. Consider Context

- Time pressure? Use automatic + quick visual review
- Publication quality? Thorough visual review essential
- Many trials? Can afford higher rejection rate
- Few trials? Be more conservative

## Comparison: Automatic vs Visual

| Aspect | Automatic | Visual |
|--------|-----------|--------|
| Speed | Fast (seconds) | Slow (minutes-hours) |
| Objectivity | High | Moderate (subjective) |
| Sensitivity | Misses subtle artifacts | Catches everything visible |
| Specificity | May reject good epochs | Can make nuanced decisions |
| Scalability | Handles any N epochs | Tedious for >200 epochs |
| Reproducibility | Perfect | Variable between reviewers |
| Documentation | Detailed metrics | Need manual notes |

**Recommendation**: Use both!
1. Automatic rejection to remove obvious artifacts
2. Visual review to catch remaining issues

## See Also

- Automatic rejection: `reject_epochs_automatic()`
- Epoch averaging: `average_epochs()`
- ERP plotting: `plot_erp()`
- Epoch plotting: `plot_epochs()`

