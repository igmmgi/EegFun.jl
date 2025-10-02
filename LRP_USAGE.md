# Lateralized Readiness Potential (LRP) - Usage Guide

## Overview

The `lrp()` function calculates the lateralized readiness potential from two ERP datasets representing left-hand and right-hand responses. This measure isolates lateralized motor preparation activity by comparing contralateral vs ipsilateral activation.

## Basic Usage

### Single pair of conditions

```julia
using eegfun

# Load your data (e.g., from JLD2 files)
using JLD2
erps = load("participant_1_erps.jld2", "erps")

# If erps is a vector where:
# - erps[1] = left-hand responses
# - erps[2] = right-hand responses

lrp_result = lrp(erps[1], erps[2])
```

### Multiple condition pairs (NEW!)

**This is probably what you want if you have many conditions!**

```julia
# Load your data with 16 conditions (odd=left, even=right)
erps = load("participant_1_erps.jld2", "erps")

# Define all condition pairs automatically
pairs = [(i, i+1) for i in 1:2:15]  # Creates [(1,2), (3,4), (5,6), ..., (15,16)]

# Calculate all LRPs at once
lrp_results = lrp(erps, pairs)

# Now lrp_results is a vector where:
# lrp_results[1] = LRP for conditions 1 vs 2
# lrp_results[2] = LRP for conditions 3 vs 4
# lrp_results[3] = LRP for conditions 5 vs 6
# ... etc.

# Save all results
save("participant_1_lrp.jld2", "lrp", lrp_results)
```

### If you have separate variables

```julia
# Load conditions separately
erp_left_hand = load("participant_1_condition_1.jld2", "erp")
erp_right_hand = load("participant_1_condition_2.jld2", "erp")

# Calculate LRP
lrp_result = lrp(erp_left_hand, erp_right_hand)
```

## What do "left hand" and "right hand" mean?

The terminology refers to **which hand the participant used to respond**:

- **First argument (erp_left)**: ERP data when participant responded with LEFT hand
- **Second argument (erp_right)**: ERP data when participant responded with RIGHT hand

The LRP calculation assumes:
- For left-hand responses: RIGHT hemisphere (e.g., C4) is more active (contralateral)
- For right-hand responses: LEFT hemisphere (e.g., C3) is more active (contralateral)

## Channel Selection

By default, the function automatically detects all lateral channel pairs based on odd/even numbering:

```julia
# Automatically detects pairs like: C3/C4, C1/C2, Fp1/Fp2, CP3/CP4, etc.
lrp_result = lrp(erp_left, erp_right)
```

You can select specific channels using the channel predicate system (specify left/odd channels only):

```julia
# Only calculate LRP for C3/C4 and CP3/CP4
lrp_result = lrp(erp_left, erp_right, 
                 channel_selection = channels([:C3, :CP3]))

# Use pattern matching for all C-channels (C1/C2, C3/C4, C5/C6, etc.)
lrp_result = lrp(erp_left, erp_right,
                 channel_selection = channels(x -> startswith.(string.(x), "C")))

# Select channels by index
lrp_result = lrp(erp_left, erp_right,
                 channel_selection = channels([1, 3, 5]))  # First 3 odd-numbered channels
```

**Important:** You only need to specify the **left/odd hemisphere channels** (e.g., C3, C1, Fp1). 
The function automatically pairs them with their right/even counterparts (C4, C2, Fp2).

## Batch Processing from Files (NEW!)

**The easiest way to process multiple participants!**

### Basic Batch Processing

```julia
using eegfun

# Process all participants in the current directory
# with 16 conditions (odd=left, even=right)
pairs = [(i, i+1) for i in 1:2:15]
lrp("erps_cleaned", pairs)

# This automatically:
# - Finds all files matching "erps_cleaned" 
# - Calculates LRP for all 8 pairs
# - Saves to new directory "lrp_erps_cleaned_1-2_3-4_5-6_..."
# - Creates a log file
```

### Batch Processing with Options

```julia
# Specific participants only
lrp("erps_cleaned", [(1,2), (3,4)], participants=[1, 2, 3])

# Only C3/C4 and CP3/CP4 pairs
lrp("erps_cleaned", [(1,2), (3,4)], 
    channel_selection=channels([:C3, :CP3]))

# Custom input/output directories
lrp("erps_cleaned", [(1,2)], 
    input_dir="/data/study1/erps",
    output_dir="/data/study1/lrp")

# Full example: 20 participants, all C-channels
pairs = [(i, i+1) for i in 1:2:15]
lrp("erps_cleaned", pairs,
    input_dir="/data/study1",
    participants=1:20,
    channel_selection=channels(x -> startswith.(string.(x), "C")))
```

### What Gets Created

After running batch processing:

```
your_data_directory/
├── 1_erps_cleaned.jld2
├── 2_erps_cleaned.jld2
├── 3_erps_cleaned.jld2
└── lrp_1-2_3-4_5-6_7-8_9-10_11-12_13-14_15-16/  ← NEW!
    ├── 1_erps_cleaned.jld2     ← Contains "lrp" variable
    ├── 2_erps_cleaned.jld2     ← Contains "lrp" variable  
    ├── 3_erps_cleaned.jld2     ← Contains "lrp" variable
    └── lrp.log                 ← Processing log
```

### Loading Batch Results

```julia
using JLD2

# Load LRP results for one participant
lrp_data = load("lrp_1-2_3-4_5-6_7-8_9-10_11-12_13-14_15-16/1_erps_cleaned.jld2", "lrp")

# lrp_data is a Vector{ErpData}:
# lrp_data[1] = LRP for conditions 1 vs 2
# lrp_data[2] = LRP for conditions 3 vs 4
# etc.
```

## Complete Example Workflow

### Single Participant (In-Memory)

```julia
using eegfun
using JLD2

# Load ERP data for one participant
erps = load("participant_05_erps_cleaned.jld2", "erps")

# Assuming conditions are organized as:
# erps[1] = left-hand responses
# erps[2] = right-hand responses

# Calculate LRP
lrp_data = lrp(erps[1], erps[2])

# The result contains LRP values for each lateral channel pair
println("LRP channels: ", channel_labels(lrp_data))
# Output: [:C3, :C4, :C1, :C2, :Fp1, :Fp2] (if these pairs exist in your data)

# Save the result
save("participant_05_lrp.jld2", "lrp", lrp_data)
```

### Multiple Participants - Manual Loop (Alternative)

If you prefer more control, you can manually loop through participants:

```julia
using eegfun
using JLD2

participants = 1:20
condition_pairs = [(i, i+1) for i in 1:2:15]

for participant in participants
    filename = "$(participant)_erps_cleaned.jld2"
    
    if !isfile(filename)
        @warn "Skipping participant $participant (file not found)"
        continue
    end
    
    erps = load(filename, "erps")
    lrp_results = lrp(erps, condition_pairs)
    save("$(participant)_lrp.jld2", "lrp", lrp_results)
    
    println("✓ Processed participant $participant")
end
```

**Note:** Batch processing with `lrp("erps_cleaned", pairs)` is usually easier!

### Multiple Participants (Simple Case - 2 conditions)

```julia
using eegfun
using JLD2

participants = 1:20

for participant in participants
    # Load ERP data
    filename = "participant_$(participant)_erps_cleaned.jld2"
    
    if !isfile(filename)
        println("Skipping participant $participant (file not found)")
        continue
    end
    
    erps = load(filename, "erps")
    
    # Calculate LRP (assuming conditions 1 and 2 are left/right hand)
    lrp_data = lrp(erps[1], erps[2])
    
    # Save result
    save("participant_$(participant)_lrp.jld2", "lrp", lrp_data)
    
    println("✓ Processed participant $participant")
end
```

### With Different Condition Numbering

If your left/right hand conditions are not the first two conditions:

```julia
# Example: conditions 3 and 4 are your left/right hand responses
lrp_data = lrp(erps[3], erps[4])

# Or if they're not consecutive:
# Condition 2 = left hand, Condition 5 = right hand
lrp_data = lrp(erps[2], erps[5])
```

## Interpreting the Results

The LRP result is an `ErpData` object containing lateralized activity:

```julia
lrp_data = lrp(erp_left, erp_right)

# Positive values indicate greater contralateral activation
# (the expected pattern for motor preparation)

# You can plot the LRP
using CairoMakie
plot_erp(lrp_data, channels = [:C3, :C4])

# Or analyze specific time windows
# Find peak between 100-300 ms
time_window = (lrp_data.data.time .>= 0.1) .& (lrp_data.data.time .<= 0.3)
peak_lrp = maximum(lrp_data.data.C3[time_window])
```

## Common Patterns

### Grand Average LRP

```julia
# Calculate LRP for each participant first
lrp_data_all = []
for participant in 1:20
    erps = load("participant_$(participant)_erps.jld2", "erps")
    lrp_data = lrp(erps[1], erps[2])
    push!(lrp_data_all, lrp_data)
end

# Then create grand average
grand_avg_lrp = grandaverage(lrp_data_all)
```

### Comparing Different Conditions

```julia
# If you have two task conditions with left/right hand responses in each:
# Task A: conditions 1 (left) and 2 (right)
# Task B: conditions 3 (left) and 4 (right)

lrp_task_a = lrp(erps[1], erps[2])
lrp_task_b = lrp(erps[3], erps[4])

# Compare LRPs between tasks
# (both are now ErpData objects with the same channels)
```

## Troubleshooting

### Error: "No lateral channel pairs detected"

This means your data doesn't have the expected odd/even numbered electrode pairs. Specify them manually:

```julia
lrp_data = lrp(erp_left, erp_right, 
               channel_pairs = [(:C3, :C4)])
```

### Error: "Sample rates differ"

Both ERP datasets must have the same sample rate. Check your preprocessing pipeline.

### Error: "Number of time points differ"

Both ERP datasets must have the same time window. Make sure you epoched them with the same parameters.

### Wrong order of arguments?

If your LRP values are negative when you expect positive (or vice versa), you may have swapped the left/right conditions:

```julia
# If you accidentally did:
lrp_data = lrp(erp_right, erp_left)  # WRONG ORDER!

# Just swap them:
lrp_data = lrp(erp_left, erp_right)  # CORRECT ORDER
```

## Formula Reference

For each lateral channel pair (e.g., C3/C4), the LRP is calculated as:

```
LRP_C3 = 0.5 × ((C3_right - C4_right) + (C4_left - C3_left))
LRP_C4 = 0.5 × ((C4_right - C3_right) + (C3_left - C4_left))
```

Where:
- `C3_left`, `C4_left`: Activity at C3/C4 during left-hand responses
- `C3_right`, `C4_right`: Activity at C3/C4 during right-hand responses

This double-subtraction formula isolates lateralized activity while canceling non-lateralized components.

## References

- Coles, M. G. H. (1989). Modern mind-brain reading: Psychophysiology, physiology, and cognition. *Psychophysiology*, 26(3), 251-269.
- de Jong, R., Wierda, M., Mulder, G., & Mulder, L. J. (1988). Use of partial stimulus information in response processing. *Journal of Experimental Psychology: Human Perception and Performance*, 14(4), 682-692.
- Oostenveld, R., Stegeman, D. F., Praamstra, P., & van Oosterom, A. (2003). Brain symmetry and topographic analysis of lateralized event-related potentials. *Clinical Neurophysiology*, 114(7), 1194-1202.

