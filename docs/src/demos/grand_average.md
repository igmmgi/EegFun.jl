# Grand Average

Grand averaging is the process of averaging ERP data across multiple participants to isolate reliable experimental effects and reduce individual noise. `EegFun.jl` provides two primary ways to calculate grand averages: using in-memory data or processing files directly from disk.

## 1. In-Memory Grand Averaging

If you already have a collection of `ErpData` objects loaded in your session, you can average them directly.

```julia
using EegFun

# Create grand averages for all conditions found in the list
# results = [erp_p1, erp_p2, erp_p3]
grand_avgs = grand_average(results)
```

## 2. Batch Averaging (from Disk)

For large studies, it is often more efficient to process files directly from their storage directory. This approach automatically handles file discovery and condition grouping.

```julia
# Averages all JLD2 files in the folder that match the "erps_cleaned" pattern
grand_average("erps_cleaned", input_dir = "derivatives/erp_analysis/")
```

## 3. How it Works

When `grand_average` is called:
1.  **Grouping**: Data is automatically grouped by condition number across all participants.
2.  **Intersection**: The system finds the intersection of channels available across all participants.
3.  **Averaging**: Channel voltages are averaged point-by-point across the cohort.
4.  **Metadata**: The resulting `ErpData` object tracks the total `n_epochs`.

## 4. Visualizing the Result

Once calculated, grand averages can be visualized using the same tools as individual ERPs.

```julia
# Plot with standard error shading (if n > 1)
plot_erp(ga_results, labels = ["Standard", "Target"])
```

---

> [!TIP]
> Always ensure your individual ERPs are baseline-corrected and filtered before grand averaging.


## Code Examples

::: details Show Code

```julia
# Demo: Grand Average
# Shows how to create grand averages across multiple participants from disk.

using EegFun

# Suppose we have ERPs from three participants
# results = [erp_p1, erp_p2, erp_p3]
# For this demo, let's assume 'results' is already defined

# Create grand averages for all conditions found in the list
# grand_avgs = grand_average(results)

# --- Batch Averaging (from Disk) ---

# Averages all JLD2 files in the folder that match the "erps_cleaned" pattern
# This creates a new directory "grand_average_erps_cleaned" containing the results
# grand_average("erps_cleaned", input_dir = "derivatives/erp_analysis/")

# Filtering Participants and Conditions
# grand_average("erps_cleaned",
#     participant_selection = participants(1:10),  # Only first 10 participants
#     condition_selection = conditions([1, 2]),    # Only "Standard" and "Target"
#     output_dir = "results/grand_averages"
# )

# --- Visualizing the Result ---

# Load the batch results
# ga_results = read_data("results/grand_averages/grand_average_erps_cleaned.jld2")

# Plot with standard error shading (if n > 1)
# plot_erp(ga_results, labels = ["Standard", "Target"])
```

:::

## See Also

- [API Reference](../reference/index.md)
