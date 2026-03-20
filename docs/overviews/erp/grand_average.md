Grand averaging is the process of averaging ERP data across multiple participants. `EegFun.jl` provides two primary ways to calculate grand averages: using in-memory data or processing files directly from disk.

## In-Memory Grand Averaging

If you already have a collection of `ErpData` objects loaded in your session, you can average them directly.

```julia
using EegFun

# Create grand averages for all conditions found in the list
# results = [erp_p1, erp_p2, erp_p3]
grand_avgs = grand_average(results)
```

## Batch Averaging (from Disk)

For large studies, it is often more efficient to process files directly from their storage directory. This approach automatically handles file discovery and condition grouping.

```julia
# Averages all JLD2 files in the folder that match the "erps_cleaned" pattern
grand_average("erps_cleaned", input_dir = "derivatives/erp_analysis/")
```

## How it Works

When `grand_average` is called:
1.  **Grouping**: Data is automatically grouped by condition number across all participants.
2.  **Intersection**: The system finds the intersection of channels available across all participants.
3.  **Averaging**: Channel voltages are averaged point-by-point across the cohort.
4.  **Metadata**: The resulting `ErpData` object tracks the total `n_epochs`.

## Visualizing the Result

Once calculated, grand averages can be visualized using the same tools as individual ERPs.

```julia
# Plot with standard error shading (if n > 1)
plot_erp(ga_results, labels = ["Standard", "Target"])
```
