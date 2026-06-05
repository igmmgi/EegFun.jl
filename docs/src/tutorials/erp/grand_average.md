# Grand Average

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


## Code Examples

::: details Show Code

```julia
# Demo: Grand Average
# Shows how to compute grand averages across participants, both from
# in-memory ERP data and from saved JLD2 files on disk.

# Note: EegFun.example_path() resolves bundled example data paths.
# When using your own data, simply pass the file path directly, e.g.:
# dat = EegFun.read_raw_data("/path/to/your/data.bdf")

using EegFun

const DEMO_OUTPUT = "./tutorials/output/"
mkpath(DEMO_OUTPUT)

#######################################################################
# LOAD MULTI-PARTICIPANT ERP DATA
#######################################################################

# Each file is one participant's averaged ERPs (saved by batch_erp or manually)
erps_p1 = EegFun.read_data(EegFun.example_path("data/julia/erps/example1_erps_good.jld2"))
erps_p2 = EegFun.read_data(EegFun.example_path("data/julia/erps/example2_erps_good.jld2"))
erps_p3 = EegFun.read_data(EegFun.example_path("data/julia/erps/example3_erps_good.jld2"))

# Flatten into a single list — one ErpData per participant per condition
all_erps = vcat(erps_p1, erps_p2, erps_p3)


#######################################################################
# GRAND AVERAGE (IN-MEMORY)
#######################################################################

# Average across participants for every condition found in the list
grand_avgs = EegFun.grand_average(all_erps)

# Only average specific conditions
grand_avgs_12 = EegFun.grand_average(all_erps, condition_selection = EegFun.conditions([1, 2]))


#######################################################################
# VISUALIZE
#######################################################################

# Plot grand average ERPs (SE shading is automatic when n > 1)
EegFun.plot_erp(grand_avgs, channel_selection = EegFun.channels([:Cz, :Pz]))

# Single condition
EegFun.plot_erp(grand_avgs, condition_selection = EegFun.conditions([1]), channel_selection = EegFun.channels(:Pz))


#######################################################################
# BATCH GRAND AVERAGE (FROM DISK)
#######################################################################

# Averages all JLD2 files whose name contains "erps_good" in the given directory.
# Saves the result to tutorials/output/grand_average_erps_good.jld2
EegFun.grand_average("erps_good", input_dir  = EegFun.example_path("data/julia/erps/"), output_dir = DEMO_OUTPUT)

# With participant and condition filtering
EegFun.grand_average(
    "erps_good",
    input_dir = EegFun.example_path("data/julia/erps/"),
    participant_selection = EegFun.participants(1:10),
    condition_selection = EegFun.conditions([1, 2]),
    output_dir = DEMO_OUTPUT,
)


#######################################################################
# LOAD AND PLOT BATCH RESULT
#######################################################################

ga = EegFun.read_data(joinpath(DEMO_OUTPUT, "grand_average_erps_good.jld2"))
EegFun.plot_erp(ga, channel_selection = EegFun.channels([:Cz, :Pz]))
```

:::

## See Also

- [API Reference](../../reference/index.md)
