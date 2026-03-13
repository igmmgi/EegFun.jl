# Batch Processing

This demo shows how to load and process multiple EEG files, manage file collections, and group results by experimental condition.

### Batch Loading

EegFun provides utilities for working with multiple data files:

- **`read_all_data`**: Load all matching files from a directory using pattern matching
- **`read_data`**: Load a single saved data file (JLD2 format)
- **`get_files`**: Find files matching a pattern in a directory

### Grouping by Condition

After loading data from multiple participants, `group_by_condition` organises the results into an `OrderedDict` keyed by condition number. This is essential for:

- Computing grand averages per condition
- Running group-level statistics
- Comparing conditions across participants

### File Management

Helper functions for working with data files:

- **`check_files_exist`**: Verify a list of file paths exist
- **`get_files`**: Find files matching a regex pattern
- **`find_file`**: Recursively search for a file by name

## Workflow Summary

This demo covers:

- Loading multiple saved data files from a directory
- Grouping loaded data by condition number
- Computing grand averages from batch data
- File existence checking and pattern-based file finding


## Code Examples

::: details Show Code

```julia
# Demo: Batch Processing
# Shows how to load multi-participant ERP data, group by condition,
# and compute grand averages.

using EegFun
# Note: EegFun.example_path() resolves bundled example data paths.
# When using your own data, simply pass the file path directly, e.g.:
# dat = EegFun.read_raw_data("/path/to/your/data.bdf")

const ERP_DIR = EegFun.example_path("data/julia/erps/")


#######################################################################
# LOAD PROCESSED DATA
#######################################################################

# Load a single participant's ERPs
erps_p1 = EegFun.read_data(joinpath(ERP_DIR, "example1_erps_good.jld2"))

# Load all 12 participants from the ERP directory
erps_all = EegFun.read_all_data(joinpath(ERP_DIR, "erps_good"))

# Subset to first 6 participants
erps_subset = EegFun.read_all_data(joinpath(ERP_DIR, "erps_good"), EegFun.participants(1:6))


#######################################################################
# GRAND AVERAGE
#######################################################################

# Average across all participants, one ErpData per condition
grand_avgs = EegFun.grand_average(erps_all)

# Specific conditions only
grand_avgs_12 = EegFun.grand_average(erps_all, condition_selection = EegFun.conditions([1, 2]))

# From disk directly — saves to output_dir automatically
# EegFun.grand_average("erps_good", input_dir = ERP_DIR, output_dir = EegFun.example_path("data/julia/grand_average/"))


#######################################################################
# FILE UTILITIES
#######################################################################

# List all .bdf files in a directory (regex pattern)
bdf_files = EegFun.get_files(EegFun.example_path("data/bdf"), "\\.bdf")

# Verify all files exist before starting a long pipeline run
EegFun.check_files_exist(bdf_files)

# Find a single file by searching recursively
layout_path = EegFun.find_file("biosemi72.csv", EegFun.example_path("layouts"))
```

:::

## See Also

- [API Reference](../../reference/index.md)
