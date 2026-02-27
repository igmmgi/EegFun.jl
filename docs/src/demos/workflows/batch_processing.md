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
# Shows how to load and process multiple files, read saved data,
# and group results by condition.

using EegFun

#######################################################################
# LOADING MULTIPLE FILES
#######################################################################

# read_all_data loads all matching files from a directory
# Pattern matching: looks for files containing the pattern string
# erps = EegFun.read_all_data("_erps_good", input_dir = "/path/to/output")

# With participant selection (e.g., only participants 1-10)
# erps = EegFun.read_all_data("_erps_good", "/path/to/output", EegFun.participants(1:10))


#######################################################################
# READING SINGLE FILES
#######################################################################

# read_data automatically handles JLD2 files with EegFun data types
# data = EegFun.read_data("/path/to/output/participant_01_erps_good.jld2")


#######################################################################
# GROUP BY CONDITION
#######################################################################

# Create some test data to demonstrate grouping
erps = EegFun.create_batch_test_erp_data(n_participants = 5)

# Group by condition — returns OrderedDict{Int, Vector{ErpData}}
grouped = EegFun.group_by_condition(erps)

# Access specific condition groups
for (cond_num, cond_erps) in grouped
    println("Condition $cond_num: $(length(cond_erps)) ERPs")
end


#######################################################################
# GRAND AVERAGE FROM BATCH DATA
#######################################################################

# Average across participants for each condition
grand_avg = EegFun.grand_average(erps)


#######################################################################
# FILE MANAGEMENT
#######################################################################

# Check if files exist
EegFun.check_files_exist(["./resources/data/bdf/example1.bdf"])

# Get files matching a pattern in a directory
files = EegFun.get_files("./resources/data/bdf", ".*\\.bdf")

# Find a file by searching directories recursively
found = EegFun.find_file("biosemi72.csv", "./resources/layouts")
```

:::

## See Also

- [API Reference](../../reference/index.md)
