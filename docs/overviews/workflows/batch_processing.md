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
