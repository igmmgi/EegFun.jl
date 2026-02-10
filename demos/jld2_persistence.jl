# Demo: Data Persistence with JLD2
# Shows how to save and load processed EEG data using JLD2 format,
# including batch loading patterns and file organization best practices.

using EegFun
using JLD2

#######################################################################
# 1. BASIC SAVING AND LOADING
#######################################################################

# Read and prepare some data
dat = EegFun.read_raw_data("./resources/data/bdf/example1.bdf")
layout = EegFun.read_layout("./resources/layouts/biosemi/biosemi72.csv")
EegFun.polar_to_cartesian_xy!(layout)
dat = EegFun.create_eegfun_data(dat, layout)

# Some preprocessing
EegFun.rereference!(dat, :avg)
EegFun.highpass_filter!(dat, 0.1)

# Save continuous data
jldsave("continuous_data.jld2"; data = dat)

# Load it back
loaded_dat = load("continuous_data.jld2", "data")

# Verify it's the same
EegFun.all_data(loaded_dat) == EegFun.all_data(dat)


#######################################################################
# 2. SAVING EPOCHED DATA
#######################################################################

# Extract epochs
epoch_cfg = [
    EegFun.EpochCondition(name = "Condition1", trigger_sequences = [[1]]),
    EegFun.EpochCondition(name = "Condition2", trigger_sequences = [[2]]),
]
epochs = EegFun.extract_epochs(dat, epoch_cfg, (-0.2, 1.0))

# Save all conditions
jldsave("epochs_all.jld2"; data = epochs)

# Save individual conditions
jldsave("epochs_condition1.jld2"; data = epochs[1])
jldsave("epochs_condition2.jld2"; data = epochs[2])

# Load epochs back
loaded_epochs = load("epochs_all.jld2", "data")


#######################################################################
# 3. SAVING ERPs
#######################################################################

# Create ERPs
erps = EegFun.average_epochs(epochs)

# Save ERPs
jldsave("erps.jld2"; data = erps)

# Load ERPs
loaded_erps = load("erps.jld2", "data")


#######################################################################
# 4. SAVING ICA RESULTS
#######################################################################

# Run ICA
EegFun.is_extreme_value!(dat, 100)
ica_result = EegFun.run_ica(dat, sample_selection = EegFun.samples_not(:is_extreme_value_100))

# Save ICA decomposition
jldsave("ica_decomposition.jld2"; data = ica_result)

# Load ICA
loaded_ica = load("ica_decomposition.jld2", "data")

# Can use loaded ICA with databrowser
EegFun.plot_databrowser(dat, loaded_ica)


#######################################################################
# 5. SAVING MULTIPLE ITEMS IN ONE FILE
#######################################################################

# Save multiple related objects together
jldsave("analysis_results.jld2"; epochs = epochs, erps = erps, ica = ica_result)

# Load specific items
my_epochs = load("analysis_results.jld2", "epochs")
my_erps = load("analysis_results.jld2", "erps")

# Or load everything
results = load("analysis_results.jld2")
results["epochs"]
results["erps"]
results["ica"]


#######################################################################
# 6. BATCH LOADING FROM DIRECTORY
#######################################################################

# Suppose you have multiple participant files:
# - participant_01_epochs.jld2
# - participant_02_epochs.jld2
# - participant_03_epochs.jld2

# List all matching files
using Glob
epoch_files = glob("participant_*_epochs.jld2", "my_output_directory/")

# Load all participants
all_participant_epochs = [load(file, "data") for file in epoch_files]

# Now you can do group analysis
# grand_avg = EegFun.grand_average(all_participant_epochs)


#######################################################################
# 7. ORGANIZED FILE STRUCTURE
#######################################################################

# Best practice: organize outputs by processing stage
# 
# project/
# ├── data/
# │   └── raw/              # Original BDF/VHDR files
# ├── derivatives/
# │   ├── continuous/       # Filtered, rereferenced continuous data
# │   ├── ica/              # ICA decompositions
# │   ├── epochs/           # Epoched data
# │   ├── epochs_clean/     # Epochs after artifact rejection
# │   └── erps/             # Averaged ERPs
# └── results/
#     └── grand_averages/   # Group-level results

# Example: saving with organized structure
output_dir = "derivatives/epochs"
participant_id = "sub-01"
jldsave(joinpath(output_dir, "$(participant_id)_epochs.jld2"); data = epochs)


#######################################################################
# 8. NAMING CONVENTIONS
#######################################################################

# Use consistent naming for easy batch processing:
# - Descriptive prefixes: "sub-", "participant_", etc.
# - Stage indicators: "_continuous", "_epochs", "_erps", "_clean"
# - Zero-padding for sorting: "sub-01" not "sub-1"

# Examples:
# jldsave("sub-01_continuous.jld2"; data = dat)
# jldsave("sub-01_epochs_raw.jld2"; data = epochs)
# jldsave("sub-01_epochs_clean.jld2"; data = epochs_cleaned)
# jldsave("sub-01_erps.jld2"; data = erps)


#######################################################################
# 9. READING DATA IN PLOT FUNCTIONS
#######################################################################

# Many EegFun plot functions accept file paths directly
EegFun.plot_databrowser("continuous_data.jld2")
EegFun.plot_databrowser("epochs_all.jld2")
EegFun.plot_databrowser("continuous_data.jld2", "ica_decomposition.jld2")

# This is convenient for quick visualization without loading into memory


#######################################################################
# 10. MEMORY CONSIDERATIONS
#######################################################################

# For large datasets, load only what you need
# Instead of: results = load("big_file.jld2")
# Use: epochs = load("big_file.jld2", "epochs")

# For very large cohorts, process one participant at a time
# rather than loading all into memory at once
