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
