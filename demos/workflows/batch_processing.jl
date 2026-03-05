# Demo: Batch Processing
# Shows how to load multi-participant ERP data, group by condition,
# and compute grand averages.

using EegFun

const ERP_DIR = "./resources/data/julia/erps/"


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
# EegFun.grand_average("erps_good", input_dir = ERP_DIR, output_dir = "./resources/data/julia/grand_average/")


#######################################################################
# FILE UTILITIES
#######################################################################

# List all .bdf files in a directory (regex pattern)
bdf_files = EegFun.get_files("./resources/data/bdf", "\\.bdf")

# Verify all files exist before starting a long pipeline run
EegFun.check_files_exist(bdf_files)

# Find a single file by searching recursively
layout_path = EegFun.find_file("biosemi72.csv", "./resources/layouts")
