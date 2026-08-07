# Demo: Lateralised Readiness Potential (LRP)
# Shows how to calculate the LRP from left-hand and right-hand ERP conditions
# and perform batch processing across participants.

# Note: EegFun.example_path() resolves bundled example data paths.
# When using your own data, simply pass the file path directly, e.g.:
# dat = EegFun.read_raw_data("/path/to/your/data.bdf")

using EegFun
using GLMakie

const DEMO_OUTPUT = "./tutorials/output/"
mkpath(DEMO_OUTPUT)

#######################################################################
# LOAD ERP DATA
#######################################################################

# Load a participant's ERPs
erps = EegFun.read_data(EegFun.example_path("data/julia/erps/example1_erps_final.jld2"))

#######################################################################
# LATERALISED READINESS POTENTIAL (IN-MEMORY)
#######################################################################

# Calculate LRP from two conditions (e.g., Left response vs Right response)
# Assuming condition 1 is Left Hand and condition 2 is Right Hand
# The default channel_selection = EegFun.channels() auto-detects all 
# odd/even lateral pairs (e.g., C3/C4, C1/C2).
lrp_data = EegFun.lrp(erps[1], erps[2])

# Plot the LRP (e.g., at C3/C4)
EegFun.plot_erp(lrp_data, channel_selection = EegFun.channels([:C3, :C4]))

# Calculate LRP for multiple condition pairs simultaneously from an array of ERPs
# e.g., Left vs Right (1 vs 2) and Left vs Right under a different condition (3 vs 4)
pairs = [(1, 2), (3, 4)]
lrp_results = EegFun.lrp(erps, pairs)


#######################################################################
# BATCH LRP (FROM DISK)
#######################################################################

# Batch calculate LRPs for a list of condition pairs across all participant files
EegFun.lrp(
    "erps_final",
    [(1, 2)],
    input_dir = EegFun.example_path("data/julia/erps/"),
    output_dir = DEMO_OUTPUT
)
