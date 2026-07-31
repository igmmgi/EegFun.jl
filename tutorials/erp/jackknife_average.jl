# Demo: Jackknife Averaging
# Shows how to create leave-one-out (jackknife) averages for statistical
# testing of ERP latency measures.

# Note: EegFun.example_path() resolves bundled example data paths.
# When using your own data, simply pass the file path directly, e.g.:
# dat = EegFun.read_raw_data("/path/to/your/data.bdf")

using EegFun

const DEMO_OUTPUT = "./tutorials/output/"
mkpath(DEMO_OUTPUT)

#######################################################################
# LOAD MULTI-PARTICIPANT ERP DATA
#######################################################################

# We need multiple participants for jackknifing
erps_p1 = EegFun.read_data(EegFun.example_path("data/julia/erps/example1_erps_final.jld2"))
erps_p2 = EegFun.read_data(EegFun.example_path("data/julia/erps/example2_erps_final.jld2"))
erps_p3 = EegFun.read_data(EegFun.example_path("data/julia/erps/example3_erps_final.jld2"))

# Create a list of one condition across participants (e.g., Condition 1)
cond_1_erps = [erps_p1[1], erps_p2[1], erps_p3[1]]

#######################################################################
# JACKKNIFE AVERAGING (IN-MEMORY)
#######################################################################

# Calculate jackknife averages: For N participants, this creates N averages 
# where the i-th average includes all participants EXCEPT participant i.
jackknife_erps = EegFun.jackknife_average(cond_1_erps)

# Plot the jackknife averages
EegFun.plot_erp(jackknife_erps, channel_selection = EegFun.channels(:Pz))


#######################################################################
# BATCH JACKKNIFE AVERAGING (FROM DISK)
#######################################################################

# Averages all JLD2 files whose name contains "erps_final" in the given directory.
# This automatically handles multiple conditions and saves the leave-one-out 
# datasets for each participant.
EegFun.jackknife_average("erps_final", input_dir = EegFun.example_path("data/julia/erps/"), output_dir = DEMO_OUTPUT)
