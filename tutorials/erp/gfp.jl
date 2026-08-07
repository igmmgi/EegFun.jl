# Demo: Global Field Power (GFP) and Global Dissimilarity
# Shows how to compute reference-independent measures of response strength
# and topographic change from ERP data.

# Note: EegFun.example_path() resolves bundled example data paths.
# When using your own data, simply pass the file path directly, e.g.:
# dat = EegFun.read_raw_data("/path/to/your/data.bdf")

using EegFun
using GLMakie

#######################################################################
# LOAD ERP DATA
#######################################################################
dat = EegFun.read_data(EegFun.example_path("data/julia/erps/example1_erps_final.jld2"))

#######################################################################
# GLOBAL FIELD POWER (GFP)
#######################################################################

# GFP = standard deviation across all channels at each time point
# High GFP = strong, synchronized activity

# All channels
gfp_result = EegFun.gfp(dat)

# Specific channels
gfp_result = EegFun.gfp(dat, channel_selection = EegFun.channels([:Cz, :Pz, :Fz]))

# Normalized to 0-100% (for comparing across datasets)
gfp_result = EegFun.gfp(dat, normalize = true)


#######################################################################
# GLOBAL DISSIMILARITY
#######################################################################

# Global dissimilarity = rate of topographic change over time
# Peaks indicate transitions between brain states or ERP components

gd_result = EegFun.global_dissimilarity(dat)

# With normalization
gd_result = EegFun.global_dissimilarity(dat, normalize = true)


#######################################################################
# COMBINED CALCULATION
#######################################################################

# Calculate both GFP and global dissimilarity in one call
result = EegFun.gfp_and_dissimilarity(dat)

# Access the values (one result per condition)
result[1].time
result[1].gfp
result[1].dissimilarity


#######################################################################
# MULTIPLE CONDITIONS
#######################################################################

# GFP for all conditions
gfp_all = EegFun.gfp(dat)

# GFP for a specific condition
gfp_cond1 = EegFun.gfp(dat, condition_selection = EegFun.conditions([1]))
