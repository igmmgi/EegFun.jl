# Demo: Jackknife Averaging
# Shows how to create leave-one-out (jackknife) averages for statistical
# testing of ERP latency measures.

using EegFun

#######################################################################
# SINGLE-PARTICIPANT JACKKNIFE
#######################################################################

# Jackknife averaging creates N averages from N participants,
# where each average excludes one participant.
# This is used for statistical testing of peak latency measures.

# Load multiple participant ERPs
# erps = [EegFun.read_data("$(i)_erps.jld2") for i in 1:20]

# Create jackknife averages for each condition
# jackknife_erps = EegFun.jackknife_average(erps)


#######################################################################
# BATCH JACKKNIFE
#######################################################################

# Process all participant ERP/LRP files in a directory

# EegFun.jackknife_average("erps_cleaned")

# LRP data (common use case)
# EegFun.jackknife_average("lrp")

# Specific conditions only
# EegFun.jackknife_average("lrp", condition_selection = EegFun.conditions([1]))


#######################################################################
# TYPICAL WORKFLOW
#######################################################################

# Jackknife method for onset latency:
#
# Calculate LRP:               lrp("erps", [(1, 2)])
# Jackknife average:           jackknife_average("lrp")
# Measure onset latency in each jackknife average
# Apply jackknife correction to the t-test:
#      t_corrected = t_original / (n - 1)
#
# Reference: Miller, Patterson, & Ulrich (1998)
