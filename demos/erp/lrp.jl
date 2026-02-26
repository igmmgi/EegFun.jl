# Demo: Lateralised Readiness Potential (LRP)
# Shows how to calculate the LRP from left-hand and right-hand ERP conditions
# and perform batch processing across participants.

using EegFun

#######################################################################
# LOAD ERP DATA
#######################################################################

dat = EegFun.read_data("./resources/data/julia/erps/example1_erps_good.jld2")


#######################################################################
# SINGLE-PARTICIPANT LRP
#######################################################################

# Calculate LRP from two conditions:
#   - dat[1] = left-hand response ERPs
#   - dat[2] = right-hand response ERPs
# Automatically pairs odd/even channels (e.g., C3/C4, CP3/CP4)

lrp_data = EegFun.lrp(dat[1], dat[2])

# Select specific channels (left hemisphere only, auto-pairs with right)
lrp_data = EegFun.lrp(dat[1], dat[2], channel_selection = EegFun.channels([:C3, :CP3]))


#######################################################################
# BATCH LRP (MULTIPLE PARTICIPANTS)
#######################################################################

# Process all participant ERP files in a directory.
# condition_pairs specifies (left_condition, right_condition)

# EegFun.lrp("erps_cleaned", [(1, 2)])

# Multiple condition pairs
# EegFun.lrp("erps_cleaned", [(1, 2), (3, 4)])

# Specific channels only
# EegFun.lrp("erps_cleaned", [(1, 2)], channel_selection = EegFun.channels([:C3]))


#######################################################################
# VISUALIZE
#######################################################################

EegFun.plot_erp(lrp_data)
