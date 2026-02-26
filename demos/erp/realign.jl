# Demo: Response-Locked Realignment
# Shows how to realign stimulus-locked epochs to a different time point
# (e.g., response time) for response-locked ERP analysis.

using EegFun

#######################################################################
# SINGLE-PARTICIPANT REALIGNMENT
#######################################################################

# Load stimulus-locked epoched data
# epochs = EegFun.read_data("participant_1_epochs.jld2")

# Realign to response times stored in the :rt column
# EegFun.realign!(epochs, :rt)

# Non-mutating version (returns a new copy)
# realigned = EegFun.realign(epochs, :rt)


#######################################################################
# BATCH REALIGNMENT
#######################################################################

# Process all participant epoch files in a directory.
# Each epoch must have a column with the realignment time
# (e.g., response time relative to stimulus onset).

# EegFun.realign("epochs_cleaned", :rt)

# Specify input directory
# EegFun.realign("epochs_cleaned", :rt, input_dir = "/path/to/epochs/")

# Specific participants
# EegFun.realign("epochs_cleaned", :rt, participant_selection = EegFun.participants([1, 2, 3]))


#######################################################################
# TYPICAL WORKFLOW
#######################################################################

# Response-locked LRP analysis:
#
# Extract stimulus-locked epochs with RT column
# Realign to response time:
#      realign("epochs_cleaned", :rt)
# Average to ERPs:
#      average_epochs in realigned directory
# Calculate LRP:
#      lrp("realigned_erps", [(1, 2)])
# Jackknife average:
#      jackknife_average("lrp")
