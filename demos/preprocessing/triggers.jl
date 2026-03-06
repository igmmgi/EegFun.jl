# Demo: Triggers
# Shows how to inspect, count, and search for trigger sequences in EEG data.

using EegFun

#######################################################################
# LOAD DATA
#######################################################################

dat = EegFun.read_raw_data("./resources/data/bdf/example1.bdf")
layout = EegFun.read_layout("./resources/layouts/biosemi/biosemi72.csv")
EegFun.polar_to_cartesian_xy!(layout)
dat = EegFun.create_test_eegfun_data(dat, layout)


#######################################################################
# TRIGGER OVERVIEW
#######################################################################

# Count how often each trigger occurs (raw and cleaned)
trigger_info = EegFun.trigger_count(dat)

# View the trigger table
trigger_info

# Access the underlying DataFrame
trigger_info.data


#######################################################################
# SEARCH FOR TRIGGER SEQUENCES
#######################################################################

# Find all onsets of trigger value 1
idx = EegFun.search_sequence(dat.data.trigger, 1)

# Find a two-trigger sequence (e.g., warning cue followed by target)
idx = EegFun.search_sequence(dat.data.trigger, [1, 2])

# Find a sequence while ignoring zeros between triggers
idx = EegFun.search_sequence(dat.data.trigger, [1, 2], ignore_values = [0])


#######################################################################
# WILDCARDS AND RANGES
#######################################################################

# Use :any wildcard to match any trigger between two specific triggers
idx = EegFun.search_sequence(dat.data.trigger, [1, :any, 2])

# Use ranges to match a set of trigger values
idx = EegFun.search_sequence(dat.data.trigger, [1:5])     # triggers 1 through 5
idx = EegFun.search_sequence(dat.data.trigger, [1:3, 10:12])  # multiple ranges


#######################################################################
# MULTIPLE SEQUENCES (OR LOGIC)
#######################################################################

# Find onsets where any of several sequences occurs
idx = EegFun.search_sequence(dat.data.trigger, [[1, 2], [3, 4]])


#######################################################################
# VISUALIZE TRIGGERS
#######################################################################

# Plot trigger timing overview
EegFun.plot_trigger_overview(dat)
EegFun.plot_trigger_timing(dat)
