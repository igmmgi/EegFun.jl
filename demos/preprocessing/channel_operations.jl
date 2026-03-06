# Demo: Channel Operations
# Shows how to average, delete, and compute differences between channels.
# Covers channel_average!, channel_delete!, and channel_difference!.

using EegFun

#######################################################################
# LOAD DATA
#######################################################################

dat = EegFun.read_raw_data("./resources/data/bdf/example1.bdf")
layout = EegFun.read_layout("./resources/layouts/biosemi/biosemi72.csv")
EegFun.polar_to_cartesian_xy!(layout)

dat = EegFun.create_test_eegfun_data(dat, layout)

# minimal preprocessing
EegFun.rereference!(dat, :avg)
EegFun.highpass_filter!(dat, 0.1)


#######################################################################
# CHANNEL DIFFERENCE (EOG CALCULATION)
#######################################################################

# Vertical EOG: average of upper channels minus average of lower channels
EegFun.channel_difference!(
    dat,
    channel_selection1 = EegFun.channels([:Fp1, :Fp2]),
    channel_selection2 = EegFun.channels([:IO1, :IO2]),
    channel_out = :vEOG,
)
# Horizontal EOG: left minus right
EegFun.channel_difference!(
    dat,
    channel_selection1 = EegFun.channels([:F9]),
    channel_selection2 = EegFun.channels([:F10]),
    channel_out = :hEOG,
)



#######################################################################
# EOG WITH CONFIGURATION
#######################################################################

# Calculate both EOG channels from a configuration
eog_cfg = EegFun.EogConfig(
    vEOG_criterion = 50.0,
    hEOG_criterion = 30.0,
    vEOG_channels = [["Fp1", "Fp2"], ["IO1", "IO2"], ["vEOG"]],
    hEOG_channels = [["F9"], ["F10"], ["hEOG"]],
)
EegFun.calculate_eog_channels!(dat, eog_cfg)


#######################################################################
# CHANNEL AVERAGING
#######################################################################

# Average all channels into a single column
EegFun.channel_average!(dat)

dat.data[!, :avg] # average of all channels

# Average specific channel groups (average of Fp1 and Fp2, average of O1 and O2 with 
# default output labels of :Fp1_Fp2 and :O1_O2)
EegFun.channel_average!(dat, channel_selections = [EegFun.channels([:Fp1, :Fp2]), EegFun.channels([:O1, :O2])])

# Average with custom output labels
EegFun.channel_average!(
    dat,
    channel_selections = [EegFun.channels([:Fp1, :Fp2]), EegFun.channels([:O1, :O2])],
    output_labels = [:avg1, :avg2],
)

# Reduce to just the averaged channels (drops originals)
dat_reduced = EegFun.channel_average(dat, channel_selections = [EegFun.channels([:Fp1, :Fp2])], reduce = true)


#######################################################################
# CHANNEL DELETION
#######################################################################

# Delete a single channel
EegFun.channel_delete!(dat, :avg)
EegFun.channel_delete!(dat, [:Fp1_Fp2, :O1_O2])

# Delete multiple channels
EegFun.channel_delete!(dat, [:vEOG, :hEOG])

# Non-mutating version (returns a copy)
dat_subset = EegFun.channel_delete(dat, [:Fp1, :Fp2])
