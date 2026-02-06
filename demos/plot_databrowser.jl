using EegFun
using JLD2

# read raw data
dat = EegFun.read_raw_data("./resources/data/bdf/example1.bdf");

# read and preprate layout file
layout_file = EegFun.read_layout("./resources/layouts/biosemi/biosemi72.csv");
EegFun.polar_to_cartesian_xy!(layout_file)

# create EegFun data structure (EegFun.ContinuousData)
dat = EegFun.create_eegfun_data(dat, layout_file);

dat = EegFun.create_eegfun_data(dat);


# Basic databrowser
EegFun.plot_databrowser(dat);

# return some analysis settings
# check highpass filter, lowpass filter, rereference and mark some regions
# these are returned in analysis_settings
_, _, analysis_settings = EegFun.plot_databrowser(dat)

# we can now apply these settings to the data
dat_new = EegFun.apply_analysis_settings(dat, analysis_settings)

# settings applied
EegFun.plot_databrowser(dat_new)

#######################################################################
# We can add new columns to our data
dat = EegFun.read_raw_data("./resources/data/bdf/example1.bdf");
dat = EegFun.create_eegfun_data(dat, layout_file);

# Some minimal preprocessing (average reference and highpass filter)
EegFun.rereference!(dat, :avg)
EegFun.highpass_filter!(dat, 0.5)

# E.g. mark data sections (some extreme values and intervals around triggers; Boolean)
EegFun.is_extreme_value!(dat, 100);
EegFun.mark_epoch_windows!(dat, [1, 2], [-0.2, 1.0]) # simple epoch marking with trigger 1 and 3

# Add new channel with a 10 Hz sine wave (50 μV amplitude)
dat.data[:, :sinewave] = 50.0 .* sin.(2π .* 10.0 .* dat.data[:, :time])

# these channels can now be viewed in the databrowser under the Extra channels tab
EegFun.plot_databrowser(dat)

EegFun.all_labels(dat)
EegFun.meta_labels(dat)
EegFun.channel_labels(dat)
EegFun.extra_labels(dat)

# try some custom styling
EegFun.plot_databrowser(dat; :channel_line_width => 2, :selection_color => (:green, 0.1))


#######################################################################
# run ica
ica_result = EegFun.run_ica(dat; sample_selection = EegFun.samples_not(:is_extreme_value_100), percentage_of_data = 25)
EegFun.plot_databrowser(dat, ica_result)


#######################################################################
# EPOCHS
epoch_cfg = [EegFun.EpochCondition(name = "ExampleEpoch1", trigger_sequences = [[1]])]
epochs = EegFun.extract_epochs(dat, epoch_cfg, (-0.5, 1.0))

# databrowser for epochs
EegFun.plot_databrowser(epochs[1])  # can now browse epochs
EegFun.plot_databrowser(epochs[1], ica_result)


# We can save our EegFun data to a file and pass a file name to plot_databrowser 
# plot_databrowser epochs from file
jldsave("dat.jld2", data = dat)
jldsave("epochs.jld2", data = epochs)
jldsave("ica.jld2", data = ica_result)

# Plot by loading from file
EegFun.plot_databrowser("dat.jld2")
EegFun.plot_databrowser("epochs.jld2")
EegFun.plot_databrowser("dat.jld2", "ica.jld2")
EegFun.plot_databrowser("epochs.jld2", "ica.jld2")
