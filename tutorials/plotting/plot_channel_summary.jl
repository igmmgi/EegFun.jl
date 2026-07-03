# Demo: Channel Summary Plots
# Shows visualization of channel-wise summary statistics.

# Note: EegFun.example_path() resolves bundled example data paths.
# When using your own data, simply pass the file path directly, e.g.:
# dat = EegFun.read_raw_data("/path/to/your/data.bdf")

using EegFun

# read raw data
dat = EegFun.read_raw_data(EegFun.example_path("data/bdf/example1.bdf"));
dat = EegFun.create_eegfun_data(dat);

# Some minimal preprocessing (average reference and highpass filter)
EegFun.rereference!(dat, :avg)
EegFun.highpass_filter!(dat, 1)

# basic channel summary statistics
cs = EegFun.channel_summary(dat)
cs = EegFun.channel_summary(dat, channel_selection = EegFun.channels([:Fp1, :Fp2]))
cs = EegFun.channel_summary(dat, channel_selection = EegFun.channels([:Fp1, :Fp2]), sample_selection = x -> x.sample .< 2000)
cs = EegFun.channel_summary(dat, channel_selection = x -> endswith.(string.(x), "z")) # all midline channels 
cs = EegFun.channel_summary(dat, channel_selection = x -> .!(endswith.(string.(x), "z"))) # all non-midline channels 

# Plotting Channel Summaries
EegFun.plot_channel_summary(cs, :range)
EegFun.plot_channel_summary(cs, :min)
EegFun.plot_channel_summary(cs, :min, bar_color = :red)
EegFun.plot_channel_summary(cs, [:min, :max, :std, :range, :var, :zvar])
EegFun.plot_channel_summary(cs, [:range, :var])

#################################
# Epoched DataFrameEeg
#################################
# some epoched data
epoch_cfg = EegFun.EpochCondition(name = "ExampleEpoch1", trigger_sequences = [[1]])
epochs = EegFun.extract_epochs(dat, epoch_cfg, (-2, 4))
cs = EegFun.channel_summary(epochs)

EegFun.plot_channel_summary(cs, :range, average_over = :epoch)
EegFun.plot_channel_summary(cs, [:min, :max, :std, :range, :var, :zvar], average_over = :epoch)

# some epoched data
epoch_cfg = [
    EegFun.EpochCondition(name = "ExampleEpoch1", trigger_sequences = [[1]]),
    EegFun.EpochCondition(name = "ExampleEpoch2", trigger_sequences = [[2]]),
]
epochs = EegFun.extract_epochs(dat, epoch_cfg, (-2, 4))
cs = EegFun.channel_summary(epochs)

EegFun.plot_channel_summary(cs, :range, average_over = :epoch)
EegFun.plot_channel_summary(cs, [:min, :max, :std, :range, :var, :zvar], average_over = :epoch)

