using EegFun
using JLD2
# using GLMakie

# Get some basic data with initial preprocessing steps (high-pass filter, epoch)
data_file = joinpath(@__DIR__, "..", "..", "..", "AttentionExp", "Flank_C_3.bdf")

# data_file = joinpath(@__DIR__, "..", "..", "..", "AttentionExp", "recoded", "Flank_C_3.bdf")
layout_file = EegFun.read_layout("./data/layouts/biosemi/biosemi72.csv");
EegFun.polar_to_cartesian_xy!(layout_file)
dat = EegFun.read_raw_data(data_file);

dat = EegFun.create_eeg_dataframe(dat, layout_file);
# EegFun.highpass_filter!(dat, 1)
EegFun.rereference!(dat, :avg)

EegFun.plot_databrowser(dat);
# fig, ax = EegFun.plot_databrowser!(dat)

# return some analysis settings
fig, ax, analysis_settings = EegFun.plot_databrowser(dat)
dat_new = EegFun.apply_analysis_settings(dat, analysis_settings)
EegFun.plot_databrowser(dat_new)

EegFun.highpass_filter!(dat, 1)
fig, ax = EegFun.plot_databrowser(dat)

EegFun.lowpass_filter!(dat, 20)
EegFun.plot_databrowser(dat)

EegFun.is_extreme_value!(dat, 500);
EegFun.mark_epoch_windows!(dat, [1, 3, 4, 5], [-1, 1.0]) # simple epoch marking with trigger 1 and 3
EegFun.plot_databrowser(dat)

# try some custom styling
EegFun.plot_databrowser(dat; :channel_line_width => 10, :selection_color => (:red, 0.5))

# with ica 
data_file = joinpath(@__DIR__, "..", "..", "..", "AttentionExp", "Flank_C_3.bdf")
layout_file = EegFun.read_layout("./data/layouts/biosemi/biosemi72.csv");
EegFun.polar_to_cartesian_xy!(layout_file)
dat = EegFun.read_raw_data(data_file)
dat = EegFun.create_eeg_dataframe(dat, layout_file)
EegFun.rereference!(dat, :avg)
EegFun.highpass_filter!(dat, 1)
EegFun.is_extreme_value!(dat, 100);

ica_result = EegFun.run_ica(dat; sample_selection = EegFun.samples_not(:is_extreme_value_100), percentage_of_data = 25)
EegFun.plot_databrowser(dat, ica_result)

# EPOCHS
epoch_cfg = [EegFun.EpochCondition(name = "ExampleEpoch1", trigger_sequences = [[1]])]
epochs = EegFun.EpochData[]
for (idx, epoch) in enumerate(epoch_cfg)
    push!(epochs, EegFun.extract_epochs(dat, idx, epoch, -2, 4))
end

EegFun.plot_databrowser(epochs[1])
EegFun.plot_databrowser(epochs[1], ica_result)

# plot_databrowser epochs from file
jldsave("dat.jld2", data = dat)
jldsave("epochs.jld2", data = epochs)
jldsave("ica.jld2", data = ica_result)

EegFun.plot_databrowser("dat.jld2")
EegFun.plot_databrowser("epochs.jld2")
EegFun.plot_databrowser("dat.jld2", "ica.jld2")
EegFun.plot_databrowser("epochs.jld2", "ica.jld2")
