using eegfun
using JLD2
# using GLMakie

# Get some basic data with initial preprocessing steps (high-pass filter, epoch)
data_file = joinpath(@__DIR__, "..", "..", "..", "AttentionExp", "Flank_C_3.bdf")
# data_file = joinpath(@__DIR__, "..", "..", "..", "AttentionExp", "recoded", "Flank_C_3.bdf")
layout_file = eegfun.read_layout("./data/layouts/biosemi/biosemi72.csv");
eegfun.polar_to_cartesian_xy!(layout_file)
dat = eegfun.read_bdf(data_file);
dat = eegfun.create_eeg_dataframe(dat, layout_file);
eegfun.filter_data!(dat, "hp", 1)
eegfun.rereference!(dat, :avg)

eegfun.plot_databrowser(dat);
# fig, ax = eegfun.plot_databrowser!(dat)


# return some analysis settings
fig, ax, analysis_settings = eegfun.plot_databrowser(dat)
dat_new = eegfun.apply_analysis_settings(dat, analysis_settings)
eegfun.plot_databrowser(dat_new)


eegfun.filter_data!(dat, "hp", 1)
fig, ax = eegfun.plot_databrowser(dat)

eegfun.filter_data!(dat, "lp", 20)
eegfun.plot_databrowser(dat)

eegfun.is_extreme_value!(dat, 500);
eegfun.mark_epoch_windows!(dat, [1, 3, 4, 5], [-1, 1.0]) # simple epoch marking with trigger 1 and 3
eegfun.plot_databrowser(dat)

# try some custom styling
eegfun.plot_databrowser(dat; :channel_line_width => 10, :selection_color => (:red, 0.5))

# with ica 
data_file = joinpath(@__DIR__, "..", "..", "..", "AttentionExp", "Flank_C_3.bdf")
layout_file = eegfun.read_layout("./data/layouts/biosemi/biosemi72.csv");
eegfun.polar_to_cartesian_xy!(layout_file)
dat = eegfun.read_bdf(data_file);
dat = eegfun.create_eeg_dataframe(dat, layout_file);
eegfun.rereference!(dat, :avg)
eegfun.filter_data!(dat, "hp", 1)
eegfun.is_extreme_value!(dat, 100);

ica_result = eegfun.run_ica(dat; sample_selection = eegfun.samples_not(:is_extreme_value_100), percentage_of_data = 25)
eegfun.plot_databrowser(dat, ica_result)

# EPOCHS
epoch_cfg = [eegfun.EpochCondition(name = "ExampleEpoch1", trigger_sequences = [[1]])]
epochs = eegfun.EpochData[]
for (idx, epoch) in enumerate(epoch_cfg)
    push!(epochs, eegfun.extract_epochs(dat, idx, epoch, -2, 4))
end

eegfun.plot_databrowser(epochs[1])
eegfun.plot_databrowser(epochs[1], ica_result)



# plot_databrowser epochs from file
jldsave( "dat.jld2", data = dat)
jldsave( "epochs.jld2", data = epochs)
jldsave( "ica.jld2", data = ica_result)


eegfun.plot_databrowser("dat.jld2")
eegfun.plot_databrowser("epochs.jld2")
eegfun.plot_databrowser("dat.jld2", "ica.jld2")
eegfun.plot_databrowser("epochs.jld2", "ica.jld2")





