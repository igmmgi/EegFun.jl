using eegfun
using GLMakie

# Get some basic data with initial preprocessing steps (high-pass filter, epoch)
data_file = joinpath(@__DIR__, "..", "..", "..",  "Flank_C_3.bdf")
layout_file = eegfun.read_layout("./data/layouts/biosemi/biosemi72.csv");

dat = eegfun.read_bdf(data_file);
dat = eegfun.create_eeg_dataframe(dat, layout_file);

eegfun.rereference!(dat, :avg)

eegfun.filter_data!(dat, "hp", 1)

# run ICA
ica_result = eegfun.run_ica(dat);

# plot ICA components
eegfun.plot_ica_topoplot(ica_result);


eegfun.plot_ica_topoplot(ica_result, dat.layout, component_selection = eegfun.components(1:15));
eegfun.plot_ica_topoplot(ica_result, dat.layout, component_selection = eegfun.components(1:15), use_global_scale = true);


