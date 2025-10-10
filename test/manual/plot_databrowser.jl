using eegfun
using GLMakie

# Get some basic data with initial preprocessing steps (high-pass filter, epoch)
data_file = joinpath(@__DIR__, "..", "..", "..",  "Flank_C_3.bdf")
layout_file = eegfun.read_layout("./data/layouts/biosemi/biosemi72.csv");
eegfun.polar_to_cartesian_xy!(layout_file)

dat = eegfun.read_bdf(data_file);
dat = eegfun.create_eeg_dataframe(dat, layout_file);

eegfun.rereference!(dat, :avg)
eegfun.plot_databrowser(dat)

eegfun.filter_data!(dat, "hp", 1)
eegfun.plot_databrowser(dat)

eegfun.filter_data!(dat, "lp", 20)
eegfun.plot_databrowser(dat)

# try some custom styling
eegfun.plot_databrowser(dat; :line_width => 10, :selection_color => (:red, 0.5))


# with ica 
data_file = joinpath(@__DIR__, "..", "..", "..",  "Flank_C_3.bdf")
layout_file = eegfun.read_layout("./data/layouts/biosemi/biosemi72.csv");
eegfun.polar_to_cartesian_xy!(layout_file)
dat = eegfun.read_bdf(data_file);
dat = eegfun.create_eeg_dataframe(dat, layout_file);
eegfun.rereference!(dat, :avg)
eegfun.filter_data!(dat, "hp", 1)
eegfun.is_extreme_value!(dat, 100);


ica_result = eegfun.run_ica(dat; sample_selection = eegfun.samples_not(:is_extreme_value_100), percentage_of_data = 25)
eegfun.plot_databrowser(dat, ica_result)
