# package
using eegfun
using GLMakie
# using CairoMakie
# using BenchmarkTools

function get_data() 
dat = eegfun.read_bdf("../../Flank_C_3.bdf");
layout = eegfun.read_layout("../../data/layouts/biosemi72.csv");
dat = eegfun.create_eeg_dataframe(dat, layout);
return dat, layout
end


# Test plot_channel_summary
dat, layout = get_data()

cs = eegfun.channel_summary(dat) 
eegfun.plot_channel_summary(cs, :max)

