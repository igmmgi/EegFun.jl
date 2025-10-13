
eegfun.all_data(dat)
eegfun.meta_data(dat)
eegfun.channel_data(dat)
eegfun.extra_data(dat)

dat_subset = eegfun.subset(dat, channel_selection = eegfun.channels([:Fp1, :Fp2]))
dat_subset = eegfun.subset(dat, sample_selection = x -> x.sample .<= 10_000) # first 10000 samples
dat_subset = eegfun.subset(dat, sample_selection = x -> x.time .<= 10) # first 10 seconds

dat_subset = eegfun.subset(dat, channel_selection = eegfun.channels([:Fp1, :Fp2, :vEOG, :hEOG]))

# # save / load
# save_object("$(subject)_continuous.jld2", dat)
# dat1 = load_object("3_continuous.jld2")

