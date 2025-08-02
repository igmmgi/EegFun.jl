# using Logging # Show all messages
# global_logger(ConsoleLogger(stderr, Logging.Debug))
# # Show only info and above
# global_logger(ConsoleLogger(stderr, Logging.Info))
# # Show only warnings and errors
# global_logger(ConsoleLogger(stderr, Logging.Warn))

# package
using eegfun
using GLMakie
using DataFrames
using BenchmarkTools
using AlgebraOfGraphics

# using CairoMakie
# eegfun.preprocess_eeg_data("pipeline.toml")
# load data
dat = eegfun.read_bdf("../Flank_C_3.bdf");
layout = eegfun.read_layout("./data/layouts/biosemi72.csv");
eegfun.plot_trigger_timing(dat)

# create our eeg ContinuousData type
dat = eegfun.create_eeg_dataframe(dat, layout);
eegfun.plot_trigger_overview(dat)

eegfun.plot_trigger_timing(dat)



eegfun.filter_data!(dat, "hp", 1)
eegfun.rereference!(dat, :avg)
eegfun.channel_difference!(dat, channel_selection1 = eegfun.channels([:Fp1, :Fp2]), channel_selection2 = eegfun.channels([:IO1, :IO2]), channel_out = :vEOG); # vertical EOG = mean(Fp1, Fp2) - mean(IO1, I02)
eegfun.mark_epoch_windows!(dat, [1, 3], [-0.5, 1.0]) # simple epoch marking with trigger 1 and 3;

eegfun.plot_databrowser(dat)

# set_aog_theme!()
# fig = Figure()
# all_data = eegfun.all_data(epochs[1])
# mydata = stack(all_data, [:Fp1, :Fp2], variable_name = :channel, value_name = :value)
# plt =
#     data(mydata) *
#     mapping(:time => "Time [ms]", :value => "Amplitude [μV]", color = :channel => nonnumeric) *
#     visual(Lines) *
#     mapping(layout = :epoch => nonnumeric) 
# # plt = paginate(plt, layout = 4)
# # draw(plt, 2)
# draw(plt)


df1 = eegfun.channel_summary(dat)
df2 = eegfun.channel_summary(epochs[1])

eegfun.plot_channel_summary(df1, :range)
eegfun.plot_channel_summary(df1, :min)

eegfun.plot_channel_summary(df2, :range, average_over=:epoch)
eegfun.plot_channel_summary(df2, :min, average_over=:epoch)






# how to filter data
eegfun.filter_data!(dat, "hp", 1)
eegfun.filter_data!(dat, "lp", 20)

# how to rereference data
eegfun.rereference!(dat, :avg)

# how to create channel differences (e.g. EOG calculations)
eegfun.channel_difference!(dat, channel_selection1 = eegfun.channels([:Fp1, :Fp2]), channel_selection2 = eegfun.channels([:IO1, :IO2]), channel_out = :vEOG); # vertical EOG = mean(Fp1, Fp2) - mean(IO1, I02)
eegfun.channel_difference!(dat, channel_selection1 = eegfun.channels([:F9]),        channel_selection2 = eegfun.channels([:F10]),       channel_out = :hEOG); # vertical EOG = mean(Fp1, Fp2) - mean(IO1, I02)

# how to detect EOG onsets (a new Bool column is added to the data frame: :is_vEOG, :is_hEOG)
eegfun.detect_eog_onsets!(dat, 50, :vEOG, :is_vEOG)
eegfun.detect_eog_onsets!(dat, 30, :hEOG, :is_hEOG)

# how to mark epoch windows (a new Bool column is added to the data frame: :epoch_window)
eegfun.mark_epoch_windows!(dat, [1, 3], [-0.5, 1.0]) # simple epoch marking with trigger 1 and 3

# how to plot the databrowser
eegfun.plot_databrowser(dat)



# define neighbours 2D/3D defined by distance (mm)
eegfun.polar_to_cartesian_xy!(layout);
eegfun.get_layout_neighbours_xy!(layout, 40);
eegfun.polar_to_cartesian_xyz!(layout);
eegfun.get_layout_neighbours_xyz!(layout, 40);

# how to plot layout without/with neighbours
eegfun.plot_layout_2d(layout)
eegfun.plot_layout_2d(layout, neighbours = true)
eegfun.plot_layout_3d(layout, neighbours= true)

eegfun.channel_labels(layout)
eegfun.positions_polar(layout)
eegfun.positions_2D(layout)
eegfun.positions_3D(layout)


# how to subset layout a layout
subset_layout = eegfun.subset_layout(layout, channel_selection = x -> .!endswith.(string.(x), "z"));
eegfun.plot_layout_2d(subset_layout);


# how to print neighbours to a file
eegfun.print_layout_neighbours(layout, "electrode_neighbours_1.toml")
eegfun.print_layout_neighbours(layout.neighbours, "electrode_neighbours_2.toml")

# we can get raw trigger info 
eegfun.trigger_count(dat);
eegfun.plot_trigger_overview(dat)
eegfun.plot_trigger_timing(dat)


eegfun.all_data(dat)
eegfun.meta_data(dat)
eegfun.channel_data(dat)
eegfun.extra_data(dat)

eegfun.all_labels(dat)
eegfun.meta_labels(dat)
eegfun.channel_labels(dat)
eegfun.extra_labels(dat)

# we can get raw trigger info 
eegfun.trigger_count(dat);

eegfun.plot_trigger_overview(dat)
eegfun.plot_trigger_timing(dat)

# how to filter data
eegfun.filter_data!(dat, "hp", 1)

# how to rereference data
eegfun.rereference!(dat, :avg)
# databrowser
eegfun.plot_databrowser(dat)

# Subset channels/samples 
dat_subset = eegfun.subset(dat, channel_selection = eegfun.channels([:Fp1, :Fp2]))
dat_subset = eegfun.subset(dat, sample_selection = x -> x.sample .<= 10_000) # first 10000 samples
dat_subset = eegfun.subset(dat, sample_selection = x -> x.time .<= 10) # first 10 seconds

# or subset data and plot
eegfun.plot_databrowser_subset(dat, channel_selection = eegfun.channels([:Fp1, :Fp2]), sample_selection = x -> x.time .< 20)

# how to calculate channel differences
eegfun.channel_difference!(dat, channels_in1 = eegfun.channels([:Fp1, :Fp2]), channels_in2 = eegfun.channels([:IO1, :IO2]), channel_out = :vEOG); # vertical EOG = mean(Fp1, Fp2) - mean(IO1, I02)

eegfun.channel_difference!(dat, channels_in1 = eegfun.channels([:F9]), channels_in2 = eegfun.channels([:F10]), channel_out = :hEOG); # vertical EOG = mean(Fp1, Fp2) - mean(IO1, I02)

# how to detect EOG onsets
eegfun.detect_eog_onsets!(dat, 50, :vEOG, :is_vEOG)
eegfun.detect_eog_onsets!(dat, 30, :hEOG, :is_hEOG)

dat_subset = eegfun.subset(dat, channel_selection = eegfun.channels([:Fp1, :Fp2, :vEOG, :hEOG]))

eegfun.plot_databrowser_subset(dat, channel_selection = eegfun.channels([:Fp1,:vEOG, :hEOG]))


eegfun.channels(dat)     # original channels in the layout
eegfun.all_channels(dat) # all channels in the data

# basic channel summary statistics
eegfun.channel_summary(dat) 
eegfun.channel_summary(dat, channel_selection = eegfun.channels([:Fp1, :Fp2])) 
eegfun.channel_summary(dat, channel_selection = eegfun.channels([:Fp1, :Fp2]), sample_selection = x -> x.sample .< 2000)
eegfun.channel_summary(dat, channel_selection = x -> endswith.(string.(x), "z")) # all midline channels 
eegfun.channel_summary(dat, channel_selection = x -> .!(endswith.(string.(x), "z"))) # all non-midline channels 
eegfun.channel_summary(dat, include_extra = true) # include additional channels (e.g. vEOG, hEOG)

# add bool columns to the data frame
eegfun.is_extreme_value!(dat, 100);
eegfun.is_extreme_value!(dat, 100; include_extra = true);
eegfun.is_extreme_value!(dat, 100; channel_selection = eegfun.channels_not([:Fp1, :Fp2]));
eegfun.is_extreme_value!(dat, 100; channel_selection = x -> endswith.(string.(x), "z"));
eegfun.is_extreme_value!(dat, 100; channel_selection = x -> .!(endswith.(string.(x), "z")));
eegfun.is_extreme_value!(dat, 100; channel_selection = x -> .!(endswith.(string.(x), "z")), sample_selection = x -> x.sample .> 1000);

# retrun count of extreme values at specific electrodes at different thresholds
eegfun.n_extreme_value(dat, 100)
eegfun.n_extreme_value(dat, 100, include_extra = true)
eegfun.n_extreme_value(dat, 100, channel_selection = eegfun.channels([:Fp1, :Fp2])) # count extreme values at Fp1 at 100 uV threshold
eegfun.n_extreme_value(dat, 100, channel_selection = x -> endswith.(string.(x), "z"))
eegfun.n_extreme_value(dat, 100, channel_selection = x -> .!(endswith.(string.(x), "z")), sample_selection = x -> x.sample .< 10)



# while preprocessing routine
eegfun.preprocess_eeg_data("pipeline.toml")

# config test
config = eegfun.load_config("pipeline.toml");
eegfun.print_config(config)
eegfun.print_config(config, "config_output.toml")

# read/load layout
layout = eegfun.read_layout("./data/layouts/biosemi72.csv");

# 2D layout
eegfun.polar_to_cartesian_xy!(layout)
fig, ax = eegfun.plot_layout_2d(layout);

eegfun.print_neighbours_dict(neighbours, "electrode_neighbours.toml")
eegfun.plot_layout_2d(layout, neighbours)

# Can set theme for plots
# set_theme!(figure_padding=0)

# Layout with Regions of Interest (ROIs) highlighted
fig, ax = eegfun.plot_layout_2d(layout)
eegfun.add_topo_rois!(ax, layout, [[:PO7, :PO3, :P1], [:PO8, :PO4, :P2]], border_size = 10)
eegfun.add_topo_rois!(ax, layout, [[:PO7, :PO3, :P1], [:PO8, :PO4, :P2]], border_size = 5)
eegfun.add_topo_rois!(ax, layout, [[:Fp1]], border_size = 15, roi_kwargs = Dict(:fill => [true], :fillcolor => [:red], :fillalpha => [0.2]))
eegfun.add_topo_rois!(ax, layout, [[:CPz, :C2, :FCz,  :C1]], border_size = 15, roi_kwargs = Dict(:fill => [true], :fillcolor => [:red], :fillalpha => [0.2]))

# save a basic figure
# NB. for vector graphics, use CairoMakie
save("topo_roi.png", fig)

# 3D layout
eegfun.polar_to_cartesian_xyz!(layout)
neighbours = eegfun.get_electrode_neighbours_xyz(layout, 40);
fig, ax = eegfun.plot_layout_3d(layout)
fig, ax = eegfun.plot_layout_3d(layout, neighbours)

# read data giving BioSemiBDF type
dat = eegfun.read_bdf("../Flank_C_3.bdf");

# plot events
eegfun.plot_events(dat)
eegfun.plot_events_timing(dat)

# create the ContinuousData eeg data type
dat = eegfun.create_eeg_dataframe(dat, layout);
dat.data # DataFrame

# above plot_events_* also works for ContinuousData type
fig, ax = eegfun.plot_events(dat)
fig, ax = eegfun.plot_events_timing(dat)

# Uses vscode viewer if available
# TODO: IDE/code editor agnostic viewer!
eegfun.viewer(dat)
eegfun.head(dat)
eegfun.tail(dat)

# rereference
eegfun.rereference!(dat, :avg)
# rereference!(dat, :mastoid)

# initial high-pass filter to remove slow drifts
eegfun.filter_data!(dat, "hp", 1)

# autodetect vEOG/hEOG signals using simple step detection
# a new Bool column is added to the data frame: :is_vEOG, :is_hEOG
eegfun.detect_eog_onsets!(dat, 50, :vEOG, :is_vEOG)
eegfun.detect_eog_onsets!(dat, 30, :hEOG, :is_hEOG)

# detect extreme values
# a new Bool column is added to the data frame: :is_extreme_value, :is_extreme_value500, :is_extreme_value1000
eegfun.is_extreme_value!(dat, 100);
eegfun.is_extreme_value!(dat, 500, channels = eegfun.channels([:Fp1]), channel_out = :is_extreme_value500);
eegfun.is_extreme_value!(dat, 1000, channels = eegfun.channels_not([:Fp1, :AF3]), channel_out = :is_extreme_value1000);

# count extreme values at specific electrodes at different thresholds
eegfun.n_extreme_value(dat, 100) # count extreme values across all electrodes
eegfun.n_extreme_value(dat, 100, channels = eegfun.channels([:Fp1, :vEOG, :is_extreme_value500])) # count extreme values at Fp1 at 100 uV threshold

# mark trigger windows
eegfun.mark_epoch_windows!(dat, [1, 3], [-0.5, 1.0]) # simple epoch marking with trigger 1 and 3
eegfun.plot_databrowser(dat) # epoch window within extra_channel menu

# or epoch marking using EpochCondition data type (allows potential for more complex epoch marking)
epoch1 = eegfun.EpochCondition(name = "ExampleEpoch1", trigger_sequences = [[1, 3], [3, :any]]) # 1 -> 3 or 3 -> any sequences
eegfun.mark_epoch_windows!(dat, [epoch1], [-2.0, 2.0]) # epoch window within extra_channel menu

# plot databrowser
eegfun.plot_databrowser(dat);
eegfun.plot_databrowser(dat, [dat.layout.label; :vEOG; :hEOG])

   
    




# plot channel summary
fig, ax = eegfun.plot_channel_summary(summary, :range)

# viewer(summary)
# fig, ax = plot_channel_summary(summary, :range)
# # channel_summary(dat, channels(dat)[1:66])
# # channel_summary(dat, 1:5)

# # bad channels
# jp = channel_joint_probability(dat, threshold=5.0, normval=2)
# # jp = channel_joint_probability(dat, threshold=5.0, normval=2, filter_samples = :epoch_window)
# fig, ax = plot_joint_probability(jp)

# cm = correlation_matrix(dat)
# fig, ax = plot_correlation_heatmap(cm)
# cm = correlation_matrix(dat, filter_samples = :epoch_window)
# fig, ax = plot_correlation_heatmap(cm)


# # save / load
# save_object("$(subject)_continuous.jld2", dat)
# dat1 = load_object("3_continuous.jld2")

# ICA "continuous" data
# ica_result = run_ica(dat; exclude_samples = [:is_extreme_value], include_samples = [:epoch_window])
ica_result = eegfun.run_ica(dat; exclude_samples = [:is_extreme_value])


# plot ICA components
eegfun.plot_ica_topoplot(ica_result, dat.layout)
eegfun.plot_ica_topoplot(ica_result, dat.layout; use_global_scale = true)
eegfun.plot_ica_topoplot(ica_result, dat.layout, component_selection = eegfun.components(1:15))
eegfun.plot_ica_topoplot(ica_result, dat.layout, component_selection = eegfun.components(1:15); use_global_scale = true)
eegfun.plot_ica_topoplot(ica_result, dat.layout, component_selection = eegfun.components([1,3]))
eegfun.plot_ica_topoplot(ica_result, dat.layout, component_selection = eegfun.components([1,3]);  use_global_scale = true)

eegfun.plot_ica_topoplot(ica_result, dat.layout, component_selection = eegfun.components([1, 3, 5]);
                  use_global_scale = true,
                  colorbar_kwargs = Dict(:colorbar_plot_numbers => [ 2]))
eegfun.plot_ica_topoplot(ica_result, dat.layout, component_selection = eegfun.components([1, 3, 5, 7, 9]); dims = (2, 3),
                  use_global_scale = true,
                  colorbar_kwargs = Dict(:colorbar_plot_numbers => [ 5]))


eegfun.plot_ica_component_activation(dat, ica_result)

eegfun.plot_databrowser(dat, ica_result)










# EPOCHS
epoch_cfg = [eegfun.EpochCondition(name = "ExampleEpoch1", trigger_sequences = [[1]])]
epochs = []
for (idx, epoch) in enumerate(epoch_cfg)
     push!(epochs, eegfun.extract_epochs(dat, idx, epoch, -2, 4))
end

plot_databrowser(epochs[1])
plot_databrowser(epochs[2])
plot_databrowser(epochs[2], ica_result)
plot_epochs(epochs[1], :Cz)

# average epochs
erps = []
for (idx, epoch) in enumerate(epochs)
    push!(erps, eegfun.average_epochs(epochs[idx]))
end

# ERP Plot
plot_erp(erps[1])
plot_erp(erps[1], :Fp1)
plot_erp(erps[1], [:Fp1, :Fp2])
plot_erp(erps[2], [:Fp1, :Fp2, :Cz])
plot_erp(erps[1], [:Fp1, :Fp2], kwargs = Dict(:average_channels => true))

plot_erp(erps[1], [:Fp1, :Fp2], kwargs = Dict(:average_channels => true))



plot_erp(erps[1], [:Fp1, :Fp2], kwargs = Dict(:average_channels => false, :add_topoplot => true))

plot_erp(erps[1], erps[3], [:PO7, :PO8])

plot_erp([erps[1], erps[2],erps[3]], [:PO7, :PO8], kwargs = Dict(:average_channels => true, :add_topoplot => true))


# bad_chans, opt_params = find_bad_channels(epochs[1], AutoRejectParams())
# 
# # Visualize results for a specific epoch
# fig = plot_interpolation_comparison(epochs[1], bad_chans, 11)  # Show epoch 1
 
# epochs_cleaned = remove_bad_epochs(epochs)
# 
# # Epoch Data Browser
# plot_databrowser(epochs[1])
# plot_databrowser(epochs[1], [epochs[1].layout.label; "hEOG"; "vEOG"])
# 
# plot_epochs(epochs[1], [:Fp1])
# plot_epochs(epochs[1], :Fp1)
# # plot_epochs(epochs, "Fp1")
# plot_epochs(epochs[1], ["PO7", "PO8"])
# 
# # ERP Plot
# # f, ax = plot_erp(erp, :Fp1)
# plot_erp(erps[1], ["Fp1", "Fp2"])
# # plot_erp(erp, [:Fp1, :Fp2])
# # plot_erp(erp, ["Fp1", "Fp2"])
# plot_erp(erp, [:Fp1, :Fp2], kwargs = Dict(:yreversed => true))
# 
# # Topoplot
# include("topo.jl")
# plot_topoplot(erp, xlim = [0.3, 0.3], ylim = [-10, 10])
# plot_topoplot(erp; method = :spherical_splines)
# 
# # ERP Image
# plot_erp_image(epochs, :Fp1)
# plot_erp_image(epochs, "Fp1")
# plot_erp_image(epochs, [:Fp1, :Fp2], colorrange = [-50, 50])
# 
# 
# save_object("$(subject)_$(cond)_epochs.jld2", epochs)
# save_object("$(subject)_$(cond)_erp.jld2", erp)
# 
# test_plot_eog_detection(dat, 1000:14000, "vEOG", "is_vEOG")
# test_plot_eog_detection(dat, 1000:4000, "hEOG", "is_hEOG")
# 
# 
# function grand_average_erps(subjects, conditions)
# 
#     file_problem = check_files_exist(subjects, conditions, "erp")
#     if file_problem
#         return
#     end
# 
#     sample_rate = nothing
#     layout = nothing
# 
#     for condition in conditions
# 
#         ind_subject_condition_array = []
# 
#         for subject in subjects
# 
#             # load individual subject data
#             erp = load_object("$(subject)_$(condition)_erp.jld2")
# 
#             # basic data checks to make sure sample layout and sample rate
#             if isnothing(sample_rate)
#                 sample_rate = erp.sample_rate
#             end
#             if isnothing(layout)
#                 layout = erp.layout
#             end
#             if !isnothing(sample_rate)
#                 if sample_rate != erp.sample_rate
#                     throw(DomainError([sample_rate, erp.sample_rate], "sample rates across files do not match!"))
#                 end
#             end
#             if !isnothing(layout)
#                 if layout != erp.layout
#                     throw(DomainError([layout, erp.layout], "layout across files do not match!"))
#                 end
#             end
#             # update sample_rate/layout from current file
#             sample_rate = erp.sample_rate
#             layout = erp.layout
# 
#             # perform subject average
#             for subject in subjects
#                 push!(ind_subject_condition_array, erp.data)
#             end
#             grand_average = reduce(.+, ind_subject_condition_array) ./ length(ind_subject_condition_array)
#             save_object("$(cond)_ga_erp.jld2", grand_average)
# 
#         end
#     end
# end
# grand_average_erps([3, 4], 1)




########################################################################

# # Add this to sandbox.jl
# function debug_analysis_info(dat)
#     println("Analysis Info:")
#     println("  Reference: $(dat.analysis_info.reference)")
#     println("  HP Filter: $(dat.analysis_info.hp_filter)")
#     println("  LP Filter: $(dat.analysis_info.lp_filter)")
# end

# # Then use it:
# debug_analysis_info(dat)
# filter_data!(dat, "hp", "fir", 1)
# debug_analysis_info(dat)





# @btime dat1 = repair_bad_channels(dat, [:Fp1], neighbours)
# polar_to_cartesian_xyz!(layout)
# @btime dat2 = repair_channels_spherical_spline(dat, [:Fp1])
# lines(dat.data[:, :Fp1])
# lines!(dat1.data[:, :Fp1], color = :red)
# lines!(dat2.data[:, :Fp1], color = :green)



# fig, ax = plot_channel_spectrum(dat )
# fig, ax = plot_channel_spectrum(dat,:P2)
# fig, ax = plot_channel_spectrum(dat,[:P2,:P1])

using eegfun
using GLMakie
dat = eegfun.read_bdf("../Flank_C_3.bdf");
layout = eegfun.read_layout("./data/layouts/biosemi72.csv");
dat = eegfun.create_eeg_dataframe(dat, layout);
# preprocessing steps
eegfun.rereference!(dat, :avg)
eegfun.filter_data!(dat, "hp", "fir", 1, order=1)
eegfun.diff_channel!(dat, [:Fp1, :Fp2], [:IO1, :IO2], :vEOG); # vertical EOG = mean(Fp1, Fp2) - mean(IO1, I02)
eegfun.diff_channel!(dat, :F9, :F10, :hEOG);                  # horizontal EOG = F9 - F10
eegfun.detect_eog_onsets!(dat, 50, :vEOG, :is_vEOG)
eegfun.detect_eog_onsets!(dat, 30, :hEOG, :is_hEOG)
eegfun.is_extreme_value!(dat, 500);
# mark trigger windows
eegfun.mark_epoch_windows!(dat, [1, 3, 4, 5], [-1, 1.0]) # simple epoch marking with trigger 1 and 3
# eegfun.plot_databrowser(dat) # epoch window within extra_channel menu
# ICA "continuous" data
ica_result = eegfun.run_ica(dat; sample_selection = eegfun.samples_not(:is_extreme_value))
# eegfun.plot_databrowser(dat)


# TODO: legend
# fig, ax = eegfun.plot_channel_spectrum(dat, channel_selection = eegfun.channels([:Fp1, :Fp2]))
# plot ICA components

eegfun.plot_ica_topoplot(ica_result, dat.layout)

# eegfun.plot_ica_topoplot(ica_result, dat.layout, component_selection = eegfun.components(1:10))
# eegfun.plot_ica_topoplot(ica_result, dat.layout; use_global_scale = true)
# eegfun.plot_ica_topoplot(ica_result, dat.layout, component_selection = eegfun.components(1:15))
# eegfun.plot_ica_topoplot(ica_result, dat.layout, component_selection = eegfun.components(1:15); use_global_scale = true)
# eegfun.plot_ica_topoplot(ica_result, dat.layout, component_selection = eegfun.components([1,3]))
# eegfun.plot_ica_topoplot(ica_result, dat.layout, component_selection = eegfun.components([1,3]);  use_global_scale = true)
# eegfun.plot_ica_topoplot(ica_result, dat.layout, component_selection = eegfun.components([1, 3, 5]); use_global_scale = true, colorbar_kwargs = Dict(:colorbar_plot_numbers => [ 2]))
# eegfun.plot_ica_topoplot(ica_result, dat.layout, component_selection = eegfun.components([1, 3, 5, 7, 9]); dims = (2, 3), use_global_scale = true, colorbar_kwargs = Dict(:colorbar_plot_numbers => [ 5]))

fig, ax = eegfun.plot_ica_component_activation(dat, ica_result)

fig, ax = eegfun.plot_ica_component_spectrum(ica_result, dat, component_selection = eegfun.components(1))

eegfun.plot_databrowser(dat, ica_result)

# dat_ica_removed, removed_activations = remove_ica_components(dat, ica_result, [1])
# dat_ica_reconstructed =  restore_original_data(dat_ica_removed, ica_result, [1], removed_activations)

eog_comps, eog_comps_metrics_df = eegfun.identify_eog_components(ica_result, dat, sample_selection = eegfun.samples_not(:is_extreme_value))

ecg_comps, ecg_comps_metrics_df = eegfun.identify_ecg_components(ica_result, dat, sample_selection = eegfun.samples_not(:is_extreme_value))

line_noise_comps, line_noise_comps_metrics_df = eegfun.identify_line_noise_components(ica_result, dat)
channel_noise_comps, channel_noise_comps_metrics_df = eegfun.identify_spatial_kurtosis_components(ica_result, dat)


# Method 1: Combine existing results
artifacts = eegfun.combine_artifact_components(
    eog_comps,
    ecg_comps,
    line_noise_comps,
    channel_noise_comps
)

all_comps = eegfun.get_all_ica_components(artifacts)

dat_ica_removed, ica_result_updated = eegfun.remove_ica_components(dat, ica_result, component_selection = eegfun.components(all_comps))
dat_ica_reconstructed, ica_result_restored = eegfun.restore_ica_components(dat_ica_removed, ica_result_updated, component_selection = eegfun.components(all_comps))
dat.data ≈ dat_ica_reconstructed.data


eegfun.plot_databrowser(dat_ica_removed)


fig, ax = eegfun.plot_eog_component_features(eog_comps, eog_comps_metrics_df)
fig, ax = eegfun.plot_ecg_component_features_(ecg_comps, ecg_comps_metrics_df)
fig, ax = eegfun.plot_line_noise_components(line_noise_comps, line_noise_comps_metrics_df)
fig, ax = eegfun.plot_spatial_kurtosis_components(channel_noise_comps, channel_noise_comps_metrics_df)

