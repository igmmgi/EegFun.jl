using BioSemiBDF
using CSV
using DSP
using DataFrames
using GLMakie
using JLD2
using LibGEOS
using LinearAlgebra
using OrderedCollections
using Random
using ScatteredInterpolation
using StatsBase
using Printf


include("types.jl")
include("utils.jl")
include("analyse.jl")
include("baseline.jl")
include("channel_difference.jl")
include("epochs.jl")
include("filter.jl")
include("layout.jl")
include("ica.jl")
include("plot.jl")
include("plot_databrowser.jl")
include("plot_events.jl")
include("plot_layout.jl")
include("rereference.jl")
include("topo.jl")
include("utils.jl")
include("viewer.jl")


layout = read_layout("./layouts/biosemi72.csv");
subject = 3
dat = read_bdf("../Flank_C_$(subject).bdf");
dat = create_eeg_dataframe(dat, layout);
filter_data!(dat, "hp", "iir", 1, order=1)
rereference!(dat, :avg)
diff_channel!(dat, [:Fp1, :Fp2], [:IO1, :IO2], :vEOG);
diff_channel!(dat, :F9, :F10, :hEOG);
# autodetect EOG signals
detect_eog_onsets!(dat, 50, :vEOG, :is_vEOG)
detect_eog_onsets!(dat, 30, :hEOG, :is_hEOG)
is_extreme_value!(dat, dat.layout.label, 50);
 plot_databrowser(dat);
# extract epochs

epochs = []
for (idx, epoch) in enumerate([1, 4, 5, 3])
     push!(epochs, extract_epochs(dat, idx, epoch, -2, 4))
end
plot_databrowser(epochs[1])


# plot_databrowser(dat)
# plot_databrowser(dat, [dat.layout.label; :vEOG; :hEOG])

# include("test/runtests.jl")
# test_baseline()
# test_filter()

# using Logging
# # Show all messages
# global_logger(ConsoleLogger(stderr, Logging.Debug))
# # Show only info and above
# global_logger(ConsoleLogger(stderr, Logging.Info))
# # Show only warnings and errors
# global_logger(ConsoleLogger(stderr, Logging.Warn))

# basic layouts
layout = read_layout("./layouts/biosemi72.csv");
layout = read_layout("./layouts/biosemi64.csv");

# 2D layout
polar_to_cartesian_xy!(layout)
plot_layout_2d(layout);

neighbours, nneighbours = get_electrode_neighbours_xy(layout, 80);
plot_layout_2d(layout, neighbours)

# fig, ax = plot_layout_2d(layout)
# add_topo_rois!(ax, layout, [[:PO7, :PO3, :P1], [:PO8, :PO4, :P2]], border_size = 10)
# add_topo_rois!(ax, layout, [[:PO7, :PO3, :P1], [:PO8, :PO4, :P2]], border_size = 5)
# add_topo_rois!(ax, layout, [[:Fp1]], border_size = 5, roi_kwargs = Dict(:fill => [true], :fillcolor => [:red], :fillalpha => [0.2]))
# add_topo_rois!(ax, layout, [[:CPz, :C2, :FCz,  :C1]], border_size = 5, roi_kwargs = Dict(:fill => [true], :fillcolor => [:blue], :fillalpha => [0.2]))

# 3D layout
polar_to_cartesian_xyz!(layout)
neighbours, nneighbours = get_electrode_neighbours_xyz(layout, 40);
plot_layout_3d(layout)
plot_layout_3d(layout, neighbours)


# read bdf file
subject = 3
dat = read_bdf("../Flank_C_$(subject).bdf");
# viewer(dat)

# save / load
# save_object("$(subject)_continuous_raw.jld2", dat)
# dat = load_object("3_continuous_raw.jld2")

# plot_events(dat)

# preprocess the eeg data
dat = create_eeg_dataframe(dat, layout);
plot_databrowser(dat)

plot_events(dat);
# viewer(dat) # requires vscode
# head(dat) # requires vscode
# save_object("$(subject)_continuous_raw_eegfun.jld2", dat)
# dat = load_object("3_continuous_raw.jld2")

# rereference!(dat.data, dat.layout.label, :mastoid)
# rereference!(dat.data, dat.layout.label, [:Fp1])


subject = 3
dat = read_bdf("../Flank_C_$(subject).bdf");
dat = create_eeg_dataframe(dat, layout);
# rereference!(dat, :Fp1)
# filter_data!(dat, "hp", "iir", 1, order=1)
# filter_data!(dat, "lp", "iir", 10, order=2)
# is_extreme_value!(dat, dat.layout.label, 500,  channel_out = :is_extreme_value500);
# is_extreme_value!(dat, dat.layout.label, 1000, channel_out = :is_extreme_value1000);
# plot_databrowser(dat)
# # search for some bad channels
# channel_summary(dat)
# channel_summary(dat, channels(dat))
# channel_summary(dat, channels(dat)[1:66])
# channel_summary(dat, 1:5)
# bad channels
# channel_joint_probability(dat, threshold=5.0, normval=2)
# cm = correlation_matrix(dat)
# plot_correlation_heatmap(cm)
filter_data!(dat, "hp", "fir", 1)
# calculate EOG channels
diff_channel!(dat, [:Fp1, :Fp2], [:IO1, :IO2], :vEOG);
diff_channel!(dat, :F9, :F10, :hEOG);
## # autodetect EOG signals
detect_eog_onsets!(dat, 50, :vEOG, :is_vEOG)
detect_eog_onsets!(dat, 30, :hEOG, :is_hEOG)
is_extreme_value!(dat, dat.layout.label, 500);
plot_databrowser(dat)



# ICA "continuous" data
dat_ica = filter_data(dat, "hp", "iir", 1, order=1)
good_samples = findall(dat_ica.data[!, :is_extreme_value] .== false)
good_channels = setdiff(dat_ica.layout.label, [:PO9])
dat_for_ica = create_ica_data_matrix(dat_ica.data, good_channels, samples_to_include = good_samples)
ica_result = infomax_ica(dat_for_ica, good_channels, n_components = length(good_channels) - 1, params=IcaPrms())

# TODO: head shape size
plot_ica_topoplot(ica_result, dat.layout)
plot_ica_topoplot(ica_result, dat.layout, comps = 1:15)
plot_ica_topoplot(ica_result, dat.layout, comps = [1,3])

# TODO: subplot y size
# TODO: data resetting?
plot_ica_component_activation(dat, ica_result)

dat_ica_removed, removed_activations = remove_ica_components(dat, ica_result, [1])
dat_ica_reconstructed =  restore_original_data(dat_ica_removed, ica_result, [1], removed_activations)


plot_databrowser(dat)
plot_databrowser(dat, :Fp1)
plot_databrowser(dat, [:Fp1, :Fp2])
plot_databrowser(dat, dat.layout.label[1:3])
plot_databrowser(dat, dat.layout.label[[1,3,5]])
plot_databrowser(dat, ica_result)


# extract epochs
epochs = EpochData[]
for (idx, epoch) in enumerate([1, 4, 5, 3])
     push!(epochs, extract_epochs(dat, idx, epoch, -2, 4))
end
plot_databrowser(epochs[1])

# average epochs
erps = []
for (idx, epoch) in enumerate(epochs)
    push!(erps, average_epochs(epochs[idx]))
end

# ERP Plot
plot_erp(erps[1])
plot_erp(erps[1], :Fp1)
plot_erp(erps[1], [:Fp1, :Fp2])
plot_erp(erps[2], [:Fp1, :Fp2, :Cz])
plot_erp(erps[1], [:Fp1, :Fp2], average_channels = true)
plot_erp(erps[1], erps[2], [:PO7, :Fp2])


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
# diff_channel!(dat, "F9", "F10", "hEOG");
# diff_channel!(dat, ["Fp1", "Fp2"], ["IO1", "IO2"], "vEOG");
# 
# 
# detect_eog_onsets!(dat, 50, :vEOG, :is_vEOG)
# detect_eog_onsets!(dat, 30, :hEOG, :is_hEOG)
# 
# test_plot_eog_detection(dat, 1000:14000, "vEOG", "is_vEOG")
# test_plot_eog_detection(dat, 1000:4000, "hEOG", "is_hEOG")
# 
# 
# 
# function test_analysis()
#     for subject = 3:4
#         for condition in [1, 3]
#             println("Reading file: Flank_C_$(subject).bdf")
#             dat = read_bdf("../test_data/Flank_C_$(subject).bdf")
#             dat = create_eeg_dataframe(dat, "/home/ian/Documents/Julia/EEGfun/layouts/biosemi72.csv")
#             filter_data!(dat, "hp", 1, 2)
#             epochs = extract_epochs(dat, condition, -0.5, 2)
#             erp = average_epochs(epochs)
#             save_object("$(subject)_$(condition)_epochs.jld2", epochs)
#             save_object("$(subject)_$(condition)_erp.jld2", erp)
#         end
#     end
# end
# 
# @time test_analysis()
# 
# 
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

# Add this to sandbox.jl
function debug_analysis_info(dat)
    println("Analysis Info:")
    println("  Reference: $(dat.analysis_info.reference)")
    println("  HP Filter: $(dat.analysis_info.hp_filter)")
    println("  LP Filter: $(dat.analysis_info.lp_filter)")
end

# Then use it:
debug_analysis_info(dat)
filter_data!(dat, "hp", "fir", 1)
debug_analysis_info(dat)



















