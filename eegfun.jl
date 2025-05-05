# using Logging
# Show all messages
# global_logger(ConsoleLogger(stderr, Logging.Debug))
# # Show only info and above
# global_logger(ConsoleLogger(stderr, Logging.Info))
# # Show only warnings and errors
# global_logger(ConsoleLogger(stderr, Logging.Warn))

# package
include("src/eegfun.jl")

# load layout
layout = read_layout("./data/layouts/biosemi72.csv");

# 2D layout
polar_to_cartesian_xy!(layout)
fig, ax = plot_layout_2d(layout);
 
# neighbours, nneighbours = get_electrode_neighbours_xy(layout, 80);
# plot_layout_2d(layout, neighbours)
# 
# fig, ax = plot_layout_2d(layout)
# add_topo_rois!(ax, layout, [[:PO7, :PO3, :P1], [:PO8, :PO4, :P2]], border_size = 10)
# add_topo_rois!(ax, layout, [[:PO7, :PO3, :P1], [:PO8, :PO4, :P2]], border_size = 5)
# add_topo_rois!(ax, layout, [[:Fp1]], border_size = 5, roi_kwargs = Dict(:fill => [true], :fillcolor => [:red], :fillalpha => [0.2]))
# add_topo_rois!(ax, layout, [[:CPz, :C2, :FCz,  :C1]], border_size = 5, roi_kwargs = Dict(:fill => [true], :fillcolor => [:blue], :fillalpha => [0.2]))

# 3D layout
# polar_to_cartesian_xyz!(layout)
# neighbours, nneighbours = get_electrode_neighbours_xyz(layout, 40);
# plot_layout_3d(layout)
# plot_layout_3d(layout, neighbours)

subject = 3
dat = read_bdf("../Flank_C_$(subject).bdf");

fig, ax = plot_events(dat)


dat = create_eeg_dataframe(dat, layout);
fig, ax = plot_events(dat)
viewer(dat)
head(dat)
tail(dat)

# rereference
rereference!(dat, :avg)
# rereference!(dat, :mastoid)

# initial high-pass filter to remove slow drifts
filter_data!(dat, "hp", "iir", 0.1, order=1)

# caculate EOG channels
diff_channel!(dat, [:Fp1, :Fp2], [:IO1, :IO2], :vEOG);
diff_channel!(dat, :F9, :F10, :hEOG);

# autodetect EOG signals
detect_eog_onsets!(dat, 50, :vEOG, :is_vEOG)
detect_eog_onsets!(dat, 30, :hEOG, :is_hEOG)

# detect extreme values
is_extreme_value!(dat, dat.layout.label, 100);
# is_extreme_value!(dat, dat.layout.label, 500,  channel_out = :is_extreme_value500);
# is_extreme_value!(dat, dat.layout.label, 1000, channel_out = :is_extreme_value1000);

# n_extreme_value(dat.data, [:Fp1], 100)
# n_extreme_value(dat,  100)

# mark trigger windows
mark_epoch_windows!(dat, [1, 4, 5, 3], [-0.5, 1.0])

# plot databrowser
plot_databrowser(dat);
plot_databrowser(dat, [dat.layout.label; :vEOG; :hEOG])

summary = channel_summary(dat)
summary = channel_summary(dat, filter_samples = :epoch_window)
viewer(summary)
fig, ax = plot_channel_summary(summary, :range)
# channel_summary(dat, channels(dat)[1:66])
# channel_summary(dat, 1:5)

# bad channels
jp = channel_joint_probability(dat, threshold=5.0, normval=2)
# jp = channel_joint_probability(dat, threshold=5.0, normval=2, filter_samples = :epoch_window)
fig, ax = plot_joint_probability(jp)

cm = correlation_matrix(dat)
fig, ax = plot_correlation_heatmap(cm)
cm = correlation_matrix(dat, filter_samples = :epoch_window)
fig, ax = plot_correlation_heatmap(cm)


# # save / load
# save_object("$(subject)_continuous.jld2", dat)
# dat1 = load_object("3_continuous.jld2")

# ICA "continuous" data
# ica_result = run_ica(dat; exclude_samples = [:is_extreme_value], include_samples = [:epoch_window])
ica_result = run_ica(dat; exclude_samples = [:is_extreme_value])


# plot ICA components
plot_ica_topoplot(ica_result, dat.layout)
plot_ica_topoplot(ica_result, dat.layout; use_global_scale = true)
plot_ica_topoplot(ica_result, dat.layout, comps = 1:15)
plot_ica_topoplot(ica_result, dat.layout, comps = 1:15; use_global_scale = true)
plot_ica_topoplot(ica_result, dat.layout, comps = [1,3])
plot_ica_topoplot(ica_result, dat.layout, comps = [1,3];  use_global_scale = true)

plot_ica_topoplot(ica_result, dat.layout, comps = [1, 3, 5];
                  use_global_scale = true,
                  colorbar_kwargs = Dict(:colorbar_plot_numbers => [ 2]))
plot_ica_topoplot(ica_result, dat.layout, comps = [1, 3, 5, 7, 9]; dims = (2, 3),
                  use_global_scale = true,
                  colorbar_kwargs = Dict(:colorbar_plot_numbers => [ 5]))


plot_ica_component_activation(dat, ica_result)

plot_databrowser(dat, ica_result)

dat_ica_removed, removed_activations = remove_ica_components(dat, ica_result, [1])
# dat_ica_reconstructed =  restore_original_data(dat_ica_removed, ica_result, [1], removed_activations)

eye_components = identify_eye_components(ica_result, dat)

# To see which components are related to vertical eye movements
println("Vertical eye components: ", eye_components[:vertical_eye])

# To see which components are related to horizontal eye movements
println("Horizontal eye components: ", eye_components[:horizontal_eye])

# To visualize the correlations
fig = plot_eye_component_correlations(ica_result, dat)




# EPOCHS
epochs = []
for (idx, epoch) in enumerate([1, 4, 5, 3])
     push!(epochs, extract_epochs(dat, idx, epoch, -2, 4))
end
plot_databrowser(epochs[1])
plot_databrowser(epochs[2])
plot_databrowser(epochs[2], ica_result)

plot_epochs(epochs[1], :Cz)

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

