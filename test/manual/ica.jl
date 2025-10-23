using eegfun
using GLMakie
# using BenchmarkTools

# Get some basic data with initial preprocessing steps (high-pass filter, epoch)
data_file = joinpath(@__DIR__, "..", "..", "..", "Flank_C_3.bdf")
layout_file = eegfun.read_layout("./data/layouts/biosemi/biosemi72.csv");
eegfun.polar_to_cartesian_xy!(layout_file)
dat = eegfun.read_bdf(data_file);
dat = eegfun.create_eeg_dataframe(dat, layout_file);
eegfun.rereference!(dat, :avg)
# eegfun.filter_data!(dat, "hp", 0.5)
eegfun.filter_data!(dat, "hp", 1)
# eegfun.resample!(dat, 4)
eegfun.is_extreme_value!(dat, 100);

eegfun.channel_difference!(
    dat,
    channel_selection1 = eegfun.channels([:Fp1, :Fp2]),
    channel_selection2 = eegfun.channels([:IO1, :IO2]),
    channel_out = :vEOG,
); # vertical EOG = mean(Fp1, Fp2) - mean(IO1, I02)
eegfun.channel_difference!(
    dat,
    channel_selection1 = eegfun.channels([:F9]),
    channel_selection2 = eegfun.channels([:F10]),
    channel_out = :hEOG,
); # vertical EOG = mean(Fp1, Fp2) - mean(IO1, I02)

# ICA on continuous data
# ica_result = eegfun.run_ica(dat; sample_selection = eegfun.samples_not(:is_extreme_value_100))
ica_result = eegfun.run_ica(dat; sample_selection = eegfun.samples_not(:is_extreme_value_100))
ica_result = eegfun.run_ica(dat; sample_selection = eegfun.samples_not(:is_extreme_value_100), percentage_of_data = 25)


fig, ax, analysis_settings = eegfun.plot_databrowser(dat, ica_result)
dat_new = eegfun.apply_analysis_settings(dat, ica_result, analysis_settings)
eegfun.plot_databrowser(dat_new)

eegfun.plot_topography(
    ica_result,
    component_selection = eegfun.components(1:4),
    method = :spherical_spline,
    colorbar_plot = true,
    colorbar_position = :right,
    colorbar_vertical = true,
    colorbar_components = [1, 2],
);
eegfun.plot_topography(
    ica_result,
    component_selection = eegfun.components(1:4),
    method = :spherical_spline,
    colorbar_plot = true,
    colorbar_position = :below,
    colorbar_vertical = false,
);


# Test single component plotting (simplified approach)
eegfun.plot_topography(ica_result, method = :multiquadratic);
eegfun.plot_topography(ica_result, method = :spherical_spline);
eegfun.plot_topography(
    ica_result,
    component_selection = eegfun.components(1),
    method = :multiquadratic,
    colorbar_plot = true,
    colorbar_position = :right,
    colorbar_vertical = true,
);
eegfun.plot_topography(
    ica_result,
    component_selection = eegfun.components(1),
    method = :multiquadratic,
    colorbar_plot = true,
    colorbar_position = :below,
    colorbar_vertical = false,
);

# Test multiple components (original approach)
eegfun.plot_topography(ica_result, component_selection = eegfun.components(1), method = :multiquadratic);
eegfun.plot_topography(
    ica_result,
    component_selection = eegfun.components(1),
    method = :multiquadratic,
    point_plot = true,
);
eegfun.plot_topography(
    ica_result,
    component_selection = eegfun.components(1),
    method = :multiquadratic,
    point_plot = true,
    label_plot = true,
);
eegfun.plot_topography(
    ica_result,
    component_selection = eegfun.components(1),
    method = :multiquadratic,
    point_plot = true,
    label_plot = true,
    point_color = :red,
    point_marker = :x,
    point_markersize = 30,
    head_color = :blue,
    head_linewidth = 5,
    head_radius = 1.1,
);
eegfun.plot_topography(
    ica_result,
    component_selection = eegfun.components(1),
    method = :multiquadratic,
    point_plot = true,
    label_plot = true,
    point_color = :red,
    point_marker = :x,
    point_markersize = 30,
    head_color = :blue,
    head_linewidth = 5,
    colorbar_plot = true,
);
eegfun.plot_topography(
    ica_result,
    component_selection = eegfun.components(1),
    method = :multiquadratic,
    colorbar_plot = true,
);
eegfun.plot_topography(
    ica_result,
    component_selection = eegfun.components(1),
    method = :multiquadratic,
    colorbar_plot = true,
    colorbar_position = :below,
    colorbar_vertical = false,
);
eegfun.plot_topography(
    ica_result,
    component_selection = eegfun.components(1),
    method = :multiquadratic,
    colorbar_plot = true,
);
eegfun.plot_topography(ica_result, component_selection = eegfun.components(1:2), method = :spherical_spline);
eegfun.plot_topography(
    ica_result,
    component_selection = eegfun.components(1:2),
    method = :spherical_spline,
    colorbar_plot = true,
);
eegfun.plot_topography(
    ica_result,
    component_selection = eegfun.components(1:4),
    method = :spherical_spline,
    colorbar_plot = true,
);
eegfun.plot_topography(
    ica_result,
    component_selection = eegfun.components(1:4),
    method = :spherical_spline,
    colorbar_plot = true,
    colorbar_position = :below,
    colorbar_vertical = false,
);
eegfun.plot_topography(ica_result, component_selection = eegfun.components(1:4), method = :spherical_spline);
eegfun.plot_topography(
    ica_result,
    component_selection = eegfun.components(1:4),
    method = :spherical_spline,
    dims = (4, 1),
);
eegfun.plot_topography(
    ica_result,
    component_selection = eegfun.components(1:4),
    method = :spherical_spline,
    colorbar_plot = false,
);
eegfun.plot_topography(
    ica_result,
    component_selection = eegfun.components(1:4),
    method = :spherical_spline,
    colorbar_plot = true,
    colorbar_plot_numbers = [2],
);
eegfun.plot_topography(
    ica_result,
    component_selection = eegfun.components(1:4),
    method = :spherical_spline,
    colorbar_plot = true,
    colorbar_plot_numbers = [1, 3],
);
eegfun.plot_topography(
    ica_result,
    component_selection = eegfun.components(1:4),
    method = :spherical_spline,
    colorbar_plot = true,
    colorbar_plot_numbers = [1],
);

# plot_topography
eegfun.plot_topography(ica_result, component_selection = eegfun.components(1:10), method = :multiquadratic);
eegfun.plot_topography(ica_result, component_selection = eegfun.components(1:10), method = :multiquadratic);
eegfun.plot_topography(ica_result, component_selection = eegfun.components(1:4), method = :spherical_spline);
eegfun.plot_topography(ica_result)
eegfun.plot_topography(ica_result, colorbar_plot = true, use_global_scale = true);
eegfun.plot_topography(ica_result, component_selection = eegfun.components(1:4), plot_labels = false);
eegfun.plot_topography(
    ica_result,
    component_selection = eegfun.components(1:4),
    plot_labels = false,
    colorbar_plot = true,
);
eegfun.plot_topography(ica_result, component_selection = eegfun.components(1:10), colorbar_plot = true);
eegfun.plot_topography(
    ica_result,
    component_selection = eegfun.components(1:10),
    colorbar_plot = true,
    colorbar_plot_numbers = [5, 10],
);
eegfun.plot_topography(
    ica_result,
    component_selection = eegfun.components(1:10),
    colorbar_plot = true,
    colorbar_plot_numbers = [1, 5, 10],
);
eegfun.plot_topography(ica_result, component_selection = eegfun.components(1:10), use_global_scale = true);
eegfun.plot_topography(ica_result; component_selection = eegfun.components(1:10), method = :spherical_spline)
eegfun.plot_topography(ica_result; component_selection = eegfun.components(1:10), method = :multiquadratic)


eegfun.plot_topography(
    ica_result,
    component_selection = eegfun.components([1, 3, 5]);
    use_global_scale = true,
    colorbar_plot = true,
    colorbar_plot_numbers = [5],
)



eegfun.plot_topography(
    ica_result,
    component_selection = eegfun.components([1, 3, 5, 7, 9]);
    dims = (2, 3),
    use_global_scale = true,
    colorbar_plot = true,
    colorbar_components = [3, 5],
)


# plot_ica_component_activation
eegfun.plot_ica_component_activation(dat, ica_result)
eegfun.plot_ica_component_activation(dat, ica_result, method = :multiquadratic)
eegfun.plot_ica_component_activation(dat, ica_result, method = :spherical_spline)

#################################
# Epoched DataFrameEeg
#################################
# some epoched data
epoch_cfg = [
    eegfun.EpochCondition(name = "ExampleEpoch1", trigger_sequences = [[1]]),
    eegfun.EpochCondition(name = "ExampleEpoch2", trigger_sequences = [[3]]),
]
epochs = eegfun.extract_epochs(dat, epoch_cfg, -2, 4)

# ICA on epoched data
ica_result = eegfun.run_ica(epochs[1]; sample_selection = eegfun.samples_not(:is_extreme_value_100))
ica_result = eegfun.run_ica(epochs; sample_selection = eegfun.samples_not(:is_extreme_value_100))

# ICA Plots
eegfun.plot_topography(ica_result, component_selection = eegfun.components(1:4), method = :multiquadratic);
eegfun.plot_topography(ica_result, component_selection = eegfun.components(1:4), method = :spherical_spline);



eegfun.plot_topography(ica_result)
eegfun.plot_ica_component_activation(dat, ica_result)

eegfun.plot_component_spectrum(ica_result, dat, component_selection = eegfun.components(1))
eegfun.plot_databrowser(dat, ica_result)

# Identify ICA components
eog_comps, eog_comps_metrics_df =
    eegfun.identify_eog_components(dat, ica_result, sample_selection = eegfun.samples_not(:is_extreme_value_100))
ecg_comps, ecg_comps_metrics_df =
    eegfun.identify_ecg_components(dat, ica_result, sample_selection = eegfun.samples_not(:is_extreme_value_100))

line_noise_comps, line_noise_comps_metrics_df = eegfun.identify_line_noise_components(dat, ica_result)
channel_noise_comps, channel_noise_comps_metrics_df = eegfun.identify_spatial_kurtosis_components(dat, ica_result)

# Combine existing results and plot
artifacts = eegfun.combine_artifact_components(eog_comps, ecg_comps, line_noise_comps, channel_noise_comps)
eegfun.plot_artifact_components(ica_result, artifacts; gridscale = 500, plot_points = true)


# Get all identified component artifacts
all_comps = eegfun.get_all_ica_components(artifacts)
dat_ica_removed, ica_result_updated =
    eegfun.remove_ica_components(dat, ica_result, component_selection = eegfun.components(all_comps))

# Reconstruct for sanity check
dat_ica_reconstructed, ica_result_restored = eegfun.restore_ica_components(
    dat_ica_removed,
    ica_result_updated,
    component_selection = eegfun.components(all_comps),
)

# Original should = reconstructed
eegfun.channel_data(dat) ≈ eegfun.channel_data(dat_ica_reconstructed)

# Plot component features
fig, ax = eegfun.plot_eog_component_features(eog_comps, eog_comps_metrics_df)
fig, ax = eegfun.plot_ecg_component_features_(ecg_comps, ecg_comps_metrics_df)
fig, ax = eegfun.plot_line_noise_components(line_noise_comps, line_noise_comps_metrics_df)
fig, ax = eegfun.plot_spatial_kurtosis_components(channel_noise_comps, channel_noise_comps_metrics_df)


GLMakie.closeall()
