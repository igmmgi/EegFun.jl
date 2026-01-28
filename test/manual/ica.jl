using EegFun
using GLMakie
using BenchmarkTools
# using BenchmarkTools
# Get some basic data with initial preprocessing steps (high-pass filter, epoch)
data_file = joinpath(@__DIR__, "..", "..", "..", "AttentionExp", "Flank_C_17.bdf")
layout_file = EegFun.read_layout("./data/layouts/biosemi/biosemi72.csv");
EegFun.polar_to_cartesian_xy!(layout_file)
dat = EegFun.read_bdf(data_file);
dat = EegFun.create_eeg_dataframe(dat, layout_file);
EegFun.rereference!(dat, :avg)
# EegFun.highpass_filter(dat, "hp", 0.5)
EegFun.highpass_filter!(dat, 1)
# EegFun.resample!(dat, 4)
EegFun.is_extreme_value!(dat, 200);
EegFun.channel_difference!(
    dat,
    channel_selection1 = EegFun.channels([:Fp1, :Fp2]),
    channel_selection2 = EegFun.channels([:IO1, :IO2]),
    channel_out = :vEOG,
); # vertical EOG = mean(Fp1, Fp2) - mean(IO1, I02)
EegFun.channel_difference!(
    dat,
    channel_selection1 = EegFun.channels([:F9]),
    channel_selection2 = EegFun.channels([:F10]),
    channel_out = :hEOG,
); # horizontal EOG = F9 - F10
# ICA on continuous data
# ica_result = EegFun.run_ica(dat; sample_selection = EegFun.samples_not(:is_extreme_value_100))
# ica_result = EegFun.run_ica(dat; sample_selection = EegFun.samples_not(:is_extreme_value_100))
# ica_result = EegFun.run_ica(dat; sample_selection = EegFun.samples_not(:is_extreme_value_200), percentage_of_data = 50)


ica_result_infomax = EegFun.run_ica(dat; sample_selection = EegFun.samples_not(:is_extreme_value_200), percentage_of_data = 50)
# ica_result_infomax_extended = EegFun.run_ica(
#     dat;
#     sample_selection = EegFun.samples_not(:is_extreme_value_200),
#     percentage_of_data = 10,
#     algorithm = :infomax_extended,
# )

EegFun.plot_topography(ica_result_infomax, component_selection = EegFun.components(1:4));


EegFun.plot_ica_component_activation(dat, ica_result_infomax)
EegFun.plot_ica_component_spectrum(dat, ica_result_infomax, component_selection = EegFun.components(1:70))

EegFun.plot_ica_component_activation(dat, ica_result_infomax)

# Calculate components for valid samples
selected_samples = EegFun.get_selected_samples(dat, EegFun.samples_not(:is_extreme_value_200))
components, n_components = EegFun._prepare_ica_data_matrix(dat, ica_result, selected_samples)

crosscov(components[1, 1:2000], components[17, 1:2000])
corspearman(abs.(zscore(components[1, :])), abs.(zscore(components[17, :])))

lp_filter = EegFun.create_filter("lp", "iir", 5.0, dat.sample_rate; order = 3)
comp1 = EegFun.filtfilt(lp_filter.filter_object, components[1, :])
comp2 = EegFun.filtfilt(lp_filter.filter_object, components[17, :])

component_artifacts, component_metrics =
    EegFun.identify_components(dat, ica_result, sample_selection = EegFun.samples_not(:is_extreme_value_200))

EegFun.identify_ecg_components(dat, ica_result; sample_selection = EegFun.samples_not(:is_extreme_value_200))


fig, ax, analysis_settings = EegFun.plot_databrowser(dat, ica_result_infomax)
dat_new = EegFun.apply_analysis_settings(dat, ica_result, analysis_settings)
EegFun.plot_databrowser(dat_new)

EegFun.plot_topography(
    ica_result_infomax,
    component_selection = EegFun.components(1:4),
    method = :spherical_spline,
    colorbar_plot = true,
    colorbar_position = :below,
    colorbar_vertical = false,
    colorbar_plot_numbers = [1, 2],
);

EegFun.plot_topography(
    ica_result_infomax,
    component_selection = EegFun.components(1:4),
    method = :spherical_spline,
    colorbar_plot = true,
    colorbar_position = :below,
    colorbar_vertical = false,
);

# TODO: units and colorbar appropriateness for ICA plots???

# Test single component plotting (simplified approach)
EegFun.plot_topography(ica_result_infomax, method = :multiquadratic);
EegFun.plot_topography(ica_result, method = :spherical_spline);
EegFun.plot_topography(
    ica_result,
    component_selection = EegFun.components(1),
    colorbar_plot = true,
    colorbar_position = :right,
    colorbar_vertical = true,
);
EegFun.plot_topography(
    ica_result,
    component_selection = EegFun.components(1),
    colorbar_plot = true,
    colorbar_position = :below,
    colorbar_vertical = false,
);

# Test multiple components (original approach)
EegFun.plot_topography(ica_result, component_selection = EegFun.components(1), method = :multiquadratic);
EegFun.plot_topography(ica_result, component_selection = EegFun.components(1), point_plot = true);
EegFun.plot_topography(ica_result, component_selection = EegFun.components(1), point_plot = true, label_plot = true);
EegFun.plot_topography(
    ica_result,
    component_selection = EegFun.components(1),
    point_plot = true,
    label_plot = true,
    point_color = :red,
    point_marker = :x,
    point_markersize = 30,
    head_color = :blue,
    head_linewidth = 5,
    head_radius = 1.1,
);

EegFun.plot_topography(
    ica_result_infomax,
    component_selection = EegFun.components(1),
    point_plot = true,
    label_plot = true,
    point_color = :red,
    point_marker = :x,
    point_markersize = 30,
    head_color = :blue,
    head_linewidth = 5,
    colorbar_plot = true,
);

EegFun.plot_topography(ica_result, component_selection = EegFun.components(1), colorbar_plot = true);
EegFun.plot_topography(
    ica_result,
    component_selection = EegFun.components(1),
    colorbar_plot = true,
    colorbar_position = :below,
    colorbar_vertical = false,
);
EegFun.plot_topography(ica_result, component_selection = EegFun.components(1), colorbar_plot = true);
EegFun.plot_topography(ica_result, component_selection = EegFun.components(1:2), method = :spherical_spline);
EegFun.plot_topography(ica_result, component_selection = EegFun.components(1:2), method = :spherical_spline, colorbar_plot = true);
EegFun.plot_topography(ica_result, component_selection = EegFun.components(1:4), method = :spherical_spline, colorbar_plot = true);
EegFun.plot_topography(
    ica_result,
    component_selection = EegFun.components(1:4),
    method = :spherical_spline,
    colorbar_plot = true,
    colorbar_position = :below,
    colorbar_vertical = false,
);
EegFun.plot_topography(ica_result_infomax, component_selection = EegFun.components(1:4), method = :spherical_spline);
EegFun.plot_topography(ica_result, component_selection = EegFun.components(1:4), method = :spherical_spline, dims = (4, 1));
EegFun.plot_topography(ica_result, component_selection = EegFun.components(1:4), method = :spherical_spline, colorbar_plot = false);

EegFun.plot_topography(
    ica_result_infomax,
    component_selection = EegFun.components(1:4),
    method = :spherical_spline,
    colorbar_plot = true,
    colorbar_position = (2, 1),
    colorbar_vertical = false,
);

EegFun.plot_topography(
    ica_result_infomax,
    component_selection = EegFun.components(1:4),
    method = :spherical_spline,
    colorbar_plot = true,
    colorbar_components = [1, 3],
);


EegFun.plot_topography(
    ica_result,
    component_selection = EegFun.components(1:4),
    method = :spherical_spline,
    colorbar_plot = true,
    colorbar_plot_numbers = [1],
);

# plot_topography
EegFun.plot_topography(ica_result, component_selection = EegFun.components(1:10), method = :multiquadratic);
EegFun.plot_topography(ica_result, component_selection = EegFun.components(1:10), method = :multiquadratic);
EegFun.plot_topography(ica_result, component_selection = EegFun.components(1:4), method = :spherical_spline);
EegFun.plot_topography(ica_result)
EegFun.plot_topography(ica_result, colorbar_plot = true, use_global_scale = true);
EegFun.plot_topography(ica_result, component_selection = EegFun.components(1:4), plot_labels = false);
EegFun.plot_topography(ica_result, component_selection = EegFun.components(1:4), plot_labels = false, colorbar_plot = true);
EegFun.plot_topography(ica_result, component_selection = EegFun.components(1:10), colorbar_plot = true);
EegFun.plot_topography(ica_result, component_selection = EegFun.components(1:10), colorbar_plot = true, colorbar_plot_numbers = [5, 10]);
EegFun.plot_topography(ica_result, component_selection = EegFun.components(1:10), colorbar_plot = true, colorbar_plot_numbers = [1, 5, 10]);
EegFun.plot_topography(ica_result, component_selection = EegFun.components(1:10), use_global_scale = true);
EegFun.plot_topography(ica_result; component_selection = EegFun.components(1:10), method = :spherical_spline)
EegFun.plot_topography(ica_result; component_selection = EegFun.components(1:10), method = :multiquadratic)


EegFun.plot_topography(
    ica_result,
    component_selection = EegFun.components([1, 3, 5]);
    use_global_scale = true,
    colorbar_plot = true,
    colorbar_plot_numbers = [5],
)



EegFun.plot_topography(
    ica_result,
    component_selection = EegFun.components([1, 3, 5, 7, 9]);
    dims = (2, 3),
    use_global_scale = true,
    colorbar_plot = true,
    colorbar_components = [3, 5],
)


# plot_ica_component_activation
EegFun.plot_ica_component_activation(dat, ica_result)
EegFun.plot_ica_component_activation(dat, ica_result, method = :multiquadratic)
EegFun.plot_ica_component_activation(dat, ica_result, method = :spherical_spline)

component_artifacts, component_metrics =
    EegFun.identify_components(dat, ica_result, sample_selection = EegFun.samples_not(:is_extreme_value_100))




#################################
# Epoched DataFrameEeg
#################################
# some epoched data
epoch_cfg = [
    EegFun.EpochCondition(name = "ExampleEpoch1", trigger_sequences = [[1]]),
    EegFun.EpochCondition(name = "ExampleEpoch2", trigger_sequences = [[3]]),
]
epochs = EegFun.extract_epochs(dat, epoch_cfg, -2, 4)

# ICA on epoched data
ica_result = EegFun.run_ica(epochs[1]; sample_selection = EegFun.samples_not(:is_extreme_value_100))
ica_result = EegFun.run_ica(epochs; sample_selection = EegFun.samples_not(:is_extreme_value_100))

# ICA Plots
EegFun.plot_topography(ica_result, component_selection = EegFun.components(1:4), method = :multiquadratic);
EegFun.plot_topography(ica_result, component_selection = EegFun.components(1:4), method = :spherical_spline);



EegFun.plot_topography(ica_result)
EegFun.plot_ica_component_activation(dat, ica_result)

EegFun.plot_ica_component_spectrum(dat, ica_result, component_selection = EegFun.components(1))
EegFun.plot_databrowser(dat, ica_result)

# Identify ICA components
eog_comps, eog_comps_metrics_df =
    EegFun.identify_eog_components(dat, ica_result, sample_selection = EegFun.samples_not(:is_extreme_value_100))
ecg_comps, ecg_comps_metrics_df =
    EegFun.identify_ecg_components(dat, ica_result, sample_selection = EegFun.samples_not(:is_extreme_value_100))

line_noise_comps, line_noise_comps_metrics_df = EegFun.identify_line_noise_components(dat, ica_result)
channel_noise_comps, channel_noise_comps_metrics_df = EegFun.identify_spatial_kurtosis_components(ica_result)

# Combine existing results and plot
artifacts = EegFun.combine_artifact_components(eog_comps, ecg_comps, line_noise_comps, channel_noise_comps)
EegFun.plot_artifact_components(ica_result, artifacts; gridscale = 500, plot_points = true)


# Get all identified component artifacts
all_comps = EegFun.get_all_ica_components(artifacts)
dat_ica_removed, ica_result_updated = EegFun.remove_ica_components(dat, ica_result, component_selection = EegFun.components(all_comps))

# Reconstruct for sanity check
dat_ica_reconstructed, ica_result_restored =
    EegFun.restore_ica_components(dat_ica_removed, ica_result_updated, component_selection = EegFun.components(all_comps))

# Original should = reconstructed
EegFun.channel_data(dat) ≈ EegFun.channel_data(dat_ica_reconstructed)

# Plot component features
fig, ax = EegFun.plot_eog_component_features(eog_comps, eog_comps_metrics_df)
fig, ax = EegFun.plot_ecg_component_features_(ecg_comps, ecg_comps_metrics_df)
fig, ax = EegFun.plot_line_noise_components(line_noise_comps, line_noise_comps_metrics_df)
fig, ax = EegFun.plot_spatial_kurtosis_components(channel_noise_comps, channel_noise_comps_metrics_df)


GLMakie.closeall()
