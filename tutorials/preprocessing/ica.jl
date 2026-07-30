# Demo: Independent Component Analysis (ICA)
# Shows ICA decomposition on continuous and epoched data, component identification
# (EOG, ECG, line noise, channel noise), component removal/restoration, and visualization.

# Note: EegFun.example_path() resolves bundled example data paths.
# When using your own data, simply pass the file path directly, e.g.:
# dat = EegFun.read_raw_data("/path/to/your/data.bdf")

using EegFun
dat = EegFun.read_raw_data(EegFun.example_path("data/bdf/example1.bdf"));

# read and prepare layout file
layout = EegFun.read_layout(EegFun.example_path("layouts/biosemi/biosemi72.csv"));
EegFun.polar_to_cartesian_xy!(layout)

# create EegFun data structure (EegFun.ContinuousData)
dat = EegFun.create_eegfun_data(dat, layout);

# Some minimal preprocessing (average reference, highpass filter, and detect extreme values)
EegFun.rereference!(dat, :avg)
EegFun.highpass_filter!(dat, 1)
EegFun.is_extreme_value!(dat, 250);

# Calculate EOG signals
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
EegFun.detect_eog_onsets!(dat, 50, :vEOG, :is_vEOG)
EegFun.detect_eog_onsets!(dat, 30, :hEOG, :is_hEOG)

# TODO: 
# picard gpu broken!
# amica gpu seems slow! Actually, still seems to be running on the CPU!

# ICA on continuous data 
@time ica_result = EegFun.run_ica(dat; sample_selection = EegFun.samples_not(:is_extreme_value_250), algorithm = :infomax)
@time ica_result = EegFun.run_ica(dat; sample_selection = EegFun.samples_not(:is_extreme_value_250), algorithm = :infomax, use_gpu = true)

@time ica_result = EegFun.run_ica(dat; sample_selection = EegFun.samples_not(:is_extreme_value_250), algorithm = :infomax_extended)
@time ica_result =
    EegFun.run_ica(dat; sample_selection = EegFun.samples_not(:is_extreme_value_250), algorithm = :infomax_extended, use_gpu = true)

@time ica_result = EegFun.run_ica(dat; sample_selection = EegFun.samples_not(:is_extreme_value_250), algorithm = :picard)
@time ica_result = EegFun.run_ica(dat; sample_selection = EegFun.samples_not(:is_extreme_value_250), algorithm = :picard, use_gpu = true)

@time ica_result = EegFun.run_ica(dat; sample_selection = EegFun.samples_not(:is_extreme_value_250), algorithm = :picard_extended)
@time ica_result =
    EegFun.run_ica(dat; sample_selection = EegFun.samples_not(:is_extreme_value_250), algorithm = :picard_extended, use_gpu = true)

# @time ica_result = EegFun.run_ica(dat; sample_selection = EegFun.samples_not(:is_extreme_value_250), algorithm = :amica)
# @time ica_result = EegFun.run_ica(dat; sample_selection = EegFun.samples_not(:is_extreme_value_250), algorithm = :amica, use_gpu = true)


# Databrowser (here we can turn on/off component removal's)
EegFun.plot_databrowser(dat, ica_result)


# ICA type plots
# Basic Topoplots
EegFun.plot_topography(ica_result, component_selection = EegFun.components(1:20));

# # Component spectra
# EegFun.plot_ica_component_spectrum(dat, ica_result, component_selection = EegFun.components(1:70))
# 
# # Component data/activation
# EegFun.plot_ica_component_activation(dat, ica_result)
# 
# # Databrowser (here we can turn on/off component removal's)
# EegFun.plot_databrowser(dat)
# 
# 
# # Databrowser (here we can turn on/off component removal's)
# EegFun.plot_databrowser(dat, ica_result)
# 
# 
# # Extended ICA
# ica_result = EegFun.run_ica(
#     dat;
#     sample_selection = EegFun.samples_not(:is_extreme_value_200),
#     percentage_of_data = 20,
#     algorithm = :infomax_extended,
# )
# 
# # ICA type plots
# # Basic Topoplots
# EegFun.plot_topography(ica_result, component_selection = EegFun.components(1:20));
# 
# # Component spectra
# EegFun.plot_ica_component_spectrum(dat, ica_result, component_selection = EegFun.components(1:70))
# 
# # Component data/activation
# EegFun.plot_ica_component_activation(dat, ica_result)
# 
# # Databrowser (here we can turn on/off component removal's)
# EegFun.plot_databrowser(dat, ica_result)
# 
# # identify components (default correlation method)
# @time component_artifacts, component_metrics =
#     EegFun.identify_components(dat, ica_result, sample_selection = EegFun.samples_not(:is_extreme_value_200));
# 
# # identify components (Combined method: union of correlation and ICLabel)
# component_artifacts_comb, component_metrics_comb =
#     EegFun.identify_components(dat, ica_result, method = :combined, sample_selection = EegFun.samples_not(:is_extreme_value_200));
# 
# # or individually
# eog_comps, eog_comps_metrics_df =
#     EegFun.identify_eog_components(dat, ica_result, sample_selection = EegFun.samples_not(:is_extreme_value_200));
# 
# ecg_comps, ecg_comps_metrics_df =
#     EegFun.identify_ecg_components(dat, ica_result, sample_selection = EegFun.samples_not(:is_extreme_value_200));
# 
# line_noise_comps, line_noise_comps_metrics_df = EegFun.identify_line_noise_components(dat, ica_result);
# 
# channel_noise_comps, channel_noise_comps_metrics_df = EegFun.identify_spatial_kurtosis_components(ica_result);
# 
# 
# # Get all identified component artifacts
# all_comps = EegFun.get_all_ica_components(component_artifacts)
# dat_ica_removed, ica_result_updated =
#     EegFun.subtract_ica_components(dat, ica_result, component_selection = EegFun.components(all_comps))
# 
# # Reconstruct for sanity check (ie., add components back to data)
# dat_ica_reconstructed, ica_result_restored =
#     EegFun.add_ica_components(dat_ica_removed, ica_result_updated, component_selection = EegFun.components(all_comps))
# 
# # Original should = reconstructed
# EegFun.channel_data(dat) ≈ EegFun.channel_data(dat_ica_reconstructed)
# 
# # Plot component features
# EegFun.plot_eog_component_features(eog_comps, eog_comps_metrics_df)
# EegFun.plot_ecg_component_features(ecg_comps, ecg_comps_metrics_df)
# EegFun.plot_line_noise_components(line_noise_comps, line_noise_comps_metrics_df)
# EegFun.plot_spatial_kurtosis_components(channel_noise_comps, channel_noise_comps_metrics_df)
# 
# 
# #################################
# # Epoched DataFrameEeg
# #################################
# # Create some epoched data
# epoch_cfg = [
#     EegFun.EpochCondition(name = "ExampleEpoch1", trigger_sequences = [[1]]),
#     EegFun.EpochCondition(name = "ExampleEpoch2", trigger_sequences = [[2]]),
# ]
# epochs = EegFun.extract_epochs(dat, epoch_cfg, (-0.2, 1.0))  # -200 to 1000 ms
# 
# # ICA on epoched data
# ica_result = EegFun.run_ica(epochs; sample_selection = EegFun.samples_not(:is_extreme_value_200))
# 
# # ICA type plots
# EegFun.plot_topography(ica_result, component_selection = EegFun.components(1:4));
# EegFun.plot_ica_component_activation(dat, ica_result)
# EegFun.plot_ica_component_spectrum(dat, ica_result, component_selection = EegFun.components(1:70))
# EegFun.plot_databrowser(dat, ica_result)
