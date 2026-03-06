# Cheatsheet

Quick reference for common EEG processing tasks in EegFun.jl.

## Data Import

| I want to... | Function |
| --- | --- |
| Load a BioSemi .bdf file | `EegFun.read_raw_data("file.bdf")` |
| Load a BrainVision file | `EegFun.read_raw_data("file.vhdr")` |
| Load an EEGLAB file | `EegFun.read_eeglab("file.set")` |
| Load a FieldTrip .mat file | `EegFun.read_fieldtrip("file.mat")` |
| Load a channel layout | `EegFun.read_layout("layout.csv")` |
| Create EegFun data structure | `EegFun.create_eegfun_data(raw, layout)` |

## Data Persistence (JLD2)

| I want to... | Function |
| --- | --- |
| Save data to JLD2 | `jldsave("file.jld2"; data = dat)` |
| Load data from JLD2 | `load("file.jld2", "data")` |
| Smart-load EegFun data | `EegFun.read_data("file.jld2")` |
| Batch load all matching files | `EegFun.read_all_data("pattern", "dir/")` |
| Group loaded data by condition | `EegFun.group_by_condition(data)` |

## Preprocessing

| I want to... | Function |
| --- | --- |
| Rereference to average | `EegFun.rereference!(dat, :avg)` |
| Rereference to specific channels | `EegFun.rereference!(dat, [:M1, :M2])` |
| High-pass filter | `EegFun.highpass_filter!(dat, 0.1)` |
| Low-pass filter | `EegFun.lowpass_filter!(dat, 30.0)` |
| Bandpass filter | `EegFun.bandpass_filter!(dat, 0.1, 30.0)` |
| Resample data | `EegFun.resample!(dat, 256)` |
| Apply baseline correction | `EegFun.baseline!(dat, (-0.2, 0.0))` |
| Create channel difference | `EegFun.channel_difference!(dat, channel_selection1=channels([:C3]), channel_selection2=channels([:C4]), channel_out=:laterality)` |
| Average channels into ROI | `EegFun.channel_average!(dat, channel_selections=[channels([:Fz, :Cz, :Pz])], output_labels=[:midline])` |
| Delete channels | `EegFun.channel_delete!(dat, channels([:M1, :M2]))` |
| Mark extreme values | `EegFun.is_extreme_value!(dat, 100)` |
| View trigger counts | `EegFun.trigger_count(dat)` |
| Search trigger sequences | `EegFun.trigger_search(dat, [1, 2])` |

## ICA

| I want to... | Function |
| --- | --- |
| Run ICA decomposition | `EegFun.run_ica(dat)` |
| Identify EOG components | `EegFun.identify_eog_components(dat, ica)` |
| Identify ECG components | `EegFun.identify_ecg_components(dat, ica)` |
| Identify line noise components | `EegFun.identify_line_noise_components(dat, ica)` |
| Remove components | `EegFun.remove_components!(dat, ica, [1, 3])` |
| Restore components | `EegFun.restore_components!(dat, ica, [1])` |

## Epoching

| I want to... | Function |
| --- | --- |
| Define epoch conditions | `cfg = [EegFun.EpochCondition(name="Cond1", trigger_sequences=[[1]])]` |
| Extract epochs | `EegFun.extract_epochs(dat, cfg, (-0.2, 1.0))` |
| Baseline correct epochs | `EegFun.baseline!(epochs, (-0.2, 0.0))` |
| Detect bad epochs automatically | `EegFun.detect_bad_epochs_automatic(epochs)` |
| Detect bad epochs interactively | `EegFun.detect_bad_epochs_interactive(epochs)` |
| Reject marked epochs | `EegFun.reject_epochs(epochs)` |
| Repair artifacts (interpolation) | `EegFun.repair_artifacts(epochs)` |

## ERP Operations

| I want to... | Function |
| --- | --- |
| Average epochs into ERPs | `EegFun.average_epochs(epochs)` |
| Combine conditions | `EegFun.condition_combine(erps, [[1, 2]])` |
| Subtract conditions | `EegFun.condition_difference(erps, [(1, 2)])` |
| Average conditions | `EegFun.condition_average(erps, [[1, 2]])` |
| Compute grand average | `EegFun.grand_average(all_erps)` |
| Jackknife average | `EegFun.jackknife_average(all_erps)` |
| Compute GFP | `EegFun.gfp(erps)` |
| Realign ERPs | `EegFun.realign(erps, :peak)` |
| Compute LRP | `EegFun.lrp(erps, channel_selection1=channels([:C3]), channel_selection2=channels([:C4]))` |

## ERP Measurements & Export

| I want to... | Function |
| --- | --- |
| Extract mean amplitude | `EegFun.erp_measurements("erps", "mean_amplitude", analysis_interval=(0.3, 0.5))` |
| Extract peak amplitude | `EegFun.erp_measurements("erps", "max_peak_amplitude", analysis_interval=(0.3, 0.5))` |
| Extract peak latency | `EegFun.erp_measurements("erps", "max_peak_latency", analysis_interval=(0.3, 0.5))` |
| Explore measurements in GUI | `EegFun.plot_erp_measurement_gui(erps)` |

## Time-Frequency Analysis

| I want to... | Function |
| --- | --- |
| Morlet wavelet decomposition | `EegFun.tf_morlet(epochs, frequencies=4:1:30, n_cycles=3)` |
| Multitaper decomposition | `EegFun.tf_multitaper(epochs, frequencies=4:1:30)` |
| STFT decomposition | `EegFun.tf_stft(epochs, frequencies=4:1:30)` |
| Baseline correct TF data | `EegFun.tf_baseline!(tf, (-0.2, 0.0); method=:db)` |
| TF channel average | `EegFun.channel_average!(tf, channel_selections=[channels([:Fz, :Cz])], output_labels=[:frontal])` |
| TF channel difference | `EegFun.channel_difference!(tf, channel_selection1=channels([:C3]), channel_selection2=channels([:C4]))` |
| TF condition difference | `EegFun.condition_difference(tf_data, [(1, 2)])` |
| TF condition average | `EegFun.condition_average(tf_data, [[1, 2]])` |
| TF grand average | `EegFun.grand_average(all_tf)` |

## Statistics

| I want to... | Function |
| --- | --- |
| Prepare group data for stats | `EegFun.prepare_stats("erps", :paired, condition_selection=conditions([1, 2]))` |
| Run analytic t-test | `EegFun.analytic_test(stat_data)` |
| Run permutation test | `EegFun.permutation_test(stat_data, n_permutations=1000)` |
| Test decoding against chance | `EegFun.test_against_chance(decoded_list)` |
| Cluster permutation (decoding) | `EegFun.test_against_chance_cluster(decoded_list)` |

## Plotting

| I want to... | Function |
| --- | --- |
| Browse continuous data | `EegFun.plot_databrowser(dat)` |
| Browse with ICA | `EegFun.plot_databrowser(dat, ica)` |
| Plot ERPs | `EegFun.plot_erp(erps, channel_selection=channels(:Pz))` |
| Plot topography | `EegFun.plot_topography(erps, interval_selection=times(0.3, 0.5))` |
| Plot ERP statistics | `EegFun.plot_erp_stats(result, channel_selection=channels(:Pz))` |
| Plot TF power | `EegFun.plot_tf(tf, channel_selection=channels(:Cz))` |
| Plot TF statistics | `EegFun.plot_tf_stats(result, channel_selection=channels(:Cz))` |
| Plot TF topography stats | `EegFun.plot_topo_stats(result, freq_range=(8.0, 12.0))` |
| Plot ICA components | `EegFun.plot_topography(ica)` |
| Plot epochs grid | `EegFun.plot_epochs(epochs)` |
| Plot ERP image | `EegFun.plot_erp_image(epochs, channel_selection=channels(:Pz))` |
| Plot GFP | `EegFun.plot_gfp(erps)` |
| Plot channel layout | `EegFun.plot_layout(layout)` |
| Plot artifact detection | `EegFun.plot_artifact_detection(epochs)` |
| Plot frequency spectrum | `EegFun.plot_frequency_spectrum(dat)` |
| Plot filter response | `EegFun.plot_filter(dat)` |
| Plot decoding results | `EegFun.plot_decoding(decoded)` |
| Plot RSA results | `EegFun.plot_rsa(rsa_result)` |

## Selection Helpers

| I want to select... | Helper |
| --- | --- |
| Specific channels | `EegFun.channels([:Fz, :Cz, :Pz])` |
| All channels except | `EegFun.channels_not([:M1, :M2])` |
| Channel range | `EegFun.channels(1:32)` |
| Specific conditions | `EegFun.conditions([1, 2])` |
| All conditions except | `EegFun.conditions_not([3])` |
| Time interval | `EegFun.times((0.3, 0.5))` |
| Specific participants | `EegFun.participants([1, 2, 3])` |
| All participants except | `EegFun.participants_not([10])` |
| Specific samples | `EegFun.samples(:is_extreme_value_100)` |
| All samples except | `EegFun.samples_not(:is_extreme_value_100)` |
| Specific components | `EegFun.components([1, 3, 5])` |
