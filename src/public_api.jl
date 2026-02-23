# Public API declarations
# The `public` keyword marks names as part of the public API without
# exporting them into the user's namespace (no namespace clash risk).
# Internal functions are prefixed with _ and are not listed here.
# Alphabetically sorted within each category.

# =====================
# Types
# =====================
publicAbstractDataState
publicAnalysisData
publicAnalysisSettings
publicAnalyticResult
publicArtifactComponents
publicArtifactInfo
publicBaselineInfo
publicBatchConfig
publicBatchResult
publicChannelRepairInfo
publicChannelState
publicCluster
publicClusterInfo
publicClusters
publicComponent
publicConfigParameter
publicContinuousData
publicContinuousDataState
publicContinuousRepairInfo
publicDataBrowserState
publicDecodedData
publicDecodingParameters
publicDecodingStatisticsResult
publicEegConfig
publicEegData
publicEegFunData
publicEogConfig
publicEpochData
publicEpochedDataState
publicEpochInfo
publicEpochRejectionInfo
publicEpochRejectionState
publicEpochRepairInfo
publicErpData
publicErpMeasurementsResult
publicExtraChannelInfo
publicFilterConfig
publicFilterInfo
publicFilterSection
publicFilterState
publicIcaComponentState
publicIcaConfig
publicIcaPrms
publicInfoIca
publicInterval
publicLayout
publicLoggingState
publicMarker
publicMasks
publicMultiDataFrameEeg
publicNeighbours
publicNoiseCeiling
publicPermutationDistribution
publicPermutationResult
publicPipelineTemplateOptions
publicPlotHelpInfo
publicPlotLayout
publicPreprocessConfig
publicRejection
publicRsaData
publicSelectionState
publicSharedSelectionState
publicSingleDataFrameEeg
publicSpectrumData
publicStatisticalData
publicStatMatrix
publicStatsResult
publicTemporalCluster
publicTestInfo
publicTFAnalysisData
publicTFAnalyticResult
publicTFCluster
publicTFClusterPermutationResult
publicTFClusters
publicTFMasks
publicTFStatisticalData
publicTFStatMatrix
publicTFStatsResult
publicTimeFreqData
publicTimeFreqEpochData
publicToggleConfig
publicTopoSelectionState
publicTriggerInfo
publicTTestResult
publicUIStyle
publicViewState
publicWorkArrays
publicZScoreRejectionInfo

# =====================
# Data access and utilities
# =====================
publicall_data
publicchannel_data
publicchannel_labels
publicchannels
publicconditions
publicextra_data
publicfile_name
publicgroup_by_condition
publichead
publicmeta_data
publicn_values
publicparticipants
publicsection
publicsubset
publicsubset_bad_data
publictail
publictimes
publicto_data_frame
publicviewer

# =====================
# I/O
# =====================
publiccheck_files_exist
publicfind_file
publicget_files
publicread_all_data
publicread_config
publicread_csv
publicread_data
publicread_eeglab
publicread_fieldtrip
publicread_layout
publicread_raw_data

# =====================
# Layout
# =====================
publicaverage_number_of_neighbours
publiccheck_channel_neighbors
publicclear_neighbours
publiccreate_custom_layout
publiccreate_layout
publicget_neighbours_xy
publicget_neighbours_xyz
publichas_valid_coordinates
publicpolar_to_cartesian_xy
publicpolar_to_cartesian_xyz
publicprint_layout_neighbours
publicsubset_layout
publicsubset_layout!
publicvalidate_layout

# =====================
# Preprocessing
# =====================
publicaverage_epochs
publicbaseline
publicbaseline!
publiccreate_highpass_filter
publiccreate_lowpass_filter
publicdetrend
publicepochs_table
publicextract_epochs
publicfind_times
publicget_filter_characteristics
publichighpass_filter
publichighpass_filter!
publiclog_epochs_table
publiclowpass_filter
publiclowpass_filter!
publicmark_epoch_intervals
publicmark_epoch_intervals!
publicmirror
publicmirror!
publicprint_filter_characteristics
publicreject_epochs
publicreject_epochs!
publicrereference
publicrereference!
publicresample
publicresample!
publicunmirror
publicunmirror!

# =====================
# Channel and condition operations
# =====================
publicchannel_average
publicchannel_average!
publicchannel_delete
publicchannel_delete!
publicchannel_difference
publicchannel_difference!
publicchannel_joint_probability
publicchannel_repairable
publicchannel_repairable!
publicchannel_summary
publiccondition_average
publiccondition_combine
publiccondition_difference
publiccondition_parse_epoch
publiccorrelation_matrix
publiccorrelation_matrix_dual_selection
publiccorrelation_matrix_eog
publicrename_channel
publicrename_channel!

# =====================
# Artifact detection and repair
# =====================
publiccompare_rejections
publicconsecutive
publiccreate_continuous_repair_info
publiccreate_epoch_repair_info
publicdetect_bad_epochs_automatic
publicdetect_bad_epochs_interactive
publicget_rejected_epochs
publicidentify_bad_channels
publicis_extreme_value
publicis_extreme_value!
publicis_step_value
publicis_step_value!
publicn_extreme_value
publicn_step_value
publicrepair_artifacts
publicrepair_artifacts!
publicrepair_artifacts_neighbor
publicrepair_artifacts_neighbor!
publicrepair_artifacts_spherical_spline
publicrepair_artifacts_spherical_spline!
publicrepair_channels
publicrepair_channels_neighbor
publicrepair_channels_per_epoch
publicrepair_channels_spherical
publicsummarize_electrode_repairs
publicunique_rejections

# =====================
# ICA
# =====================
publiccalculate_eog_channels
publiccalculate_eog_channels!
publiccombine_artifact_components
publicdetect_eog_onsets
publicdetect_eog_onsets!
publicdetect_eog_signals
publicdetect_eog_signals!
publicget_all_ica_components
publicget_eog_channels
publicget_selected_components
publicidentify_components
publicidentify_ecg_components
publicidentify_eog_components
publicidentify_line_noise_components
publicidentify_spatial_kurtosis_components
publicinfomax_extended_ica
publicinfomax_ica
publicpartition_channels_by_eog_correlation
publicremove_ica_components
publicrestore_ica_components
publicrun_ica
publicsummarize_ica_components

# =====================
# ERP analysis
# =====================
publicerp_measurements
publicerp_measurements!
publicgfp
publicgfp_and_dissimilarity
publicglobal_dissimilarity
publicgrand_average
publicjackknife_average
publiclrp
publicrealign
publicrealign!

# =====================
# Statistics
# =====================
publicanalytic_test
publiccompute_probability
publicpermutation_test
publicprepare_stats

# =====================
# Time-frequency
# =====================
publicfreq_spectrum
publictf_baseline
publictf_baseline!
publictf_morlet
publictf_multitaper
publictf_stft

# =====================
# Decoding / MVPA
# =====================
publiccreate_work_arrays
publicdecode_libsvm
publiclibsvm_classifier
publicprepare_decoding
publicresample_temporal_data
publictest_against_chance
publictest_against_chance_cluster

# =====================
# RSA
# =====================
publicadd_noise_ceiling
publicadd_noise_ceiling!
publiccompare_models
publiccompute_noise_ceiling
publiccreate_model_rdms
publiccreate_rdm_from_categorical
publiccreate_rdm_from_distances
publiccreate_rdm_from_matrix
publiccreate_rdm_from_reaction_times
publiccreate_rdm_from_similarity_ratings
publiccreate_rdm_from_timeseries
publiccreate_rdm_from_vectors
publiccreate_temporal_model_rdms
publicnormalize_rdm
publicrsa
publicrsa_crossvalidated

# =====================
# Triggers
# =====================
publicsearch_sequence
publictrigger_count

# =====================
# Pipelines and config
# =====================
publicapply_analysis_settings
publicapply_analysis_settings!
publicgenerate_config_template
publicgenerate_pipeline_template
publicpreprocess_v1
publicpreprocess_v2
publicprint_config
publicshow_parameter_info

# =====================
# Plotting
# =====================
publicplot_artifact_components
publicplot_artifact_detection
publicplot_artifact_rejection
publicplot_artifact_repair
publicplot_channel_spectrum
publicplot_channel_summary
publicplot_channel_summary!
publicplot_confusion_matrix
publicplot_correlation_heatmap
publicplot_correlation_heatmap!
publicplot_databrowser
publicplot_decoding
publicplot_ecg_component_features
publicplot_eog_component_features
publicplot_epochs
publicplot_erp
publicplot_erp!
publicplot_erp_filter_gui
publicplot_erp_image
publicplot_erp_measurement_gui
publicplot_erp_measurements
publicplot_erp_stats
publicplot_filter_response
publicplot_frequency_spectrum
publicplot_gfp
publicplot_gui
publicplot_ica_component_activation
publicplot_ica_component_spectrum
publicplot_joint_probability
publicplot_joint_probability!
publicplot_layout_2d
publicplot_layout_2d!
publicplot_layout_3d
publicplot_layout_3d!
publicplot_line_noise_components
publicplot_model_correlations
publicplot_power_spectrum
publicplot_rdm_heatmap
publicplot_rdm_timecourse
publicplot_rsa
publicplot_spatial_kurtosis_components
publicplot_tf_stats
publicplot_time_frequency
publicplot_topography
publicplot_topography!
publicplot_topo_stats
publicplot_trigger_overview
publicplot_trigger_timing

# =====================
# Plot helpers
# =====================
publicadd_topo_rois
publicadd_topo_rois!
publicadd_zscore_columns
publiccombine_boolean_columns
publicget_plot_help_info
publicget_selected_channels
publicget_selected_conditions
publicget_selected_epochs
publicget_selected_samples
publicprint_plot_help
publicshow_plot_help

# =====================
# Logging
# =====================
publicclose_logging
publicsetup_logging

# =====================
# Data creation (testing and demos)
# =====================
publiccreate_batch_test_erp_data
publiccreate_eegfun_data
publiccreate_test_continuous_data
publiccreate_test_continuous_data_empty_triggers
publiccreate_test_continuous_data_with_artifacts
publiccreate_test_continuous_data_with_triggers
publiccreate_test_epoch_data
publiccreate_test_epoch_data_vector
publiccreate_test_epoch_data_with_artifacts
publiccreate_test_epoch_data_with_rt
publiccreate_test_erp_data
publiccreate_test_layout
publiccreate_test_lrp_data
publiccreate_test_summary_data
publiccreate_test_summary_data_with_epochs
publiccreate_test_tf_data

# =====================
# Demos
# =====================
publicgenerate_signal
publicsignal_example_1
publicsignal_example_2
publicsignal_to_data
publicsimulate_erp

# =====================
# Misc
# =====================
publicEegFun_version_info
publicget_package_version
