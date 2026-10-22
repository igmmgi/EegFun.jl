# Public API declarations
# The `public` keyword marks names as part of the public API without
# exporting them into the user's namespace (no namespace clash risk).
# Internal functions are prefixed with _ and are not listed here.
# Internal types (state machines, GUI state, buffers) are NOT listed here
# but remain accessible via EegFun.TypeName if needed.
# Alphabetically sorted within each category.

# TODO: seems far too many in here!
# We need to trim this!

# === Types — User-facing data, result, and config types ===
public AnalysisSettings
public AnalyticResult
public ArtifactComponents
public ArtifactInfo
public BaselineInfo
public BatchResult
public CleanLineConfig
public Component
public ContinuousData
public ContinuousRepairInfo
public DecodedData
public DecodingParameters
public DecodingStatisticsResult
public EegConfig
public EegData
public EegFunData
public EegFunError
public EogConfig
public EpochData
public EpochRejectionInfo
public EpochRepairInfo
public ErpData
public ErpMeasurementsResult
public FilesConfig
public FilterConfig
public FilterInfo
public FilterSection
public IcaConfig
public InputFilesConfig
public IcaPrms
public InfoIca
public Interval
public Layout
public LMM_FIT_STATUS
public LmmStatsResult
public MultiDataFrameEeg
public Neighbours
public NoiseCeiling
public OutputFilesConfig
public PermutationResult
public PipelineConfig
public PipelineTemplateOptions
public PreprocessConfig
public Rejection
public ResampleConfig
public RideResult
public RsaData
public SingleDataFrameEeg
public SpectrumData
public StatisticalData
public StatsResult
public TFAnalyticResult
public TFClusterPermutationResult
public TFStatisticalData
public TFStatsResult
public TimeFreqData
public TimeFreqEpochData
public TriggerInfo
public ValidationResult

# === Data access and utilities ===
public all_data
public all_labels
public channel_data
public channel_labels
public channels
public channels_not
public coef
public common_channels
public components
public components_not
public condition_name
public condition_number
public conditions
public conditions_not
public confint
public create_eegfun_data
public duration
public epoch_to_continuous
public nobs
public epochs
public epochs_not
public extra_data
public extra_labels
public file_name
public frequency_vector
public group_by_condition
public has_channels
public head
public meta_data
public meta_labels
public n_channels
public n_epochs
public n_samples
public n_values
public participants
public participants_not
public reference
public sample_rate
public samples
public samples_and
public samples_and_not
public samples_not
public samples_or
public samples_or_not
public section
public stderror
public subset
public subset_bad_data
public tail
public time_vector
public times
public to_data_frame
public viewer

# === I/O ===
public download_eegfun_datasets
public find_file
public get_files
public read_all_data
public read_config
public read_fieldtrip_csv
public read_data
public read_eeglab
public read_fieldtrip
public read_layout
public read_raw_data
public read_xdf

# === Layout ===
public average_number_of_neighbours
public check_channel_neighbors
public clear_neighbours!
public create_custom_layout
public create_layout
public get_neighbours_xy!
public get_neighbours_xyz!
public has_2d_coords
public has_3d_coords
public has_neighbours
public has_valid_coordinates
public n_layout
public neighbour_criterion
public neighbours
public polar_to_cartesian_xy!
public polar_to_cartesian_xyz!
public positions_2D
public positions_3D
public positions_polar
public print_layout_neighbours
public subset_layout
public subset_layout!
public validate_layout

# === Preprocessing ===
public average_epochs
public baseline
public baseline!
public cleanline
public cleanline!
public create_highpass_filter
public create_lowpass_filter
public detrend
public epochs_count
public epochs_table
public extract_epochs
public extract_predictor_stats
public find_times
public get_filter_characteristics
public highpass_filter
public highpass_filter!
public log_epochs_table
public lowpass_filter
public lowpass_filter!
public mark_epoch_intervals!
public mirror
public mirror!
public print_filter_characteristics
public reject_epochs
public reject_epochs!
public rereference
public rereference!
public resample
public resample!
public unmirror
public unmirror!

# === Channel and condition operations ===
public channel_average
public channel_average!
public channel_delete
public channel_delete!
public channel_difference
public channel_difference!
public channel_joint_probability
public channel_repairable
public channel_repairable!
public channel_summary
public condition_average
public condition_combine
public condition_difference
public condition_parse_epoch
public correlation_matrix
public correlation_matrix_dual_selection
public correlation_matrix_eog
public rename_channel
public rename_channel!

# === RSA (Representational Similarity Analysis) ===
public add_noise_ceiling!
public compare_models
public compute_noise_ceiling
public create_model_rdms
public create_rdm
public create_temporal_model_rdms
public normalize_rdm
public rsa
public rsa_crossvalidated

# === Artifact detection and repair ===
public compare_rejections
public compute_csd
public compute_csd!
public consecutive
public create_continuous_repair_info
public create_epoch_repair_info
public detect_bad_epochs_automatic
public detect_bad_epochs_interactive
public find_bridged_channels
public get_rejected_epochs
public identify_bad_channels
public is_extreme_value
public is_extreme_value!
public is_flatline
public is_flatline!
public is_peak_to_peak
public is_peak_to_peak!
public is_step_value
public is_step_value!
public n_extreme_value
public n_flatline
public n_peak_to_peak
public n_step_value
public repair_artifacts
public repair_artifacts!
public repair_artifacts_neighbor!
public repair_artifacts_spherical_spline!
public repair_channels
public repair_channels!
public repair_channels_neighbor!
public summarize_electrode_repairs
public unique_rejections

# === ICA ===
public add_ica_components
public add_ica_components!
public amica_ica
public calculate_eog_channels!
public combine_artifact_components
public detect_eog_onsets!
public get_all_components
public get_all_ica_components
public get_eog_channels
public get_selected_components
public identify_components
public identify_ecg_components
public identify_eog_components
public identify_line_noise_components
public identify_spatial_kurtosis_components
public infomax_extended_ica
public infomax_ica
public partition_channels_by_eog_correlation
public picard_ica
public run_ica
public subtract_ica_components
public subtract_ica_components!
public summarize_ica_components

# === ERP analysis and RIDE ===
public calculate_trigger_interval
public calculate_trigger_interval!
public erp_measurements
public erp_measurements!
public extract_ride
public find_nearest_peak
public gfp
public gfp_and_dissimilarity
public global_dissimilarity
public grand_average
public jackknife_average
public lrp
public realign
public realign!
public ride_decompose
public shift_matrix
public shift_matrix!
public shift_signal
public shift_signal!
public woody_filter

# === Statistics ===
public analytic_test
public fit_mass_lmm
public generate_permutation_matrix
public permutation_test
public prepare_stats

# === Time-frequency ===
public freq_spectrum
public tf_baseline
public tf_baseline!
public tf_morlet
public tf_multitaper
public tf_stft

# === Decoding / MVPA ===
public decode_libsvm
public libsvm_classifier
public prepare_decoding
public resample_temporal_data
public test_against_chance
public test_against_chance_cluster

# === Triggers ===
public allcomb
public search_sequence
public trigger_count
public trigger_info

# === Pipelines and config ===
public apply_analysis_settings
public apply_analysis_settings!
public batch_process
public check_raw_data
public export_bids
public from_dict
public generate_config_template
public generate_pipeline_template
public preprocess
public print_config
public read_bids
public show_parameter_info
public validate_config

# === Plotting ===
public plot_artifact_components
public plot_artifact_detection
public plot_artifact_rejection
public plot_artifact_repair
public plot_channel_spectrum
public plot_channel_summary
public plot_channel_summary!
public plot_confusion_matrix
public plot_correlation_heatmap
public plot_correlation_heatmap!
public plot_databrowser
public plot_decoding
public plot_ecg_component_features
public plot_eog_component_features
public plot_epochs
public plot_erp
public plot_erp!
public plot_erp_errorbar
public plot_erp_errorbar!
public plot_erp_filter_gui
public plot_erp_image
public plot_erp_measurement_gui
public plot_erp_measurements
public plot_erp_stats
public plot_filter_response
public plot_gfp
public plot_gui
public plot_ica_component_activation
public plot_ica_component_spectrum
public plot_joint_probability
public plot_joint_probability!
public plot_layout
public plot_layout!
public plot_layout_2d
public plot_layout_2d!
public plot_layout_3d
public plot_layout_3d!
public plot_line_noise_components
public plot_model_correlations
public plot_rdm_heatmap
public plot_rdm_timecourse
public plot_rsa
public plot_spatial_kurtosis_components
public plot_stat_heatmap
public plot_stat_heatmap!
public plot_tf
public plot_tf_stats
public plot_topography
public plot_topography!
public plot_topography_3d
public plot_topography_stats
public plot_trigger_overview
public plot_trigger_timing
public theme_eegfun

# === Plot helpers (user-facing selection helpers) ===
public add_topo_rois!
public get_selected_channels
public get_selected_conditions
public get_selected_epochs
public get_selected_samples

# === Logging ===
public close_global_logging
public close_logging
public setup_global_logging
public setup_logging

# === Synthetic signal generation ===
public generate_signal
public signal_to_data

# === Demos ===
public signal_example_composition
public signal_example_convolution
public signal_example_decoding
public signal_example_dotproduct
public signal_example_ecg_einthoven
public signal_example_ecg_hrv
public signal_example_emg
public signal_example_heog
public signal_example_ica_clt
public signal_example_ica_geometry
public signal_example_ica_mixture
public signal_example_ica_optimization
public signal_example_ica_separation
public signal_example_ica_sphering
public signal_example_sampling
public signal_example_tf
public simulate_erp

# === Misc ===
public example_path
public version_info
