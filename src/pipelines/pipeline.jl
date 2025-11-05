"""
    preprocess(config::String; log_level::Symbol = :info)

Preprocess EEG data according to the specified configuration file.

# Arguments
- `config::String`: Path to the configuration file in TOML format
- `log_level::Symbol`: Log level for preprocessing (:debug, :info, :warn, :error)
"""
function preprocess(config::String; log_level::Symbol = :info)

    # set up the global log for overall processing
    global_log = setup_global_logging("preprocess_log.txt", log_level = log_level)

    # initialize variable for outer scope
    output_directory = ""
    all_epoch_counts = DataFrame[]  # Vector to store all epoch counts

    try

        # Setup for all analyses files:
        # This involves loading the end-user config file, merging it with the default config, 
        # and creating the PreprocessConfig struct.
        # The epoch conditions files are also loaded and parsed (see XXX for example configuration files)

        @info section("Setup")
        @info "Configuration Files:"
        !isfile(config) && @minimal_error "Config file does not exist: $config"
        cfg = load_config(config)
        cfg == nothing && @minimal_error "Failed to load configuration from: $config"

        # try and merge user config above with default config
        default_config = load_config(joinpath(@__DIR__, "..", "..", "src", "config", "default.toml"))
        default_config == nothing && @minimal_error "Failed to load default configuration"
        cfg = _merge_configs(default_config, cfg)

        # Create the PreprocessConfig object for easier access
        preprocess_cfg = PreprocessConfig(cfg["preprocess"])

        # check if all requested raw data files exist
        raw_data_files = get_files(cfg["files"]["input"]["directory"], cfg["files"]["input"]["raw_data_files"])
        raw_data_files_exist = check_files_exist(raw_data_files)
        !raw_data_files_exist && @minimal_error "Missing raw data files requested within TOML file!"
        @info "Found $(length(raw_data_files)) files: $(join(raw_data_files, ", "))"

        # Read the epoch conditions defined within the toml file (See XXX for examples)
        !isfile(cfg["files"]["input"]["epoch_condition_file"]) &&
            @minimal_error "File missing: $(cfg["files"]["input"]["epoch_condition_file"])"
        epoch_cfgs = condition_parse_epoch(TOML.parsefile(cfg["files"]["input"]["epoch_condition_file"]))
        @info "Loading/parsing epoch file: $(cfg["files"]["input"]["epoch_condition_file"])"

        # Find and load layout file
        layout_file = find_file(cfg["files"]["input"]["layout_file"], joinpath(@__DIR__, "..", "..", "data", "layouts"))
        layout_file === nothing && @minimal_error "Layout file not found: $(cfg["files"]["input"]["layout_file"])"
        layout = read_layout(layout_file)

        # Check if requested output directory exists and if not, create it
        output_directory = cfg["files"]["output"]["directory"]
        !isdir(output_directory) && mkdir(output_directory)
        @info "Output directory: $output_directory"

        # print config to output directory
        print_config(cfg, joinpath(output_directory, "config.toml"))

        # Layout coordinates and calculation of channel neighbours (2D/3D)
        polar_to_cartesian_xy!(layout)
        polar_to_cartesian_xyz!(layout)
        get_layout_neighbours_xy!(layout, preprocess_cfg.neighbour_criterion)
        print_layout_neighbours(layout, joinpath(output_directory, "neighbours.toml"))

        # Actual start of preprocessing pipeline!
        # This is the main loop that processes each raw data file.
        # TODO: embarrasingly parallel? use Threads.@threads?
        # Track processing results
        processed_files = 0
        failed_files = String[]
        for data_file in raw_data_files

            try

                # Individual file processing
                @info section("Processing")
                @info "File: $data_file"

                # Set up per-file logging (temporarily replaces global logger)
                setup_logging(
                    joinpath(output_directory, "$(basename_without_ext(data_file))_log.txt"),
                    log_level = log_level,
                )

                ################### LOAD RAW DATA FILE ###################
                # TODO: update for different file types!
                @info section("Raw Data")
                dat = create_eeg_dataframe(read_bdf(data_file), layout)

                # Save the original data in Julia format
                if cfg["files"]["output"]["save_continuous_data"]
                    @info "Saving continuous data"
                    jldsave(make_output_filename(output_directory, data_file, "_continuous"); data = dat)
                end

                # Mark epoch windows
                # This is useful for x (time/sample) subsetting within the preprocessing pipeline
                @info section("Marking epoch windows")
                @info "Epoch windows: $([preprocess_cfg.epoch_start, preprocess_cfg.epoch_end])"
                mark_epoch_windows!(dat, epoch_cfgs, [preprocess_cfg.epoch_start, preprocess_cfg.epoch_end])
                
                ################### REREFERENCE DATA ###################
                @info section("Rereference")
                rereference!(dat, preprocess_cfg.reference_channel)

                ################### APPLY INITIAL FILTERS ###################
                @info section("Initial Filters")
                @info "Continuous data filters: $(_applied_filters(preprocess_cfg.filter, filter_sections = [:highpass, :lowpass]))"
                filter_data!(dat, preprocess_cfg.filter)

                #################### CALCULATE EOG CHANNELS ###################
                @info section("EOG")
                @info subsection("Calculating EOG (vEOG/hEOG) channels")
                calculate_eog_channels!(dat, preprocess_cfg.eog)

                # Autodetect EOG signals
                @info subsection("Detecting EOG (vEOG/hEOG) onsets")
                detect_eog_signals!(dat, preprocess_cfg.eog)

                # Calculate correlations between all channels and EOG channels
                @info subsection("Channel x vEOG/hEOG Correlation Matrix")
                hEOG_vEOG_cm = eegfun.correlation_matrix_eog(dat, preprocess_cfg.eog)
                eegfun.add_zscore_columns!(hEOG_vEOG_cm)
                log_pretty_table(hEOG_vEOG_cm; title = "Channel x vEOG/hEOG Correlation Matrix (whole dataset)")

                # Calculate correlations between all channels and EOG channels (epoch window)
                @info subsection("Channel x vEOG/hEOG Correlation Matrix (epoch window)")
                hEOG_vEOG_cm_epoch = eegfun.correlation_matrix_eog(dat, preprocess_cfg.eog; sample_selection = samples(:epoch_window))
                eegfun.add_zscore_columns!(hEOG_vEOG_cm_epoch)
                log_pretty_table(hEOG_vEOG_cm_epoch; title = "Channel x vEOG/hEOG Correlation Matrix (epoch window)")

                #################### INITIAL EPOCH and ERP EXTRACTION ###################
                # This is just the initial epoch extraction and not the cleaned epoched data.
                # Here, only the most basic of "preprocessing" is applied (i.e., no ICA artifact correction or 
                # subsequent channel repair)
                @info section("Initial Epoch Extraction")
                epochs_original = extract_epochs(dat, epoch_cfgs, preprocess_cfg.epoch_start, preprocess_cfg.epoch_end)
                erps_original = average_epochs(epochs_original)

                if cfg["files"]["output"]["save_epoch_data_original"]
                    @info "Saving epoch data (original)"
                    jldsave(make_output_filename(output_directory, data_file, "_epochs_all"); data = epochs_original)
                end

                if cfg["files"]["output"]["save_erp_data_original"]
                    @info "Saving ERP data (original)"
                    jldsave(make_output_filename(output_directory, data_file, "_erps_all"); data = erps_original)
                end

                ############################### INITIAL ARTIFACT DETECTION ###############################
                # This is the initial artifact detection on the continuous data and just looks for sections of 
                # data that are extreme (i.e., beyond a certain threshold and unlikely to be real data).
                # Also, try and identify bad channels based on the channel joint probability and z-score variance measures.
                @info section("Artifact Detection: Continuous Data")

                ################### INITIAL CHANNEL SUMMARY ###################
                @info subsection("Channel Summary")
                summary_whole_dataset = channel_summary(dat)
                log_pretty_table(summary_whole_dataset; title = "Channel Summary (whole dataset)")

                summary_epoch_window = channel_summary(dat, sample_selection = samples(:epoch_window))
                log_pretty_table(summary_epoch_window; title = "Channel Summary (epoch window)")

                #################### DETECT EXTREME VALUES IN CONTINUOUS DATA ###################
                @info subsection("Artifact Detection (extreme values)")
                @info "Detecting extreme values: $(preprocess_cfg.eeg.extreme_value_criterion) μV"
                is_extreme_value!(
                    dat,
                    preprocess_cfg.eeg.extreme_value_criterion,
                    channel_out = _flag_symbol("is_extreme_value", preprocess_cfg.eeg.extreme_value_criterion),
                )

                @info subsection("Artifact Detection (criterion values)")
                @info "Detecting artifact values: $(preprocess_cfg.eeg.artifact_value_criterion) μV"
                is_extreme_value!(
                    dat,
                    preprocess_cfg.eeg.artifact_value_criterion,
                    channel_out = _flag_symbol("is_artifact_value", preprocess_cfg.eeg.artifact_value_criterion),
                )

                #################### CHANNEL JOINT PROBABILITY IN CONTINUOUS DATA ###################
                @info subsection("Bad Channel Detection using Channel Joint Probability + Z-Score Variance")
                cjp_whole_dataset = channel_joint_probability(dat)
                log_pretty_table(cjp_whole_dataset; title = "Channel Joint Probability (whole dataset)")

                cjp_epoch_window = channel_joint_probability(dat, sample_selection = samples(:epoch_window))
                log_pretty_table(cjp_epoch_window; title = "Channel Joint Probability (epoch window)")

                @info subsubsection("Bad Channels")
                bad_channels_whole_dataset = identify_bad_channels(summary_whole_dataset, cjp_whole_dataset)
                bad_channels_epoch_window = identify_bad_channels(summary_epoch_window, cjp_epoch_window)
                
                # Separate identification within whole dataset and epoch windows and taking common seems more reliable
                bad_channels = intersect(bad_channels_whole_dataset, bad_channels_epoch_window)
               
                # Some channels may be classified as "bad" due to EOG-related activity.
                # Partition into non-EOG-related (retain) and EOG-related (prob. handled better via ICA later).
                bad_channels_non_eog_related, bad_channels_eog_related = partition_channels_by_eog_correlation(
                    bad_channels, hEOG_vEOG_cm_epoch;
                    eog_channels = [:hEOG, :vEOG], threshold = 0.3, use_z = false,
                )
               
                @info "Bad channels (non-EOG related): $(length(bad_channels_non_eog_related)) channels - $(bad_channels_non_eog_related)"
                @info "Bad channels (EOG related): $(length(bad_channels_eog_related)) channels - $(bad_channels_eog_related)"
                
                # Analyze which channels can be repaired (needed for ICA and repair steps)
                continuous_repair_info = nothing
                if !isempty(bad_channels_non_eog_related)
                    continuous_repair_info = create_continuous_repair_info(:neighbor_interpolation; name = "continuous_repair")
                    channel_repairable!(continuous_repair_info, bad_channels_non_eog_related, layout)
                end

                #################### Independent Component Analysis (ICA) ###################
                # We perform the ica on "continuous" data (clean sections) that usually has a 
                # more extreme high-pass filter applied. 
                # We then run ica on clean sections of "continuous" data
                if preprocess_cfg.ica.apply

                    @info section("ICA")
                    
                    dat_ica = copy(dat) # we need a copy of the data for the ICA
                    
                    if !isnothing(continuous_repair_info) && !isempty(continuous_repair_info.repaired)
                        @info subsection("Removing repairable bad channels for ICA")
                        dat_ica = subset(dat_ica, channel_selection = channels_not(continuous_repair_info.repaired), include_extra = true)
                        @info "Removed $(length(continuous_repair_info.repaired)) repairable channels for ICA: $(continuous_repair_info.repaired)"
                    end

                    # Apply ICA-specific filters
                    @info subsection("ICA filters")
                    @info "ICA data filters: $(_applied_filters(preprocess_cfg.filter, filter_sections = [:ica_highpass, :ica_lowpass]))"
                    filter_data!(dat_ica, preprocess_cfg.filter, filter_sections = [:ica_highpass, :ica_lowpass])

                    @info subsection("Running ICA")
                    ica = run_ica(
                        dat_ica;
                        sample_selection = samples_not(
                            _flag_symbol("is_extreme_value", preprocess_cfg.eeg.extreme_value_criterion),
                        ),
                        percentage_of_data = preprocess_cfg.ica.percentage_of_data,
                    )

                    # Identify all artifact components 
                    @info subsection("Component Identification")
                    component_artifacts, component_metrics = identify_components(
                        dat_ica,
                        ica,
                        sample_selection = samples_not(
                            _flag_symbol("is_extreme_value", preprocess_cfg.eeg.extreme_value_criterion),
                        ),
                    )

                    # Print component metrics to log files
                    log_pretty_table(component_metrics[:eog_metrics]; title = "EOG Component Metrics")
                    log_pretty_table(component_metrics[:ecg_metrics]; title = "ECG Component Metrics")
                    log_pretty_table(component_metrics[:line_noise_metrics]; title = "Line Noise Component Metrics")
                    log_pretty_table(component_metrics[:channel_noise_metrics]; title = "Channel Noise Component Metrics")

                    @info subsection("Removing ICA components")
                    remove_ica_components!(
                        dat,
                        ica,
                        component_selection = components(get_all_ica_components(component_artifacts)),
                    )
                    @info "Removed $(length(get_all_ica_components(component_artifacts))) ICA components"

                    # save ica results
                    if cfg["files"]["output"]["save_ica_data"]
                        @info "Saving ica data"
                        jldsave(make_output_filename(output_directory, data_file, "_ica"); data = ica)
                    end

                end

                #################### REPAIR BAD CHANNELS BEFORE EPOCHING ###################
                if !isnothing(continuous_repair_info)
                    @info section("Channel Repair")
                    
                    # Perform repairs with tracking (similar to repair_artifacts! for epochs)
                    repair_channels!(dat, continuous_repair_info; method = :neighbor_interpolation)
                    
                    @info continuous_repair_info
                end

                #################### RECALCULATE EOG CHANNELS AFTER ICA AND REPAIR ###################
                # After ICA component removal and channel repair, EOG channels need to be recalculated
                # because the underlying channel data has changed
                @info section("EOG Recalculation")
                @info subsection("Recalculating EOG (vEOG/hEOG) channels after ICA and repair")
                calculate_eog_channels!(dat, preprocess_cfg.eog)

                #################### DETECT ARTIFACT VALUES IN CONTINUOUS DATA (FOR EPOCHING) ###################
                @info subsection("Detecting artifact values in continuous data")
                is_extreme_value!(
                    dat,
                    preprocess_cfg.eeg.artifact_value_criterion,
                    channel_out = Symbol(
                        "is_artifact_value" * "_" * string(preprocess_cfg.eeg.artifact_value_criterion),
                    ),
                )

                # Save the original data in Julia format
                if cfg["files"]["output"]["save_continuous_data"]
                    @info "Saving continuous data"
                    jldsave(make_output_filename(output_directory, data_file, "_continuous_cleaned"); data = dat)
                end

                #################### EPOCH EXTRACTION ###################
                @info subsection("Extracting cleaned epoched data")
                epochs = extract_epochs(dat, epoch_cfgs, preprocess_cfg.epoch_start, preprocess_cfg.epoch_end)

                #################### DETECT BAD EPOCHS ###################
                @info subsection("Automatic epoch detection")
                rejection_step1 = detect_bad_epochs_automatic(
                    epochs;
                    z_criterion = 0.0,
                    abs_criterion = preprocess_cfg.eeg.artifact_value_criterion,
                    name = "rejection_step1",
                )
                channel_repairable!(rejection_step1, epochs[1].layout)
                @info rejection_step1
                
                #################### CHANNEL REPAIR PER EPOCH ###################
                # Repair channels identified in rejection_step1 before rejecting epochs
                # This may save epochs that would otherwise be rejected
                @info section("Channel Repair per Epoch")
                repair_artifacts!(epochs, rejection_step1)
                
                #################### RE-DETECT ARTIFACTS AFTER REPAIR ###################
                # Re-detect artifacts after repair to get updated rejection info
                @info subsection("Re-detecting artifacts after repair")
                rejection_info_step2 = detect_bad_epochs_automatic(
                    epochs;
                    z_criterion = 0.0,
                    abs_criterion = preprocess_cfg.eeg.artifact_value_criterion,
                    name = "rejection_step2",
                )
                channel_repairable!(rejection_info_step2, epochs[1].layout)
                @info rejection_info_step2
                
                #################### EPOCH REJECTION ###################
                @info subsection("Rejecting bad epochs")
                epochs = reject_epochs(epochs, rejection_info_step2)
                
                #################### SAVE ARTIFACT INFO ###################
                # Collect all artifact-related info into a single structure
                artifact_info = ArtifactInfo(
                    continuous_repair_info !== nothing ? [continuous_repair_info] : ContinuousRepairInfo[],
                    vcat(rejection_step1, rejection_info_step2),
                )
                jldsave(make_output_filename(output_directory, data_file, "_artifact_info"); data = artifact_info)
                @info "Saved artifact info: $(artifact_info)"

                #################### LOG EPOCH COUNTS AND STORE FOR SUMMARY ###################
                df = log_epochs_table(epochs_original, epochs, title = "Epoch counts per condition (after repair and rejection):")
                push!(all_epoch_counts, df)

                #################### SAVE EPOCH DATA ###################
                if cfg["files"]["output"]["save_epoch_data_cleaned"]
                    @info "Saving epoch data (cleaned)"
                    jldsave(make_output_filename(output_directory, data_file, "_epochs_cleaned"); data = epochs)
                end

                #################### SAVE ERP DATA ###################
                if cfg["files"]["output"]["save_erp_data_cleaned"]
                    erps = average_epochs(epochs)
                    @info "Saving ERP data (cleaned)"
                    jldsave(make_output_filename(output_directory, data_file, "_erps_cleaned"); data = erps)
                end

                @info section("End of Processing")
                @info "Successfully processed $data_file"
                processed_files += 1

            catch e # Log error to both console and global log
                @error "Error processing $data_file" exception=(e, catch_backtrace())
                push!(failed_files, data_file)
            finally # Close per-file logging (restores global logger)
                close_logging()
            end
        end

        #################### FINAL SUMMARY ###################
        @info section("Summary")
        @info "$processed_files success, $(length(failed_files)) fail"
        !isempty(failed_files) && @info "Failed files: $(join(failed_files, ", "))"

        # Print combined epoch counts
        if !isempty(all_epoch_counts)
            combined_counts = vcat(all_epoch_counts...)
            log_pretty_table(combined_counts, title = "Combined epoch counts across all files:")
            jldsave(joinpath(output_directory, "epoch_summary.jld2"); data = combined_counts)
        end

    finally
        close_global_logging()
        log_source = "preprocess_log.txt"
        log_destination = joinpath(output_directory, "preprocess_log.txt")
        if log_source != log_destination
            mv(log_source, log_destination, force = true)
        end
    end

end

