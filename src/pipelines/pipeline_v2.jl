"""
    preprocess_v2(config::String; base_dir::Union{String,Nothing} = nothing, log_level::Symbol = :info)

Preprocess EEG data according to the specified configuration file.

# Arguments
- `config::String`: Path to the configuration file in TOML format
- `base_dir::Union{String,Nothing}`: Base directory for resolving relative paths in TOML file.
  If `nothing` (default), uses the directory of the config file. This allows relative paths
  in your TOML to work relative to where your analysis script is located.
- `log_level::Symbol`: Log level for preprocessing (:debug, :info, :warn, :error)

# Notes
- Relative paths in the TOML config file are resolved relative to `base_dir`
- If `base_dir` is not provided, it defaults to the directory containing the config file
- Absolute paths in the TOML are used as-is
- This allows you to use relative paths in your TOML files that work regardless of where
  you run the script from, as long as the config file is in the same directory as your analysis script
"""
function preprocess_v2(config::String; base_dir::Union{String,Nothing} = nothing, log_level::Symbol = :info)

    # Use config file's directory as base_dir if not provided
    # This makes relative paths in TOML work relative to the analysis script location
    base_dir === nothing && (base_dir = abspath(dirname(abspath(config))))
    @info "Using base directory for relative paths: $base_dir"

    # Generate timestamp for unique log filename
    timestamp = Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")
    log_filename = "preprocess_log_$(timestamp).txt"

    # set up the global log for overall processing
    setup_global_logging(log_filename, log_level = log_level)

    # initialize variable for outer scope
    output_directory = ""
    all_epoch_counts = DataFrame[]  # Vector to store all epoch counts

    try

        @info section("Setup")
        @info "Configuration Files:"
        !isfile(config) && @minimal_error "Config file does not exist: $config"
        cfg = load_config(config)
        cfg === nothing && @minimal_error "Failed to load configuration from: $config"

        # try and merge user config above with default config
        default_config = load_config(joinpath(@__DIR__, "..", "..", "src", "config", "default.toml"))
        default_config === nothing && @minimal_error "Failed to load default configuration"
        cfg = _merge_configs(default_config, cfg)

        # Resolve relative/absolute paths
        resolve_path(path::String) = isabspath(path) ? path : joinpath(base_dir, path)
        
        # Resolve input directory
        input_directory = resolve_path(cfg["files"]["input"]["directory"])
        !isdir(input_directory) && @minimal_error "Input directory does not exist: $input_directory"
        
        # Resolve output directory
        output_directory = resolve_path(cfg["files"]["output"]["directory"])
        !isdir(output_directory) && mkpath(output_directory)
        
        # Resolve epoch condition file path
        epoch_condition_file = resolve_path(cfg["files"]["input"]["epoch_condition_file"])

        # Create the PreprocessConfig object for easier access
        preprocess_cfg = PreprocessConfig(cfg["preprocess"])

        # check if all requested raw data files exist
        raw_data_files = get_files(input_directory, cfg["files"]["input"]["raw_data_files"])
        raw_data_files_exist = check_files_exist(raw_data_files)
        !raw_data_files_exist && @minimal_error "Missing raw data files requested within TOML file!"
        @info "Found $(length(raw_data_files)) files: $(print_vector(basename.(raw_data_files)))"

        # Read the epoch conditions defined within the toml file (See XXX for examples)
        !isfile(epoch_condition_file) &&
            @minimal_error "File missing: $epoch_condition_file"
        epoch_cfgs = condition_parse_epoch(TOML.parsefile(epoch_condition_file))
        @info "Loading/parsing epoch file: $epoch_condition_file"

        # Output directory was already set early, just log it
        @info "Output directory: $output_directory"

        # print config to output directory
        print_config(cfg, joinpath(output_directory, "config.toml"))

        # Actual start of preprocessing pipeline!
        # This is the main loop that processes each raw data file.
        # TODO: embarrasingly parallel? use Threads.@threads?
        # Track processing results
        processed_files = 0
        failed_files = String[]


        for (file_idx, data_file) in enumerate(raw_data_files)
            @info "Processing file $file_idx/$(length(raw_data_files)): $(basename(data_file))"
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
                @info section("Raw Data")
                dat = load_data(data_file)

                # Mark epoch windows
                # This is useful for x (time/sample) subsetting within the preprocessing pipeline
                @info section("Marking epoch windows")
                @info "Epoch windows: $([preprocess_cfg.epoch_start, preprocess_cfg.epoch_end])"
                mark_epoch_windows!(dat, epoch_cfgs, [preprocess_cfg.epoch_start, preprocess_cfg.epoch_end])

                #################### CALCULATE EOG CHANNELS ###################
                @info section("EOG")
                @info subsection("Calculating EOG (vEOG/hEOG) channels")
                calculate_eog_channels!(dat, preprocess_cfg.eog)

                # Autodetect EOG signals
                @info subsection("Detecting EOG (vEOG/hEOG) onsets")
                detect_eog_signals!(dat, preprocess_cfg.eog)

                # Calculate correlations between all channels and EOG channels
                @info subsection("Channel x vEOG/hEOG Correlation Matrix")
                hEOG_vEOG_cm = correlation_matrix_eog(dat, preprocess_cfg.eog)
                add_zscore_columns!(hEOG_vEOG_cm)
                log_pretty_table(hEOG_vEOG_cm; title = "Channel x vEOG/hEOG Correlation Matrix (whole dataset)")

                # Calculate correlations between all channels and EOG channels (epoch window)
                @info subsection("Channel x vEOG/hEOG Correlation Matrix (epoch window)")
                hEOG_vEOG_cm_epoch =
                    correlation_matrix_eog(dat, preprocess_cfg.eog; sample_selection = samples(:epoch_window))
                add_zscore_columns!(hEOG_vEOG_cm_epoch)
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
                    jldsave(
                        make_output_filename(output_directory, data_file, "_epochs_original");
                        data = epochs_original,
                    )
                end

                if cfg["files"]["output"]["save_erp_data_original"]
                    @info "Saving ERP data (original)"
                    jldsave(make_output_filename(output_directory, data_file, "_erps_original"); data = erps_original)
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
                @info "Detecting extreme values: $(preprocess_cfg.eeg.extreme_value_abs_criterion) μV"
                is_extreme_value!(
                    dat,
                    preprocess_cfg.eeg.extreme_value_abs_criterion,
                    channel_out = _flag_symbol("is_extreme_value", preprocess_cfg.eeg.extreme_value_abs_criterion),
                )

                @info subsection("Artifact Detection (criterion values)")
                @info "Detecting artifact values: $(preprocess_cfg.eeg.artifact_value_abs_criterion) μV"
                is_extreme_value!(
                    dat,
                    preprocess_cfg.eeg.artifact_value_abs_criterion,
                    channel_out = _flag_symbol("is_artifact_value", preprocess_cfg.eeg.artifact_value_abs_criterion),
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
                    bad_channels,
                    hEOG_vEOG_cm_epoch;
                    eog_channels = [:hEOG, :vEOG],
                    threshold = 0.3,
                    use_z = false,
                )

                @info "Bad channels (non-EOG related): $(length(bad_channels_non_eog_related)) channels - $(bad_channels_non_eog_related)"
                @info "Bad channels (EOG related): $(length(bad_channels_eog_related)) channels - $(bad_channels_eog_related)"

                # Analyze which channels can be repaired (needed for ICA and repair steps)
                continuous_repair_info = nothing
                if !isempty(bad_channels_non_eog_related)
                    continuous_repair_info =
                        create_continuous_repair_info(:neighbor_interpolation; name = "continuous_repair")
                    channel_repairable!(continuous_repair_info, bad_channels_non_eog_related, dat.layout)
                end

                #################### Independent Component Analysis (ICA) ###################
                # We perform the ica on "continuous" data (clean sections) that usually has a 
                # more extreme high-pass filter applied. 
                # We then run ica on clean sections of "continuous" data
                component_artifacts = nothing  # Initialize in case ICA is not applied
                if preprocess_cfg.ica.apply
                    @info section("ICA")

                    dat_ica = copy(dat) # we need a copy of the data for the ICA

                    if !isnothing(continuous_repair_info) && !isempty(continuous_repair_info.repaired)
                        @info subsection("Removing repairable bad channels for ICA")
                        dat_ica = subset(
                            dat_ica,
                            channel_selection = channels_not(continuous_repair_info.repaired),
                            include_extra = true,
                        )
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
                            _flag_symbol("is_extreme_value", preprocess_cfg.eeg.extreme_value_abs_criterion),
                        ),
                        percentage_of_data = preprocess_cfg.ica.percentage_of_data,
                    )

                    # Identify all artifact components 
                    @info subsection("Component Identification")
                    component_artifacts, component_metrics = identify_components(
                        dat, # dat_ica vs. dat makes a difference here! TODO: what is going on?
                        ica,
                        sample_selection = samples_not(
                            _flag_symbol("is_extreme_value", preprocess_cfg.eeg.extreme_value_abs_criterion),
                        ),
                    )

                    # Print component metrics to log files
                    log_pretty_table(component_metrics[:eog_metrics]; title = "EOG Component Metrics")
                    log_pretty_table(component_metrics[:ecg_metrics]; title = "ECG Component Metrics")
                    log_pretty_table(component_metrics[:line_noise_metrics]; title = "Line Noise Component Metrics")
                    log_pretty_table(
                        component_metrics[:channel_noise_metrics];
                        title = "Channel Noise Component Metrics",
                    )

                    @info subsection("Removing ICA components")
                    all_removed_components = get_all_ica_components(component_artifacts)
                    @info "Removed $(length(all_removed_components)) ICA components" component_artifacts
                    remove_ica_components!(dat, ica, component_selection = components(all_removed_components))

                    # save ica results
                    if cfg["files"]["output"]["save_ica_data"]
                        @info "Saving ica data"
                        jldsave(make_output_filename(output_directory, data_file, "_ica"); data = ica)
                    end

                end

                #################### REPAIR BAD CHANNELS BEFORE EPOCHING ###################
                if !isnothing(continuous_repair_info)
                    @info section("Channel Repair")
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
                @info section("Detecting artifact values in continuous data")
                is_extreme_value!(
                    dat,
                    preprocess_cfg.eeg.artifact_value_abs_criterion,
                    channel_out = Symbol(
                        "is_artifact_value" * "_" * string(preprocess_cfg.eeg.artifact_value_abs_criterion),
                    ),
                )

                #################### EPOCH EXTRACTION ###################
                @info section("Extracting cleaned epoched data")
                epochs = extract_epochs(dat, epoch_cfgs, preprocess_cfg.epoch_start, preprocess_cfg.epoch_end)

                # Check if any epochs have empty data
                empty_epochs = [i for (i, ep) in enumerate(epochs) if isempty(ep.data)]
                if !isempty(empty_epochs)
                    @eegfun.minimal_error_throw "Epoch extraction resulted in empty epochs for conditions: $(join([epochs[i].condition_name for i in empty_epochs], ", ")). Check epoch window parameters and trigger locations."
                end

                #################### BASELINE WHOLE EPOCHS ##############
                @info section("Baseline whole epochs")
                baseline!(epochs)

                #################### DETECT BAD EPOCHS ###################
                @info section("Automatic epoch detection")
                rejection_info_step1 = detect_bad_epochs_automatic(
                    epochs;
                    z_criterion = 0.0,
                    abs_criterion = preprocess_cfg.eeg.artifact_value_abs_criterion,
                    name = "rejection_step1",
                )
                channel_repairable!(rejection_info_step1, epochs[1].layout)
                @info "" # formatting
                @info rejection_info_step1

                #################### CHANNEL REPAIR PER EPOCH ###################
                # Repair channels identified in rejection_step1 before rejecting epochs
                # This may save epochs that would otherwise be rejected
                @info section("Channel Repair per Epoch")
                repair_artifacts!(epochs, rejection_info_step1)

                #################### SAVE EPOCH DATA ###################
                if cfg["files"]["output"]["save_epoch_data_cleaned"]
                    @info "Saving epoch data (cleaned)"
                    jldsave(make_output_filename(output_directory, data_file, "_epochs_cleaned"); data = epochs)
                end

                #################### RE-DETECT ARTIFACTS AFTER REPAIR ###################
                # Re-detect artifacts after repair to get updated rejection info
                @info subsection("Re-detecting artifacts after repair")
                rejection_info_step2 = detect_bad_epochs_automatic(
                    epochs;
                    z_criterion = 0.0,
                    abs_criterion = preprocess_cfg.eeg.artifact_value_abs_criterion,
                    name = "rejection_step2",
                )
                channel_repairable!(rejection_info_step2, epochs[1].layout)
                @info "" # formatting
                @info rejection_info_step2

                #################### COMPARE REJECTION STEPS ###################
                @info subsection("Rejection Step Comparison (before vs after repair)")
                rejection_comparison = compare_rejections(rejection_info_step1, rejection_info_step2)
                log_pretty_table(
                    rejection_comparison;
                    title = "Rejection Step Comparison: Effectiveness of Channel Repair",
                )

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

                #################### EPOCH REJECTION ###################
                @info subsection("Rejecting bad epochs")
                epochs = reject_epochs(epochs, rejection_info_step2)

                #################### SAVE ARTIFACT INFO ###################
                # Collect all artifact-related info into a single structure
                @info subsection("Artifact Information")
                artifact_info = ArtifactInfo(
                    continuous_repair_info !== nothing ? [continuous_repair_info] : ContinuousRepairInfo[],
                    vcat(rejection_info_step1, rejection_info_step2),
                    component_artifacts,  # Save ICA components if ICA was applied, otherwise nothing
                )
                jldsave(make_output_filename(output_directory, data_file, "_artifact_info"); data = artifact_info)
                @info "Saved artifact info: $(artifact_info)"

                #################### LOG EPOCH COUNTS AND STORE FOR SUMMARY ###################
                df = log_epochs_table(
                    epochs_original,
                    epochs,
                    title = "Epoch counts per condition (after repair and rejection):",
                )
                push!(all_epoch_counts, df)

                #################### SAVE EPOCH DATA ###################
                if cfg["files"]["output"]["save_epoch_data_good"]
                    @info "Saving epoch data (good)"
                    jldsave(make_output_filename(output_directory, data_file, "_epochs_good"); data = epochs)
                end

                #################### SAVE ERP DATA ###################
                if cfg["files"]["output"]["save_erp_data_good"]
                    erps = average_epochs(epochs)
                    @info "Saving ERP data (good)"
                    jldsave(make_output_filename(output_directory, data_file, "_erps_good"); data = erps)
                end

                @info section("End of Processing")
                @info "Successfully processed file $file_idx/$(length(raw_data_files)): $(basename(data_file))"
                processed_files += 1

            catch e
                @minimal_stacktrace "Error processing $data_file" e 5 # avoid Julia spew!
                push!(failed_files, data_file)
            finally # Close per-file logging (restores global logger)
                close_logging()
            end
        end

        #################### FINAL SUMMARY ###################
        @info section("Summary")
        @info "$processed_files success, $(length(failed_files)) fail"
        !isempty(failed_files) && @info "Failed files: $(join(failed_files, ", "))"

        # Print electrode repair summary (load from saved artifact info files)
        @info subsection("Electrode Repair Summary Across All Participants (Continuous Level Only)")
        electrode_repair_summary = summarize_electrode_repairs("_artifact_info", input_dir = output_directory)
        if !isempty(electrode_repair_summary)
            log_pretty_table(
                electrode_repair_summary;
                title = "Electrode Repairs at Continuous Level: Number of Participants Affected",
            )
        end

        # Print ICA component summary (load from saved artifact info files)
        @info subsection("ICA Component Removal Summary")
        ica_per_file, ica_avg = summarize_ica_components("_artifact_info", input_dir = output_directory)
        if !isempty(ica_per_file)
            log_pretty_table(
                ica_per_file;
                title = "ICA Components Removed per Participant",
                alignment = [:l, :r, :r, :r, :r, :r, :r],  # First column left, rest right
            )
            log_pretty_table(
                ica_avg;
                title = "Average ICA Components Removed per Participant",
                alignment = [:l, :r],  # First column left, second right
            )
        end

        # Print combined epoch counts
        if !isempty(all_epoch_counts)
            epoch_summary, file_summary = _epoch_and_file_summary(all_epoch_counts)
            # Merge with existing summaries (replaces data for files that already exist)
            merged_epoch_summary, merged_file_summary = _merge_summaries(epoch_summary, file_summary, output_directory)
            log_pretty_table(
                merged_epoch_summary,
                title = "Combined epoch counts across all files:",
                alignment = [:l, :r, :l, :r, :r, :r],
            )
            log_pretty_table(
                merged_file_summary,
                title = "Average percentage per condition (averaged across conditions):",
                alignment = [:l, :r],
            )
            @info "Mean percentage (averaged across all conditions and files): $(round(Statistics.mean(merged_file_summary.percentage), digits = 1)) %"
            jldsave(joinpath(output_directory, "epoch_summary.jld2"); data = merged_epoch_summary)
            jldsave(joinpath(output_directory, "file_summary.jld2"); data = merged_file_summary)
        end

    finally
        close_global_logging()
        if !isempty(output_directory) && isdir(output_directory) && isfile(log_filename)
            log_destination = joinpath(output_directory, log_filename)
            if log_filename != log_destination
                mv(log_filename, log_destination, force = true)
            end
        end
    end

end
