"""
    _format_time(seconds::Real)

Format seconds into a human-readable string (e.g., "1h 5m" or "45.2s").
"""
function _format_time(seconds::Real)
    if seconds < 60
        return "$(round(seconds, digits=1))s"
    elseif seconds < 3600
        mins = floor(Int, seconds / 60)
        secs = floor(Int, seconds % 60)
        return "$(mins)m $(secs)s"
    else
        hrs = floor(Int, seconds / 3600)
        mins = floor(Int, (seconds % 3600) / 60)
        return "$(hrs)h $(mins)m"
    end
end

"""
    preprocess(config::String; base_dir=nothing, log_level=:info, skip_existing=false, dry_run=false)

Preprocess EEG data according to the specified configuration file.

# Arguments
- `config::String`: Path to the configuration file in TOML format
- `base_dir::Union{String,Nothing}`: Base directory for resolving relative paths in TOML file.
  If `nothing` (default), uses the directory of the config file. This allows relative paths
  in your TOML to work relative to where your analysis script is located.
- `log_level::Symbol`: Log level for preprocessing (:debug, :info, :warn, :error)
- `skip_existing::Bool`: If `true`, skip files whose final output (`_epochs_final.jld2`)
  already exists in the output directory. Useful for resuming after a crash (default: `false`).
- `dry_run::Bool`: If `true`, validate the configuration, resolve all paths, check that all
  input files exist, and report what *would* happen — without actually processing any data.
  This also serves as a standalone config validator (default: `false`).

# Notes
- Relative paths in the TOML config file are resolved relative to `base_dir`
- If `base_dir` is not provided, it defaults to the directory containing the config file
- Absolute paths in the TOML are used as-is
- Set `reference_channel = "none"` in the TOML to skip rereferencing (e.g., for data
  recorded with an implicit reference that has already been handled)
- Filter sections with `apply = false` are skipped automatically
- Interactive mode (`interactive_continuous`, `interactive_ica`, `interactive_epochs`) will pause execution for *each* file. It is recommended to use interactive mode only for single-file pilot runs to fine-tune parameters, rather than batch processing across a full dataset.
"""
function preprocess(config::String; base_dir::Union{String,Nothing} = nothing, log_level::Symbol = :info, skip_existing::Bool = false, dry_run::Bool = false)

    # Use config file's directory as base_dir if not provided
    # This makes relative paths in TOML work relative to the analysis script location
    isnothing(base_dir) && (base_dir = abspath(dirname(abspath(config))))
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

        # Setup for all analyses files:
        # This involves loading the end-user config file, merging it with the default config, 
        # and creating the PreprocessConfig struct.
        @info section("Setup")
        @info "Configuration Files:"
        cfg = read_config(config)
        isnothing(cfg) && @minimal_error "Failed to load configuration from: $config"

        # Resolve relative paths in config relative to base_dir
        resolve_path(path::String) = isabspath(path) ? path : joinpath(base_dir, path)

        # Resolve input directory
        input_directory = resolve_path(cfg["files"]["input"]["directory"])
        !isdir(input_directory) && @minimal_error "Input directory does not exist: $input_directory"

        # Resolve output directory
        output_directory = resolve_path(cfg["files"]["output"]["directory"])
        !isdir(output_directory) && mkpath(output_directory)

        # Resolve layout file path - try resolved path first, then fall back to package layouts
        layout_file_path = resolve_path(cfg["files"]["input"]["layout_file"])
        if !isfile(layout_file_path)
            # Fall back to searching in package layouts directory
            layout_file = find_file(cfg["files"]["input"]["layout_file"], joinpath(@__DIR__, "..", "..", "resources", "layouts"))
            if !isnothing(layout_file)
                layout_file_path = layout_file
            end
        end
        !isfile(layout_file_path) && @minimal_error "Layout file not found: $layout_file_path"
        layout = read_layout(layout_file_path)

        # Resolve epoch condition file path
        epoch_condition_file = resolve_path(cfg["files"]["input"]["epoch_condition_file"])

        # Create the PreprocessConfig object for easier access
        preprocess_cfg = PreprocessConfig(cfg["preprocess"])

        # check if all requested raw data files exist
        raw_data_files =
            get_files(input_directory, cfg["files"]["input"]["raw_data_files"]; recursive = get(cfg["files"]["input"], "recursive", false))
        raw_data_files_exist = check_files_exist(raw_data_files)
        !raw_data_files_exist && @minimal_error "Missing raw data files requested within TOML file!"
        isempty(raw_data_files) &&
            @minimal_error "No files found in '$input_directory' matching pattern '$(cfg["files"]["input"]["raw_data_files"])'. Check the 'directory' and 'raw_data_files' settings in your pipeline.toml."
        @info "Found $(length(raw_data_files)) files: $(_print_vector(basename.(raw_data_files)))"

        # Read the epoch conditions defined within the toml file
        !isfile(epoch_condition_file) && @minimal_error "File missing: $epoch_condition_file"
        epoch_cfgs = condition_parse_epoch(TOML.parsefile(epoch_condition_file))
        @info "Loading/parsing epoch file: $epoch_condition_file"

        # Output directory was already set early, just log it
        @info "Output directory: $output_directory"

        # print config to output directory
        print_config(cfg, joinpath(output_directory, "config.toml"))

        # Layout coordinates and calculation of channel neighbours (2D/3D)
        polar_to_cartesian_xy!(layout)
        polar_to_cartesian_xyz!(layout)
        get_neighbours_xy!(layout, preprocess_cfg.neighbour_criterion)
        print_layout_neighbours(layout, joinpath(output_directory, "neighbours.toml"))

        # ────────────────────────────────────────────────────────
        # DRY RUN: Report what would happen and return early
        # ────────────────────────────────────────────────────────
        if dry_run
            @info section("Dry Run Summary")
            @info "✓ Configuration loaded and validated successfully"
            @info "✓ Input directory exists: $input_directory"
            @info "✓ Output directory: $output_directory"
            @info "✓ Layout file loaded: $layout_file_path ($(length(layout.data.label)) channels)"
            @info "✓ Epoch condition file loaded: $epoch_condition_file ($(length(epoch_cfgs)) condition(s))"
            @info "✓ Found $(length(raw_data_files)) file(s) to process:"
            for (i, f) in enumerate(raw_data_files)
                existing = isfile(_make_output_filename(output_directory, f, "_epochs_final"))
                status = existing ? "(already processed)" : "(pending)"
                @info "  $i. $(basename(f)) $status"
            end
            @info ""
            @info subsection("Preprocessing Settings")
            @info "  Reference: $(preprocess_cfg.reference_channel)"
            @info "  Epoch window: $(preprocess_cfg.epoch_start)s to $(preprocess_cfg.epoch_end)s"
            @info "  Highpass filter: $(preprocess_cfg.filter.highpass.apply ? "$(preprocess_cfg.filter.highpass.freq) Hz" : "off")"
            @info "  Lowpass filter: $(preprocess_cfg.filter.lowpass.apply ? "$(preprocess_cfg.filter.lowpass.freq) Hz" : "off")"
            @info "  ICA: $(preprocess_cfg.ica.apply ? "on ($(preprocess_cfg.ica.percentage_of_data)% of data)" : "off")"
            @info "  CleanLine: $(preprocess_cfg.cleanline.apply ? "on ($(preprocess_cfg.cleanline.line_frequencies) Hz)" : "off")"
            @info "  Resampling: $(preprocess_cfg.resample.apply ? "on ($(preprocess_cfg.resample.target_rate) Hz)" : "off")"
            @info "  Artifact threshold (abs): $(preprocess_cfg.eeg.artifact_value_abs_criterion) μV"
            @info "  Artifact threshold (z): $(preprocess_cfg.eeg.artifact_value_z_criterion == 0 ? "off" : preprocess_cfg.eeg.artifact_value_z_criterion)"
            @info "  Extreme value threshold: $(preprocess_cfg.eeg.extreme_value_abs_criterion) μV"
            @info "  Interactive continuous: $(preprocess_cfg.interactive_continuous)"
            @info "  Interactive ICA: $(preprocess_cfg.interactive_ica)"
            @info "  Interactive epochs: $(preprocess_cfg.interactive_epochs)"
            @info ""
            @info subsection("Output Files")
            save_flags = [
                ("Continuous (raw)",        cfg["files"]["output"]["save_continuous_data_raw"]),
                ("Continuous (corrected)",  cfg["files"]["output"]["save_continuous_data_corrected"]),
                ("ICA data",               cfg["files"]["output"]["save_ica_data"]),
                ("Epochs (uncorrected)",    cfg["files"]["output"]["save_epoch_data_uncorrected"]),
                ("Epochs (unrejected)",     cfg["files"]["output"]["save_epoch_data_unrejected"]),
                ("Epochs (final)",          cfg["files"]["output"]["save_epoch_data_final"]),
                ("ERPs (uncorrected)",      cfg["files"]["output"]["save_erp_data_uncorrected"]),
                ("ERPs (unrejected)",       cfg["files"]["output"]["save_erp_data_unrejected"]),
                ("ERPs (final)",            cfg["files"]["output"]["save_erp_data_final"]),
            ]
            for (label, flag) in save_flags
                @info "  $(flag ? "✓" : "✗") $label"
            end
            @info ""
            @info "Dry run complete. No data was processed."
            return nothing
        end

        # Actual start of preprocessing pipeline!
        # This is the main loop that processes each raw data file.
        # Track processing results
        processed_files = 0
        skipped_files = 0
        failed_files = String[]
        
        # Track timing for ETA calculation
        total_processing_time = 0.0
        n_processed_for_eta = 0

        for (file_idx, data_file) in enumerate(raw_data_files)
            # Skip files that have already been processed (if skip_existing is enabled)
            if skip_existing
                final_output = _make_output_filename(output_directory, data_file, "_epochs_final")
                if isfile(final_output)
                    @info "Skipping file $file_idx/$(length(raw_data_files)): $(basename(data_file)) (already processed)"
                    skipped_files += 1
                    continue
                end
            end
            
            file_start_time = time()
            
            eta_str = "calculating..."
            if n_processed_for_eta > 0
                avg_time = total_processing_time / n_processed_for_eta
                remaining_files = length(raw_data_files) - file_idx + 1
                eta_str = _format_time(avg_time * remaining_files)
            end

            @info "Processing file $file_idx/$(length(raw_data_files)): $(basename(data_file)) (ETA: $eta_str)"

            try

                # Individual file processing
                @info section("Processing")
                @info "File: $data_file"

                # Set up per-file logging (temporarily replaces global logger)
                setup_logging(joinpath(output_directory, "$(basename_without_ext(data_file))_log.txt"), log_level = log_level)

                ################### LOAD RAW DATA FILE ###################
                @info section("Raw Data")
                dat = create_eegfun_data(read_raw_data(data_file), layout)

                # Save the original data in Julia format
                if cfg["files"]["output"]["save_continuous_data_raw"]
                    @info "Saving continuous data (original)"
                    jldsave(_make_output_filename(output_directory, data_file, "_continuous_raw"); data = dat)
                end

                # Mark epoch intervals
                # This is useful for x (time/sample) subsetting within the preprocessing pipeline
                @info section("Marking epoch intervals")
                @info "Epoch intervals: $([preprocess_cfg.epoch_start, preprocess_cfg.epoch_end])"
                mark_epoch_intervals!(dat, epoch_cfgs, [preprocess_cfg.epoch_start, preprocess_cfg.epoch_end])

                ################### APPLY INITIAL FILTERS ###################
                @info section("Initial Filters")
                highpass_filter!(dat, preprocess_cfg.filter)
                lowpass_filter!(dat, preprocess_cfg.filter)

                #################### RESAMPLE ###################
                if preprocess_cfg.resample.apply
                    @info section("Resampling")
                    @info "Resampling data to $(preprocess_cfg.resample.target_rate) Hz"
                    resample!(dat, preprocess_cfg.resample.target_rate)
                end

                ################### CLEANLINE ###################
                if preprocess_cfg.cleanline.apply
                    @info section("Line Noise Removal (CleanLine)")
                    cleanline!(
                        dat;
                        line_frequencies = preprocess_cfg.cleanline.line_frequencies,
                        bandwidth = preprocess_cfg.cleanline.bandwidth,
                        sliding_win_length = preprocess_cfg.cleanline.sliding_win_length,
                        sliding_win_step = preprocess_cfg.cleanline.sliding_win_step,
                        time_bandwidth = preprocess_cfg.cleanline.time_bandwidth,
                        k_tapers = preprocess_cfg.cleanline.k_tapers,
                        p_value = preprocess_cfg.cleanline.p_value,
                        pad = preprocess_cfg.cleanline.pad
                    )
                end

                #################### CALCULATE EOG CHANNELS ###################
                @info section("EOG")
                @info subsection("Calculating EOG (vEOG/hEOG) channels")
                calculate_eog_channels!(dat, preprocess_cfg.eog)

                # Autodetect EOG signals
                @info subsection("Detecting EOG (vEOG/hEOG) onsets")
                detect_eog_onsets!(dat, preprocess_cfg.eog)

                # Calculate correlations between all channels and EOG channels
                @info subsection("Channel x vEOG/hEOG Correlation Matrix")
                hEOG_vEOG_cm = correlation_matrix_eog(dat, preprocess_cfg.eog)
                add_zscore_columns!(hEOG_vEOG_cm)
                log_pretty_table(hEOG_vEOG_cm; title = "Channel x vEOG/hEOG Correlation Matrix (whole dataset)")

                # Calculate correlations between all channels and EOG channels (epoch interval)
                @info subsection("Channel x vEOG/hEOG Correlation Matrix (epoch interval)")
                hEOG_vEOG_cm_epoch = correlation_matrix_eog(dat, preprocess_cfg.eog; sample_selection = samples(:epoch_interval))
                add_zscore_columns!(hEOG_vEOG_cm_epoch)
                log_pretty_table(hEOG_vEOG_cm_epoch; title = "Channel x vEOG/hEOG Correlation Matrix (epoch interval)")



                ############################### INITIAL ARTIFACT DETECTION ###############################
                # This is the initial artifact detection on the continuous data and just looks for sections of 
                # data that are extreme (i.e., beyond a certain threshold and unlikely to be real data).
                # Also, try and identify bad channels based on the channel joint probability and z-score variance measures.
                @info section("Artifact Detection: Continuous Data")

                ################### INITIAL CHANNEL SUMMARY ###################
                @info subsection("Channel Summary")
                summary_whole_dataset = channel_summary(dat)
                log_pretty_table(summary_whole_dataset; title = "Channel Summary (whole dataset)")

                summary_epoch_interval = channel_summary(dat, sample_selection = samples(:epoch_interval))
                log_pretty_table(summary_epoch_interval; title = "Channel Summary (epoch interval)")

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

                cjp_epoch_interval = channel_joint_probability(dat, sample_selection = samples(:epoch_interval))
                log_pretty_table(cjp_epoch_interval; title = "Channel Joint Probability (epoch interval)")

                @info subsubsection("Bad Channels")
                bad_channels_whole_dataset = identify_bad_channels(summary_whole_dataset, cjp_whole_dataset)
                bad_channels_epoch_interval = identify_bad_channels(summary_epoch_interval, cjp_epoch_interval)

                # Separate identification within whole dataset and epoch intervals and taking common seems more reliable
                bad_channels = intersect(bad_channels_whole_dataset, bad_channels_epoch_interval)

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
                    continuous_repair_info = create_continuous_repair_info(:neighbor_interpolation; name = "continuous_repair")
                    channel_repairable!(continuous_repair_info, bad_channels_non_eog_related, dat.layout)
                end

                if preprocess_cfg.interactive_continuous
                    @info section("Interactive Continuous Data Review")
                    @info "Reviewing continuous data before ICA. Close the window to continue."
                    res = plot_databrowser(dat)
                    wait(res.fig.scene)
                    
                    # Apply any manual settings configured in the databrowser
                    apply_analysis_settings!(dat, res.analysis_settings)
                    
                    # Track manually repaired channels for the summary logs
                    manual_repaired = res.analysis_settings[].repaired_channels
                    if !isempty(manual_repaired)
                        @info "Manually repaired channels: $manual_repaired"
                        if isnothing(continuous_repair_info)
                            continuous_repair_info = create_continuous_repair_info(:neighbor_interpolation; name = "continuous_repair_manual")
                        end
                        # Append to the info object without re-repairing, since apply_analysis_settings! did it
                        append!(continuous_repair_info.repaired, manual_repaired)
                        unique!(continuous_repair_info.repaired)
                    end

                    @info "Continuous data review complete."
                end

                #################### REPAIR BAD CHANNELS ###################
                if !isnothing(continuous_repair_info)
                    @info section("Channel Repair")
                    repair_channels!(dat, continuous_repair_info; method = :neighbor_interpolation)
                    @info continuous_repair_info
                end

                ################### REREFERENCE DATA ###################
                if preprocess_cfg.reference_channel != :none
                    @info section("Rereference")
                    rereference!(dat, preprocess_cfg.reference_channel)
                else
                    @info section("Rereference")
                    @info "Skipping rereferencing (reference_channel = none)"
                end

                #################### INITIAL EPOCH and ERP EXTRACTION ###################
                # This is just the initial epoch extraction and not the cleaned epoched data.
                # Here, preprocessing (rereferencing and channel repair) has been applied, but
                # no ICA artifact correction has been performed.
                @info section("Initial Epoch Extraction (Pre-ICA)")
                epochs_uncorrected = extract_epochs(dat, epoch_cfgs, (preprocess_cfg.epoch_start, preprocess_cfg.epoch_end))
                erps_uncorrected = average_epochs(epochs_uncorrected)

                if cfg["files"]["output"]["save_epoch_data_uncorrected"]
                    @info "Saving epoch data (original)"
                    jldsave(_make_output_filename(output_directory, data_file, "_epochs_uncorrected"); data = epochs_uncorrected)
                end

                if cfg["files"]["output"]["save_erp_data_uncorrected"]
                    @info "Saving ERP data (original)"
                    jldsave(_make_output_filename(output_directory, data_file, "_erps_uncorrected"); data = erps_uncorrected)
                end

                #################### Independent Component Analysis (ICA) ###################
                # We perform the ica on "continuous" data (clean sections) that usually has a 
                # more extreme high-pass filter applied. 
                # We then run ica on clean sections of "continuous" data
                component_artifacts = nothing  # Initialize in case ICA is not applied
                if preprocess_cfg.ica.apply
                    @info section("ICA")

                    dat_ica = copy(dat) # we need a copy of the data for the ICA

                    # Apply ICA-specific filters
                    @info subsection("ICA filters")
                    highpass_filter!(dat_ica, preprocess_cfg.filter, section = :ica_highpass)
                    lowpass_filter!(dat_ica, preprocess_cfg.filter, section = :ica_lowpass)

                    extreme_flag = _flag_symbol("is_extreme_value", preprocess_cfg.eeg.extreme_value_abs_criterion)
                    # Exclude both automatic extreme values and any manually marked bad regions from the GUI
                    artifact_flags = [extreme_flag]
                    hasproperty(dat_ica, :selected_region) && push!(artifact_flags, :selected_region)

                    @info subsection("Running ICA")
                    
                    # Calculate rank to prevent ICA crash from rank-deficient data
                    n_repaired = isnothing(continuous_repair_info) ? 0 : length(continuous_repair_info.repaired)
                    n_ref = preprocess_cfg.reference_channel == :avg ? 1 : 0
                    n_comps = length(channel_labels(dat_ica)) - n_repaired - n_ref
                    
                    ica = run_ica(
                        dat_ica;
                        n_components = n_comps,
                        sample_selection = samples_or_not(artifact_flags),
                        percentage_of_data = preprocess_cfg.ica.percentage_of_data,
                    )

                    # Identify all artifact components 
                    @info subsection("Component Identification")
                    component_artifacts, component_metrics = identify_components(
                        dat,
                        ica,
                        sample_selection = samples_or_not(artifact_flags),
                    )

                    # Print component metrics to log files
                    log_pretty_table(component_metrics[:eog_metrics]; title = "EOG Component Metrics")
                    log_pretty_table(component_metrics[:ecg_metrics]; title = "ECG Component Metrics")
                    log_pretty_table(component_metrics[:line_noise_metrics]; title = "Line Noise Component Metrics")
                    log_pretty_table(component_metrics[:channel_noise_metrics]; title = "Channel Noise Component Metrics")

                    if preprocess_cfg.interactive_ica
                        @info section("Interactive ICA Component Review")
                        @info "Use the databrowser ICA menu to review and manually select components to remove. Close the window to continue."
                        res = plot_databrowser(dat, ica)
                        wait(res.fig.scene)
                        
                        # Add user's manually removed components to the "manual" category in component_artifacts
                        manual_removals = res.analysis_settings[].removed_ica_components
                        if !isempty(manual_removals)
                            @info "Manually removed components: $manual_removals"
                            if !haskey(component_artifacts.artifacts, :manual)
                                component_artifacts.artifacts[:manual] = Int[]
                            end
                            append!(component_artifacts.artifacts[:manual], manual_removals)
                            unique!(component_artifacts.artifacts[:manual])
                        end
                    end

                    @info subsection("Removing ICA components")
                    all_removed_components = get_all_ica_components(component_artifacts)
                    @info "Removed $(length(all_removed_components)) ICA components" component_artifacts
                    subtract_ica_components!(dat, ica, component_selection = components(all_removed_components))

                    # save ica results
                    if cfg["files"]["output"]["save_ica_data"]
                        @info "Saving ica data"
                        jldsave(_make_output_filename(output_directory, data_file, "_ica"); data = ica)
                    end

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
                    channel_out = Symbol("is_artifact_value" * "_" * string(preprocess_cfg.eeg.artifact_value_abs_criterion)),
                )

                # Save the cleaned continuous data in Julia format
                if cfg["files"]["output"]["save_continuous_data_corrected"]
                    @info "Saving continuous data"
                    jldsave(_make_output_filename(output_directory, data_file, "_continuous_corrected"); data = dat)
                end

                #################### EPOCH EXTRACTION ###################
                @info section("Extracting cleaned epoched data")
                epochs = extract_epochs(dat, epoch_cfgs, (preprocess_cfg.epoch_start, preprocess_cfg.epoch_end))

                # Check if any epochs have empty data
                empty_epochs = [i for (i, ep) in enumerate(epochs) if isempty(ep.data)]
                if !isempty(empty_epochs)
                    @minimal_error "Epoch extraction resulted in empty epochs for conditions: $(join([epochs[i].condition_name for i in empty_epochs], ", ")). Check epoch interval parameters and trigger locations."
                end

                #################### BASELINE WHOLE EPOCHS ##############
                @info section("Baseline whole epochs")
                baseline!(epochs)

                #################### DETECT BAD EPOCHS ###################
                @info section("Automatic epoch detection")
                
                # Determine interval selection for artifact rejection
                if !isnothing(preprocess_cfg.eeg.artifact_interval_start) && !isnothing(preprocess_cfg.eeg.artifact_interval_end)
                    rej_interval = times(preprocess_cfg.eeg.artifact_interval_start, preprocess_cfg.eeg.artifact_interval_end)
                    @info "Restricting artifact rejection to window: $(preprocess_cfg.eeg.artifact_interval_start)s to $(preprocess_cfg.eeg.artifact_interval_end)s"
                else
                    rej_interval = times()
                end

                rejection_info_step1 = detect_bad_epochs_automatic(
                    epochs;
                    z_criterion = preprocess_cfg.eeg.artifact_value_z_criterion,
                    abs_criterion = preprocess_cfg.eeg.artifact_value_abs_criterion,
                    interval_selection = rej_interval,
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
                if cfg["files"]["output"]["save_epoch_data_unrejected"]
                    @info "Saving epoch data (cleaned)"
                    jldsave(_make_output_filename(output_directory, data_file, "_epochs_unrejected"); data = epochs)
                end

                #################### RE-DETECT ARTIFACTS AFTER REPAIR ###################
                # Re-detect artifacts after repair to get updated rejection info
                @info subsection("Re-detecting artifacts after repair")
                rejection_info_step2 = detect_bad_epochs_automatic(
                    epochs;
                    z_criterion = preprocess_cfg.eeg.artifact_value_z_criterion,
                    abs_criterion = preprocess_cfg.eeg.artifact_value_abs_criterion,
                    interval_selection = rej_interval,
                    name = "rejection_step2",
                )
                channel_repairable!(rejection_info_step2, epochs[1].layout)
                @info "" # formatting
                @info rejection_info_step2

                #################### COMPARE REJECTION STEPS ###################
                @info subsection("Rejection Step Comparison (before vs after repair)")
                rejection_comparison = compare_rejections(rejection_info_step1, rejection_info_step2)
                log_pretty_table(rejection_comparison; title = "Rejection Step Comparison: Effectiveness of Channel Repair")

                if preprocess_cfg.interactive_epochs
                    @info section("Interactive Epoch Rejection Review")
                    @info "Review epoch artifacts. Checked epochs will be rejected. Close the window to continue."
                    state = detect_bad_epochs_interactive(epochs; artifact_info=rejection_info_step2)
                    wait(state.fig.scene)
                    
                    # Merge manual rejection state into the automatic rejection info
                    interactive_info = EegFun._to_rejection_info(state)
                    # We can just use the interactive_info directly as it contains all manual rejections!
                    # Wait, the interactive GUI starts with checkboxes active for automatically rejected epochs!
                    # So state.rejected contains BOTH automatic AND manual rejections.
                    # We can just override rejection_info_step2 with this combined info
                    rejection_info_step2 = interactive_info
                end

                #################### SAVE ERP DATA ###################
                if cfg["files"]["output"]["save_erp_data_unrejected"]
                    erps = average_epochs(epochs)
                    @info "Saving ERP data (cleaned)"
                    jldsave(_make_output_filename(output_directory, data_file, "_erps_unrejected"); data = erps)
                end

                #################### EPOCH REJECTION ###################
                @info subsection("Rejecting bad epochs")
                epochs = reject_epochs(epochs, rejection_info_step2)

                #################### SAVE ARTIFACT INFO ###################
                # Collect all artifact-related info into a single structure
                @info subsection("Artifact Information")
                artifact_info = ArtifactInfo(
                    !isnothing(continuous_repair_info) ? [continuous_repair_info] : ContinuousRepairInfo[],
                    vcat(rejection_info_step1, rejection_info_step2),
                    component_artifacts,  # Save ICA components if ICA was applied, otherwise nothing
                )
                jldsave(_make_output_filename(output_directory, data_file, "_artifact_info"); data = artifact_info)
                @info "Saved artifact info: $(artifact_info)"

                #################### LOG EPOCH COUNTS AND STORE FOR SUMMARY ###################
                df = log_epochs_table(epochs_uncorrected, epochs, title = "Epoch counts per condition (after repair and rejection):")
                push!(all_epoch_counts, df)

                #################### SAVE EPOCH DATA ###################
                if cfg["files"]["output"]["save_epoch_data_final"]
                    @info "Saving epoch data (good)"
                    jldsave(_make_output_filename(output_directory, data_file, "_epochs_final"); data = epochs)
                end

                #################### SAVE ERP DATA ###################
                if cfg["files"]["output"]["save_erp_data_final"]
                    erps = average_epochs(epochs)
                    @info "Saving ERP data (good)"
                    jldsave(_make_output_filename(output_directory, data_file, "_erps_final"); data = erps)
                end

                @info section("End of Processing")
                file_elapsed = time() - file_start_time
                total_processing_time += file_elapsed
                n_processed_for_eta += 1
                
                @info "Successfully processed file $file_idx/$(length(raw_data_files)): $(basename(data_file)) in $(_format_time(file_elapsed))"
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
        @info "$processed_files success, $(length(failed_files)) fail$(skipped_files > 0 ? ", $skipped_files skipped" : "")"
        !isempty(failed_files) && @info "Failed files: $(join(failed_files, ", "))"

        # Print electrode repair summary (load from saved artifact info files)
        @info subsection("Electrode Repair Summary Across All Participants (Continuous Level Only)")
        electrode_per_file, electrode_summary = summarize_electrode_repairs("_artifact_info", input_dir = output_directory)
        if !isempty(electrode_per_file)
            log_pretty_table(electrode_per_file; title = "Electrode Repairs per Participant:")
        end
        if !isempty(electrode_summary)
            log_pretty_table(electrode_summary; title = "Electrode Repairs: Number of Participants Affected")
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
            log_pretty_table(merged_epoch_summary, title = "Combined epoch counts across all files:", alignment = [:l, :r, :l, :r, :r, :r])
            log_pretty_table(
                merged_file_summary,
                title = "Average percentage per condition (averaged across conditions):",
                alignment = [:l, :r],
            )
            @info "Mean percentage (averaged across all conditions and files): $(round(mean(merged_file_summary.percentage), digits = 1)) %"
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
        elseif isfile(log_filename)
            rm(log_filename, force = true)
        end
    end

end
