
"""
    preprocess_eeg_data(config::String; log::Bool = true, global_log_file::String = "")

Preprocess EEG data according to the specified configuration file.

# Arguments
- `config::String`: Path to the configuration file in TOML format
"""
function preprocess_eeg_data(config::String)

    # set up the global log for overall processing
    global_log = setup_global_logging("preprocess_eeg_data.log")

    # initialize variable for outer scope
    output_directory = ""
    all_epoch_counts = DataFrame[]  # Vector to store all epoch counts

    try

        # logging start of preprocessing, load user config, default config and merge
        @info "EEG Preprocessing started at $(now()) ..."
        @info "Configuration file: $config"
        !isfile(config) && @minimal_error "Config file does not exist: $config"
        cfg = load_config(config)

        # try and merge user config above with default config
        default_config = load_config(joinpath(@__DIR__, "..", "..", "src", "config", "default.toml"))
        cfg = _merge_configs(default_config, cfg)
        @info "Configuration loaded and merged successfully ..."

        # check if all requested raw data files exist
        @info "Checking if raw data file(s) exist ..."
        raw_data_files = get_files(cfg["files"]["input"]["directory"], cfg["files"]["input"]["raw_data_files"])
        raw_data_files_exist = check_files_exist(raw_data_files)
        !raw_data_files_exist && @minimal_error "Missing raw data files requested within TOML file!"
        @info "Found $(length(raw_data_files)) raw data files: $(join(raw_data_files, ", "))"

        # Read the epoch conditions defined within the toml file (See XXX for examples)
        !isfile(cfg["files"]["input"]["epoch_condition_file"]) && @minimal_error "File missing: $(cfg["files"]["input"]["epoch_condition_file"])"
        epoch_cfgs = parse_epoch_conditions(TOML.parsefile(cfg["files"]["input"]["epoch_condition_file"]))
        @info "Epoch conditions file: $(cfg["files"]["input"]["epoch_condition_file"]) loaded successfully"

        # Find and load layout file
        layout_file = find_file(cfg["files"]["input"]["layout_file"], joinpath(@__DIR__, "..", "..", "data", "layouts"))
        layout_file === nothing && @minimal_error "Electrode configuration file not found: $layout_name"
        layout = read_layout(layout_file)
        @info "Layout file: $layout_file loaded successfully"

        # Check if requested output directory exists and if not, create it
        output_directory = cfg["files"]["output"]["directory"]
        !isdir(output_directory) && mkdir(output_directory)
        @info "Output directory: $output_directory"

        # print config to output directory
        print_config(cfg, joinpath(output_directory, "config.toml"))
        @info "Configuration printed to output directory: $output_directory"

        # Layout coordinates and calculation of channel neighbours (2D)
        polar_to_cartesian_xy!(layout)
        polar_to_cartesian_xyz!(layout)
        get_layout_neighbours_xy!(layout, cfg["preprocess"]["layout"]["neighbour_criterion"])
        print_layout_neighbours(layout, joinpath(output_directory, "neighbours_xy.toml"))
        @info "Layout neighbours printed to output directory: $output_directory"

        # Actual start of preprocessing pipeline!
        # Track processing results
        processed_files = 0
        failed_files = String[]
        # TODO: embarrasingly parallel? use Threads.@threads?
        for data_file in raw_data_files
            try
                @info "Processing $data_file"

                # Set up per-file logging (temporarily replaces global logger)
                setup_logging(joinpath(output_directory, "$(basename_without_ext(data_file)).log"))

                # read raw data file and create our Julia DataFrame
                # TODO: update for different file types!
                dat = create_eeg_dataframe(read_bdf(data_file), layout)

                # Save the original data in Julia format
                if cfg["files"]["output"]["save_continuous_data"]
                    @info "Saving continuous data"
                    jldsave(make_output_filename(output_directory, data_file, "_continuous"); dat = dat)
                end

                # Create comprehensive preprocessing configuration
                preprocess_cfg = PreprocessConfig(cfg["preprocess"])
                
                # rereference the data
                @info "Rereferencing with channel: $(preprocess_cfg.reference_channel)"
                rereference!(dat, Symbol(preprocess_cfg.reference_channel))

                # Apply filters based on configuration
                @info "Filtering data with configuration: $(preprocess_cfg.filter)"
                filter_data!(dat, preprocess_cfg.filter)

                # Calculate EOG channels based on configuration
                @info "Calculating EOG channels based on configuration: $(preprocess_cfg.eog)"
                calculate_eog_channels!(dat, preprocess_cfg.eog)

                # Autodetect EOG signals
                @info "Autodetecting EOG signals based on configuration: $(preprocess_cfg.eog)"
                detect_eog_signals!(dat, preprocess_cfg.eog)

                # detect extreme values
                @info "Detecting extreme values based on configuration: $(preprocess_cfg.eeg["artifact_value_criterion"])"
                is_extreme_value!(
                    dat,
                    preprocess_cfg.eeg["artifact_value_criterion"],
                    channel_out = Symbol( "is_extreme_value" * "_" * string(preprocess_cfg.eeg["artifact_value_criterion"])),
                )

                # We perform the ica on "continuous" data (clean sections) that usually has a 
                # more extreme high-pass filter run ica on clean sections of "continuous" data
                if cfg["ica"]["apply"]

                    dat_ica = copy(dat)
                    dat_cleaned = copy(dat)

                    # Apply ICA-specific filters
                    filter_data!(dat_ica, preprocess_cfg.filter, filter_sections = ["ica_highpass", "ica_lowpass"])

                    ica_result = run_ica(
                        dat_ica;
                        sample_selection = samples_not(
                            Symbol(
                                "is_extreme_value" * "_" * string(preprocess_cfg.eeg["artifact_value_criterion"]),
                            ),
                        ),
                    )

                    # Identify all artifact components in one unified call
                    component_artifacts, component_metrics = identify_components(
                        dat_ica,
                        ica_result,
                        sample_selection = samples_not(
                            Symbol(
                                "is_extreme_value" * "_" * string(preprocess_cfg.eeg["artifact_value_criterion"]),
                            ),
                        ),
                    )

                    # Print component metrics to log files
                    @info "EOG Component Metrics:"
                    log_pretty_table(component_metrics[:eog_metrics], title = "EOG Component Metrics")
                    @info "ECG Component Metrics:"
                    log_pretty_table(component_metrics[:ecg_metrics], title = "ECG Component Metrics")
                    @info "Line Noise Component Metrics:"
                    log_pretty_table(component_metrics[:line_noise_metrics], title = "Line Noise Component Metrics")
                    @info "Channel Noise Component Metrics:"
                    log_pretty_table(component_metrics[:channel_noise_metrics], title = "Channel Noise Component Metrics")

                    remove_ica_components!(
                        dat_cleaned,
                        ica_result,
                        component_selection = components(get_all_ica_components(component_artifacts)),
                    )

                    # save ica results
                    if cfg["files"]["output"]["save_ica_data"]
                        @info "Saving ica data"
                        jldsave(make_output_filename(output_directory, data_file, "_ica"); ica_result = ica_result)
                    end

                else
                    dat_cleaned = dat
                end

                # epoch data
                epochs_original = extract_epochs(dat, epoch_cfgs, preprocess_cfg.eeg["epoch_start"], preprocess_cfg.eeg["epoch_end"])

                # detect standard artifact values
                is_extreme_value!(
                    dat_cleaned,
                    preprocess_cfg.eeg["artifact_value_criterion"],
                    channel_out = Symbol(
                        "is_artifact_value" * "_" * string(preprocess_cfg.eeg["artifact_value_criterion"]),
                    ),
                )

                epochs_cleaned = extract_epochs(dat_cleaned, epoch_cfgs, preprocess_cfg.eeg["epoch_start"], preprocess_cfg.eeg["epoch_end"])


                epochs_cleaned = reject_epochs(
                    epochs_cleaned,
                    Symbol(
                        "is_artifact_value" * "_" * string(preprocess_cfg.eeg["artifact_value_criterion"]),
                    ),
                );

                # Log epoch counts and store for summary
                df = log_epochs_table(epochs_original, epochs_cleaned, title = "Epoch counts per condition:")
                push!(all_epoch_counts, df)

                # save epochs
                if cfg["files"]["output"]["save_epoch_data_original"]
                    @info "Saving epoch data (original)"
                    jldsave(
                        make_output_filename(output_directory, data_file, "_epochs_original");
                        epochs = epochs_original,
                    )
                end

                if cfg["files"]["output"]["save_epoch_data_cleaned"]
                    @info "Saving epoch data (cleaned)"
                    jldsave(
                        make_output_filename(output_directory, data_file, "_epochs_cleaned");
                        epochs = epochs_cleaned,
                    )
                end

                # save erp data
                if cfg["files"]["output"]["save_erp_data_original"]
                    erps_original = average_epochs(epochs_original)
                    @info "Saving erp data (original)"
                    jldsave(make_output_filename(output_directory, data_file, "_erps_original"); erps = erps_original)
                end

                if cfg["files"]["output"]["save_erp_data_cleaned"]
                    erps_cleaned = average_epochs(epochs_cleaned)
                    @info "Saving erp data (cleaned)"
                    jldsave(make_output_filename(output_directory, data_file, "_erps_cleaned"); erps = erps_cleaned)
                end

                @info "Successfully processed $data_file"
                processed_files += 1

            catch e # Log error to both console and global log
                @error "Error processing $data_file" exception=(e, catch_backtrace())
                push!(failed_files, data_file)
            finally # Close per-file logging (restores global logger)
                close_logging()
            end
        end

        # Write final summary
        @info "Preprocessing complete: $processed_files success, $(length(failed_files)) fail"
        !isempty(failed_files) && @info "Failed files: $(join(failed_files, ", "))"

        # Print combined epoch counts
        if !isempty(all_epoch_counts)
            combined_counts = vcat(all_epoch_counts...)
            log_pretty_table(combined_counts, title = "Combined epoch counts across all files:")
            jldsave(joinpath(output_directory, "epoch_summary.jld2"); df = combined_counts)
        end

    finally
        close_global_logging()
        log_source = "preprocess_eeg_data.log"
        log_destination = joinpath(output_directory, "preprocess_eeg_data.log")
        if log_source != log_destination
            mv(log_source, log_destination, force = true)
        end
    end

end
