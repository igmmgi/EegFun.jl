"""
    preprocess(config::String; log_level::String = "info")

Preprocess EEG data according to the specified configuration file.

# Arguments
- `config::String`: Path to the configuration file in TOML format
- `log_level::String`: Log level for preprocessing ("debug", "info", "warn", "error")
"""
function preprocess(config::String; log_level::String = "info")

    # set up the global log for overall processing
    global_log = setup_global_logging("preprocess_eeg_data.log", log_level = log_level)

    # initialize variable for outer scope
    output_directory = ""
    all_epoch_counts = DataFrame[]  # Vector to store all epoch counts

    try

        # Setup for all analyses files:
        # This involves loading the end-user config file, merging it with the default config, 
        # and creating the PreprocessConfig object.
        # The epoch conditions files are also loaded and parsed (see XXX for example configuration files)

        @info section("Setup")
        !isfile(config) && @minimal_error "Config file does not exist: $config"
        cfg = load_config(config)
        cfg == nothing && @minimal_error "Failed to load configuration from: $config"

        # try and merge user config above with default config
        default_config = load_config(joinpath(@__DIR__, "..", "..", "src", "config", "default.toml"))
        default_config == nothing && @minimal_error "Failed to load default configuration"
        cfg = _merge_configs(default_config, cfg)

        # Create the PreprocessConfig object for log level
        preprocess_cfg = PreprocessConfig(cfg["preprocess"])

        # check if all requested raw data files exist
        raw_data_files = get_files(cfg["files"]["input"]["directory"], cfg["files"]["input"]["raw_data_files"])
        raw_data_files_exist = check_files_exist(raw_data_files)
        !raw_data_files_exist && @minimal_error "Missing raw data files requested within TOML file!"
        @info "Found $(length(raw_data_files)) files: $(join(raw_data_files, ", "))"

        # Read the epoch conditions defined within the toml file (See XXX for examples)
        !isfile(cfg["files"]["input"]["epoch_condition_file"]) &&
            @minimal_error "File missing: $(cfg["files"]["input"]["epoch_condition_file"])"
        epoch_cfgs = parse_epoch_conditions(TOML.parsefile(cfg["files"]["input"]["epoch_condition_file"]))
        @info "Epoch file: $(cfg["files"]["input"]["epoch_condition_file"]) loaded"

        # Find and load layout file
        layout_file = find_file(cfg["files"]["input"]["layout_file"], joinpath(@__DIR__, "..", "..", "data", "layouts"))
        layout_file === nothing && @minimal_error "Layout file not found: $layout_name"
        layout = read_layout(layout_file)

        # Check if requested output directory exists and if not, create it
        output_directory = cfg["files"]["output"]["directory"]
        !isdir(output_directory) && mkdir(output_directory)
        @info "Output directory: $output_directory"

        # print config to output directory
        print_config(cfg, joinpath(output_directory, "config.toml"))

        # Layout coordinates and calculation of channel neighbours (2D)
        polar_to_cartesian_xy!(layout)
        polar_to_cartesian_xyz!(layout)
        get_layout_neighbours_xy!(layout, cfg["preprocess"]["layout"]["neighbour_criterion"])
        print_layout_neighbours(layout, joinpath(output_directory, "neighbours_xy.toml"))

        # Actual start of preprocessing pipeline!
        # This is the main loop that processes each raw data file.
        # TODO: embarrasingly parallel? use Threads.@threads?
        # Track processing results
        processed_files = 0
        failed_files = String[]
        for data_file in raw_data_files

            try
                @info section("Processing")
                @info "File: $data_file"

                # Set up per-file logging (temporarily replaces global logger)
                setup_logging(
                    joinpath(output_directory, "$(basename_without_ext(data_file)).log"),
                    log_level = log_level,
                )

                ################### LOAD RAW DATA FILE ###################
                # TODO: update for different file types!
                @info section("Raw Data")
                dat = create_eeg_dataframe(read_bdf(data_file), layout)

                # Save the original data in Julia format
                if cfg["files"]["output"]["save_continuous_data"]
                    @info "Saving continuous data"
                    jldsave(make_output_filename(output_directory, data_file, "_continuous"); dat = dat)
                end

                ################### REREFERENCE DATA ###################
                @info section("Rereference")
                rereference!(dat, preprocess_cfg.reference_channel)

                ################### APPLY FILTERS ###################
                @info section("Initial Filters")
                @info "Continuous data filters: $(_applied_filters(preprocess_cfg.filter, filter_sections = [:highpass, :lowpass]))"
                filter_data!(dat, preprocess_cfg.filter)

                ################### INITIAL CHANNEL SUMMARY ###################
                @debug section("Channel Summary")
                summary = channel_summary(dat)
                log_pretty_table(summary; title = "Channel Summary (whole dataset)", log_level = :debug)

                # Mark epoch windows
                mark_epoch_windows!(dat, epoch_cfgs, [preprocess_cfg.epoch_start, preprocess_cfg.epoch_end])
                summary = channel_summary(dat, sample_selection = samples(:epoch_window))
                log_pretty_table(summary; title = "Channel Summary (epoch window)", log_level = :debug)

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
                log_pretty_table(hEOG_vEOG_cm; title = "Channel x vEOG/hEOG Correlation Matrix (whole dataset)", log_level = :debug)

                # Calculate correlations between all channels and EOG channels
                @info subsection("Channel x vEOG/hEOG Correlation Matrix")
                hEOG_vEOG_cm = eegfun.correlation_matrix_eog(dat, sample_selection = samples(:epoch_window), preprocess_cfg.eog)
                eegfun.add_zscore_columns!(hEOG_vEOG_cm)
                log_pretty_table(hEOG_vEOG_cm; title = "Channel x vEOG/hEOG Correlation Matrix (epoch window)", log_level = :debug)

                #################### DETECT BAD ELECTRODES IN CONTINUOUS DATA ###################
                @info section("Bad Channel Detection")
                cjp = channel_joint_probability(dat)
                log_pretty_table(cjp; title = "Channel Joint Probability (whole dataset)", log_level = :debug)

                cjp = channel_joint_probability(dat, sample_selection = samples(:epoch_window))
                log_pretty_table(cjp; title = "Channel Joint Probability (epoch window)", log_level = :debug)


                #################### DETECT EXTREME VALUES IN CONTINUOUS DATA ###################
                @info section("Artifact Detection (extreme values)")
                @info "Detecting extreme values: $(preprocess_cfg.eeg.artifact_value_criterion) criterion"
                is_extreme_value!(
                    dat,
                    Int(preprocess_cfg.eeg.artifact_value_criterion),
                    channel_out = Symbol(
                        "is_extreme_value" * "_" * string(preprocess_cfg.eeg.artifact_value_criterion),
                    ),
                )

                #################### Independent Component Analysis (ICA) ###################
                # We perform the ica on "continuous" data (clean sections) that usually has a 
                # more extreme high-pass filter applied. 
                # We then run ica on clean sections of "continuous" data
                if preprocess_cfg.ica.apply

                    @info section("ICA")
                    dat_ica = copy(dat)
                    dat_cleaned = copy(dat)

                    # Apply ICA-specific filters
                    @info subsection("ICA filters")
                    @info "ICA data filters: $(_applied_filters(preprocess_cfg.filter, filter_sections = [:ica_highpass, :ica_lowpass]))"
                    filter_data!(dat_ica, preprocess_cfg.filter, filter_sections = [:ica_highpass, :ica_lowpass])

                    @info subsection("Running ICA")
                    ica_result = run_ica(
                        dat_ica;
                        sample_selection = samples_not(
                            Symbol("is_extreme_value" * "_" * string(preprocess_cfg.eeg.artifact_value_criterion)),
                        ),
                        percentage_of_data = preprocess_cfg.ica.percentage_of_data,
                    )

                    # Identify all artifact components 
                    @info subsection("Component Identification")
                    component_artifacts, component_metrics = identify_components(
                        dat_ica,
                        ica_result,
                        sample_selection = samples_not(
                            Symbol("is_extreme_value" * "_" * string(preprocess_cfg.eeg.artifact_value_criterion)),
                        ),
                    )

                    # Print component metrics to log files
                    log_pretty_table(component_metrics[:eog_metrics]; title = "EOG Component Metrics", log_level = :debug)
                    log_pretty_table(component_metrics[:ecg_metrics], title = "ECG Component Metrics", log_level = :debug)
                    log_pretty_table(component_metrics[:line_noise_metrics], title = "Line Noise Component Metrics", log_level = :debug)
                    log_pretty_table(component_metrics[:channel_noise_metrics], title = "Channel Noise Component Metrics", log_level = :debug)

                    @info subsection("Removing ICA components")
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

                #################### EPOCHING DATA ###################
                @info section("Epoching")
                epochs_original = extract_epochs(dat, epoch_cfgs, preprocess_cfg.epoch_start, preprocess_cfg.epoch_end)

                #################### DETECT STANDARD ARTIFACT VALUES IN EPOCHED DATA ###################
                @info subsection("Detecting artifact values in epoched data")
                is_extreme_value!(
                    dat_cleaned,
                    Int(preprocess_cfg.eeg.artifact_value_criterion),
                    channel_out = Symbol(
                        "is_artifact_value" * "_" * string(preprocess_cfg.eeg.artifact_value_criterion),
                    ),
                )

                @info subsection("Extracting cleaned epoched data")
                epochs_cleaned =
                    extract_epochs(dat_cleaned, epoch_cfgs, preprocess_cfg.epoch_start, preprocess_cfg.epoch_end)
                epochs_cleaned = reject_epochs(
                    epochs_cleaned,
                    Symbol("is_artifact_value" * "_" * string(preprocess_cfg.eeg.artifact_value_criterion)),
                );

                #################### LOG EPOCH COUNTS AND STORE FOR SUMMARY ###################
                df = log_epochs_table(epochs_original, epochs_cleaned, title = "Epoch counts per condition:")
                push!(all_epoch_counts, df)

                # save epochs
                @info section("Saving data")
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

                #################### SAVE ERP DATA ###################
                if cfg["files"]["output"]["save_erp_data_original"]
                    erps_original = average_epochs(epochs_original)
                    @info "Saving ERP data (original)"
                    jldsave(make_output_filename(output_directory, data_file, "_erps_original"); erps = erps_original)
                end

                if cfg["files"]["output"]["save_erp_data_cleaned"]
                    erps_cleaned = average_epochs(epochs_cleaned)
                    @info "Saving ERP data (cleaned)"
                    jldsave(make_output_filename(output_directory, data_file, "_erps_cleaned"); erps = erps_cleaned)
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

# Helper function to return applied filters string during preprocessing
function _applied_filters(filter_cfg::FilterConfig; filter_sections::Vector{Symbol} = [:highpass, :lowpass])
    sections = [getfield(filter_cfg, section) for section in filter_sections]
    applied_filters = []
    for section in sections
        if section.apply
            push!(
                applied_filters,
                "$(section.freq)Hz $(section.type), $(section.method), $(section.func), order $(section.order)",
            )
        end
    end
    return join(applied_filters, "; ")
end

# Helper function to return EOG configuration string for printing
function _eog_config_string(eog_cfg::EogConfig)
    vEOG_channels =
        join(eog_cfg.vEOG_channels[1], ", ") *
        " - " *
        join(eog_cfg.vEOG_channels[2], ", ") *
        " → " *
        eog_cfg.vEOG_channels[3][1]
    hEOG_channels =
        join(eog_cfg.hEOG_channels[1], ", ") *
        " - " *
        join(eog_cfg.hEOG_channels[2], ", ") *
        " → " *
        eog_cfg.hEOG_channels[3][1]
    return "vEOG: $vEOG_channels ($(eog_cfg.vEOG_criterion)μV); hEOG: $hEOG_channels ($(eog_cfg.hEOG_criterion)μV)"
end

# Helper functions for section headers
function _center_title(title::String, width::Int)
    title_length = length(title)
    total_dashes = width - title_length - 2  # 2 spaces around title
    left_dashes = div(total_dashes, 2)
    right_dashes = total_dashes - left_dashes
    return "-" ^ left_dashes * " $title " * "-" ^ right_dashes
end

function section(title::String; width::Int = 80)
    dash_line = "-" ^ width
    middle_line = _center_title(title, width)
    return "\n$dash_line\n$middle_line\n$dash_line"
end

function subsection(title::String; width::Int = 80)
    return "\n" * _center_title(title, width)
end
