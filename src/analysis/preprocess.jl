
"""
    preprocess_eeg_data(config::String; log::Bool = true, global_log_file::String = "")

Preprocess EEG data according to the specified configuration file.

# Arguments
- `config::String`: Path to the configuration file in TOML format

# Suggested improvements:
# TODO: Break into smaller functions:
#   - setup_preprocessing(config) -> (config_data, layout, epoch_conditions, output_dir)
#   - process_single_file(file, config_data, layout, epoch_conditions, output_dir)
#   - generate_preprocessing_report(results, output_dir)
"""
function preprocess_eeg_data(config::String)

    # set up the global log for overall processing
    global_log = setup_global_logging("preprocess_eeg_data.log")

    # initialize variable for outer scope
    output_directory = ""
    all_epoch_counts = DataFrame[]  # Vector to store all epoch counts

    try
        @info "EEG Preprocessing started at $(now()) ..."
        @info "Configuration file: $config"

        # check if config file exists and load
        if !isfile(config)
            error_msg = "Config file does not exist: $config"
            @minimal_error error_msg
        end
        cfg = load_config(config)

        # try and merge user config above with default config
        default_config = load_config(joinpath(@__DIR__, "..", "..", "src", "config", "default.toml"))
        cfg = _merge_configs(default_config, cfg)
        @info "Configuration loaded and merged successfully ..."

        # check if all requested raw data files exist
        @info "Checking if raw data file(s) exist ..."
        raw_data_files = get_files(cfg["files"]["input"]["directory"], cfg["files"]["input"]["raw_data_files"])
        raw_data_files_exist = check_files_exist(raw_data_files)
        if !raw_data_files_exist
            @minimal_error "Missing raw data files requested within TOML file!"
        end
        @info "Found $(length(raw_data_files)) raw data files: $(join(raw_data_files, ", "))"

        # We now need some code to read the epoch conditions defined within the toml file
        # For ease, these are defined in a separate toml file
        # See XXX for examples
        @info "Reading epoch conditions from $(cfg["files"]["input"]["epoch_condition_file"])"
        if !isfile(cfg["files"]["input"]["epoch_condition_file"])
            @minimal_error "Epoch definition file does not exist: $(cfg["files"]["input"]["epoch_condition_file"])"
        end
        epoch_cfgs = parse_epoch_conditions(TOML.parsefile(cfg["files"]["input"]["epoch_condition_file"]))
        @info "Epoch conditions loaded successfully"

        # Find and load layout file
        layout_file = find_file(cfg["files"]["input"]["layout_file"], joinpath(@__DIR__, "..", "..", "data", "layouts"))
        if layout_file === nothing
            @minimal_error "Electrode configuration file not found: $layout_name"
        end
        layout = read_layout(layout_file)
        @info "Electrode layout loaded from $layout_file"

        # Check if requested output directory exists and if not, create it
        output_directory = cfg["files"]["output"]["directory"]
        if !isdir(output_directory)
            mkdir(output_directory)
            @info "Created output directory: $output_directory"
        end

        # print config to output directory
        print_config(cfg, joinpath(output_directory, "config.toml"))

        # take layout file and add in 2D and 3D coordinates
        polar_to_cartesian_xy!(layout)
        polar_to_cartesian_xyz!(layout)

        get_layout_neighbours_xy!(layout, cfg["preprocess"]["layout"]["neighbour_criterion"])
        print_layout_neighbours(layout, joinpath(output_directory, "neighbours_xy.toml"))

        # Track processing results
        processed_files = 0
        failed_files = String[]

        # Process files in parallel
        for data_file in raw_data_files
            try
                @info "Processing $data_file"

                # Set up per-file logging (temporarily replaces global logger)
                setup_logging(joinpath(output_directory, "$(basename_without_ext(data_file)).log"))

                # read raw data file and create our Julia DataFrame
                dat = create_eeg_dataframe(read_bdf(data_file), layout)

                # rereference the data
                rereference!(dat, Symbol(cfg["preprocess"]["reference_channel"]))

                # Save the results
                if cfg["files"]["output"]["save_continuous_data"]
                    @info "Saving continuous data"
                    jldsave(make_output_filename(output_directory, data_file, "_continuous"); dat = dat)
                end

                # initial high-pass filter to remove DC offset/slow drifts
                if cfg["filter"]["highpass"]["apply"]
                    filter_data!(
                        dat,
                        "hp",
                        cfg["filter"]["highpass"]["freq"];
                        order = cfg["filter"]["highpass"]["order"],
                        filter_method = cfg["filter"]["highpass"]["method"],
                        filter_func = cfg["filter"]["highpass"]["func"],
                    )
                end

                # initial low-pass filter 
                if cfg["filter"]["lowpass"]["apply"]
                    filter_data!(
                        dat,
                        "lp",
                        cfg["filter"]["lowpass"]["freq"];
                        order = cfg["filter"]["lowpass"]["order"],
                        filter_method = cfg["filter"]["lowpass"]["method"],
                        filter_func = cfg["filter"]["lowpass"]["func"],
                    )
                end

                # caculate vEOG and hEOG channels
                channel_difference!(
                    dat,
                    channel_selection1 = channels(Symbol.(cfg["preprocess"]["eog"]["vEOG_channels"][1])),
                    channel_selection2 = channels(Symbol.(cfg["preprocess"]["eog"]["vEOG_channels"][2])),
                    channel_out        = Symbol.(cfg["preprocess"]["eog"]["vEOG_channels"][3][1]),
                )
                channel_difference!(
                    dat,
                    channel_selection1 = channels(Symbol.(cfg["preprocess"]["eog"]["hEOG_channels"][1])),
                    channel_selection2 = channels(Symbol.(cfg["preprocess"]["eog"]["hEOG_channels"][2])),
                    channel_out        = Symbol.(cfg["preprocess"]["eog"]["hEOG_channels"][3][1]),
                )

                # autodetect EOG signals
                detect_eog_onsets!(
                    dat,
                    cfg["preprocess"]["eog"]["vEOG_criterion"],
                    Symbol(cfg["preprocess"]["eog"]["vEOG_channels"][3][1]),
                    Symbol("is_" * cfg["preprocess"]["eog"]["vEOG_channels"][3][1]),
                )
                detect_eog_onsets!(
                    dat,
                    cfg["preprocess"]["eog"]["hEOG_criterion"],
                    Symbol(cfg["preprocess"]["eog"]["hEOG_channels"][3][1]),
                    Symbol("is_" * cfg["preprocess"]["eog"]["hEOG_channels"][3][1]),
                )

                # detect extreme values
                is_extreme_value!(
                    dat,
                    cfg["preprocess"]["eeg"]["extreme_value_criterion"],
                    mode = :combined,
                    channel_out = Symbol(
                        "is_extreme_value" * "_" * string(cfg["preprocess"]["eeg"]["extreme_value_criterion"]),
                    ),
                )

                # We perform the ica on "continuous" data (clean sections) that usually has a 
                # more extreme high-pass filter run ica on clean sections of "continuous" data
                if cfg["ica"]["apply"]

                    dat_ica = copy(dat)
                    dat_cleaned = copy(dat)

                    if cfg["filter"]["ica_highpass"]["apply"]
                        # apply high-pass filter to data
                        filter_data!(
                            dat_ica,
                            "hp",
                            cfg["filter"]["ica_highpass"]["freq"];
                            order = cfg["filter"]["ica_highpass"]["order"],
                            filter_method = cfg["filter"]["ica_highpass"]["method"],
                            filter_func = cfg["filter"]["ica_highpass"]["func"],
                        )
                    end

                    if cfg["filter"]["ica_lowpass"]["apply"]
                        # apply low-pass filter to data
                        filter_data!(
                            dat_ica,
                            "lp",
                            cfg["filter"]["ica_lowpass"]["freq"];
                            order = cfg["filter"]["ica_lowpass"]["order"],
                            filter_method = cfg["filter"]["ica_lowpass"]["method"],
                            filter_func = cfg["filter"]["ica_lowpass"]["func"],
                        )
                    end

                    ica_result = run_ica(
                        dat_ica;
                        sample_selection = samples_not(
                            Symbol(
                                "is_extreme_value" * "_" * string(cfg["preprocess"]["eeg"]["extreme_value_criterion"]),
                            ),
                        ),
                    )

                    # automatically identify components that are likely to be artifacts
                    eog_comps, eog_comps_metrics_df = identify_eog_components(
                        ica_result,
                        dat_ica,
                        sample_selection = samples_not(
                            Symbol(
                                "is_extreme_value" * "_" * string(cfg["preprocess"]["eeg"]["extreme_value_criterion"]),
                            ),
                        ),
                    )
                    ecg_comps, ecg_comps_metrics_df = identify_ecg_components(
                        ica_result,
                        dat_ica,
                        sample_selection = samples_not(
                            Symbol(
                                "is_extreme_value" * "_" * string(cfg["preprocess"]["eeg"]["extreme_value_criterion"]),
                            ),
                        ),
                    )
                    line_noise_comps, line_noise_comps_metrics_df = identify_line_noise_components(
                        ica_result,
                        dat_ica,
                        sample_selection = samples_not(
                            Symbol(
                                "is_extreme_value" * "_" * string(cfg["preprocess"]["eeg"]["extreme_value_criterion"]),
                            ),
                        ),
                    )
                    channel_noise_comps, channel_noise_comps_metrics_df =
                        identify_spatial_kurtosis_components(ica_result, dat_ica)

                    # Combine above component artifact results into a single structure
                    component_artifacts =
                        combine_artifact_components(eog_comps, ecg_comps, line_noise_comps, channel_noise_comps)
                    println(component_artifacts)
                    # plot_topography(ica_result)
                    plot_ica_component_activation(dat, ica_result)

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
                epochs_original = [
                    extract_epochs(
                        dat,
                        idx,
                        epoch_cfg,
                        cfg["preprocess"]["epoch_start"],
                        cfg["preprocess"]["epoch_end"],
                    ) for (idx, epoch_cfg) in enumerate(epoch_cfgs)
                ]

                # detect standard artifact values
                is_extreme_value!(
                    dat_cleaned,
                    cfg["preprocess"]["eeg"]["artifact_value_criterion"],
                    mode = :combined,
                    channel_out = Symbol(
                        "is_artifact_value" * "_" * string(cfg["preprocess"]["eeg"]["artifact_value_criterion"]),
                    ),
                )

                epochs_cleaned = [
                    extract_epochs(
                        dat_cleaned,
                        idx,
                        epoch_cfg,
                        cfg["preprocess"]["epoch_start"],
                        cfg["preprocess"]["epoch_end"],
                    ) for (idx, epoch_cfg) in enumerate(epoch_cfgs)
                ]

                epochs_cleaned = reject_epochs.(
                    epochs_cleaned,
                    Ref(
                        Symbol(
                            "is_artifact_value" * "_" * string(cfg["preprocess"]["eeg"]["artifact_value_criterion"]),
                        ),
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
                    erps_original = [average_epochs(epoch) for epoch in epochs_original]
                    @info "Saving erp data (original)"
                    jldsave(make_output_filename(output_directory, data_file, "_erps_original"); erps = erps_original)
                end

                if cfg["files"]["output"]["save_erp_data_cleaned"]
                    erps_cleaned = [average_epochs(epoch) for epoch in epochs_cleaned]
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
        @info "Preprocessing complete: $processed_files files processed successfully, $(length(failed_files)) files failed"
        if !isempty(failed_files)
            @info "Failed files: $(join(failed_files, ", "))"
        end

        # Print combined epoch counts
        if !isempty(all_epoch_counts)
            combined_counts = vcat(all_epoch_counts...)
            log_pretty_table(combined_counts, title = "Combined epoch counts across all files:")
            jldsave(joinpath(output_directory, "epoch_summary.jld2"); df = combined_counts)
        end

    finally
        # close global logging and move log file to output directory
        close_global_logging()
        log_source = "preprocess_eeg_data.log"
        log_dest = joinpath(output_directory, "preprocess_eeg_data.log")
        if log_source != log_dest
            mv(log_source, log_dest, force = true)
        end
    end

end
