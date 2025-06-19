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

        # Check if electrode configuration file exists and load
        layout_file = joinpath(@__DIR__, "..", "..", "data", "layouts", cfg["files"]["input"]["layout_file"])
        if !isfile(layout_file)
            # not default layout file so check if custom file specified
            layout_file = cfg["files"]["input"]["layout_file"]
            if !isfile(layout_file)
                @minimal_error "Electrode configuration file does not exist: $layout_file"
            end
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

        # take layout file and add in 2D and 3D coordinated
        polar_to_cartesian_xy!(layout)
        polar_to_cartesian_xyz!(layout)

        neighbours_xy = get_electrode_neighbours_xy(layout, cfg["preprocess"]["layout"]["neighbour_criterion"])
        print_neighbours_dict(neighbours_xy, joinpath(output_directory, "neighbours_xy.toml"))
        
        # Track processing results
        processed_files = 0
        failed_files = String[]

        # Process files in parallel
        @threads for file in raw_data_files
            
            try
                @info "Processing $file"
                
                # Set up per-file logging (temporarily replaces global logger)
                setup_logging(joinpath(output_directory, "$(basename_without_ext(file)).log"))

                # read raw data file and create our Julia DataFrame
                dat = create_eeg_dataframe(read_bdf(file), layout)

                # rereference the data
                rereference!(dat, Symbol(cfg["preprocess"]["reference_channel"]))

                # Save the results
                if cfg["files"]["output"]["save_continuous_data"]
                    @info "Saving continuous data"
                    jldsave(make_output_filename(output_directory, file, "_continuous"); dat=dat)
                end


                # initial high-pass filter to remove DC offset/slow drifts
                if cfg["filter"]["highpass"]["on"]
                    filter_data!(
                        dat,
                        "hp",
                        cfg["filter"]["highpass"]["type"],
                        cfg["filter"]["highpass"]["cutoff"],
                        order = cfg["filter"]["highpass"]["order"],
                    )
                end

                # initial low-pass filter 
                if cfg["filter"]["lowpass"]["on"]
                    filter_data!(
                        dat,
                        "lp",
                        cfg["filter"]["lowpass"]["type"],
                        cfg["filter"]["lowpass"]["cutoff"],
                        order = cfg["filter"]["lowpass"]["order"],
                    )
                end

                # TODO: do not hard code these channels!
                # caculate EOG channels
                diff_channel!(dat, [:Fp1, :Fp2], [:IO1, :IO2], :vEOG)
                diff_channel!(dat, :F9, :F10, :hEOG)

                # autodetect EOG signals
                # TODO: do not hard code output channel names!
                detect_eog_onsets!(dat, cfg["preprocess"]["eog"]["vEOG_criterion"], :vEOG, :is_vEOG)
                detect_eog_onsets!(dat, cfg["preprocess"]["eog"]["hEOG_criterion"], :hEOG, :is_hEOG)

                # detect extreme values
                # TODO: output channel names should be configurable
                is_extreme_value!(dat, dat.layout.label, cfg["preprocess"]["eog"]["extreme_value_criterion"])

                # We perform the ica on "continuous" data (clean sections) that usually has a 
                # more extreme high-pass filter run ica on clean sections of "continuous" data
                if cfg["ica"]["run"]

                    if cfg["ica"]["filter"]["highpass"]["on"]
                        # apply high-pass filter to data
                        dat_ica = filter_data(
                            dat,
                            "hp",
                            cfg["ica"]["filter"]["highpass"]["type"],
                            cfg["ica"]["filter"]["highpass"]["cutoff"],
                            order = cfg["ica"]["filter"]["highpass"]["order"],
                        )
                    else
                        dat_ica = copy(dat)
                    end

                    ica_result = run_ica(dat_ica; exclude_samples = [:is_extreme_value])

                    # save ica results
                    if cfg["files"]["output"]["save_ica_data"]
                        @info "Saving ica data"
                        jldsave(make_output_filename(output_directory, file, "_ica"); ica_result=ica_result)
                    end

                end

                # epoch data
                epochs = []
                for (idx, epoch_cfg) in enumerate(epoch_cfgs)
                    # Use the new extract_epochs function that accepts EpochCondition objects
                    push!(
                        epochs,
                        extract_epochs(
                            dat,
                            idx,
                            epoch_cfg,
                            cfg["preprocess"]["epoch_start"],
                            cfg["preprocess"]["epoch_end"],
                        ),
                    )
                end

                # Log epoch counts
                df = DataFrame(file = basename(file), condition = 1:length(epochs), n_epochs_total = n_epochs.(epochs))
                push!(all_epoch_counts, df)  # Store the DataFrame
                @info "Epoch counts per condition:\n$(pretty_table(String, df, show_row_number=false, show_subheader=false))"

                # save epochs
                if cfg["files"]["output"]["save_epoch_data"]
                    @info "Saving epoch data"
                    jldsave(make_output_filename(output_directory, file, "_epochs"); epochs=epochs)
                end

                # average epochs
                erps = []
                for (idx, epoch) in enumerate(epochs)
                    push!(erps, average_epochs(epochs[idx]))
                end

                df.n_epochs_erp = n_average.(erps)
                # calculate the percentage of epochs that went into the ERP
                df.n_epochs_erp_percentage = (df.n_epochs_erp ./ df.n_epochs_total) .* 100
                @info "Epoch counts per condition:\n$(pretty_table(String, df, show_row_number=false, show_subheader=false))"

                # save erps
                if cfg["files"]["output"]["save_erp_data"]
                    @info "Saving erp data"
                    jldsave(make_output_filename(output_directory, file, "_erps"); erps=erps)
                end

                @info "Successfully processed $file"
                processed_files += 1

            catch e # Log error to both console and global log
                @error "Error processing $file" exception=(e, catch_backtrace())
                push!(failed_files, file)
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
            @info "Combined epoch counts across all files:\n$(pretty_table(String, combined_counts, show_row_number=false, show_subheader=false))"
            jldsave(joinpath(output_directory, "epoch_summary.jld2"); df=combined_counts)
        end
        
    finally
        # close global logging and move log file to output directory
        close_global_logging()
        mv("preprocess_eeg_data.log", joinpath(output_directory, "preprocess_eeg_data.log"), force=true)
    end

end
