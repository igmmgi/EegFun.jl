using Base.Threads
using JLD2
using PrettyTables

"""
    make_output_filename(output_dir::String, input_file::String, suffix::String)

Create an output filename from input file path with given suffix.

# Arguments
- `output_dir::String`: Output directory path
- `input_file::String`: Input file path
- `suffix::String`: Suffix to add (e.g., "_ica", "_continuous")

# Returns
- `String`: Full output filename path

# Example
```julia
filename = make_output_filename("/output", "data/file.bdf", "_ica")
# Returns: "/output/file_ica.jld2"
```
"""
function make_output_filename(output_dir::String, input_file::String, suffix::String)
    base_name = splitext(basename(input_file))[1]
    return joinpath(output_dir, "$(base_name)$(suffix).jld2")
end

"""
    preprocess_eeg_data(config::String; log::Bool = true, global_log_file::String = "")

Preprocess EEG data according to the specified configuration file.

# Arguments
- `config::String`: Path to the configuration file in TOML format
"""
function preprocess_eeg_data(config::String)
    # Set up the global log
    global_log = setup_global_logging("preprocess_eeg_data.log")
    
    # Initialize variable for outer scope
    output_data_directory = ""
    all_epoch_counts = DataFrame[]  # Vector to store all epoch counts
    
    try
        @info "EEG Preprocessing started at $(now())"
        @info "Configuration file: $config"
        
        # check if config file exists and load
        if !isfile(config)
            error_msg = "Config file does not exist: $config"
            @minimal_error error_msg
        end
        config_data = load_config(config)

        # merge with the default config
        default_config = load_config(joinpath(@__DIR__, "..", "..", "src", "config", "default.toml"))
        config_data = _merge_configs(default_config, config_data)

        print_config(config_data)
        @info "Configuration loaded successfully"

        # check if all requested raw data files exist
        @info "Checking if raw data files exist"
        raw_data_files =
            get_files(config_data["files"]["input"]["raw_data_directory"], config_data["files"]["input"]["raw_data_files"])
        raw_data_files_exist = check_files_exist(raw_data_files)
        if !raw_data_files_exist
            error_msg = "Missing raw data files requested within TOML file"
            @minimal_error error_msg
        end
        @info "Found $(length(raw_data_files)) raw data files: $(join(raw_data_files, ", "))"

        # Check if electrode configuration file exists and load
        layout_file = joinpath(@__DIR__, "..", "..", "data", "layouts", config_data["files"]["input"]["layout_file"])
        if !isfile(layout_file)
            # not default layout file so check if custom file specified
            layout_file = config_data["files"]["input"]["layout_file"]
            if !isfile(layout_file)
                error_msg = "Electrode configuration file does not exist: $layout_file"
                @minimal_error error_msg
            end
        end
        layout = read_layout(layout_file)
        @info "Electrode layout loaded from $layout_file"

        # Check if requested output directory exists and if not, create it
        output_data_directory = config_data["files"]["output"]["output_data_directory"]
        if !isdir(output_data_directory)
            mkdir(output_data_directory)
            @info "Created output directory: $output_data_directory"
        end
        
        # delete current files in output directory
        if config_data["files"]["output"]["delete_current_files"]
            @info "Deleting current files in output directory: $output_data_directory"
            rm.(readdir(output_data_directory, join = true), recursive = true)
        end

        # print config to output directory
        print_config(config_data, joinpath(output_data_directory, "config.toml"))

        # take layout file and add in 2D and 3D coordinated
        polar_to_cartesian_xy!(layout)
        polar_to_cartesian_xyz!(layout)
        
        # Track processing results
        processed_files = 0
        failed_files = String[]

        # Process files in parallel
        @threads for file in raw_data_files
            
            try
                # Set up per-file logging (temporarily replaces global logger)
                setup_logging(joinpath(output_data_directory, splitext(basename(file))[1] * ".log"))
            
                @info "Processing $file"

                # read raw data file and create our Julia DataFrame
                dat = read_bdf(file)
                dat = create_eeg_dataframe(dat, layout)

                # rereference the data
                # Why does this not work?
                rereference!(dat, Symbol(config_data["preprocessing"]["reference_electrode"]))

                # Save the results
                if config_data["files"]["output"]["save_continuous_data"]
                    @info "Saving continuous data"
                    jldsave(make_output_filename(output_data_directory, file, "_continuous"); dat=dat)
                end

                # initial high-pass filter to remove slow drifts
                filter_data!(
                    dat,
                    "hp",
                    config_data["filtering"]["highpass"]["type"],
                    config_data["filtering"]["highpass"]["cutoff"],
                    order = config_data["filtering"]["highpass"]["order"],
                )

                # caculate EOG channels
                diff_channel!(dat, [:Fp1, :Fp2], [:IO1, :IO2], :vEOG)
                diff_channel!(dat, :F9, :F10, :hEOG)

                # autodetect EOG signals
                detect_eog_onsets!(dat, 50, :vEOG, :is_vEOG)
                detect_eog_onsets!(dat, 30, :hEOG, :is_hEOG)

                # detect extreme values
                is_extreme_value!(dat, dat.layout.label, 100)


                # run ica on clean sections of "continuous" data
                if config_data["ica"]["ica_run"]
                    ica_result = run_ica(dat; exclude_samples = [:is_extreme_value])
                    # save ica results
                    if config_data["files"]["output"]["save_ica_data"]
                        @info "Saving ica data"
                        jldsave(make_output_filename(output_data_directory, file, "_ica"); ica_result=ica_result)
                    end
                end


                # epoch data
                epochs = []
                for (idx, epoch) in enumerate([1, 4, 5, 3])
                    push!(epochs, extract_epochs(dat, idx, epoch, -2, 4))
                end

                # Log epoch counts
                df = DataFrame(file = basename(file), condition = 1:length(epochs), n_epochs_total = n_epochs.(epochs))
                push!(all_epoch_counts, df)  # Store the DataFrame
                @info "Epoch counts per condition:\n$(pretty_table(String, df, show_row_number=false, show_subheader=false))"

                # save epochs
                if config_data["files"]["output"]["save_epoch_data"]
                    @info "Saving epoch data"
                    jldsave(make_output_filename(output_data_directory, file, "_epochs"); epochs=epochs)
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
                if config_data["files"]["output"]["save_erp_data"]
                    @info "Saving erp data"
                    jldsave(make_output_filename(output_data_directory, file, "_erps"); erps=erps)
                end

                @info "Successfully processed $file"
                processed_files += 1

            catch e
                # Log error to both console and global log
                @error "Error processing $file" exception=(e, catch_backtrace())
                push!(failed_files, file)
            finally
                # Close per-file logging (restores global logger)
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
            jldsave(joinpath(output_data_directory, "epoch_summary.jld2"); df=combined_counts)
        end
        
    finally
       
        # Close global logging
        close_global_logging()
        
        # Move log file if needed
        mv("preprocess_eeg_data.log", joinpath(output_data_directory, "preprocess_eeg_data.log"), force=true)
    end
end
