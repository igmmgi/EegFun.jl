function preprocess_eeg_data(config::String)

    # check if config file exists and load
    if !isfile(config)
        @minimal_error "Config file does not exist: $config"
    end
    config = load_config(config)

    # merge with the default config
    default_config = load_config(joinpath(@__DIR__, "..", "..", "src", "config", "default.toml"))
    config = merge_configs(default_config, config)

    print_config(config)

    # check if all requested raw data files exist
    raw_data_files =
        get_files(config["files"]["input"]["raw_data_directory"], config["files"]["input"]["raw_data_files"])
    raw_data_files_exist = check_files_exist(raw_data_files)
    if !raw_data_files_exist
        @minimal_error "Missing raw data files requested within TOML file"
    end
    @info "Found $(length(raw_data_files)) raw data files: $(join(raw_data_files, ", "))"

    # Check if electrode configuration file exists and load
    layout_file = joinpath(@__DIR__, "..", "..", "data", "layouts", config["files"]["input"]["layout_file"])
    if !isfile(layout_file)
        # not default layout file so check if custom file specified
        layout_file = config["files"]["input"]["layout_file"]
        if !isfile(layout_file)
            @minimal_error "Electrode configuration file does not exist: $layout_file"
        end
    end
    layout = read_layout(layout_file)

    # Check if requested output directory exists and if not, create it
    output_data_directory = config["files"]["output"]["output_data_directory"]
    if !isdir(output_data_directory)
        mkdir(output_data_directory)
    end
    # delete current files in output directory
    if config["files"]["output"]["delete_current_files"]
        @info "Deleting current files in output directory: $output_data_directory"
        rm.(readdir(output_data_directory, join=true), recursive=true)
    end

    # print config to output directory
    print_config(config, joinpath(output_data_directory, "config.txt"))

    # take layout file and add in 2D and 3D coordinated
    polar_to_cartesian_xy!(layout)
    polar_to_cartesian_xyz!(layout)

    for file in raw_data_files
        @info "Processing $file"

        # read raw data file and create our Julia DataFrame
        dat = read_bdf(file)
        dat = create_eeg_dataframe(dat, layout)

        # rereference the data
        rereference!(dat, Symbol(config["preprocessing"]["reference_electrode"]))

        # initial high-pass filter to remove slow drifts
        filter_data!(
            dat,
            "hp",
            config["filtering"]["highpass"]["type"]["value"],
            config["filtering"]["highpass"]["cutoff"]["value"],
            order = config["filtering"]["highpass"]["order"]["value"],
        )

        # caculate EOG channels
        diff_channel!(dat, [:Fp1, :Fp2], [:IO1, :IO2], :vEOG)
        diff_channel!(dat, :F9, :F10, :hEOG)

        # autodetect EOG signals
        detect_eog_onsets!(dat, 50, :vEOG, :is_vEOG)
        detect_eog_onsets!(dat, 30, :hEOG, :is_hEOG)

        # detect extreme values
        is_extreme_value!(dat, dat.layout.label, 100)

        # run ica
        if config["ica"]["ica_run"]
            ica_result = run_ica(dat; exclude_samples = [:is_extreme_value])

            # save ica results
            if config["files"]["output"]["save_ica_data"]
                save_object("$(output_data_directory)/$(splitext(basename(file))[1])_ica.jld2", ica_result)
            end

        end

        # Save the results
        if config["files"]["output"]["save_continuous_data"]
            save_object("$(output_data_directory)/$(splitext(basename(file))[1])_continuous.jld2", dat)
        end

        # epoch data
        epochs = []
        for (idx, epoch) in enumerate([1, 4, 5, 3])
            push!(epochs, extract_epochs(dat, idx, epoch, -2, 4))
        end

        # save epochs
        if config["files"]["output"]["save_epoch_data"]
            save_object("$(output_data_directory)/$(splitext(basename(file))[1])_epochs.jld2", epochs)
        end

        # average epochs
        erps = []
        for (idx, epoch) in enumerate(epochs)
            push!(erps, average_epochs(epochs[idx]))
        end

        # save erps
        if config["files"]["output"]["save_erp_data"]
            save_object("$(output_data_directory)/$(splitext(basename(file))[1])_erps.jld2", erps)
        end

    end

end

