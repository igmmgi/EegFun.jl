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


function basename_without_ext(path::String) 
    return splitext(basename(path))[1]
end

function make_output_filename(output_dir::String, input_file::String, suffix::String)
    base_name = basename_without_ext(input_file)
    return joinpath(output_dir, "$(base_name)$(suffix).jld2")
end


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
    # Set up the global log
    global_log = setup_global_logging("preprocess_eeg_data.log")
    
    # Initialize variable for outer scope
    output_directory = ""
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
        @info "Configuration loaded successfully"

        # check if all requested raw data files exist
        @info "Checking if raw data files exist"
        raw_data_files =
            get_files(config_data["files"]["input"]["directory"], config_data["files"]["input"]["raw_data_files"])
        raw_data_files_exist = check_files_exist(raw_data_files)
        if !raw_data_files_exist
            error_msg = "Missing raw data files requested within TOML file"
            @minimal_error error_msg
        end
        @info "Found $(length(raw_data_files)) raw data files: $(join(raw_data_files, ", "))"

        # We now need some code to read the epoch conditions defined within the toml file
        @info "Reading epoch conditions from $(config_data["files"]["input"]["epoch_condition_file"])"
        if !isfile(config_data["files"]["input"]["epoch_condition_file"])
            error_msg = "Epoch definition file does not exist: $(config_data["files"]["input"]["epoch_condition_file"])"
            @minimal_error error_msg
        end
        epoch_definition_file = TOML.parsefile(config_data["files"]["input"]["epoch_condition_file"])
        epoch_conditions = parse_epoch_conditions(epoch_definition_file)
        @info "Epoch conditions loaded successfully"

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
        output_directory = config_data["files"]["output"]["directory"]
        if !isdir(output_directory)
            mkdir(output_directory)
            @info "Created output directory: $output_directory"
        end

        # print config to output directory
        print_config(config_data, joinpath(output_directory, "config.toml"))

        # take layout file and add in 2D and 3D coordinated
        polar_to_cartesian_xy!(layout)
        polar_to_cartesian_xyz!(layout)

        neighbours_xy = get_electrode_neighbours_xy(layout, config_data["preprocess"]["layout"]["neighbour_criterion"])
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
                dat = read_bdf(file)
                dat = create_eeg_dataframe(dat, layout)

                # rereference the data
                rereference!(dat, Symbol(config_data["preprocess"]["reference_channel"]))

                # Save the results
                if config_data["files"]["output"]["save_continuous_data"]
                    @info "Saving continuous data"
                    jldsave(make_output_filename(output_directory, file, "_continuous"); dat=dat)
                end

                # initial high-pass filter 
                if config_data["filter"]["highpass"]["on"]
                    filter_data!(
                        dat,
                        "hp",
                        config_data["filter"]["highpass"]["type"],
                        config_data["filter"]["highpass"]["cutoff"],
                        order = config_data["filter"]["highpass"]["order"],
                    )
                end

                # initial low-pass filter 
                if config_data["filter"]["lowpass"]["on"]
                    filter_data!(
                        dat,
                        "lp",
                        config_data["filter"]["lowpass"]["type"],
                        config_data["filter"]["lowpass"]["cutoff"],
                        order = config_data["filter"]["lowpass"]["order"],
                    )
                end


                # caculate EOG channels
                diff_channel!(dat, [:Fp1, :Fp2], [:IO1, :IO2], :vEOG)
                diff_channel!(dat, :F9, :F10, :hEOG)

                # autodetect EOG signals
                detect_eog_onsets!(dat, config_data["preprocess"]["eog"]["vEOG_criterion"], :vEOG, :is_vEOG)
                detect_eog_onsets!(dat, config_data["preprocess"]["eog"]["hEOG_criterion"], :hEOG, :is_hEOG)

                # detect extreme values
                is_extreme_value!(dat, dat.layout.label, config_data["preprocess"]["eog"]["extreme_value_criterion"])

                # We perform the ica on "continuous" data (clean sections) that usually has a more extreme high-pass filter
                # run ica on clean sections of "continuous" data
                if config_data["ica"]["run"]

                    if config_data["ica"]["filter"]["highpass"]["on"]
                        # apply high-pass filter to data
                        dat_ica = filter_data(
                            dat,
                            "hp",
                            config_data["ica"]["filter"]["highpass"]["type"],
                            config_data["ica"]["filter"]["highpass"]["cutoff"],
                            order = config_data["ica"]["filter"]["highpass"]["order"],
                        )
                    else
                        dat_ica = copy(dat)
                    end

                    ica_result = run_ica(dat_ica; exclude_samples = [:is_extreme_value])

                    # save ica results
                    if config_data["files"]["output"]["save_ica_data"]
                        @info "Saving ica data"
                        jldsave(make_output_filename(output_directory, file, "_ica"); ica_result=ica_result)
                    end

                end


                # epoch data
                epochs = []
                for (idx, epoch_condition) in enumerate(epoch_conditions)
                    # TODO: Implement proper sequence matching
                    push!(
                        epochs,
                        extract_epochs(
                            dat,
                            idx,
                            epoch_condition.trigger_sequence,
                            config_data["preprocess"]["epoch_start"],
                            config_data["preprocess"]["epoch_end"];
                            zero_position = epoch_condition.reference_index,
                        ),
                    )
                end

                # Log epoch counts
                df = DataFrame(file = basename(file), condition = 1:length(epochs), n_epochs_total = n_epochs.(epochs))
                push!(all_epoch_counts, df)  # Store the DataFrame
                @info "Epoch counts per condition:\n$(pretty_table(String, df, show_row_number=false, show_subheader=false))"

                # save epochs
                if config_data["files"]["output"]["save_epoch_data"]
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
                if config_data["files"]["output"]["save_erp_data"]
                    @info "Saving erp data"
                    jldsave(make_output_filename(output_directory, file, "_erps"); erps=erps)
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
            jldsave(joinpath(output_directory, "epoch_summary.jld2"); df=combined_counts)
        end
        
    finally
       
        # Close global logging
        close_global_logging()
        
        # Move log file if needed
        mv("preprocess_eeg_data.log", joinpath(output_directory, "preprocess_eeg_data.log"), force=true)
    end
end

"""
    EpochCondition

Defines parameters for extracting epochs for a specific experimental condition.

# Fields
- `name::String`: Descriptive condition name
- `trigger_sequence::Vector{Int}`: Trigger sequence to match (e.g., [1, 2, 3])
- `reference_index::Int`: Which trigger position is t=0 (1-based, default: 1)
- `timing_pairs::Union{Nothing,Vector{Tuple{Int,Int}}}`: Which trigger pairs to apply min/max intervals to (optional, default: nothing)
- `min_interval::Union{Nothing,Float64}`: Minimum time between specified trigger pairs (optional, default: nothing)
- `max_interval::Union{Nothing,Float64}`: Maximum time between specified trigger pairs (optional, default: nothing)
"""
@kwdef struct EpochCondition
    name::String
    trigger_sequence::Vector{Int}
    reference_index::Int = 1
    timing_pairs::Union{Nothing,Vector{Tuple{Int,Int}}} = nothing
    min_interval::Union{Nothing,Float64} = nothing
    max_interval::Union{Nothing,Float64} = nothing
end

"""
    parse_epoch_conditions(config::Dict) -> Vector{EpochCondition}

Parse epoch conditions from configuration dictionary.
"""
function parse_epoch_conditions(config::Dict)
    defaults = get(config, "epochs", Dict())
    default_max_interval = get(defaults, "default_max_interval", 2.0)
    default_min_interval = get(defaults, "default_min_interval", 0.0)
    
    conditions = EpochCondition[]
    condition_configs = get(defaults, "conditions", [])
    
    for condition_config in condition_configs
        name = condition_config["name"]
        
        # Parse trigger sequence (single sequence per condition)
        trigger_sequence = Vector{Int}(condition_config["trigger_sequence"])
        
        # Parse reference index (single index) - default to 1
        reference_index = get(condition_config, "reference_index", 1)
        
        # Parse timing pairs (optional - if not specified, no timing constraints)
        timing_pairs_raw = get(condition_config, "timing_pairs", nothing)
        if timing_pairs_raw === nothing
            # No timing constraints
            timing_pairs = nothing
            min_interval = nothing
            max_interval = nothing
        else
            # Parse timing constraints
            timing_pairs = [(pair[1], pair[2]) for pair in timing_pairs_raw]
            min_interval = get(condition_config, "min_interval", default_min_interval)
            max_interval = get(condition_config, "max_interval", default_max_interval)
        end
        
        # Validation
        if reference_index < 1 || reference_index > length(trigger_sequence)
            error("reference_index must be between 1 and $(length(trigger_sequence)) for condition '$name'")
        end
        
        # Only validate timing constraints if they're specified
        if timing_pairs !== nothing
            if min_interval !== nothing && max_interval !== nothing && min_interval >= max_interval
                error("min_interval must be < max_interval for condition '$name'")
            end
            
            # Validate timing pairs
            for (start_idx, end_idx) in timing_pairs
                if start_idx < 1 || start_idx > length(trigger_sequence) || end_idx < 1 || end_idx > length(trigger_sequence)
                    error("timing_pairs contains invalid indices for sequence of length $(length(trigger_sequence)) in condition '$name'")
                end
                if start_idx >= end_idx
                    error("timing_pairs must have start_idx < end_idx for condition '$name'")
                end
            end
        end
        
        push!(conditions, EpochCondition(
            name, trigger_sequence, reference_index, timing_pairs, min_interval, max_interval
        ))
    end
    
    return conditions
end
