function _format_time(seconds::Real)
    seconds < 60 && return "$(round(seconds; digits = 1))s"
    if seconds < 3600
        return "$(floor(Int, seconds / 60))m $(floor(Int, seconds % 60))s"
    end
    return "$(floor(Int, seconds / 3600))h $(floor(Int, (seconds % 3600) / 60))m"
end

@inline function _save_pipeline_jld2(output_directory::String, data_file::String, suffix::String, obj, msg::String)
    @info msg
    jldsave(_make_output_filename(output_directory, data_file, suffix); data = obj, format_version = EEGFUN_FORMAT_VERSION)
end

"""
    preprocess(config::String; base_dir = nothing, log_level = :info, skip_existing = false, dry_run = false)

Run the EEG preprocessing pipeline defined by TOML file `config`.

# Keyword Arguments
- `base_dir::Union{String,Nothing} = nothing`: Base directory for relative paths in `config` (defaults to `dirname(abspath(config))`).
- `log_level::Symbol = :info`: Logging verbosity (`:debug`, `:info`, `:warn`, `:error`).
- `skip_existing::Bool = false`: Skip files whose final `_epochs.jld2` already exists.
- `dry_run::Bool = false`: Validate configuration, paths, and input files without writing outputs or processing data.
"""
function preprocess(
    config::String;
    base_dir::Union{String,Nothing} = nothing,
    log_level::Symbol = :info,
    skip_existing::Bool = false,
    dry_run::Bool = false,
)
    isnothing(base_dir) && (base_dir = abspath(dirname(abspath(config))))
    @info "Using base directory for relative paths: $base_dir"

    log_filename = "preprocess_log_$(Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")).txt"
    setup_global_logging(log_filename; log_level = log_level)

    output_directory_ref = Ref("")
    all_epoch_counts = DataFrame[]

    try
        setup = _pipeline_setup(config, base_dir, dry_run, output_directory_ref)
        output_directory = setup.output_directory
        raw_data_files = setup.raw_data_files

        if dry_run
            _pipeline_dry_run_summary(
                setup.input_directory,
                output_directory,
                setup.layout_file_path,
                setup.layout,
                setup.epoch_condition_file,
                setup.epoch_cfgs,
                raw_data_files,
                setup.preprocess_cfg,
                setup.files_cfg,
            )
            return nothing
        end

        processed_files = 0
        skipped_files = 0
        failed_files = String[]
        total_processing_time = 0.0
        n_processed_for_eta = 0
        n_total = length(raw_data_files)

        for (file_idx, data_file) in enumerate(raw_data_files)
            if skip_existing && isfile(_make_output_filename(output_directory, data_file, "_epochs"))
                @info "Skipping file $file_idx/$n_total: $(basename(data_file)) (already processed)"
                skipped_files += 1
                continue
            end

            file_start_time = time()
            eta_str =
                n_processed_for_eta > 0 ? _format_time((total_processing_time / n_processed_for_eta) * (n_total - file_idx + 1)) :
                "calculating..."
            @info "Processing file $file_idx/$n_total: $(basename(data_file)) (ETA: $eta_str)"

            try
                _pipeline_process_file!(
                    data_file,
                    setup.layout,
                    setup.preprocess_cfg,
                    setup.files_cfg,
                    setup.epoch_cfgs,
                    output_directory,
                    log_level,
                    all_epoch_counts,
                )
                file_elapsed = time() - file_start_time
                total_processing_time += file_elapsed
                n_processed_for_eta += 1
                @info "Successfully processed file $file_idx/$n_total: $(basename(data_file)) in $(_format_time(file_elapsed))"
                processed_files += 1
            catch e
                @minimal_stacktrace "Error processing $data_file" e 5
                push!(failed_files, data_file)
            finally
                close_logging()
            end
        end

        _pipeline_final_summary(processed_files, skipped_files, failed_files, output_directory, all_epoch_counts)
    finally
        close_global_logging()
        _pipeline_finalize_log(log_filename, output_directory_ref[])
    end
end

function _pipeline_setup(config::String, base_dir::String, dry_run::Bool, output_directory_ref::Ref{String})
    @info section("Setup")
    @info "Configuration Files:"
    cfg = read_config(config)
    isnothing(cfg) && @minimal_error "Failed to load configuration from: $config"

    resolve_path(path::String) = isabspath(path) ? path : joinpath(base_dir, path)

    pipeline_cfg = PipelineConfig(cfg)
    preprocess_cfg = pipeline_cfg.preprocess
    files_cfg = pipeline_cfg.files

    input_directory = resolve_path(files_cfg.input.directory)
    !isdir(input_directory) && @minimal_error "Input directory does not exist: $input_directory"

    output_directory = resolve_path(files_cfg.output.directory)
    output_directory_ref[] = output_directory
    !dry_run && !isdir(output_directory) && mkpath(output_directory)

    layout_file_path = resolve_path(files_cfg.input.layout_file)
    if !isfile(layout_file_path)
        layout_file = find_file(files_cfg.input.layout_file, joinpath(@__DIR__, "..", "..", "resources", "layouts"))
        !isnothing(layout_file) && (layout_file_path = layout_file)
    end
    !isfile(layout_file_path) && @minimal_error "Layout file not found: $layout_file_path"
    layout = read_layout(layout_file_path)

    epoch_condition_file = resolve_path(files_cfg.input.epoch_condition_file)

    raw_data_files = get_files(input_directory, files_cfg.input.raw_data_files; recursive = files_cfg.input.recursive)
    !check_files_exist(raw_data_files) && @minimal_error "Missing raw data files requested within TOML file!"
    isempty(raw_data_files) &&
        @minimal_error "No files found in '$input_directory' matching pattern '$(files_cfg.input.raw_data_files)'. Check the 'directory' and 'raw_data_files' settings in your pipeline.toml."
    @info "Found $(length(raw_data_files)) files: $(_print_vector(basename.(raw_data_files)))"

    !isfile(epoch_condition_file) && @minimal_error "File missing: $epoch_condition_file"
    epoch_cfgs = condition_parse_epoch(TOML.parsefile(epoch_condition_file))
    @info "Loading/parsing epoch file: $epoch_condition_file"
    @info "Output directory: $output_directory"

    !dry_run && print_config(cfg, joinpath(output_directory, "config.toml"))

    polar_to_cartesian_xy!(layout)
    polar_to_cartesian_xyz!(layout)
    get_neighbors_xy!(layout, preprocess_cfg.neighbor_criterion)
    !dry_run && print_layout_neighbors(layout, joinpath(output_directory, "neighbors.toml"))

    return (;
        cfg,
        preprocess_cfg,
        files_cfg,
        input_directory,
        output_directory,
        layout_file_path,
        layout,
        epoch_condition_file,
        epoch_cfgs,
        raw_data_files,
    )
end

function _pipeline_process_file!(
    data_file::String,
    layout::Layout,
    preprocess_cfg::PreprocessConfig,
    files_cfg::FilesConfig,
    epoch_cfgs,
    output_directory::String,
    log_level::Symbol,
    all_epoch_counts::Vector{DataFrame},
)
    @info section("Processing")
    @info "File: $data_file"
    setup_logging(joinpath(output_directory, "$(basename_without_ext(data_file))_log.txt"); log_level = log_level)

    @info section("Raw Data")
    dat = create_eegfun_data(read_raw_data(data_file), layout)

    missing_channels = setdiff(layout.data.label, channel_labels(dat))
    !isempty(missing_channels) && @minimal_error "Data file is missing channels defined in layout: $(join(missing_channels, ", "))"

    files_cfg.output.save_continuous_data_raw &&
        _save_pipeline_jld2(output_directory, data_file, "_continuous_raw", dat, "Saving continuous data (original)")

    @info section("Marking epoch intervals")
    @info "Epoch intervals: $([preprocess_cfg.epoch_start, preprocess_cfg.epoch_end])"
    mark_epoch_intervals!(dat, epoch_cfgs, [preprocess_cfg.epoch_start, preprocess_cfg.epoch_end])

    _pipeline_initial_filters!(dat, preprocess_cfg)
    hEOG_vEOG_cm_epoch = _pipeline_eog!(dat, preprocess_cfg)
    continuous_repair_info = _pipeline_continuous_artifacts!(dat, preprocess_cfg, hEOG_vEOG_cm_epoch)

    @info section("Rereference")
    if preprocess_cfg.reference_channel != :none
        rereference!(dat, preprocess_cfg.reference_channel)
    else
        @info "Skipping rereferencing (reference_channel = none)"
    end

    @info section("Initial Epoch Extraction (Pre-ICA)")
    epochs_raw = extract_epochs(dat, epoch_cfgs, (preprocess_cfg.epoch_start, preprocess_cfg.epoch_end))
    erps_raw = average_epochs(epochs_raw)

    files_cfg.output.save_epoch_data_raw &&
        _save_pipeline_jld2(output_directory, data_file, "_epochs_raw", epochs_raw, "Saving epoch data (original)")
    files_cfg.output.save_erp_data_raw &&
        _save_pipeline_jld2(output_directory, data_file, "_erps_raw", erps_raw, "Saving ERP data (original)")

    component_artifacts = _pipeline_ica!(dat, preprocess_cfg, continuous_repair_info, output_directory, data_file, files_cfg)

    _pipeline_epoch_and_reject!(
        dat,
        preprocess_cfg,
        epoch_cfgs,
        continuous_repair_info,
        component_artifacts,
        output_directory,
        data_file,
        files_cfg,
        epochs_raw,
        all_epoch_counts,
    )
    return nothing
end

function _pipeline_final_summary(
    processed_files::Int,
    skipped_files::Int,
    failed_files::Vector{String},
    output_directory::String,
    all_epoch_counts::Vector{DataFrame},
)
    @info section("Summary")
    @info "$processed_files success, $(length(failed_files)) fail$(skipped_files > 0 ? ", $skipped_files skipped" : "")"
    !isempty(failed_files) && @info "Failed files: $(join(failed_files, ", "))"

    @info subsection("Electrode Repair Summary Across All Participants (Continuous Level Only)")
    electrode_per_file, electrode_summary = summarize_electrode_repairs("_artifact_info"; input_dir = output_directory)
    !isempty(electrode_per_file) && log_pretty_table(electrode_per_file; title = "Electrode Repairs per Participant:")
    !isempty(electrode_summary) && log_pretty_table(electrode_summary; title = "Electrode Repairs: Number of Participants Affected")

    @info subsection("ICA Component Removal Summary")
    ica_per_file, ica_avg = summarize_ica_components("_artifact_info"; input_dir = output_directory)
    if !isempty(ica_per_file)
        log_pretty_table(ica_per_file; title = "ICA Components Removed per Participant", alignment = [:l, :r, :r, :r, :r, :r, :r])
        log_pretty_table(ica_avg; title = "Average ICA Components Removed per Participant", alignment = [:l, :r])
    end

    if !isempty(all_epoch_counts)
        epoch_summary, file_summary = _epoch_and_file_summary(all_epoch_counts)
        merged_epoch_summary, merged_file_summary = _merge_summaries(epoch_summary, file_summary, output_directory)
        log_pretty_table(merged_epoch_summary; title = "Combined epoch counts across all files:", alignment = [:l, :r, :l, :r, :r, :r])
        log_pretty_table(
            merged_file_summary;
            title = "Average percentage per condition (averaged across conditions):",
            alignment = [:l, :r],
        )
        @info "Mean percentage (averaged across all conditions and files): $(round(mean(merged_file_summary.percentage); digits = 1)) %"
        jldsave(joinpath(output_directory, "epoch_summary.jld2"); data = merged_epoch_summary, format_version = EEGFUN_FORMAT_VERSION)
        jldsave(joinpath(output_directory, "file_summary.jld2"); data = merged_file_summary, format_version = EEGFUN_FORMAT_VERSION)
    end
    return nothing
end

function _pipeline_finalize_log(log_filename::String, output_directory::String)
    if !isempty(output_directory) && isdir(output_directory) && isfile(log_filename)
        log_destination = joinpath(output_directory, log_filename)
        log_filename != log_destination && mv(log_filename, log_destination; force = true)
    elseif isfile(log_filename)
        rm(log_filename; force = true)
    end
    return nothing
end

function _pipeline_dry_run_summary(
    input_directory,
    output_directory,
    layout_file_path,
    layout,
    epoch_condition_file,
    epoch_cfgs,
    raw_data_files,
    preprocess_cfg::PreprocessConfig,
    files_cfg::FilesConfig,
)
    @info section("Dry Run Summary")
    @info "✓ Configuration loaded and validated successfully"
    @info "✓ Input directory exists: $input_directory"
    @info "✓ Output directory: $output_directory"
    @info "✓ Layout file loaded: $layout_file_path ($(length(layout.data.label)) channels)"
    @info "✓ Epoch condition file loaded: $epoch_condition_file ($(length(epoch_cfgs)) condition(s))"
    @info "✓ Found $(length(raw_data_files)) file(s) to process:"
    for (i, f) in enumerate(raw_data_files)
        status = isfile(_make_output_filename(output_directory, f, "_epochs")) ? "(already processed)" : "(pending)"
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
    for (label, flag) in (
        ("Continuous (raw)", files_cfg.output.save_continuous_data_raw),
        ("Continuous (corrected)", files_cfg.output.save_continuous_data_corrected),
        ("ICA data", files_cfg.output.save_ica_data),
        ("Epochs (raw)", files_cfg.output.save_epoch_data_raw),
        ("Epochs (corrected)", files_cfg.output.save_epoch_data_corrected),
        ("Epochs (final)", files_cfg.output.save_epoch_data),
        ("ERPs (raw)", files_cfg.output.save_erp_data_raw),
        ("ERPs (corrected)", files_cfg.output.save_erp_data_corrected),
        ("ERPs (final)", files_cfg.output.save_erp_data),
    )
        @info "  $(flag ? "✓" : "✗") $label"
    end
    @info ""
    @info "Dry run complete. No data was processed."
    return nothing
end

function _pipeline_initial_filters!(dat::ContinuousData, preprocess_cfg::PreprocessConfig)
    @info section("Initial Filters")
    highpass_filter!(dat, preprocess_cfg.filter)
    lowpass_filter!(dat, preprocess_cfg.filter)

    if preprocess_cfg.resample.apply
        @info section("Resampling")
        @info "Resampling data to $(preprocess_cfg.resample.target_rate) Hz"
        resample!(dat, preprocess_cfg.resample.target_rate)
    end

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
            pad = preprocess_cfg.cleanline.pad,
        )
    end
end

function _pipeline_eog!(dat::ContinuousData, preprocess_cfg::PreprocessConfig)
    @info section("EOG")
    @info subsection("Calculating EOG (vEOG/hEOG) channels")
    calculate_eog_channels!(dat, preprocess_cfg.eog)

    @info subsection("Detecting EOG (vEOG/hEOG) onsets")
    detect_eog_onsets!(dat, preprocess_cfg.eog)

    @info subsection("Channel x vEOG/hEOG Correlation Matrix")
    hEOG_vEOG_cm = correlation_matrix_eog(dat, preprocess_cfg.eog)
    add_zscore_columns!(hEOG_vEOG_cm)
    log_pretty_table(hEOG_vEOG_cm; title = "Channel x vEOG/hEOG Correlation Matrix (whole dataset)")

    @info subsection("Channel x vEOG/hEOG Correlation Matrix (epoch interval)")
    hEOG_vEOG_cm_epoch = correlation_matrix_eog(dat, preprocess_cfg.eog; sample_selection = samples(:epoch_interval))
    add_zscore_columns!(hEOG_vEOG_cm_epoch)
    log_pretty_table(hEOG_vEOG_cm_epoch; title = "Channel x vEOG/hEOG Correlation Matrix (epoch interval)")
    return hEOG_vEOG_cm_epoch
end

function _wait_interactive_figure(fig::Figure, screen = nothing)
    if !isnothing(screen) && applicable(wait, screen)
        wait(screen)
        return nothing
    end
    for scr in fig.scene.current_screens
        applicable(wait, scr) && wait(scr)
    end
    return nothing
end

function _pipeline_continuous_artifacts!(dat::ContinuousData, preprocess_cfg::PreprocessConfig, hEOG_vEOG_cm_epoch::DataFrame)
    @info section("Artifact Detection: Continuous Data")

    @info subsection("Channel Summary")
    summary_whole = channel_summary(dat)
    log_pretty_table(summary_whole; title = "Channel Summary (whole dataset)")
    summary_epoch = channel_summary(dat; sample_selection = samples(:epoch_interval))
    log_pretty_table(summary_epoch; title = "Channel Summary (epoch interval)")

    @info subsection("Artifact Detection (extreme values)")
    @info "Detecting extreme values: $(preprocess_cfg.eeg.extreme_value_abs_criterion) μV"
    is_extreme_value!(
        dat,
        preprocess_cfg.eeg.extreme_value_abs_criterion;
        channel_out = _flag_symbol("is_extreme_value", preprocess_cfg.eeg.extreme_value_abs_criterion),
    )

    @info subsection("Artifact Detection (criterion values)")
    @info "Detecting artifact values: $(preprocess_cfg.eeg.artifact_value_abs_criterion) μV"
    is_extreme_value!(
        dat,
        preprocess_cfg.eeg.artifact_value_abs_criterion;
        channel_out = _flag_symbol("is_artifact_value", preprocess_cfg.eeg.artifact_value_abs_criterion),
    )

    @info subsection("Bad Channel Detection using Channel Joint Probability + Z-Score Variance")
    cjp_whole = channel_joint_probability(dat)
    log_pretty_table(cjp_whole; title = "Channel Joint Probability (whole dataset)")
    cjp_epoch = channel_joint_probability(dat; sample_selection = samples(:epoch_interval))
    log_pretty_table(cjp_epoch; title = "Channel Joint Probability (epoch interval)")

    @info subsubsection("Bad Channels")
    bad_channels = intersect(identify_bad_channels(summary_whole, cjp_whole), identify_bad_channels(summary_epoch, cjp_epoch))

    bad_non_eog, bad_eog = partition_channels_by_eog_correlation(
        bad_channels,
        hEOG_vEOG_cm_epoch;
        eog_channels = [:hEOG, :vEOG],
        threshold = 0.3,
        use_z = false,
    )
    @info "Bad channels (non-EOG related): $(length(bad_non_eog)) channels - $(bad_non_eog)"
    @info "Bad channels (EOG related): $(length(bad_eog)) channels - $(bad_eog)"

    repair_method = preprocess_cfg.channel_repair_method
    channels_to_repair = copy(bad_non_eog)
    repair_info_name = "continuous_repair"

    if preprocess_cfg.interactive_continuous
        @info section("Interactive Continuous Data Review")
        @info "Reviewing continuous data before ICA. Close the window to continue."
        res = plot_databrowser(dat)
        _wait_interactive_figure(res.fig, res.screen)

        settings = res.analysis_settings[]
        apply_analysis_settings!(
            dat,
            AnalysisSettings(
                hp_filter = settings.hp_filter,
                lp_filter = settings.lp_filter,
                reference = settings.reference,
                repaired_channels = Symbol[],
                repair_method = :none,
                selected_regions = settings.selected_regions,
                removed_ica_components = settings.removed_ica_components,
            ),
        )

        if !isempty(settings.repaired_channels)
            @info "Manually marked channels for repair: $(settings.repaired_channels)"
            isempty(channels_to_repair) && (repair_info_name = "continuous_repair_manual")
            settings.repair_method != :none && (repair_method = settings.repair_method)
            append!(channels_to_repair, settings.repaired_channels)
            unique!(channels_to_repair)
        end
        @info "Continuous data review complete."
    end

    isempty(channels_to_repair) && return nothing

    continuous_repair_info = create_continuous_repair_info(repair_method; name = repair_info_name)
    channel_repairable!(continuous_repair_info, channels_to_repair, dat.layout)

    @info section("Channel Repair")
    repair_channels!(dat, continuous_repair_info; method = repair_method)
    @info continuous_repair_info
    return continuous_repair_info
end

function _pipeline_ica!(
    dat::ContinuousData,
    preprocess_cfg::PreprocessConfig,
    continuous_repair_info::Union{ContinuousRepairInfo,Nothing},
    output_directory::String,
    data_file::String,
    files_cfg::FilesConfig,
)
    preprocess_cfg.ica.apply || return nothing
    @info section("ICA")

    dat_ica = copy(dat)
    @info subsection("ICA filters")
    highpass_filter!(dat_ica, preprocess_cfg.filter; section = :ica_highpass)
    lowpass_filter!(dat_ica, preprocess_cfg.filter; section = :ica_lowpass)

    artifact_flags = [_flag_symbol("is_extreme_value", preprocess_cfg.eeg.extreme_value_abs_criterion)]
    hasproperty(dat_ica.data, :selected_region) && push!(artifact_flags, :selected_region)

    @info subsection("Running ICA")
    n_repaired = isnothing(continuous_repair_info) ? 0 : length(continuous_repair_info.repaired)
    n_ref = preprocess_cfg.reference_channel == :avg ? 1 : 0
    n_comps = max(1, n_channels(dat_ica) - n_repaired - n_ref)

    ica = run_ica(
        dat_ica;
        n_components = n_comps,
        sample_selection = samples_or_not(artifact_flags),
        percentage_of_data = preprocess_cfg.ica.percentage_of_data,
    )

    @info subsection("Component Identification")
    component_artifacts, component_metrics =
        identify_components(dat, ica; method = preprocess_cfg.ica.component_method, sample_selection = samples_or_not(artifact_flags))

    log_pretty_table(component_metrics[:eog_metrics]; title = "EOG Component Metrics")
    log_pretty_table(component_metrics[:ecg_metrics]; title = "ECG Component Metrics")
    log_pretty_table(component_metrics[:line_noise_metrics]; title = "Line Noise Component Metrics")
    log_pretty_table(component_metrics[:channel_noise_metrics]; title = "Channel Noise Component Metrics")

    if preprocess_cfg.interactive_ica
        @info section("Interactive ICA Component Review")
        @info "Use the databrowser ICA menu to review and manually select components to remove. Close the window to continue."
        res = plot_databrowser(dat, ica)
        _wait_interactive_figure(res.fig, res.screen)

        manual_removals = res.analysis_settings[].removed_ica_components
        if !isempty(manual_removals)
            @info "Manually removed components: $manual_removals"
            append!(component_artifacts.manual, manual_removals)
            unique!(component_artifacts.manual)
        end
    end

    @info subsection("Removing ICA components")
    all_removed_components = get_all_ica_components(component_artifacts)
    @info "Removed $(length(all_removed_components)) ICA components" component_artifacts
    !isempty(all_removed_components) && subtract_ica_components!(dat, ica; component_selection = components(all_removed_components))

    files_cfg.output.save_ica_data && _save_pipeline_jld2(output_directory, data_file, "_ica", ica, "Saving ica data")
    return component_artifacts
end

function _pipeline_epoch_and_reject!(
    dat::ContinuousData,
    preprocess_cfg::PreprocessConfig,
    epoch_cfgs::Vector{EpochCondition},
    continuous_repair_info::Union{ContinuousRepairInfo,Nothing},
    component_artifacts::Union{ArtifactComponents,Nothing},
    output_directory::String,
    data_file::String,
    files_cfg::FilesConfig,
    epochs_raw::Vector{EpochData},
    all_epoch_counts::Vector{DataFrame},
)
    @info section("EOG Recalculation")
    @info subsection("Recalculating EOG (vEOG/hEOG) channels after ICA and repair")
    calculate_eog_channels!(dat, preprocess_cfg.eog)

    @info section("Detecting artifact values in continuous data")
    is_extreme_value!(
        dat,
        preprocess_cfg.eeg.artifact_value_abs_criterion;
        channel_out = _flag_symbol("is_artifact_value", preprocess_cfg.eeg.artifact_value_abs_criterion),
    )

    files_cfg.output.save_continuous_data_corrected &&
        _save_pipeline_jld2(output_directory, data_file, "_continuous_corrected", dat, "Saving continuous data")

    @info section("Extracting cleaned epoched data")
    epochs = extract_epochs(dat, epoch_cfgs, (preprocess_cfg.epoch_start, preprocess_cfg.epoch_end))

    empty_conds = [ep.condition_name for ep in epochs if isempty(ep.data)]
    !isempty(empty_conds) && @minimal_error(
        "Epoch extraction resulted in empty epochs for conditions: $(join(empty_conds, ", ")). Check epoch interval parameters and trigger locations.",
    )

    @info section("Baseline whole epochs")
    baseline!(epochs)

    rej_interval = if !isnothing(preprocess_cfg.eeg.artifact_interval_start) && !isnothing(preprocess_cfg.eeg.artifact_interval_end)
        @info "Restricting artifact rejection to window: $(preprocess_cfg.eeg.artifact_interval_start)s to $(preprocess_cfg.eeg.artifact_interval_end)s"
        times(preprocess_cfg.eeg.artifact_interval_start, preprocess_cfg.eeg.artifact_interval_end)
    else
        times()
    end

    local rejection_info_step1::Vector{EpochRejectionInfo}
    local rejection_info_step2::Vector{EpochRejectionInfo}

    if preprocess_cfg.autoreject.apply
        @info section("Adaptive Artifact Handling (Autoreject)")
        autoreject_res = autoreject!(
            epochs;
            mode = preprocess_cfg.autoreject.mode,
            interval_selection = rej_interval,
            n_folds = preprocess_cfg.autoreject.n_folds,
            n_candidates = preprocess_cfg.autoreject.n_candidates,
            consensus_perc = preprocess_cfg.autoreject.consensus_perc,
            max_repair_perc = preprocess_cfg.autoreject.max_repair_perc,
            repair_method = preprocess_cfg.autoreject.repair_method,
            seed = preprocess_cfg.autoreject.seed,
            reject = false,
        )
        log_pretty_table(autoreject_res.summary_table; title = "Autoreject: Summary")
        log_pretty_table(autoreject_res.channel_table; title = "Autoreject: Learned Channel Thresholds")
        rejection_info_step1 = autoreject_res.rejection_info
        rejection_info_step2 = autoreject_res.rejection_info
    else
        @info section("Automatic epoch detection")
        rejection_info_step1 = detect_bad_epochs_automatic(
            epochs;
            z_criterion = preprocess_cfg.eeg.artifact_value_z_criterion,
            abs_criterion = preprocess_cfg.eeg.artifact_value_abs_criterion,
            interval_selection = rej_interval,
            name = "rejection_step1",
        )
        channel_repairable!(rejection_info_step1, epochs[1].layout)
        @info ""
        @info rejection_info_step1

        @info section("Channel Repair per Epoch")
        repair_artifacts!(epochs, rejection_info_step1; method = preprocess_cfg.channel_repair_method)

        @info subsection("Re-detecting artifacts after repair")
        rejection_info_step2 = detect_bad_epochs_automatic(
            epochs;
            z_criterion = preprocess_cfg.eeg.artifact_value_z_criterion,
            abs_criterion = preprocess_cfg.eeg.artifact_value_abs_criterion,
            interval_selection = rej_interval,
            name = "rejection_step2",
        )
        channel_repairable!(rejection_info_step2, epochs[1].layout)
        @info ""
        @info rejection_info_step2

        @info subsection("Rejection Step Comparison (before vs after repair)")
        log_pretty_table(
            compare_rejections(rejection_info_step1, rejection_info_step2);
            title = "Rejection Step Comparison: Effectiveness of Channel Repair",
        )
    end

    files_cfg.output.save_epoch_data_corrected &&
        _save_pipeline_jld2(output_directory, data_file, "_epochs_corrected", epochs, "Saving epoch data (cleaned)")

    if preprocess_cfg.interactive_epochs
        @info section("Interactive Epoch Rejection Review")
        @info "Review epoch artifacts. Checked epochs will be rejected. Close the window to continue."
        states = detect_bad_epochs_interactive(epochs; artifact_info = rejection_info_step2)
        for st in (states isa AbstractVector ? states : (states,))
            _wait_interactive_figure(st.fig)
        end
        rejection_info_step2 = states isa AbstractVector ? _to_rejection_info(states) : [_to_rejection_info(states)]
    end

    files_cfg.output.save_erp_data_corrected &&
        _save_pipeline_jld2(output_directory, data_file, "_erps_corrected", average_epochs(epochs), "Saving ERP data (cleaned)")

    @info subsection("Rejecting bad epochs")
    epochs = reject_epochs(epochs, rejection_info_step2)

    @info subsection("Artifact Information")
    artifact_info = ArtifactInfo(
        !isnothing(continuous_repair_info) ? [continuous_repair_info] : ContinuousRepairInfo[],
        vcat(rejection_info_step1, rejection_info_step2),
        component_artifacts,
    )
    jldsave(
        _make_output_filename(output_directory, data_file, "_artifact_info");
        data = artifact_info,
        format_version = EEGFUN_FORMAT_VERSION,
    )
    @info "Saved artifact info: $(artifact_info)"

    push!(all_epoch_counts, log_epochs_table(epochs_raw, epochs; title = "Epoch counts per condition (after repair and rejection):"))

    files_cfg.output.save_epoch_data && _save_pipeline_jld2(output_directory, data_file, "_epochs", epochs, "Saving epoch data (good)")
    files_cfg.output.save_erp_data &&
        _save_pipeline_jld2(output_directory, data_file, "_erps", average_epochs(epochs), "Saving ERP data (good)")

    @info section("End of Processing")
end

"""
    validate_config(config::String; base_dir::Union{String,Nothing} = nothing)

Validate a preprocessing pipeline configuration file (`dry_run = true`) without processing data or writing output files.
"""
function validate_config(config::String; base_dir::Union{String,Nothing} = nothing)
    preprocess(config; base_dir = base_dir, log_level = :info, dry_run = true)
end
