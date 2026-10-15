"""Build a Symbol like `:is_extreme_value_100` from a base string and criterion."""
_flag_symbol(base::AbstractString, criterion) = Symbol("$(base)_$(criterion)")

"""Return `true` if all paths in `files` exist, warning for each missing file."""
function check_files_exist(files::Vector{String})
    all_files_exist = true
    for fname in files
        if !isfile(fname)
            @minimal_warning "File not found: $(fname)"
            all_files_exist = false
        end
    end
    return all_files_exist
end


"""Center a title between dashes for log formatting."""
function _center_title(title::String, width::Int)
    total_dashes = width - length(title) - 2  # 2 spaces around title
    left_dashes = div(total_dashes, 2)
    right_dashes = total_dashes - left_dashes
    return "-"^left_dashes * " $title " * "-"^right_dashes
end

"""Return a triple-line section header with `title` centered between dashes."""
function section(title::String; width::Int = 80)
    dash_line = "-"^width
    middle_line = _center_title(title, width)
    return "\n$dash_line\n$middle_line\n$dash_line"
end

"""Return a single-line sub-section header."""
subsection(title::String; width::Int = 80) = "\n" * _center_title(title, width)
"""Return a `# Title` sub-sub-section header."""
subsubsection(title::String) = "\n# " * title

"""
    check_raw_data(directory::String; raw_data_files::String="\\\\.(bdf|edf|vhdr|eeg|fif|xdf)\$", recursive::Bool=false)

Load raw EEG files (or headers where supported) in `directory` and display a summary of channels and sampling rate, warning if there are discrepancies.
"""
function check_raw_data(directory::String; raw_data_files::String = "\\.(bdf|edf|vhdr|eeg|fif|xdf)\$", recursive::Bool = false)
    files = get_files(directory, raw_data_files, recursive = recursive)
    isempty(files) && @minimal_error "No files found in $directory matching pattern $raw_data_files"

    @info "Checking $(length(files)) raw data files..."

    file_names = String[]
    n_chans = Int[]
    sampling_rates = Float64[]
    chan_labels_hash = UInt64[]
    file_sizes_mb = Float64[]
    durations_min = Float64[]

    for f in files
        push!(file_names, basename(f))

        # File size in MB
        push!(file_sizes_mb, round(filesize(f) / (1024 * 1024), digits = 2))

        ext = get_file_extension(f)
        if ext == ".bdf"
            hdr = read_raw_data(f, header_only = true)
            chans = Symbol.(hdr.channel_labels[1:(end-1)])
            sr = Float64(hdr.sample_rate[1])
            num_records = hdr.num_data_records
            if num_records <= 0
                n_chan_total = length(hdr.channel_labels)
                hdr_size = 256 * (n_chan_total + 1) # BDF header is 256 bytes per channel + 256 bytes for general info
                data_size = filesize(f) - hdr_size
                bps = n_chan_total * sr * 3
                num_records = floor(Int, data_size / bps)
                @minimal_warning "File $(basename(f)) has invalid duration in header. Estimated $(num_records) records from file size."
            end
            duration_sec = num_records * hdr.duration_data_records
        else
            # NOTE: Reads full data file. In the future, implementing header-only 
            # reading for these formats would significantly speed up data_check.
            dat_raw = read_raw_data(f)
            # The 1-argument create_eegfun_data builds a dummy layout matching the data channels
            dat = create_eegfun_data(dat_raw)
            chans = channel_labels(dat)
            sr = Float64(sample_rate(dat))
            # Fallback duration calculation from number of samples
            duration_sec = size(dat.data, 1) / sr
        end

        push!(n_chans, length(chans))
        push!(sampling_rates, sr)
        push!(chan_labels_hash, hash(chans))
        push!(durations_min, round(duration_sec / 60.0, digits = 2))
    end

    df = DataFrame(
        file = file_names,
        size_MB = file_sizes_mb,
        duration_min = durations_min,
        n_channels = n_chans,
        sampling_rate = sampling_rates,
        labels_hash = chan_labels_hash,
    )

    # Identify common baselines
    mode_chans = mode(n_chans)
    mode_sr = mode(sampling_rates)
    mode_hash = mode(chan_labels_hash)
    med_dur = isempty(durations_min) ? 0.0 : median(durations_min)

    # Add problem flag column
    issues = String[]
    for i = 1:nrow(df)
        has_issue = false
        if df.n_channels[i] != mode_chans || df.sampling_rate[i] != mode_sr || df.labels_hash[i] != mode_hash
            has_issue = true
        elseif med_dur > 0 && df.duration_min[i] < med_dur * 0.5
            has_issue = true
        end
        push!(issues, has_issue ? "*" : "")
    end

    df.issue = issues

    # Highlight rows with issues in red in the terminal output
    hl = PrettyTables.TextHighlighter((data, i, j) -> data[i, :issue] == "*", PrettyTables.Crayon(foreground = :red, bold = true))

    log_pretty_table(df; title = "Raw Data Files Metadata Check", highlighters = [hl])

    # Check for discrepancies
    if length(unique(n_chans)) > 1
        bad_files = df.file[df.n_channels .!= mode_chans]
        @minimal_warning "Discrepancy detected! Files have different numbers of channels: $(join(bad_files, ", "))"
    end
    if length(unique(sampling_rates)) > 1
        bad_files = df.file[df.sampling_rate .!= mode_sr]
        @minimal_warning "Discrepancy detected! Files have different sampling rates: $(join(bad_files, ", "))"
    end
    if length(unique(chan_labels_hash)) > 1
        bad_files = df.file[df.labels_hash .!= mode_hash]
        @minimal_warning "Discrepancy detected! Files have different channel labels: $(join(bad_files, ", "))"
    end

    # Flag files that are unusually short or small (e.g. < 50% of median)
    if !isempty(durations_min)
        med_dur = median(durations_min)
        for i = 1:nrow(df)
            if df.duration_min[i] < med_dur * 0.5
                @minimal_warning "File $(df.file[i]) is unusually short ($(df.duration_min[i]) min compared to median $(med_dur) min)."
            end
        end
    end

    if length(unique(n_chans)) == 1 && length(unique(sampling_rates)) == 1 && length(unique(chan_labels_hash)) == 1
        @info "All $(length(files)) files have consistent metadata!"
    end

    # Write to log file
    log_file = joinpath(directory, "check_raw_data.log")
    open(log_file, "w") do io
        println(io, "Raw Data Files Metadata Check")
        println(io, "=============================")
        show(io, df, allrows = true, allcols = true)
    end
    @info "Saved check results to: $log_file"

    return df
end
