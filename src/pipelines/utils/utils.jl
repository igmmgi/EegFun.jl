_flag_symbol(base::AbstractString, criterion) = Symbol(base, '_', criterion)

function check_files_exist(files::Vector{String})
    all_exist = true
    for fname in files
        if !isfile(fname)
            @minimal_warning "File not found: $(fname)"
            all_exist = false
        end
    end
    return all_exist
end

function _center_title(title::String, width::Int)
    total_dashes = max(0, width - length(title) - 2)
    left_dashes = div(total_dashes, 2)
    return "-"^left_dashes * " $title " * "-"^(total_dashes - left_dashes)
end

"""
    section(title::String; width::Int = 80) -> String

Format a triple-line section header with `title` centered between dashes.
"""
function section(title::String; width::Int = 80)
    dash_line = "-"^width
    return "\n$dash_line\n$(_center_title(title, width))\n$dash_line"
end

subsection(title::String; width::Int = 80) = "\n" * _center_title(title, width)
subsubsection(title::String) = "\n# " * title

"""
    check_raw_data(directory::String; raw_data_files::String = "\\\\.(bdf|edf|vhdr|eeg|fif|xdf)\$", recursive::Bool = false) -> DataFrame

Inspect raw EEG files in `directory`, summarize channel counts, sampling rates, and durations,
and warn on metadata discrepancies across files.
"""
function check_raw_data(directory::String; raw_data_files::String = "\\.(bdf|edf|vhdr|eeg|fif|xdf)\$", recursive::Bool = false)
    files = get_files(directory, raw_data_files; recursive = recursive)
    isempty(files) && @minimal_error "No files found in $directory matching pattern $raw_data_files"

    n_files = length(files)
    @info "Checking $n_files raw data files..."

    file_names = Vector{String}(undef, n_files)
    file_sizes_mb = Vector{Float64}(undef, n_files)
    durations_min = Vector{Float64}(undef, n_files)
    n_chans = Vector{Int}(undef, n_files)
    sampling_rates = Vector{Float64}(undef, n_files)
    chan_labels_hash = Vector{UInt64}(undef, n_files)

    for (idx, f) in enumerate(files)
        file_names[idx] = basename(f)
        file_sizes_mb[idx] = round(filesize(f) / (1024 * 1024); digits = 2)

        if get_file_extension(f) == ".bdf"
            hdr = read_raw_data(f; header_only = true)
            chans = Symbol.(hdr.channel_labels[1:(end-1)])
            sr = Float64(hdr.sample_rate[1])
            num_records = hdr.num_data_records
            if num_records <= 0
                n_chan_total = length(hdr.channel_labels)
                hdr_size = 256 * (n_chan_total + 1)
                num_records = floor(Int, (filesize(f) - hdr_size) / (n_chan_total * sr * 3))
                @minimal_warning "File $(basename(f)) has invalid duration in header. Estimated $(num_records) records from file size."
            end
            duration_sec = num_records * hdr.duration_data_records
        else
            dat = create_eegfun_data(read_raw_data(f))
            chans = channel_labels(dat)
            sr = Float64(sample_rate(dat))
            duration_sec = size(dat.data, 1) / sr
        end

        n_chans[idx] = length(chans)
        sampling_rates[idx] = sr
        chan_labels_hash[idx] = hash(chans)
        durations_min[idx] = round(duration_sec / 60.0; digits = 2)
    end

    mode_chans = mode(n_chans)
    mode_sr = mode(sampling_rates)
    mode_hash = mode(chan_labels_hash)
    med_dur = isempty(durations_min) ? 0.0 : median(durations_min)

    issues = [
        (
            n_chans[i] != mode_chans ||
            sampling_rates[i] != mode_sr ||
            chan_labels_hash[i] != mode_hash ||
            (med_dur > 0 && durations_min[i] < med_dur * 0.5)
        ) ? "*" : "" for i = 1:n_files
    ]

    df = DataFrame(
        file = file_names,
        size_MB = file_sizes_mb,
        duration_min = durations_min,
        n_channels = n_chans,
        sampling_rate = sampling_rates,
        labels_hash = chan_labels_hash,
        issue = issues,
    )

    hl = PrettyTables.TextHighlighter((data, i, j) -> data[i, :issue] == "*", PrettyTables.Crayon(foreground = :red, bold = true))
    log_pretty_table(df; title = "Raw Data Files Metadata Check", highlighters = [hl])

    ch_mismatch = any(!=(mode_chans), n_chans)
    sr_mismatch = any(!=(mode_sr), sampling_rates)
    hash_mismatch = any(!=(mode_hash), chan_labels_hash)

    ch_mismatch &&
        @minimal_warning "Discrepancy detected! Files have different numbers of channels: $(join(file_names[n_chans .!= mode_chans], ", "))"
    sr_mismatch &&
        @minimal_warning "Discrepancy detected! Files have different sampling rates: $(join(file_names[sampling_rates .!= mode_sr], ", "))"
    hash_mismatch &&
        @minimal_warning "Discrepancy detected! Files have different channel labels: $(join(file_names[chan_labels_hash .!= mode_hash], ", "))"

    if med_dur > 0
        for i = 1:n_files
            if durations_min[i] < med_dur * 0.5
                @minimal_warning "File $(file_names[i]) is unusually short ($(durations_min[i]) min compared to median $(med_dur) min)."
            end
        end
    end

    if !ch_mismatch && !sr_mismatch && !hash_mismatch
        @info "All $n_files files have consistent metadata!"
    end

    log_file = joinpath(directory, "check_raw_data.log")
    open(log_file, "w") do io
        println(io, "Raw Data Files Metadata Check")
        println(io, "=============================")
        show(io, df; allrows = true, allcols = true)
    end
    @info "Saved check results to: $log_file"
    return df
end
