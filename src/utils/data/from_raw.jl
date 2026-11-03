# === DATA CREATION FUNCTIONS ===

"""
    _create_layout_from_labels(labels::Vector{Symbol})::Layout

Creates a fake layout with channel labels and zero positions for quick visualization.
Useful when no proper electrode layout is available.
"""
function _create_layout_from_labels(labels::Vector{Symbol})::Layout
    n_channels = length(labels)
    df = DataFrame(label = labels, inc = zeros(n_channels), azi = zeros(n_channels))
    return Layout(df, nothing, nothing, nothing)
end

"""Create a time+sample+trigger+channel DataFrame from raw BiosemiData."""
function _create_eegfun_dataframe(dat::BiosemiDataFormat.BiosemiData)::DataFrame
    df = DataFrame(time = dat.time, sample = 1:length(dat.time), trigger = _clean_triggers(dat.triggers.raw))

    channels = Symbol.(dat.header.channel_labels[1:(end-1)])
    for (i, ch) in enumerate(channels)
        # Allocate Float64 vector directly into the DataFrame column
        df[!, ch] = Float64.(@view dat.data[:, i])
    end

    return df
end


function create_eegfun_data(dat::BiosemiDataFormat.BiosemiData, layout::Layout)::ContinuousData
    @info "Creating EEG DataFrame from Biosemi data"
    file_name = filename(dat)
    df = _create_eegfun_dataframe(dat)
    sr = Int(dat.header.sample_rate[1])
    return ContinuousData(file_name, df, copy(layout), sr, AnalysisInfo(sample_rate = Float64(sr)))
end

# Optional Layout variant for quick databrowser visualization
function create_eegfun_data(dat::BiosemiDataFormat.BiosemiData)
    channel_labels = Symbol.(dat.header.channel_labels[1:(end-1)])
    layout = _create_layout_from_labels(channel_labels)
    return create_eegfun_data(dat, layout)
end

_is_edf_status_label(lbl::AbstractString) = lowercase(strip(lbl)) in ("status", "edf annotations", "tal")

function _unit_to_microvolt_scale(unit_str::AbstractString)::Float64
    u = strip(unit_str)
    isempty(u) && return 1.0
    ul = lowercase(u)
    if u == "V" || ul in ("v", "volt", "volts")
        return 1e6
    elseif u == "mV" || ul in ("mv", "millivolt", "millivolts")
        return 1e3
    else
        return 1.0
    end
end

function _edf_data_channels(dat::EuropeanDataFormat.EdfData)
    hdr_labels = dat.header.channel_labels
    hdr_units = dat.header.channel_unit
    n_data_cols = size(dat.data, 2)

    if n_data_cols == length(hdr_labels)
        keep_data = [i for (i, lbl) in enumerate(hdr_labels) if !_is_edf_status_label(lbl)]
        labels = Symbol.(strip.(hdr_labels[keep_data]))
        scales = [i <= length(hdr_units) ? _unit_to_microvolt_scale(hdr_units[i]) : 1.0 for i in keep_data]
        return keep_data, labels, scales
    else
        non_status_hdr = [i for (i, lbl) in enumerate(hdr_labels) if !_is_edf_status_label(lbl)]
        if length(non_status_hdr) == n_data_cols
            labels = Symbol.(strip.(hdr_labels[non_status_hdr]))
            scales = [i <= length(hdr_units) ? _unit_to_microvolt_scale(hdr_units[i]) : 1.0 for i in non_status_hdr]
            return collect(1:n_data_cols), labels, scales
        else
            sel_hdr = length(non_status_hdr) >= n_data_cols ? non_status_hdr[1:n_data_cols] : collect(1:n_data_cols)
            labels = Symbol.(strip.(hdr_labels[sel_hdr]))
            scales = [i <= length(hdr_units) ? _unit_to_microvolt_scale(hdr_units[i]) : 1.0 for i in sel_hdr]
            return collect(1:n_data_cols), labels, scales
        end
    end
end

"""Create a time+sample+trigger+channel DataFrame from raw EdfData."""
function _create_eegfun_dataframe(dat::EuropeanDataFormat.EdfData)::DataFrame

    n_samples = size(dat.data, 1)
    sample_rate = dat.header.sample_rate[1]

    time = collect(0:(n_samples-1)) ./ sample_rate
    sample = 1:n_samples

    trigger = zeros(Int, n_samples)
    trigger_info = fill("", n_samples)

    if !isnothing(dat.triggers) && length(dat.triggers.idx) > 0
        for i = 1:length(dat.triggers.idx)
            sample_idx = dat.triggers.idx[i]
            if 1 <= sample_idx <= n_samples
                trigger[sample_idx] = dat.triggers.val[i]
                trigger_info[sample_idx] = string(dat.triggers.val[i])
            end
        end
    end

    data_cols, channel_labels, scales = _edf_data_channels(dat)
    scaled_data = Float64.(@view dat.data[:, data_cols]) .* scales'

    df = hcat(
        DataFrame(time = time, sample = sample, trigger = trigger, trigger_info = trigger_info),
        DataFrame(scaled_data, channel_labels),
    )
    return df
end

"""
    download_eegfun_datasets()

Automatically download and cache built-in tutorial datasets using `DataDeps.jl`.
Returns the absolute file path to the directory containing the downloaded datasets.

# Examples
```julia
data_dir = download_eegfun_datasets()
file_path = joinpath(data_dir, "participant1.bdf")
dat = read_raw_data(file_path)
```
"""
function download_eegfun_datasets()
    return datadep"TutorialDataSets"
end

function create_eegfun_data(dat::EuropeanDataFormat.EdfData, layout::Layout)::ContinuousData
    @info "Creating EEG DataFrame (*.edf)"
    file_name = basename_without_ext(dat.filename)
    df = _create_eegfun_dataframe(dat)
    sr = Int(dat.header.sample_rate[1])
    return ContinuousData(file_name, df, copy(layout), sr, AnalysisInfo(sample_rate = Float64(sr)))
end

# Optional Layout variant for quick databrowser visualization
function create_eegfun_data(dat::EuropeanDataFormat.EdfData)
    _, channel_labels, _ = _edf_data_channels(dat)
    layout = _create_layout_from_labels(channel_labels)
    return create_eegfun_data(dat, layout)
end

"""Create a time+sample+trigger+channel DataFrame from raw BrainVision data."""
function _create_eegfun_dataframe(dat::BrainVisionDataFormat.BrainVisionData)::DataFrame

    # Check if data is available
    if isnothing(dat.data)
        @minimal_error "BrainVision data is empty (data field is nothing)"
    end

    if isnothing(dat.header)
        @minimal_error "BrainVision header is empty (header field is nothing)"
    end

    # Extract basic information
    # BrainVision data is now samples × channels (after package modification)
    n_samples = size(dat.data, 1)
    n_channels = size(dat.data, 2)
    sample_rate = dat.header.Fs

    # Create time vector
    time = collect(0:(n_samples-1)) ./ sample_rate

    # Create sample vector
    sample = 1:n_samples

    # Extract channel labels from header 
    channel_labels = dat.header.label[1:n_channels]

    # Verify channel count matches
    if length(channel_labels) != n_channels
        @minimal_error "Number of channel labels ($(length(channel_labels))) does not match number of channels ($n_channels)"
    end

    # Create trigger and marker strings columns from markers
    if isnothing(dat.markers) || isempty(dat.markers)
        trigger = zeros(Int, n_samples)
        trigger_info = fill("", n_samples)
    else
        trigger, trigger_info = _extract_triggers_from_markers(dat.markers, n_samples)
    end

    # Create the DataFrame 
    df = hcat(DataFrame(time = time, sample = sample, trigger = trigger, trigger_info = trigger_info), DataFrame(dat.data, channel_labels))

    return df
end

"""
    create_eegfun_data(dat::BrainVisionDataFormat.BrainVisionData, layout::Layout)::ContinuousData

Creates a ContinuousData object from a BrainVisionDataFormat data structure and a layout.

# Arguments
- `dat::BrainVisionDataFormat.BrainVisionData`: The BrainVisionDataFormat data structure containing EEG data.
- `layout::Layout`: The layout object containing electrode information.

# Returns
- `ContinuousData`: ContinuousData object containing the EEG data and layout information.

# Examples
```julia
# Create ContinuousData from BrainVisionDataFormat data and layout
eeg_data = create_eegfun_data(brainvision_data, layout)

# Quick visualization without layout (for databrowser only)
eeg_data = create_eegfun_data(brainvision_data)
```
"""
function create_eegfun_data(dat::BrainVisionDataFormat.BrainVisionData, layout::Layout)::ContinuousData
    @info "Creating EEG DataFrame (*.vhdr/*.vmrk/*.eeg)"
    file_name = basename_without_ext(dat.filename)
    df = _create_eegfun_dataframe(dat)
    sr = round(Int, dat.header.Fs)
    return ContinuousData(file_name, df, copy(layout), sr, AnalysisInfo(sample_rate = Float64(dat.header.Fs)))
end

# Optional Layout variant for quick databrowser visualization
function create_eegfun_data(dat::BrainVisionDataFormat.BrainVisionData)
    channel_labels = Symbol.(dat.header.label[1:size(dat.data, 2)])
    layout = _create_layout_from_labels(channel_labels)
    return create_eegfun_data(dat, layout)
end

function _parse_marker_code(val::AbstractString)::Int
    s = strip(val)
    isempty(s) && return 0
    direct = tryparse(Int, s)
    !isnothing(direct) && return direct
    m = match(r"-?\d+", s)
    return !isnothing(m) ? something(tryparse(Int, m.match), 0) : 0
end

"""
    _extract_triggers_from_markers(markers::Vector{BrainVisionDataFormat.BrainVisionMarker}, n_samples::Int)::Tuple{Vector{Int}, Vector{String}}

Extract trigger values from BrainVision markers and create trigger and marker string vectors.
Numeric codes are parsed directly from marker values (e.g. `"S 10"` -> `10`, `"101"` -> `101`);
non-numeric markers (e.g. `"New Segment"`) receive `0` in `:trigger` while preserving their text in `:trigger_info`.
"""
function _extract_triggers_from_markers(
    markers::Vector{BrainVisionDataFormat.BrainVisionMarker},
    n_samples::Int,
)::Tuple{Vector{Int},Vector{String}}
    trigger = zeros(Int, n_samples)
    trigger_info = fill("", n_samples)

    @debug "Processing $(length(markers)) markers for $n_samples samples"

    if !isempty(markers)
        first_few_samples = [marker.sample for marker in markers[1:min(5, length(markers))]]
        @debug "First few marker sample indices: $first_few_samples"
    end

    for marker in markers
        sample_idx = marker.sample
        # BrainVision markers are 1-based; sample == 0 denotes global/non-positioned markers
        sample_idx == 0 && continue
        if 1 <= sample_idx <= n_samples
            code = _parse_marker_code(marker.value)
            info_str = isempty(strip(marker.value)) ? strip(marker.type) : marker.value
            if code != 0 || trigger[sample_idx] == 0
                trigger[sample_idx] = code
            end
            if code != 0 || isempty(trigger_info[sample_idx])
                trigger_info[sample_idx] = info_str
            end
        else
            @minimal_warning "Marker sample $sample_idx out of bounds (1:$n_samples), skipping"
        end
    end

    non_zero_triggers = count(x -> x != 0, trigger)
    non_empty_strings = count(x -> x != "", trigger_info)
    @debug "Created trigger with $non_zero_triggers non-zero values and $non_empty_strings non-empty marker strings"

    return trigger, trigger_info
end

_is_stim_channel(ch) = ch.kind == 3 || startswith(uppercase(ch.ch_name), "STI") || uppercase(ch.ch_name) == "STATUS"

"""Create a time+sample+trigger+channel DataFrame from raw FifData."""
function _create_eegfun_dataframe(dat::FunctionalImageFormat.FifData)::DataFrame
    _, n_samples = size(dat.data)
    sample_rate = Float64(dat.header.sfreq)

    time = collect(0:(n_samples-1)) ./ sample_rate
    sample = 1:n_samples

    trigger = zeros(Int, n_samples)
    trigger_info = fill("", n_samples)

    if !isnothing(dat.triggers) && length(dat.triggers.idx) > 0
        for i = 1:length(dat.triggers.idx)
            sample_idx = dat.triggers.idx[i]
            if 1 <= sample_idx <= n_samples
                trigger[sample_idx] = dat.triggers.val[i]
                trigger_info[sample_idx] = string(dat.triggers.val[i])
            end
        end
    end

    # Filter out stimulus channels (kind == 3 or named Status/STI)
    data_ch_indices = Int[]
    scales = Float64[]
    channel_labels = Symbol[]

    for (i, ch) in enumerate(dat.header.ch_info)
        if !_is_stim_channel(ch)
            push!(data_ch_indices, i)
            push!(channel_labels, Symbol(ch.ch_name))
            if ch.unit == 107 # Volts (with FIFF unit_mul) -> Microvolts
                push!(scales, 10.0^(6 + Int(ch.unit_mul)))
            else
                push!(scales, 1.0)
            end
        end
    end

    scaled_data = Float64.(dat.data[data_ch_indices, :]') .* scales'

    df = hcat(
        DataFrame(time = time, sample = sample, trigger = trigger, trigger_info = trigger_info),
        DataFrame(scaled_data, channel_labels),
    )
    return df
end

function create_eegfun_data(dat::FunctionalImageFormat.FifData, layout::Layout)::ContinuousData
    @info "Creating EEG DataFrame (*.fif)"
    file_name = basename_without_ext(dat.filename)
    df = _create_eegfun_dataframe(dat)
    sr = round(Int, dat.header.sfreq)
    return ContinuousData(file_name, df, copy(layout), sr, AnalysisInfo(sample_rate = Float64(dat.header.sfreq)))
end

# Optional Layout variant for quick databrowser visualization
function create_eegfun_data(dat::FunctionalImageFormat.FifData)
    channel_labels = Symbol.([ch.ch_name for ch in dat.header.ch_info if !_is_stim_channel(ch)])
    layout = _create_layout_from_labels(channel_labels)
    return create_eegfun_data(dat, layout)
end

function _xdf_srate(eeg_stream)::Float64
    return eeg_stream.header.nominal_srate > 0 ? Float64(eeg_stream.header.nominal_srate) :
           (eeg_stream.header.effective_srate > 0 ? Float64(eeg_stream.header.effective_srate) : 1.0)
end

function _xdf_channel_metadata(eeg_stream)
    n_channels = eeg_stream.header.channel_count
    channel_labels = String[]
    channel_scales = Float64[]

    info_root =
        haskey(eeg_stream.header.info, "info") && eeg_stream.header.info["info"] isa AbstractDict ? eeg_stream.header.info["info"] :
        eeg_stream.header.info

    if haskey(info_root, "desc") &&
       info_root["desc"] isa AbstractDict &&
       haskey(info_root["desc"], "channels") &&
       info_root["desc"]["channels"] isa AbstractDict &&
       haskey(info_root["desc"]["channels"], "channel")
        channels_info = info_root["desc"]["channels"]["channel"]
        if channels_info isa AbstractVector
            for ch in channels_info
                if ch isa AbstractDict
                    push!(channel_labels, string(get(ch, "label", "Unknown")))
                    push!(channel_scales, _unit_to_microvolt_scale(string(get(ch, "unit", ""))))
                end
            end
        elseif channels_info isa AbstractDict
            push!(channel_labels, string(get(channels_info, "label", "Unknown")))
            push!(channel_scales, _unit_to_microvolt_scale(string(get(channels_info, "unit", ""))))
        end
    end

    if length(channel_labels) != n_channels
        channel_labels = ["ch$i" for i = 1:n_channels]
        channel_scales = ones(Float64, n_channels)
    end

    return Symbol.(channel_labels[1:n_channels]), channel_scales[1:n_channels]
end

"""Create a time+sample+trigger+channel DataFrame from raw XdfData."""
function _create_eegfun_dataframe(dat::ExtensibleDataFormat.XdfData)::DataFrame
    eeg_streams = [s for s in values(dat.streams) if s.header.type == "EEG"]
    if isempty(eeg_streams)
        @minimal_error "No EEG stream found in XDF data"
    end
    eeg_stream = eeg_streams[1]

    n_samples = size(eeg_stream.time_series, 1)
    srate = _xdf_srate(eeg_stream)

    # Enforce uniform time array for ContinuousData (required for DSP filters)
    time = collect((0:(n_samples-1)) ./ srate)
    sample = 1:n_samples

    channel_labels, channel_scales = _xdf_channel_metadata(eeg_stream)

    trigger = zeros(Int, n_samples)
    trigger_info = fill("", n_samples)

    marker_streams = [s for s in values(dat.streams) if s.header.type == "Markers"]
    if !isempty(marker_streams) && !isempty(eeg_stream.timestamps)
        marker_stream = marker_streams[1]
        for (i, t) in enumerate(marker_stream.timestamps)
            idx = find_closest_time_index(eeg_stream.timestamps, t)
            marker_val = marker_stream.time_series[i, 1]
            int_val = _parse_marker_code(string(marker_val))
            trigger[idx] = int_val != 0 ? int_val : 1
            trigger_info[idx] = string(marker_val)
        end
    end

    scaled_data = Float64.(eeg_stream.time_series) .* channel_scales'

    df = hcat(
        DataFrame(time = time, sample = sample, trigger = trigger, trigger_info = trigger_info),
        DataFrame(scaled_data, channel_labels),
    )
    return df
end

function create_eegfun_data(dat::ExtensibleDataFormat.XdfData, layout::Layout)::ContinuousData
    @info "Creating EEG DataFrame (*.xdf)"
    file_name = basename_without_ext(dat.filename)
    df = _create_eegfun_dataframe(dat)

    eeg_streams = [s for s in values(dat.streams) if s.header.type == "EEG"]
    sample_rate = !isempty(eeg_streams) ? _xdf_srate(eeg_streams[1]) : 1.0
    sr = round(Int, sample_rate)

    return ContinuousData(file_name, df, copy(layout), sr, AnalysisInfo(sample_rate = sample_rate))
end

function create_eegfun_data(dat::ExtensibleDataFormat.XdfData)
    eeg_streams = [s for s in values(dat.streams) if s.header.type == "EEG"]
    if isempty(eeg_streams)
        @minimal_error "No EEG stream found in XDF data"
    end
    channel_labels, _ = _xdf_channel_metadata(eeg_streams[1])
    layout = _create_layout_from_labels(channel_labels)
    return create_eegfun_data(dat, layout)
end

function _create_eegfun_epoch_dataframes(dat::FunctionalImageFormat.FifEpochs)::Vector{DataFrame}
    _, n_samples, n_epochs = size(dat.data)
    sample_rate = Float64(dat.header.sfreq)

    # Filter out stimulus channels
    data_ch_indices = Int[]
    scales = Float64[]
    channel_labels = Symbol[]

    stim_idx = nothing
    stim_info = nothing

    for (i, ch) in enumerate(dat.header.ch_info)
        if _is_stim_channel(ch) && stim_idx === nothing
            stim_idx = i
            stim_info = ch
        end
        if !_is_stim_channel(ch)
            push!(data_ch_indices, i)
            push!(channel_labels, Symbol(ch.ch_name))
            if ch.unit == 107 # Volts (with FIFF unit_mul) -> Microvolts
                push!(scales, 10.0^(6 + Int(ch.unit_mul)))
            else
                push!(scales, 1.0)
            end
        end
    end

    dfs = Vector{DataFrame}(undef, n_epochs)

    for e = 1:n_epochs
        first_samp = 0
        if !isempty(dat.first_sample)
            first_samp = length(dat.first_sample) >= e ? dat.first_sample[e] : dat.first_sample[1]
        end
        time = collect(first_samp:(first_samp+n_samples-1)) ./ sample_rate
        sample = 1:n_samples
        epoch_col = fill(e, n_samples)

        trigger = zeros(Int, n_samples)
        trigger_info = fill("", n_samples)

        if stim_idx !== nothing
            raw_float = @views dat.data[stim_idx, :, e] ./ (stim_info.range * stim_info.cal)
            raw = Int32.(round.(raw_float))

            if raw[1] > 0
                trigger[1] = raw[1]
                trigger_info[1] = string(raw[1])
            end
            for i = 2:n_samples
                if raw[i] != raw[i-1] && raw[i] > 0
                    trigger[i] = raw[i]
                    trigger_info[i] = string(raw[i])
                end
            end
        end

        scaled_data = Float64.(dat.data[data_ch_indices, :, e]') .* scales'

        dfs[e] = hcat(
            DataFrame(time = time, sample = sample, epoch = epoch_col, trigger = trigger, trigger_info = trigger_info),
            DataFrame(scaled_data, channel_labels),
        )
    end

    return dfs
end

function create_eegfun_data(dat::FunctionalImageFormat.FifEpochs, layout::Layout)::Union{EpochData,ErpData,Vector{ErpData}}
    file_name = basename_without_ext(dat.filename)
    dfs = _create_eegfun_epoch_dataframes(dat)
    sr = round(Int, dat.header.sfreq)
    sfreq = Float64(dat.header.sfreq)

    condition_name = isempty(dat.comments) ? "FIF Epochs" : dat.comments[1]

    if !isempty(dat.nave) && length(dat.nave) == length(dfs) && all(>(1), dat.nave)
        @info "Creating EEG DataFrame (*.fif evoked)"
        for df in dfs
            hasproperty(df, :epoch) && select!(df, Not(:epoch))
        end
        if length(dfs) == 1
            return ErpData(file_name, 1, condition_name, dfs[1], copy(layout), sr, AnalysisInfo(sample_rate = sfreq), Int(dat.nave[1]))
        else
            return ErpData[
                ErpData(
                    file_name,
                    i,
                    i <= length(dat.comments) && !isempty(dat.comments[i]) ? dat.comments[i] : "Condition $i",
                    dfs[i],
                    copy(layout),
                    sr,
                    AnalysisInfo(sample_rate = sfreq),
                    Int(dat.nave[i]),
                ) for i in eachindex(dfs)
            ]
        end
    else
        @info "Creating EEG DataFrame (*.fif epochs)"
        return EpochData(file_name, 1, condition_name, dfs, copy(layout), sr, AnalysisInfo(sample_rate = sfreq))
    end
end

# Optional Layout variant for quick databrowser visualization
function create_eegfun_data(dat::FunctionalImageFormat.FifEpochs)
    channel_labels = Symbol.([ch.ch_name for ch in dat.header.ch_info if !_is_stim_channel(ch)])
    layout = _create_layout_from_labels(channel_labels)
    return create_eegfun_data(dat, layout)
end
