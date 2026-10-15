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
    return ContinuousData(file_name, df, layout, dat.header.sample_rate[1], AnalysisInfo(sample_rate = dat.header.sample_rate[1]))
end

# Optional Layout variant for quick databrowser visualization
function create_eegfun_data(dat::BiosemiDataFormat.BiosemiData)
    channel_labels = Symbol.(dat.header.channel_labels[1:(end-1)])
    layout = _create_layout_from_labels(channel_labels)
    return create_eegfun_data(dat, layout)
end

"""Create a time+sample+trigger+channel DataFrame from raw EdfData."""
function _create_eegfun_dataframe(dat::EuropeanDataFormat.EdfData)::DataFrame

    n_samples = size(dat.data, 1)
    n_channels = size(dat.data, 2)
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

    df = hcat(
        DataFrame(time = time, sample = sample, trigger = trigger, trigger_info = trigger_info),
        DataFrame(Float64.(dat.data), Symbol.(dat.header.channel_labels[1:size(dat.data, 2)])),
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
    return ContinuousData(file_name, df, layout, dat.header.sample_rate[1], AnalysisInfo(sample_rate = dat.header.sample_rate[1]))
end

# Optional Layout variant for quick databrowser visualization
function create_eegfun_data(dat::EuropeanDataFormat.EdfData)
    channel_labels = Symbol.(dat.header.channel_labels[1:size(dat.data, 2)])
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
    return ContinuousData(file_name, df, layout, dat.header.Fs, AnalysisInfo(sample_rate = dat.header.Fs))
end

# Optional Layout variant for quick databrowser visualization
function create_eegfun_data(dat::BrainVisionDataFormat.BrainVisionData)
    channel_labels = Symbol.(dat.header.label[1:size(dat.data, 2)])
    layout = _create_layout_from_labels(channel_labels)
    return create_eegfun_data(dat, layout)
end

"""
    _extract_triggers_from_markers(markers::Vector{BrainVisionDataFormat.BrainVisionMarker}, n_samples::Int)::Tuple{Vector{Int}, Vector{String}}

Extract trigger values from BrainVision markers and create trigger and marker string vectors.

# Arguments
- `markers::Vector{BrainVisionDataFormat.BrainVisionMarker}`: Vector of BrainVision markers
- `n_samples::Int`: Number of samples in the data

# Returns
- `Tuple{Vector{Int}, Vector{String}}`: (trigger vector, marker string vector) with values at appropriate sample positions
"""
function _extract_triggers_from_markers(
    markers::Vector{BrainVisionDataFormat.BrainVisionMarker},
    n_samples::Int,
)::Tuple{Vector{Int},Vector{String}}
    trigger = zeros(Int, n_samples)
    trigger_info = fill("", n_samples)

    @debug "Processing $(length(markers)) markers for $n_samples samples"

    # Debug: Show first few marker sample indices
    if length(markers) > 0
        first_few_samples = [marker.sample for marker in markers[1:min(5, length(markers))]]
        @debug "First few marker sample indices: $first_few_samples"
    end

    # First pass: extract all unique trigger values (including empty strings for system markers)
    # Check if markers are 0-based or 1-based by looking at the first marker
    is_zero_based = length(markers) > 0 && markers[1].sample == 0
    @debug "Detected $(is_zero_based ? "0-based" : "1-based") indexing for marker samples"

    unique_values = Set{String}()
    valid_markers = 0
    for marker in markers
        # Convert to 1-based if needed
        sample_idx = is_zero_based ? marker.sample + 1 : marker.sample

        if 1 <= sample_idx <= n_samples
            # Include all markers, even those with empty values (like "New Segment")
            push!(unique_values, marker.value)
            valid_markers += 1
        else
            @minimal_warning "Marker sample $sample_idx out of bounds (1:$n_samples), skipping"
        end
    end

    @debug "Found $valid_markers valid markers with $(length(unique_values)) unique values: $(collect(unique_values))"

    # Create mapping from original values to sequential integers (1, 2, 3, ...)
    value_to_trigger = Dict{String,Int}()
    for (i, value) in enumerate(sort(collect(unique_values)))
        value_to_trigger[value] = i
    end

    # Second pass: assign sequential trigger values and original marker strings
    for marker in markers
        # Convert to 1-based if needed
        sample_idx = is_zero_based ? marker.sample + 1 : marker.sample

        if 1 <= sample_idx <= n_samples
            # Include all markers, even those with empty values
            trigger[sample_idx] = value_to_trigger[marker.value]
            trigger_info[sample_idx] = marker.value
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
            if ch.unit == 107 # Volts -> Microvolts
                push!(scales, 1e6)
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
    return ContinuousData(file_name, df, layout, Int(round(dat.header.sfreq)), AnalysisInfo(sample_rate = Int(round(dat.header.sfreq))))
end

# Optional Layout variant for quick databrowser visualization
function create_eegfun_data(dat::FunctionalImageFormat.FifData)
    channel_labels = Symbol.([ch.ch_name for ch in dat.header.ch_info if !_is_stim_channel(ch)])
    layout = _create_layout_from_labels(channel_labels)
    return create_eegfun_data(dat, layout)
end

"""Create a time+sample+trigger+channel DataFrame from raw XdfData."""
function _create_eegfun_dataframe(dat::ExtensibleDataFormat.XdfData)::DataFrame
    eeg_streams = [s for s in values(dat.streams) if s.header.type == "EEG"]
    if isempty(eeg_streams)
        @minimal_error "No EEG stream found in XDF data"
    end
    eeg_stream = eeg_streams[1]

    n_samples = size(eeg_stream.time_series, 1)
    n_channels = eeg_stream.header.channel_count

    srate =
        eeg_stream.header.nominal_srate > 0 ? eeg_stream.header.nominal_srate :
        (eeg_stream.header.effective_srate > 0 ? eeg_stream.header.effective_srate : 1.0)

    # Enforce uniform time array for ContinuousData (required for DSP filters)
    time = collect((0:(n_samples-1)) ./ srate)
    sample = 1:n_samples

    channel_labels = String[]
    if haskey(eeg_stream.header.info, "desc") &&
       haskey(eeg_stream.header.info["desc"], "channels") &&
       haskey(eeg_stream.header.info["desc"]["channels"], "channel")
        channels_info = eeg_stream.header.info["desc"]["channels"]["channel"]
        if channels_info isa Vector
            for ch in channels_info
                push!(channel_labels, get(ch, "label", "Unknown"))
            end
        elseif channels_info isa Dict
            push!(channel_labels, get(channels_info, "label", "Unknown"))
        end
    end

    if length(channel_labels) != n_channels
        channel_labels = ["ch$i" for i = 1:n_channels]
    end

    trigger = zeros(Int, n_samples)
    trigger_info = fill("", n_samples)

    marker_streams = [s for s in values(dat.streams) if s.header.type == "Markers"]
    if !isempty(marker_streams)
        marker_stream = marker_streams[1]
        for (i, t) in enumerate(marker_stream.timestamps)
            idx = find_closest_time_index(eeg_stream.timestamps, t)
            marker_val = marker_stream.time_series[i, 1]
            int_val = tryparse(Int, string(marker_val))
            trigger[idx] = int_val !== nothing ? int_val : 1
            trigger_info[idx] = string(marker_val)
        end
    end

    df = hcat(
        DataFrame(time = time, sample = sample, trigger = trigger, trigger_info = trigger_info),
        DataFrame(Float64.(eeg_stream.time_series), Symbol.(channel_labels[1:n_channels])),
    )
    return df
end

function create_eegfun_data(dat::ExtensibleDataFormat.XdfData, layout::Layout)::ContinuousData
    @info "Creating EEG DataFrame (*.xdf)"
    file_name = basename_without_ext(dat.filename)
    df = _create_eegfun_dataframe(dat)

    eeg_streams = [s for s in values(dat.streams) if s.header.type == "EEG"]
    sample_rate = !isempty(eeg_streams) ? eeg_streams[1].header.nominal_srate : 1.0

    return ContinuousData(file_name, df, layout, Int(round(sample_rate)), AnalysisInfo(sample_rate = Int(round(sample_rate))))
end

function create_eegfun_data(dat::ExtensibleDataFormat.XdfData)
    eeg_streams = [s for s in values(dat.streams) if s.header.type == "EEG"]
    if isempty(eeg_streams)
        @minimal_error "No EEG stream found in XDF data"
    end
    eeg_stream = eeg_streams[1]
    n_channels = eeg_stream.header.channel_count

    channel_labels = String[]
    if haskey(eeg_stream.header.info, "desc") &&
       haskey(eeg_stream.header.info["desc"], "channels") &&
       haskey(eeg_stream.header.info["desc"]["channels"], "channel")
        channels_info = eeg_stream.header.info["desc"]["channels"]["channel"]
        if channels_info isa Vector
            for ch in channels_info
                push!(channel_labels, get(ch, "label", "Unknown"))
            end
        elseif channels_info isa Dict
            push!(channel_labels, get(channels_info, "label", "Unknown"))
        end
    end
    if length(channel_labels) != n_channels
        channel_labels = ["ch$i" for i = 1:n_channels]
    end

    layout = _create_layout_from_labels(Symbol.(channel_labels))
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
            if ch.unit == 107 # Volts -> Microvolts
                push!(scales, 1e6)
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

function create_eegfun_data(dat::FunctionalImageFormat.FifEpochs, layout::Layout)::Union{EpochData,ErpData}
    file_name = basename_without_ext(dat.filename)
    dfs = _create_eegfun_epoch_dataframes(dat)

    condition_name = isempty(dat.comments) ? "FIF Epochs" : dat.comments[1]

    if length(dfs) == 1 && !isempty(dat.nave) && dat.nave[1] > 1
        @info "Creating EEG DataFrame (*.fif evoked)"
        return ErpData(
            file_name,
            1,
            condition_name,
            dfs[1],
            layout,
            Int(round(dat.header.sfreq)),
            AnalysisInfo(sample_rate = Int(round(dat.header.sfreq))),
            Int(dat.nave[1]),
        )
    else
        @info "Creating EEG DataFrame (*.fif epochs)"
        return EpochData(
            file_name,
            1,
            condition_name,
            dfs,
            layout,
            Int(round(dat.header.sfreq)),
            AnalysisInfo(sample_rate = Int(round(dat.header.sfreq))),
        )
    end
end

# Optional Layout variant for quick databrowser visualization
function create_eegfun_data(dat::FunctionalImageFormat.FifEpochs)
    channel_labels = Symbol.([ch.ch_name for ch in dat.header.ch_info if !_is_stim_channel(ch)])
    layout = _create_layout_from_labels(channel_labels)
    return create_eegfun_data(dat, layout)
end
