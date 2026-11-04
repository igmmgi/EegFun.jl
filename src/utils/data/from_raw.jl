function _create_layout_from_labels(labels::Vector{Symbol})::Layout
    n = length(labels)
    return Layout(DataFrame(label = labels, inc = zeros(n), azi = zeros(n)), nothing, nothing, nothing)
end

function _build_samples_x_channels_df(
    meta_cols::Vector{Any},
    meta_names::Vector{Symbol},
    raw_mat::AbstractMatrix,
    col_indices::AbstractVector{Int},
    ch_labels::Vector{Symbol},
    scales::Union{Nothing,Vector{Float64}} = nothing,
)::DataFrame
    n_samples = size(raw_mat, 1)
    n_meta = length(meta_cols)
    n_ch = length(col_indices)
    cols = Vector{Any}(undef, n_meta + n_ch)
    names_out = Vector{Symbol}(undef, n_meta + n_ch)
    @inbounds for i = 1:n_meta
        cols[i] = meta_cols[i]
        names_out[i] = meta_names[i]
    end
    for k = 1:n_ch
        c_idx = col_indices[k]
        dst = Vector{Float64}(undef, n_samples)
        if isnothing(scales) || scales[k] == 1.0
            @inbounds @simd for t = 1:n_samples
                dst[t] = Float64(raw_mat[t, c_idx])
            end
        else
            sc = scales[k]
            @inbounds @simd for t = 1:n_samples
                dst[t] = Float64(raw_mat[t, c_idx]) * sc
            end
        end
        cols[n_meta+k] = dst
        names_out[n_meta+k] = ch_labels[k]
    end
    return DataFrame(cols, names_out; copycols = false)
end

function _build_channels_x_samples_df(
    meta_cols::Vector{Any},
    meta_names::Vector{Symbol},
    raw_mat::AbstractMatrix,
    ch_indices::Vector{Int},
    ch_labels::Vector{Symbol},
    scales::Vector{Float64},
)::DataFrame
    n_samples = size(raw_mat, 2)
    n_meta = length(meta_cols)
    n_ch = length(ch_indices)
    cols = Vector{Any}(undef, n_meta + n_ch)
    names_out = Vector{Symbol}(undef, n_meta + n_ch)
    @inbounds for i = 1:n_meta
        cols[i] = meta_cols[i]
        names_out[i] = meta_names[i]
    end
    ch_vecs = Vector{Vector{Float64}}(undef, n_ch)
    for k = 1:n_ch
        v = Vector{Float64}(undef, n_samples)
        ch_vecs[k] = v
        cols[n_meta+k] = v
        names_out[n_meta+k] = ch_labels[k]
    end

    tb = 64
    for t0 = 1:tb:n_samples
        t1 = min(t0 + tb - 1, n_samples)
        for k = 1:n_ch
            c_idx = ch_indices[k]
            sc = scales[k]
            dst = ch_vecs[k]
            @inbounds @simd for t = t0:t1
                dst[t] = Float64(raw_mat[c_idx, t]) * sc
            end
        end
    end
    return DataFrame(cols, names_out; copycols = false)
end

function _create_eegfun_dataframe(dat::BiosemiDataFormat.BiosemiData)::DataFrame
    n_samples = length(dat.time)
    time_col = dat.time isa Vector{Float64} ? copy(dat.time) : collect(Float64, dat.time)
    ch_labels = Symbol.(dat.header.channel_labels[1:(end-1)])
    return _build_samples_x_channels_df(
        Any[time_col, collect(1:n_samples), _clean_triggers(dat.triggers.raw)],
        [:time, :sample, :trigger],
        dat.data,
        1:length(ch_labels),
        ch_labels,
    )
end

"""
    create_eegfun_data(dat, [layout::Layout])

Convert raw imported EEG data (`BiosemiData`, `EdfData`, `BrainVisionData`, `FifData`, `FifEpochs`, or `XdfData`)
into EegFun `ContinuousData`, `EpochData`, or `ErpData` (scaling V/mV channels to μV).
"""
function create_eegfun_data(dat::BiosemiDataFormat.BiosemiData, layout::Layout)::ContinuousData
    @info "Creating EEG DataFrame from Biosemi data"
    sr = Int(dat.header.sample_rate[1])
    return ContinuousData(filename(dat), _create_eegfun_dataframe(dat), copy(layout), sr, AnalysisInfo(sample_rate = Float64(sr)))
end

create_eegfun_data(dat::BiosemiDataFormat.BiosemiData) =
    create_eegfun_data(dat, _create_layout_from_labels(Symbol.(dat.header.channel_labels[1:(end-1)])))

_is_edf_status_label(lbl::AbstractString) = lowercase(strip(lbl)) in ("status", "edf annotations", "tal")

function _unit_to_microvolt_scale(unit_str::AbstractString)::Float64
    u = strip(unit_str)
    isempty(u) && return 1.0
    ul = lowercase(u)
    (u == "V" || ul in ("v", "volt", "volts")) && return 1e6
    (u == "mV" || ul in ("mv", "millivolt", "millivolts")) && return 1e3
    return 1.0
end

function _edf_data_channels(dat::EuropeanDataFormat.EdfData)
    hdr_labels = dat.header.channel_labels
    hdr_units = dat.header.channel_unit
    n_data_cols = size(dat.data, 2)
    non_status = [i for (i, lbl) in enumerate(hdr_labels) if !_is_edf_status_label(lbl)]

    if n_data_cols == length(hdr_labels)
        sel_hdr = non_status
        data_cols = non_status
    elseif length(non_status) == n_data_cols
        sel_hdr = non_status
        data_cols = collect(1:n_data_cols)
    else
        sel_hdr = length(non_status) >= n_data_cols ? non_status[1:n_data_cols] : collect(1:n_data_cols)
        data_cols = collect(1:n_data_cols)
    end
    labels = Symbol[Symbol(strip(hdr_labels[i])) for i in sel_hdr]
    scales = Float64[i <= length(hdr_units) ? _unit_to_microvolt_scale(hdr_units[i]) : 1.0 for i in sel_hdr]
    return data_cols, labels, scales
end

function _fill_idx_triggers!(trigger::Vector{Int}, trigger_info::Vector{String}, trigs, n_samples::Int)
    if !isnothing(trigs) && !isempty(trigs.idx)
        for i in eachindex(trigs.idx)
            s_idx = trigs.idx[i]
            if 1 <= s_idx <= n_samples
                v = trigs.val[i]
                trigger[s_idx] = v
                trigger_info[s_idx] = string(v)
            end
        end
    end
end

function _create_eegfun_dataframe(dat::EuropeanDataFormat.EdfData)::DataFrame
    n_samples = size(dat.data, 1)
    inv_sr = 1.0 / dat.header.sample_rate[1]
    time = Float64[(i - 1) * inv_sr for i = 1:n_samples]
    trigger = zeros(Int, n_samples)
    trigger_info = fill("", n_samples)
    _fill_idx_triggers!(trigger, trigger_info, dat.triggers, n_samples)

    data_cols, ch_labels, scales = _edf_data_channels(dat)
    return _build_samples_x_channels_df(
        Any[time, collect(1:n_samples), trigger, trigger_info],
        [:time, :sample, :trigger, :trigger_info],
        dat.data,
        data_cols,
        ch_labels,
        scales,
    )
end

"""
    download_eegfun_datasets() -> String

Download and cache built-in tutorial datasets via `DataDeps.jl`, returning the local directory path.
"""
download_eegfun_datasets() = datadep"TutorialDataSets"

function create_eegfun_data(dat::EuropeanDataFormat.EdfData, layout::Layout)::ContinuousData
    @info "Creating EEG DataFrame (*.edf)"
    sr = Int(dat.header.sample_rate[1])
    return ContinuousData(
        basename_without_ext(dat.filename),
        _create_eegfun_dataframe(dat),
        copy(layout),
        sr,
        AnalysisInfo(sample_rate = Float64(sr)),
    )
end

function create_eegfun_data(dat::EuropeanDataFormat.EdfData)
    _, ch_labels, _ = _edf_data_channels(dat)
    return create_eegfun_data(dat, _create_layout_from_labels(ch_labels))
end

function _create_eegfun_dataframe(dat::BrainVisionDataFormat.BrainVisionData)::DataFrame
    isnothing(dat.data) && @minimal_error "BrainVision data is empty (data field is nothing)"
    isnothing(dat.header) && @minimal_error "BrainVision header is empty (header field is nothing)"

    n_samples, n_ch = size(dat.data)
    inv_sr = 1.0 / dat.header.Fs
    time = Float64[(i - 1) * inv_sr for i = 1:n_samples]
    ch_labels = Symbol.(dat.header.label[1:n_ch])
    length(ch_labels) == n_ch || @minimal_error "Number of channel labels ($(length(ch_labels))) does not match number of channels ($n_ch)"

    if isnothing(dat.markers) || isempty(dat.markers)
        trigger = zeros(Int, n_samples)
        trigger_info = fill("", n_samples)
    else
        trigger, trigger_info = _extract_triggers_from_markers(dat.markers, n_samples)
    end

    return _build_samples_x_channels_df(
        Any[time, collect(1:n_samples), trigger, trigger_info],
        [:time, :sample, :trigger, :trigger_info],
        dat.data,
        1:n_ch,
        ch_labels,
    )
end

function create_eegfun_data(dat::BrainVisionDataFormat.BrainVisionData, layout::Layout)::ContinuousData
    @info "Creating EEG DataFrame (*.vhdr/*.vmrk/*.eeg)"
    sr = round(Int, dat.header.Fs)
    return ContinuousData(
        basename_without_ext(dat.filename),
        _create_eegfun_dataframe(dat),
        copy(layout),
        sr,
        AnalysisInfo(sample_rate = Float64(dat.header.Fs)),
    )
end

create_eegfun_data(dat::BrainVisionDataFormat.BrainVisionData) =
    create_eegfun_data(dat, _create_layout_from_labels(Symbol.(dat.header.label[1:size(dat.data, 2)])))

function _parse_marker_code(val::AbstractString)::Int
    s = strip(val)
    isempty(s) && return 0
    direct = tryparse(Int, s)
    !isnothing(direct) && return direct
    m = match(r"-?\d+", s)
    return !isnothing(m) ? something(tryparse(Int, m.match), 0) : 0
end

function _extract_triggers_from_markers(
    markers::Vector{BrainVisionDataFormat.BrainVisionMarker},
    n_samples::Int,
)::Tuple{Vector{Int},Vector{String}}
    trigger = zeros(Int, n_samples)
    trigger_info = fill("", n_samples)

    for marker in markers
        s_idx = marker.sample
        s_idx == 0 && continue
        if 1 <= s_idx <= n_samples
            code = _parse_marker_code(marker.value)
            info_str = isempty(strip(marker.value)) ? strip(marker.type) : marker.value
            (code != 0 || trigger[s_idx] == 0) && (trigger[s_idx] = code)
            (code != 0 || isempty(trigger_info[s_idx])) && (trigger_info[s_idx] = info_str)
        else
            @minimal_warning "Marker sample $s_idx out of bounds (1:$n_samples), skipping"
        end
    end

    return trigger, trigger_info
end

_is_stim_channel(ch) = ch.kind == 3 || startswith(uppercase(ch.ch_name), "STI") || uppercase(ch.ch_name) == "STATUS"

function _fif_channel_selection(ch_info)
    ch_indices = Int[]
    ch_labels = Symbol[]
    scales = Float64[]
    stim_idx = nothing
    stim_info = nothing

    for (i, ch) in enumerate(ch_info)
        if _is_stim_channel(ch)
            if isnothing(stim_idx)
                stim_idx = i
                stim_info = ch
            end
        else
            push!(ch_indices, i)
            push!(ch_labels, Symbol(ch.ch_name))
            push!(scales, ch.unit == 107 ? 10.0^(6 + Int(ch.unit_mul)) : 1.0)
        end
    end
    return ch_indices, ch_labels, scales, stim_idx, stim_info
end

function _create_eegfun_dataframe(dat::FunctionalImageFormat.FifData)::DataFrame
    _, n_samples = size(dat.data)
    inv_sr = 1.0 / Float64(dat.header.sfreq)
    time = Float64[(i - 1) * inv_sr for i = 1:n_samples]
    trigger = zeros(Int, n_samples)
    trigger_info = fill("", n_samples)
    _fill_idx_triggers!(trigger, trigger_info, dat.triggers, n_samples)

    ch_indices, ch_labels, scales, _, _ = _fif_channel_selection(dat.header.ch_info)
    return _build_channels_x_samples_df(
        Any[time, collect(1:n_samples), trigger, trigger_info],
        [:time, :sample, :trigger, :trigger_info],
        dat.data,
        ch_indices,
        ch_labels,
        scales,
    )
end

function create_eegfun_data(dat::FunctionalImageFormat.FifData, layout::Layout)::ContinuousData
    @info "Creating EEG DataFrame (*.fif)"
    sr = round(Int, dat.header.sfreq)
    return ContinuousData(
        basename_without_ext(dat.filename),
        _create_eegfun_dataframe(dat),
        copy(layout),
        sr,
        AnalysisInfo(sample_rate = Float64(dat.header.sfreq)),
    )
end

create_eegfun_data(dat::FunctionalImageFormat.FifData) =
    create_eegfun_data(dat, _create_layout_from_labels(Symbol[Symbol(ch.ch_name) for ch in dat.header.ch_info if !_is_stim_channel(ch)]))

_xdf_srate(eeg_stream)::Float64 =
    eeg_stream.header.nominal_srate > 0 ? Float64(eeg_stream.header.nominal_srate) :
    (eeg_stream.header.effective_srate > 0 ? Float64(eeg_stream.header.effective_srate) : 1.0)

function _xdf_channel_metadata(eeg_stream)
    n_channels = eeg_stream.header.channel_count
    ch_labels = String[]
    ch_scales = Float64[]

    info_root =
        haskey(eeg_stream.header.info, "info") && eeg_stream.header.info["info"] isa AbstractDict ? eeg_stream.header.info["info"] :
        eeg_stream.header.info

    desc = get(info_root, "desc", nothing)
    if desc isa AbstractDict
        ch_block = get(desc, "channels", nothing)
        if ch_block isa AbstractDict && haskey(ch_block, "channel")
            channels_info = ch_block["channel"]
            if channels_info isa AbstractVector
                for ch in channels_info
                    if ch isa AbstractDict
                        push!(ch_labels, string(get(ch, "label", "Unknown")))
                        push!(ch_scales, _unit_to_microvolt_scale(string(get(ch, "unit", ""))))
                    end
                end
            elseif channels_info isa AbstractDict
                push!(ch_labels, string(get(channels_info, "label", "Unknown")))
                push!(ch_scales, _unit_to_microvolt_scale(string(get(channels_info, "unit", ""))))
            end
        end
    end

    if length(ch_labels) != n_channels
        ch_labels = ["ch$i" for i = 1:n_channels]
        ch_scales = ones(Float64, n_channels)
    end
    return Symbol.(ch_labels[1:n_channels]), ch_scales[1:n_channels]
end

function _create_eegfun_dataframe(dat::ExtensibleDataFormat.XdfData)::DataFrame
    eeg_streams = [s for s in values(dat.streams) if s.header.type == "EEG"]
    isempty(eeg_streams) && @minimal_error "No EEG stream found in XDF data"
    eeg_stream = eeg_streams[1]

    n_samples = size(eeg_stream.time_series, 1)
    inv_sr = 1.0 / _xdf_srate(eeg_stream)
    time = Float64[(i - 1) * inv_sr for i = 1:n_samples]
    ch_labels, ch_scales = _xdf_channel_metadata(eeg_stream)

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

    return _build_samples_x_channels_df(
        Any[time, collect(1:n_samples), trigger, trigger_info],
        [:time, :sample, :trigger, :trigger_info],
        eeg_stream.time_series,
        1:length(ch_labels),
        ch_labels,
        ch_scales,
    )
end

function create_eegfun_data(dat::ExtensibleDataFormat.XdfData, layout::Layout)::ContinuousData
    @info "Creating EEG DataFrame (*.xdf)"
    df = _create_eegfun_dataframe(dat)
    eeg_streams = [s for s in values(dat.streams) if s.header.type == "EEG"]
    srate = !isempty(eeg_streams) ? _xdf_srate(eeg_streams[1]) : 1.0
    return ContinuousData(basename_without_ext(dat.filename), df, copy(layout), round(Int, srate), AnalysisInfo(sample_rate = srate))
end

function create_eegfun_data(dat::ExtensibleDataFormat.XdfData)
    eeg_streams = [s for s in values(dat.streams) if s.header.type == "EEG"]
    isempty(eeg_streams) && @minimal_error "No EEG stream found in XDF data"
    ch_labels, _ = _xdf_channel_metadata(eeg_streams[1])
    return create_eegfun_data(dat, _create_layout_from_labels(ch_labels))
end

function _create_eegfun_epoch_dataframes(dat::FunctionalImageFormat.FifEpochs)::Vector{DataFrame}
    _, n_samples, n_epochs = size(dat.data)
    inv_sr = 1.0 / Float64(dat.header.sfreq)
    ch_indices, ch_labels, scales, stim_idx, stim_info = _fif_channel_selection(dat.header.ch_info)
    dfs = Vector{DataFrame}(undef, n_epochs)

    for e = 1:n_epochs
        first_samp = !isempty(dat.first_sample) ? (length(dat.first_sample) >= e ? dat.first_sample[e] : dat.first_sample[1]) : 0
        time = Float64[(first_samp + i - 1) * inv_sr for i = 1:n_samples]
        trigger = zeros(Int, n_samples)
        trigger_info = fill("", n_samples)

        if !isnothing(stim_idx) && n_samples > 0
            inv_stim = 1.0 / (stim_info.range * stim_info.cal)
            prev_raw = round(Int32, dat.data[stim_idx, 1, e] * inv_stim)
            if prev_raw > 0
                trigger[1] = prev_raw
                trigger_info[1] = string(prev_raw)
            end
            @inbounds for i = 2:n_samples
                cur = round(Int32, dat.data[stim_idx, i, e] * inv_stim)
                if cur != prev_raw && cur > 0
                    trigger[i] = cur
                    trigger_info[i] = string(cur)
                end
                prev_raw = cur
            end
        end

        dfs[e] = _build_channels_x_samples_df(
            Any[time, collect(1:n_samples), fill(e, n_samples), trigger, trigger_info],
            [:time, :sample, :epoch, :trigger, :trigger_info],
            @view(dat.data[:, :, e]),
            ch_indices,
            ch_labels,
            scales,
        )
    end
    return dfs
end

function create_eegfun_data(dat::FunctionalImageFormat.FifEpochs, layout::Layout)::Union{EpochData,ErpData,Vector{ErpData}}
    fname = basename_without_ext(dat.filename)
    dfs = _create_eegfun_epoch_dataframes(dat)
    sr = round(Int, dat.header.sfreq)
    sfreq = Float64(dat.header.sfreq)
    cond_name = isempty(dat.comments) ? "FIF Epochs" : dat.comments[1]

    if !isempty(dat.nave) && length(dat.nave) == length(dfs) && all(>(1), dat.nave)
        @info "Creating EEG DataFrame (*.fif evoked)"
        for df in dfs
            hasproperty(df, :epoch) && select!(df, Not(:epoch))
        end
        length(dfs) == 1 &&
            return ErpData(fname, 1, cond_name, dfs[1], copy(layout), sr, AnalysisInfo(sample_rate = sfreq), Int(dat.nave[1]))
        return ErpData[
            ErpData(
                fname,
                i,
                i <= length(dat.comments) && !isempty(dat.comments[i]) ? dat.comments[i] : "Condition $i",
                dfs[i],
                copy(layout),
                sr,
                AnalysisInfo(sample_rate = sfreq),
                Int(dat.nave[i]),
            ) for i in eachindex(dfs)
        ]
    else
        @info "Creating EEG DataFrame (*.fif epochs)"
        return EpochData(fname, 1, cond_name, dfs, copy(layout), sr, AnalysisInfo(sample_rate = sfreq))
    end
end

create_eegfun_data(dat::FunctionalImageFormat.FifEpochs) =
    create_eegfun_data(dat, _create_layout_from_labels(Symbol[Symbol(ch.ch_name) for ch in dat.header.ch_info if !_is_stim_channel(ch)]))
