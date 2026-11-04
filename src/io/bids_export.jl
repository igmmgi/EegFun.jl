"""
    export_bids(; task, raw_dir = ".", raw_pattern = "\\\\.bdf", derivatives_dir = "",
                  output_dir = "./bids_dataset", layout_file = "", participant_map = nothing,
                  dataset_name = "", dataset_authors = String[], dataset_license = "CC0",
                  power_line_frequency = 50, overwrite = false, ...)

Export an EegFun project (raw data and optional preprocessed derivatives) to a
[BIDS-compliant](https://bids.neuroimaging.io/index.html) directory structure.

# Keyword Arguments
- `task::String`: **Required.** BIDS task label (e.g., `"posner"`).
- `raw_dir::String`: Directory containing raw data files (default: `"."`).
- `raw_pattern::String`: Regex pattern for raw files (default: `"\\\\.bdf"`).
- `derivatives_dir::String`: Pipeline output directory containing JLD2 files (default: `""`).
- `output_dir::String`: Output BIDS directory (default: `"./bids_dataset"`).
- `layout_file::String`: Electrode layout CSV for 3D electrode positions (default: `""`).
- `participant_map::Union{Nothing,Dict{String,String}}`: Optional mapping from raw basenames (without extension) to BIDS subject IDs (e.g., `"sub-01"`).
- `dataset_name::String`: Name for `dataset_description.json` (default: `""`).
- `dataset_authors::Vector{String}`: Authors for `dataset_description.json` (default: `String[]`).
- `dataset_license::String`: Dataset license (default: `"CC0"`).
- `power_line_frequency::Int`: Power line frequency in Hz (default: `50`).
- `overwrite::Bool`: Overwrite existing output directory (default: `false`).
- Recommended sidecar metadata keywords: `task_description`, `instructions`, `manufacturer`,
  `manufacturer_model`, `cap_manufacturer`, `cap_model`, `institution_name`,
  `institution_address`, `institution_department`, `eeg_ground`, `eeg_placement_scheme`.
"""
function export_bids(;
    task::String,
    raw_dir::String = ".",
    raw_pattern::String = "\\.bdf",
    derivatives_dir::String = "",
    output_dir::String = "./bids_dataset",
    layout_file::String = "",
    participant_map::Union{Nothing,Dict{String,String}} = nothing,
    dataset_name::String = "",
    dataset_authors::Vector{String} = String[],
    dataset_license::String = "CC0",
    power_line_frequency::Int = 50,
    overwrite::Bool = false,
    task_description::String = "",
    instructions::String = "",
    manufacturer::String = "",
    manufacturer_model::String = "",
    cap_manufacturer::String = "",
    cap_model::String = "",
    institution_name::String = "",
    institution_address::String = "",
    institution_department::String = "",
    eeg_ground::String = "",
    eeg_placement_scheme::String = "",
)
    isempty(task) && @minimal_error "task is required (e.g., \"posner\")"
    !isdir(raw_dir) && @minimal_error "raw_dir does not exist: $raw_dir"
    !isempty(derivatives_dir) && !isdir(derivatives_dir) && @minimal_error "derivatives_dir does not exist: $derivatives_dir"
    (isdir(output_dir) && !overwrite) && @minimal_error "Output directory already exists: $output_dir. Use overwrite=true to replace."

    pattern = Regex(raw_pattern)
    raw_files = sort!(filter(f -> occursin(pattern, f), readdir(raw_dir)))
    isempty(raw_files) && @minimal_error "No raw files matching '$raw_pattern' found in $raw_dir"
    @info "Found $(length(raw_files)) raw file(s)"

    pmap = _bids_build_participant_map(raw_files, participant_map)

    layout = nothing
    if !isempty(layout_file)
        if isfile(layout_file)
            layout = read_layout(layout_file)
            @info "Loaded layout: $(layout_file) ($(nrow(layout.data)) electrodes)"
        else
            @minimal_warning "Layout file not found: $layout_file — electrode files will be skipped"
        end
    end

    mkpath(output_dir)
    @info "Creating BIDS dataset at: $output_dir"

    recommended = Dict{String,String}()
    for (bids_key, val) in (
        ("TaskDescription", task_description),
        ("Instructions", instructions),
        ("Manufacturer", manufacturer),
        ("ManufacturersModelName", manufacturer_model),
        ("CapManufacturer", cap_manufacturer),
        ("CapManufacturersModelName", cap_model),
        ("InstitutionName", institution_name),
        ("InstitutionAddress", institution_address),
        ("InstitutionalDepartmentName", institution_department),
        ("EEGGround", eeg_ground),
        ("EEGPlacementScheme", eeg_placement_scheme),
    )
        !isempty(val) && (recommended[bids_key] = val)
    end

    _bids_write_dataset_description(output_dir, dataset_name, dataset_authors, dataset_license)
    _bids_write_participants_tsv(output_dir, pmap)
    _bids_write_readme(output_dir, dataset_name)

    for raw_file in raw_files
        base, raw_ext = splitext(raw_file)
        sub_id = pmap[base]
        sub_eeg_dir = joinpath(output_dir, sub_id, "eeg")
        mkpath(sub_eeg_dir)

        bids_prefix = "$(sub_id)_task-$(task)"
        raw_src = joinpath(raw_dir, raw_file)
        cp(raw_src, joinpath(sub_eeg_dir, "$(bids_prefix)_eeg$(raw_ext)"); force = overwrite)
        @info "  $sub_id: copied raw data"

        try
            raw_data = read_raw_data(raw_src)
            dat = isnothing(layout) ? create_eegfun_data(raw_data) : create_eegfun_data(raw_data, layout)

            _bids_write_eeg_sidecar(sub_eeg_dir, bids_prefix, dat, task, power_line_frequency, recommended)
            _bids_write_channels_tsv(sub_eeg_dir, bids_prefix, dat)
            _bids_write_events_tsv(sub_eeg_dir, bids_prefix, dat)
            _bids_write_events_json(sub_eeg_dir, bids_prefix)

            if !isnothing(layout) && has_valid_coordinates(layout)
                _bids_write_electrodes_tsv(sub_eeg_dir, sub_id, layout)
                _bids_write_coordsystem_json(sub_eeg_dir, sub_id)
            end
            @info "  $sub_id: wrote sidecar files"
        catch e
            @minimal_warning "  $sub_id: could not read raw data for metadata — $(sprint(showerror, e))"
        end
    end

    !isempty(derivatives_dir) && _bids_copy_derivatives(output_dir, derivatives_dir, raw_files, pmap, task, overwrite)

    @info "BIDS export complete: $output_dir"
    @info "Validate your dataset at: https://bids-standard.github.io/bids-validator/"
end

function _bids_write_readme(output_dir::String, dataset_name::String)
    name = isempty(dataset_name) ? "Untitled Dataset" : dataset_name
    open(joinpath(output_dir, "README"), "w") do io
        println(io, name)
        println(io, "="^length(name))
        println(io)
        println(io, "This dataset was generated using EegFun.jl (https://github.com/igmmgi/EegFun.jl)")
        println(io)
        println(io, "This dataset is formatted according to the Brain Imaging Data Structure (BIDS).")
        println(io, "For more information, see https://bids.neuroimaging.io/")
    end
end

function _bids_build_participant_map(raw_files::Vector{String}, user_map::Union{Nothing,Dict{String,String}})
    pmap = Dict{String,String}()
    for raw_file in raw_files
        base = splitext(raw_file)[1]
        pmap[base] =
            (!isnothing(user_map) && haskey(user_map, base)) ? user_map[base] : @sprintf("sub-%02d", _extract_participant_id(raw_file))
    end

    sub_ids = collect(values(pmap))
    if length(unique(sub_ids)) != length(sub_ids)
        duplicates = filter(id -> count(==(id), sub_ids) > 1, unique(sub_ids))
        @minimal_warning "Duplicate BIDS subject IDs detected: $(join(duplicates, ", ")). Provide a participant_map to fix."
    end
    return pmap
end

function _bids_write_dataset_description(output_dir::String, name::String, authors::Vector{String}, license::String)
    desc = OrderedDict{String,Any}(
        "Name" => isempty(name) ? "Untitled Dataset" : name,
        "BIDSVersion" => "1.9.0",
        "DatasetType" => "raw",
        "License" => isempty(license) ? "CC0" : license,
        "GeneratedBy" => [OrderedDict("Name" => "EegFun.jl")],
    )
    !isempty(authors) && (desc["Authors"] = authors)
    _bids_write_json(joinpath(output_dir, "dataset_description.json"), desc)
end

function _bids_write_participants_tsv(output_dir::String, pmap::Dict{String,String})
    sub_ids = sort!(collect(values(pmap)))
    open(joinpath(output_dir, "participants.tsv"), "w") do io
        println(io, "participant_id")
        for sub_id in sub_ids
            println(io, sub_id)
        end
    end
    _bids_write_json(
        joinpath(output_dir, "participants.json"),
        OrderedDict("participant_id" => OrderedDict("Description" => "Unique participant identifier")),
    )
end

function _bids_count_channel_types(ch_lists...)
    n_eeg = n_eog = n_ecg = n_emg = n_misc = n_trig = 0
    for list in ch_lists, ch in list
        t = _bids_classify_channel(ch)
        if t == "EEG"
            n_eeg += 1
        elseif t == "EOG"
            n_eog += 1
        elseif t == "ECG"
            n_ecg += 1
        elseif t == "EMG"
            n_emg += 1
        elseif t == "TRIG"
            n_trig += 1
        else
            n_misc += 1
        end
    end
    return n_eeg, n_eog, n_ecg, n_emg, n_misc, n_trig
end

function _bids_write_eeg_sidecar(
    dir::String,
    prefix::String,
    dat::ContinuousData,
    task::String,
    power_line_freq::Int,
    recommended::Dict{String,String},
)
    ref = dat.analysis_info.reference
    sidecar = OrderedDict{String,Any}(
        "TaskName" => task,
        "SamplingFrequency" => dat.sample_rate,
        "EEGReference" => ref == :none ? "n/a" : string(ref),
        "PowerLineFrequency" => power_line_freq,
        "SoftwareFilters" => "n/a",
        "HardwareFilters" => "n/a",
        "RecordingType" => "continuous",
        "RecordingDuration" => round(n_samples(dat) / dat.sample_rate; digits = 2),
        "SubjectArtefactDescription" => "n/a",
    )

    hp, lp = dat.analysis_info.hp_filter, dat.analysis_info.lp_filter
    if hp > 0 || lp > 0
        filters = OrderedDict{String,Any}()
        hp > 0 && (filters["HighpassFilter"] = OrderedDict("CutoffFrequency" => hp))
        lp > 0 && (filters["LowpassFilter"] = OrderedDict("CutoffFrequency" => lp))
        sidecar["SoftwareFilters"] = filters
    end

    n_eeg, n_eog, n_ecg, n_emg, n_misc, n_trig = _bids_count_channel_types(channel_labels(dat), extra_labels(dat))
    sidecar["EEGChannelCount"] = n_eeg
    sidecar["EOGChannelCount"] = n_eog
    sidecar["ECGChannelCount"] = n_ecg
    sidecar["EMGChannelCount"] = n_emg
    sidecar["MISCChannelCount"] = n_misc
    sidecar["TriggerChannelCount"] = n_trig

    for (key, value) in recommended
        sidecar[key] = value
    end

    _bids_write_json(joinpath(dir, "$(prefix)_eeg.json"), sidecar)
end

function _bids_write_channels_tsv(dir::String, prefix::String, dat::ContinuousData)
    ch_labels = channel_labels(dat)
    extra = extra_labels(dat)
    open(joinpath(dir, "$(prefix)_channels.tsv"), "w") do io
        println(io, "name\ttype\tunits\tstatus")
        for list in (ch_labels, extra), ch in list
            println(io, ch, '\t', _bids_classify_channel(ch), "\tµV\tgood")
        end
    end
end

function _write_events_rows(io::IO, triggers::AbstractVector, time_col::AbstractVector, trig_info)
    @inbounds for i in eachindex(triggers)
        trig = triggers[i]
        trig == 0 && continue
        onset = round(time_col[i]; digits = 6)
        tt = if !isnothing(trig_info) && !ismissing(trig_info[i])
            s = strip(string(trig_info[i]))
            isempty(s) ? "n/a" : replace(s, '\t' => ' ', '\n' => ' ', '\r' => ' ')
        else
            "n/a"
        end
        println(io, onset, "\t0\t", tt, '\t', Int(trig))
    end
end

function _bids_write_events_tsv(dir::String, prefix::String, dat::ContinuousData)
    df = dat.data
    hasproperty(df, :trigger) || return
    trig_info = hasproperty(df, :trigger_info) ? df.trigger_info : nothing
    open(joinpath(dir, "$(prefix)_events.tsv"), "w") do io
        println(io, "onset\tduration\ttrial_type\tvalue")
        _write_events_rows(io, df.trigger, df.time, trig_info)
    end
end

function _bids_write_events_json(dir::String, prefix::String)
    events_desc =
        OrderedDict{String,Any}("value" => OrderedDict("Description" => "Trigger value sent by the stimulus presentation software"))
    _bids_write_json(joinpath(dir, "$(prefix)_events.json"), events_desc)
end

function _write_electrode_rows(io::IO, labels::AbstractVector, xs::AbstractVector, ys::AbstractVector, zs::AbstractVector)
    @inbounds for i in eachindex(labels)
        println(io, labels[i], '\t', round(xs[i]; digits = 6), '\t', round(ys[i]; digits = 6), '\t', round(zs[i]; digits = 6))
    end
end

function _bids_write_electrodes_tsv(dir::String, prefix::String, layout::Layout)
    _ensure_coordinates_3d!(layout)
    df = layout.data
    open(joinpath(dir, "$(prefix)_electrodes.tsv"), "w") do io
        println(io, "name\tx\ty\tz")
        _write_electrode_rows(io, df.label, df.x3, df.y3, df.z3)
    end
end

function _bids_write_coordsystem_json(dir::String, prefix::String)
    coord = OrderedDict{String,Any}(
        "EEGCoordinateSystem" => "Other",
        "EEGCoordinateUnits" => "n/a",
        "EEGCoordinateSystemDescription" => "Spherical coordinates converted from inclination/azimuth via EegFun.jl polar_to_cartesian_xyz!",
    )
    _bids_write_json(joinpath(dir, "$(prefix)_coordsystem.json"), coord)
end

const _BIDS_DERIVATIVE_SUFFIXES = OrderedDict(
    "_continuous_raw" => "continuousRaw",
    "_continuous_corrected" => "continuousCorrected",
    "_continuous" => "continuous",
    "_epochs_raw" => "epochsRaw",
    "_epochs_corrected" => "epochsCorrected",
    "_epochs" => "epochs",
    "_erps_raw" => "erpsRaw",
    "_erps_corrected" => "erpsCorrected",
    "_erps" => "erps",
    "_ica" => "ica",
    "_artifact_info" => "artifactInfo",
)

function _bids_copy_derivatives(
    output_dir::String,
    deriv_dir::String,
    raw_files::Vector{String},
    pmap::Dict{String,String},
    task::String,
    overwrite::Bool,
)
    deriv_out = joinpath(output_dir, "derivatives", "EegFun")
    mkpath(deriv_out)

    desc = OrderedDict{String,Any}(
        "Name" => "EegFun Preprocessed Data",
        "BIDSVersion" => "1.9.0",
        "DatasetType" => "derivative",
        "GeneratedBy" => [OrderedDict("Name" => "EegFun.jl")],
    )
    _bids_write_json(joinpath(deriv_out, "dataset_description.json"), desc)

    jld2_files = filter(f -> endswith(f, ".jld2"), readdir(deriv_dir))
    isempty(jld2_files) && return

    n_copied = 0
    for raw_file in raw_files
        base = basename_without_ext(raw_file)
        sub_id = pmap[splitext(raw_file)[1]]
        sub_dir = joinpath(deriv_out, sub_id, "eeg")
        mkpath(sub_dir)
        bids_prefix = "$(sub_id)_task-$(task)"

        for jld2_file in jld2_files
            startswith(jld2_file, "$(base)_") || continue
            bids_name = nothing
            for (suffix, desc_label) in _BIDS_DERIVATIVE_SUFFIXES
                if jld2_file == "$(base)$(suffix).jld2"
                    bids_name = "$(bids_prefix)_desc-$(desc_label)_eeg.jld2"
                    break
                end
            end
            if isnothing(bids_name)
                custom_desc = splitext(jld2_file[(ncodeunits(base)+2):end])[1]
                bids_name = "$(bids_prefix)_desc-$(custom_desc)_eeg.jld2"
            end
            cp(joinpath(deriv_dir, jld2_file), joinpath(sub_dir, bids_name); force = overwrite)
            n_copied += 1
        end
    end
    @info "Copied $n_copied derivative file(s) to: $deriv_out"
end

function _bids_classify_channel(label::Symbol)::String
    s = uppercase(string(label))
    (contains(s, "EOG") || s in ("IO1", "IO2", "LO1", "LO2", "SO1", "SO2")) && return "EOG"
    contains(s, "EMG") && return "EMG"
    (contains(s, "ECG") || contains(s, "EKG")) && return "ECG"
    s in ("M1", "M2", "A1", "A2", "TP9", "TP10") && return "EEG"
    (contains(s, "TRIGGER") || contains(s, "STATUS") || s == "STI") && return "TRIG"
    return "EEG"
end

function _bids_write_json(path::String, data::OrderedDict)
    open(path, "w") do io
        _json_print(io, data, 0)
        println(io)
    end
end

@inline function _print_json_indent(io::IO, indent::Int)
    for _ = 1:indent
        print(io, "  ")
    end
end

function _json_print(io::IO, data::AbstractDict, indent::Int)
    print(io, "{\n")
    n = length(data)
    i = 0
    for (key, value) in pairs(data)
        i += 1
        _print_json_indent(io, indent + 1)
        print(io, '"', key, "\": ")
        _json_print(io, value, indent + 1)
        i < n && print(io, ',')
        print(io, '\n')
    end
    _print_json_indent(io, indent)
    print(io, '}')
end

function _json_print(io::IO, data::AbstractVector, indent::Int)
    if isempty(data)
        print(io, "[]")
        return
    end
    if all(x -> !(x isa AbstractDict || x isa AbstractVector), data)
        print(io, '[')
        for (i, item) in enumerate(data)
            _json_print(io, item, indent)
            i < length(data) && print(io, ", ")
        end
        print(io, ']')
    else
        print(io, "[\n")
        for (i, item) in enumerate(data)
            _print_json_indent(io, indent + 1)
            _json_print(io, item, indent + 1)
            i < length(data) && print(io, ',')
            print(io, '\n')
        end
        _print_json_indent(io, indent)
        print(io, ']')
    end
end

function _json_print(io::IO, data::AbstractString, _::Int)
    escaped = replace(data, "\\" => "\\\\", "\"" => "\\\"", "\n" => "\\n", "\t" => "\\t", "\r" => "\\r")
    print(io, '"', escaped, '"')
end
_json_print(io::IO, data::Bool, _::Int) = print(io, data ? "true" : "false")
_json_print(io::IO, data::Number, _::Int) = print(io, data)
_json_print(io::IO, data::Symbol, _::Int) = print(io, '"', data, '"')
_json_print(io::IO, ::Nothing, _::Int) = print(io, "null")
