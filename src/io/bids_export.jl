# =============================================================================
# BIDS EXPORT
# =============================================================================

"""
    export_bids(; task, raw_dir, raw_pattern, derivatives_dir, output_dir,
                  layout_file, participant_map, dataset_name, dataset_authors,
                  power_line_frequency, overwrite, ...)

Export an EegFun project (raw data + preprocessed outputs) to a
[BIDS-compliant](https://bids.neuroimaging.io/index.html) directory structure.

The generated dataset can be validated using the
[BIDS Validator](https://bids-standard.github.io/bids-validator/).

# Keyword Arguments
- `task::String`: **Required.** BIDS task label (e.g., `"posner"`).
- `raw_dir::String`: Directory containing raw data files (default: `"."`).
- `raw_pattern::String`: Regex pattern for raw files (default: `"\\\\.bdf"`).
- `derivatives_dir::String`: Pipeline output directory containing JLD2 files (default: `""`). If empty, only raw data is exported.
- `output_dir::String`: Where to create the BIDS dataset (default: `"./bids_dataset"`).
- `layout_file::String`: Electrode layout CSV for electrode positions (default: `""`). If empty, electrode/coordsystem files are skipped.
- `participant_map::Union{Nothing,Dict{String,String}}`: Optional mapping from raw filenames (without extension) to BIDS subject IDs (e.g., `"sub-01"`). If `nothing`, IDs are auto-extracted from filenames.
- `dataset_name::String`: Name for `dataset_description.json` (default: `""`).
- `dataset_authors::Vector{String}`: Authors for `dataset_description.json` (default: `String[]`).
- `dataset_license::String`: License for dataset (default: `"CC0"`).
- `power_line_frequency::Int`: Power line frequency in Hz, 50 (EU) or 60 (US) (default: `50`).
- `overwrite::Bool`: Overwrite existing output directory (default: `false`).

## Recommended metadata (suppresses BIDS validator warnings)
- `task_description::String`: Longer description of the task (default: `""`).
- `instructions::String`: Instructions given to participants (default: `""`).
- `manufacturer::String`: EEG equipment manufacturer, e.g., `"BioSemi"` (default: `""`).
- `manufacturer_model::String`: Equipment model name, e.g., `"ActiveTwo"` (default: `""`).
- `cap_manufacturer::String`: EEG cap manufacturer, e.g., `"BioSemi"` (default: `""`).
- `cap_model::String`: Cap model name, e.g., `"headcap with 72 active electrodes"` (default: `""`).
- `institution_name::String`: Name of the recording institution (default: `""`).
- `institution_address::String`: Address of the institution (default: `""`).
- `institution_department::String`: Department name (default: `""`).
- `eeg_ground::String`: Description of ground electrode location (default: `""`).
- `eeg_placement_scheme::String`: Electrode placement scheme, e.g., `"10-20"` (default: `""`).

# Example
```julia
EegFun.export_bids(
    task = "posner",
    raw_dir = ".",
    derivatives_dir = "./preprocessed_files",
    output_dir = "./bids_dataset",
    layout_file = "biosemi72.csv",
    dataset_name = "Visual Attention (Posner Cueing)",
    dataset_authors = ["Author A", "Author B"],
    manufacturer = "BioSemi",
    manufacturer_model = "ActiveTwo",
    cap_manufacturer = "BioSemi",
    task_description = "Visual attention task using Posner cueing paradigm",
    power_line_frequency = 50,
)
```
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
    # Recommended metadata
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

    # --- Validation ---
    isempty(task) && @minimal_error "task is required (e.g., \"posner\")"
    !isdir(raw_dir) && @minimal_error "raw_dir does not exist: $raw_dir"
    !isempty(derivatives_dir) && !isdir(derivatives_dir) && @minimal_error "derivatives_dir does not exist: $derivatives_dir"

    if isdir(output_dir) && !overwrite
        @minimal_error "Output directory already exists: $output_dir. Use overwrite=true to replace."
    end

    # --- Discover raw files ---
    pattern = Regex(raw_pattern)
    raw_files = sort(filter(f -> occursin(pattern, f), readdir(raw_dir)))
    isempty(raw_files) && @minimal_error "No raw files matching '$raw_pattern' found in $raw_dir"
    @info "Found $(length(raw_files)) raw file(s)"

    # --- Build participant map ---
    pmap = _bids_build_participant_map(raw_files, participant_map)

    # --- Load layout (optional) ---
    layout = nothing
    if !isempty(layout_file) && isfile(layout_file)
        layout = read_layout(layout_file)
        @info "Loaded layout: $(layout_file) ($(nrow(layout.data)) electrodes)"
    elseif !isempty(layout_file)
        @minimal_warning "Layout file not found: $layout_file — electrode files will be skipped"
    end

    # --- Create directory tree ---
    mkpath(output_dir)
    @info "Creating BIDS dataset at: $output_dir"

    # --- Collect recommended metadata into a Dict for sidecar generation ---
    recommended = Dict{String,String}()
    !isempty(task_description) && (recommended["TaskDescription"] = task_description)
    !isempty(instructions) && (recommended["Instructions"] = instructions)
    !isempty(manufacturer) && (recommended["Manufacturer"] = manufacturer)
    !isempty(manufacturer_model) && (recommended["ManufacturersModelName"] = manufacturer_model)
    !isempty(cap_manufacturer) && (recommended["CapManufacturer"] = cap_manufacturer)
    !isempty(cap_model) && (recommended["CapManufacturersModelName"] = cap_model)
    !isempty(institution_name) && (recommended["InstitutionName"] = institution_name)
    !isempty(institution_address) && (recommended["InstitutionAddress"] = institution_address)
    !isempty(institution_department) && (recommended["InstitutionalDepartmentName"] = institution_department)
    !isempty(eeg_ground) && (recommended["EEGGround"] = eeg_ground)
    !isempty(eeg_placement_scheme) && (recommended["EEGPlacementScheme"] = eeg_placement_scheme)

    # --- Top-level files ---
    _bids_write_dataset_description(output_dir, dataset_name, dataset_authors, dataset_license)
    _bids_write_participants_tsv(output_dir, pmap)
    _bids_write_readme(output_dir, dataset_name)

    # --- Per-subject raw data ---
    for raw_file in raw_files
        base = splitext(raw_file)[1]
        sub_id = pmap[base]
        sub_eeg_dir = joinpath(output_dir, sub_id, "eeg")
        mkpath(sub_eeg_dir)

        bids_prefix = "$(sub_id)_task-$(task)"

        # Copy raw file
        raw_src = joinpath(raw_dir, raw_file)
        raw_ext = splitext(raw_file)[2]
        raw_dest = joinpath(sub_eeg_dir, "$(bids_prefix)_eeg$(raw_ext)")
        cp(raw_src, raw_dest, force = overwrite)
        @info "  $sub_id: copied raw data"

        # Read raw data to extract metadata
        try
            raw_data = read_raw_data(raw_src)
            dat = isnothing(layout) ? create_eegfun_data(raw_data) : create_eegfun_data(raw_data, layout)

            # EEG sidecar JSON
            _bids_write_eeg_sidecar(sub_eeg_dir, bids_prefix, dat, task, power_line_frequency, recommended)

            # Channels TSV
            _bids_write_channels_tsv(sub_eeg_dir, bids_prefix, dat)

            # Events TSV
            _bids_write_events_tsv(sub_eeg_dir, bids_prefix, dat)
            _bids_write_events_json(sub_eeg_dir, bids_prefix)

            # Electrode positions (only if layout has real coordinates)
            # BIDS: electrodes/coordsystem are session-level (no task entity)
            if !isnothing(layout) && has_valid_coordinates(layout)
                _bids_write_electrodes_tsv(sub_eeg_dir, sub_id, layout)
                _bids_write_coordsystem_json(sub_eeg_dir, sub_id)
            end

            @info "  $sub_id: wrote sidecar files"
        catch e
            @minimal_warning "  $sub_id: could not read raw data for metadata — $(sprint(showerror, e))"
        end
    end

    # --- Derivatives ---
    if !isempty(derivatives_dir)
        _bids_copy_derivatives(output_dir, derivatives_dir, raw_files, pmap, task, overwrite)
    end

    @info "BIDS export complete: $output_dir"
    @info "Validate your dataset at: https://bids-standard.github.io/bids-validator/"
end


# =============================================================================
# TOP-LEVEL FILES
# =============================================================================

"""Write a README file for the BIDS dataset."""
function _bids_write_readme(output_dir::String, dataset_name::String)
    name = isempty(dataset_name) ? "Untitled Dataset" : dataset_name
    readme_path = joinpath(output_dir, "README")
    open(readme_path, "w") do io
        println(io, name)
        println(io, "=" ^ length(name))
        println(io)
        println(io, "This dataset was generated using EegFun.jl (https://github.com/igmmgi/EegFun.jl)")
        println(io)
        println(io, "This dataset is formatted according to the Brain Imaging Data Structure (BIDS).")
        println(io, "For more information, see https://bids.neuroimaging.io/")
    end
end

# =============================================================================
# PARTICIPANT MAP
# =============================================================================

"""Build a mapping from raw filename basenames to BIDS subject IDs."""
function _bids_build_participant_map(raw_files::Vector{String}, user_map::Union{Nothing,Dict{String,String}})
    pmap = Dict{String,String}()

    for raw_file in raw_files
        base = splitext(raw_file)[1]

        if !isnothing(user_map) && haskey(user_map, base)
            pmap[base] = user_map[base]
        else
            pid = _extract_participant_id(raw_file)
            pmap[base] = @sprintf("sub-%02d", pid)
        end
    end

    # Check for duplicate subject IDs
    sub_ids = collect(values(pmap))
    if length(unique(sub_ids)) != length(sub_ids)
        duplicates = filter(id -> count(==(id), sub_ids) > 1, unique(sub_ids))
        @minimal_warning "Duplicate BIDS subject IDs detected: $(join(duplicates, ", ")). Provide a participant_map to fix."
    end

    return pmap
end


# =============================================================================
# DATASET-LEVEL FILES
# =============================================================================

"""Write dataset_description.json."""
function _bids_write_dataset_description(output_dir::String, name::String, authors::Vector{String}, license::String)
    desc = OrderedDict{String,Any}(
        "Name" => isempty(name) ? "Untitled Dataset" : name,
        "BIDSVersion" => "1.9.0",
        "DatasetType" => "raw",
        "License" => isempty(license) ? "CC0" : license,
        "GeneratedBy" => [OrderedDict("Name" => "EegFun.jl")],
    )
    if !isempty(authors)
        desc["Authors"] = authors
    end

    path = joinpath(output_dir, "dataset_description.json")
    _bids_write_json(path, desc)
end

"""Write participants.tsv and participants.json."""
function _bids_write_participants_tsv(output_dir::String, pmap::Dict{String,String})
    sub_ids = sort(collect(values(pmap)))

    # participants.tsv
    tsv_path = joinpath(output_dir, "participants.tsv")
    open(tsv_path, "w") do io
        println(io, "participant_id")
        for sub_id in sub_ids
            println(io, sub_id)
        end
    end

    # participants.json
    json_path = joinpath(output_dir, "participants.json")
    _bids_write_json(json_path, OrderedDict("participant_id" => OrderedDict("Description" => "Unique participant identifier")))
end


# =============================================================================
# PER-SUBJECT SIDECAR FILES
# =============================================================================

"""Write *_eeg.json sidecar."""
function _bids_write_eeg_sidecar(
    dir::String,
    prefix::String,
    dat::ContinuousData,
    task::String,
    power_line_freq::Int,
    recommended::Dict{String,String},
)
    ref = dat.analysis_info.reference
    ref_str = ref == :none ? "n/a" : string(ref)

    sidecar = OrderedDict{String,Any}(
        "TaskName" => task,
        "SamplingFrequency" => dat.sample_rate,
        "EEGReference" => ref_str,
        "PowerLineFrequency" => power_line_freq,
        "SoftwareFilters" => "n/a",
        "HardwareFilters" => "n/a",
        "RecordingType" => "continuous",
        "RecordingDuration" => round(nrow(dat.data) / dat.sample_rate, digits = 2),
        "SubjectArtefactDescription" => "n/a",
    )

    # Add filter info if available
    hp = dat.analysis_info.hp_filter
    lp = dat.analysis_info.lp_filter
    if hp > 0 || lp > 0
        filters = OrderedDict{String,Any}()
        hp > 0 && (filters["HighpassFilter"] = OrderedDict("CutoffFrequency" => hp))
        lp > 0 && (filters["LowpassFilter"] = OrderedDict("CutoffFrequency" => lp))
        sidecar["SoftwareFilters"] = filters
    end

    # Channel counts (auto-computed from channels.tsv classification)
    all_channels = vcat(channel_labels(dat), extra_labels(dat))
    ch_types = [_bids_classify_channel(ch) for ch in all_channels]
    sidecar["EEGChannelCount"] = count(==("EEG"), ch_types)
    sidecar["EOGChannelCount"] = count(==("EOG"), ch_types)
    sidecar["ECGChannelCount"] = count(==("ECG"), ch_types)
    sidecar["EMGChannelCount"] = count(==("EMG"), ch_types)
    sidecar["MISCChannelCount"] = count(==("MISC"), ch_types)
    sidecar["TriggerChannelCount"] = count(==("TRIG"), ch_types)

    # Add user-provided recommended fields
    for (key, value) in recommended
        sidecar[key] = value
    end

    _bids_write_json(joinpath(dir, "$(prefix)_eeg.json"), sidecar)
end

"""Write *_channels.tsv."""
function _bids_write_channels_tsv(dir::String, prefix::String, dat::ContinuousData)
    ch_labels = channel_labels(dat)
    extra = extra_labels(dat)

    path = joinpath(dir, "$(prefix)_channels.tsv")
    open(path, "w") do io
        println(io, "name\ttype\tunits\tstatus")
        for ch in vcat(ch_labels, extra)
            ch_type = _bids_classify_channel(ch)
            println(io, "$(ch)\t$(ch_type)\tµV\tgood")
        end
    end
end

"""Write *_events.tsv from trigger column."""
function _bids_write_events_tsv(dir::String, prefix::String, dat::ContinuousData)
    df = dat.data
    !hasproperty(df, :trigger) && return  # no triggers to write

    triggers = df.trigger
    time_col = df.time

    path = joinpath(dir, "$(prefix)_events.tsv")
    open(path, "w") do io
        println(io, "onset\tduration\ttrial_type\tvalue")
        for i in eachindex(triggers)
            trig = triggers[i]
            trig == 0 && continue  # skip non-events
            onset = round(time_col[i], digits = 6)
            println(io, "$(onset)\t0\tn/a\t$(Int(trig))")
        end
    end
end

"""Write *_events.json sidecar to describe additional columns."""
function _bids_write_events_json(dir::String, prefix::String)
    # Only define columns not already specified by the BIDS schema
    # (onset, duration, trial_type are built-in — do NOT redefine them)
    events_desc =
        OrderedDict{String,Any}("value" => OrderedDict("Description" => "Trigger value sent by the stimulus presentation software"))
    _bids_write_json(joinpath(dir, "$(prefix)_events.json"), events_desc)
end

"""Write *_electrodes.tsv."""
function _bids_write_electrodes_tsv(dir::String, prefix::String, layout::Layout)
    _ensure_coordinates_3d!(layout)

    path = joinpath(dir, "$(prefix)_electrodes.tsv")
    open(path, "w") do io
        println(io, "name\tx\ty\tz")
        for row in eachrow(layout.data)
            x = round(row.x3, digits = 6)
            y = round(row.y3, digits = 6)
            z = round(row.z3, digits = 6)
            println(io, "$(row.label)\t$(x)\t$(y)\t$(z)")
        end
    end
end

"""Write *_coordsystem.json."""
function _bids_write_coordsystem_json(dir::String, prefix::String)
    coord = OrderedDict{String,Any}(
        "EEGCoordinateSystem" => "Other",
        "EEGCoordinateUnits" => "n/a",
        "EEGCoordinateSystemDescription" => "Spherical coordinates converted from inclination/azimuth via EegFun.jl polar_to_cartesian_xyz!",
    )
    _bids_write_json(joinpath(dir, "$(prefix)_coordsystem.json"), coord)
end


# =============================================================================
# DERIVATIVES
# =============================================================================

# Mapping from pipeline file suffixes to BIDS desc- labels
const _BIDS_DERIVATIVE_SUFFIXES = OrderedDict(
    "_continuous_original" => "continuousOriginal",
    "_continuous_cleaned"  => "continuousCleaned",
    "_epochs_original"     => "epochsOriginal",
    "_epochs_cleaned"      => "epochsCleaned",
    "_epochs_good"         => "epochsGood",
    "_erps_original"       => "erpsOriginal",
    "_erps_cleaned"        => "erpsCleaned",
    "_erps_good"           => "erpsGood",
    "_ica"                 => "ica",
    "_artifact_info"       => "artifactInfo",
)

"""Copy preprocessed JLD2 files into derivatives/EegFun/sub-XX/eeg/ with BIDS naming."""
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

    # Derivatives dataset_description.json
    desc = OrderedDict{String,Any}(
        "Name" => "EegFun Preprocessed Data",
        "BIDSVersion" => "1.9.0",
        "DatasetType" => "derivative",
        "GeneratedBy" => [OrderedDict("Name" => "EegFun.jl")],
    )
    _bids_write_json(joinpath(deriv_out, "dataset_description.json"), desc)

    # Find all JLD2 files in derivatives
    jld2_files = filter(f -> endswith(f, ".jld2"), readdir(deriv_dir))
    isempty(jld2_files) && return

    n_copied = 0
    for raw_file in raw_files
        base = splitext(raw_file)[1]
        sub_id = pmap[base]
        sub_dir = joinpath(deriv_out, sub_id, "eeg")
        mkpath(sub_dir)

        bids_prefix = "$(sub_id)_task-$(task)"

        # Match JLD2 files belonging to this participant
        participant_files = filter(f -> startswith(f, base), jld2_files)

        for jld2_file in participant_files
            # Try to match a known suffix
            bids_name = nothing
            for (suffix, desc_label) in _BIDS_DERIVATIVE_SUFFIXES
                if startswith(jld2_file, "$(base)$(suffix)")
                    bids_name = "$(bids_prefix)_desc-$(desc_label)_eeg.jld2"
                    break
                end
            end

            # Fallback: keep original suffix
            if isnothing(bids_name)
                remaining = replace(jld2_file, base => "")
                remaining = replace(remaining, r"^_" => "")
                bids_name = "$(bids_prefix)_desc-$(splitext(remaining)[1])_eeg.jld2"
            end

            src = joinpath(deriv_dir, jld2_file)
            dest = joinpath(sub_dir, bids_name)
            cp(src, dest, force = overwrite)
            n_copied += 1
        end
    end

    @info "Copied $n_copied derivative file(s) to: $deriv_out"
end


# =============================================================================
# UTILITY HELPERS
# =============================================================================

"""Classify a channel label as EEG, EOG, EMG, ECG, or TRIG for BIDS channels.tsv."""
function _bids_classify_channel(label::Symbol)::String
    s = uppercase(string(label))

    # EOG channels
    (contains(s, "EOG") || s in ("IO1", "IO2", "LO1", "LO2", "SO1", "SO2")) && return "EOG"

    # EMG channels
    contains(s, "EMG") && return "EMG"

    # ECG channels
    (contains(s, "ECG") || contains(s, "EKG")) && return "ECG"

    # Mastoid / reference channels
    (s in ("M1", "M2", "A1", "A2", "TP9", "TP10")) && return "EEG"

    # Trigger / status channels
    (contains(s, "TRIGGER") || contains(s, "STATUS") || s == "STI") && return "TRIG"

    # Default to EEG
    return "EEG"
end

"""Write a dictionary as pretty-printed JSON to a file."""
function _bids_write_json(path::String, data::OrderedDict)
    open(path, "w") do io
        _json_print(io, data, 0)
        println(io)  # trailing newline
    end
end

"""Simple JSON serialiser (avoids adding JSON.jl as a dependency)."""
function _json_print(io::IO, data::AbstractDict, indent::Int)
    print(io, "{\n")
    entries = collect(data)
    for (i, (key, value)) in enumerate(entries)
        print(io, "  "^(indent + 1), "\"", key, "\": ")
        _json_print(io, value, indent + 1)
        i < length(entries) && print(io, ",")
        print(io, "\n")
    end
    print(io, "  "^indent, "}")
end

function _json_print(io::IO, data::AbstractVector, indent::Int)
    if isempty(data)
        print(io, "[]")
        return
    end
    # Check if all elements are simple (not dicts/arrays)
    all_simple = all(x -> !(x isa AbstractDict || x isa AbstractVector), data)
    if all_simple
        print(io, "[")
        for (i, item) in enumerate(data)
            _json_print(io, item, indent)
            i < length(data) && print(io, ", ")
        end
        print(io, "]")
    else
        print(io, "[\n")
        for (i, item) in enumerate(data)
            print(io, "  "^(indent + 1))
            _json_print(io, item, indent + 1)
            i < length(data) && print(io, ",")
            print(io, "\n")
        end
        print(io, "  "^indent, "]")
    end
end

function _json_print(io::IO, data::AbstractString, _::Int)
    escaped = replace(data, "\\" => "\\\\", "\"" => "\\\"", "\n" => "\\n", "\t" => "\\t", "\r" => "\\r")
    print(io, "\"", escaped, "\"")
end
_json_print(io::IO, data::Bool, _::Int) = print(io, data ? "true" : "false")
_json_print(io::IO, data::Number, _::Int) = print(io, data)
_json_print(io::IO, data::Symbol, _::Int) = print(io, "\"", data, "\"")
_json_print(io::IO, data::Nothing, _::Int) = print(io, "null")
