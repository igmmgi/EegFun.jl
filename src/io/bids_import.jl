# === BIDS IMPORT ===

const _BIDS_RAW_EXTENSIONS = (".bdf", ".edf", ".vhdr", ".fif", ".fiff", ".set", ".xdf", ".xdfz", ".xdf.gz")

_bids_regex_escape(s::AbstractString) = replace(s, r"([.^$*+?()\[\]{}|\\])" => s"\\\1")

function _bids_parse_coord(v)::Float64
    (ismissing(v) || v === nothing) && return 0.0
    v isa Number && return isnan(v) ? 0.0 : Float64(v)
    s = strip(string(v))
    (isempty(s) || lowercase(s) == "n/a") && return 0.0
    parsed = tryparse(Float64, s)
    return isnothing(parsed) ? 0.0 : parsed
end

function _bids_electrodes_to_layout(elec_df::DataFrame)::Layout
    n = nrow(elec_df)
    x3 = [_bids_parse_coord(v) for v in elec_df.x]
    y3 = [_bids_parse_coord(v) for v in elec_df.y]
    z3 = [_bids_parse_coord(v) for v in elec_df.z]
    inc = Vector{Float64}(undef, n)
    azi = Vector{Float64}(undef, n)

    for i = 1:n
        r = sqrt(x3[i]^2 + y3[i]^2 + z3[i]^2)
        if r > 0.0
            inc[i] = acosd(clamp(z3[i] / r, -1.0, 1.0))
            azi[i] = (x3[i] == 0.0 && y3[i] == 0.0) ? 0.0 : atand(y3[i], x3[i])
        else
            inc[i] = 0.0
            azi[i] = 0.0
        end
    end

    layout_data = DataFrame(label = Symbol.(elec_df.name), inc = inc, azi = azi, x3 = x3, y3 = y3, z3 = z3)
    layout = Layout(layout_data, nothing, nothing, nothing)
    polar_to_cartesian_xy!(layout)
    return layout
end

"""
    list_bids_runs(dataset_dir::String; datatype::String="eeg") -> DataFrame

Discover all BIDS recordings in `dataset_dir` for the given `datatype` and return a `DataFrame`
with columns `subject`, `session`, `task`, `run`, `datatype`, and `file`.
"""
function list_bids_runs(dataset_dir::String; datatype::String = "eeg")::DataFrame
    !isdir(dataset_dir) && @minimal_error "BIDS dataset directory not found: $dataset_dir"

    subjects = String[]
    sessions = String[]
    tasks = String[]
    runs = Int[]
    datatypes = String[]
    files_out = String[]

    sub_dirs = sort(filter(d -> startswith(d, "sub-") && isdir(joinpath(dataset_dir, d)), readdir(dataset_dir)); by = _natural_sort_key)

    dt_escaped = _bids_regex_escape(datatype)
    entity_rx = Regex("^sub-([^_]+)(?:_ses-([^_]+))?_task-([^_]+)(?:_.*?_run-(\\d+)|_run-(\\d+))?(?:_.*)?_$(dt_escaped)\$")

    for sub_dir in sub_dirs
        sub_path = joinpath(dataset_dir, sub_dir)
        candidate_dirs = Tuple{String,String}[]

        direct_dt = joinpath(sub_path, datatype)
        isdir(direct_dt) && push!(candidate_dirs, ("", direct_dt))

        for entry in sort(readdir(sub_path); by = _natural_sort_key)
            if startswith(entry, "ses-") && isdir(joinpath(sub_path, entry))
                ses_dt = joinpath(sub_path, entry, datatype)
                isdir(ses_dt) && push!(candidate_dirs, (replace(entry, r"^ses-" => ""), ses_dt))
            end
        end

        for (ses_default, dt_dir) in candidate_dirs
            for fname in sort(readdir(dt_dir); by = _natural_sort_key)
                ext = get_file_extension(fname)
                ext in _BIDS_RAW_EXTENSIONS || continue
                stem = fname[1:(end-ncodeunits(ext))]
                m = match(entity_rx, stem)
                isnothing(m) && continue

                sub_val = String(m.captures[1])
                ses_val = isnothing(m.captures[2]) ? ses_default : String(m.captures[2])
                task_val = String(m.captures[3])
                run_cap = !isnothing(m.captures[4]) ? m.captures[4] : m.captures[5]
                run_val = isnothing(run_cap) ? 0 : something(tryparse(Int, run_cap), 0)

                push!(subjects, sub_val)
                push!(sessions, ses_val)
                push!(tasks, task_val)
                push!(runs, run_val)
                push!(datatypes, datatype)
                push!(files_out, joinpath(dt_dir, fname))
            end
        end
    end

    return DataFrame(subject = subjects, session = sessions, task = tasks, run = runs, datatype = datatypes, file = files_out)
end

"""
    read_bids(dataset_dir::String; subject::String, task::String, session::String="", run::Int=0, datatype::String="eeg") -> ContinuousData

Read an EEG recording from a BIDS dataset, automatically integrating curated events and 3D layout coordinates if available.

# Keyword Arguments
- `subject::String`: The participant ID (e.g., `"01"` or `"sub-01"`).
- `task::String`: The task label (e.g., `"posner"`).
- `session::String`: The session label (optional).
- `run::Int`: The run number (optional).
- `datatype::String`: The BIDS data type folder (default: `"eeg"`).
"""
function read_bids(dataset_dir::String; subject::String, task::String, session::String = "", run::Int = 0, datatype::String = "eeg")
    # Strip 'sub-' or 'ses-' prefixes if the user provided them accidentally
    sub = replace(subject, r"^sub-" => "")
    ses = replace(session, r"^ses-" => "")

    # 1. Build the BIDS prefix
    sub_prefix = "sub-$(sub)"
    ses_prefix = isempty(ses) ? "" : "_ses-$(ses)"
    task_prefix = "_task-$(task)"
    base_prefix = "$(sub_prefix)$(ses_prefix)$(task_prefix)"

    # 2. Find the target directory
    eeg_dir = isempty(ses) ? joinpath(dataset_dir, sub_prefix, datatype) : joinpath(dataset_dir, sub_prefix, "ses-$(ses)", datatype)

    if !isdir(eeg_dir)
        @minimal_error "BIDS directory not found: $eeg_dir"
    end

    # 3. Find the raw binary file (supporting both _run-1 and zero-padded _run-01)
    files = sort(readdir(eeg_dir); by = _natural_sort_key)
    dt_escaped = _bids_regex_escape(datatype)
    base_escaped = _bids_regex_escape(base_prefix)
    run_rx = run > 0 ? Regex("^$(base_escaped)_run-0*$(run)_$(dt_escaped)\$") : Regex("^$(base_escaped)_$(dt_escaped)\$")

    binary_files = filter(files) do f
        ext = get_file_extension(f)
        ext in _BIDS_RAW_EXTENSIONS || return false
        stem = f[1:(end-ncodeunits(ext))]
        return occursin(run_rx, stem)
    end

    if isempty(binary_files)
        expected_prefix = run > 0 ? "$(base_prefix)_run-$(run)_$(datatype)" : "$(base_prefix)_$(datatype)"
        @minimal_error "No raw EEG file found matching prefix '$expected_prefix' in $eeg_dir"
    end

    matched_file = binary_files[1]
    matched_ext = get_file_extension(matched_file)
    matched_stem = matched_file[1:(end-ncodeunits(matched_ext))]
    recording_prefix = replace(matched_stem, Regex("_$(dt_escaped)\$") => "")

    binary_path = joinpath(eeg_dir, matched_file)
    @info "Reading BIDS data: $(matched_file)"

    # 4. Read electrodes.tsv if available (for exact 3D coordinates)
    layout = nothing
    electrodes_file = joinpath(eeg_dir, "$(sub_prefix)$(ses_prefix)_electrodes.tsv")
    if isfile(electrodes_file)
        elec_df = DataFrame(CSV.File(electrodes_file))
        if all(col -> col in propertynames(elec_df), [:name, :x, :y, :z])
            @info "  ↳ Injecting 3D Coordinates from electrodes.tsv"
            layout = _bids_electrodes_to_layout(elec_df)
        end
    end

    # 5. Read the raw data & create ContinuousData
    dat = if matched_ext == ".set"
        loaded = read_eeglab(binary_path)
        eeg_dat = loaded isa Tuple ? loaded[1] : loaded
        if !isnothing(layout)
            eeg_dat.layout = copy(layout)
        end
        eeg_dat
    else
        raw_data = read_raw_data(binary_path)
        isnothing(layout) ? create_eegfun_data(raw_data) : create_eegfun_data(raw_data, layout)
    end

    # 6. Inject events.tsv if available (overrides raw trigger column)
    events_file = joinpath(eeg_dir, "$(recording_prefix)_events.tsv")
    if !isfile(events_file) && run > 0
        alt_events = joinpath(eeg_dir, "$(base_prefix)_run-$(run)_events.tsv")
        isfile(alt_events) && (events_file = alt_events)
    end

    if isfile(events_file)
        @info "  ↳ Injecting Annotated Events from events.tsv"
        events_df = DataFrame(CSV.File(events_file))

        # Ensure trigger and trigger_info columns exist and erase raw triggers
        if !hasproperty(dat.data, :trigger)
            dat.data[!, :trigger] = zeros(Int, nrow(dat.data))
        else
            dat.data.trigger .= 0
        end
        if !hasproperty(dat.data, :trigger_info)
            dat.data[!, :trigger_info] = fill("", nrow(dat.data))
        else
            dat.data.trigger_info .= ""
        end

        sr = sample_rate(dat)
        for row in eachrow(events_df)
            if hasproperty(row, :onset) && !ismissing(row.onset)
                onset_val = row.onset isa Number ? Float64(row.onset) : tryparse(Float64, strip(string(row.onset)))
                (isnothing(onset_val) || isnan(onset_val)) && continue
                samp_idx = Int(round(onset_val * sr)) + 1
                if 1 <= samp_idx <= n_samples(dat)
                    val = 1
                    info = "event"

                    if hasproperty(row, :value) && !ismissing(row.value)
                        if row.value isa Integer
                            val = Int(row.value)
                        elseif row.value isa Real && !isnan(row.value)
                            val = round(Int, row.value)
                        else
                            s_val = strip(string(row.value))
                            if !isempty(s_val) && lowercase(s_val) != "n/a"
                                parsed = tryparse(Int, s_val)
                                if !isnothing(parsed)
                                    val = parsed
                                else
                                    m = match(r"-?\d+", s_val)
                                    val = !isnothing(m) ? something(tryparse(Int, m.match), 1) : 1
                                end
                            end
                        end
                    end

                    # Prefer trial_type for the string label if available and not "n/a"
                    if hasproperty(row, :trial_type) && !ismissing(row.trial_type)
                        tt_str = strip(string(row.trial_type))
                        if !isempty(tt_str) && lowercase(tt_str) != "n/a"
                            info = tt_str
                        elseif hasproperty(row, :value) && !ismissing(row.value) && lowercase(strip(string(row.value))) != "n/a"
                            info = strip(string(row.value))
                        end
                    elseif hasproperty(row, :value) && !ismissing(row.value)
                        v_str = strip(string(row.value))
                        if !isempty(v_str) && lowercase(v_str) != "n/a"
                            info = v_str
                        end
                    end

                    dat.data.trigger[samp_idx] = val
                    dat.data.trigger_info[samp_idx] = info
                end
            end
        end
    end

    return dat
end
