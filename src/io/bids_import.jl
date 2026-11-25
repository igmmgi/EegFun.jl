# Brain Imaging Data Structure (BIDS-EEG) dataset discovery and import.

######################################
######### Types & Constants ##########
######################################

const _BIDS_RAW_EXTENSIONS = (".bdf", ".edf", ".vhdr", ".fif", ".fiff", ".set", ".xdf", ".xdfz", ".xdf.gz")

######################################
########### User Functions ###########
######################################

"""
    list_bids_runs(dataset_dir::String; datatype::String = "eeg") -> DataFrame

Discover all BIDS recordings in `dataset_dir` for `datatype` and return a `DataFrame`
with columns `subject`, `session`, `task`, `run`, `datatype`, and `file`.
"""
function list_bids_runs(dataset_dir::AbstractString; datatype::AbstractString = "eeg")::DataFrame
    dataset_dir_s = String(dataset_dir)
    datatype_s = String(datatype)
    !isdir(dataset_dir_s) && @minimal_error "BIDS dataset directory not found: $dataset_dir_s"

    subjects = String[]
    sessions = String[]
    tasks = String[]
    runs = Int[]
    datatypes = String[]
    files_out = String[]

    sub_dirs =
        sort!(filter(d -> startswith(d, "sub-") && isdir(joinpath(dataset_dir_s, d)), readdir(dataset_dir_s)); by = _natural_sort_key)
    dt_escaped = _bids_regex_escape(datatype_s)
    entity_rx = Regex("^sub-([^_]+)(?:_ses-([^_]+))?_task-([^_]+)(?:_.*?_run-(\\d+)|_run-(\\d+))?(?:_.*)?_$(dt_escaped)\$")

    for sub_dir in sub_dirs
        sub_path = joinpath(dataset_dir_s, sub_dir)
        candidate_dirs = Tuple{String,String}[]

        direct_dt = joinpath(sub_path, datatype_s)
        isdir(direct_dt) && push!(candidate_dirs, ("", direct_dt))

        for entry in sort!(readdir(sub_path); by = _natural_sort_key)
            if startswith(entry, "ses-") && isdir(joinpath(sub_path, entry))
                ses_dt = joinpath(sub_path, entry, datatype_s)
                isdir(ses_dt) && push!(candidate_dirs, (replace(entry, r"^ses-" => ""), ses_dt))
            end
        end

        for (ses_default, dt_dir) in candidate_dirs
            for fname in sort!(readdir(dt_dir); by = _natural_sort_key)
                ext = _get_file_extension(fname)
                ext in _BIDS_RAW_EXTENSIONS || continue
                m = match(entity_rx, fname[1:(end-ncodeunits(ext))])
                isnothing(m) && continue

                run_cap = !isnothing(m.captures[4]) ? m.captures[4] : m.captures[5]
                push!(subjects, String(m.captures[1]))
                push!(sessions, isnothing(m.captures[2]) ? ses_default : String(m.captures[2]))
                push!(tasks, String(m.captures[3]))
                push!(runs, isnothing(run_cap) ? 0 : something(tryparse(Int, run_cap), 0))
                push!(datatypes, datatype_s)
                push!(files_out, joinpath(dt_dir, fname))
            end
        end
    end

    return DataFrame(
        AbstractVector[subjects, sessions, tasks, runs, datatypes, files_out],
        [:subject, :session, :task, :run, :datatype, :file];
        copycols = false,
    )
end

"""
    read_bids(dataset_dir::AbstractString; subject::AbstractString, task::AbstractString, session::AbstractString = "", run::Int = 0, datatype::AbstractString = "eeg") -> ContinuousData

Read an EEG recording from a BIDS dataset, integrating `electrodes.tsv` 3D coordinates and `events.tsv` annotations when present.

# Keyword Arguments
- `subject::AbstractString`: Participant ID (e.g., `"01"` or `"sub-01"`).
- `task::AbstractString`: Task label (e.g., `"posner"`).
- `session::AbstractString`: Optional session label.
- `run::Int`: Optional run number (`0` omits the run entity).
- `datatype::AbstractString`: BIDS data type folder (default: `"eeg"`).
"""
function read_bids(
    dataset_dir::AbstractString;
    subject::AbstractString,
    task::AbstractString,
    session::AbstractString = "",
    run::Int = 0,
    datatype::AbstractString = "eeg",
)
    dataset_dir_s = String(dataset_dir)
    datatype_s = String(datatype)
    sub = replace(String(subject), r"^sub-" => "")
    ses = replace(String(session), r"^ses-" => "")

    sub_prefix = "sub-$(sub)"
    ses_prefix = isempty(ses) ? "" : "_ses-$(ses)"
    base_prefix = "$(sub_prefix)$(ses_prefix)_task-$(task)"

    eeg_dir = isempty(ses) ? joinpath(dataset_dir_s, sub_prefix, datatype_s) : joinpath(dataset_dir_s, sub_prefix, "ses-$(ses)", datatype_s)
    !isdir(eeg_dir) && @minimal_error "BIDS directory not found: $eeg_dir"

    files = sort!(readdir(eeg_dir); by = _natural_sort_key)
    dt_escaped = _bids_regex_escape(datatype_s)
    base_escaped = _bids_regex_escape(base_prefix)
    run_rx =
        run > 0 ? Regex("^$(base_escaped)(?:_.*)?_run-0*$(run)(?:_.*)?_$(dt_escaped)\$") :
        Regex("^$(base_escaped)(?:_(?!run-)[^_]+)*_$(dt_escaped)\$")

    binary_files = filter(files) do f
        ext = _get_file_extension(f)
        ext in _BIDS_RAW_EXTENSIONS && occursin(run_rx, f[1:(end-ncodeunits(ext))])
    end
    if isempty(binary_files)
        expected_prefix = run > 0 ? "$(base_prefix)_run-$(run)_$(datatype)" : "$(base_prefix)_$(datatype)"
        @minimal_error "No raw EEG file found matching prefix '$expected_prefix' in $eeg_dir"
    end

    matched_file = binary_files[1]
    matched_ext = _get_file_extension(matched_file)
    recording_prefix = replace(matched_file[1:(end-ncodeunits(matched_ext))], Regex("_$(dt_escaped)\$") => "")
    binary_path = joinpath(eeg_dir, matched_file)
    @info "Reading BIDS data: $(matched_file)"

    layout = nothing
    electrodes_file = joinpath(eeg_dir, "$(sub_prefix)$(ses_prefix)_electrodes.tsv")
    if !isfile(electrodes_file)
        sub_ses_rx = Regex("^$(_bids_regex_escape("$(sub_prefix)$(ses_prefix)"))(?:_[^_]+)*_electrodes\\.tsv\$")
        elec_idx = findfirst(f -> occursin(sub_ses_rx, f), files)
        !isnothing(elec_idx) && (electrodes_file = joinpath(eeg_dir, files[elec_idx]))
    end
    if isfile(electrodes_file)
        elec_df = DataFrame(CSV.File(electrodes_file))
        if all(col -> hasproperty(elec_df, col), (:name, :x, :y, :z))
            @info "  ↳ Injecting 3D Coordinates from electrodes.tsv"
            layout = _bids_electrodes_to_layout(elec_df)
        end
    end

    dat = if matched_ext == ".set"
        loaded = read_eeglab(binary_path)
        eeg_dat = loaded isa Tuple ? loaded[1] : loaded
        !isnothing(layout) && _align_bids_layout!(eeg_dat, layout)
        eeg_dat
    else
        raw_data = read_raw_data(binary_path)
        d = create_eegfun_data(raw_data)
        !isnothing(layout) && _align_bids_layout!(d, layout)
        d
    end

    events_file = joinpath(eeg_dir, "$(recording_prefix)_events.tsv")
    if !isfile(events_file) && run > 0
        alt_events = joinpath(eeg_dir, "$(base_prefix)_run-$(run)_events.tsv")
        isfile(alt_events) && (events_file = alt_events)
    end
    if isfile(events_file)
        @info "  ↳ Injecting Annotated Events from events.tsv"
        _inject_bids_events!(dat, DataFrame(CSV.File(events_file)))
    end

    return dat
end

######################################
######### Internal Functions #########
######################################

_bids_regex_escape(s::AbstractString) = replace(s, r"([.^$*+?()\[\]{}|\\])" => s"\\\1")

function _bids_parse_coord(v)::Float64
    (ismissing(v) || v === nothing) && return 0.0
    v isa Number && return isnan(v) ? 0.0 : Float64(v)
    s = strip(string(v))
    (isempty(s) || lowercase(s) == "n/a") && return 0.0
    return something(tryparse(Float64, s), 0.0)
end

function _build_bids_layout_df(names_col::AbstractVector, xs::AbstractVector, ys::AbstractVector, zs::AbstractVector)
    n = length(names_col)
    labels = Vector{Symbol}(undef, n)
    x3 = Vector{Float64}(undef, n)
    y3 = Vector{Float64}(undef, n)
    z3 = Vector{Float64}(undef, n)
    inc = Vector{Float64}(undef, n)
    azi = Vector{Float64}(undef, n)

    @inbounds for i = 1:n
        labels[i] = Symbol(names_col[i])
        xi = _bids_parse_coord(xs[i])
        yi = _bids_parse_coord(ys[i])
        zi = _bids_parse_coord(zs[i])
        x3[i] = xi
        y3[i] = yi
        z3[i] = zi
        r = sqrt(xi * xi + yi * yi + zi * zi)
        if r > 0.0
            inc[i] = acosd(clamp(zi / r, -1.0, 1.0))
            azi[i] = (xi == 0.0 && yi == 0.0) ? 0.0 : atand(yi, xi)
        else
            inc[i] = 0.0
            azi[i] = 0.0
        end
    end
    return DataFrame(AbstractVector[labels, inc, azi, x3, y3, z3], [:label, :inc, :azi, :x3, :y3, :z3]; copycols = false)
end

function _bids_electrodes_to_layout(elec_df::DataFrame)::Layout
    layout = Layout(_build_bids_layout_df(elec_df.name, elec_df.x, elec_df.y, elec_df.z), nothing, nothing, nothing)
    polar_to_cartesian_xy!(layout)
    return layout
end

function _try_parse_bids_event_int(v)::Union{Int,Nothing}
    (v === nothing || ismissing(v)) && return nothing
    v isa Integer && return Int(v) != 0 ? Int(v) : nothing
    (v isa Real && !isnan(v)) && (r = round(Int, v); return r != 0 ? r : nothing)
    s_val = strip(string(v))
    (isempty(s_val) || lowercase(s_val) == "n/a") && return nothing
    parsed = tryparse(Int, s_val)
    !isnothing(parsed) && return parsed != 0 ? parsed : nothing
    m = match(r"-?\d+", s_val)
    if !isnothing(m)
        mc = tryparse(Int, m.match)
        !isnothing(mc) && mc != 0 && return mc
    end
    return nothing
end

function _parse_bids_event_value(v)::Int
    c = _try_parse_bids_event_int(v)
    return isnothing(c) ? 1 : c
end

function _parse_bids_event_label(tt_raw, val_raw)::String
    if tt_raw !== nothing && !ismissing(tt_raw)
        tt_str = strip(string(tt_raw))
        (!isempty(tt_str) && lowercase(tt_str) != "n/a") && return tt_str
    end
    if val_raw !== nothing && !ismissing(val_raw)
        v_str = strip(string(val_raw))
        (!isempty(v_str) && lowercase(v_str) != "n/a") && return v_str
    end
    return "event"
end

_inject_bids_events!(dat, ::DataFrame) = nothing

function _inject_bids_events!(dat::ContinuousData, events_df::DataFrame)
    hasproperty(events_df, :onset) || return
    n = nrow(dat.data)
    if !hasproperty(dat.data, :trigger)
        dat.data[!, :trigger] = zeros(Int, n)
    else
        fill!(dat.data.trigger, 0)
    end
    if !hasproperty(dat.data, :trigger_info)
        dat.data[!, :trigger_info] = fill("", n)
    else
        fill!(dat.data.trigger_info, "")
    end

    trig_col = dat.data.trigger
    info_col = dat.data.trigger_info
    onsets = events_df.onset
    vals = hasproperty(events_df, :value) ? events_df.value : nothing
    tts = hasproperty(events_df, :trial_type) ? events_df.trial_type : nothing
    sr = sample_rate(dat)

    n_ev = length(onsets)
    explicit_codes = Vector{Union{Int,Nothing}}(undef, n_ev)
    labels = Vector{String}(undef, n_ev)
    used_codes = Set{Int}()
    unmapped_labels = String[]

    for i = 1:n_ev
        val_raw = isnothing(vals) ? nothing : vals[i]
        tt_raw = isnothing(tts) ? nothing : tts[i]
        c = _try_parse_bids_event_int(val_raw)
        lbl = _parse_bids_event_label(tt_raw, val_raw)
        explicit_codes[i] = c
        labels[i] = lbl
        if !isnothing(c)
            push!(used_codes, c)
        elseif !(lbl in unmapped_labels)
            push!(unmapped_labels, lbl)
        end
    end

    sort!(unmapped_labels)
    label_to_code = Dict{String,Int}()
    next_code = 1
    for lbl in unmapped_labels
        while next_code in used_codes
            next_code += 1
        end
        label_to_code[lbl] = next_code
        push!(used_codes, next_code)
        next_code += 1
    end

    for i = 1:n_ev
        raw_onset = onsets[i]
        ismissing(raw_onset) && continue
        onset_val = raw_onset isa Number ? Float64(raw_onset) : tryparse(Float64, strip(string(raw_onset)))
        (isnothing(onset_val) || isnan(onset_val)) && continue
        samp_idx = round(Int, onset_val * sr) + 1
        (1 <= samp_idx <= n) || continue

        c = explicit_codes[i]
        lbl = labels[i]
        trig_col[samp_idx] = isnothing(c) ? label_to_code[lbl] : c
        info_col[samp_idx] = lbl
    end
end
