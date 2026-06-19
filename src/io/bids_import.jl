# =============================================================================
# BIDS IMPORT
# =============================================================================

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
    run_prefix = run > 0 ? "_run-$(run)" : ""

    bids_prefix = "$(sub_prefix)$(ses_prefix)$(task_prefix)$(run_prefix)_$(datatype)"

    # 2. Find the target directory
    eeg_dir = isempty(ses) ? joinpath(dataset_dir, sub_prefix, datatype) : joinpath(dataset_dir, sub_prefix, "ses-$(ses)", datatype)

    if !isdir(eeg_dir)
        @minimal_error "BIDS directory not found: $eeg_dir"
    end

    # 3. Find the raw binary file
    files = readdir(eeg_dir)
    binary_files = filter(f -> startswith(f, bids_prefix) && get_file_extension(f) in [".bdf", ".edf", ".vhdr", ".fif", ".fiff"], files)

    if isempty(binary_files)
        @minimal_error "No raw EEG file found matching prefix '$bids_prefix' in $eeg_dir"
    end

    binary_path = joinpath(eeg_dir, binary_files[1])
    @info "Reading BIDS data: $(binary_files[1])"

    # 4. Read the raw data (matrix only)
    raw_data = read_raw_data(binary_path)

    # 5. Read electrodes.tsv if available (for exact 3D coordinates)
    layout = nothing
    electrodes_file = joinpath(eeg_dir, "$(sub_prefix)$(ses_prefix)_electrodes.tsv")
    if isfile(electrodes_file)
        elec_df = DataFrame(CSV.File(electrodes_file))
        if all(col -> col in propertynames(elec_df), [:name, :x, :y, :z])
            @info "  ↳ Injecting 3D Coordinates from electrodes.tsv"
            layout_data = DataFrame(
                label = Symbol.(elec_df.name),
                x2 = zeros(nrow(elec_df)),
                y2 = zeros(nrow(elec_df)),
                x3 = elec_df.x,
                y3 = elec_df.y,
                z3 = elec_df.z,
            )
            layout = Layout("BIDS_$(sub)", layout_data)
        end
    end

    # 6. Create ContinuousData
    dat = isnothing(layout) ? create_eegfun_data(raw_data) : create_eegfun_data(raw_data, layout)

    # 7. Inject events.tsv if available (overrides raw trigger column)
    events_file = joinpath(eeg_dir, "$(bids_prefix)_events.tsv")
    if isfile(events_file)
        @info "  ↳ Injecting Annotated Events from events.tsv"
        events_df = DataFrame(CSV.File(events_file))

        # Erase raw triggers
        dat.data.trigger .= 0
        dat.data.trigger_info .= ""

        sr = sample_rate(dat)
        for row in eachrow(events_df)
            if hasproperty(row, :onset) && !ismissing(row.onset)
                samp_idx = Int(round(row.onset * sr)) + 1
                if 1 <= samp_idx <= nrow(dat.data)
                    val = 1
                    info = "event"

                    if hasproperty(row, :value) && !ismissing(row.value)
                        val = row.value isa Number ? Int(row.value) : parse(Int, string(row.value))
                    end

                    # Prefer trial_type for the string label if available
                    if hasproperty(row, :trial_type) && !ismissing(row.trial_type)
                        info = string(row.trial_type)
                    elseif hasproperty(row, :value) && !ismissing(row.value)
                        info = string(row.value)
                    end

                    dat.data.trigger[samp_idx] = val
                    dat.data.trigger_info[samp_idx] = info
                end
            end
        end
    end

    return dat
end
