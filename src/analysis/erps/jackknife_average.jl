"""
Jackknife averaging for ERP/LRP data.

The jackknife technique creates leave-one-out averages: for each participant,
compute the average of all other participants (excluding that participant).
This is commonly used with LRP data to reduce variance for statistical testing.
"""

# === JACKKNIFE-SPECIFIC VALIDATION ===

"""Validate jackknife parameters."""
function _validate_jackknife_params(erps::AbstractVector{<:ErpData})
    length(erps) < 2 && return "Need at least 2 participants for jackknife averaging"

    # Validate that all ERPs have the same structure
    if !_have_same_structure(erps)
        return "ERPs have inconsistent structure (sample rate, number of samples, or channel labels)"
    end

    return nothing
end

"""Generate default output directory name for jackknife operation."""
function _default_jackknife_output_dir(input_dir::String, pattern::String)
    joinpath(input_dir, "jackknife_$(clean_pattern(pattern))")
end

# === JACKKNIFE-SPECIFIC HELPERS ===

"""
Create jackknife averages: for each participant, average all other participants.

Returns a vector where element `i` is the leave-one-out average of all participants except
participant `i` (`X̄_{-i}`), or when `pseudovalues = true`, the Smulders (2010) pseudovalue
`J_i = N * X̄ - (N - 1) * X̄_{-i}`.
"""
function _create_jackknife_averages(erps::AbstractVector{<:ErpData}; pseudovalues::Bool = false)::Vector{ErpData}
    n_participants = length(erps)

    if n_participants < 2
        @minimal_error("Need at least 2 participants for jackknife averaging, got $(n_participants)")
    end

    # Get metadata columns and EEG channels from first ERP
    first_erp = erps[1]
    metadata_cols = meta_labels(first_erp)
    eeg_channels = setdiff(propertynames(first_erp.data), metadata_cols)
    n_timepoints = n_samples(first_erp)

    # Pre-allocate result vector
    jackknife_erps = ErpData[]

    @info "Creating jackknife $(pseudovalues ? "pseudovalues" : "averages") for $n_participants participants"

    # Step 1: Pre-compute the total sum and grand mean for each channel to enable O(N) jackknifing
    total_sums = Dict{Symbol,Vector{Float64}}()
    grand_means = Dict{Symbol,Vector{Float64}}()
    for ch in eeg_channels
        total_sums[ch] = zeros(Float64, n_timepoints)
        for erp in erps
            total_sums[ch] .+= erp.data[!, ch]
        end
        if pseudovalues
            grand_means[ch] = total_sums[ch] ./ n_participants
        end
    end

    # Pre-calculate total epochs
    total_epochs_all = sum(erp.n_epochs for erp in erps)

    # Step 2: Create jackknife average (or pseudovalue) for each participant
    for excluded_idx = 1:n_participants
        excluded_erp = erps[excluded_idx]

        @info "  Participant $excluded_idx: $(pseudovalues ? "computing pseudovalue from" : "averaging") $(n_participants - 1) other participants"

        # Create a copy of the first ERP's data as the base (for structure and time vector)
        jackknife_data = copy(first_erp.data)

        # Preserve excluded participant/file metadata columns if present in DataFrame
        if hasproperty(jackknife_data, :file)
            jackknife_data[!, :file] .= excluded_erp.file
        end
        if hasproperty(jackknife_data, :participant) && hasproperty(excluded_erp.data, :participant) && !isempty(excluded_erp.data)
            jackknife_data[!, :participant] .= excluded_erp.data[1, :participant]
        end

        # Calculate leave-one-out jackknife average X̄_{-i} = (Total - Excluded) / (N - 1)
        # If pseudovalues = true, apply Smulders (2010): J_i = N * X̄ - (N - 1) * X̄_{-i}
        for ch in eeg_channels
            avg_col = (total_sums[ch] .- excluded_erp.data[!, ch]) ./ (n_participants - 1)
            if pseudovalues
                jackknife_data[!, ch] = n_participants .* grand_means[ch] .- (n_participants - 1) .* avg_col
            else
                jackknife_data[!, ch] = avg_col
            end
        end

        # Calculate epochs for this jackknife average
        total_epochs = total_epochs_all - excluded_erp.n_epochs

        # Get condition info
        cond_name = excluded_erp.condition_name
        condition = excluded_erp.condition
        jackknife_cond_name = pseudovalues ? "$(cond_name)_jackknife_pseudovalue_$(excluded_idx)" : "$(cond_name)_jackknife_$(excluded_idx)"

        # Stamp excluded_erp.file onto the leave-one-out ErpData
        jackknife_erp = ErpData(
            excluded_erp.file,
            condition,
            jackknife_cond_name,
            jackknife_data,
            copy(excluded_erp.layout),
            excluded_erp.sample_rate,
            copy(excluded_erp.analysis_info),
            total_epochs,
        )

        push!(jackknife_erps, jackknife_erp)
    end

    @info "Created $(length(jackknife_erps)) jackknife $(pseudovalues ? "pseudovalues" : "averages")"
    return jackknife_erps
end

"""Select condition indices matching `.condition` numbers (with positional fallback when unset/non-unique)."""
function _select_jackknife_conditions(data::AbstractVector{<:ErpData}, condition_selection)
    isempty(data) && return Int[]
    isnothing(condition_selection) && return collect(1:length(data))
    condition_selection === Colon() && return collect(1:length(data))

    cond_nums = [d.condition for d in data]
    unique_positive = allunique(cond_nums) && all(>(0), cond_nums)

    if condition_selection isa Integer
        return unique_positive ? findall(==(condition_selection), cond_nums) : intersect([Int(condition_selection)], 1:length(data))
    elseif condition_selection isa AbstractVector{<:Integer} || condition_selection isa AbstractRange{<:Integer}
        targets = collect(Int, condition_selection)
        return unique_positive ? findall(in(targets), cond_nums) : intersect(targets, 1:length(data))
    elseif condition_selection isa Function
        sel_type = typeof(condition_selection)
        if hasfield(sel_type, :condition_indices)
            targets = getfield(condition_selection, :condition_indices)
            return unique_positive ? findall(in(targets), cond_nums) : intersect(targets, 1:length(data))
        elseif hasfield(sel_type, :f) && hasfield(typeof(getfield(condition_selection, :f)), :condition_indices)
            targets = getfield(getfield(condition_selection, :f), :condition_indices)
            return unique_positive ? findall(!in(targets), cond_nums) : [i for i = 1:length(data) if i ∉ targets]
        end
    end
    return get_selected_conditions(data, condition_selection)
end

"""
Load ERP/LRP data from multiple files and organize by condition while preserving `(participant_id, file)` alignment.
Returns `(Dict{Int, Vector{Tuple{Int, String, ErpData}}}, Vector{Tuple{Int, String}})`.
"""
function _load_and_group_for_jackknife(files::Vector{String}, input_dir::String, condition_selection)
    entries_by_condition = Dict{Int,Vector{Tuple{Int,String,ErpData}}}()
    file_records = Tuple{Int,String}[]

    for (i, file) in enumerate(files)
        input_path = joinpath(input_dir, file)
        @info "Loading: $file ($i/$(length(files)))"

        participant_id = _extract_participant_id(file)
        push!(file_records, (participant_id, file))

        # Read data (using read_data which finds by type)
        data = read_data(input_path)

        if isnothing(data)
            @minimal_warning "No data variables found in $file. Skipping."
            continue
        end

        # Validate that data is Vector{ErpData} or ErpData
        if !(data isa Union{AbstractVector{<:ErpData},ErpData})
            @minimal_warning "Invalid data type in $file: expected Vector{ErpData} or ErpData, got $(typeof(data)). Skipping."
            continue
        end

        # Handle both single ErpData and Vector{ErpData}
        data_vec = data isa ErpData ? [data] : data

        # Select conditions if specified
        selected_idx = _select_jackknife_conditions(data_vec, condition_selection)
        data_selected = data_vec[selected_idx]

        # Group by condition while tracking exact (participant_id, file) origin
        for erp in data_selected
            cond_num = erp.condition
            entries = get!(() -> Tuple{Int,String,ErpData}[], entries_by_condition, cond_num)
            push!(entries, (participant_id, file, erp))
        end
    end

    return entries_by_condition, file_records
end

# === MAIN API FUNCTIONS ===

"""
    jackknife_average(erps::AbstractVector{<:ErpData}; condition_selection = conditions(), pseudovalues::Bool = false)::Vector{ErpData}
    jackknife_average(measurements::ErpMeasurementsResult; pseudovalues::Bool = true)::ErpMeasurementsResult
    jackknife_average(file_pattern::String; input_dir = pwd(), participant_selection = participants(),
                      condition_selection = conditions(), output_dir = nothing, pseudovalues::Bool = false)

Create leave-one-out jackknife subaverages (or Smulders, 2010 pseudovalues) from ERP/LRP data.

For each participant `i`, computes the leave-one-out average `X̄_{-i}` of all other `N - 1`
participants (excluding `i`). This technique is standard for nonlinear ERP/LRP latency measures
(e.g., fractional peak or fractional area onset latency) where single-participant waveforms have
low signal-to-noise ratio.

# Statistical Testing & Variance Correction
Because leave-one-out subaverages `X̄_{-i}` share `N - 2` participants, their sample variance is
compressed by `(N - 1)^2` relative to single-participant variance:
- **Ulrich & Miller (2001) correction (`pseudovalues = false`)**: When running standard `t`-tests
  or repeated-measures ANOVAs directly on measurements extracted from leave-one-out subaverages
  `X̄_{-i}`, divide `t`-statistics by `N - 1` (`t_corrected = t / (N - 1)`) and `F`-statistics by
  `(N - 1)^2` (`F_corrected = F / (N - 1)^2`).
- **Smulders (2010) pseudovalues (`pseudovalues = true`)**: Converts leave-one-out values `X̄_{-i}`
  into single-participant jackknife pseudovalues:
  ```math
  J_i = N\\bar{X} - (N - 1)\\bar{X}_{-i}
  ```
  where ``\\bar{X}`` is the mean across all ``N`` participants. Standard `t`-tests, ANOVAs,
  correlations, and mixed-effects models can then be applied directly to ``J_i`` **without**
  manual ``t/(N-1)`` or ``F/(N-1)^2`` rescaling (they yield the exact Ulrich & Miller corrected
  test statistics automatically). Note: for nonlinear measures such as onset latency, extract
  measurements on the `pseudovalues = false` leave-one-out ERPs first and then convert the
  resulting `ErpMeasurementsResult` via `jackknife_average(meas; pseudovalues = true)`.

# Arguments
- `erps::AbstractVector{<:ErpData}`: Vector of ERP/LRP data, one per participant
- `condition_selection`: Condition selector predicate (default: `conditions()`)
- `pseudovalues::Bool`: If `true`, return Smulders (2010) pseudovalues `J_i = N*X̄ - (N-1)*X̄_{-i}` (default: `false`)

# References
- Miller, J., Patterson, T., & Ulrich, R. (1998). Jackknife-based method for measuring LRP onset latency differences. *Psychophysiology*, 35(1), 99-115.
- Ulrich, R., & Miller, J. (2001). Using the jackknife-based scoring method for measuring LRP onset effects in factorial designs. *Psychophysiology*, 38(5), 816-827.
- Smulders, F. T. Y. (2010). Simplifying jackknifing of ERPs and getting more out of it. *Psychophysiology*, 47(2), 387-392.
"""
function jackknife_average(erps::AbstractVector{<:ErpData}; condition_selection = conditions(), pseudovalues::Bool = false)::Vector{ErpData}
    @info "Starting jackknife averaging"

    # Apply condition_selection first
    erps_filtered = erps[_select_jackknife_conditions(erps, condition_selection)]

    # Validate inputs
    error_msg = _validate_jackknife_params(erps_filtered)
    if !isnothing(error_msg)
        @minimal_error(error_msg)
    end

    # Create jackknife averages (or pseudovalues)
    jackknife_results = _create_jackknife_averages(erps_filtered; pseudovalues = pseudovalues)

    @info "Jackknife averaging complete"
    return jackknife_results
end

"""
Convert leave-one-out jackknife measurements `X̄_{-i}` in an `ErpMeasurementsResult` into
Smulders (2010) pseudovalues `J_i = N*X̄ - (N - 1)*X̄_{-i}` within each condition.
"""
function jackknife_average(measurements::ErpMeasurementsResult; pseudovalues::Bool = true)::ErpMeasurementsResult
    df = copy(measurements.data)
    !pseudovalues &&
        return ErpMeasurementsResult(df, measurements.analysis_type, measurements.analysis_interval, measurements.baseline_interval)

    meta_cols = Set([
        :file,
        :participant,
        :condition,
        :condition_name,
        :epoch,
        :analysis_type,
        :baseline_interval_start,
        :baseline_interval_end,
        :analysis_interval_start,
        :analysis_interval_end,
    ])
    meas_cols = [col for col in propertynames(df) if col ∉ meta_cols && eltype(df[!, col]) <: Real]

    group_indices = if hasproperty(df, :condition)
        [findall(==(c), df.condition) for c in unique(df.condition)]
    else
        [collect(1:nrow(df))]
    end

    for row_idxs in group_indices
        n = length(row_idxs)
        n < 2 && @minimal_error("Need at least 2 participants per condition for jackknife pseudovalues, got $n")
        for col in meas_cols
            vals = Float64.(df[row_idxs, col])
            x_bar = mean(vals)
            df[row_idxs, col] = n .* x_bar .- (n - 1) .* vals
        end
    end

    return ErpMeasurementsResult(df, measurements.analysis_type, measurements.analysis_interval, measurements.baseline_interval)
end

function jackknife_average(
    file_pattern::String;
    input_dir::String = pwd(),
    participant_selection = participants(),
    condition_selection = conditions(),
    output_dir::Union{String,Nothing} = nothing,
    pseudovalues::Bool = false,
)
    _with_batch_logging("jackknife averaging", "jackknife.log", input_dir, output_dir) do ctx
        @log_call "jackknife_average"

        # Setup directories
        output_dir = _prepare_output_dir!(ctx, _default_jackknife_output_dir(input_dir, file_pattern))

        # Find files
        files = _find_batch_files(file_pattern, input_dir, participant_selection)

        if isempty(files)
            @minimal_warning "No JLD2 files found matching pattern '$file_pattern' in $input_dir"
            return nothing
        end

        if length(files) < 2
            @minimal_warning "Need at least 2 participants for jackknife averaging, found $(length(files))"
            return nothing
        end

        @info "Found $(length(files)) JLD2 files matching pattern '$file_pattern'"

        # Load and group data by condition while preserving (participant_id, file) identity
        entries_by_condition, file_records = _load_and_group_for_jackknife(files, input_dir, condition_selection)

        if isempty(entries_by_condition)
            @minimal_warning "No valid data found in any files"
            return nothing
        end

        sorted_conds = sort(collect(keys(entries_by_condition)))
        @info "Found conditions: $sorted_conds"

        # Create jackknife averages for each condition and map back to exact (participant_id, file)
        jackknife_by_file = Dict{Tuple{Int,String},Vector{ErpData}}()

        for cond_num in sorted_conds
            entries = entries_by_condition[cond_num]
            @info "Creating jackknife averages for condition $cond_num (n=$(length(entries)) participants)"

            if length(entries) < 2
                @minimal_warning "Only $(length(entries)) participant(s) for condition $cond_num. Skipping jackknife."
                continue
            end

            erps = ErpData[entry[3] for entry in entries]
            jackknife_erps = _create_jackknife_averages(erps; pseudovalues = pseudovalues)

            for ((participant_id, file, _), jk_erp) in zip(entries, jackknife_erps)
                vec_for_file = get!(() -> ErpData[], jackknife_by_file, (participant_id, file))
                push!(vec_for_file, jk_erp)
            end
        end

        if isempty(jackknife_by_file)
            @minimal_warning "No jackknife averages created"
            return nothing
        end

        # Save jackknife data: one file per participant, strictly aligned by (participant_id, file)
        for (participant_id, file) in file_records
            participant_jackknife = get(jackknife_by_file, (participant_id, file), ErpData[])

            if !isempty(participant_jackknife)
                output_file = file  # Keep original filename
                output_path = joinpath(output_dir, output_file)

                # If single condition, save as single ErpData, otherwise as Vector
                data_to_save = length(participant_jackknife) == 1 ? participant_jackknife[1] : participant_jackknife
                jldsave(output_path; data = data_to_save, format_version = EEGFUN_FORMAT_VERSION)

                @info "  Saved participant $participant_id: $output_file"
            end
        end

        @info "Jackknife averaging complete!"
        @info "Output saved to: $output_dir"
    end
end
