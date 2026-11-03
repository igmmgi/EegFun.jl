const _JACKKNIFE_META_COLS = (
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
)

_default_jackknife_output_dir(input_dir::String, pattern::String) = joinpath(input_dir, "jackknife_$(clean_pattern(pattern))")

function _create_jackknife_averages(erps::AbstractVector{<:ErpData}; pseudovalues::Bool = false)::Vector{ErpData}
    n_participants = length(erps)
    n_participants < 2 && @minimal_error("Need at least 2 participants for jackknife averaging, got $(n_participants)")

    first_erp = erps[1]
    meta_set = Set(meta_labels(first_erp))
    all_cols = propertynames(first_erp.data)
    eeg_channels = Symbol[c for c in all_cols if c ∉ meta_set]
    ch_map = Dict{Symbol,Int}(ch => i for (i, ch) in enumerate(eeg_channels))
    n_tp = n_samples(first_erp)

    @info "Creating jackknife $(pseudovalues ? "pseudovalues" : "averages") for $n_participants participants"

    sums = [zeros(Float64, n_tp) for _ in eeg_channels]
    for (c, ch) in enumerate(eeg_channels)
        s = sums[c]
        for erp in erps
            s .+= erp.data[!, ch]
        end
    end

    total_epochs_all = sum(erp.n_epochs for erp in erps)
    inv_n = 1.0 / n_participants
    inv_nm1 = 1.0 / (n_participants - 1)
    jackknife_erps = Vector{ErpData}(undef, n_participants)

    for excluded_idx = 1:n_participants
        excluded_erp = erps[excluded_idx]
        @info "  Participant $excluded_idx: $(pseudovalues ? "computing pseudovalue from" : "averaging") $(n_participants - 1) other participants"

        jackknife_data = DataFrame(; copycols = false)
        for col in all_cols
            c_idx = get(ch_map, col, 0)
            if c_idx > 0
                s = sums[c_idx]
                x = excluded_erp.data[!, col]
                out = Vector{Float64}(undef, n_tp)
                if pseudovalues
                    @. out = n_participants * (s * inv_n) - (n_participants - 1) * ((s - x) * inv_nm1)
                else
                    @. out = (s - x) * inv_nm1
                end
                jackknife_data[!, col] = out
            elseif col === :file
                jackknife_data[!, :file] = fill(excluded_erp.file, n_tp)
            elseif col === :participant && hasproperty(excluded_erp.data, :participant) && !isempty(excluded_erp.data)
                jackknife_data[!, :participant] = fill(excluded_erp.data[1, :participant], n_tp)
            else
                jackknife_data[!, col] = copy(first_erp.data[!, col])
            end
        end

        cond_name = excluded_erp.condition_name
        jk_name = pseudovalues ? "$(cond_name)_jackknife_pseudovalue_$(excluded_idx)" : "$(cond_name)_jackknife_$(excluded_idx)"
        jackknife_erps[excluded_idx] = ErpData(
            excluded_erp.file,
            excluded_erp.condition,
            jk_name,
            jackknife_data,
            copy(excluded_erp.layout),
            excluded_erp.sample_rate,
            copy(excluded_erp.analysis_info),
            total_epochs_all - excluded_erp.n_epochs,
        )
    end

    @info "Created $(length(jackknife_erps)) jackknife $(pseudovalues ? "pseudovalues" : "averages")"
    return jackknife_erps
end

function _load_and_group_for_jackknife(files::Vector{String}, input_dir::String, condition_selection)
    entries_by_condition = Dict{Int,Vector{Tuple{Int,String,ErpData}}}()
    file_records = Tuple{Int,String}[]

    for (i, file) in enumerate(files)
        @info "Loading: $file ($i/$(length(files)))"
        participant_id = _extract_participant_id(file)
        push!(file_records, (participant_id, file))

        data = read_data(joinpath(input_dir, file))
        if isnothing(data)
            @minimal_warning "No data variables found in $file. Skipping."
            continue
        elseif !(data isa Union{AbstractVector{<:ErpData},ErpData})
            @minimal_warning "Invalid data type in $file: expected Vector{ErpData} or ErpData, got $(typeof(data)). Skipping."
            continue
        end

        data_vec = data isa ErpData ? [data] : data
        for idx in _select_conditions_by_field(data_vec, condition_selection)
            erp = data_vec[idx]
            push!(get!(() -> Tuple{Int,String,ErpData}[], entries_by_condition, erp.condition), (participant_id, file, erp))
        end
    end

    return entries_by_condition, file_records
end

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
- **Ulrich & Miller (2001) correction (`pseudovalues = false`)**: Divide `t`-statistics by `N - 1`
  and `F`-statistics by `(N - 1)^2`.
- **Smulders (2010) pseudovalues (`pseudovalues = true`)**: Converts leave-one-out values `X̄_{-i}`
  into single-participant pseudovalues ``J_i = N\\bar{X} - (N - 1)\\bar{X}_{-i}``, allowing standard
  `t`-tests, ANOVAs, correlations, and mixed-effects models without manual rescaling. For nonlinear
  measures (e.g., onset latency), extract measurements on `pseudovalues = false` subaverages first,
  then convert via `jackknife_average(meas; pseudovalues = true)`.

# References
- Miller, J., Patterson, T., & Ulrich, R. (1998). Jackknife-based method for measuring LRP onset latency differences. *Psychophysiology*, 35(1), 99-115.
- Ulrich, R., & Miller, J. (2001). Using the jackknife-based scoring method for measuring LRP onset effects in factorial designs. *Psychophysiology*, 38(5), 816-827.
- Smulders, F. T. Y. (2010). Simplifying jackknifing of ERPs and getting more out of it. *Psychophysiology*, 47(2), 387-392.
"""
function jackknife_average(erps::AbstractVector{<:ErpData}; condition_selection = conditions(), pseudovalues::Bool = false)::Vector{ErpData}
    @info "Starting jackknife averaging"
    erps_filtered = erps[_select_conditions_by_field(erps, condition_selection)]
    length(erps_filtered) < 2 && @minimal_error("Need at least 2 participants for jackknife averaging")
    _have_same_structure(erps_filtered) ||
        @minimal_error("ERPs have inconsistent structure (sample rate, number of samples, or channel labels)")

    jackknife_results = _create_jackknife_averages(erps_filtered; pseudovalues = pseudovalues)
    @info "Jackknife averaging complete"
    return jackknife_results
end

function jackknife_average(measurements::ErpMeasurementsResult; pseudovalues::Bool = true)::ErpMeasurementsResult
    df = copy(measurements.data)
    !pseudovalues &&
        return ErpMeasurementsResult(df, measurements.analysis_type, measurements.analysis_interval, measurements.baseline_interval)

    meas_cols = Symbol[col for col in propertynames(df) if col ∉ _JACKKNIFE_META_COLS && eltype(df[!, col]) <: Real]
    group_indices = hasproperty(df, :condition) ? [findall(==(c), df.condition) for c in unique(df.condition)] : [collect(1:nrow(df))]

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
        output_dir = _prepare_output_dir!(ctx, _default_jackknife_output_dir(input_dir, file_pattern))

        files = _find_batch_files(file_pattern, input_dir, participant_selection)
        if isempty(files)
            @minimal_warning "No JLD2 files found matching pattern '$file_pattern' in $input_dir"
            return nothing
        elseif length(files) < 2
            @minimal_warning "Need at least 2 participants for jackknife averaging, found $(length(files))"
            return nothing
        end
        @info "Found $(length(files)) JLD2 files matching pattern '$file_pattern'"

        entries_by_condition, file_records = _load_and_group_for_jackknife(files, input_dir, condition_selection)
        if isempty(entries_by_condition)
            @minimal_warning "No valid data found in any files"
            return nothing
        end

        sorted_conds = sort!(collect(keys(entries_by_condition)))
        @info "Found conditions: $sorted_conds"

        jackknife_by_file = Dict{Tuple{Int,String},Vector{ErpData}}()
        for cond_num in sorted_conds
            entries = entries_by_condition[cond_num]
            @info "Creating jackknife averages for condition $cond_num (n=$(length(entries)) participants)"
            if length(entries) < 2
                @minimal_warning "Only $(length(entries)) participant(s) for condition $cond_num. Skipping jackknife."
                continue
            end
            erps = ErpData[entry[3] for entry in entries]
            _have_same_structure(erps) || @minimal_error("ERPs for condition $cond_num have inconsistent structure")
            for ((participant_id, file, _), jk_erp) in zip(entries, _create_jackknife_averages(erps; pseudovalues = pseudovalues))
                push!(get!(() -> ErpData[], jackknife_by_file, (participant_id, file)), jk_erp)
            end
        end

        if isempty(jackknife_by_file)
            @minimal_warning "No jackknife averages created"
            return nothing
        end

        for (participant_id, file) in file_records
            participant_jackknife = get(jackknife_by_file, (participant_id, file), ErpData[])
            isempty(participant_jackknife) && continue
            output_path = joinpath(output_dir, file)
            data_to_save = length(participant_jackknife) == 1 ? participant_jackknife[1] : participant_jackknife
            jldsave(output_path; data = data_to_save, format_version = EEGFUN_FORMAT_VERSION)
            @info "  Saved participant $participant_id: $file"
        end

        @info "Jackknife averaging complete!"
        @info "Output saved to: $output_dir"
    end
end
