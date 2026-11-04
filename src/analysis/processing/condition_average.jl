function _condition_average_default_output_dir(input_dir::String, pattern::String, groups::AbstractVector{<:AbstractVector{<:Integer}})
    groups_str = join([join(group, "-") for group in groups], "_")
    return joinpath(input_dir, "averages_$(clean_pattern(pattern))_$(groups_str)")
end

function _linear_average_df(dfs::AbstractVector{DataFrame}, cond_chs::Set{Symbol})::DataFrame
    df1 = dfs[1]
    col_names = propertynames(df1)
    n_rows = nrow(df1)
    n_dfs = length(dfs)
    inv_n = 1.0 / n_dfs
    out_cols = Vector{Any}(undef, length(col_names))

    for (j, col) in enumerate(col_names)
        c1 = df1[!, j]
        if col in cond_chs
            out = Vector{Float64}(undef, n_rows)
            if c1 isa Vector{Float64}
                v1 = c1::Vector{Float64}
                @inbounds @simd for i = 1:n_rows
                    out[i] = v1[i]
                end
            else
                @inbounds @simd for i = 1:n_rows
                    out[i] = Float64(c1[i])
                end
            end
            for m = 2:n_dfs
                cm = dfs[m][!, col]
                if cm isa Vector{Float64}
                    vm = cm::Vector{Float64}
                    @inbounds @simd for i = 1:n_rows
                        out[i] += vm[i]
                    end
                else
                    @inbounds @simd for i = 1:n_rows
                        out[i] += Float64(cm[i])
                    end
                end
            end
            if n_dfs > 1
                @inbounds @simd for i = 1:n_rows
                    out[i] *= inv_n
                end
            end
            out_cols[j] = out
        else
            out_cols[j] = copy(c1)
        end
    end
    return DataFrame(out_cols, col_names; copycols = false)
end

function _circular_average_df(dfs::AbstractVector{DataFrame}, cond_chs::Set{Symbol})::DataFrame
    df1 = dfs[1]
    col_names = propertynames(df1)
    n_rows = nrow(df1)
    out_cols = Vector{Any}(undef, length(col_names))
    cos_sum = Vector{Float64}(undef, n_rows)

    for (j, col) in enumerate(col_names)
        if col in cond_chs
            sin_sum = zeros(Float64, n_rows)
            fill!(cos_sum, 0.0)
            for df in dfs
                p = df[!, col]
                @inbounds @simd for i = 1:n_rows
                    s, c = sincos(Float64(p[i]))
                    sin_sum[i] += s
                    cos_sum[i] += c
                end
            end
            @inbounds @simd for i = 1:n_rows
                sin_sum[i] = atan(sin_sum[i], cos_sum[i])
            end
            out_cols[j] = sin_sum
        else
            out_cols[j] = copy(df1[!, j])
        end
    end
    return DataFrame(out_cols, col_names; copycols = false)
end

function _create_average_wave(dats::AbstractVector{T}, conditions::AbstractVector{<:Integer}, avg_cond::Int) where {T<:ConditionAverageData}
    for i = 2:length(dats)
        _have_same_structure(dats[1], dats[i]) || @minimal_error("Data have inconsistent structure")
    end

    cond_chs = Set{Symbol}(_condition_channels(dats[1]))
    for ch in _condition_channels(dats[1])
        if !all(d -> all(df -> hasproperty(df, ch), _data_frames(d)), dats)
            @minimal_warning "Channel $ch not found in all conditions, keeping values from first condition"
            delete!(cond_chs, ch)
        end
    end

    avg_frames = if dats[1] isa TimeFreqData
        pow_dfs = DataFrame[d.data_power for d in dats]
        pha_dfs = DataFrame[d.data_phase for d in dats]
        (_linear_average_df(pow_dfs, cond_chs), _circular_average_df(pha_dfs, cond_chs))
    else
        n_containers = length(_data_frames(dats[1]))
        ntuple(k -> _linear_average_df(DataFrame[_data_frames(d)[k] for d in dats], cond_chs), n_containers)
    end

    n_epochs = dats[1] isa ErpData ? sum(d.n_epochs for d in dats) : 0
    return _derive_condition(dats[1], avg_cond, "average_$(join(conditions, "_"))", avg_frames; n_epochs = n_epochs)
end

function _condition_average_process_file(
    filepath::String,
    output_path::String,
    condition_groups::AbstractVector{<:AbstractVector{<:Integer}},
)
    filename = basename(filepath)
    file_data = read_data(filepath)
    isnothing(file_data) && return BatchResult(false, filename, "No data variables found")

    if file_data isa AbstractVector{<:ErpData} || file_data isa AbstractVector{<:TimeFreqData}
        average_waves = condition_average(file_data, condition_groups)
        jldsave(output_path; data = average_waves, format_version = EEGFUN_FORMAT_VERSION)
        return BatchResult(true, filename, "Created $(length(average_waves)) average(s)")
    end
    return BatchResult(false, filename, "Invalid data type: expected Vector{ErpData} or Vector{TimeFreqData}, got $(typeof(file_data))")
end

"""
    condition_average(data::AbstractVector{<:Union{ErpData, TimeFreqData}}, condition_groups) -> Vector
    condition_average(file_pattern::String, condition_groups; input_dir = pwd(),
                      participant_selection = participants(), output_dir = nothing)

Average ERP waveforms or time-frequency representations (linear power, circular phase) across each
group of condition indices in `condition_groups`, or batch-process matching JLD2 files.
"""
function condition_average(data::AbstractVector{T}, condition_groups) where {T<:ConditionAverageData}
    groups = isempty(condition_groups) ? condition_groups : (first(condition_groups) isa Integer ? [condition_groups] : condition_groups)

    tag = _kind_tag(T)
    @info "Creating $(tag)condition average waves..."

    cond_nums = Int[d.condition for d in data]
    use_cond_field = allunique(cond_nums) && all(>(0), cond_nums)
    average_waves = T[]

    for (group_idx, conditions) in enumerate(groups)
        group_data = T[]
        all_found = true
        for cond in conditions
            found_idx = use_cond_field ? findfirst(d -> d.condition == cond, data) : ((1 <= cond <= length(data)) ? cond : nothing)
            if isnothing(found_idx)
                @minimal_warning "Condition $cond not found. Skipping group $conditions."
                all_found = false
                break
            end
            push!(group_data, data[found_idx])
        end
        all_found || continue
        push!(average_waves, _create_average_wave(group_data, conditions, group_idx))
    end

    isempty(average_waves) && @minimal_error("No valid condition groups found")
    @info "Created $(length(average_waves)) $(tag)average wave(s)"
    return average_waves
end

function condition_average(
    file_pattern::String,
    condition_groups::AbstractVector{<:AbstractVector{<:Integer}};
    input_dir::String = pwd(),
    participant_selection = participants(),
    output_dir::Union{String,Nothing} = nothing,
)
    _with_batch_logging("condition averaging", "condition_average.log", input_dir, output_dir) do ctx
        @log_call "condition_average"
        error_msg = _validate_condition_groups(condition_groups)
        !isnothing(error_msg) && @minimal_error(error_msg)

        output_dir = _prepare_output_dir!(ctx, _condition_average_default_output_dir(input_dir, file_pattern, condition_groups))
        process_fn = (input_path, output_path) -> _condition_average_process_file(input_path, output_path, condition_groups)
        return batch_process(process_fn, file_pattern, input_dir, output_dir, participant_selection, "Creating averages")
    end
end

condition_average(file_pattern::String, condition_group::AbstractVector{<:Integer}; kwargs...) =
    condition_average(file_pattern, [collect(Int, condition_group)]; kwargs...)
