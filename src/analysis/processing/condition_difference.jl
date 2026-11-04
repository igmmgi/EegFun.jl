function _condition_difference_default_output_dir(
    input_dir::String,
    pattern::String,
    pairs::Union{AbstractVector{<:Tuple{Integer,Integer}},AbstractVector{<:AbstractVector{<:Integer}}},
)
    pairs_str = join([join(pair, "-") for pair in pairs], "_")
    return joinpath(input_dir, "differences_$(clean_pattern(pattern))_$(pairs_str)")
end

function _subtract_df(df1::DataFrame, df2::DataFrame, cond_chs::Set{Symbol})::DataFrame
    col_names = propertynames(df1)
    n_rows = nrow(df1)
    out_cols = Vector{Any}(undef, length(col_names))
    for (j, col) in enumerate(col_names)
        c1 = df1[!, j]
        if col in cond_chs
            c2 = df2[!, col]
            out = Vector{Float64}(undef, n_rows)
            if c1 isa Vector{Float64} && c2 isa Vector{Float64}
                v1 = c1::Vector{Float64}
                v2 = c2::Vector{Float64}
                @inbounds @simd for i = 1:n_rows
                    out[i] = v1[i] - v2[i]
                end
            else
                @inbounds @simd for i = 1:n_rows
                    out[i] = Float64(c1[i]) - Float64(c2[i])
                end
            end
            out_cols[j] = out
        else
            out_cols[j] = copy(c1)
        end
    end
    return DataFrame(out_cols, col_names; copycols = false)
end

function _create_difference_wave(dat1::T, dat2::T, cond1::Int, cond2::Int, diff_cond::Int) where {T<:ConditionAverageData}
    _have_same_structure(dat1, dat2) || @minimal_error("Data have inconsistent structure")

    frames1 = _data_frames(dat1)
    frames2 = _data_frames(dat2)
    cond_chs = Set{Symbol}(_condition_channels(dat1))
    for ch in _condition_channels(dat1)
        if !all(df2 -> hasproperty(df2, ch), frames2)
            @minimal_warning "Channel $ch not found in condition $cond2, keeping original values"
            delete!(cond_chs, ch)
        end
    end

    diff_frames = map((df1, df2) -> _subtract_df(df1, df2, cond_chs), frames1, frames2)
    n_epochs = dat1 isa ErpData ? min(dat1.n_epochs, dat2.n_epochs) : 0
    return _derive_condition(dat1, diff_cond, "difference_$(cond1)_$(cond2)", diff_frames; n_epochs = n_epochs)
end

function _condition_difference_process_file(
    filepath::String,
    output_path::String,
    condition_pairs::Union{AbstractVector{<:Tuple{Integer,Integer}},AbstractVector{<:AbstractVector{<:Integer}}},
)
    filename = basename(filepath)
    file_data = read_data(filepath)
    isnothing(file_data) && return BatchResult(false, filename, "No data variables found")

    if file_data isa AbstractVector{<:ErpData} || file_data isa AbstractVector{<:TimeFreqData}
        difference_waves = condition_difference(file_data, condition_pairs)
        jldsave(output_path; data = difference_waves, format_version = EEGFUN_FORMAT_VERSION)
        return BatchResult(true, filename, "Created $(length(difference_waves)) difference(s)")
    end
    return BatchResult(false, filename, "Invalid data type: expected Vector{ErpData} or Vector{TimeFreqData}, got $(typeof(file_data))")
end

"""
    condition_difference(data::AbstractVector{<:Union{ErpData, TimeFreqData}}, condition_pairs) -> Vector
    condition_difference(file_pattern::String, condition_pairs; input_dir = pwd(),
                         participant_selection = participants(), output_dir = nothing)

Create condition difference waves (`a - b` for each `(a, b)` in `condition_pairs`, wrapping TF phase
differences onto `[-π, π]`), or batch-process matching JLD2 files.
"""
function condition_difference(data::AbstractVector{T}, condition_pairs) where {T<:ConditionAverageData}
    pairs =
        if condition_pairs isa Tuple ||
           (condition_pairs isa AbstractVector && !isempty(condition_pairs) && first(condition_pairs) isa Integer)
            [condition_pairs]
        else
            condition_pairs
        end

    tag = _kind_tag(T)
    @info "Creating $(tag)condition difference waves..."

    cond_nums = Int[d.condition for d in data]
    use_cond_field = allunique(cond_nums) && all(>(0), cond_nums)
    difference_waves = T[]

    for (pair_idx, pair) in enumerate(pairs)
        cond1, cond2 = pair[1], pair[2]
        idx1 = use_cond_field ? findfirst(d -> d.condition == cond1, data) : ((1 <= cond1 <= length(data)) ? cond1 : nothing)
        idx2 = use_cond_field ? findfirst(d -> d.condition == cond2, data) : ((1 <= cond2 <= length(data)) ? cond2 : nothing)

        if isnothing(idx1)
            @minimal_warning "Condition $cond1 not found. Skipping pair ($cond1, $cond2)."
            continue
        end
        if isnothing(idx2)
            @minimal_warning "Condition $cond2 not found. Skipping pair ($cond1, $cond2)."
            continue
        end

        push!(difference_waves, _create_difference_wave(data[idx1], data[idx2], cond1, cond2, pair_idx))
    end

    isempty(difference_waves) && @minimal_error("No valid condition pairs found")
    @info "Created $(length(difference_waves)) $(tag)difference wave(s)"
    return difference_waves
end

function condition_difference(
    file_pattern::String,
    condition_pairs::Union{AbstractVector{<:Tuple{Integer,Integer}},AbstractVector{<:AbstractVector{<:Integer}}};
    input_dir::String = pwd(),
    participant_selection = participants(),
    output_dir::Union{String,Nothing} = nothing,
)
    _with_batch_logging("condition differencing", "condition_difference.log", input_dir, output_dir) do ctx
        @log_call "condition_difference"
        error_msg = _validate_condition_pairs(condition_pairs)
        !isnothing(error_msg) && @minimal_error(error_msg)

        output_dir = _prepare_output_dir!(ctx, _condition_difference_default_output_dir(input_dir, file_pattern, condition_pairs))
        process_fn = (input_path, output_path) -> _condition_difference_process_file(input_path, output_path, condition_pairs)
        return batch_process(process_fn, file_pattern, input_dir, output_dir, participant_selection, "Creating differences")
    end
end
