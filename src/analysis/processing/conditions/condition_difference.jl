# Pairwise condition difference waves for ERP and time-frequency data.

######################################
########### User Functions ###########
######################################

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
    error_msg = _validate_condition_pairs(pairs)
    !isnothing(error_msg) && @minimal_error(error_msg)

    tag = _kind_tag(T)
    @info "Creating $(tag)condition difference waves..."

    difference_waves = T[]
    part_groups = _group_by_participant_file(data)
    if !isnothing(part_groups)
        for pg in part_groups
            _condition_difference_single!(difference_waves, pg, pairs)
        end
    else
        _condition_difference_single!(difference_waves, data, pairs)
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
        _ = output_dir
        @log_call "condition_difference"
        error_msg = _validate_condition_pairs(condition_pairs)
        !isnothing(error_msg) && @minimal_error(error_msg)

        pairs_str = join([join(pair, "-") for pair in condition_pairs], "_")
        default_dir = joinpath(input_dir, "differences_$(_clean_pattern(file_pattern))_$(pairs_str)")
        return _batch_jld2_files(
            data -> (diffs = condition_difference(data, condition_pairs); (diffs, "Created $(length(diffs)) difference(s)")),
            Union{AbstractVector{<:ErpData},AbstractVector{<:TimeFreqData}},
            ctx,
            default_dir,
            file_pattern,
            input_dir,
            participant_selection,
            "Creating differences";
            expected_label = "Vector{ErpData} or Vector{TimeFreqData}",
        )
    end
end

condition_difference(file_pattern::String, condition_pair::Tuple{Integer,Integer}; kwargs...) =
    condition_difference(file_pattern, [condition_pair]; kwargs...)

condition_difference(file_pattern::String, condition_pair::AbstractVector{<:Integer}; kwargs...) =
    condition_difference(file_pattern, [collect(Int, condition_pair)]; kwargs...)

######################################
######### Internal Functions #########
######################################

function _subtract_df(df1::DataFrame, df2::DataFrame, cond_chs::Set{Symbol})::DataFrame
    col_names = propertynames(df1)
    n_rows = nrow(df1)
    out_cols = Vector{Any}(undef, length(col_names))
    for (j, col) in enumerate(col_names)
        c1 = df1[!, j]
        out_cols[j] = col in cond_chs ? _sub_col!(Vector{Float64}(undef, n_rows), c1, df2[!, col]) : copy(c1)
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

function _condition_difference_single!(difference_waves::Vector{T}, data::AbstractVector{T}, pairs) where {T<:ConditionAverageData}
    cond_nums = Int[d.condition for d in data]
    use_cond_field = allunique(cond_nums) && all(>(0), cond_nums)
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
    return difference_waves
end
