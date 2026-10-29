"""
Condition difference waves for ERP and time-frequency data (in-memory and batch).
"""

# === DIFFERENCE-SPECIFIC HELPERS ===

"""Generate default output directory name for difference operation."""
function _condition_difference_default_output_dir(
    input_dir::String,
    pattern::String,
    pairs::Union{Vector{Tuple{Int,Int}},Vector{Vector{Int}}},
)
    pairs_str = join([join(pair, "-") for pair in pairs], "_")
    joinpath(input_dir, "differences_$(clean_pattern(pattern))_$(pairs_str)")
end

# === DIFFERENCE-SPECIFIC PROCESSING ===

"""
Create a difference wave by subtracting `dat2` from `dat1` (ERP data, or TF power and phase).
"""
function _create_difference_wave(dat1::T, dat2::T, cond1::Int, cond2::Int, diff_cond::Int) where {T<:ConditionAverageData}
    # Validate that both datasets have the same structure
    _have_same_structure(dat1, dat2) || @minimal_error("Data have inconsistent structure")

    # Copy dat1's channel data (all containers, e.g. TF power and phase) and subtract dat2 in place
    diff_frames = map(copy, _data_frames(dat1))
    frames2 = _data_frames(dat2)
    for ch in _condition_channels(dat1)
        if hasproperty(dat2.data, ch)
            for (df, df2) in zip(diff_frames, frames2)
                df[!, ch] .-= df2[!, ch]
            end
        else
            @minimal_warning "Channel $ch not found in condition $cond2, keeping original values"
        end
    end

    # ERP: n_epochs reflects the minimum (conservative estimate)
    n_epochs = dat1 isa ErpData ? min(dat1.n_epochs, dat2.n_epochs) : 0
    return _derive_condition(dat1, diff_cond, "difference_$(cond1)_$(cond2)", diff_frames; n_epochs = n_epochs)
end

"""
Process a single ERP file through difference wave creation.
Returns BatchResult with success/failure info.
"""
function _condition_difference_process_file(
    filepath::String,
    output_path::String,
    condition_pairs::Union{AbstractVector{<:Tuple{Integer,Integer}},AbstractVector{<:AbstractVector{<:Integer}}},
)
    filename = basename(filepath)

    # Read data (using read_data which finds by type)
    file_data = read_data(filepath)

    if isnothing(file_data)
        return BatchResult(false, filename, "No data variables found")
    end

    # Dispatch to the appropriate in-memory condition_difference based on data type
    if file_data isa AbstractVector{<:ErpData} || file_data isa AbstractVector{<:TimeFreqData}
        difference_waves = condition_difference(file_data, condition_pairs)
        jldsave(output_path; data = difference_waves)
        return BatchResult(true, filename, "Created $(length(difference_waves)) difference(s)")
    else
        return BatchResult(false, filename, "Invalid data type: expected Vector{ErpData} or Vector{TimeFreqData}, got $(typeof(file_data))")
    end
end

# === IN-MEMORY API FUNCTION ===

"""
    condition_difference(data::Vector{<:ErpData}, condition_pairs) -> Vector{ErpData}
    condition_difference(data::Vector{TimeFreqData}, condition_pairs) -> Vector{TimeFreqData}
    condition_difference(file_pattern::String, condition_pairs;
    input_dir::String = pwd(), participant_selection::Function = participants(),
    output_dir = nothing)

Create condition difference waves (or time-frequency differences), or batch-process JLD2 files.

For each pair `(a, b)`, subtracts condition `b` from condition `a`. `condition_pairs` can be
`Vector{Tuple{Int,Int}}` or `Vector{Vector{Int}}`. The `file_pattern` method loads JLD2 files.

# Example
```julia
erps = average_epochs(epochs)
diff_waves = condition_difference(erps, [(1, 2)])

# Batch
condition_difference("erps_unrejected", [(1, 2)])
```
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

    difference_waves = T[]

    for (pair_idx, pair) in enumerate(pairs)
        cond1, cond2 = pair[1], pair[2]

        # Find data for each condition (last match wins)
        dat1 = nothing
        dat2 = nothing
        for d in data
            cond_num = d.condition
            if cond_num == cond1
                dat1 = d
            elseif cond_num == cond2
                dat2 = d
            end
        end

        if isnothing(dat1)
            @minimal_warning "Condition $cond1 not found. Skipping pair ($cond1, $cond2)."
            continue
        end
        if isnothing(dat2)
            @minimal_warning "Condition $cond2 not found. Skipping pair ($cond1, $cond2)."
            continue
        end

        diff_wave = _create_difference_wave(dat1, dat2, cond1, cond2, pair_idx)
        push!(difference_waves, diff_wave)
    end

    if isempty(difference_waves)
        @minimal_error("No valid condition pairs found")
    end

    @info "Created $(length(difference_waves)) $(tag)difference wave(s)"
    return difference_waves
end

# === BATCH API FUNCTION ===

function condition_difference(
    file_pattern::String,
    condition_pairs::Union{AbstractVector{<:Tuple{Integer,Integer}},AbstractVector{<:AbstractVector{<:Integer}}};
    input_dir::String = pwd(),
    participant_selection = participants(),
    output_dir::Union{String,Nothing} = nothing,
)

    _with_batch_logging("condition differencing", "condition_difference.log", input_dir, output_dir) do ctx
        @log_call "condition_difference"

        # Validation

        error_msg = _validate_condition_pairs(condition_pairs)
        if !isnothing(error_msg)
            @minimal_error(error_msg)
        end

        # Setup directories
        output_dir = _prepare_output_dir!(ctx, _condition_difference_default_output_dir(input_dir, file_pattern, condition_pairs))

        process_fn = (input_path, output_path) -> _condition_difference_process_file(input_path, output_path, condition_pairs)
        result = batch_process(process_fn, file_pattern, input_dir, output_dir, participant_selection, "Creating differences")
    end
end
