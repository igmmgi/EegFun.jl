function _parse_channel_label(label::Symbol)
    s = String(label)
    digits_only = filter(isdigit, s)
    letters_only = filter(isletter, s)
    return letters_only, isempty(digits_only) ? nothing : parse(Int, digits_only)
end

function _get_channel_pairs_from_selection(erp_left::ErpData, erp_right::ErpData, channel_selection)::Vector{Tuple{Symbol,Symbol}}
    left_labels = channel_labels(erp_left)
    right_labels = channel_labels(erp_right)
    selected = get_selected_channels(erp_left, channel_selection, include_meta = false, include_extra = false)
    auto_detect = length(selected) == length(left_labels) && selected == left_labels

    candidates = auto_detect ? intersect(left_labels, right_labels) : selected
    pairs = Tuple{Symbol,Symbol}[]

    for ch_left in candidates
        letters, digit = _parse_channel_label(ch_left)
        if isnothing(digit)
            auto_detect || @minimal_warning "Channel $ch_left has no number, cannot pair. Skipping."
            continue
        elseif iseven(digit)
            auto_detect || @minimal_warning "Channel $ch_left is even-numbered (expected odd). Skipping."
            continue
        end

        ch_right = Symbol(letters, digit + 1)
        if ch_right ∉ left_labels
            auto_detect || @minimal_warning "Right pair $ch_right for $ch_left not found in left dataset. Skipping."
            continue
        elseif ch_right ∉ right_labels
            auto_detect || @minimal_warning "Right pair $ch_right for $ch_left not found in right dataset. Skipping."
            continue
        end
        push!(pairs, (ch_left, ch_right))
    end

    return pairs
end

function _create_lrp_layout(original_layout::Layout, lrp_channels::AbstractVector{Symbol})::Layout
    keep_set = Set(lrp_channels)
    mask = Bool[label ∈ keep_set for label in original_layout.data.label]
    filtered_layout_df = original_layout.data[mask, :]

    filtered_neighbors = if !isnothing(original_layout.neighbors)
        filtered_dict = Dict{Symbol,Neighbors}()
        for channel in lrp_channels
            if haskey(original_layout.neighbors, channel)
                filtered_dict[channel] = _subset_neighbors(original_layout.neighbors[channel], keep_set)
            end
        end
        isempty(filtered_dict) ? nothing : filtered_dict
    else
        nothing
    end

    return Layout(filtered_layout_df, filtered_neighbors, original_layout.criterion, original_layout.criterion_type)
end

function _calculate_lrp(erp_left::ErpData, erp_right::ErpData, pairs::Vector{Tuple{Symbol,Symbol}})::ErpData
    lrp_df = DataFrame(; copycols = false)
    for col in meta_labels(erp_left)
        lrp_df[!, col] = copy(erp_left.data[!, col])
    end

    lrp_labels = Vector{Symbol}(undef, 2 * length(pairs))
    for (idx, (ch_left, ch_right)) in enumerate(pairs)
        l_in_l = erp_left.data[!, ch_left]
        r_in_l = erp_left.data[!, ch_right]
        l_in_r = erp_right.data[!, ch_left]
        r_in_r = erp_right.data[!, ch_right]

        col_l = @. 0.5 * ((l_in_r - r_in_r) + (r_in_l - l_in_l))
        lrp_df[!, ch_left] = col_l
        lrp_df[!, ch_right] = .-col_l
        lrp_labels[2*idx-1] = ch_left
        lrp_labels[2*idx] = ch_right
    end

    return ErpData(
        erp_left.file,
        erp_left.condition,
        "lrp",
        lrp_df,
        _create_lrp_layout(erp_left.layout, lrp_labels),
        erp_left.sample_rate,
        copy(erp_left.analysis_info),
        min(erp_left.n_epochs, erp_right.n_epochs),
    )
end

function _validate_lrp_pairs(condition_pairs::AbstractVector{<:Tuple{Integer,Integer}})
    isempty(condition_pairs) && @minimal_error("Condition pairs cannot be empty")
    for (left, right) in condition_pairs
        left < 1 && @minimal_error("Condition indices must be positive, got left=$left")
        right < 1 && @minimal_error("Condition indices must be positive, got right=$right")
    end
end

function _resolve_lrp_cond_idx(erps::AbstractVector{<:ErpData}, cond_nums::Vector{Int}, use_cond_field::Bool, cond::Int, side::String)::Int
    idx = use_cond_field ? findfirst(d -> d.condition == cond, erps) : (1 <= cond <= length(erps) ? cond : nothing)
    isnothing(idx) &&
        @minimal_error("$side condition $cond not found (available conditions: $(use_cond_field ? cond_nums : collect(1:length(erps))))")
    return idx
end

_default_lrp_output_dir(input_dir::String, condition_pairs::AbstractVector{<:Tuple{Integer,Integer}}) =
    joinpath(input_dir, "lrp_" * join(["$(l)-$(r)" for (l, r) in condition_pairs], "_"))

function _process_lrp_file(
    filepath::String,
    output_path::String,
    condition_pairs::AbstractVector{<:Tuple{Integer,Integer}},
    channel_selection,
)
    filename = basename(filepath)
    erps_data = read_data(filepath)
    isnothing(erps_data) && return BatchResult(false, filename, "No data variables found")
    erps_data isa Vector{<:ErpData} ||
        return BatchResult(false, filename, "Invalid data type: expected Vector{ErpData}, got $(typeof(erps_data))")
    isempty(erps_data) && return BatchResult(false, filename, "Empty erps array")

    lrp_results = lrp(erps_data, condition_pairs; channel_selection = channel_selection)
    jldsave(output_path; data = lrp_results, format_version = EEGFUN_FORMAT_VERSION)
    return BatchResult(true, filename, "Calculated LRP for $(length(lrp_results)) pair(s)")
end

"""
    lrp(erp_left::ErpData, erp_right::ErpData; channel_selection = channels()) -> ErpData
    lrp(erps::AbstractVector{<:ErpData}, condition_pairs; condition_selection = conditions(), channel_selection = channels()) -> Vector{ErpData}
    lrp(file_pattern::String, condition_pairs; input_dir = pwd(), channel_selection = channels(), participant_selection = participants(), output_dir = nothing)

Calculate the lateralized readiness potential (LRP) from ERP data using the double-subtraction formula:
- `LRP_odd  = 0.5 * ((odd_right - even_right) + (even_left - odd_left))`
- `LRP_even = -LRP_odd`

# Arguments
- `erp_left`, `erp_right`: Paired left- and right-response `ErpData`
- `channel_selection`: Odd/left hemisphere channel selector; default `channels()` auto-detects all odd/even pairs (e.g., `C3`/`C4`, `CP3`/`CP4`)
- `condition_pairs`: Vector of `(left_cond, right_cond)` integer tuples

# References
- Coles, M. G. H. (1989). Modern mind-brain reading. *Psychophysiology*, 26(3), 251-269.
- de Jong, R. et al. (1988). Use of partial stimulus information. *JEP:HPP*, 14(4), 682-692.
- Oostenveld, R. et al. (2003). Brain symmetry and topographic analysis. *Clin. Neurophysiology*, 114(7), 1194-1202.
"""
function lrp(erp_left::ErpData, erp_right::ErpData; channel_selection = channels())::ErpData
    @info "Calculating lateralized readiness potential (LRP)"
    _have_same_structure(erp_left, erp_right) || @minimal_error("Left and right ERPs have inconsistent structure")

    pairs = _get_channel_pairs_from_selection(erp_left, erp_right, channel_selection)
    isempty(pairs) && @minimal_error(
        "No valid lateral channel pairs found. Check that selected channels have odd numbers and their even counterparts exist."
    )

    @info "Using $(length(pairs)) channel pair(s): $(pairs)"
    lrp_data = _calculate_lrp(erp_left, erp_right, pairs)
    @info "LRP calculation complete"
    return lrp_data
end

function lrp(
    erps::AbstractVector{<:ErpData},
    condition_pairs::AbstractVector{<:Tuple{Integer,Integer}};
    condition_selection = conditions(),
    channel_selection = channels(),
)::Vector{ErpData}
    _validate_lrp_pairs(condition_pairs)

    erps_filtered = erps[_select_conditions_by_field(erps, condition_selection)]
    isempty(erps_filtered) && @minimal_error("No ERP conditions available after applying condition_selection")

    @info "Calculating LRP for $(length(condition_pairs)) condition pair(s) from $(length(erps_filtered)) condition(s)"
    cond_nums = Int[d.condition for d in erps_filtered]
    use_cond_field = allunique(cond_nums) && all(>(0), cond_nums)
    lrp_results = Vector{ErpData}(undef, length(condition_pairs))

    for (idx, (left_raw, right_raw)) in enumerate(condition_pairs)
        left_cond, right_cond = Int(left_raw), Int(right_raw)
        left_idx = _resolve_lrp_cond_idx(erps_filtered, cond_nums, use_cond_field, left_cond, "Left")
        right_idx = _resolve_lrp_cond_idx(erps_filtered, cond_nums, use_cond_field, right_cond, "Right")

        @info "  Processing pair $idx: condition $left_cond (left) vs $right_cond (right)"
        lrp_data = lrp(erps_filtered[left_idx], erps_filtered[right_idx]; channel_selection = channel_selection)
        lrp_data.condition = idx
        lrp_data.condition_name = "lrp_$(left_cond)_$(right_cond)"
        lrp_results[idx] = lrp_data
    end

    @info "LRP calculation complete for all $(length(condition_pairs)) pair(s)"
    return lrp_results
end

function lrp(
    file_pattern::String,
    condition_pairs::AbstractVector{<:Tuple{Integer,Integer}};
    input_dir::String = pwd(),
    channel_selection = channels(),
    participant_selection = participants(),
    output_dir::Union{String,Nothing} = nothing,
)
    _with_batch_logging("LRP calculation", "lrp.log", input_dir, output_dir) do ctx
        @log_call "lrp"
        _validate_lrp_pairs(condition_pairs)
        output_dir = _prepare_output_dir!(ctx, _default_lrp_output_dir(input_dir, condition_pairs))
        process_fn = (input_path, output_path) -> _process_lrp_file(input_path, output_path, condition_pairs, channel_selection)
        batch_process(process_fn, file_pattern, input_dir, output_dir, participant_selection, "Calculating LRP")
    end
end
