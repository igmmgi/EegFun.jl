function _linear_colmeans_idx(df::DataFrame, idxs::Vector{Int})::Vector{Float64}
    n = nrow(df)
    isempty(idxs) && return Float64[]
    out = Vector{Float64}(undef, n)
    c1 = df[!, idxs[1]]
    if c1 isa Vector{Float64}
        v1 = c1::Vector{Float64}
        @inbounds @simd for i = 1:n
            out[i] = v1[i]
        end
    else
        @inbounds @simd for i = 1:n
            out[i] = Float64(c1[i])
        end
    end
    for m = 2:length(idxs)
        cm = df[!, idxs[m]]
        if cm isa Vector{Float64}
            vm = cm::Vector{Float64}
            @inbounds @simd for i = 1:n
                out[i] += vm[i]
            end
        else
            @inbounds @simd for i = 1:n
                out[i] += Float64(cm[i])
            end
        end
    end
    if length(idxs) > 1
        k_f64 = Float64(length(idxs))
        @inbounds @simd for i = 1:n
            out[i] /= k_f64
        end
    end
    return out
end

function _circular_colmeans_idx(df::DataFrame, idxs::Vector{Int})::Vector{Float64}
    n = nrow(df)
    sin_sum = zeros(Float64, n)
    cos_sum = zeros(Float64, n)
    for idx in idxs
        vals = df[!, idx]
        @inbounds @simd for i = 1:n
            s, c = sincos(Float64(vals[i]))
            sin_sum[i] += s
            cos_sum[i] += c
        end
    end
    @inbounds @simd for i = 1:n
        sin_sum[i] = atan(sin_sum[i], cos_sum[i])
    end
    return sin_sum
end

_circular_colmeans(df::DataFrame, cols::AbstractVector{Symbol})::Vector{Float64} =
    _circular_colmeans_idx(df, Int[columnindex(df, c) for c in cols])

function _assign_group_mean!(df::DataFrame, grp::Vector{Symbol}, lbl::Symbol, circular::Bool)
    idxs = Int[columnindex(df, c) for c in grp]
    df[!, lbl] = circular ? _circular_colmeans_idx(df, idxs) : _linear_colmeans_idx(df, idxs)
    return nothing
end

function _assign_group_mean!(dfs::AbstractVector{<:DataFrame}, grp::Vector{Symbol}, lbl::Symbol, circular::Bool)
    isempty(dfs) && return nothing
    idxs = Int[columnindex(dfs[1], c) for c in grp]
    work = length(dfs) * nrow(dfs[1]) * max(1, length(idxs))
    min_work = circular ? 65_536 : 1_000_000
    if length(dfs) >= 8 && Threads.nthreads() > 1 && work >= min_work
        Threads.@threads :static for e in eachindex(dfs)
            dfs[e][!, lbl] = circular ? _circular_colmeans_idx(dfs[e], idxs) : _linear_colmeans_idx(dfs[e], idxs)
        end
    else
        for df in dfs
            df[!, lbl] = circular ? _circular_colmeans_idx(df, idxs) : _linear_colmeans_idx(df, idxs)
        end
    end
    return nothing
end

"""
    channel_average!(dat::EegData; channel_selections = [channels()], output_labels = nothing,
                     include_extra::Bool = false, reduce::Bool = false) -> EegData

Average each group of selected channels in `dat` (linearly for amplitude/power, circularly for TF phase).
Output labels default to `:avg` when all channels are selected, or underscore-joined channel names otherwise.
When `reduce = true`, keeps only metadata and the averaged channels in both data and layout.
"""
function channel_average!(
    dat::EegData;
    channel_selections = nothing,
    channel_selection = nothing,
    output_labels = nothing,
    include_extra::Bool = false,
    reduce::Bool = false,
)
    selections = if !isnothing(channel_selection)
        [channel_selection]
    elseif !isnothing(channel_selections)
        if channel_selections isa AbstractVector{Symbol} ||
           channel_selections isa AbstractVector{<:AbstractString} ||
           channel_selections isa Function
            [channel_selections]
        else
            channel_selections
        end
    else
        [channels()]
    end

    if !isnothing(output_labels)
        output_labels = output_labels isa Union{Symbol,AbstractString} ? [Symbol(output_labels)] : Symbol.(output_labels)
    end

    selected_channel_groups =
        Vector{Symbol}[get_selected_channels(dat, sel; include_meta = false, include_extra = include_extra) for sel in selections]
    if any(isempty, selected_channel_groups)
        @minimal_error "Channel selections at indices $(findall(isempty, selected_channel_groups)) produced no channels"
    end

    if isnothing(output_labels)
        original_channels = channel_labels(dat)
        labels = Vector{Symbol}(undef, length(selected_channel_groups))
        for (i, grp) in enumerate(selected_channel_groups)
            is_all = !isempty(original_channels) && length(grp) == length(original_channels) && all(in(original_channels), grp)
            labels[i] = is_all ? :avg : Symbol(join(string.(grp), "_"))
        end
    else
        length(output_labels) == length(selected_channel_groups) ||
            @minimal_error "N output_labels ($(length(output_labels))) must match N channel selections ($(length(selected_channel_groups)))"
        labels = Symbol.(output_labels)
    end

    frames = _data_frames(dat)
    is_tf = dat isa Union{TimeFreqData,TimeFreqEpochData}
    existing = Set(all_labels(dat))
    for (grp, lbl) in zip(selected_channel_groups, labels)
        lbl in existing && @minimal_warning "Overwriting existing channel '$(lbl)'"
        push!(existing, lbl)
        @info "channel_average!: $(_print_vector(grp)) → :$(lbl)"
        if is_tf
            _assign_group_mean!(frames[1], grp, lbl, false)
            _assign_group_mean!(frames[2], grp, lbl, true)
        else
            for container in frames
                _assign_group_mean!(container, grp, lbl, false)
            end
        end
    end

    avg_layout = _layout_from_groups(dat.layout, selected_channel_groups, labels)
    if reduce
        meta_cols = meta_labels(dat)
        _set_data_frames!(dat, map(container -> _map_df(df -> _reduced_df(df, meta_cols, labels), container), frames))
        dat.layout = avg_layout
    else
        dat.layout = _append_layouts(dat.layout, avg_layout)
    end
    return dat
end

function channel_average!(dat::AbstractVector{<:EegData}; kwargs...)
    for d in dat
        channel_average!(d; kwargs...)
    end
    return dat
end

"""
    channel_average(dat, args...; kwargs...)
    channel_average(file_pattern::String, channel_selections::AbstractVector;
                    output_labels = nothing, input_dir::String = pwd(),
                    participant_selection = participants(), condition_selection = conditions(),
                    output_dir = nothing, reduce::Bool = false)

Non-mutating version of [`channel_average!`](@ref), or batch-process matching JLD2 files into `output_dir`.
"""
channel_average(dat::EegData, args...; kwargs...) = channel_average!(copy(dat), args...; kwargs...)
channel_average(dat::AbstractVector{<:EegData}; kwargs...) = [channel_average(d; kwargs...) for d in dat]

_default_channel_average_output_dir(input_dir::String, pattern::String) = joinpath(input_dir, "averaged_channels_$(clean_pattern(pattern))")

function _process_channel_average_file(
    filepath::String,
    output_path::String,
    channel_selections::AbstractVector,
    output_labels::AbstractVector{Symbol},
    condition_selection,
    reduce::Bool,
)
    filename = basename(filepath)
    data = read_data(filepath)
    isnothing(data) && return BatchResult(false, filename, "No data variables found")
    if !(data isa AbstractVector{<:Union{ErpData,EpochData}})
        return BatchResult(false, filename, "Invalid data type: expected Vector{ErpData} or Vector{EpochData}")
    end

    data = _condition_select(data, condition_selection)
    for item in data
        channel_average!(item; channel_selections = channel_selections, output_labels = output_labels, reduce = reduce)
    end
    jldsave(output_path; data = data, format_version = EEGFUN_FORMAT_VERSION)
    return BatchResult(true, filename, "Averaged $(length(channel_selections)) channel group(s)")
end

function channel_average(
    file_pattern::String,
    channel_selections::AbstractVector;
    output_labels::Union{AbstractVector{Symbol},Nothing} = nothing,
    input_dir::String = pwd(),
    participant_selection = participants(),
    condition_selection = conditions(),
    output_dir::Union{String,Nothing} = nothing,
    reduce::Bool = false,
)
    _with_batch_logging("channel averaging", "channel_average.log", input_dir, output_dir) do ctx
        @log_call "channel_average"
        if isnothing(output_labels)
            output_labels = [Symbol("avg_$i") for i = 1:length(channel_selections)]
        end
        length(output_labels) == length(channel_selections) || @minimal_error(
            "Number of output_labels ($(length(output_labels))) must match number of channel_selections ($(length(channel_selections)))"
        )

        output_dir = _prepare_output_dir!(ctx, _default_channel_average_output_dir(input_dir, file_pattern))
        process_fn =
            (input_path, output_path) ->
                _process_channel_average_file(input_path, output_path, channel_selections, output_labels, condition_selection, reduce)
        return batch_process(process_fn, file_pattern, input_dir, output_dir, participant_selection, "Averaging channels")
    end
end

function _reduced_df(df::DataFrame, meta_cols::AbstractVector{Symbol}, labels::AbstractVector{Symbol})
    n_meta = length(meta_cols)
    n_lbl = length(labels)
    col_names = Vector{Symbol}(undef, n_meta + n_lbl)
    out_cols = Vector{Any}(undef, n_meta + n_lbl)
    for (i, c) in enumerate(meta_cols)
        col_names[i] = c
        out_cols[i] = copy(df[!, c])
    end
    for (j, lbl) in enumerate(labels)
        col_names[n_meta+j] = lbl
        out_cols[n_meta+j] = df[!, lbl]
    end
    return DataFrame(out_cols, col_names; copycols = false)
end

function _layout_from_groups(layout::Layout, channel_groups::Vector{Vector{Symbol}}, labels::Vector{Symbol})::Layout
    df = layout.data
    isempty(df) && return Layout(DataFrame(), nothing, nothing, nothing)

    work_layout = copy(layout)
    _ensure_all_coordinates!(work_layout)
    wdf = work_layout.data
    labels_in_layout = df[!, :label]
    x3 = wdf[!, :x3]
    y3 = wdf[!, :y3]
    z3 = wdf[!, :z3]

    out_labels = Symbol[]
    out_incs = Float64[]
    out_azis = Float64[]
    for (grp, lbl) in zip(channel_groups, labels)
        mx = 0.0
        my = 0.0
        mz = 0.0
        cnt = 0
        for i in eachindex(labels_in_layout)
            if labels_in_layout[i] in grp
                mx += Float64(x3[i])
                my += Float64(y3[i])
                mz += Float64(z3[i])
                cnt += 1
            end
        end
        cnt == 0 && continue
        inv_cnt = 1.0 / cnt
        mx *= inv_cnt
        my *= inv_cnt
        mz *= inv_cnt
        norm_m = sqrt(mx^2 + my^2 + mz^2)
        norm_m == 0 && continue

        push!(out_labels, lbl)
        push!(out_incs, acos(clamp(mz / norm_m, -1.0, 1.0)) * (180 / pi))
        push!(out_azis, atan(my, mx) * (180 / pi))
    end

    if isempty(out_labels)
        new_df = df[1:0, :]
    else
        new_df = DataFrame()
        :label in propertynames(df) && (new_df[!, :label] = out_labels)
        :inc in propertynames(df) && (new_df[!, :inc] = out_incs)
        :azi in propertynames(df) && (new_df[!, :azi] = out_azis)
    end
    new_layout = Layout(new_df, nothing, nothing, nothing)
    _ensure_all_coordinates!(new_layout)
    return new_layout
end

function _ensure_all_coordinates!(layout::Layout)
    _ensure_coordinates_3d!(layout)
    _ensure_coordinates_2d!(layout)
    return layout
end

function _append_layouts(base_layout::Layout, avg_layout::Layout)::Layout
    base_copy = copy(base_layout)
    _ensure_all_coordinates!(base_copy)
    if isempty(avg_layout.data)
        return Layout(base_copy.data, nothing, nothing, nothing)
    end
    if !hasproperty(avg_layout.data, :x3) || !hasproperty(avg_layout.data, :x2)
        avg_copy = copy(avg_layout)
        _ensure_all_coordinates!(avg_copy)
        return Layout(vcat(base_copy.data, avg_copy.data), nothing, nothing, nothing)
    end
    return Layout(vcat(base_copy.data, avg_layout.data), nothing, nothing, nothing)
end
