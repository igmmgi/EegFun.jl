"""
Channel averaging for time-frequency data.

Provides TF-specific `channel_average!` and `channel_average` that apply
averaging to both `data_power` and `data_phase` DataFrames.
"""

"""
    channel_average!(dat::TimeFreqData; channel_selections::AbstractVector = [channels()],
                     output_labels = nothing, reduce::Bool = false) -> Nothing

Create averaged channels for TimeFreqData, applying to both power and phase DataFrames.

See `channel_average!` for full documentation. This method ensures both `data_power`
and `data_phase` are updated consistently.

# Examples
```julia
# Average all channels
channel_average!(tf_data)

# Average specific channel groups
channel_average!(tf_data, channel_selections = [channels([:Ch1, :Ch2])], output_labels = [:frontal])

# Reduce to only averaged channels
channel_average!(tf_data, channel_selections = [channels([:Ch1, :Ch2])]; reduce = true)
```
"""
function channel_average!(
    dat::TimeFreqData;
    channel_selections::AbstractVector = [channels()],
    output_labels = nothing,
    reduce::Bool = false,
)::Nothing

    # Validate that all channel_selections are callable
    for (i, sel) in enumerate(channel_selections)
        if !(sel isa Function)
            @minimal_error_throw "Channel selection at index $i must be a Function, got $(typeof(sel))"
        end
    end

    # Resolve actual channel groups
    selected_channel_groups::Vector{Vector{Symbol}} = if length(channel_selections) == 1 && channel_selections[1] == channels()
        all_cols = all_labels(dat)
        meta_cols = meta_labels(dat)
        original_channels = [col for col in all_cols if col ∉ meta_cols && !contains(string(col), "_")]
        [original_channels]
    else
        [get_selected_channels(dat, sel; include_meta = false, include_extra = false) for sel in channel_selections]
    end

    # Validate channel groups
    if any(isempty, selected_channel_groups)
        empty_selections = findall(isempty, selected_channel_groups)
        @minimal_error_throw "Channel selections at indices $empty_selections produced no channels"
    end

    # Determine labels
    if isnothing(output_labels)
        all_cols = all_labels(dat)
        meta_cols = meta_labels(dat)
        original_channels = [col for col in all_cols if col ∉ meta_cols && !contains(string(col), "_")]

        labels = Vector{Symbol}(undef, length(selected_channel_groups))
        for (i, grp) in enumerate(selected_channel_groups)
            is_all = length(grp) == length(original_channels) && all(in(original_channels), grp)
            labels[i] = is_all ? :avg : Symbol(join(string.(grp), "_"))
        end
    else
        if length(output_labels) != length(selected_channel_groups)
            @minimal_error_throw "N output_labels ($(length(output_labels))) must match N channel selections ($(length(selected_channel_groups)))"
        end
        labels = Symbol.(output_labels)
    end

    # Apply averaging to BOTH power and phase DataFrames
    for (grp, lbl) in zip(selected_channel_groups, labels)
        @info "channel_average! (TF): $(_print_vector(grp)) → :$(lbl)"
        dat.data_power[!, lbl] = _colmeans(dat.data_power, grp)
        dat.data_phase[!, lbl] = _colmeans(dat.data_phase, grp)
    end

    # Handle layout
    if reduce
        # Build reduced DataFrames (meta + averaged columns only)
        meta_cols = meta_labels(dat)

        new_power = isempty(meta_cols) ? DataFrame() : dat.data_power[:, meta_cols]
        new_phase = isempty(meta_cols) ? DataFrame() : dat.data_phase[:, meta_cols]
        for (grp, lbl) in zip(selected_channel_groups, labels)
            new_power[!, lbl] = _colmeans(dat.data_power, grp)
            new_phase[!, lbl] = _colmeans(dat.data_phase, grp)
        end
        dat.data_power = new_power
        dat.data_phase = new_phase

        dat.layout = _layout_from_groups(dat.layout, selected_channel_groups, labels)
    else
        avg_layout = _layout_from_groups(dat.layout, selected_channel_groups, labels)
        dat.layout = _append_layouts(dat.layout, avg_layout)
    end

    return nothing
end


"""
    channel_average(dat::TimeFreqData; kwargs...) -> TimeFreqData

Non-mutating version of `channel_average!` for TimeFreqData.
Creates a copy and applies the operation.

# Examples
```julia
# Average all channels (non-mutating)
result = EegFun.channel_average(tf_data)

# Average specific channel groups
result = EegFun.channel_average(tf_data,
    channel_selections = [EegFun.channels([:Fz, :Cz, :Pz])],
    output_labels = [:midline])
```
"""
function channel_average(dat::TimeFreqData; kwargs...)
    data_copy = copy(dat)
    channel_average!(data_copy; kwargs...)
    return data_copy
end

# Handle Vector{TimeFreqData}
function channel_average(dat::Vector{TimeFreqData}; kwargs...)
    return [channel_average(d; kwargs...) for d in dat]
end
