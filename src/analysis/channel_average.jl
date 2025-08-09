"""
    channel_average!(dat::EegData, group_selections::AbstractVector{<:Function};
                     output_labels = nothing,
                     include_extra::Bool = false,
                     reduce::Bool = false,
                     update_layout::Bool = false) -> Nothing

Create averaged channels for each provided channel-selection predicate and add them to `dat`.

- Each function in `group_selections` is a channel predicate (see `channels(...)`).
- For every predicate, the selected channels are averaged sample-wise to produce a new column.
- Output labels default to the underscore-joined contributing channel labels (e.g., :Fp1_Fp2).

Layout handling
- `reduce = false, update_layout = false` (default): append averaged columns only; layout unchanged.
- `reduce = true`: replace data with metadata + averaged columns and build a matching averaged layout (positions computed from averaged 3D → polar → regenerated 2D/3D).
- `update_layout = true`: keep original data but merge averaged rows into the layout (replace existing labels if present; recompute 2D/3D via helpers).

Examples
```julia
# Average specific pairs (append columns only)
channel_average!(dat, [channels([:Fp1, :Fp2]), channels([:O1, :O2])])

# Average by patterns with custom labels and update layout
channel_average!(dat,
    [x -> startswith.(string.(x), "F"), x -> startswith.(string.(x), "P")];
    output_labels = [:F_avg, :P_avg],
    update_layout = true,
)

# Produce a reduced dataset (meta + averages only) with averaged layout
channel_average!(dat, [channels([:Fp1, :Fp2])]; reduce = true)
```
"""
function channel_average!(
    dat::EegData,
    group_selections::AbstractVector{<:Function};
    output_labels = nothing,
    include_extra::Bool = false,
    reduce::Bool = false,
    update_layout::Bool = false,
)::Nothing

    # Resolve actual channel groups (exclude metadata; optionally include extra)
    selected_groups = [
        get_selected_channels(dat, sel; include_meta = false, include_extra = include_extra)
        for sel in group_selections
    ]

    # Validate groups
    any(isempty, selected_groups) && @minimal_error_throw "One or more group selections produced no channels"

    # Determine labels (explicit if/else for clarity)
    if output_labels === nothing
        all_ch = channel_labels(dat)
        labels = Vector{Symbol}(undef, length(selected_groups))
        @inbounds for (i, grp) in enumerate(selected_groups)
            is_all = length(grp) == length(all_ch) && all(in(all_ch), grp)
            labels[i] = is_all ? :avg : Symbol(join(string.(grp), "_"))
        end
    else
        # Accept Vector of Symbol/String; coerce to Symbol
        if length(output_labels) != length(selected_groups)
            @minimal_error_throw "Length of output_labels ($(length(output_labels))) must match number of groups ($(length(selected_groups)))"
        end
        labels = Symbol.(output_labels)
    end

    # Internal helpers (single DataFrame vs vector of DataFrames)
    function _add_avg!(df::DataFrame, cols::Vector{Symbol}, outlbl::Symbol)
        df[!, outlbl] = colmeans(df, cols)
        return nothing
    end
    function _add_avg!(dfs::Vector{DataFrame}, cols::Vector{Symbol}, outlbl::Symbol)
        _add_avg!.(dfs, Ref(cols), Ref(outlbl))
        return nothing
    end

    # Apply for each group (append columns)
    for (grp, lbl) in zip(selected_groups, labels)
        if lbl ∈ all_labels(dat)
            @minimal_warning "Overwriting existing channel '$(lbl)'"
        end
        @info "channel_average!: $(join(string.(grp), ", ")) -> $(lbl)"
        _add_avg!(dat.data, grp, lbl)
    end

    # If reducing to only averages, rebuild data and layout to match
    if reduce
        # Build reduced data: meta + averaged columns only
        dat.data = _build_reduced_df(dat, selected_groups, labels)
        # Always produce a matching layout for the averaged channels
        dat.layout = _layout_from_groups(dat.layout, selected_groups, labels)
    elseif update_layout
        # Merge averaged rows into existing layout and ensure coordinates
        avg_layout = _layout_from_groups(dat.layout, selected_groups, labels)
        dat.layout = _merge_layout_with_averages(dat.layout, avg_layout)
    end

    return nothing
end

"""
    channel_average!(dat::EegData, groups;
                     output_labels = nothing,
                     include_extra::Bool = false,
                     reduce::Bool = false,
                     update_layout::Bool = false) -> Nothing

Convenience wrapper to average explicit channel groups (e.g., pairs). `groups` can be vectors of Symbols or Strings.

Examples
```julia
# [["a","b"], ["c","d"]] → a_b and c_d (append only)
channel_average!(dat, [["a","b"], ["c","d"]])

# Reduce with custom labels
channel_average!(dat, [[:Fp1,:Fp2]]; output_labels = [:F_avg], reduce = true)
```
"""
function channel_average!(
    dat::EegData,
    groups::AbstractVector{<:AbstractVector};
    output_labels = nothing,
    include_extra::Bool = false,
    reduce::Bool = false,
    update_layout::Bool = false,
)::Nothing
    # Coerce group labels to Symbols once
    sym_groups = [Symbol.(g) for g in groups]
    preds = [channels(g) for g in sym_groups]
    return channel_average!(
        dat,
        preds;
        output_labels = output_labels,
        include_extra = include_extra,
        reduce = reduce,
        update_layout = update_layout,
    )
end

@add_nonmutating channel_average!

# Internal: create a reduced DataFrame(s) with only meta + averaged columns
function _build_reduced_df(dat::SingleDataFrameEeg, groups::Vector{Vector{Symbol}}, labels::Vector{Symbol})
    meta_cols = meta_labels(dat)
    new_df = isempty(meta_cols) ? DataFrame() : dat.data[:, meta_cols]
    for (grp, lbl) in zip(groups, labels)
        new_df[!, lbl] = colmeans(dat.data, grp)
    end
    return new_df
end
function _build_reduced_df(dat::MultiDataFrameEeg, groups::Vector{Vector{Symbol}}, labels::Vector{Symbol})
    meta_cols = meta_labels(dat)
    new_epochs = Vector{DataFrame}(undef, length(dat.data))
    for (i, df) in pairs(dat.data)
        new_df = isempty(meta_cols) ? DataFrame() : df[:, meta_cols]
        for (grp, lbl) in zip(groups, labels)
            new_df[!, lbl] = colmeans(df, grp)
        end
        new_epochs[i] = new_df
    end
    return new_epochs
end

# Internal: build a layout with one row per averaged group
function _layout_from_groups(layout::Layout, groups::Vector{Vector{Symbol}}, labels::Vector{Symbol})::Layout
    df = layout.data
    if isempty(df)
        return Layout(DataFrame(), nothing, nothing)
    end
    # Ensure we have coordinates via existing helpers
    work_layout = copy(layout)
    _ensure_coordinates_3d!(work_layout)
    _ensure_coordinates_2d!(work_layout)
    wdf = work_layout.data

    # Prepare output rows (preserve schema from original df)
    new_df = df[1:0, :]
    labels_in_layout = df[!, :label]

    for (grp, lbl) in zip(groups, labels)
        sel_idx = findall(in(grp), labels_in_layout)
        isempty(sel_idx) && continue

        # Mean 3D vector from ensured coords
        mx = mean(Float64.(wdf[sel_idx, :x3]))
        my = mean(Float64.(wdf[sel_idx, :y3]))
        mz = mean(Float64.(wdf[sel_idx, :z3]))
        norm_m = sqrt(mx^2 + my^2 + mz^2)
        norm_m == 0 && continue

        # Convert to polar (degrees)
        inc_rad = acos(clamp(mz / norm_m, -1.0, 1.0))
        azi_rad = atan(my, mx)
        inc_deg = inc_rad * (180 / pi)
        azi_deg = azi_rad * (180 / pi)

        # Build row with label and polar
        row = df[1:1, :]
        row[1, :label] = lbl
        if :inc in names(row); row[1, :inc] = inc_deg; end
        if :azi in names(row); row[1, :azi] = azi_deg; end
        append!(new_df, row)
    end

    # Derive coordinates via existing helpers
    new_layout = Layout(new_df, nothing, nothing)
    _ensure_coordinates_3d!(new_layout)
    _ensure_coordinates_2d!(new_layout)
    return new_layout
end

# Internal: merge averaged rows into an existing layout (replace if label exists), ensure coords
function _merge_layout_with_averages(layout::Layout, avg_layout::Layout)::Layout
    # Ensure both sides have full coordinate columns to align schemas
    base_layout = copy(layout)
    _ensure_coordinates_3d!(base_layout)
    _ensure_coordinates_2d!(base_layout)
    base_df = base_layout.data

    add_layout = copy(avg_layout)
    _ensure_coordinates_3d!(add_layout)
    _ensure_coordinates_2d!(add_layout)
    add_df = add_layout.data

    if isempty(add_df)
        return base_layout
    end

    common_cols = intersect(names(base_df), names(add_df))

    for i in 1:nrow(add_df)
        lbl = add_df[i, :label]
        idx = findfirst(==(lbl), base_df.label)
        if isnothing(idx)
            # Create a one-row DataFrame with base schema and fill from add_df where available
            template = base_df[1:1, :]
            template[1, :label] = lbl
            for col in common_cols
                template[1, col] = add_df[i, col]
            end
            append!(base_df, template; cols = :setequal)
        else
            # Replace available columns on existing row
            for col in common_cols
                base_df[idx, col] = add_df[i, col]
            end
        end
    end

    new_layout = Layout(base_df, nothing, nothing)
    # Ensure consistency
    _ensure_coordinates_3d!(new_layout)
    _ensure_coordinates_2d!(new_layout)
    return new_layout
end


