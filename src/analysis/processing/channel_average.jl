"""
    channel_average!(dat::EegData; channel_selections::AbstractVector{<:Function} = [channels()],
                     output_labels = nothing,
                     include_extra::Bool = false,
                     reduce::Bool = false) -> Nothing

Create averaged channels for each provided channel-selection predicate and add them to `dat`.

- `channel_selections` is a vector of channel predicates (see `channels(...)`) that defaults to `[channels()]` (all channels).
- For every predicate, the selected channels are averaged to produce a new column.
- Output labels default to the underscore-joined contributing channel labels (e.g., :Fp1_Fp2).
- Can be used for single or multiple channel selections by passing a vector with one or more functions.

Layout handling
- `reduce = false` (default): append averaged columns to data and averaged channels to layout.
- `reduce = true`: replace data with metadata + averaged columns and create new layout with only averaged channels.

Examples
```julia
# Average all channels (default)
channel_average!(dat)

# Single channel selection
channel_average!(dat, channel_selections = [channels([:Fp1, :Fp2])])

# Multiple channel selections (append columns and layout)
channel_average!(dat, channel_selections = [channels([:Fp1, :Fp2]), channels([:O1, :O2])])

# Average by patterns with custom labels
channel_average!(dat,
    channel_selections = [x -> startswith.(string.(x), "F"), x -> startswith.(string.(x), "P")];
    output_labels = [:F_avg, :P_avg],
)

# Produce a reduced dataset (meta + averages only) with new layout
channel_average!(dat, channel_selections = [channels([:Fp1, :Fp2])]; reduce = true)
```
"""
function channel_average!(
    dat::EegData;
    channel_selections::AbstractVector = [channels()],
    output_labels = nothing,
    include_extra::Bool = false,
    reduce::Bool = false,
)::Nothing

    # Validate that all channel_selections are callable
    for (i, sel) in enumerate(channel_selections)
        if !(sel isa Function)
            @minimal_error_throw "Channel selection at index $i must be a Function, got $(typeof(sel))"
        end
    end

    # Resolve actual channel groups (exclude metadata; optionally include extra)
    # For the default case (channels()), we want to select only original channel columns, not derived averaged columns
    selected_channel_groups::Vector{Vector{Symbol}} = if length(channel_selections) == 1 && channel_selections[1] == channels()
        # Default case: select only original channel columns
        all_cols = all_labels(dat)
        meta_cols = meta_labels(dat)
        original_channels = [col for col in all_cols if col ∉ meta_cols && !contains(string(col), "_")]
        [original_channels]
    else
        # Custom selections: use the provided functions
        [get_selected_channels(dat, sel; include_meta = false, include_extra = include_extra) for sel in channel_selections]
    end

    # Validate channel groups
    if any(isempty, selected_channel_groups)
        empty_selections = findall(isempty, selected_channel_groups)
        @minimal_error_throw "Channel selections at indices $empty_selections produced no channels"
    end

    # Determine labels (explicit if/else for clarity)
    if output_labels === nothing
        # Get original channel columns (excluding any previously added averaged columns)
        # We need to identify which columns are the original channels vs derived averaged columns
        all_cols = all_labels(dat)
        meta_cols = meta_labels(dat)

        # Original channels are those that are not metadata and don't contain underscores (averaged columns)
        original_channels = [col for col in all_cols if col ∉ meta_cols && !contains(string(col), "_")]

        labels = Vector{Symbol}(undef, length(selected_channel_groups))
        for (i, grp) in enumerate(selected_channel_groups)
            # Check if this group covers all original channel columns (excluding metadata)
            is_all = length(grp) == length(original_channels) && all(in(original_channels), grp)
            labels[i] = is_all ? :avg : Symbol(join(string.(grp), "_"))
        end
    else
        # Accept Vector of Symbol/String; coerce to Symbol
        if length(output_labels) != length(selected_channel_groups)
            @minimal_error_throw "N output_labels ($(length(output_labels))) must match N channel selections ($(length(selected_channel_groups)))"
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

    # Apply for each channel selection (append columns)
    for (grp, lbl) in zip(selected_channel_groups, labels)
        if lbl ∈ all_labels(dat)
            @minimal_warning "Overwriting existing channel '$(lbl)'"
        end
        @info "channel_average!: $(print_vector(grp)) → :$(lbl)"
        _add_avg!(dat.data, grp, lbl)
    end

    # Handle data and layout based on reduce parameter
    if reduce
        # Replace data with metadata + averaged columns only
        dat.data = _build_reduced_df(dat, selected_channel_groups, labels)
        # Create new layout with only averaged channels
        dat.layout = _layout_from_groups(dat.layout, selected_channel_groups, labels)
    else
        # Append averaged channels to existing layout
        avg_layout = _layout_from_groups(dat.layout, selected_channel_groups, labels)
        dat.layout = _append_layouts(dat.layout, avg_layout)
    end

    return nothing
end

@add_nonmutating channel_average!

# Handle collections of EEG data by broadcasting
function channel_average(
    dat::Vector{<:EegData};
    channel_selections::AbstractVector = [channels()],
    output_labels = nothing,
    include_extra::Bool = false,
    reduce::Bool = false,
)
    return channel_average.(
        dat;
        channel_selections = channel_selections,
        output_labels = output_labels,
        include_extra = include_extra,
        reduce = reduce,
    )
end

#=============================================================================
    BATCH PROCESSING (FILE PATTERN)
=============================================================================#

"""Generate default output directory name for channel averaging."""
function _default_channel_average_output_dir(input_dir::String, pattern::String)
    joinpath(input_dir, "averaged_channels_$(pattern)")
end

"""
Process a single file through channel averaging pipeline.
Returns BatchResult with success/failure info.
"""
function _process_channel_average_file(
    filepath::String,
    output_path::String,
    channel_selections::Vector{<:Function},
    output_labels::Vector{Symbol},
    condition_selection::Function,
    reduce::Bool,
)
    filename = basename(filepath)

    # Load data
    data = read_data(filepath)
    if isnothing(data)
        return BatchResult(false, filename, "No data variables found")
    end

    # Validate that data is valid EEG data (Vector of ErpData or EpochData)
    if !(data isa Vector{<:Union{ErpData,EpochData}})
        return BatchResult(false, filename, "Invalid data type: expected Vector{ErpData} or Vector{EpochData}")
    end

    # Select conditions
    data = _condition_select(data, condition_selection)

    # Apply channel averaging to each data item
    foreach(data) do item
        channel_average!(item, channel_selections = channel_selections, output_labels = output_labels, reduce = reduce)
    end

    # Save (always use "data" as variable name since load_data finds by type)
    jldsave(output_path; data = data)

    n_groups = length(channel_selections)
    return BatchResult(true, filename, "Averaged $n_groups channel group(s)")
end

"""
    channel_average(file_pattern::String, channel_selections::Vector{Function}; 
                    output_labels::Union{Vector{Symbol}, Nothing} = nothing,
                    input_dir::String = pwd(), 
                    participant_selection::Function = participants(),
                    condition_selection::Function = conditions(),
                    output_dir::Union{String, Nothing} = nothing,
                    reduce::Bool = false)

Batch process EEG data files to average specified channels using predicates.

This function loads JLD2 files containing EEG data, averages specified channel groups,
and saves the resulting data with new averaged channels to a new directory.

# Arguments
- `file_pattern::String`: Pattern to match JLD2 files (e.g., "erps_cleaned", "epochs_cleaned")
- `channel_selections::Vector{Function}`: Channel selection predicates (e.g., `[channels([:Fp1, :Fp2]), channels([:PO7, :PO8])]`)
- `output_labels::Union{Vector{Symbol}, Nothing}`: Labels for averaged channels (default: auto-generated from selection)
- `input_dir::String`: Input directory containing JLD2 files (default: current directory)
- `participant_selection::Function`: Participant selection predicate (default: `participants()` for all)
- `condition_selection::Function`: Condition selection predicate (default: `conditions()` for all)
- `output_dir::Union{String, Nothing}`: Output directory (default: auto-generated)
- `reduce::Bool`: Whether to keep only averaged channels (true) or append to existing (false, default)

# Examples
```julia
# Average frontal and parietal channels using predicates
channel_average("erps_cleaned", [channels([:Fp1, :Fp2]), channels([:PO7, :PO8])])

# With custom output labels
channel_average("erps_cleaned", 
                [channels([:Fp1, :Fp2]), channels([:PO7, :PO8])],
                output_labels = [:frontal, :parietal])

# Process specific participants and conditions
channel_average("erps_cleaned", [channels([:Fp1, :Fp2])], 
                input_dir = "/path/to/data", 
                participants = [1, 2, 3], 
                conditions = [1, 2])

# Create reduced dataset with only averaged channels
channel_average("erps_cleaned", 
                [channels([:Fp1, :Fp2]), channels([:PO7, :PO8])], 
                reduce = true)
```
"""
function channel_average(
    file_pattern::String,
    channel_selections::Vector{<:Function};
    output_labels::Union{Vector{Symbol},Nothing} = nothing,
    input_dir::String = pwd(),
    participant_selection::Function = participants(),
    condition_selection::Function = conditions(),
    output_dir::Union{String,Nothing} = nothing,
    reduce::Bool = false,
)

    # Setup logging
    log_file = "channel_average.log"
    setup_global_logging(log_file)

    try
        @info "Batch channel averaging started at $(now())"
        @log_call "channel_average"

        # Validation (early return on error)
        if (error_msg = _validate_input_dir(input_dir)) !== nothing
            @minimal_error_throw(error_msg)
        end

        # Generate output labels if not provided
        if isnothing(output_labels)
            output_labels = [Symbol("avg_$i") for i = 1:length(channel_selections)]
        end

        # Validate output labels length matches channel selections
        if length(output_labels) != length(channel_selections)
            @minimal_error_throw(
                "Number of output_labels ($(length(output_labels))) must match number of channel_selections ($(length(channel_selections)))"
            )
        end

        # Setup directories
        output_dir = something(output_dir, _default_channel_average_output_dir(input_dir, file_pattern))
        mkpath(output_dir)

        # Find files
        files = _find_batch_files(file_pattern, input_dir, participant_selection)

        if isempty(files)
            @minimal_warning "No JLD2 files found matching pattern '$file_pattern' in $input_dir"
            return nothing
        end

        @info "Found $(length(files)) JLD2 files matching pattern '$file_pattern'"
        @info "Channel selections: $(length(channel_selections)) group(s)"
        @info "Output labels: $output_labels"
        @info "Reduce mode: $reduce"

        # Create processing function with captured parameters
        process_fn =
            (input_path, output_path) ->
                _process_channel_average_file(input_path, output_path, channel_selections, output_labels, condition_selection, reduce)

        # Execute batch operation
        results = _run_batch_operation(process_fn, files, input_dir, output_dir; operation_name = "Averaging channels")

        return _log_batch_summary(results, output_dir)

    finally
        _cleanup_logging(log_file, output_dir)
    end
end


# Internal: create a reduced DataFrame(s) with only meta + averaged columns
function _build_reduced_df(dat::SingleDataFrameEeg, channel_groups::Vector{Vector{Symbol}}, labels::Vector{Symbol})
    meta_cols = meta_labels(dat)
    new_df = isempty(meta_cols) ? DataFrame() : dat.data[:, meta_cols]
    for (grp, lbl) in zip(channel_groups, labels)
        new_df[!, lbl] = colmeans(dat.data, grp)
    end
    return new_df
end

function _build_reduced_df(dat::MultiDataFrameEeg, channel_groups::Vector{Vector{Symbol}}, labels::Vector{Symbol})
    meta_cols = meta_labels(dat)
    new_epochs = Vector{DataFrame}(undef, length(dat.data))
    for (i, df) in pairs(dat.data)
        new_df = isempty(meta_cols) ? DataFrame() : df[:, meta_cols]
        for (grp, lbl) in zip(channel_groups, labels)
            new_df[!, lbl] = colmeans(df, grp)
        end
        new_epochs[i] = new_df
    end
    return new_epochs
end



# Internal: build a layout with one row per averaged group
function _layout_from_groups(layout::Layout, channel_groups::Vector{Vector{Symbol}}, labels::Vector{Symbol})::Layout
    df = layout.data
    isempty(df) && return Layout(DataFrame(), nothing, nothing)

    # Ensure coordinates once
    work_layout = copy(layout)
    _ensure_all_coordinates!(work_layout)
    wdf = work_layout.data
    labels_in_layout = df[!, :label]

    # Build new rows
    new_rows = DataFrame[]
    for (grp, lbl) in zip(channel_groups, labels)
        sel_idx = findall(in(grp), labels_in_layout)
        isempty(sel_idx) && continue

        # Calculate averaged coordinates
        coords = _average_coordinates(wdf[sel_idx, :])
        isnothing(coords) && continue

        # Create row with averaged coordinates
        row = _create_layout_row(df, lbl, coords)
        push!(new_rows, row)
    end

    # Build new layout and ensure coordinates
    new_df = isempty(new_rows) ? df[1:0, :] : reduce(vcat, new_rows)
    new_layout = Layout(new_df, nothing, nothing)
    _ensure_all_coordinates!(new_layout)
    return new_layout
end

# Helper: ensure both 3D and 2D coordinates
function _ensure_all_coordinates!(layout::Layout)
    _ensure_coordinates_3d!(layout)
    _ensure_coordinates_2d!(layout)
    return nothing
end

# Helper: average 3D coordinates and convert to polar
function _average_coordinates(coord_df)

    # Average 3D Cartesian coordinates first 
    mx, my, mz = mean(Float64.(coord_df[!, :x3])), mean(Float64.(coord_df[!, :y3])), mean(Float64.(coord_df[!, :z3]))
    norm_m = sqrt(mx^2 + my^2 + mz^2)
    norm_m == 0 && return nothing

    inc_deg = acos(clamp(mz / norm_m, -1.0, 1.0)) * (180 / pi)
    azi_deg = atan(my, mx) * (180 / pi)

    return (inc = inc_deg, azi = azi_deg)
end

# Helper: create a layout row with averaged coordinates
function _create_layout_row(template_df, label, coords)
    new_df = DataFrame()
    if :label in propertynames(template_df)
        new_df[!, :label] = [label]
    end
    if :inc in propertynames(template_df)
        new_df[!, :inc] = [coords.inc]
    end
    if :azi in propertynames(template_df)
        new_df[!, :azi] = [coords.azi]
    end
    return new_df
end

# Internal: append averaged channels to existing layout
function _append_layouts(base_layout::Layout, avg_layout::Layout)::Layout
    # Ensure both layouts have coordinates
    base_copy = copy(base_layout)
    _ensure_all_coordinates!(base_copy)
    avg_copy = copy(avg_layout)
    _ensure_all_coordinates!(avg_copy)

    # Simple append - no complex merging needed
    combined_df = vcat(base_copy.data, avg_copy.data)
    return Layout(combined_df, nothing, nothing)
end

