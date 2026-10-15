function _default_subset_output_dir(input_dir::String, pattern::String)
    joinpath(input_dir, "subset_$(clean_pattern(pattern))")
end

function _process_subset_file(filepath::String, output_path::String; kwargs...)
    filename = basename(filepath)
    dat = read_data(filepath)

    if isnothing(dat)
        return BatchResult(false, filename, "No data found in file")
    end

    try
        subset_data = subset(dat; kwargs...)
        jldsave(output_path; data = subset_data)

        # Determine stats for logging
        if subset_data isa Vector && !isempty(subset_data) && hasproperty(subset_data[1], :data)
            n_items = sum(length(cond.data) for cond in subset_data)
            return BatchResult(true, filename, "Processed subset. Elements remaining: $(n_items)")
        elseif hasproperty(subset_data, :data) && subset_data.data isa Vector
            n_items = length(subset_data.data)
            return BatchResult(true, filename, "Processed subset. Elements remaining: $(n_items)")
        else
            return BatchResult(true, filename, "Processed subset.")
        end
    catch e
        return BatchResult(false, filename, "Error: $(sprint(showerror, e))")
    end
end

"""
    subset(file_pattern::String; kwargs...)

Batch process JLD2 files that match `file_pattern`, applying the `subset` operation.
Supports all `subset` kwargs (`epoch_selection`, `channel_selection`, `sample_selection`, `interval_selection`).
"""
function subset(
    file_pattern::String;
    input_dir::String = pwd(),
    participant_selection::Function = participants(),
    output_dir::Union{String,Nothing} = nothing,
    kwargs...,
)
    log_file = "subset.log"
    setup_global_logging(log_file)

    try
        @info "Batch subset started at $(now())"

        error_msg = _validate_input_dir(input_dir)
        if !isnothing(error_msg)
            @minimal_error(error_msg)
        end

        output_dir = something(output_dir, _default_subset_output_dir(input_dir, file_pattern))
        mkpath(output_dir)

        process_fn = (input_path, output_path) -> _process_subset_file(input_path, output_path; kwargs...)

        batch_process(process_fn, file_pattern, input_dir, output_dir, participant_selection, "Subsetting parameters")

    finally
        _cleanup_logging(log_file, output_dir)
    end
end




"""
    log_pretty_table(df::DataFrame; log_level::Symbol = :info, kwargs...)

Log a pretty table with specified log level. For general DataFrame logging.
Sets show_row_number=false and show_subheader=false by default for cleaner logs.

# Arguments
- `df::DataFrame`: The DataFrame to log
- `log_level::Symbol`: Log level (:debug, :info, :warn, :error) (default: :info)
- `kwargs...`: Additional arguments passed to pretty_table

# Examples
```julia
# Log with default info level
log_pretty_table(df; title = "My Table")

# Log with debug level
log_pretty_table(df; log_level = :debug, title = "Debug Table")

# Log with warn level
log_pretty_table(df; log_level = :warn, title = "Warning Table")
```
"""
function log_pretty_table(df::DataFrame; log_level::Symbol = :info, kwargs...)

    table_output = sprint() do output_io
        # Set a large display size to avoid terminal limitations (i.e., cropping)
        # when capturing the table string output.
        io_context = IOContext(output_io, :displaysize => (2000, 2000), :color => true)
        # Default column_labels to column names only — this suppresses the auto-generated
        # type subheader row that PrettyTables adds when given a DataFrame directly.
        # Callers can still override by passing their own column_labels kwarg.
        kw = Dict{Symbol,Any}(kwargs)
        get!(kw, :column_labels, names(df))
        pretty_table(io_context, df; kw...)
    end

    # Map symbol to standard LogLevel
    levels = Dict(:debug => Logging.Debug, :info => Logging.Info, :warn => Logging.Warn, :error => Logging.Error)

    lvl = get(levels, log_level, Logging.Info)
    if !haskey(levels, log_level)
        @minimal_warning "Unknown log level: $log_level, using :info instead"
    end

    # Log with minimal metadata (no file/line noise)
    @logmsg lvl "\n\n$table_output\n" _module=nothing _file=nothing _line=nothing

    return nothing
end


"""
    channels()
    channels(name::Symbol)
    channels(names::Symbol...)
    channels(names::Vector{Symbol})
    channels(index::Int)
    channels(indices::Union{Vector{Int}, UnitRange})
    channels(ranges::Vector{UnitRange{Int}})
    channels(predicate::Function)

Predicate generators for channel selection in subsetting and plotting operations.

Returns a function that, given a vector of channel labels, produces a boolean mask
selecting the matching channels. Used with the `channel_selection` keyword argument.

# Examples
```julia
channels()                          # Select all channels (default)
channels(:Cz)                       # Select a single channel by name
channels(:Fz, :Cz, :Pz)            # Select multiple channels by name
channels([:Fz, :Cz, :Pz])          # Same, from a Vector
channels(1:10)                      # Select first 10 channels by index
channels(x -> startswith.(string.(x), "F"))  # Custom predicate
```
"""
channels() = x -> fill(true, length(x))  # Default: select all channels given
channels(channel_names::Vector{Symbol}) = x -> x .∈ Ref(channel_names)
channels(channel_name::Symbol) = channels([channel_name])
channels(channel_names::Symbol...) = channels(collect(channel_names))
channels(channel_number::Int) = channels([channel_number])
channels(channel_numbers::Union{Vector{Int},UnitRange}) = x -> [i in channel_numbers for i = 1:length(x)]
channels(channel_ranges::Vector{UnitRange{Int}}) = x -> [i in union(channel_ranges...) for i = 1:length(x)]
channels(predicate::Function) = predicate  # Allow custom function predicates
"""
    channels_not(name::Symbol)
    channels_not(names::Symbol...)
    channels_not(names::Vector{Symbol})
    channels_not(indices::Union{Vector{Int}, UnitRange})

Predicate generators that **exclude** the specified channels.
Inverse of `channels(...)`.

# Examples
```julia
channels_not(:Fp1, :Fp2)           # Exclude two frontal channels
channels_not(1:4)                   # Exclude first 4 channels by index
```
"""
channels_not(channel_names::Vector{Symbol}) = x -> .!(x .∈ Ref(channel_names))
channels_not(channel_name::Symbol) = channels_not([channel_name])
channels_not(channel_names::Symbol...) = channels_not(collect(channel_names))
channels_not(channel_numbers::Union{Vector{Int},UnitRange}) = x -> .!([i in channel_numbers for i = 1:length(x)])

"""
    components()
    components(index::Int)
    components(indices::Int...)
    components(indices::Union{Vector{Int}, UnitRange})
    components(predicate::Function)

Predicate generators for ICA component selection.

# Examples
```julia
components()                        # Select all components (default)
components(1, 2, 3)                 # Select components 1, 2, and 3
components(1:5)                     # Select first 5 components
```
"""
components() = x -> fill(true, length(x))  # Default: select all components given
components(component_numbers::Union{Vector{Int},UnitRange}) = x -> [i in component_numbers for i = 1:length(x)]
components(component_number::Int) = x -> x .== component_number
components(component_numbers::Int...) = components(collect(component_numbers))
components(predicate::Function) = predicate  # Allow custom function predicates
"""
    components_not(index::Int)
    components_not(indices::Int...)
    components_not(indices::Union{Vector{Int}, UnitRange})

Predicate generators that **exclude** the specified ICA components.

# Examples
```julia
components_not(1)                   # Exclude component 1
components_not(1:3)                 # Exclude first 3 components
```
"""
components_not(component_numbers::Union{Vector{Int},UnitRange}) = x -> .!([i in component_numbers for i = 1:length(x)])
components_not(component_number::Int) = x -> .!(x .== component_number)
components_not(component_numbers::Int...) = components_not(collect(component_numbers))

"""
    samples()
    samples(column::Symbol)
    samples(predicate::Function)

Predicate generators for sample (row) selection. Returns a function that,
given a DataFrame, produces a boolean mask over rows.

# Examples
```julia
samples()                           # Select all samples (default)
samples(:is_clean)                  # Select rows where :is_clean is true
samples(df -> df.time .> 0)         # Custom predicate on the DataFrame
```
"""
samples() = x -> fill(true, nrow(x))
samples(column::Symbol) = x -> x[!, column]
# Allow custom function predicates
samples(predicate::Function) = predicate
"""
    samples_or(columns::Vector{Symbol})
    samples_and(columns::Vector{Symbol})
    samples_not(column::Symbol)
    samples_or_not(columns::Vector{Symbol})
    samples_and_not(columns::Vector{Symbol})

Combine or negate sample predicates across boolean columns.

# Examples
```julia
samples_or([:is_extreme, :is_drift])    # Rows where ANY flag is true
samples_and([:is_clean, :is_valid])     # Rows where ALL flags are true
samples_not(:is_bad)                    # Rows where :is_bad is false
samples_or_not([:is_extreme, :is_bad])  # Rows where NO flag is true
```
"""
samples_or(columns::Vector{Symbol}) =
    x -> isempty(columns) ? fill(false, nrow(x)) : reduce((a, b) -> a .| b, (x[!, col] for col in columns))
"""See `samples_or`."""
samples_and(columns::Vector{Symbol}) =
    x -> isempty(columns) ? fill(true, nrow(x)) : reduce((a, b) -> a .& b, (x[!, col] for col in columns))

"""See `samples_or`."""
samples_not(column::Symbol) = x -> .!(x[!, column])

"""See `samples_or`."""
samples_or_not(columns::Vector{Symbol}) =
    x -> isempty(columns) ? fill(true, nrow(x)) : .!(reduce((a, b) -> a .| b, (x[!, col] for col in columns)))

"""See `samples_or`."""
samples_and_not(columns::Vector{Symbol}) =
    x -> isempty(columns) ? fill(true, nrow(x)) : .!(reduce((a, b) -> a .& b, (x[!, col] for col in columns)))

# Helper functions for time interval selection (Interval type)
# These return Interval values (tuples, ranges, or nothing), not predicates
"""
    times()
    times(t::Real)
    times(start::Real, stop::Real)
    times(interval::Tuple{Real,Real})

Create time interval values for the `interval_selection` keyword argument.
Unlike other predicate generators, `times()` returns an `AbstractSelection` value
(a `TimeSelection` or `AllSelection`), not a predicate function.

# Examples
```julia
times()                             # Select all time points (returns AllSelection)
times(0.3)                          # Select single time point (snaps to nearest sample)
times(-0.2, 0.5)                    # Select time window from -200ms to 500ms
times((-0.2, 0.5))                  # Same, from a Tuple
```
"""
times() = AllSelection()
times(t::Real) = TimeSelection(Float64(t), Float64(t))
times(start::Real, stop::Real) = TimeSelection(Float64(start), Float64(stop))
times(interval::Tuple{Real,Real}) = TimeSelection(Float64(interval[1]), Float64(interval[2]))




"""
    epochs()
    epochs(index::Int)
    epochs(indices::Int...)
    epochs(indices::Union{Vector{Int}, UnitRange})
    epochs(predicate::Function)

Predicate generators for epoch selection in subsetting operations.
The `predicate` variant receives individual epoch DataFrames and should return `Bool`.

# Examples
```julia
epochs()                              # Select all epochs (default)
epochs(1:5)                           # Select first 5 epochs by numerical index
epochs(2, 4, 6)                       # Select epochs 2, 4, and 6
epochs(df -> nrow(df) > 100)          # Keep epochs with >100 samples
epochs(df -> df.interval[1] < 1.0)    # Keep epochs by metadata value
```
"""
epochs() = x -> fill(true, length(x))  # Default: select all epochs given
epochs(epoch_numbers::Union{Vector{Int},UnitRange}) = x -> [i in epoch_numbers for i = 1:length(x)]
epochs(epoch_number::Int) = x -> [i == epoch_number for i = 1:length(x)]
epochs(epoch_numbers::Int...) = epochs(collect(epoch_numbers))
epochs(predicate::Function) = x -> map(predicate, x)  # Allow custom function predicates over dataframes

"""
    epochs_not(index::Int)
    epochs_not(indices::Int...)
    epochs_not(indices::Union{Vector{Int}, UnitRange})

Predicate generators that **exclude** the specified epochs. Operates inversely to `epochs`.

# Examples
```julia
epochs_not(1)                         # Exclude epoch 1
epochs_not(1:10)                      # Exclude first 10 epochs
```
"""
epochs_not(epoch_numbers::Union{Vector{Int},UnitRange}) = x -> .!([i in epoch_numbers for i = 1:length(x)])
epochs_not(epoch_number::Int) = x -> .!([i == epoch_number for i = 1:length(x)])
epochs_not(epoch_numbers::Int...) = epochs_not(collect(epoch_numbers))

"""
    participants()
    participants(id::Int)
    participants(ids::Int...)
    participants(ids::Union{Vector{Int}, UnitRange})
    participants(predicate::Function)

Predicate generators for participant ID selection in batch operations.

# Examples
```julia
participants()                        # Select all participants (default)
participants(1, 2, 3)                 # Select participants 1, 2, and 3
participants(1:10)                    # Select participants 1 through 10
```
"""
participants() = x -> fill(true, length(x))  # Default: select all participants given
participants(participant_ids::Union{Vector{Int},UnitRange}) = x -> [id in participant_ids for id in x]
participants(participant_id::Int) = x -> x .== participant_id
participants(participant_ids::Int...) = participants(collect(participant_ids))
participants(predicate::Function) = predicate  # Allow custom function predicates
"""
    participants_not(id::Int)
    participants_not(ids::Int...)
    participants_not(ids::Union{Vector{Int}, UnitRange})

Predicate generators that **exclude** the specified participant IDs.

# Examples
```julia
participants_not(5)                   # Exclude participant 5
participants_not(1:3)                 # Exclude participants 1, 2, and 3
```
"""
participants_not(participant_ids::Union{Vector{Int},UnitRange}) = x -> .!([id in participant_ids for id in x])
participants_not(participant_id::Int) = x -> .!(x .== participant_id)
participants_not(participant_ids::Int...) = participants_not(collect(participant_ids))

"""
    conditions()
    conditions(index::Int)
    conditions(indices::Int...)
    conditions(indices::Union{Vector{Int}, UnitRange})
    conditions(name::String)
    conditions(names::Vector{String})
    conditions(predicate::Function)

Predicate generators for condition selection by index or name.
Used with the `condition_selection` keyword argument on `Vector{ErpData}`,
`Vector{EpochData}`, and similar multi-condition containers.

# Examples
```julia
conditions()                          # Select all conditions (default)
conditions(1, 3)                      # Select conditions 1 and 3 by index
conditions("congruent")               # Select condition by name
conditions(["congruent", "neutral"])   # Select multiple conditions by name
```
"""
conditions() = x -> fill(true, length(x))  # Default: select all conditions given
conditions(condition_indices::Union{Vector{Int},UnitRange}) = x -> intersect(condition_indices, 1:length(x))
conditions(condition_index::Int) = x -> intersect([condition_index], 1:length(x))
conditions(condition_indices::Int...) = conditions(collect(condition_indices))
conditions(condition_names::Vector{String}) =
    x -> Vector{Int}(filter(!isnothing, [findfirst(dat -> condition_name(dat) == n, x) for n in condition_names]))
conditions(cond_name::String) = x -> Vector{Int}(filter(!isnothing, [findfirst(dat -> condition_name(dat) == cond_name, x)]))
conditions(predicate::Function) = predicate  # Allow custom function predicates
"""
    conditions_not(index::Int)
    conditions_not(indices::Int...)
    conditions_not(indices::Union{Vector{Int}, UnitRange})
    conditions_not(name::String)
    conditions_not(names::String...)
    conditions_not(names::Vector{String})

Predicate generators that **exclude** the specified conditions.

# Examples
```julia
conditions_not(1)                     # Exclude condition 1
conditions_not("incongruent")         # Exclude condition by name
```
"""
conditions_not() = x -> fill(true, length(x))
conditions_not(condition_indices::Union{Vector{Int},UnitRange}) = x -> .!([i in condition_indices for i = 1:length(x)])
conditions_not(condition_index::Int) = x -> .!([i == condition_index for i = 1:length(x)])
conditions_not(condition_indices::Int...) = conditions_not(collect(condition_indices))
conditions_not(condition_names::Vector{String}) = x -> .!([condition_name(dat) in condition_names for dat in x])
conditions_not(cond_name::String) = x -> .!([condition_name(dat) == cond_name for dat in x])
conditions_not(condition_names::String...) = conditions_not(collect(condition_names))

"""Preserve user-specified channel name ordering if all names are in the selection."""
function _handle_channel_names_order(user_order::Vector{Symbol}, selectable_cols::Vector{Symbol}, selected::Vector{Symbol})
    # Validate: check for missing channels and duplicates
    seen = Set{Symbol}()
    for ch in user_order
        if ch ∉ selectable_cols
            @minimal_warning "Channel $(ch) not found in data!"
        elseif ch ∈ seen
            @minimal_warning "Channel $(ch) already specified!"
        else
            push!(seen, ch)
        end
    end

    # Check if this is channels([...]) or channels_not([...])
    existing_in_order = [ch for ch in user_order if ch in selectable_cols]
    existing_in_selected = [ch for ch in existing_in_order if ch in selected]

    # If all existing channels are in selected, it's channels([...]) - preserve order
    if !isempty(existing_in_order) && length(existing_in_selected) == length(existing_in_order)
        result = Symbol[]
        for ch in user_order
            if ch ∉ selectable_cols || ch ∈ result
                continue
            end
            push!(result, ch)
        end
        return result
    end
    return nothing  # Use default order
end

"""Preserve user-specified channel number ordering if all indices are in the selection."""
function _handle_channel_numbers_order(user_order_numbers, selectable_cols::Vector{Symbol}, selected::Vector{Symbol})
    # Validate: check for invalid indices and duplicates
    seen = Set{Int}()
    for i in user_order_numbers
        if i < 1 || i > length(selectable_cols)
            @minimal_warning "Channel index $(i) out of range (valid: 1:$(length(selectable_cols)))!"
        elseif i ∈ seen
            @minimal_warning "Channel index $(i) already specified!"
        else
            push!(seen, i)
        end
    end

    # Check if this is channels([...]) or channels_not([...])
    valid_indices = [i for i in user_order_numbers if 1 <= i <= length(selectable_cols)]
    selected_indices = [i for i in valid_indices if selectable_cols[i] in selected]

    # If all valid indices are in selected, it's channels([...]) - preserve order
    if !isempty(valid_indices) && length(selected_indices) == length(valid_indices)
        result = Symbol[]
        for i in user_order_numbers
            if i < 1 || i > length(selectable_cols) || selectable_cols[i] ∈ result
                continue
            end
            push!(result, selectable_cols[i])
        end
        return result
    end
    return nothing  # Use default order
end

"""Apply a channel predicate to data and return the matching column names (with optional metadata/extra)."""
function get_selected_channels(dat, channel_selection::Function; include_meta::Bool = true, include_extra::Bool = true)
    # Columns/channels in dataframe to include
    metadata_cols = include_meta ? meta_labels(dat) : Symbol[]
    selectable_cols = include_extra ? vcat(channel_labels(dat), extra_labels(dat)) : channel_labels(dat)

    # Apply channel selection to non-metadata columns
    selection_mask = channel_selection(selectable_cols)
    selected = selectable_cols[selection_mask]

    # Preserve user-specified order if available
    selection_type = typeof(channel_selection)
    if hasfield(selection_type, :channel_names)
        user_order_names = getfield(channel_selection, :channel_names)
        ordered = _handle_channel_names_order(user_order_names, selectable_cols, selected)
        !isnothing(ordered) && (selected = ordered)
    elseif hasfield(selection_type, :channel_numbers)
        user_order_numbers = getfield(channel_selection, :channel_numbers)
        ordered = _handle_channel_numbers_order(user_order_numbers, selectable_cols, selected)
        !isnothing(ordered) && (selected = ordered)
    end

    return vcat(metadata_cols, selected)
end


"""Apply a component predicate and return matching indices."""
function get_selected_components(ica_result::InfoIca, component_selection::Function)
    all_components = 1:length(ica_result.ica_label)
    return all_components[component_selection(all_components)]
end

"""Apply a sample predicate and return matching row indices."""
function get_selected_samples(dat::SingleDataFrameEeg, sample_selection::Function)
    return findall(sample_selection(dat.data))
end
function get_selected_samples(dat::SingleDataFrameEeg, ::Union{Nothing,AllSelection})
    return 1:nrow(dat.data)
end

# Helper to select samples based on a predicate
function get_selected_samples(dat::MultiDataFrameEeg, sample_selection::Function)
    return findall(sample_selection(dat.data[1])) # assume all data have the same samples
end
function get_selected_samples(dat::MultiDataFrameEeg, ::Union{Nothing,AllSelection})
    return 1:nrow(dat.data[1])
end

# Helper to select samples from a DataFrame
function get_selected_samples(dat::DataFrame, sample_selection::Function)
    return findall(sample_selection(dat))
end
function get_selected_samples(dat::DataFrame, ::Union{Nothing,AllSelection})
    return 1:nrow(dat)
end

# Helper to convert TimeSelection and Tuple into sample functions
function get_selected_samples(dat::Union{SingleDataFrameEeg,MultiDataFrameEeg,DataFrame}, sel::Union{TimeSelection,Tuple{Real,Real}})
    func = _interval_to_samples(sel)
    return get_selected_samples(dat, func)
end

"""Apply an epoch predicate and return matching epoch indices."""
function get_selected_epochs(dat::MultiDataFrameEeg, epoch_selection::Function)
    # The dat.data accessor maps to data_power for TimeFreqEpochData
    return findall(epoch_selection(dat.data))
end

"""Apply a condition predicate and return matching dataset indices."""
function get_selected_conditions(datasets::Vector{<:EegFunData}, condition_selection::Function)
    res = condition_selection(datasets)
    return res isa AbstractVector{Bool} ? findall(res) : res
end




"""
    epoch_to_continuous(dat::MultiDataFrameEeg, epoch_idx::Int) -> ContinuousData

Extract a single epoch from multi-epoch EEG data and return it as ContinuousData.

# Arguments
- `dat::MultiDataFrameEeg`: The multi-DataFrame EEG data
- `epoch_idx::Int`: Index of the epoch to extract (1-based)

# Returns
- `ContinuousData`: Single DataFrame containing only the specified epoch

# Examples
```julia
# Extract epoch 3 as continuous data
single_dat = epoch_to_continuous(dat, 3)

# Now you can use single DataFrame functions
summary = channel_summary(single_dat)
```
"""
function epoch_to_continuous(dat::T, epoch_idx::Int)::ContinuousData where {T<:MultiDataFrameEeg}
    # Validate epoch index
    if epoch_idx < 1 || epoch_idx > length(dat.data)
        @minimal_error "Epoch index $epoch_idx out of range (1:$(length(dat.data)))"
    end
    return ContinuousData(dat.file, dat.data[epoch_idx], dat.layout, dat.sample_rate, dat.analysis_info)
end
