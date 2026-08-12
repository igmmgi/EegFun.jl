"""
    _make_output_filename(output_dir::String, input_file::String, suffix::String)

Create an output filename from input file path with given suffix.

# Arguments
- `output_dir::String`: Output directory path
- `input_file::String`: Input file path
- `suffix::String`: Suffix to add (e.g., "_ica", "_continuous")

# Returns
- `String`: Full output filename path

# Example
```julia
filename = _make_output_filename("/output", "data/file.bdf", "_ica")
# Returns: "/output/file_ica.jld2"
```
"""
function _make_output_filename(output_dir::String, input_file::String, suffix::String)
    base_name = basename_without_ext(input_file)
    return joinpath(output_dir, "$(base_name)$(suffix).jld2")
end


"""
    basename_without_ext(path::String)

Extract the base filename without extension from a file path.

# Arguments
- `path::String`: File path to process

# Returns
- `String`: Base filename without extension

# Example
```julia
filename = basename_without_ext("data/file.bdf")
```
"""
basename_without_ext(path::String) = splitext(basename(path))[1]



"""
    consecutive(f::Function, A::AbstractVector; step::Int=1) -> Vector

Apply function f to consecutive pairs of elements in vector A.

# Arguments
- `f::Function`: Function to apply to pairs
- `A::AbstractVector`: Input vector
- `step::Int=1`: Step size between pairs

# Returns
- `Vector`: Results of applying f to consecutive pairs
"""
function consecutive(f::Function, A::AbstractVector; step::Int = 1)
    step < 1 && @minimal_error "Step must be positive"
    (length(A) < step + 1) && @minimal_error "Vector too short for given step size"
    return [f(A[i+step], A[i]) for i = 1:(length(A)-step)]
end

"""
    _splitgroups(v::AbstractVector) -> Tuple{Vector{Int64},Vector{Int64}}

Split vector into groups based on consecutive numbers.
Refined to a single pass to avoid intermediate allocations.

# Returns
- `Tuple{Vector{Int64},Vector{Int64}}`: Start and end indices of groups
"""
function _splitgroups(v::AbstractVector{<:Integer})
    isempty(v) && return Int64[], Int64[]

    n = length(v)
    start_vals = [Int64(v[1])]
    end_vals = Int64[]

    @inbounds for i = 2:n
        if v[i] > v[i-1] + 1
            push!(end_vals, Int64(v[i-1]))
            push!(start_vals, Int64(v[i]))
        end
    end
    push!(end_vals, Int64(v[end]))

    return start_vals, end_vals
end





"""
    _find_idx_start_end(time::AbstractVector, start_time::Real, end_time::Real) -> Tuple{Int,Int}
    _find_idx_start_end(time::AbstractVector, limits::AbstractVector) -> Tuple{Int,Int}

Find start and end indices corresponding to time interval.
Assumes time vector is sorted in ascending order.

# Returns
- `Tuple{Int,Int}`: Start and end indices
"""
_find_idx_start_end(time::AbstractVector, start_time::Real, end_time::Real) =
    searchsortedfirst(time, start_time), searchsortedlast(time, end_time)
_find_idx_start_end(time::AbstractVector, limits::AbstractVector) = searchsortedfirst(time, limits[1]), searchsortedlast(time, limits[end])


"""
    find_times(time::AbstractVector, requested_times::AbstractVector) -> (indices::Vector{Int}, times::Vector{Float64})

Find nearest time points in a sorted time vector for each requested time.

For each requested time, finds the nearest matching time point in the sorted time vector
and returns both the indices and the actual time values.

# Arguments
- `time::AbstractVector`: Sorted time vector (assumed to be in ascending order)
- `requested_times::AbstractVector`: Vector of requested time points

# Returns
- `indices::Vector{Int}`: Indices of nearest matching time points
- `times::Vector{Float64}`: Actual time values at those indices

# Example
```julia
time_vec = [0.0, 0.01, 0.02, 0.03, 0.04]
requested = [0.005, 0.015, 0.025]
indices, times = find_times(time_vec, requested)
# indices = [1, 2, 3]
# times = [0.0, 0.01, 0.02]
```
"""
function find_times(time::AbstractVector, requested_times::AbstractVector)::Tuple{Vector{Int},Vector{Float64}}
    # Pre-allocate for performance if we have many requested times
    n_req = length(requested_times)
    indices = sizehint!(Int[], n_req)
    times_out = sizehint!(Float64[], n_req)

    isempty(time) && return indices, times_out

    time_min, time_max = time[1], time[end] # Assumes sorted

    for t_req in requested_times
        # Only include if within data range
        if t_req >= time_min && t_req <= time_max
            idx = searchsortedfirst(time, t_req)
            # Find nearest (check previous index if closer)
            if idx > 1 && abs(time[idx-1] - t_req) <= abs(time[min(idx, length(time))] - t_req)
                idx = idx - 1
            end
            idx = min(idx, length(time))

            # Avoid duplicates of the same index
            if isempty(indices) || indices[end] != idx
                push!(indices, idx)
                push!(times_out, time[idx])
            end
        end
    end

    return indices, times_out
end



"""
    detrend(x::AbstractVector, y::AbstractVector) -> Vector{Float64}

Remove linear trend from data using closed-form least squares regression.
More efficient than matrix-based regression `(X \\ y)`.

# Arguments
- `x::AbstractVector`: Independent variable (e.g., time points)
- `y::AbstractVector`: Dependent variable to detrend

# Returns
- `Vector{Float64}`: Detrended data with linear trend removed
"""
function detrend(x::AbstractVector, y::AbstractVector)::Vector{Float64}
    n = length(x)
    n == length(y) || @minimal_error "x and y must have the same length"
    n < 2 && @minimal_error "Need at least 2 points for detrending"

    # Use closed-form formulas for slope (m) and intercept (b)
    # m = (nΣxy - ΣxΣy) / (nΣx² - (Σx)²)
    # b = (Σy - mΣx) / n

    sum_x = sum(x)
    sum_y = sum(y)
    sum_xy = dot(x, y)
    sum_x2 = dot(x, x)

    denom = (n * sum_x2 - sum_x^2)

    # Handle zero denominator (vertical line or identical points)
    if abs(denom) < 1e-15
        return y .- (sum_y / n) # Just remove mean
    end

    m = (n * sum_xy - sum_x * sum_y) / denom
    b = (sum_y - m * sum_x) / n

    # Return detrended: y - (mx + b)
    return [y_val - (m * x_val + b) for (x_val, y_val) in zip(x, y)]
end


"""
    _extract_int(s::String) -> Union{Int, Nothing}

Extract the first integer found in a string.

# Arguments
- `s::String`: Input string

# Returns
- `Union{Int, Nothing}`: First integer found, or `nothing` if no digits

# Example
```julia
_extract_int("channel_123_data")  # Returns: 123
_extract_int("no_numbers_here")   # Returns: nothing
```
"""
function _extract_int(s::String)::Union{Int,Nothing}
    m = match(r"\d+", s)
    return isnothing(m) ? nothing : parse(Int, m.match)
end


"""
    _natural_sort_key(s::String) -> String

Generate a sort key for natural/numeric sorting of strings.
Pads numeric parts with zeros so that "file_3" sorts before "file_12".

# Arguments
- `s::String`: Input string to generate sort key for

# Returns
- `String`: Transformed string with numeric parts zero-padded

# Example
```julia
# Use with sort
files = ["file_10.jld2", "file_2.jld2", "file_1.jld2"]
sort(files, by=_natural_sort_key)  # Returns: ["file_1.jld2", "file_2.jld2", "file_10.jld2"]

# Use with DataFrame sort
sort(df, :file, by=_natural_sort_key)
```
"""
_natural_sort_key(s::String)::String = replace(s, r"\d+" => m -> lpad(String(m), 10, '0'))




"""
    Base.copy(dat::ContinuousData) -> ContinuousData
    Base.copy(dat::EpochData) -> EpochData
    Base.copy(dat::ErpData) -> ErpData
    Base.copy(info::AnalysisInfo) -> AnalysisInfo
    Base.copy(layout::Layout) -> Layout
    Base.copy(tf_data::TimeFreqData) -> TimeFreqData
    Base.copy(spectrum_data::SpectrumData) -> SpectrumData
    Base.copy(ica::InfoIca) -> InfoIca

Create a shallow copy of an EegFun data object. DataFrames are copied with
`copycols=true` for independence; small immutable fields are shared.
"""
function Base.copy(dat::ContinuousData)::ContinuousData
    return ContinuousData(dat.file, copy(dat.data, copycols = true), copy(dat.layout), dat.sample_rate, copy(dat.analysis_info))
end

function Base.copy(dat::EpochData)::EpochData
    return EpochData(
        dat.file,
        dat.condition,
        dat.condition_name,
        [copy(epoch, copycols = true) for epoch in dat.data],
        copy(dat.layout),
        dat.sample_rate,
        copy(dat.analysis_info),
    )
end

function Base.copy(dat::ErpData)::ErpData
    return ErpData(
        dat.file,
        dat.condition,
        dat.condition_name,
        copy(dat.data, copycols = true),
        copy(dat.layout),
        dat.sample_rate,
        copy(dat.analysis_info),
        dat.n_epochs,
    )
end

function Base.copy(info::AnalysisInfo)::AnalysisInfo
    return AnalysisInfo(
        reference = info.reference,
        hp_filter = info.hp_filter,
        lp_filter = info.lp_filter,
        sample_rate = info.sample_rate,
        n_ica_components_removed = info.n_ica_components_removed,
        n_channels_repaired = info.n_channels_repaired,
        repaired_channels = copy(info.repaired_channels),
    )
end

# Function to parse component input text into a list of component indices
"""
    _parse_string_to_ints(text::String) -> Vector{Int}
    _parse_string_to_ints(text::String, max_count::Int) -> Vector{Int}

Parse a comma/semicolon-separated string potentially containing ranges (e.g., "1,3:5,8")
into a sorted, unique vector of integers.

# Arguments
- `text::String`: The input string (supports `,`, `;`, and `:` for ranges)
- `max_count::Int`: Optional maximum number of results to return

# Returns
- `Vector{Int}`: Sorted, unique vector of integers found in the text.
                 Returns an empty vector if input is empty or invalid.
"""
function _parse_string_to_ints(text::String)
    isempty(text) && return Int[]

    # Check for decimal points and error (EEG channel/component indices are always integers)
    occursin('.', text) && throw(ArgumentError("Decimal points not allowed in integer selection: '$text'"))

    components = Int[]

    # Split by comma or semicolon
    for part in split(text, r"[,;]")
        s_part = strip(part)
        isempty(s_part) && continue

        if occursin(':', s_part)
            range_splits = split(s_part, ':')
            if length(range_splits) == 2
                start_str, end_str = strip(range_splits[1]), strip(range_splits[2])
                if all(isdigit, start_str) && all(isdigit, end_str)
                    s_idx, e_idx = parse(Int, start_str), parse(Int, end_str)
                    s_idx <= e_idx && append!(components, s_idx:e_idx)
                end
            end
        elseif all(isdigit, s_part)
            push!(components, parse(Int, s_part))
        else
            @minimal_warning "Skipping invalid integer selection part: '$s_part'"
        end
    end

    return unique!(sort!(components))
end

function _parse_string_to_ints(text::String, max_count::Int)
    all_indices = _parse_string_to_ints(text)
    if length(all_indices) > max_count
        resize!(all_indices, max_count)
    end
    return all_indices
end


"""
    _get_defaults(kwargs_dict::Dict{Symbol,Tuple{Any,String}}) -> Dict{Symbol,Any}

Extract default values from a keyword arguments dictionary, discarding descriptions.
"""
function _get_defaults(kwargs_dict::Dict{Symbol,Tuple{Any,String}})::Dict{Symbol,Any}
    return Dict(key => value[1] for (key, value) in kwargs_dict)
end

"""
    _merge_plot_kwargs(defaults_dict, user_kwargs; validate=true)

Helper function to merge user keyword arguments with defaults and optionally validate parameter names.

# Arguments
- `defaults_dict`: Dictionary with (default_value, description) tuples
- `user_kwargs`: User-provided keyword arguments (NamedTuple or Dict)
- `validate`: Whether to validate unknown parameters (default: true)

# Returns
- `Dict{Symbol,Any}`: Merged and validated keyword arguments

# Example
```julia
# With validation (default):
plot_kwargs = _merge_plot_kwargs(PLOT_KWARGS, kwargs)

# Without validation (for multi-component functions):
plot_kwargs = _merge_plot_kwargs(PLOT_KWARGS, kwargs; validate=false)
```
"""
function _merge_plot_kwargs(defaults_dict::Dict{Symbol,Tuple{Any,String}}, user_kwargs::NamedTuple; validate::Bool = true)::Dict{Symbol,Any}

    # Get default and user kwargs
    defaults = _get_defaults(defaults_dict)
    user_dict = Dict{Symbol,Any}(pairs(user_kwargs))

    # Check for unknown parameters (only if validation is enabled)
    if validate
        valid_keys = keys(defaults)
        unknown_keys = setdiff(keys(user_dict), valid_keys)
        if !isempty(unknown_keys)
            @minimal_error "Unknown keyword arguments: $(join(unknown_keys, ", ")). Valid arguments: $(join(valid_keys, ", "))"
        end
    end

    # Merge defaults with user kwargs
    merged_kwargs = merge(defaults, user_dict)

    return merged_kwargs
end

# Convenience function for the common pattern
function _merge_plot_kwargs(defaults_dict::Dict{Symbol,Tuple{Any,String}}, user_kwargs::Dict; validate::Bool = true)::Dict{Symbol,Any}
    return _merge_plot_kwargs(defaults_dict, NamedTuple(user_kwargs); validate = validate)
end

# Handle empty keyword arguments (Base.Pairs)
function _merge_plot_kwargs(defaults_dict::Dict{Symbol,Tuple{Any,String}}, user_kwargs::Base.Pairs; validate::Bool = true)::Dict{Symbol,Any}
    return _merge_plot_kwargs(defaults_dict, NamedTuple(user_kwargs); validate = validate)
end



"""
    _orientation(p1::Vector{Float64}, p2::Vector{Float64}, p3::Vector{Float64}) -> Int

Calculate the orientation of three points (p1, p2, p3).
Returns 0 if collinear, 1 if clockwise, 2 if counterclockwise.

This is a computational geometry utility function used in convex hull algorithms.
"""
function _orientation(p1::Vector{Float64}, p2::Vector{Float64}, p3::Vector{Float64})::Int
    val = (p2[2] - p1[2]) * (p3[1] - p2[1]) - (p2[1] - p1[1]) * (p3[2] - p2[2])
    if abs(val) < 1e-10
        return 0  # collinear
    elseif val > 0
        return 1  # clockwise
    else
        return 2  # counterclockwise
    end
end

# Helper function to generate documentation
function _generate_kwargs_doc(kwargs_dict::Dict{Symbol,Tuple{Any,String}})::String
    doc_lines = ["# Keyword Arguments"]
    push!(
        doc_lines,
        "All keyword arguments below have sensible defaults. You can override any by passing the corresponding keyword argument.",
    )
    push!(doc_lines, "")

    # Define which legend arguments are actually useful to show in the main docs
    allowed_legend_keys = [:legend, :legend_position, :legend_label, :legend_channel, :legend_nbanks]
    
    legend_keys = Symbol[]
    layout_keys = Symbol[]
    main_keys = Symbol[]

    for k in keys(kwargs_dict)
        str_k = string(k)
        if startswith(str_k, "legend")
            if k in allowed_legend_keys
                push!(legend_keys, k)
            end
        elseif startswith(str_k, "layout_")
            push!(layout_keys, k)
        else
            push!(main_keys, k)
        end
    end

    sort!(legend_keys)
    sort!(layout_keys)
    sort!(main_keys)

    function _format_default(val)
        str = repr(val)
        if length(str) > 80
            return str[1:77] * "..."
        end
        return str
    end

    function _add_section(title, keys_list)
        if !isempty(keys_list)
            push!(doc_lines, "### $(title)")
            for k in keys_list
                default_val, desc = kwargs_dict[k]
                type_info = typeof(default_val)
                formatted_val = _format_default(default_val)
                push!(doc_lines, "- `$(k)::$(type_info) = $(formatted_val)`: $(desc)")
            end
            push!(doc_lines, "")
        end
    end

    _add_section("General Plot Settings", main_keys)
    _add_section("Layout Options", layout_keys)
    _add_section("Legend Options", legend_keys)

    # Append note for the hidden Makie kwargs
    push!(doc_lines, "> **Further Legend Options**")
    push!(doc_lines, "> `EegFun.jl` passes any argument prefixed with `legend_` directly to Makie's legend system.")
    push!(doc_lines, "> While only the most common options are listed above, you can pass **any** standard Makie legend attribute (e.g., `legend_bgcolor`, `legend_patchsize`).")
    push!(doc_lines, "> See the Makie.jl documentation for a full list of discoverable arguments.")

    return join(doc_lines, "\n")
end


"""
    _combine_boolean_columns!(dat::ContinuousData, columns::Vector{Symbol}, operation::Symbol; output_column::Symbol = :combined_flags)

Combine multiple boolean columns using logical operations.

# Arguments
- `dat::ContinuousData`: The continuous EEG data object
- `columns::Vector{Symbol}`: Vector of column names to combine
- `operation::Symbol`: Logical operation to apply (:and, :or, :nand, :nor)
- `output_column::Symbol`: Name of the output column (default: :combined_flags)

# Modifies
- `dat`: Adds the combined boolean column to the data

# Examples
```julia
# Combine multiple artifact detection columns with AND operation
_combine_boolean_columns!(dat, [:is_extreme_value_100, :is_eog_onset], :and)

# Combine with OR operation
_combine_boolean_columns!(dat, [:is_extreme_value_100, :is_eog_onset], :or, output_column = :any_artifact)
```
"""
function _combine_boolean_columns!(dat::ContinuousData, columns::Vector{Symbol}, operation::Symbol; output_column::Symbol = :combined_flags)
    # Input validation
    isempty(columns) && @minimal_error "Must specify at least one column to combine"
    missing_cols = filter(col -> !hasproperty(dat.data, col), columns)
    !isempty(missing_cols) && @minimal_error "Columns not found in data: $(join(missing_cols, ", "))"
    operation ∉ [:and, :or, :nand, :nor] && @minimal_error "Invalid operation :$operation. Must be one of: :and, :or, :nand, :nor"

    # Get the boolean columns
    bool_columns = [dat.data[!, col] for col in columns]

    # Apply the logical operation
    result = if operation == :and
        all.(zip(bool_columns...))
    elseif operation == :or
        any.(zip(bool_columns...))
    elseif operation == :nand
        .!(all.(zip(bool_columns...)))
    elseif operation == :nor
        .!(any.(zip(bool_columns...)))
    end

    # Store the result
    dat.data[!, output_column] = result

    @info "_combine_boolean_columns!: Combined $(length(columns)) columns using :$operation operation into column :$output_column"
end


"""
    example_path(relative_path::String) -> String

Return the absolute path to a bundled example resource file.

This resolves paths relative to the package's `resources/` directory using `pkgdir`,
so it works regardless of the current working directory.

# Arguments
- `relative_path::String`: Path relative to the `resources/` directory

# Returns
- `String`: Absolute path to the resource file

# Examples
```julia
# Load example BDF data
dat = read_raw_data(EegFun.example_path("data/bdf/example1.bdf"))

# Load a layout
layout = read_layout(EegFun.example_path("layouts/biosemi/biosemi72.csv"))
```
"""
function example_path(relative_path::String)
    path = joinpath(pkgdir(EegFun), "resources", relative_path)
    isfile(path) || isdir(path) || error("Example resource not found: $relative_path\nExpected at: $path")
    return path
end
