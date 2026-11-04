_make_output_filename(output_dir::String, input_file::String, suffix::String) =
    joinpath(output_dir, "$(basename_without_ext(input_file))$(suffix).jld2")

function basename_without_ext(path::String)
    base = basename(path)
    endswith(lowercase(base), ".xdf.gz") && return base[1:(end-7)]
    return splitext(base)[1]
end

function consecutive(f::Function, A::AbstractVector; step::Int = 1)
    step < 1 && @minimal_error "Step must be positive"
    length(A) < step + 1 && @minimal_error "Vector too short for given step size"
    return [f(A[i+step], A[i]) for i = 1:(length(A)-step)]
end

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

_find_idx_start_end(time::AbstractVector, start_time::Real, end_time::Real) =
    searchsortedfirst(time, start_time), searchsortedlast(time, end_time)
_find_idx_start_end(time::AbstractVector, limits::AbstractVector) = searchsortedfirst(time, limits[1]), searchsortedlast(time, limits[end])

"""
    find_times(time::AbstractVector, requested_times::AbstractVector) -> (indices::Vector{Int}, times::Vector{Float64})

Find the nearest deduplicated sample indices and time values in sorted `time` for each value in
`requested_times` falling within `[first(time), last(time)]`.
"""
function find_times(time::AbstractVector, requested_times::AbstractVector)::Tuple{Vector{Int},Vector{Float64}}
    n_req = length(requested_times)
    indices = sizehint!(Int[], n_req)
    times_out = sizehint!(Float64[], n_req)
    n = length(time)
    n == 0 && return indices, times_out

    @inbounds begin
        time_min = time[1]
        time_max = time[n]
        last_idx = 0
        for t_req in requested_times
            if t_req >= time_min && t_req <= time_max
                idx = searchsortedfirst(time, t_req)
                if idx > 1 && abs(time[idx-1] - t_req) <= abs(time[idx] - t_req)
                    idx -= 1
                end
                if idx != last_idx
                    push!(indices, idx)
                    push!(times_out, Float64(time[idx]))
                    last_idx = idx
                end
            end
        end
    end
    return indices, times_out
end

"""
    detrend(x::AbstractVector, y::AbstractVector) -> Vector{Float64}

Remove the least-squares linear trend of `y` with respect to `x`.
"""
function detrend(x::AbstractVector, y::AbstractVector)::Vector{Float64}
    n = length(x)
    n == length(y) || @minimal_error "x and y must have the same length"
    n < 2 && @minimal_error "Need at least 2 points for detrending"

    sum_x = 0.0
    sum_y = 0.0
    sum_xy = 0.0
    sum_x2 = 0.0
    @inbounds @simd for i in eachindex(x, y)
        xi = Float64(x[i])
        yi = Float64(y[i])
        sum_x += xi
        sum_y += yi
        sum_xy = muladd(xi, yi, sum_xy)
        sum_x2 = muladd(xi, xi, sum_x2)
    end

    denom = n * sum_x2 - sum_x * sum_x
    out = Vector{Float64}(undef, n)
    if abs(denom) < 1e-15
        mean_y = sum_y / n
        @inbounds @simd for i = 1:n
            out[i] = Float64(y[i]) - mean_y
        end
        return out
    end

    m = (n * sum_xy - sum_x * sum_y) / denom
    b = (sum_y - m * sum_x) / n
    @inbounds @simd for i = 1:n
        out[i] = Float64(y[i]) - muladd(m, Float64(x[i]), b)
    end
    return out
end

function _extract_int(s::String)::Union{Int,Nothing}
    cu = codeunits(s)
    n = length(cu)
    i = 1
    @inbounds while i <= n && !(0x30 <= cu[i] <= 0x39)
        i += 1
    end
    i > n && return nothing
    j = i + 1
    @inbounds while j <= n && (0x30 <= cu[j] <= 0x39)
        j += 1
    end
    return parse(Int, SubString(s, i, j - 1))
end

function _natural_sort_key(s::AbstractString)::String
    cu = codeunits(s)
    n = length(cu)
    extra = 0
    has_short_run = false
    i = 1
    @inbounds while i <= n
        if 0x30 <= cu[i] <= 0x39
            j = i + 1
            while j <= n && (0x30 <= cu[j] <= 0x39)
                j += 1
            end
            run_len = j - i
            if run_len < 10
                extra += 10 - run_len
                has_short_run = true
            end
            i = j
        else
            i += 1
        end
    end
    !has_short_run && return String(s)

    out = Base.StringVector(n + extra)
    pos = 1
    i = 1
    @inbounds while i <= n
        b = cu[i]
        if 0x30 <= b <= 0x39
            j = i + 1
            while j <= n && (0x30 <= cu[j] <= 0x39)
                j += 1
            end
            pad = 10 - (j - i)
            while pad > 0
                out[pos] = 0x30
                pos += 1
                pad -= 1
            end
            while i < j
                out[pos] = cu[i]
                pos += 1
                i += 1
            end
        else
            out[pos] = b
            pos += 1
            i += 1
        end
    end
    return String(out)
end

Base.copy(dat::ContinuousData)::ContinuousData =
    ContinuousData(dat.file, copy(dat.data, copycols = true), copy(dat.layout), dat.sample_rate, copy(dat.analysis_info))

Base.copy(dat::EpochData)::EpochData = EpochData(
    dat.file,
    dat.condition,
    dat.condition_name,
    [copy(epoch, copycols = true) for epoch in dat],
    copy(dat.layout),
    dat.sample_rate,
    copy(dat.analysis_info),
)

Base.copy(dat::ErpData)::ErpData = ErpData(
    dat.file,
    dat.condition,
    dat.condition_name,
    copy(dat.data, copycols = true),
    copy(dat.layout),
    dat.sample_rate,
    copy(dat.analysis_info),
    dat.n_epochs,
)

Base.copy(info::AnalysisInfo)::AnalysisInfo = AnalysisInfo(
    reference = info.reference,
    hp_filter = info.hp_filter,
    lp_filter = info.lp_filter,
    sample_rate = info.sample_rate,
    n_ica_components_removed = info.n_ica_components_removed,
    n_channels_repaired = info.n_channels_repaired,
    repaired_channels = copy(info.repaired_channels),
)

function _parse_string_to_ints(text::String)
    isempty(text) && return Int[]
    occursin('.', text) && @minimal_error("Decimal points not allowed in integer selection: '$text'")

    components = Int[]
    for part in eachsplit(text, (',', ';'))
        s_part = strip(part)
        isempty(s_part) && continue
        if occursin(':', s_part)
            range_splits = split(s_part, ':')
            if length(range_splits) == 2
                start_str, end_str = strip(range_splits[1]), strip(range_splits[2])
                if !isempty(start_str) && !isempty(end_str) && all(isdigit, start_str) && all(isdigit, end_str)
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
    length(all_indices) > max_count && resize!(all_indices, max_count)
    return all_indices
end

_get_defaults(kwargs_dict::Dict{Symbol,Tuple{Any,String}})::Dict{Symbol,Any} = Dict{Symbol,Any}(k => v[1] for (k, v) in kwargs_dict)

function _merge_plot_kwargs(defaults_dict::Dict{Symbol,Tuple{Any,String}}, user_kwargs::NamedTuple; validate::Bool = true)::Dict{Symbol,Any}
    defaults = _get_defaults(defaults_dict)
    if validate
        unknown_keys = [k for k in keys(user_kwargs) if !haskey(defaults, k)]
        if !isempty(unknown_keys)
            @minimal_error "Unknown keyword arguments: $(join(unknown_keys, ", ")). Valid arguments: $(join(keys(defaults), ", "))"
        end
    end
    for (k, v) in pairs(user_kwargs)
        defaults[k] = v
    end
    return defaults
end

_merge_plot_kwargs(defaults_dict::Dict{Symbol,Tuple{Any,String}}, user_kwargs::Dict; validate::Bool = true)::Dict{Symbol,Any} =
    _merge_plot_kwargs(defaults_dict, NamedTuple(user_kwargs); validate = validate)
_merge_plot_kwargs(defaults_dict::Dict{Symbol,Tuple{Any,String}}, user_kwargs::Base.Pairs; validate::Bool = true)::Dict{Symbol,Any} =
    _merge_plot_kwargs(defaults_dict, NamedTuple(user_kwargs); validate = validate)

function _orientation(p1::Vector{Float64}, p2::Vector{Float64}, p3::Vector{Float64})::Int
    val = (p2[2] - p1[2]) * (p3[1] - p2[1]) - (p2[1] - p1[1]) * (p3[2] - p2[2])
    abs(val) < 1e-10 && return 0
    return val > 0 ? 1 : 2
end

function _generate_kwargs_doc(kwargs_dict::Dict{Symbol,Tuple{Any,String}})::String
    doc_lines = [
        "# Keyword Arguments",
        "All keyword arguments below have sensible defaults. You can override any by passing the corresponding keyword argument.",
        "",
    ]

    allowed_legend_keys = (:legend, :legend_position, :legend_label, :legend_channel, :legend_nbanks)
    legend_keys = Symbol[]
    layout_keys = Symbol[]
    main_keys = Symbol[]

    for k in keys(kwargs_dict)
        str_k = string(k)
        if startswith(str_k, "legend")
            k in allowed_legend_keys && push!(legend_keys, k)
        elseif startswith(str_k, "layout_")
            push!(layout_keys, k)
        else
            push!(main_keys, k)
        end
    end
    sort!(legend_keys)
    sort!(layout_keys)
    sort!(main_keys)

    function _add_section(title, keys_list)
        isempty(keys_list) && return
        push!(doc_lines, "### $(title)")
        for k in keys_list
            default_val, desc = kwargs_dict[k]
            str = repr(default_val)
            formatted_val = length(str) > 80 ? first(str, 77) * "..." : str
            push!(doc_lines, "- `$(k)::$(typeof(default_val)) = $(formatted_val)`: $(desc)")
        end
        push!(doc_lines, "")
    end

    _add_section("General Plot Settings", main_keys)
    _add_section("Layout Options", layout_keys)
    _add_section("Legend Options", legend_keys)

    push!(
        doc_lines,
        "> **Further Legend Options**",
        "> `EegFun.jl` passes any argument prefixed with `legend_` directly to Makie's legend system.",
        "> While only the most common options are listed above, you can pass **any** standard Makie legend attribute (e.g., `legend_bgcolor`, `legend_patchsize`).",
        "> See the Makie.jl documentation for a full list of discoverable arguments.",
    )
    return join(doc_lines, "\n")
end

function _combine_boolean_columns!(dat::ContinuousData, columns::Vector{Symbol}, operation::Symbol; output_column::Symbol = :combined_flags)
    isempty(columns) && @minimal_error "Must specify at least one column to combine"
    for col in columns
        if !hasproperty(dat.data, col)
            missing_cols = filter(c -> !hasproperty(dat.data, c), columns)
            @minimal_error "Columns not found in data: $(join(missing_cols, ", "))"
        end
    end
    operation in (:and, :or, :nand, :nor) || @minimal_error "Invalid operation :$operation. Must be one of: :and, :or, :nand, :nor"

    n = nrow(dat.data)
    result = Vector{Bool}(undef, n)
    _copy_bools!(result, dat.data[!, columns[1]])
    if operation === :and || operation === :nand
        for k = 2:length(columns)
            _and_bools!(result, dat.data[!, columns[k]])
        end
        if operation === :nand
            @inbounds for i = 1:n
                result[i] = !result[i]
            end
        end
    else
        for k = 2:length(columns)
            _or_bools!(result, dat.data[!, columns[k]])
        end
        if operation === :nor
            @inbounds for i = 1:n
                result[i] = !result[i]
            end
        end
    end

    dat.data[!, output_column] = result
    @info "_combine_boolean_columns!: Combined $(length(columns)) columns using :$operation operation into column :$output_column"
end

"""
    example_path(relative_path::String) -> String

Return the absolute path to a bundled resource file under `pkgdir(EegFun)/resources/`.
"""
function example_path(relative_path::String)
    path = joinpath(pkgdir(EegFun), "resources", relative_path)
    isfile(path) || isdir(path) || @minimal_error("Example resource not found: $relative_path\nExpected at: $path")
    return path
end

function _foreach_chunk(f, n::Integer; ntasks::Integer = Threads.nthreads())
    n <= 0 && return nothing
    ntasks = clamp(ntasks, 1, n)
    if ntasks == 1
        f(1:n)
        return nothing
    end
    chunk_size = cld(n, ntasks)
    tasks = [Threads.@spawn f(chunk) for chunk in Iterators.partition(1:n, chunk_size)]
    foreach(wait, tasks)
    return nothing
end
