# Per-channel descriptive statistics and variance summaries.

######################################
########### User Functions ###########
######################################

"""
    channel_summary(dat::SingleDataFrameEeg; sample_selection = samples(), interval_selection::Interval = times(),
                    channel_selection = channels(), include_meta::Bool = false, include_extra::Bool = false) -> DataFrame
    channel_summary(dat::MultiDataFrameEeg; sample_selection = samples(), interval_selection::Interval = times(),
                    channel_selection = channels(), include_meta::Bool = false, include_extra::Bool = false) -> DataFrame
    channel_summary(file_pattern::String; input_dir = pwd(), participant_selection = participants(),
                    condition_selection = conditions(), sample_selection = samples(), interval_selection = times(),
                    channel_selection = channels(), include_extra::Bool = false, output_dir = nothing,
                    output_file::String = "channel_summary")

Compute summary statistics (`min`, `max`, `std`, `range`, `var`, `zvar`) per channel (and per epoch for
`MultiDataFrameEeg`, or saved to CSV across JLD2 files in batch mode).
"""
function channel_summary(
    dat::SingleDataFrameEeg;
    sample_selection = samples(),
    interval_selection::Interval = times(),
    channel_selection = channels(),
    include_meta::Bool = false,
    include_extra::Bool = false,
)::DataFrame
    n_samples(dat) == 0 && @minimal_error("Cannot compute channel summary: data is empty")
    selected_channels = get_selected_channels(dat, channel_selection; include_meta = include_meta, include_extra = include_extra)
    selected_samples = get_selected_samples(dat, _combine_interval_sample(interval_selection, sample_selection))
    return _channel_summary_impl(dat.data, selected_samples, selected_channels)
end

function channel_summary(
    dat::MultiDataFrameEeg;
    sample_selection = samples(),
    interval_selection::Interval = times(),
    channel_selection = channels(),
    include_meta::Bool = false,
    include_extra::Bool = false,
)::DataFrame
    isempty(dat.data) && @minimal_error("Cannot compute channel summary: no epochs in data")

    valid_epochs = Tuple{Int,DataFrame}[]
    sizehint!(valid_epochs, length(dat.data))
    for (epoch_idx, epoch_df) in enumerate(dat.data)
        if nrow(epoch_df) == 0
            @minimal_warning("Skipping empty epoch $(epoch_idx)")
            continue
        end
        orig_ep = hasproperty(epoch_df, :epoch) && !isempty(epoch_df.epoch) ? Int(epoch_df.epoch[1]) : epoch_idx
        push!(valid_epochs, (orig_ep, epoch_df))
    end
    isempty(valid_epochs) && @minimal_error("No valid epochs found for channel summary")

    first_df = valid_epochs[1][2]
    first_single = ContinuousData(dat.file, first_df, dat.layout, dat.sample_rate, dat.analysis_info)
    selected_channels = get_selected_channels(first_single, channel_selection; include_meta = include_meta, include_extra = include_extra)
    isempty(selected_channels) && @minimal_error("No channels selected for channel summary")
    col_indices = Int[columnindex(first_df, ch) for ch in selected_channels]

    all_samples_fast =
        (sample_selection isa AllSelection || sample_selection isa _AllSamples) &&
        (isnothing(interval_selection) || interval_selection isa AllSelection)
    combined_sel = all_samples_fast ? nothing : _combine_interval_sample(interval_selection, sample_selection)

    n_valid = length(valid_epochs)
    n_ch = length(selected_channels)
    total_rows = n_valid * n_ch

    epoch_col = Vector{Int}(undef, total_rows)
    channel_col = Vector{Symbol}(undef, total_rows)
    min_col = Vector{Float64}(undef, total_rows)
    max_col = Vector{Float64}(undef, total_rows)
    std_col = Vector{Float64}(undef, total_rows)
    range_col = Vector{Float64}(undef, total_rows)
    var_col = Vector{Float64}(undef, total_rows)
    zvar_col = Vector{Float64}(undef, total_rows)

    empty_idx = Int[]
    process_epoch = function (k::Int)
        orig_ep, epoch_df = valid_epochs[k]
        r_start = (k - 1) * n_ch + 1
        r_end = k * n_ch
        rng = r_start:r_end
        @inbounds for i = 1:n_ch
            epoch_col[r_start+i-1] = orig_ep
            channel_col[r_start+i-1] = selected_channels[i]
        end

        if all_samples_fast
            _write_channel_summary_slice!(
                epoch_df,
                col_indices,
                1:nrow(epoch_df),
                empty_idx,
                @view(min_col[rng]),
                @view(max_col[rng]),
                @view(std_col[rng]),
                @view(range_col[rng]),
                @view(var_col[rng]),
                @view(zvar_col[rng]),
            )
        else
            sel_samps = get_selected_samples(epoch_df, combined_sel)
            isempty(sel_samps) && @minimal_error("No samples selected for channel summary")
            contig = _contiguous_sample_range(sel_samps)
            _write_channel_summary_slice!(
                epoch_df,
                col_indices,
                contig,
                sel_samps,
                @view(min_col[rng]),
                @view(max_col[rng]),
                @view(std_col[rng]),
                @view(range_col[rng]),
                @view(var_col[rng]),
                @view(zvar_col[rng]),
            )
        end
    end

    if all_samples_fast && n_valid > 1 && total_rows * nrow(first_df) >= 16_384
        Threads.@threads :static for k = 1:n_valid
            process_epoch(k)
        end
    else
        for k = 1:n_valid
            process_epoch(k)
        end
    end

    return DataFrame(
        epoch = epoch_col,
        channel = channel_col,
        min = min_col,
        max = max_col,
        std = std_col,
        range = range_col,
        var = var_col,
        zvar = zvar_col;
        copycols = false,
    )
end

channel_summary(dat::AbstractVector{<:MultiDataFrameEeg}; kwargs...) = channel_summary.(dat; kwargs...)

function channel_summary(
    file_pattern::String;
    input_dir::String = pwd(),
    participant_selection = participants(),
    condition_selection = conditions(),
    sample_selection = samples(),
    interval_selection::Interval = times(),
    channel_selection = channels(),
    include_extra::Bool = false,
    output_dir::Union{String,Nothing} = nothing,
    output_file::String = "channel_summary",
)
    _with_batch_logging("channel summary", "$(output_file).log", input_dir, output_dir) do ctx
        @log_call "channel_summary"
        out_dir = _prepare_output_dir!(ctx, _default_channel_summary_output_dir(input_dir, file_pattern))
        files = _find_batch_files(file_pattern, input_dir, participant_selection)
        if isempty(files)
            @minimal_warning "No JLD2 files found matching pattern '$file_pattern' in $input_dir"
            return nothing
        end

        @info "Found $(length(files)) JLD2 files to process"
        all_summaries = DataFrame[]
        n_success, n_error = 0, 0

        for (i, file) in enumerate(files)
            @info "Channel summary: $file ($i/$(length(files)))"
            result, summary_dfs = try
                _process_channel_summary_file(
                    joinpath(input_dir, file),
                    condition_selection,
                    sample_selection,
                    interval_selection,
                    channel_selection,
                    include_extra,
                )
            catch e
                @error "Error processing $file" exception = (e, catch_backtrace())
                (BatchResult(false, file, "Exception: $(sprint(showerror, e))"), DataFrame[])
            end

            if result.success
                @info "  ✓ $(result.message)"
                append!(all_summaries, summary_dfs)
                n_success += 1
            else
                @minimal_warning "  ✗ $(result.message)"
                n_error += 1
            end
        end

        if !isempty(all_summaries)
            output_path = joinpath(out_dir, "$(output_file).csv")
            CSV.write(output_path, reduce(vcat, all_summaries))
            @info "Combined results saved to: $output_path"
        else
            @minimal_warning "No results to save"
        end

        @info "Batch operation complete! Processed $n_success files successfully, $n_error errors"
        @info "Output saved to: $out_dir"
    end
end

######################################
######### Internal Functions #########
######################################

@inline function _col_summary_stats(src::AbstractVector)
    n = length(src)
    v1 = Float64(@inbounds src[1])
    n == 1 && return (v1, v1, NaN, 0.0, NaN)

    mn = mx = v1
    s = 0.0
    @inbounds for i = 1:n
        x = Float64(src[i])
        x < mn && (mn = x)
        x > mx && (mx = x)
        s += x
    end
    μ = s / n
    ssq = 0.0
    @inbounds @simd for i = 1:n
        d = Float64(src[i]) - μ
        ssq = muladd(d, d, ssq)
    end
    v = ssq / (n - 1)
    return (mn, mx, sqrt(v), mx - mn, v)
end

function _fill_zvar!(zvar_slice::AbstractVector{Float64}, var_slice::AbstractVector{Float64})
    n = length(var_slice)
    if all(==(0.0), var_slice)
        fill!(zvar_slice, 0.0)
        return nothing
    end
    n == 1 && (zvar_slice[1] = NaN; return nothing)
    μ = mean(var_slice)
    s = std(var_slice)
    s == 0.0 && (fill!(zvar_slice, 0.0); return nothing)
    inv_s = 1.0 / s
    @inbounds @simd for i = 1:n
        zvar_slice[i] = (var_slice[i] - μ) * inv_s
    end
    return nothing
end

function _write_channel_summary_slice!(
    data::DataFrame,
    col_indices::Vector{Int},
    contig::Union{Nothing,UnitRange{Int}},
    sample_selection::Vector{Int},
    min_s::AbstractVector{Float64},
    max_s::AbstractVector{Float64},
    std_s::AbstractVector{Float64},
    range_s::AbstractVector{Float64},
    var_s::AbstractVector{Float64},
    zvar_s::AbstractVector{Float64},
)
    @inbounds for c in eachindex(col_indices)
        col = data[!, col_indices[c]]
        mn, mx, sd, rg, vr = _col_summary_stats(!isnothing(contig) ? @view(col[contig]) : @view(col[sample_selection]))
        min_s[c] = mn
        max_s[c] = mx
        std_s[c] = sd
        range_s[c] = rg
        var_s[c] = vr
    end
    _fill_zvar!(zvar_s, var_s)
    return nothing
end

function _channel_summary_impl(data::DataFrame, sample_selection::Vector{Int}, channel_selection::Vector{Symbol})::DataFrame
    isempty(sample_selection) && @minimal_error("No samples selected for channel summary")
    isempty(channel_selection) && @minimal_error("No channels selected for channel summary")

    col_indices = Vector{Int}(undef, length(channel_selection))
    for (i, ch) in enumerate(channel_selection)
        idx = columnindex(data, ch)
        idx == 0 && @minimal_error("Channels not found in data: $(setdiff(channel_selection, propertynames(data)))")
        col_indices[i] = idx
    end

    n_rows = nrow(data)
    min_s, max_s = extrema(sample_selection)
    if min_s < 1 || max_s > n_rows
        invalid_samples = sample_selection[(sample_selection .< 1) .| (sample_selection .> n_rows)]
        @minimal_error("Invalid sample indices: $(invalid_samples)")
    end

    contig = _contiguous_sample_range(sample_selection)
    n_ch = length(channel_selection)
    min_vec = Vector{Float64}(undef, n_ch)
    max_vec = Vector{Float64}(undef, n_ch)
    std_vec = Vector{Float64}(undef, n_ch)
    range_vec = Vector{Float64}(undef, n_ch)
    var_vec = Vector{Float64}(undef, n_ch)
    zvar_vec = Vector{Float64}(undef, n_ch)

    _write_channel_summary_slice!(data, col_indices, contig, sample_selection, min_vec, max_vec, std_vec, range_vec, var_vec, zvar_vec)

    return DataFrame(
        channel = copy(channel_selection),
        min = min_vec,
        max = max_vec,
        std = std_vec,
        range = range_vec,
        var = var_vec,
        zvar = zvar_vec;
        copycols = false,
    )
end

_default_channel_summary_output_dir(input_dir::String, pattern::String) = joinpath(input_dir, "channel_summary_$(_clean_pattern(pattern))")

function _process_channel_summary_file(
    filepath::String,
    condition_selection,
    sample_selection,
    interval_selection::Interval,
    channel_selection,
    include_extra::Bool,
)
    filename = basename(filepath)
    data_var = read_data(filepath)
    isnothing(data_var) && return (BatchResult(false, filename, "No recognized data variable"), DataFrame[])
    (data_var isa EegData || data_var isa AbstractVector{<:EegData}) || return (
        BatchResult(false, filename, "Invalid data type: expected EegData or Vector{<:EegData}, got $(typeof(data_var))"),
        DataFrame[],
    )

    items = data_var isa AbstractVector ? data_var : [data_var]
    selected_indices = get_selected_conditions(items, condition_selection)
    selected_items = items[selected_indices]

    summary_dfs = DataFrame[]
    file_stem = splitext(filename)[1]
    for (cond_idx, data) in enumerate(selected_items)
        condition = hasfield(typeof(data), :condition) ? data.condition : selected_indices[cond_idx]
        summary_df = channel_summary(
            data;
            sample_selection = sample_selection,
            interval_selection = interval_selection,
            channel_selection = channel_selection,
            include_extra = include_extra,
        )
        insertcols!(summary_df, 1, :file => file_stem)
        insertcols!(summary_df, 2, :condition => condition)
        push!(summary_dfs, summary_df)
    end

    return (BatchResult(true, filename, "Processed $(length(summary_dfs)) condition(s)"), summary_dfs)
end
