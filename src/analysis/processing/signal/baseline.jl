# Time-domain baseline subtraction for continuous, epoched, and ERP data.

######################################
########### User Functions ###########
######################################

"""
    baseline(dat::Union{ContinuousData,ErpData,EpochData}, baseline_interval::Interval = times(); channel_selection = channels())
    baseline(dat::AbstractVector{<:Union{ContinuousData,ErpData,EpochData}}, baseline_interval::Interval = times(); channel_selection = channels())
    baseline(file_pattern::String, baseline_interval; input_dir = pwd(),
             participant_selection = participants(), condition_selection = conditions(), output_dir = nothing)

Subtract the mean amplitude over `baseline_interval` from each selected channel without modifying `dat`
(or batch-process matching `.jld2` files when given `file_pattern::String`).
Omitting `baseline_interval` uses the full time range; passing `nothing` returns a copy unchanged.
"""
baseline(dat::_TimeDomainEegData, args...; kwargs...) = baseline!(copy(dat), args...; kwargs...)

baseline(dat::AbstractVector{<:_TimeDomainEegData}, args...; kwargs...) = [baseline(d, args...; kwargs...) for d in dat]

"""
    baseline!(dat::Union{ContinuousData,ErpData,EpochData}, baseline_interval::Interval = times(); channel_selection = channels())
    baseline!(dat::AbstractVector{<:Union{ContinuousData,ErpData,EpochData}}, baseline_interval::Interval = times(); channel_selection = channels())

In-place mutating version of [`baseline`](@ref).
"""
function baseline!(dat::_TimeDomainEegData, baseline_interval::Interval = times(); channel_selection = channels())
    isnothing(baseline_interval) && return dat
    bl_idx = _validate_baseline_interval(dat, baseline_interval)
    selected_channels = get_selected_channels(dat, channel_selection, include_meta = false, include_extra = false)
    if isempty(selected_channels)
        @minimal_warning "No channels selected for baseline correction"
        return dat
    end

    time = time_vector(dat)
    t_start = round(time[bl_idx[1]], digits = 3)
    t_end = round(time[bl_idx[2]], digits = 3)
    @info "Applying baseline correction to $(length(selected_channels)) channels over interval $(bl_idx[1]):$(bl_idx[2]) ($(t_start) to $(t_end) s)"
    _apply_baseline!(dat.data, selected_channels, bl_idx)
    return dat
end

baseline!(dat::AbstractVector{<:_TimeDomainEegData}, baseline_interval::Interval = times(); channel_selection = channels()) =
    isnothing(baseline_interval) ? dat : (foreach(d -> baseline!(d, baseline_interval; channel_selection = channel_selection), dat); dat)

function baseline(
    file_pattern::String,
    baseline_interval::Interval;
    input_dir::String = pwd(),
    participant_selection = participants(),
    condition_selection = conditions(),
    output_dir::Union{String,Nothing} = nothing,
)
    _with_batch_logging("baseline correction", "baseline.log", input_dir, output_dir) do ctx
        _ = output_dir
        @log_call "baseline"
        bl_tuple = _validate_baseline_interval(baseline_interval)
        default_dir = joinpath(input_dir, "baseline_$(_clean_pattern(file_pattern))_$(bl_tuple[1])_to_$(bl_tuple[2])")
        return _batch_jld2_files(
            data -> (baseline!(data, bl_tuple), "Baseline corrected successfully"),
            AbstractVector{<:Union{ErpData,EpochData}},
            ctx,
            default_dir,
            file_pattern,
            input_dir,
            participant_selection,
            "Baseline correction";
            condition_selection = condition_selection,
        )
    end
end

######################################
######### Internal Functions #########
######################################

function _apply_baseline!(dat::DataFrame, ch_indices::Vector{Int}, baseline_interval::Tuple{Int,Int})
    start_idx, end_idx = baseline_interval
    inv_len = 1.0 / (end_idx - start_idx + 1)
    for col_idx in ch_indices
        col = dat[!, col_idx]
        if col isa Vector{Float64}
            col_f64 = col::Vector{Float64}
            s = 0.0
            @inbounds @simd for i = start_idx:end_idx
                s += col_f64[i]
            end
            mean_val = s * inv_len
            @inbounds @simd for i in eachindex(col_f64)
                col_f64[i] -= mean_val
            end
        else
            mean_val = mean(@view col[start_idx:end_idx])
            dat[!, col_idx] = col .- mean_val
        end
    end
    return nothing
end

_apply_baseline!(dat::DataFrame, channels::AbstractVector{Symbol}, baseline_interval::Tuple{Int,Int}) =
    _apply_baseline!(dat, Int[columnindex(dat, ch) for ch in channels], baseline_interval)

function _apply_baseline!(dat::AbstractVector{DataFrame}, channels::AbstractVector{Symbol}, baseline_interval::Tuple{Int,Int})
    isempty(dat) && return nothing
    ch_indices = Int[columnindex(dat[1], ch) for ch in channels]
    if length(dat) > 1 && Threads.nthreads() > 1 && length(dat) * length(ch_indices) * nrow(dat[1]) >= 65_536
        Threads.@threads for df in dat
            _apply_baseline!(df, ch_indices, baseline_interval)
        end
    else
        for df in dat
            _apply_baseline!(df, ch_indices, baseline_interval)
        end
    end
    return nothing
end

function _to_interval(interval::Tuple{Real,Real})
    start, stop = Float64(interval[1]), Float64(interval[2])
    start > stop && @minimal_error("Baseline start ($start) must be <= stop ($stop)")
    return (start, stop)
end

function _to_interval(interval::TimeSelection)
    start, stop = interval.start, interval.stop
    start > stop && @minimal_error("Baseline start ($start) must be <= stop ($stop)")
    return (start, stop)
end

function _to_interval(interval::AbstractVector{<:Real})
    length(interval) == 2 || @minimal_error("Baseline interval vector must have 2 elements [start, stop], got: $interval")
    start, stop = Float64(interval[1]), Float64(interval[2])
    start > stop && @minimal_error("Baseline start ($start) must be <= stop ($stop)")
    return (start, stop)
end

_validate_baseline_interval(baseline_interval::Interval)::Tuple{Real,Real} = _to_interval(baseline_interval)

_validate_baseline_interval(time::AbstractVector, ::AllSelection)::Tuple{Int,Int} = (1, length(time))

function _validate_baseline_interval(time::AbstractVector, baseline_interval::Interval)::Tuple{Int,Int}
    baseline_tuple = _to_interval(baseline_interval)
    time_min, time_max = first(time), last(time)
    start_time, stop_time = baseline_tuple
    if start_time < time_min || stop_time > time_max
        @minimal_error "Invalid baseline_interval: $(baseline_tuple) - time values must be within data range ($(time_min) to $(time_max))"
    end

    result = _find_idx_start_end(time, start_time, stop_time)
    isnothing(result) && @minimal_error(
        "Invalid baseline_interval: $(baseline_tuple) - interval not found in time vector (range: $(time_min) to $(time_max))"
    )

    start_idx, stop_idx = result
    if !(1 <= start_idx <= length(time)) || !(1 <= stop_idx <= length(time)) || !(start_idx <= stop_idx)
        @minimal_error "Invalid baseline_interval: ($start_idx, $stop_idx) (time vector length: $(length(time)))"
    end
    return (start_idx, stop_idx)
end

function _validate_baseline_interval(dat::MultiDataFrameEeg, baseline_interval::Interval)::Tuple{Int,Int}
    if baseline_interval isa Function
        indices = findall(baseline_interval(dat.data[1]))
        isempty(indices) && @minimal_error("Baseline window selection returned no samples!")
        return (first(indices), last(indices))
    end
    return _validate_baseline_interval(dat.data[1].time, baseline_interval)
end

function _validate_baseline_interval(dat::SingleDataFrameEeg, baseline_interval::Interval)::Tuple{Int,Int}
    if baseline_interval isa Function
        indices = findall(baseline_interval(dat.data))
        isempty(indices) && @minimal_error("Baseline window selection returned no samples!")
        if length(indices) > 1 && indices[end] - indices[1] + 1 != length(indices)
            @minimal_warning "Selected baseline window is not contiguous!"
        end
        return (first(indices), last(indices))
    end
    return _validate_baseline_interval(dat.data.time, baseline_interval)
end
