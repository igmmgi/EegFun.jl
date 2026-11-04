struct PaddedVector{T,V<:AbstractVector{T}} <: AbstractVector{T}
    data::V
    pad_len::Int
end
Base.size(p::PaddedVector) = (length(p.data) + 2 * p.pad_len,)
Base.IndexStyle(::Type{<:PaddedVector}) = IndexLinear()
@inline function Base.getindex(p::PaddedVector, i::Int)
    n = length(p.data)
    pad = p.pad_len
    @inbounds if i <= pad
        return p.data[pad-i+1]
    elseif i > n + pad
        return p.data[2*n+pad-i+1]
    else
        return p.data[i-pad]
    end
end

const _ResampleData = Union{ContinuousData,ErpData,EpochData}

function _validate_target_rate(target_rate::Real)::Int
    (!isfinite(target_rate) || target_rate <= 0 || !isinteger(target_rate)) &&
        @minimal_error("Target rate must be a positive integer (Hz), got $target_rate")
    return Int(target_rate)
end

function _prepare_resample_spec(old_rate::Real, target_rate::Real)
    rate_spec = (isinteger(old_rate) && isinteger(target_rate)) ? (Int(target_rate) // Int(old_rate)) : Float64(target_rate / old_rate)
    rate = Float64(rate_spec)
    h = DSP.resample_filter(rate_spec)
    return rate_spec, rate, h
end

function _resample_padding(N_old::Int, old_rate::Real, rate_spec, rate::Float64)
    if rate_spec isa Rational
        q = denominator(rate_spec)
        target_pad = min(N_old ÷ 2, Int(old_rate))
        pad_len = (target_pad ÷ q) * q
        (pad_len == 0 && q <= N_old) && (pad_len = q)
        return pad_len, Int(pad_len * rate_spec)
    else
        pad_len = min(N_old ÷ 2, round(Int, old_rate))
        return pad_len, round(Int, pad_len * rate)
    end
end

function _resample_signal_vec(sig::Vector{Float64}, rate_spec, h, pad_len::Int, pad_new::Int, N_new::Int)::Vector{Float64}
    if pad_len > 0
        resampled_padded = DSP.resample(PaddedVector(sig, pad_len), rate_spec, h)
        return resampled_padded[(pad_new+1):(pad_new+N_new)]
    else
        res = DSP.resample(sig, rate_spec, h)
        return length(res) == N_new ? res : res[1:min(end, N_new)]
    end
end

function _resample_df_col(df::DataFrame, pos::Int, rate_spec, h, pad_len::Int, pad_new::Int, N_new::Int)::Vector{Float64}
    col_data = df[!, pos]
    sig = col_data isa Vector{Float64} ? col_data : convert(Vector{Float64}, col_data)
    return _resample_signal_vec(sig, rate_spec, h, pad_len, pad_new, N_new)
end

function _resample_dataframe_prepared(
    df::DataFrame,
    old_rate::Real,
    target_rate::Real,
    rate_spec,
    rate::Float64,
    h,
    channel_set::Set{Symbol},
    trigger_col::Symbol,
    extra_set::Set{Symbol},
    thread_channels::Bool,
)
    N_old = nrow(df)
    N_new = ceil(Int, N_old * rate_spec)
    pad_len, pad_new = _resample_padding(N_old, old_rate, rate_spec, rate)

    col_names = propertynames(df)
    n_cols = length(col_names)
    out_cols = Vector{AbstractVector}(undef, n_cols)

    sig_positions = Int[]
    other_positions = Int[]
    for (pos, col) in enumerate(col_names)
        if col in channel_set
            push!(sig_positions, pos)
        elseif col == :time
            out_cols[pos] = collect(range(df.time[1], step = 1 / target_rate, length = N_new))
        elseif col == trigger_col
            trig_src = df[!, trigger_col]
            trig_out = zeros(eltype(trig_src), N_new)
            @inbounds for orig_idx = 1:N_old
                v = trig_src[orig_idx]
                if v != 0
                    new_idx = clamp(round(Int, (orig_idx - 1) * rate) + 1, 1, N_new)
                    trig_out[new_idx] = v
                end
            end
            out_cols[pos] = trig_out
        elseif col in extra_set && eltype(df[!, col]) <: AbstractFloat
            col_data = df[!, col]
            if allequal(col_data)
                out_cols[pos] = fill(col_data[1], N_new)
            else
                push!(sig_positions, pos)
            end
        else
            push!(other_positions, pos)
        end
    end

    if thread_channels && length(sig_positions) > 1 && length(sig_positions) * N_old >= 16_384
        Threads.@threads :static for idx in eachindex(sig_positions)
            out_cols[sig_positions[idx]] = _resample_df_col(df, sig_positions[idx], rate_spec, h, pad_len, pad_new, N_new)
        end
    else
        for pos in sig_positions
            out_cols[pos] = _resample_df_col(df, pos, rate_spec, h, pad_len, pad_new, N_new)
        end
    end

    if !isempty(other_positions)
        old_indices = Int[clamp(round(Int, (i - 1) / rate) + 1, 1, N_old) for i = 1:N_new]
        for pos in other_positions
            col_data = df[!, pos]
            if eltype(col_data) == Bool
                out_bool = rate > 1.0 ? Bool[col_data[old_indices[i]] for i = 1:N_new] : fill(false, N_new)
                @inbounds for orig_idx = 1:N_old
                    if col_data[orig_idx]
                        new_idx = clamp(round(Int, (orig_idx - 1) * rate) + 1, 1, N_new)
                        out_bool[new_idx] = true
                    end
                end
                out_cols[pos] = out_bool
            else
                out_cols[pos] = col_data[old_indices]
            end
        end
    end

    return DataFrame(out_cols, col_names; copycols = false)
end

function _resample_dataframe(
    df::DataFrame,
    old_rate::Real,
    target_rate::Real,
    channel_cols::Vector{Symbol},
    trigger_col::Symbol,
    extra_cols::Vector{Symbol} = Symbol[],
)
    rate_spec, rate, h = _prepare_resample_spec(old_rate, target_rate)
    return _resample_dataframe_prepared(
        df,
        old_rate,
        target_rate,
        rate_spec,
        rate,
        h,
        Set(channel_cols),
        trigger_col,
        Set(extra_cols),
        true,
    )
end

"""
    resample!(dat::Union{ContinuousData, ErpData, EpochData}, target_rate::Real)
    resample!(data_vec::AbstractVector{<:_ResampleData}, target_rate::Real)

Resample EEG data in-place to integer `target_rate` (Hz) using polyphase anti-alias filtering.
"""
function resample!(dat::Union{ContinuousData,ErpData}, target_rate::Real)
    target_rate_int = _validate_target_rate(target_rate)
    if target_rate_int == dat.sample_rate
        @info "Target rate matches current sample rate, no resampling needed"
        return dat
    end

    @info "Resampling data from $(dat.sample_rate) Hz to $(target_rate_int) Hz"
    dat.data = _resample_dataframe(dat.data, dat.sample_rate, target_rate_int, channel_labels(dat), :trigger, extra_labels(dat))
    hasproperty(dat.data, :sample) && (dat.data.sample = collect(1:n_samples(dat)))

    dat.sample_rate = target_rate_int
    dat.analysis_info.sample_rate = Float64(target_rate_int)
    @info "Resampling complete. New sample rate: $(dat.sample_rate) Hz, $(n_samples(dat)) samples"
    return dat
end

function _resample_single_epoch(
    epoch::DataFrame,
    old_rate::Real,
    target_rate_int::Int,
    rate_spec,
    rate::Float64,
    h,
    ch_set::Set{Symbol},
    ex_set::Set{Symbol},
    thread_channels::Bool,
)::DataFrame
    new_df = _resample_dataframe_prepared(epoch, old_rate, target_rate_int, rate_spec, rate, h, ch_set, :trigger, ex_set, thread_channels)
    if hasproperty(new_df, :sample)
        first_sample = max(1, round(Int, (epoch.sample[1] - 1) * rate) + 1)
        n_samp = nrow(new_df)
        new_df.sample = collect(first_sample:(first_sample+n_samp-1))
    end
    return new_df
end

function resample!(dat::EpochData, target_rate::Real)
    target_rate_int = _validate_target_rate(target_rate)
    if target_rate_int == dat.sample_rate
        @info "Target rate matches current sample rate, no resampling needed"
        return dat
    end

    @info "Resampling $(n_epochs(dat)) epochs from $(dat.sample_rate) Hz to $(target_rate_int) Hz"
    old_rate = dat.sample_rate
    rate_spec, rate, h = _prepare_resample_spec(old_rate, target_rate_int)
    ch_set = Set(channel_labels(dat))
    ex_set = Set(extra_labels(dat))
    n_ep = length(dat.data)

    if n_ep > 1
        Threads.@threads :static for i = 1:n_ep
            dat.data[i] = _resample_single_epoch(dat.data[i], old_rate, target_rate_int, rate_spec, rate, h, ch_set, ex_set, false)
        end
    elseif n_ep == 1
        dat.data[1] = _resample_single_epoch(dat.data[1], old_rate, target_rate_int, rate_spec, rate, h, ch_set, ex_set, true)
    end

    dat.sample_rate = target_rate_int
    dat.analysis_info.sample_rate = Float64(target_rate_int)
    !isempty(dat.data) && @info "Resampling complete. New sample rate: $(dat.sample_rate) Hz, $(nrow(dat.data[1])) samples per epoch"
    return dat
end

"""
    resample(dat::Union{ContinuousData,ErpData,EpochData}, target_rate::Real)
    resample(data_vec::AbstractVector{<:Union{ContinuousData,ErpData,EpochData}}, target_rate::Real)
    resample(file_pattern::String, target_rate::Real; input_dir = pwd(), participant_selection = participants(), output_dir = nothing)

Non-mutating or batch JLD2 version of [`resample!`](@ref).
"""
function resample(dat::_ResampleData, target_rate::Real)
    _validate_target_rate(target_rate)
    return resample!(copy(dat), target_rate)
end

function resample(data_vec::AbstractVector{<:_ResampleData}, target_rate::Real)
    _validate_target_rate(target_rate)
    return [resample(dat, target_rate) for dat in data_vec]
end

function resample!(data_vec::AbstractVector{<:_ResampleData}, target_rate::Real)
    _validate_target_rate(target_rate)
    for dat in data_vec
        resample!(dat, target_rate)
    end
    return data_vec
end

_default_resample_output_dir(input_dir::String, pattern::String, target_rate::Real) =
    joinpath(input_dir, "resampled_by_$(target_rate)_$(clean_pattern(pattern))")

function _process_resample_file(filepath::String, output_path::String, target_rate::Real)
    filename = basename(filepath)
    loaded_data = read_data(filepath)
    isnothing(loaded_data) && return BatchResult(false, filename, "No data found in file")

    is_single = loaded_data isa _ResampleData
    is_vec = loaded_data isa AbstractVector{<:_ResampleData}
    (is_single || is_vec) || return BatchResult(false, filename, "Data is not a recognized EEG data type")
    (is_vec && isempty(loaded_data)) && return BatchResult(false, filename, "Empty data vector in file")

    old_rate = is_vec ? first(loaded_data).sample_rate : loaded_data.sample_rate
    resampled_data = resample(loaded_data, target_rate)
    new_rate = is_vec ? first(resampled_data).sample_rate : resampled_data.sample_rate
    jldsave(output_path; data = resampled_data, format_version = EEGFUN_FORMAT_VERSION)
    return BatchResult(true, filename, "Resampled from $old_rate Hz to $new_rate Hz (target rate: $target_rate)")
end

function resample(
    file_pattern::String,
    target_rate::Real;
    input_dir::String = pwd(),
    participant_selection = participants(),
    output_dir::Union{String,Nothing} = nothing,
)
    _with_batch_logging("resampling", "resample.log", input_dir, output_dir) do ctx
        @log_call "resample"
        target_rate_int = _validate_target_rate(target_rate)
        output_dir = _prepare_output_dir!(ctx, _default_resample_output_dir(input_dir, file_pattern, target_rate_int))
        process_fn = (input_path, output_path) -> _process_resample_file(input_path, output_path, target_rate_int)
        batch_process(process_fn, file_pattern, input_dir, output_dir, participant_selection, "Resampling")
    end
end
