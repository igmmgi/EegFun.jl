# Polyphase anti-aliased resampling and downsampling for EEG data.

######################################
######### Types & Constants ##########
######################################

struct PaddedVector{T,V<:AbstractVector{T}} <: AbstractVector{T}
    data::V
    pad_len::Int
end

const _ResampleData = Union{ContinuousData,ErpData,EpochData}

######################################
########### User Functions ###########
######################################

"""
    resample(dat::Union{ContinuousData,ErpData,EpochData}, target_rate::Real)
    resample(data_vec::AbstractVector{<:Union{ContinuousData,ErpData,EpochData}}, target_rate::Real)
    resample(file_pattern::String, target_rate::Real; input_dir = pwd(), participant_selection = participants(), output_dir = nothing)

Non-mutating or batch JLD2 version of [`resample!`](@ref).
"""
function resample(dat::_ResampleData, target_rate::Real)
    target_rate_int = _validate_target_rate(target_rate)
    if target_rate_int == dat.sample_rate
        @info "Target rate matches current sample rate, no resampling needed"
        return copy(dat)
    end
    out = _create_subset(dat, _compute_resampled_data(dat, target_rate_int), copy(dat.layout))
    return _finalize_resample!(out, target_rate_int)
end

function resample(data_vec::AbstractVector{<:_ResampleData}, target_rate::Real)
    _validate_target_rate(target_rate)
    return [resample(dat, target_rate) for dat in data_vec]
end

"""
    resample!(dat::Union{ContinuousData, ErpData, EpochData}, target_rate::Real)
    resample!(data_vec::AbstractVector{<:_ResampleData}, target_rate::Real)

Resample EEG data in-place to integer `target_rate` (Hz) using polyphase anti-alias filtering.
"""
function resample!(dat::_ResampleData, target_rate::Real)
    target_rate_int = _validate_target_rate(target_rate)
    if target_rate_int == dat.sample_rate
        @info "Target rate matches current sample rate, no resampling needed"
        return dat
    end
    dat.data = _compute_resampled_data(dat, target_rate_int)
    return _finalize_resample!(dat, target_rate_int)
end

function resample!(data_vec::AbstractVector{<:_ResampleData}, target_rate::Real)
    _validate_target_rate(target_rate)
    foreach(dat -> resample!(dat, target_rate), data_vec)
    return data_vec
end

function resample(
    file_pattern::String,
    target_rate::Real;
    input_dir::String = pwd(),
    participant_selection = participants(),
    output_dir::Union{String,Nothing} = nothing,
)
    _with_batch_logging("resampling", "resample.log", input_dir, output_dir) do ctx
        _ = output_dir
        @log_call "resample"
        target_rate_int = _validate_target_rate(target_rate)
        default_dir = joinpath(input_dir, "resampled_by_$(target_rate_int)_$(_clean_pattern(file_pattern))")
        return _batch_jld2_files(
            Union{_ResampleData,AbstractVector{<:_ResampleData}},
            ctx,
            default_dir,
            file_pattern,
            input_dir,
            participant_selection,
            "Resampling",
        ) do loaded_data
            (loaded_data isa AbstractVector && isempty(loaded_data)) && return BatchResult(false, "", "Empty data vector in file")
            old_rate = loaded_data isa AbstractVector ? first(loaded_data).sample_rate : loaded_data.sample_rate
            resampled_data = resample(loaded_data, target_rate_int)
            new_rate = resampled_data isa AbstractVector ? first(resampled_data).sample_rate : resampled_data.sample_rate
            return (resampled_data, "Resampled from $old_rate Hz to $new_rate Hz (target rate: $target_rate_int)")
        end
    end
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

######################################
######### Internal Functions #########
######################################

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
        elseif col == :sample && eltype(df[!, col]) <: Integer
            first_sample = isempty(df) ? 1 : max(1, round(Int, (df[1, col] - 1) * rate) + 1)
            out_cols[pos] = collect(first_sample:(first_sample+N_new-1))
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
        elseif col == :trigger_info && eltype(df[!, col]) <: AbstractString
            info_src = df[!, col]
            info_out = fill("", N_new)
            has_trig = hasproperty(df, trigger_col)
            trig_vec = has_trig ? df[!, trigger_col] : nothing
            @inbounds for orig_idx = 1:N_old
                s = info_src[orig_idx]
                if !isempty(s) || (has_trig && trig_vec[orig_idx] != 0)
                    new_idx = clamp(round(Int, (orig_idx - 1) * rate) + 1, 1, N_new)
                    if !isempty(s)
                        info_out[new_idx] = s
                    end
                end
            end
            out_cols[pos] = info_out
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

function _compute_resampled_data(dat::Union{ContinuousData,ErpData}, target_rate_int::Int)::DataFrame
    @info "Resampling data from $(dat.sample_rate) Hz to $(target_rate_int) Hz"
    new_df = _resample_dataframe(dat.data, dat.sample_rate, target_rate_int, channel_labels(dat), :trigger, extra_labels(dat))
    hasproperty(new_df, :sample) && (new_df.sample = collect(1:nrow(new_df)))
    return new_df
end

function _compute_resampled_data(dat::EpochData, target_rate_int::Int)::Vector{DataFrame}
    n_ep = length(dat.data)
    @info "Resampling $n_ep epochs from $(dat.sample_rate) Hz to $(target_rate_int) Hz"
    old_rate = dat.sample_rate
    rate_spec, rate, h = _prepare_resample_spec(old_rate, target_rate_int)
    ch_set = Set(channel_labels(dat))
    ex_set = Set(extra_labels(dat))
    out_dfs = Vector{DataFrame}(undef, n_ep)
    if n_ep > 1
        Threads.@threads :static for i = 1:n_ep
            out_dfs[i] =
                _resample_dataframe_prepared(dat.data[i], old_rate, target_rate_int, rate_spec, rate, h, ch_set, :trigger, ex_set, false)
        end
    elseif n_ep == 1
        out_dfs[1] =
            _resample_dataframe_prepared(dat.data[1], old_rate, target_rate_int, rate_spec, rate, h, ch_set, :trigger, ex_set, true)
    end
    return out_dfs
end

function _finalize_resample!(dat::_ResampleData, target_rate_int::Int)
    dat.sample_rate = target_rate_int
    dat.analysis_info.sample_rate = Float64(target_rate_int)
    if dat isa EpochData
        !isempty(dat.data) && @info "Resampling complete. New sample rate: $(dat.sample_rate) Hz, $(nrow(dat.data[1])) samples per epoch"
    else
        @info "Resampling complete. New sample rate: $(dat.sample_rate) Hz, $(n_samples(dat)) samples"
    end
    return dat
end
