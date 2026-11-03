
struct PaddedVector{T,V<:AbstractVector{T}} <: AbstractVector{T}
    data::V
    pad_len::Int
end
Base.size(p::PaddedVector) = (length(p.data) + 2*p.pad_len,)
Base.IndexStyle(::Type{<:PaddedVector}) = IndexLinear()
function Base.getindex(p::PaddedVector, i::Int)
    if i <= p.pad_len
        p.data[p.pad_len-i+1]
    elseif i > length(p.data) + p.pad_len
        p.data[end-(i-length(p.data)-p.pad_len)+1]
    else
        p.data[i-p.pad_len]
    end
end

const _ResampleData = Union{ContinuousData,ErpData,EpochData}

function _validate_target_rate(target_rate::Real)::Int
    if !isfinite(target_rate) || target_rate <= 0 || !isinteger(target_rate)
        @minimal_error("Target rate must be a positive integer (Hz), got $target_rate")
    end
    return Int(target_rate)
end

function _resample_signal_vec(sig::Vector{Float64}, rate_spec, h, pad_len::Int, pad_new::Int, N_new::Int)::Vector{Float64}
    if pad_len > 0
        padded_sig = PaddedVector(sig, pad_len)
        resampled_padded = DSP.resample(padded_sig, rate_spec, h)
        return resampled_padded[(pad_new+1):(pad_new+N_new)]
    else
        res = DSP.resample(sig, rate_spec, h)
        return length(res) == N_new ? res : res[1:min(end, N_new)]
    end
end

"""
Helper function to resample a dataframe using DSP.jl polyphase filtering.
"""
function _resample_dataframe(
    df::DataFrame,
    old_rate::Real,
    target_rate::Real,
    channel_cols::Vector{Symbol},
    trigger_col::Symbol,
    extra_cols::Vector{Symbol} = Symbol[],
)
    rate_spec = (isinteger(old_rate) && isinteger(target_rate)) ? (Int(target_rate) // Int(old_rate)) : Float64(target_rate / old_rate)
    rate = Float64(rate_spec)
    N_old = nrow(df)
    N_new = ceil(Int, N_old * rate_spec)

    df_new = DataFrame()

    # 1. Resample signal channels with integer-aligned mirror padding
    h = DSP.resample_filter(rate_spec)
    if rate_spec isa Rational
        q = denominator(rate_spec)
        target_pad = min(N_old ÷ 2, Int(old_rate))
        pad_len = (target_pad ÷ q) * q
        if pad_len == 0 && q <= N_old
            pad_len = q
        end
        pad_new = Int(pad_len * rate_spec)
    else
        pad_len = min(N_old ÷ 2, round(Int, old_rate))
        pad_new = round(Int, pad_len * rate)
    end

    for ch in channel_cols
        if hasproperty(df, ch)
            sig = convert(Vector{Float64}, df[!, ch])
            df_new[!, ch] = _resample_signal_vec(sig, rate_spec, h, pad_len, pad_new, N_new)
        end
    end

    # 2. Resample time column if it exists
    if hasproperty(df, :time)
        t0 = df.time[1]
        df_new.time = range(t0, step = 1 / target_rate, length = N_new)
    end

    # 3. Handle other metadata and extra columns
    other_cols = setdiff(propertynames(df), [channel_cols..., :time])

    if hasproperty(df, trigger_col)
        trigger_type = eltype(df[!, trigger_col])
        df_new[!, trigger_col] = zeros(trigger_type, N_new)

        trigger_indices = findall(df[!, trigger_col] .!= 0)
        trigger_values = df[!, trigger_col][trigger_indices]

        for (orig_idx, trig_val) in zip(trigger_indices, trigger_values)
            new_idx = round(Int, (orig_idx - 1) * rate) + 1
            new_idx = clamp(new_idx, 1, N_new)
            df_new[!, trigger_col][new_idx] = trig_val
        end

        other_cols = setdiff(other_cols, [trigger_col])
    end

    if !isempty(other_cols)
        old_indices = [clamp(round(Int, (i - 1) / rate) + 1, 1, N_old) for i = 1:N_new]
        extra_set = Set(extra_cols)
        for col in other_cols
            col_data = df[!, col]
            if eltype(col_data) == Bool
                # OR-reduce boolean flags per output sample bin so transient flags are never lost
                out_bool = rate > 1.0 ? collect(Bool, col_data[old_indices]) : falses(N_new)
                for orig_idx = 1:N_old
                    if col_data[orig_idx]
                        new_idx = clamp(round(Int, (orig_idx - 1) * rate) + 1, 1, N_new)
                        out_bool[new_idx] = true
                    end
                end
                df_new[!, col] = out_bool
            elseif col in extra_set && eltype(col_data) <: AbstractFloat
                # Anti-alias continuous float extra channels (e.g. EOG/bipolar), preserving exact constants (e.g. RT)
                if allequal(col_data)
                    df_new[!, col] = fill(col_data[1], N_new)
                else
                    sig = convert(Vector{Float64}, col_data)
                    df_new[!, col] = _resample_signal_vec(sig, rate_spec, h, pad_len, pad_new, N_new)
                end
            else
                df_new[!, col] = col_data[old_indices]
            end
        end
    end

    return df_new[:, propertynames(df)]
end

"""
    resample!(dat::Union{ContinuousData, ErpData}, target_rate::Real)
    resample!(dat::EpochData, target_rate::Real)
    resample!(data_vec::AbstractVector{<:_ResampleData}, target_rate::Real)
    resample(dat, target_rate)
    resample(data_vec::AbstractVector{<:_ResampleData}, target_rate)

Resample EEG data to a specific integer `target_rate` (Hz) using DSP.jl's polyphase filtering.
Automatically applies an anti-aliasing filter and supports arbitrary integer target rates (including non-integer resampling ratios).
"""
function resample!(dat::Union{ContinuousData,ErpData}, target_rate::Real)
    target_rate_int = _validate_target_rate(target_rate)

    if target_rate_int == dat.sample_rate
        @info "Target rate matches current sample rate, no resampling needed"
        return dat
    end

    @info "Resampling data from $(dat.sample_rate) Hz to $(target_rate_int) Hz"

    chans = channel_labels(dat)
    extras = extra_labels(dat)
    dat.data = _resample_dataframe(dat.data, dat.sample_rate, target_rate_int, chans, :trigger, extras)

    if hasproperty(dat.data, :sample)
        dat.data.sample = 1:n_samples(dat)
    end

    dat.sample_rate = target_rate_int
    dat.analysis_info.sample_rate = Float64(target_rate_int)
    @info "Resampling complete. New sample rate: $(dat.sample_rate) Hz, $(n_samples(dat)) samples"

    return dat
end

function resample!(dat::EpochData, target_rate::Real)
    target_rate_int = _validate_target_rate(target_rate)

    if target_rate_int == dat.sample_rate
        @info "Target rate matches current sample rate, no resampling needed"
        return dat
    end

    @info "Resampling $(n_epochs(dat)) epochs from $(dat.sample_rate) Hz to $(target_rate_int) Hz"

    chans = channel_labels(dat)
    extras = extra_labels(dat)
    rate = target_rate_int / dat.sample_rate
    for (i, epoch) in enumerate(dat.data)
        dat.data[i] = _resample_dataframe(epoch, dat.sample_rate, target_rate_int, chans, :trigger, extras)

        if hasproperty(dat.data[i], :sample)
            first_sample = max(1, round(Int, (epoch.sample[1] - 1) * rate) + 1)
            n_samp = nrow(dat.data[i])
            dat.data[i].sample = first_sample:(first_sample+n_samp-1)
        end
    end

    dat.sample_rate = target_rate_int
    dat.analysis_info.sample_rate = Float64(target_rate_int)
    if !isempty(dat.data)
        @info "Resampling complete. New sample rate: $(dat.sample_rate) Hz, $(nrow(dat.data[1])) samples per epoch"
    end

    return dat
end

"""
    resample(dat::Union{ContinuousData,ErpData,EpochData}, target_rate::Real)
    resample(data_vec::AbstractVector{<:Union{ContinuousData,ErpData,EpochData}}, target_rate::Real)

Non-mutating version of `resample!`.
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


# === BATCH PROCESSING FUNCTIONS ===

"""Generate default output directory name for resampling operation."""
function _default_resample_output_dir(input_dir::String, pattern::String, target_rate::Real)
    new_rate_str = "resampled_by_$(target_rate)"
    joinpath(input_dir, "$(new_rate_str)_$(clean_pattern(pattern))")
end


"""
Process a single data file through resampling pipeline.
Returns BatchResult with success/failure info.
"""
function _process_resample_file(filepath::String, output_path::String, target_rate::Real)
    filename = basename(filepath)

    # Read data using read_data (handles single variable files automatically)
    loaded_data = read_data(filepath)

    if isnothing(loaded_data)
        return BatchResult(false, filename, "No data found in file")
    end

    is_single = loaded_data isa _ResampleData
    is_vec = loaded_data isa AbstractVector{<:_ResampleData}
    if !(is_single || is_vec)
        return BatchResult(false, filename, "Data is not a recognized EEG data type")
    end
    if is_vec && isempty(loaded_data)
        return BatchResult(false, filename, "Empty data vector in file")
    end

    # Resample (exceptions are handled by _run_batch_operation)
    old_rate = is_vec ? first(loaded_data).sample_rate : loaded_data.sample_rate
    resampled_data = resample(loaded_data, target_rate)
    new_rate = is_vec ? first(resampled_data).sample_rate : resampled_data.sample_rate

    # Save results (always use "data" as variable name since read_data finds by type)
    jldsave(output_path; data = resampled_data, format_version = EEGFUN_FORMAT_VERSION)

    message = "Resampled from $old_rate Hz to $new_rate Hz (target rate: $target_rate)"
    return BatchResult(true, filename, message)
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

        # Setup directories
        output_dir = _prepare_output_dir!(ctx, _default_resample_output_dir(input_dir, file_pattern, target_rate_int))

        process_fn = (input_path, output_path) -> _process_resample_file(input_path, output_path, target_rate_int)
        result = batch_process(process_fn, file_pattern, input_dir, output_dir, participant_selection, "Resampling")
    end
end
