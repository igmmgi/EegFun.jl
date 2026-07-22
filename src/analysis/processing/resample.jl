
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

"""
Helper function to resample a dataframe using DSP.jl polyphase filtering.
"""
function _resample_dataframe(df::DataFrame, old_rate::Real, target_rate::Real, channel_cols::Vector{Symbol}, trigger_col::Symbol)
    rate = target_rate / old_rate
    N_old = nrow(df)
    N_new = ceil(Int, N_old * rate)

    df_new = DataFrame()

    # 1. Resample signal channels
    # Pre-compute the polyphase filter to avoid redesigning it for every channel
    h = DSP.resample_filter(rate)
    pad_len = min(N_old ÷ 2, round(Int, old_rate))
    pad_new = round(Int, pad_len * rate)

    for ch in channel_cols
        if hasproperty(df, ch)
            sig = convert(Vector{Float64}, df[!, ch])
            if pad_len > 0
                padded_sig = PaddedVector(sig, pad_len)
                resampled_padded = DSP.resample(padded_sig, rate, h)
                df_new[!, ch] = resampled_padded[(pad_new+1):(pad_new+N_new)]
            else
                df_new[!, ch] = DSP.resample(sig, rate, h)
            end
        end
    end

    # 2. Resample time column if it exists
    if hasproperty(df, :time)
        t0 = df.time[1]
        df_new.time = range(t0, step = 1/target_rate, length = N_new)
    end

    # 3. Handle other metadata columns
    other_cols = setdiff(propertynames(df), [channel_cols..., :time])

    if hasproperty(df, trigger_col)
        trigger_type = eltype(df[!, trigger_col])
        df_new[!, trigger_col] = zeros(trigger_type, N_new)

        trigger_indices = findall(df[!, trigger_col] .!= 0)
        trigger_values = df[!, trigger_col][trigger_indices]

        for (orig_idx, trig_val) in zip(trigger_indices, trigger_values)
            new_idx = round(Int, orig_idx * rate)
            new_idx = clamp(new_idx, 1, N_new)
            df_new[!, trigger_col][new_idx] = trig_val
        end

        other_cols = setdiff(other_cols, [trigger_col])
    end

    if !isempty(other_cols)
        old_indices = [clamp(round(Int, i / rate), 1, N_old) for i = 1:N_new]
        for col in other_cols
            df_new[!, col] = df[!, col][old_indices]
        end
    end

    return df_new[:, propertynames(df)]
end

"""
    resample!(dat::Union{ContinuousData, ErpData}, target_rate::Real)
    resample!(dat::EpochData, target_rate::Real)
    resample!(data_vec::Vector{T}, target_rate::Real) where {T<:EegData}
    resample(dat, target_rate)
    resample(data_vec::Vector{T}, target_rate) where {T<:EegData}

Resample EEG data to a specific `target_rate` using DSP.jl's polyphase filtering.
Automatically applies an anti-aliasing filter. Supports fractional resampling.
"""
function resample!(dat::SingleDataFrameEeg, target_rate::Real)::Nothing
    if target_rate <= 0
        @minimal_error("Target rate must be positive, got $target_rate")
    end

    if target_rate == dat.sample_rate
        @info "Target rate matches current sample rate, no resampling needed"
        return nothing
    end

    @info "Resampling data from $(dat.sample_rate) Hz to $(target_rate) Hz"

    chans = channel_labels(dat)
    dat.data = _resample_dataframe(dat.data, dat.sample_rate, target_rate, chans, :trigger)

    if hasproperty(dat.data, :sample)
        dat.data.sample = 1:nrow(dat.data)
    end

    dat.sample_rate = Float64(target_rate)
    @info "Resampling complete. New sample rate: $(dat.sample_rate) Hz, $(nrow(dat.data)) samples"

    return nothing
end

function resample!(dat::EpochData, target_rate::Real)::Nothing
    if target_rate <= 0
        @minimal_error("Target rate must be positive, got $target_rate")
    end

    if target_rate == dat.sample_rate
        @info "Target rate matches current sample rate, no resampling needed"
        return nothing
    end

    @info "Resampling $(length(dat.data)) epochs from $(dat.sample_rate) Hz to $(target_rate) Hz"

    chans = channel_labels(dat)
    for (i, epoch) in enumerate(dat.data)
        dat.data[i] = _resample_dataframe(epoch, dat.sample_rate, target_rate, chans, :trigger)

        if hasproperty(dat.data[i], :sample)
            first_sample = round(Int, epoch.sample[1] * (target_rate / dat.sample_rate))
            n_samples = nrow(dat.data[i])
            dat.data[i].sample = first_sample:(first_sample+n_samples-1)
        end
    end

    dat.sample_rate = Float64(target_rate)
    @info "Resampling complete. New sample rate: $(dat.sample_rate) Hz, $(nrow(dat.data[1])) samples per epoch"

    return nothing
end

function resample(dat::T, target_rate::Real)::T where {T<:EegData}
    dat_copy = copy(dat)
    resample!(dat_copy, target_rate)
    return dat_copy
end

function resample(data_vec::Vector{T}, target_rate::Real)::Vector{T} where {T<:EegData}
    return [resample(dat, target_rate) for dat in data_vec]
end

function resample!(data_vec::Vector{T}, target_rate::Real)::Nothing where {T<:EegData}
    for dat in data_vec
        resample!(dat, target_rate)
    end
    return nothing
end


# =============================================================================
#     BATCH PROCESSING FUNCTIONS
# =============================================================================

"""Generate default output directory name for resampling operation."""
function _default_resample_output_dir(input_dir::String, pattern::String, target_rate::Real)
    new_rate_str = "resampled_by_$(target_rate)"
    joinpath(input_dir, "$(new_rate_str)_$(pattern)")
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

    if !(loaded_data isa Union{ContinuousData,EpochData,ErpData})
        return BatchResult(false, filename, "Data is not a recognized EEG data type")
    end

    # Resample
    try
        old_rate = loaded_data.sample_rate
        resampled_data = resample(loaded_data, target_rate)
        new_rate = resampled_data.sample_rate

        # Save results (always use "data" as variable name since read_data finds by type)
        jldsave(output_path; data = resampled_data)

        message = "Resampled from $old_rate Hz to $new_rate Hz (target rate: $target_rate)"
        return BatchResult(true, filename, message)
    catch e
        return BatchResult(false, filename, "Error: $(sprint(showerror, e))")
    end
end



function resample(
    file_pattern::String,
    target_rate::Real;
    input_dir::String = pwd(),
    participant_selection::Function = participants(),
    output_dir::Union{String,Nothing} = nothing,
)

    # Setup logging
    log_file = "resample.log"
    setup_global_logging(log_file)

    try
        @info "Batch resampling started at $(now())"
        @log_call "resample"

        # Validation
        error_msg = _validate_input_dir(input_dir)
        if !isnothing(error_msg)
            @minimal_error(error_msg)
        end

        if target_rate <= 0
            @minimal_error("Target rate must be positive, got $target_rate")
        end

        # Setup directories
        output_dir = something(output_dir, _default_resample_output_dir(input_dir, file_pattern, target_rate))
        mkpath(output_dir)

        # Find files
        files = _find_batch_files(file_pattern, input_dir, participant_selection)

        if isempty(files)
            @minimal_warning "No JLD2 files found matching pattern '$file_pattern' in $input_dir"
            return nothing
        end

        @info "Found $(length(files)) JLD2 files matching pattern '$file_pattern'"
        @info "Downsampling target rate: $target_rate"

        # Create processing function with captured parameters
        process_fn = (input_path, output_path) -> _process_resample_file(input_path, output_path, target_rate)

        # Execute batch operation
        results = _run_batch_operation(process_fn, files, input_dir, output_dir; operation_name = "Resampling")

        _log_batch_summary(results, output_dir)

    finally
        _cleanup_logging(log_file, output_dir)
    end
end
