function _calculate_reference!(reference::Vector{Float64}, dat::DataFrame, ref_indices::AbstractVector{<:Integer})
    n = length(reference)
    if isempty(ref_indices)
        fill!(reference, 0.0)
        return reference
    end
    c1 = dat[!, ref_indices[1]]
    if c1 isa Vector{Float64}
        v1 = c1::Vector{Float64}
        @inbounds @simd for i = 1:n
            reference[i] = v1[i]
        end
    else
        @inbounds for i = 1:n
            reference[i] = Float64(c1[i])
        end
    end
    @inbounds for k = 2:length(ref_indices)
        col = dat[!, ref_indices[k]]
        if col isa Vector{Float64}
            col_f64 = col::Vector{Float64}
            @simd for i = 1:n
                reference[i] += col_f64[i]
            end
        else
            reference .+= col
        end
    end
    inv_refs = 1.0 / length(ref_indices)
    @inbounds @simd for i = 1:n
        reference[i] *= inv_refs
    end
    return reference
end

_calculate_reference!(reference::Vector{Float64}, dat::DataFrame, reference_channels::AbstractVector{Symbol}) =
    _calculate_reference!(reference, dat, Int[columnindex(dat, ch) for ch in reference_channels])

@inline function _subtract_reference_col!(df::DataFrame, col_idx::Int, reference::Vector{Float64}, n::Int)
    col = df[!, col_idx]
    if col isa Vector{Float64}
        col_f64 = col::Vector{Float64}
        @inbounds @simd for i = 1:n
            col_f64[i] -= reference[i]
        end
    else
        df[!, col_idx] .-= reference
    end
    return nothing
end

function _apply_rereference!(dat::DataFrame, channel_selection::AbstractVector{Symbol}, reference_selection::AbstractVector{Symbol})
    n = n_samples(dat)
    reference = Vector{Float64}(undef, n)
    ref_indices = Int[columnindex(dat, ch) for ch in reference_selection]
    ch_indices = Int[columnindex(dat, ch) for ch in channel_selection]
    _calculate_reference!(reference, dat, ref_indices)
    if length(ch_indices) > 1 && Threads.nthreads() > 1 && length(ch_indices) * n >= 2_000_000
        Threads.@threads :static for idx in eachindex(ch_indices)
            _subtract_reference_col!(dat, ch_indices[idx], reference, n)
        end
    else
        for col_idx in ch_indices
            _subtract_reference_col!(dat, col_idx, reference, n)
        end
    end
    return nothing
end

function _apply_rereference!(
    dat::AbstractVector{DataFrame},
    channel_selection::AbstractVector{Symbol},
    reference_selection::AbstractVector{Symbol},
)
    isempty(dat) && return nothing
    first_df = dat[1]
    n = n_samples(first_df)
    ref_indices = Int[columnindex(first_df, ch) for ch in reference_selection]
    ch_indices = Int[columnindex(first_df, ch) for ch in channel_selection]

    _foreach_chunk(length(dat)) do df_range
        reference = Vector{Float64}(undef, n)
        for df_idx in df_range
            df = dat[df_idx]
            _calculate_reference!(reference, df, ref_indices)
            for col_idx in ch_indices
                _subtract_reference_col!(df, col_idx, reference, n)
            end
        end
    end
    return nothing
end

function _get_reference_channels(dat::EegData, reference_channel::Symbol)
    reference_channel == :none && return Symbol[]
    reference_channel == :avg && return channel_labels(dat)
    reference_channel == :mastoid && return [:M1, :M2]
    return [reference_channel]
end

function _get_reference_channels(dat, reference_channel::AbstractVector{Symbol})
    isempty(reference_channel) && return Symbol[]
    return length(reference_channel) == 1 ? _get_reference_channels(dat, reference_channel[1]) : reference_channel
end

_get_reference_channels(dat, ref::AbstractString) = _get_reference_channels(dat, Symbol(ref))
_get_reference_channels(dat, refs::AbstractVector{<:AbstractString}) = _get_reference_channels(dat, Symbol.(refs))

"""
    rereference!(dat::EegData, reference_selection; channel_selection = channels())
    rereference!(dat::AbstractVector{<:EegData}, reference_selection; channel_selection = channels())

Rereference EEG data in-place. `reference_selection` can be `:avg`, `:mastoid`, `:none`,
or a channel name / vector of channel names.
"""
function rereference!(dat::EegData, reference_selection; channel_selection = channels())
    reference_channels = _get_reference_channels(dat, reference_selection)
    if isempty(reference_channels)
        @info "No rereferencing applied (reference: $(reference_selection))"
        return dat
    end

    selected_channels = get_selected_channels(dat, channel_selection, include_meta = false, include_extra = false)
    missing_channels = [ch for ch in reference_channels if ch ∉ dat]
    !isempty(missing_channels) && @minimal_error("Missing reference channels in data: $(missing_channels)")

    @info "Rereference channels: $(reference_selection) ($(_print_vector(reference_channels; n_ends=3)))"
    _apply_rereference!(dat.data, selected_channels, reference_channels)
    dat.analysis_info.reference = reference_selection isa Symbol ? reference_selection : Symbol(join(reference_selection, '_'))
    return dat
end

function rereference!(dat::AbstractVector{<:EegData}, reference_selection; channel_selection = channels())
    for d in dat
        rereference!(d, reference_selection; channel_selection = channel_selection)
    end
    return dat
end

"""
    rereference(dat, reference_selection; channel_selection = channels())
    rereference(file_pattern::String; input_dir = pwd(), reference_selection = :avg,
                participant_selection = participants(), condition_selection = conditions(), output_dir = nothing)

Non-mutating or batch JLD2 version of [`rereference!`](@ref).
"""
rereference(dat::EegData, args...; kwargs...) = rereference!(copy(dat), args...; kwargs...)
rereference(dat::AbstractVector{<:EegData}, args...; kwargs...) = [rereference(d, args...; kwargs...) for d in dat]

function _default_rereference_output_dir(input_dir::String, pattern::String, ref_selection)
    ref_str = ref_selection isa Symbol ? string(ref_selection) : join(ref_selection, "_")
    return joinpath(input_dir, "rereferenced_$(clean_pattern(pattern))_$(ref_str)")
end

function _process_rereference_file(filepath::String, output_path::String, reference_selection, condition_selection::Function)
    fname = basename(filepath)
    data = read_data(filepath)
    isnothing(data) && return BatchResult(false, fname, "No data variables found")
    data isa AbstractVector{<:Union{ErpData,EpochData}} ||
        return BatchResult(false, fname, "Invalid data type: expected Vector{ErpData} or Vector{EpochData}")

    data = _condition_select(data, condition_selection)
    ref_str = reference_selection isa Symbol ? string(reference_selection) : join(reference_selection, ", ")
    if isempty(data)
        jldsave(output_path; data = data, format_version = EEGFUN_FORMAT_VERSION)
        return BatchResult(true, fname, "Rereferenced to $ref_str (empty data)")
    end

    rereference!(data, reference_selection)
    jldsave(output_path; data = data, format_version = EEGFUN_FORMAT_VERSION)
    return BatchResult(true, fname, "Rereferenced to $ref_str")
end

function rereference(
    file_pattern::String;
    input_dir::String = pwd(),
    reference_selection::Union{Symbol,AbstractVector{Symbol}} = :avg,
    participant_selection = participants(),
    condition_selection = conditions(),
    output_dir::Union{String,Nothing} = nothing,
)
    _with_batch_logging("rereferencing", "rereference.log", input_dir, output_dir) do ctx
        @log_call "rereference"
        output_dir = _prepare_output_dir!(ctx, _default_rereference_output_dir(input_dir, file_pattern, reference_selection))
        ref_str = reference_selection isa Symbol ? string(reference_selection) : join(reference_selection, ", ")
        @info "Reference settings: $ref_str"
        process_fn =
            (input_path, output_path) -> _process_rereference_file(input_path, output_path, reference_selection, condition_selection)
        batch_process(process_fn, file_pattern, input_dir, output_dir, participant_selection, "Rereferencing")
    end
end
