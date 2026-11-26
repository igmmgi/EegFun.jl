# Average, mastoid, and custom-channel EEG re-referencing operations.

######################################
########### User Functions ###########
######################################

"""
    rereference(dat, reference_selection; channel_selection = channels())
    rereference(file_pattern::String; input_dir = pwd(), reference_selection = :avg,
                participant_selection = participants(), condition_selection = conditions(), output_dir = nothing)

Non-mutating or batch JLD2 version of [`rereference!`](@ref).
"""
rereference(dat::EegData, args...; kwargs...) = rereference!(copy(dat), args...; kwargs...)

rereference(dat::AbstractVector{<:EegData}, args...; kwargs...) = [rereference(d, args...; kwargs...) for d in dat]

"""
    rereference!(dat::EegData, reference_selection; channel_selection = channels())
    rereference!(dat::AbstractVector{<:EegData}, reference_selection; channel_selection = channels())

Rereference EEG data in-place. `reference_selection` can be `:avg`, `:mastoid`, `:none`,
or a channel name / vector of channel names.
"""
function rereference!(dat::EegData, reference_selection; channel_selection = channels())
    selected_channels = get_selected_channels(dat, channel_selection, include_meta = false, include_extra = false)
    reference_channels = _get_reference_channels(dat, reference_selection; avg_channels = selected_channels)
    if isempty(reference_channels)
        @info "No rereferencing applied (reference: $(reference_selection))"
        return dat
    end

    missing_channels = [ch for ch in reference_channels if ch ∉ dat]
    !isempty(missing_channels) && @minimal_error("Missing reference channels in data: $(missing_channels)")

    @info "Rereference channels: $(reference_selection) ($(_print_vector(reference_channels; n_ends=3)))"
    _apply_rereference!(dat.data, selected_channels, reference_channels)
    dat.analysis_info.reference = reference_selection isa Symbol ? reference_selection : Symbol(join(reference_selection, '_'))
    return dat
end

rereference!(dat::AbstractVector{<:EegData}, reference_selection; channel_selection = channels()) =
    (foreach(d -> rereference!(d, reference_selection; channel_selection = channel_selection), dat); dat)

function rereference(
    file_pattern::String;
    input_dir::String = pwd(),
    reference_selection::Union{Symbol,AbstractVector{Symbol}} = :avg,
    participant_selection = participants(),
    condition_selection = conditions(),
    output_dir::Union{String,Nothing} = nothing,
)
    _with_batch_logging("rereferencing", "rereference.log", input_dir, output_dir) do ctx
        _ = output_dir
        @log_call "rereference"
        ref_dir_str = reference_selection isa Symbol ? string(reference_selection) : join(reference_selection, "_")
        ref_msg_str = reference_selection isa Symbol ? string(reference_selection) : join(reference_selection, ", ")
        @info "Reference settings: $ref_msg_str"
        default_dir = joinpath(input_dir, "rereferenced_$(_clean_pattern(file_pattern))_$(ref_dir_str)")
        return _batch_jld2_files(
            data ->
                isempty(data) ? (data, "Rereferenced to $ref_msg_str (empty data)") :
                (rereference!(data, reference_selection), "Rereferenced to $ref_msg_str"),
            AbstractVector{<:Union{ErpData,EpochData}},
            ctx,
            default_dir,
            file_pattern,
            input_dir,
            participant_selection,
            "Rereferencing";
            condition_selection = condition_selection,
            expected_label = "Vector{ErpData} or Vector{EpochData}",
        )
    end
end

######################################
######### Internal Functions #########
######################################

_calculate_reference!(reference::Vector{Float64}, dat::DataFrame, ref_indices::AbstractVector{<:Integer}) =
    _fill_colmeans_idx!(reference, dat, ref_indices)

_calculate_reference!(reference::Vector{Float64}, dat::DataFrame, reference_channels::AbstractVector{Symbol}) =
    _fill_colmeans_idx!(reference, dat, Int[columnindex(dat, ch) for ch in reference_channels])

@inline function _sub_ref_vec!(col::Vector{Float64}, reference::Vector{Float64}, n::Int)
    @inbounds @simd for i = 1:n
        col[i] -= reference[i]
    end
    return nothing
end

@inline _sub_ref_vec!(col::AbstractVector, reference::Vector{Float64}, ::Int) = (col .-= reference; nothing)

@inline _subtract_reference_col!(df::DataFrame, col_idx::Int, reference::Vector{Float64}, n::Int) =
    _sub_ref_vec!(df[!, col_idx], reference, n)

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

function _get_reference_channels(dat::EegData, reference_channel::Symbol; avg_channels::Union{Nothing,Vector{Symbol}} = nothing)
    reference_channel == :none && return Symbol[]
    reference_channel == :avg && return isnothing(avg_channels) ? channel_labels(dat) : avg_channels
    reference_channel == :mastoid && return [:M1, :M2]
    return [reference_channel]
end

function _get_reference_channels(dat, reference_channel::AbstractVector{Symbol}; avg_channels::Union{Nothing,Vector{Symbol}} = nothing)
    isempty(reference_channel) && return Symbol[]
    return length(reference_channel) == 1 ? _get_reference_channels(dat, reference_channel[1]; avg_channels = avg_channels) :
           reference_channel
end

_get_reference_channels(dat, ref::AbstractString; avg_channels::Union{Nothing,Vector{Symbol}} = nothing) =
    _get_reference_channels(dat, Symbol(ref); avg_channels = avg_channels)

_get_reference_channels(dat, refs::AbstractVector{<:AbstractString}; avg_channels::Union{Nothing,Vector{Symbol}} = nothing) =
    _get_reference_channels(dat, Symbol.(refs); avg_channels = avg_channels)
