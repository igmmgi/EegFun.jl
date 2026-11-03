_to_repair_channels(ch::Symbol) = [ch]
_to_repair_channels(ch::AbstractString) = [Symbol(ch)]
_to_repair_channels(chans::AbstractVector{Symbol}) = collect(chans)
_to_repair_channels(chans::AbstractVector{<:AbstractString}) = Symbol.(chans)
_to_repair_channels(chans) = collect(chans)

"""
    repair_channels!(data, channels_to_repair; method::Symbol=:neighbor_interpolation, kwargs...)

Repair bad EEG channels in `ContinuousData`, `ErpData`, `EpochData`, or an `AbstractMatrix` in-place.

# Keyword Arguments
- `method::Symbol`: `:neighbor_interpolation` (default; uses `neighbors_dict` or `data.layout.neighbors`)
  or `:spherical_spline` (accepts `m::Int=4`, `lambda::Real=1e-5`).
- `epoch_selection`: Epoch selector for `EpochData` (default: `epochs()`).
"""
function repair_channels!(data, channels_to_repair; method::Symbol = :neighbor_interpolation, repair_info = nothing, kwargs...)
    chs = _to_repair_channels(channels_to_repair)
    actually_repaired = if method == :neighbor_interpolation
        _repair_channels_neighbor!(data, chs; repair_info = repair_info, kwargs...)
    elseif method == :spherical_spline
        data isa EpochData ? _repair_channels_spherical!(data, chs; kwargs...) :
        _repair_channels_spherical!(data, chs; repair_info = repair_info, kwargs...)
    else
        @minimal_error("Unknown repair method: $method. Use :neighbor_interpolation or :spherical_spline")
    end

    if !(data isa EpochData) && hasfield(typeof(data), :analysis_info) && !isnothing(actually_repaired)
        for ch in actually_repaired
            ch ∉ data.analysis_info.repaired_channels && push!(data.analysis_info.repaired_channels, ch)
        end
        data.analysis_info.n_channels_repaired = length(data.analysis_info.repaired_channels)
    end
    return data
end

"""
    repair_channels(dat, args...; kwargs...)

Non-mutating version of [`repair_channels!`](@ref).
"""
repair_channels(dat, args...; kwargs...) = repair_channels!(copy(dat), args...; kwargs...)

function repair_channels!(dat::AbstractVector{<:EegData}, args...; kwargs...)
    for d in dat
        repair_channels!(d, args...; kwargs...)
    end
    return dat
end

repair_channels(dat::AbstractVector{<:EegData}, args...; kwargs...) = [repair_channels(d, args...; kwargs...) for d in dat]

# --- Neighbor interpolation ---

function _repair_channels_neighbor!(data, channels_to_repair; neighbors_dict = nothing, repair_info = nothing, kwargs...)
    if isnothing(neighbors_dict)
        hasfield(typeof(data), :layout) || @minimal_error("neighbors_dict must be provided for AbstractMatrix data")
        isnothing(data.layout.neighbors) && get_neighbors_xyz!(data.layout, 0.5)
        neighbors_dict = data.layout.neighbors
    end
    return data isa EpochData ? _repair_channels_neighbor!(data, channels_to_repair, neighbors_dict; kwargs...) :
           _repair_channels_neighbor!(data, channels_to_repair, neighbors_dict; repair_info = repair_info, kwargs...)
end

function _resolve_neighbor_plans(
    channels::AbstractVector{Symbol},
    channels_to_repair::AbstractVector{Symbol},
    neighbors_dict::AbstractDict{Symbol,Neighbors};
    repair_info = nothing,
    bad_channels::Union{Nothing,AbstractVector{Symbol}} = nothing,
    log_repairs::Bool = true,
)
    if !isnothing(repair_info)
        repair_info.repaired = Symbol[]
        repair_info.skipped = Symbol[]
    end

    plans = Tuple{Int,Symbol,Vector{Int},Vector{Symbol},Vector{Float64}}[]
    unique_to_repair = unique(channels_to_repair)
    isempty(unique_to_repair) && return plans

    excluded = Set{Symbol}(unique_to_repair)
    !isnothing(bad_channels) && union!(excluded, bad_channels)
    ch_indices = Dict{Symbol,Int}(ch => i for (i, ch) in enumerate(channels))

    for bad_ch in unique_to_repair
        bad_idx = get(ch_indices, bad_ch, nothing)
        if isnothing(bad_idx)
            log_repairs && @minimal_warning "Channel $bad_ch not found in channel list, skipping"
            !isnothing(repair_info) && push!(repair_info.skipped, bad_ch)
            continue
        end

        nb = get(neighbors_dict, bad_ch, nothing)
        if isnothing(nb) || isempty(nb.channels)
            log_repairs && @minimal_warning "No neighbors found for channel $bad_ch, skipping"
            !isnothing(repair_info) && push!(repair_info.skipped, bad_ch)
            continue
        end

        clean_chs = Symbol[]
        clean_idxs = Int[]
        clean_weights = Float64[]
        for (i, nch) in enumerate(nb.channels)
            nch ∈ excluded && continue
            nidx = get(ch_indices, nch, nothing)
            if !isnothing(nidx)
                push!(clean_chs, nch)
                push!(clean_idxs, nidx)
                push!(clean_weights, Float64(nb.weights[i]))
            end
        end

        if isempty(clean_idxs)
            log_repairs && @minimal_warning "No clean neighbors available for channel $bad_ch, skipping"
            !isnothing(repair_info) && push!(repair_info.skipped, bad_ch)
            continue
        end

        wsum = sum(clean_weights)
        wsum > 0.0 ? (clean_weights ./= wsum) : fill!(clean_weights, 1.0 / length(clean_weights))

        push!(plans, (bad_idx, bad_ch, clean_idxs, clean_chs, clean_weights))
        !isnothing(repair_info) && push!(repair_info.repaired, bad_ch)
        log_repairs && @info "Repaired channel $bad_ch using weighted neighbor interpolation with neighbors: $(clean_chs)"
    end
    return plans
end

@inline function _scale_col!(dst::Vector{Float64}, src::AbstractVector{<:Real}, w::Float64)
    @inbounds @simd for i in eachindex(dst)
        dst[i] = Float64(src[i]) * w
    end
    return dst
end

@inline function _muladd_col!(dst::Vector{Float64}, src::AbstractVector{<:Real}, w::Float64)
    @inbounds @simd for i in eachindex(dst)
        dst[i] = muladd(Float64(src[i]), w, dst[i])
    end
    return dst
end

function _apply_neighbor_plans!(df::DataFrame, plans::Vector{Tuple{Int,Symbol,Vector{Int},Vector{Symbol},Vector{Float64}}})
    isempty(plans) && return nothing
    n = nrow(df)
    new_cols = Vector{Vector{Float64}}(undef, length(plans))
    for (p, (_, _, _, clean_chs, weights)) in enumerate(plans)
        acc = Vector{Float64}(undef, n)
        _scale_col!(acc, df[!, clean_chs[1]], weights[1])
        for k = 2:length(clean_chs)
            _muladd_col!(acc, df[!, clean_chs[k]], weights[k])
        end
        new_cols[p] = acc
    end
    for (p, (_, bad_ch, _, _, _)) in enumerate(plans)
        df[!, bad_ch] = new_cols[p]
    end
    return nothing
end

function _repair_channels_neighbor!(
    data::AbstractMatrix,
    channels::AbstractVector{Symbol},
    channels_to_repair::AbstractVector{Symbol},
    neighbors_dict::AbstractDict{Symbol,Neighbors};
    repair_info = nothing,
    bad_channels::Union{Nothing,AbstractVector{Symbol}} = nothing,
)
    size(data, 1) == length(channels) ||
        @minimal_error("Number of channels ($(length(channels))) doesn't match data dimensions ($(size(data, 1)))")

    plans = _resolve_neighbor_plans(channels, channels_to_repair, neighbors_dict; repair_info = repair_info, bad_channels = bad_channels)
    isempty(plans) && return Symbol[]

    n_tp = size(data, 2)
    temp_traces = Matrix{Float64}(undef, length(plans), n_tp)
    for (slot, (_, _, clean_idxs, _, weights)) in enumerate(plans)
        @inbounds for t = 1:n_tp
            acc = 0.0
            for k in eachindex(clean_idxs)
                acc = muladd(Float64(data[clean_idxs[k], t]), weights[k], acc)
            end
            temp_traces[slot, t] = acc
        end
    end
    for (slot, (bad_idx, _, _, _, _)) in enumerate(plans)
        @inbounds for t = 1:n_tp
            data[bad_idx, t] = temp_traces[slot, t]
        end
    end
    return [bad_ch for (_, bad_ch, _, _, _) in plans]
end

function _repair_channels_neighbor!(
    data::Union{ContinuousData,ErpData},
    channels_to_repair::AbstractVector{Symbol},
    neighbors_dict::AbstractDict{Symbol,Neighbors};
    repair_info = nothing,
    bad_channels::Union{Nothing,AbstractVector{Symbol}} = nothing,
)
    channels = Symbol.(data.layout.data.label)
    plans = _resolve_neighbor_plans(channels, channels_to_repair, neighbors_dict; repair_info = repair_info, bad_channels = bad_channels)
    _apply_neighbor_plans!(data.data, plans)
    return [bad_ch for (_, bad_ch, _, _, _) in plans]
end

function _repair_channels_neighbor!(
    data::EpochData,
    channels_to_repair::AbstractVector{Symbol},
    neighbors_dict::AbstractDict{Symbol,Neighbors};
    epoch_selection = epochs(),
    bad_channels::Union{Nothing,AbstractVector{Symbol}} = nothing,
)
    selected_epochs = get_selected_epochs(data, epoch_selection)
    isempty(selected_epochs) && return Symbol[]
    channels = Symbol.(data.layout.data.label)
    plans = _resolve_neighbor_plans(channels, channels_to_repair, neighbors_dict; bad_channels = bad_channels)
    for epoch_idx in selected_epochs
        _apply_neighbor_plans!(data.data[epoch_idx], plans)
    end
    return [bad_ch for (_, bad_ch, _, _, _) in plans]
end

# --- Spherical spline interpolation (Perrin et al., 1989) ---

_g_factors(m::Int, nterms::Int) = [(2n + 1) / (float(n)^m * float(n + 1)^m * 4π) for n = 1:nterms]

@inline function _legendre_series(x::Float64, f::Vector{Float64})
    x = clamp(x, -1.0, 1.0)
    p_prev, p = 1.0, x
    acc = f[1] * x
    @inbounds for n = 2:length(f)
        p_next = ((2n - 1) * x * p - (n - 1) * p_prev) / n
        acc += f[n] * p_next
        p_prev, p = p, p_next
    end
    return acc
end

function _g_from_cos!(C::Matrix{Float64}, f::Vector{Float64})
    @inbounds for i in eachindex(C)
        C[i] = _legendre_series(C[i], f)
    end
    return C
end

function _spline_geometry(d, ch::Vector{Int})
    ρe = hypot.(d.x2[ch], d.y2[ch])
    if hasproperty(d, :inc) && hasproperty(d, :azi)
        rad_factor = π / 180.0
        inc = Float64.(d.inc[ch]) .* rad_factor
        azi = Float64.(d.azi[ch]) .* rad_factor
        U = hcat(sin.(inc) .* cos.(azi), sin.(inc) .* sin.(azi), cos.(inc))
        θ = abs.(inc)
        denom = dot(ρe, ρe)
        k = denom > 0 ? dot(ρe, θ) / denom : (π / 2)
        if maximum(abs, θ .- k .* ρe) > 1e-3
            @warn "x2/y2 are not an azimuthal-equidistant image of inc/azi; spline map may be distorted." maxlog = 1
        end
    else
        @warn "Layout has no :inc/:azi; assuming azimuthal-equidistant x2/y2 with ρ=1 at the equator." maxlog = 1
        k = π / 2
        θ = min.(k .* ρe, π)
        sx = [ρe[i] > 0 ? sin(θ[i]) / ρe[i] : 0.0 for i in eachindex(ρe)]
        U = hcat(sx .* d.x2[ch], sx .* d.y2[ch], cos.(θ))
    end
    return U, k
end

function _spline_good_to_bad(U_good::Matrix{Float64}, U_bad::Matrix{Float64}, f::Vector{Float64}, λ::Real)
    n = size(U_good, 1)
    A = zeros(Float64, n + 1, n + 1)
    A[1:n, 1:n] .= _g_from_cos!(U_good * U_good', f)
    @inbounds for i = 1:n
        A[i, i] += λ
        A[i, n+1] = 1.0
        A[n+1, i] = 1.0
    end
    X = lu(A) \ [Matrix{Float64}(I, n, n); zeros(Float64, 1, n)]
    Gb = _g_from_cos!(U_bad * U_good', f)
    return Gb * view(X, 1:n, :) .+ view(X, (n+1):(n+1), :)
end

_calc_g(cosang::AbstractArray, m::Int = 4, n_legendre_terms::Int = 50) = _g_from_cos!(Float64.(cosang), _g_factors(m, n_legendre_terms))

function _repair_dfs_spherical!(
    dfs::AbstractVector{DataFrame},
    layout::Layout,
    channels_to_repair::AbstractVector{Symbol};
    m::Int = 4,
    lambda::Real = 1e-5,
    n_legendre_terms::Int = 50,
    repair_info = nothing,
)
    if !isnothing(repair_info)
        repair_info.repaired = Symbol[]
        repair_info.skipped = Symbol[]
    end
    isempty(dfs) && return Symbol[]

    _ensure_coordinates_3d!(layout)
    channels = Symbol.(layout.data.label)
    unique_to_repair = unique(channels_to_repair)
    isempty(unique_to_repair) && return Symbol[]

    ch_map = Dict{Symbol,Int}(ch => i for (i, ch) in enumerate(channels))
    bad_indices = Int[]
    found_bad = Symbol[]
    for bad_ch in unique_to_repair
        bad_idx = get(ch_map, bad_ch, nothing)
        if isnothing(bad_idx)
            @minimal_warning "Channel $bad_ch not found in channel list, skipping"
            !isnothing(repair_info) && push!(repair_info.skipped, bad_ch)
            continue
        end
        push!(bad_indices, bad_idx)
        push!(found_bad, bad_ch)
    end
    isempty(bad_indices) && return Symbol[]

    good_indices = setdiff(1:length(channels), bad_indices)
    if isempty(good_indices)
        @minimal_warning "Cannot repair channels via spherical spline: no clean channels remain, skipping"
        !isnothing(repair_info) && append!(repair_info.skipped, found_bad)
        return Symbol[]
    end

    coords = _spherical_coords(layout, channels; warn_zero = false)
    op = _spline_good_to_bad(coords[good_indices, :], coords[bad_indices, :], _g_factors(m, n_legendre_terms), lambda)
    op_T = collect(transpose(op))

    good_chs = channels[good_indices]
    n_good = length(good_chs)
    n_bad = length(found_bad)
    X_good = Matrix{Float64}(undef, nrow(dfs[1]), n_good)
    Y_bad = Matrix{Float64}(undef, nrow(dfs[1]), n_bad)

    for df in dfs
        n_samp = nrow(df)
        if size(X_good, 1) != n_samp
            X_good = Matrix{Float64}(undef, n_samp, n_good)
            Y_bad = Matrix{Float64}(undef, n_samp, n_bad)
        end
        for (k, gch) in enumerate(good_chs)
            copyto!(view(X_good, :, k), df[!, gch])
        end
        mul!(Y_bad, X_good, op_T)
        for (b, bch) in enumerate(found_bad)
            df[!, bch] = Y_bad[:, b]
        end
    end

    for bad_ch in found_bad
        @info "Repaired channel $bad_ch using spherical spline interpolation"
        !isnothing(repair_info) && push!(repair_info.repaired, bad_ch)
    end
    return found_bad
end

_repair_channels_spherical!(
    data::Union{ContinuousData,ErpData},
    channels_to_repair::AbstractVector{Symbol};
    m::Int = 4,
    lambda::Real = 1e-5,
    repair_info = nothing,
) = _repair_dfs_spherical!([data.data], data.layout, channels_to_repair; m = m, lambda = lambda, repair_info = repair_info)

function _repair_channels_spherical!(
    data::EpochData,
    channels_to_repair::AbstractVector{Symbol};
    epoch_selection = epochs(),
    m::Int = 4,
    lambda::Real = 1e-5,
)
    selected_epochs = get_selected_epochs(data, epoch_selection)
    isempty(selected_epochs) && return Symbol[]
    return _repair_dfs_spherical!(view(data.data, selected_epochs), data.layout, channels_to_repair; m = m, lambda = lambda)
end
