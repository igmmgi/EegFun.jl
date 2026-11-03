struct _RideIterWorkspace{T<:Real}
    trial_residues::Vector{Vector{T}}
    sort_buffers::Vector{Vector{T}}
    aligned_buffers::Vector{Matrix{T}}
    prev_templates::Array{T,3}
end

function _RideIterWorkspace(::Type{T}, n_samples::Int, n_channels::Int, n_trials::Int, n_components::Int, method::Symbol) where {T<:Real}
    n_workers = max(1, min(Threads.nthreads(), n_channels))
    trial_residues = [Vector{T}(undef, n_samples) for _ = 1:n_workers]
    sort_buffers = method === :median ? [Vector{T}(undef, n_trials) for _ = 1:n_workers] : Vector{T}[]
    aligned_buffers = method === :median ? [Matrix{T}(undef, n_samples, n_trials) for _ = 1:n_workers] : Matrix{T}[]
    prev_templates = Array{T,3}(undef, n_samples, n_channels, n_components)
    return _RideIterWorkspace{T}(trial_residues, sort_buffers, aligned_buffers, prev_templates)
end

function _ride_tukey_window(::Type{T}, L::Int, alpha::Real) where {T<:Real}
    L <= 0 && return T[]
    w = ones(T, L)
    (L <= 2 || alpha <= 0) && return w
    a = clamp(Float64(alpha), 0.0, 1.0)
    m = round(Int, a * (L - 1) / 2.0)
    m <= 0 && return w
    @inbounds for i = 0:(m-1)
        val = T(0.5 * (1.0 - cos(π * i / m)))
        w[i+1] = val
        w[L-i] = val
    end
    return w
end

function _apply_ride_boundary_taper!(
    vec::AbstractVector{T},
    win_range::UnitRange{Int},
    taper_weights::Vector{T};
    baseline_range::Union{Nothing,UnitRange{Int}} = nothing,
    detrend_edges::Bool = false,
) where {T<:Real}
    n = length(vec)
    w_start = first(win_range)
    w_end = last(win_range)
    L = length(win_range)

    @inbounds for s = 1:(w_start-1)
        vec[s] = zero(T)
    end
    @inbounds for s = (w_end+1):n
        vec[s] = zero(T)
    end
    L == 0 && return vec

    if !isnothing(baseline_range) && !isempty(baseline_range)
        b_start = clamp(first(baseline_range), 1, n)
        b_end = clamp(last(baseline_range), b_start, n)
        base_val = mean(view(vec, b_start:b_end))
        @inbounds @simd for s = w_start:w_end
            vec[s] -= base_val
        end
    elseif detrend_edges && L >= 4
        n_edge = max(1, min(div(L, 4), round(Int, 0.05 * L)))
        left_mean = mean(view(vec, w_start:(w_start+n_edge-1)))
        right_mean = mean(view(vec, (w_end-n_edge+1):w_end))
        slope = (right_mean - left_mean) / T(max(1, L - 1))
        @inbounds @simd for idx = 1:L
            s = w_start + idx - 1
            vec[s] -= left_mean + slope * T(idx - 1)
        end
    end

    if length(taper_weights) == L
        @inbounds @simd for idx = 1:L
            vec[w_start+idx-1] *= taper_weights[idx]
        end
    end
    return vec
end

@inline function _add_shifted!(dest::AbstractVector{T}, src::AbstractVector{T}, shift::Int) where {T}
    n = length(dest)
    abs(shift) >= n && return dest
    if shift > 0
        @inbounds @simd for s = 1:(n-shift)
            dest[s+shift] += src[s]
        end
    elseif shift < 0
        s_val = -shift
        @inbounds @simd for s = 1:(n-s_val)
            dest[s] += src[s+s_val]
        end
    else
        @inbounds @simd for s = 1:n
            dest[s] += src[s]
        end
    end
    return dest
end

"""
    ride_iteration!(templates, data, latencies; max_iter=50, tol=1e-3, method=:median,
                    taper_alpha=0.2, component_windows=nothing, baseline_range=nothing)

Run the inner L1 (`:median`) or L2 (`:mean`) RIDE template estimation loop in-place with
boundary drift stabilization and Tukey tapering, returning `(iterations, converged)`.
"""
function ride_iteration!(
    templates::AbstractArray{T,3},
    data::AbstractArray{T,3},
    latencies::AbstractMatrix{Int};
    max_iter::Int = 50,
    tol::Real = 1e-3,
    method::Symbol = :median,
    taper_alpha::Real = 0.2,
    component_windows::Union{Nothing,Vector{UnitRange{Int}}} = nothing,
    baseline_range::Union{Nothing,UnitRange{Int}} = nothing,
    _workspace::Union{Nothing,_RideIterWorkspace{T}} = nothing,
) where {T<:Real}
    method in (:median, :mean) || @minimal_error("Unknown method: $method. Use :median or :mean")

    n_samples, n_channels, n_trials = size(data)
    _, _, n_components = size(templates)
    (n_samples == 0 || n_channels == 0 || n_trials == 0 || n_components == 0) && return (0, true)

    comp_wins = if isnothing(component_windows)
        [1:n_samples for _ = 1:n_components]
    else
        length(component_windows) == n_components ||
            @minimal_error("component_windows length ($(length(component_windows))) must match n_components ($n_components)")
        [clamp(first(w), 1, n_samples):clamp(last(w), clamp(first(w), 1, n_samples), n_samples) for w in component_windows]
    end
    comp_tapers = [_ride_tukey_window(T, length(w), taper_alpha) for w in comp_wins]
    detrend_multi = n_components > 1 && isnothing(baseline_range)

    ws = isnothing(_workspace) ? _RideIterWorkspace(T, n_samples, n_channels, n_trials, n_components, method) : _workspace
    prev_templates = ws.prev_templates
    copyto!(prev_templates, templates)

    n_workers = length(ws.trial_residues)
    ch_chunks = collect(Iterators.partition(1:n_channels, cld(n_channels, n_workers)))
    n_chunks = length(ch_chunks)
    inv_trials = one(T) / T(n_trials)
    converged = false
    iter = 0

    for i = 1:max_iter
        iter = i

        for k = 1:n_components
            win_k = comp_wins[k]
            w_start, w_end = first(win_k), last(win_k)
            taper_k = comp_tapers[k]

            _foreach_chunk(n_chunks; ntasks = n_chunks) do chunk_ids
                for chunk_id in chunk_ids
                    ch_range = ch_chunks[chunk_id]
                    trial_residue = ws.trial_residues[chunk_id]

                    if method === :median
                        sort_buffer = ws.sort_buffers[chunk_id]
                        aligned_buffer = ws.aligned_buffers[chunk_id]
                        for c in ch_range
                            @inbounds for t = 1:n_trials
                                copyto!(trial_residue, view(data, :, c, t))
                                for j = 1:n_components
                                    j == k && continue
                                    _subtract_shifted!(trial_residue, view(templates, :, c, j), latencies[t, j])
                                end
                                shift_signal!(view(aligned_buffer, :, t), trial_residue, -latencies[t, k])
                            end

                            @inbounds for s = 1:(w_start-1)
                                templates[s, c, k] = zero(T)
                            end
                            @inbounds for s = w_start:w_end
                                @simd for t = 1:n_trials
                                    sort_buffer[t] = aligned_buffer[s, t]
                                end
                                templates[s, c, k] = median!(sort_buffer)
                            end
                            @inbounds for s = (w_end+1):n_samples
                                templates[s, c, k] = zero(T)
                            end

                            _apply_ride_boundary_taper!(
                                view(templates, :, c, k),
                                win_k,
                                taper_k;
                                baseline_range = baseline_range,
                                detrend_edges = detrend_multi,
                            )
                        end
                    else
                        for c in ch_range
                            tmpl_ck = view(templates, :, c, k)
                            fill!(tmpl_ck, zero(T))
                            @inbounds for t = 1:n_trials
                                copyto!(trial_residue, view(data, :, c, t))
                                for j = 1:n_components
                                    j == k && continue
                                    _subtract_shifted!(trial_residue, view(templates, :, c, j), latencies[t, j])
                                end
                                _add_shifted!(tmpl_ck, trial_residue, -latencies[t, k])
                            end
                            @inbounds @simd for s = w_start:w_end
                                tmpl_ck[s] *= inv_trials
                            end
                            _apply_ride_boundary_taper!(
                                tmpl_ck,
                                win_k,
                                taper_k;
                                baseline_range = baseline_range,
                                detrend_edges = detrend_multi,
                            )
                        end
                    end
                end
            end
        end

        max_diff = zero(T)
        total_abs = zero(T)
        @inbounds for k = 1:n_components
            diff_k = zero(T)
            for c = 1:n_channels
                @simd for s = 1:n_samples
                    prev_v = prev_templates[s, c, k]
                    diff_k += abs(templates[s, c, k] - prev_v)
                    total_abs += abs(prev_v)
                end
            end
            diff_k > max_diff && (max_diff = diff_k)
        end

        total_abs > 0 && (max_diff /= total_abs)
        if max_diff < tol
            converged = true
            break
        end
        copyto!(prev_templates, templates)
    end

    return iter, converged
end
