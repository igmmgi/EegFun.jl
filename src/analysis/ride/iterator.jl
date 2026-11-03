"""
    RIDE Iteration Core

Native Julia implementation of the L1-norm / L2-norm Expectation-Maximization
loop for Residue Iteration DEcomposition (RIDE).
"""

using Statistics

"""
Construct a 1D Tukey cosine-tapered window of length `L` with taper fraction `alpha` in `[0, 1]`.
"""
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

"""
Apply baseline/edge-drift removal and Tukey boundary tapering to `vec[win_range]`,
zeroing samples outside `win_range`.
"""
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

    # Zero outside component active window
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
        # Remove linear drift anchored to left/right boundary means (prevents DC/drift leakage across RIDE components)
        n_edge = max(1, min(div(L, 4), round(Int, 0.05 * L)))
        left_mean = mean(view(vec, w_start:(w_start+n_edge-1)))
        right_mean = mean(view(vec, (w_end-n_edge+1):w_end))
        denom = T(max(1, L - 1))
        slope = (right_mean - left_mean) / denom
        @inbounds @simd for idx = 1:L
            s = w_start + idx - 1
            vec[s] -= (left_mean + slope * T(idx - 1))
        end
    end

    # Multiply by Tukey boundary taper so shifted templates smoothly vanish at boundaries
    if length(taper_weights) == L
        @inbounds @simd for idx = 1:L
            s = w_start + idx - 1
            vec[s] *= taper_weights[idx]
        end
    end

    return vec
end

"""
    ride_iteration!(
        templates::AbstractArray{T, 3},
        data::AbstractArray{T, 3},
        latencies::AbstractMatrix{Int};
        max_iter::Int = 50,
        tol::Float64 = 1e-3,
        method::Symbol = :median,
        taper_alpha::Real = 0.2,
        component_windows::Union{Nothing,Vector{UnitRange{Int}}} = nothing,
        baseline_range::Union{Nothing,UnitRange{Int}} = nothing,
    ) where {T<:Real}

Perform the RIDE iteration loop to extract component templates, applying baseline/linear-drift
stabilization and Tukey boundary tapering on each component update.

# Arguments
- `templates`: 3D array of `(samples, channels, components)`. Will be updated in-place.
- `data`: 3D array of `(samples, channels, trials)`.
- `latencies`: 2D matrix of `(trials, components)` containing the sample shifts.
- `max_iter`: Maximum number of iterations.
- `tol`: Convergence tolerance based on L1-norm of template changes.
- `method`: `:median` (L1-norm, standard RIDE) or `:mean` (L2-norm).
- `taper_alpha`: Tukey cosine taper fraction (`0.0` disables boundary tapering; default `0.2`).
- `component_windows`: Optional active sample range per component (defaults to `1:n_samples` for all).
- `baseline_range`: Optional sample range for explicit baseline subtraction per iteration.

# Returns
- `iter`: Number of iterations performed.
- `converged`: Boolean indicating if convergence was reached.
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

    # Pre-allocate channel-specific buffers outside the loop (channels run in parallel without contention)
    trial_residues = [zeros(T, n_samples) for _ = 1:n_channels]
    sort_buffers = [zeros(T, n_trials) for _ = 1:n_channels]
    aligned_buffers = [zeros(T, n_samples, n_trials) for _ = 1:n_channels]

    # Track convergence
    prev_templates = copy(templates)
    converged = false
    iter = 0

    for i = 1:max_iter
        iter = i
        max_diff = zero(T)

        for k = 1:n_components
            win_k = comp_wins[k]
            taper_k = comp_tapers[k]

            # Channels are entirely independent, perfectly parallelizable
            Threads.@threads for c = 1:n_channels

                # Use channel-specific buffers
                trial_residue = trial_residues[c]
                sort_buffer = sort_buffers[c]
                aligned_buffer = aligned_buffers[c]

                # 1. Compute residues and align them for this channel
                @inbounds for t = 1:n_trials

                    # Start with the raw data for this channel and trial
                    copyto!(trial_residue, view(data, :, c, t))

                    # Subtract all OTHER components
                    for j = 1:n_components
                        if j == k
                            continue
                        end

                        comp_j = view(templates, :, c, j)
                        shift_val = latencies[t, j]

                        # Apply shift directly via SIMD without an intermediate buffer
                        if shift_val >= n_samples || shift_val <= -n_samples
                            continue
                        elseif shift_val > 0
                            @simd for s = 1:(n_samples-shift_val)
                                trial_residue[s+shift_val] -= comp_j[s]
                            end
                        elseif shift_val < 0
                            s_val = abs(shift_val)
                            @simd for s = 1:(n_samples-s_val)
                                trial_residue[s] -= comp_j[s+s_val]
                            end
                        else
                            @simd for s = 1:n_samples
                                trial_residue[s] -= comp_j[s]
                            end
                        end
                    end

                    # Now we have the residue. Align it by shifting backwards by latencies[t, k]
                    align_val = -latencies[t, k]

                    # Apply alignment directly to aligned_buffer
                    if align_val >= n_samples || align_val <= -n_samples
                        for s = 1:n_samples
                            aligned_buffer[s, t] = zero(T)
                        end
                    elseif align_val > 0
                        for s = 1:align_val
                            aligned_buffer[s, t] = zero(T)
                        end
                        @simd for s = 1:(n_samples-align_val)
                            aligned_buffer[s+align_val, t] = trial_residue[s]
                        end
                    elseif align_val < 0
                        s_val = abs(align_val)
                        @simd for s = 1:(n_samples-s_val)
                            aligned_buffer[s, t] = trial_residue[s+s_val]
                        end
                        for s = (n_samples-s_val+1):n_samples
                            aligned_buffer[s, t] = zero(T)
                        end
                    else
                        for s = 1:n_samples
                            aligned_buffer[s, t] = trial_residue[s]
                        end
                    end
                end

                # 2. Update the template for component k, channel c
                if method == :median
                    @inbounds for s = 1:n_samples
                        @simd for t = 1:n_trials
                            sort_buffer[t] = aligned_buffer[s, t]
                        end
                        templates[s, c, k] = median!(sort_buffer)
                    end
                else
                    @inbounds for s = 1:n_samples
                        val = zero(T)
                        @simd for t = 1:n_trials
                            val += aligned_buffer[s, t]
                        end
                        templates[s, c, k] = val / n_trials
                    end
                end

                # Apply baseline/linear-drift stabilization and Tukey boundary tapering
                _apply_ride_boundary_taper!(
                    view(templates, :, c, k),
                    win_k,
                    taper_k;
                    baseline_range = baseline_range,
                    detrend_edges = detrend_multi,
                )
            end
        end

        # 3. Check convergence
        for k = 1:n_components
            diff = sum(abs, view(templates, :, :, k) .- view(prev_templates, :, :, k))
            max_diff = max(max_diff, diff)
        end

        # Normalize diff by total absolute sum to make tolerance scale-invariant
        total_abs = sum(abs, prev_templates)
        if total_abs > 0
            max_diff = max_diff / total_abs
        end

        if max_diff < tol
            converged = true
            break
        end

        copyto!(prev_templates, templates)
    end

    return iter, converged
end
