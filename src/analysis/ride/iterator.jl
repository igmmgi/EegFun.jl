"""
    RIDE Iteration Core

Native Julia implementation of the L1-norm / L2-norm Expectation-Maximization
loop for Residue Iteration DEcomposition (RIDE).
"""

using Statistics

"""
    ride_iteration!(
        templates::AbstractArray{T, 3},
        data::AbstractArray{T, 3},
        latencies::AbstractMatrix{Int};
        max_iter::Int = 50,
        tol::Float64 = 1e-3,
        method::Symbol = :median
    ) where {T<:Real}

Perform the RIDE iteration loop to extract component templates.

# Arguments
- `templates`: 3D array of `(samples, channels, components)`. Will be updated in-place.
- `data`: 3D array of `(samples, channels, trials)`.
- `latencies`: 2D matrix of `(trials, components)` containing the sample shifts.
- `max_iter`: Maximum number of iterations.
- `tol`: Convergence tolerance based on L1-norm of template changes.
- `method`: `:median` (L1-norm, standard RIDE) or `:mean` (L2-norm).

# Returns
- `iter`: Number of iterations performed.
- `converged`: Boolean indicating if convergence was reached.
"""
function ride_iteration!(
    templates::AbstractArray{T, 3},
    data::AbstractArray{T, 3},
    latencies::AbstractMatrix{Int};
    max_iter::Int = 50,
    tol::Float64 = 1e-3,
    method::Symbol = :median
) where {T<:Real}

    n_samples, n_channels, n_trials = size(data)
    _, _, n_components = size(templates)

    # Pre-allocate thread-local buffers per channel outside the loop
    trial_residues = [zeros(T, n_samples) for _ in 1:n_channels]
    sort_buffers   = [zeros(T, n_trials) for _ in 1:n_channels]
    aligned_buffers = [zeros(T, n_samples, n_trials) for _ in 1:n_channels]
    
    # Track convergence
    prev_templates = copy(templates)
    converged = false
    iter = 0

    for i in 1:max_iter
        iter = i
        max_diff = zero(T)

        for k in 1:n_components
            
            # Channels are entirely independent, perfectly parallelizable
            Threads.@threads for c in 1:n_channels
                
                # Use thread-local buffers for this specific channel
                trial_residue = trial_residues[c]
                sort_buffer = sort_buffers[c]
                aligned_buffer = aligned_buffers[c]
                
                # 1. Compute residues and align them for this channel
                @inbounds for t in 1:n_trials
                    
                    # Start with the raw data for this channel and trial
                    copyto!(trial_residue, view(data, :, c, t))
                    
                    # Subtract all OTHER components
                    for j in 1:n_components
                        if j == k
                            continue
                        end
                        
                        comp_j = view(templates, :, c, j)
                        shift_val = latencies[t, j]
                        
                        # Apply shift directly via SIMD without an intermediate buffer
                        if shift_val >= n_samples || shift_val <= -n_samples
                            continue
                        elseif shift_val > 0
                            @simd for s in 1:(n_samples - shift_val)
                                trial_residue[s + shift_val] -= comp_j[s]
                            end
                        elseif shift_val < 0
                            s_val = abs(shift_val)
                            @simd for s in 1:(n_samples - s_val)
                                trial_residue[s] -= comp_j[s + s_val]
                            end
                        else
                            @simd for s in 1:n_samples
                                trial_residue[s] -= comp_j[s]
                            end
                        end
                    end
                    
                    # Now we have the residue. Align it by shifting backwards by latencies[t, k]
                    align_val = -latencies[t, k]
                    
                    # Apply alignment directly to aligned_buffer
                    if align_val >= n_samples || align_val <= -n_samples
                        for s in 1:n_samples
                            aligned_buffer[s, t] = zero(T)
                        end
                    elseif align_val > 0
                        for s in 1:align_val
                            aligned_buffer[s, t] = zero(T)
                        end
                        @simd for s in 1:(n_samples - align_val)
                            aligned_buffer[s + align_val, t] = trial_residue[s]
                        end
                    elseif align_val < 0
                        s_val = abs(align_val)
                        @simd for s in 1:(n_samples - s_val)
                            aligned_buffer[s, t] = trial_residue[s + s_val]
                        end
                        for s in (n_samples - s_val + 1):n_samples
                            aligned_buffer[s, t] = zero(T)
                        end
                    else
                        for s in 1:n_samples
                            aligned_buffer[s, t] = trial_residue[s]
                        end
                    end
                end
                
                # 2. Update the template for component k, channel c
                if method == :median
                    # Zero-allocation median: copy to sort_buffer and mutate
                    @inbounds for s in 1:n_samples
                        @simd for t in 1:n_trials
                            sort_buffer[t] = aligned_buffer[s, t]
                        end
                        templates[s, c, k] = median!(sort_buffer)
                    end
                elseif method == :mean
                    @inbounds for s in 1:n_samples
                        val = zero(T)
                        @simd for t in 1:n_trials
                            val += aligned_buffer[s, t]
                        end
                        templates[s, c, k] = val / n_trials
                    end
                else
                    error("Unknown method: $method. Use :median or :mean")
                end
                
                # Apply Tukey window/detrending here if needed (as per MATLAB implementation)
                # For this core engine, we keep it pure math. Windows can be applied 
                # externally or added as a post-iteration step.
            end
        end

        # 3. Check convergence
        for k in 1:n_components
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
