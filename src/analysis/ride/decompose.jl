"""
    RIDE Decomposition Orchestrator

High-level interface to run the Residue Iteration DEcomposition (RIDE) on EEG data.
"""

export RideComponent, ride_decompose

"""
    RideComponent

Defines a single RIDE component to be extracted.

# Fields
- `name::Symbol`: Component name (e.g., :S, :R, :C).
- `latencies::Union{Vector{Int}, Nothing}`: Pre-defined sample latencies (relative to epoch zero) 
  for each trial. For Stimulus-locked, this is usually a vector of zeros. For Response-locked, 
  this is a vector of reaction times (in samples). If `nothing`, latencies will be estimated via Woody filter.
- `window_radius::Int`: Search window radius (in samples) for the Woody filter if latencies are estimated.
"""
struct RideComponent
    name::Symbol
    latencies::Union{Vector{Int},Nothing}
    window_radius::Int
end

# Default constructor for pre-defined components
RideComponent(name::Symbol, latencies::Vector{Int}) = RideComponent(name, latencies, 0)
# Default constructor for estimated components
RideComponent(name::Symbol, window_radius::Int) = RideComponent(name, nothing, window_radius)


"""
    RideResult{T<:Real}

Contains the decomposed templates, final latencies, and residues.
"""
struct RideResult{T<:Real}
    components::Vector{Symbol}
    templates::Array{T,3}      # (samples, channels, components)
    latencies::Matrix{Int}      # (trials, components)
    residue::Array{T,3}        # (samples, channels, trials)
    converged::Bool
    iterations::Int
end


"""
    ride_decompose(
        data::AbstractArray{T, 3}, 
        comps::Vector{RideComponent};
        outer_iters::Int = 4,
        inner_iters::Int = 50,
        inner_tol::Float64 = 1e-3,
        method::Symbol = :median
    ) where {T<:Real}

Execute the full RIDE decomposition pipeline on a 3D data tensor `(samples, channels, trials)`.
"""
function ride_decompose(
    data::AbstractArray{T,3},
    comps::Vector{RideComponent};
    outer_iters::Int = 4,
    inner_iters::Int = 50,
    inner_tol::Float64 = 1e-3,
    method::Symbol = :median,
) where {T<:Real}

    n_samples, n_channels, n_trials = size(data)
    n_components = length(comps)

    # Initialize latencies matrix (trials, components)
    latencies = zeros(Int, n_trials, n_components)

    # Track which components need estimation
    needs_estimation = falses(n_components)

    for k = 1:n_components
        if !isnothing(comps[k].latencies)
            if length(comps[k].latencies) != n_trials
                error(
                    "Component $(comps[k].name) latencies length ($(length(comps[k].latencies))) does not match number of trials ($n_trials)",
                )
            end
            latencies[:, k] .= comps[k].latencies
        else
            needs_estimation[k] = true
        end
    end

    # Initialize templates
    templates = zeros(T, n_samples, n_channels, n_components)

    converged_outer = false
    actual_outer_iters = 0

    prev_latencies = copy(latencies)

    # Outer loop: Alternate between EM template extraction and Woody latency estimation
    for outer = 1:outer_iters
        actual_outer_iters = outer

        # 1. EM Template Iteration
        # Modifies `templates` in-place
        _, inner_converged = ride_iteration!(templates, data, latencies; max_iter = inner_iters, tol = inner_tol, method = method)

        # 2. Re-estimate Latencies for unknown components
        if any(needs_estimation)
            for k = 1:n_components
                if needs_estimation[k]
                    # Compute data for Woody filter: Raw data minus all OTHER components
                    # (This prevents stimulus or response components from confusing the C-component estimator)
                    woody_data = copy(data)
                    shifted_comp = zeros(T, n_samples)

                    for j = 1:n_components
                        if j == k
                            continue
                        end

                        comp_j = view(templates, :, :, j)
                        @inbounds for t = 1:n_trials
                            shift_val = latencies[t, j]
                            for c = 1:n_channels
                                shift_signal!(shifted_comp, view(comp_j, :, c), shift_val)
                                for s = 1:n_samples
                                    woody_data[s, c, t] -= shifted_comp[s]
                                end
                            end
                        end
                    end

                    # Run Woody filter using the current template for this component
                    # And use current latencies as initial guesses
                    current_template = view(templates, :, :, k)
                    initial_lats = latencies[:, k]

                    new_lats =
                        woody_filter(woody_data, comps[k].window_radius; template = current_template, initial_latencies = initial_lats)

                    latencies[:, k] .= new_lats
                end
            end

            # Check latency convergence
            # If 99% of trials haven't changed latency, we converge
            changed_trials = sum(latencies .!= prev_latencies)
            max_allowed_changes = ceil(Int, 0.01 * n_trials * sum(needs_estimation))

            if changed_trials <= max_allowed_changes
                converged_outer = true
                break
            end

            copyto!(prev_latencies, latencies)
        else
            # If no components need estimation, we only need 1 outer iteration
            converged_outer = true
            break
        end
    end

    # 3. Calculate final residue
    residue = copy(data)
    shifted_comp = zeros(T, n_samples)

    for k = 1:n_components
        comp_k = view(templates, :, :, k)
        @inbounds for t = 1:n_trials
            shift_val = latencies[t, k]
            for c = 1:n_channels
                shift_signal!(shifted_comp, view(comp_k, :, c), shift_val)
                for s = 1:n_samples
                    residue[s, c, t] -= shifted_comp[s]
                end
            end
        end
    end

    return RideResult([c.name for c in comps], templates, latencies, residue, converged_outer, actual_outer_iters)
end
