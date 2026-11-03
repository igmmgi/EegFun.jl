"""
    RIDE Decomposition Orchestrator

High-level interface to run the Residue Iteration DEcomposition (RIDE) on EEG data.
"""

using Random
using Statistics

"""
    RideComponent

Defines a single RIDE component to be extracted.

# Fields
- `name::Symbol`: Component name (e.g., `:S`, `:R`, `:C`).
- `latencies::Union{Vector{Int}, Nothing}`: Pre-defined sample latencies (relative to epoch zero
  or trial onset) for each trial. For Stimulus-locked (`:S`), this is usually a vector of zeros.
  For Response-locked (`:R`), this is a vector of reaction times (in samples), which `ride_decompose`
  centers around `round(Int, median(rt))` (recording the median offset in `RideResult.latency_offsets`).
  If `nothing`, latencies are initialized and refined via the Woody filter.
- `window_radius::Int`: Search window radius (in samples) for the Woody filter if latencies are estimated.
- `time_window::Union{Nothing,Tuple{Float64,Float64}}`: Optional component time window `(t_start, t_end)`
  in seconds (when `times` or `sample_rate` is provided) or sample indices (for raw tensors without
  time metadata). Used for initial `:C` Woody latency estimation and template tapering.
"""
struct RideComponent
    name::Symbol
    latencies::Union{Vector{Int},Nothing}
    window_radius::Int
    time_window::Union{Nothing,Tuple{Float64,Float64}}

    function RideComponent(
        name::Symbol,
        latencies::Union{Vector{Int},Nothing},
        window_radius::Int,
        time_window::Union{Nothing,Tuple{Float64,Float64}} = nothing,
    )
        window_radius < 0 && @minimal_error("RideComponent window_radius must be >= 0, got $window_radius")
        if !isnothing(time_window) && !(time_window[1] < time_window[2])
            @minimal_error("RideComponent time_window must have start < end, got $time_window")
        end
        return new(name, latencies, window_radius, time_window)
    end
end

_normalize_ride_time_window(::Nothing) = nothing
_normalize_ride_time_window(tw::Tuple{<:Real,<:Real}) = (Float64(tw[1]), Float64(tw[2]))

# Constructor for 3 positional args + optional time_window keyword
RideComponent(
    name::Symbol,
    latencies::Union{AbstractVector{<:Integer},Nothing},
    window_radius::Int;
    time_window::Union{Nothing,Tuple{<:Real,<:Real}} = nothing,
) = RideComponent(name, isnothing(latencies) ? nothing : collect(Int, latencies), window_radius, _normalize_ride_time_window(time_window))

# Default constructor for pre-defined components (e.g., :S or :R)
RideComponent(name::Symbol, latencies::AbstractVector{<:Integer}; time_window::Union{Nothing,Tuple{<:Real,<:Real}} = nothing) =
    RideComponent(name, collect(Int, latencies), 0, _normalize_ride_time_window(time_window))

# Default constructor for estimated components (e.g., :C)
RideComponent(name::Symbol, window_radius::Int; time_window::Union{Nothing,Tuple{<:Real,<:Real}} = nothing) =
    RideComponent(name, nothing, window_radius, _normalize_ride_time_window(time_window))

# Keyword constructor
RideComponent(;
    name::Symbol,
    latencies::Union{AbstractVector{<:Integer},Nothing} = nothing,
    window_radius::Int = 0,
    time_window::Union{Nothing,Tuple{<:Real,<:Real}} = nothing,
) = RideComponent(name, isnothing(latencies) ? nothing : collect(Int, latencies), window_radius, _normalize_ride_time_window(time_window))


"""
    RideResult{T<:Real}

Contains the decomposed templates, centered single-trial latencies, median latency offsets, and residues.
"""
struct RideResult{T<:Real}
    components::Vector{Symbol}
    templates::Array{T,3}        # (samples, channels, components)
    latencies::Matrix{Int}       # (trials, components) — centered relative sample shifts
    latency_offsets::Vector{Int} # median sample offset per component (e.g., round(Int, median(rt)) for :R)
    residue::Array{T,3}          # (samples, channels, trials)
    converged::Bool
    iterations::Int
end

RideResult(
    components::Vector{Symbol},
    templates::Array{T,3},
    latencies::Matrix{Int},
    residue::Array{T,3},
    converged::Bool,
    iterations::Int,
) where {T<:Real} = RideResult(components, templates, latencies, zeros(Int, length(components)), converged, iterations)

"""
Map a component's `time_window` tuple to a valid `UnitRange{Int}` within `1:n_samples`.
"""
function _resolve_ride_time_window(
    time_window::Union{Nothing,Tuple{Float64,Float64}},
    n_samples::Int,
    times::Union{Nothing,AbstractVector{<:Real}},
    sample_rate::Union{Nothing,Real},
)::UnitRange{Int}
    n_samples <= 0 && return 1:0
    isnothing(time_window) && return 1:n_samples

    t_start, t_end = time_window
    if !isnothing(times) && length(times) == n_samples
        s_start = findfirst(>=(t_start), times)
        s_end = findlast(<=(t_end), times)
        isnothing(s_start) && (s_start = argmin(abs.(times .- t_start)))
        isnothing(s_end) && (s_end = argmin(abs.(times .- t_end)))
        s_start = clamp(s_start, 1, n_samples)
        s_end = clamp(s_end, s_start, n_samples)
        return s_start:s_end
    elseif !isnothing(sample_rate) && sample_rate > 0
        s_start = clamp(round(Int, t_start * Float64(sample_rate)) + 1, 1, n_samples)
        s_end = clamp(round(Int, t_end * Float64(sample_rate)) + 1, s_start, n_samples)
        return s_start:s_end
    elseif 0.0 <= t_start < t_end <= 1.0 && n_samples > 1
        s_start = clamp(round(Int, t_start * (n_samples - 1)) + 1, 1, n_samples)
        s_end = clamp(round(Int, t_end * (n_samples - 1)) + 1, s_start, n_samples)
        return s_start:s_end
    else
        s_start = clamp(round(Int, t_start), 1, n_samples)
        s_end = clamp(round(Int, t_end), s_start, n_samples)
        return s_start:s_end
    end
end


"""
    ride_decompose(
        data::AbstractArray{T, 3}, 
        comps::Vector{RideComponent};
        outer_iters::Int = 4,
        inner_iters::Int = 50,
        inner_tol::Real = 1e-3,
        method::Symbol = :median,
        taper_alpha::Real = 0.2,
        times::Union{Nothing,AbstractVector{<:Real}} = nothing,
        sample_rate::Union{Nothing,Real} = nothing,
        baseline_range::Union{Nothing,UnitRange{Int}} = nothing,
        rng::AbstractRNG = Random.default_rng(),
    ) where {T<:Real}

Execute the full RIDE decomposition pipeline on a 3D data tensor `(samples, channels, trials)`.
"""
function ride_decompose(
    data::AbstractArray{T,3},
    comps::Vector{RideComponent};
    outer_iters::Int = 4,
    inner_iters::Int = 50,
    inner_tol::Real = 1e-3,
    method::Symbol = :median,
    taper_alpha::Real = 0.2,
    times::Union{Nothing,AbstractVector{<:Real}} = nothing,
    sample_rate::Union{Nothing,Real} = nothing,
    baseline_range::Union{Nothing,UnitRange{Int}} = nothing,
    rng::AbstractRNG = Random.default_rng(),
) where {T<:Real}

    n_samples, n_channels, n_trials = size(data)
    n_components = length(comps)

    # Initialize centered latencies matrix (trials, components) and median offsets
    latencies = zeros(Int, n_trials, n_components)
    latency_offsets = zeros(Int, n_components)
    comp_windows = Vector{UnitRange{Int}}(undef, n_components)

    # Track which components need estimation
    needs_estimation = falses(n_components)

    for k = 1:n_components
        comp_windows[k] = _resolve_ride_time_window(comps[k].time_window, n_samples, times, sample_rate)

        if !isnothing(comps[k].latencies)
            if length(comps[k].latencies) != n_trials
                @minimal_error(
                    "Component $(comps[k].name) latencies length ($(length(comps[k].latencies))) does not match number of trials ($n_trials)",
                )
            end
            # Center pre-defined latencies (e.g., reaction times for :R) around their median
            med_offset = n_trials > 0 ? round(Int, median(comps[k].latencies)) : 0
            latency_offsets[k] = med_offset
            latencies[:, k] .= comps[k].latencies .- med_offset
        else
            needs_estimation[k] = true
        end
    end

    # Initialize unknown component latencies (e.g. :C) BEFORE the first ride_iteration!
    # using leave-one-out Woody filter within each component's time_window so :C never collapses to 0.0
    if any(needs_estimation) && n_trials > 1 && n_samples > 2
        for k = 1:n_components
            if needs_estimation[k]
                win_k = comp_windows[k]
                init_slice = view(data, win_k, :, :)
                latencies[:, k] .= woody_filter(init_slice, comps[k].window_radius; template = nothing, rng = rng)
            end
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

        # 1. EM Template Iteration (modifies `templates` in-place, with Tukey boundary tapering)
        _, inner_converged = ride_iteration!(
            templates,
            data,
            latencies;
            max_iter = inner_iters,
            tol = inner_tol,
            method = method,
            taper_alpha = taper_alpha,
            component_windows = comp_windows,
            baseline_range = baseline_range,
        )

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

                    # Run Woody filter within component time_window using current template & latency guesses
                    win_k = comp_windows[k]
                    current_template = view(templates, win_k, :, k)
                    initial_lats = latencies[:, k]

                    new_lats = woody_filter(
                        view(woody_data, win_k, :, :),
                        comps[k].window_radius;
                        template = current_template,
                        initial_latencies = initial_lats,
                        rng = rng,
                    )

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

    # Final pass of template extraction with converged latencies if latencies changed on the final outer step
    if any(needs_estimation) && !converged_outer && actual_outer_iters > 0
        ride_iteration!(
            templates,
            data,
            latencies;
            max_iter = inner_iters,
            tol = inner_tol,
            method = method,
            taper_alpha = taper_alpha,
            component_windows = comp_windows,
            baseline_range = baseline_range,
        )
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

    return RideResult([c.name for c in comps], templates, latencies, latency_offsets, residue, converged_outer, actual_outer_iters)
end
