"""
    RideComponent

Specification for a single RIDE component (`:S`, `:C`, `:R`, etc.).

# Fields
- `name::Symbol`: Component name (e.g. `:S`, `:C`, `:R`).
- `latencies::Union{Vector{Int}, Nothing}`: Known trial sample latencies (centered around their median by
  [`ride_decompose`](@ref)), or `nothing` when latencies are estimated via [`woody_filter`](@ref).
- `window_radius::Int`: Woody-filter search radius in samples when `latencies === nothing`.
- `time_window::Union{Nothing, Tuple{Float64, Float64}}`: Optional active component window `(t_start, t_end)`.
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

RideComponent(
    name::Symbol,
    latencies::Union{AbstractVector{<:Integer},Nothing},
    window_radius::Int;
    time_window::Union{Nothing,Tuple{<:Real,<:Real}} = nothing,
) = RideComponent(name, isnothing(latencies) ? nothing : collect(Int, latencies), window_radius, _normalize_ride_time_window(time_window))

RideComponent(name::Symbol, latencies::AbstractVector{<:Integer}; time_window::Union{Nothing,Tuple{<:Real,<:Real}} = nothing) =
    RideComponent(name, collect(Int, latencies), 0, _normalize_ride_time_window(time_window))

RideComponent(name::Symbol, window_radius::Int; time_window::Union{Nothing,Tuple{<:Real,<:Real}} = nothing) =
    RideComponent(name, nothing, window_radius, _normalize_ride_time_window(time_window))

RideComponent(;
    name::Symbol,
    latencies::Union{AbstractVector{<:Integer},Nothing} = nothing,
    window_radius::Int = 0,
    time_window::Union{Nothing,Tuple{<:Real,<:Real}} = nothing,
) = RideComponent(name, isnothing(latencies) ? nothing : collect(Int, latencies), window_radius, _normalize_ride_time_window(time_window))

"""
    RideResult{T<:Real}

Decomposed RIDE templates (`samples × channels × components`), centered single-trial `latencies`,
`latency_offsets`, residual tensor `residue`, convergence flag, and outer iteration count.
"""
struct RideResult{T<:Real}
    components::Vector{Symbol}
    templates::Array{T,3}
    latencies::Matrix{Int}
    latency_offsets::Vector{Int}
    residue::Array{T,3}
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
) where {T<:Real} = RideResult(components, templates, latencies, zeros(Int, length(components)), residue, converged, iterations)

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
        return s_start:clamp(s_end, s_start, n_samples)
    elseif !isnothing(sample_rate) && sample_rate > 0
        s_start = clamp(round(Int, t_start * Float64(sample_rate)) + 1, 1, n_samples)
        return s_start:clamp(round(Int, t_end*Float64(sample_rate))+1, s_start, n_samples)
    elseif 0.0 <= t_start < t_end <= 1.0 && n_samples > 1
        s_start = clamp(round(Int, t_start * (n_samples - 1)) + 1, 1, n_samples)
        return s_start:clamp(round(Int, t_end*(n_samples-1))+1, s_start, n_samples)
    else
        s_start = clamp(round(Int, t_start), 1, n_samples)
        return s_start:clamp(round(Int, t_end), s_start, n_samples)
    end
end

function _compute_ride_residue!(
    dest::Array{T,3},
    data::AbstractArray{T,3},
    templates::Array{T,3},
    latencies::Matrix{Int};
    skip_comp::Int = 0,
) where {T<:Real}
    _, n_channels, n_trials = size(data)
    n_components = size(templates, 3)

    _foreach_chunk(n_trials) do trial_range
        @inbounds for t in trial_range, c = 1:n_channels
            dest_col = view(dest, :, c, t)
            copyto!(dest_col, view(data, :, c, t))
            for j = 1:n_components
                j == skip_comp && continue
                _subtract_shifted!(dest_col, view(templates, :, c, j), latencies[t, j])
            end
        end
    end
    return dest
end

"""
    ride_decompose(data, comps; outer_iters=4, inner_iters=50, inner_tol=1e-3,
                   method=:median, taper_alpha=0.2, times=nothing, sample_rate=nothing,
                   baseline_range=nothing, rng=Random.default_rng()) -> RideResult

Run Residue Iteration DEcomposition (RIDE) on a 3D tensor `(samples, channels, trials)`.
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

    latencies = zeros(Int, n_trials, n_components)
    latency_offsets = zeros(Int, n_components)
    comp_windows = Vector{UnitRange{Int}}(undef, n_components)
    needs_estimation = falses(n_components)

    for k = 1:n_components
        comp_windows[k] = _resolve_ride_time_window(comps[k].time_window, n_samples, times, sample_rate)
        comp_lats = comps[k].latencies
        if !isnothing(comp_lats)
            length(comp_lats) == n_trials || @minimal_error(
                "Component $(comps[k].name) latencies length ($(length(comp_lats))) does not match number of trials ($n_trials)",
            )
            med_offset = n_trials > 0 ? round(Int, median(comp_lats)) : 0
            latency_offsets[k] = med_offset
            latencies[:, k] .= comp_lats .- med_offset
        else
            needs_estimation[k] = true
        end
    end

    has_estimated = any(needs_estimation)
    if has_estimated && n_trials > 1 && n_samples > 2
        for k = 1:n_components
            if needs_estimation[k]
                latencies[:, k] .= woody_filter(view(data, comp_windows[k], :, :), comps[k].window_radius; template = nothing, rng = rng)
            end
        end
    end

    templates = zeros(T, n_samples, n_channels, n_components)
    residue = Array{T,3}(undef, n_samples, n_channels, n_trials)
    ws =
        (n_samples > 0 && n_channels > 0 && n_trials > 0 && n_components > 0) ?
        _RideIterWorkspace(T, n_samples, n_channels, n_trials, n_components, method) : nothing

    converged_outer = false
    actual_outer_iters = 0
    prev_latencies = copy(latencies)
    max_allowed_changes = ceil(Int, 0.01 * n_trials * count(needs_estimation))

    for outer = 1:outer_iters
        actual_outer_iters = outer

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
            _workspace = ws,
        )

        if has_estimated
            for k = 1:n_components
                needs_estimation[k] || continue
                _compute_ride_residue!(residue, data, templates, latencies; skip_comp = k)
                win_k = comp_windows[k]
                latencies[:, k] .= woody_filter(
                    view(residue, win_k, :, :),
                    comps[k].window_radius;
                    template = view(templates, win_k, :, k),
                    initial_latencies = @view(latencies[:, k]),
                    rng = rng,
                )
            end

            changed_trials = count(i -> @inbounds(latencies[i] != prev_latencies[i]), eachindex(latencies, prev_latencies))
            if changed_trials <= max_allowed_changes
                converged_outer = true
                break
            end
            copyto!(prev_latencies, latencies)
        else
            converged_outer = true
            break
        end
    end

    if has_estimated && !converged_outer && actual_outer_iters > 0
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
            _workspace = ws,
        )
    end

    _compute_ride_residue!(residue, data, templates, latencies; skip_comp = 0)
    return RideResult([c.name for c in comps], templates, latencies, latency_offsets, residue, converged_outer, actual_outer_iters)
end
