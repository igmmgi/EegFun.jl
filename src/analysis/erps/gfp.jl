function _normalize_to_percent!(vals::Vector{Float64})
    isempty(vals) && return vals
    v_min, v_max = extrema(vals)
    span = v_max - v_min
    if span > eps()
        scale = 100.0 / span
        @inbounds @simd for i in eachindex(vals)
            vals[i] = clamp((vals[i] - v_min) * scale, 0.0, 100.0)
        end
    end
    return vals
end

function _gfp_and_gd_cols(cols::Vector{Vector{Float64}}, n_tp::Int, need_gd::Bool)
    n_ch = length(cols)
    inv_ch = 1.0 / n_ch

    means = zeros(Float64, n_tp)
    for col in cols
        @inbounds @simd for t = 1:n_tp
            means[t] += col[t]
        end
    end
    @inbounds @simd for t = 1:n_tp
        means[t] *= inv_ch
    end

    gfp_vals = zeros(Float64, n_tp)
    for col in cols
        @inbounds @simd for t = 1:n_tp
            d = col[t] - means[t]
            gfp_vals[t] += d * d
        end
    end
    @inbounds @simd for t = 1:n_tp
        gfp_vals[t] = sqrt(gfp_vals[t] * inv_ch)
    end

    !need_gd && return gfp_vals, Float64[]

    gd_vals = zeros(Float64, n_tp)
    n_tp <= 1 && return gfp_vals, gd_vals

    inv_gfp = Vector{Float64}(undef, n_tp)
    @inbounds for t = 1:n_tp
        g = gfp_vals[t]
        inv_gfp[t] = g > eps(Float64) ? 1.0 / g : 0.0
    end

    for col in cols
        @inbounds @simd for t = 2:n_tp
            d = (col[t] - means[t]) * inv_gfp[t] - (col[t-1] - means[t-1]) * inv_gfp[t-1]
            gd_vals[t] += d * d
        end
    end
    @inbounds for t = 2:n_tp
        gd_vals[t] = (inv_gfp[t] == 0.0 || inv_gfp[t-1] == 0.0) ? 0.0 : sqrt(gd_vals[t] * inv_ch)
    end
    gd_vals[1] = gd_vals[2]
    return gfp_vals, gd_vals
end

function _compute_gfp_metrics(dat::ErpData, channel_selection, operation::String; need_gd::Bool)
    selected_channels = get_selected_channels(dat, channel_selection, include_meta = false, include_extra = false)
    isempty(selected_channels) && @minimal_error("No channels selected for $operation")
    @info "Using $(length(selected_channels)) channels for $operation"

    cols = Vector{Float64}[
        let c = dat.data[!, ch]
            c isa Vector{Float64} ? c : Float64.(c)
        end for ch in selected_channels
    ]
    return _gfp_and_gd_cols(cols, nrow(dat.data), need_gd)
end

function _build_gfp_result_df(dat::ErpData, result_pairs::Pair{Symbol,Vector{Float64}}...)
    result_df = DataFrame(; copycols = false)
    for col in meta_labels(dat)
        result_df[!, col] = copy(dat.data[!, col])
    end
    n = nrow(dat.data)
    hasproperty(result_df, :condition) || (result_df[!, :condition] = fill(dat.condition, n))
    hasproperty(result_df, :condition_name) || (result_df[!, :condition_name] = fill(dat.condition_name, n))
    for (name, values) in result_pairs
        result_df[!, name] = values
    end
    return result_df
end

function _apply_to_conditions(
    fn::Function,
    dat::AbstractVector{<:ErpData};
    condition_selection = conditions(),
    operation::String = "calculation",
    kwargs...,
)::Vector{DataFrame}
    dat_filtered = dat[get_selected_conditions(dat, condition_selection)]
    @info "$operation for $(length(dat_filtered)) dataset(s)"
    return [fn(erp_data; kwargs...) for erp_data in dat_filtered]
end

"""
    gfp(dat::ErpData; channel_selection = channels(), normalize::Bool = false) -> DataFrame
    gfp(dat::AbstractVector{<:ErpData}; condition_selection = conditions(), channel_selection = channels(), normalize::Bool = false) -> Vector{DataFrame}

Calculate Global Field Power (GFP) as the population standard deviation across selected channels
at each time point (Lehmann & Skrandies, 1980).

# Arguments
- `dat`: `ErpData` or `Vector{ErpData}`
- `channel_selection`: Channel predicate or symbols (default: `channels()`)
- `normalize::Bool`: Normalize GFP to `0–100%` range (default: `false`)

# Returns
- `DataFrame` (or `Vector{DataFrame}`) with `:time`, metadata columns, and `:gfp`.
"""
function gfp(dat::ErpData; channel_selection = channels(), normalize::Bool = false)::DataFrame
    @info "Calculating Global Field Power (GFP)"
    gfp_vals, _ = _compute_gfp_metrics(dat, channel_selection, "GFP calculation"; need_gd = false)
    normalize && _normalize_to_percent!(gfp_vals)
    @info "GFP calculation complete"
    return _build_gfp_result_df(dat, :gfp => gfp_vals)
end

function gfp(
    dat::AbstractVector{<:ErpData};
    condition_selection = conditions(),
    channel_selection = channels(),
    normalize::Bool = false,
)::Vector{DataFrame}
    return _apply_to_conditions(gfp, dat; condition_selection, operation = "Calculating GFP", channel_selection, normalize)
end

"""
    global_dissimilarity(dat::ErpData; channel_selection = channels(), normalize::Bool = false) -> DataFrame
    global_dissimilarity(dat::AbstractVector{<:ErpData}; condition_selection = conditions(), channel_selection = channels(), normalize::Bool = false) -> Vector{DataFrame}

Calculate reference-independent Global Dissimilarity (GD) across adjacent time points
after subtracting the spatial channel mean and normalizing by GFP (Lehmann & Skrandies, 1980).

# Arguments
- `dat`: `ErpData` or `Vector{ErpData}`
- `channel_selection`: Channel predicate or symbols (default: `channels()`)
- `normalize::Bool`: Normalize GD to `0–100%` range (default: `false`)

# Returns
- `DataFrame` (or `Vector{DataFrame}`) with `:time`, metadata columns, and `:dissimilarity`.
"""
function global_dissimilarity(dat::ErpData; channel_selection = channels(), normalize::Bool = false)::DataFrame
    @info "Calculating Global Dissimilarity (GD)"
    _, gd_vals = _compute_gfp_metrics(dat, channel_selection, "dissimilarity calculation"; need_gd = true)
    normalize && _normalize_to_percent!(gd_vals)
    @info "Global Dissimilarity calculation complete"
    return _build_gfp_result_df(dat, :dissimilarity => gd_vals)
end

function global_dissimilarity(
    dat::AbstractVector{<:ErpData};
    condition_selection = conditions(),
    channel_selection = channels(),
    normalize::Bool = false,
)::Vector{DataFrame}
    return _apply_to_conditions(
        global_dissimilarity,
        dat;
        condition_selection,
        operation = "Calculating Global Dissimilarity",
        channel_selection,
        normalize,
    )
end

"""
    gfp_and_dissimilarity(dat::ErpData; channel_selection = channels(), normalize::Bool = false) -> DataFrame
    gfp_and_dissimilarity(dat::AbstractVector{<:ErpData}; condition_selection = conditions(), channel_selection = channels(), normalize::Bool = false) -> Vector{DataFrame}

Calculate both Global Field Power (`:gfp`) and Global Dissimilarity (`:dissimilarity`) in a single pass.
"""
function gfp_and_dissimilarity(dat::ErpData; channel_selection = channels(), normalize::Bool = false)::DataFrame
    @info "Calculating Global Field Power and Global Dissimilarity"
    gfp_vals, gd_vals = _compute_gfp_metrics(dat, channel_selection, "GFP/dissimilarity calculation"; need_gd = true)
    if normalize
        _normalize_to_percent!(gfp_vals)
        _normalize_to_percent!(gd_vals)
    end
    @info "GFP and Global Dissimilarity calculation complete"
    return _build_gfp_result_df(dat, :gfp => gfp_vals, :dissimilarity => gd_vals)
end

function gfp_and_dissimilarity(
    dat::AbstractVector{<:ErpData};
    condition_selection = conditions(),
    channel_selection = channels(),
    normalize::Bool = false,
)::Vector{DataFrame}
    return _apply_to_conditions(
        gfp_and_dissimilarity,
        dat;
        condition_selection,
        operation = "Calculating GFP and Global Dissimilarity",
        channel_selection,
        normalize,
    )
end
