# === Artifact Identification Helpers ===

"""
    _robust_z(x::AbstractVector{<:Real})

Compute robust z-scores using median and Median Absolute Deviation (MAD),
scaled by 1.4826 to match standard deviation under normality.
Falls back to standard z-score if MAD == 0, or zeros if variance is zero.
"""
function _robust_z(x::AbstractVector{<:Real})
    n = length(x)
    n == 0 && return Float64[]
    n == 1 && return [0.0]
    med = median(x)
    mad = median(abs.(x .- med))
    s = 1.482602218505602 * mad
    if s > eps(Float64)
        return (x .- med) ./ s
    else
        s_std = std(x)
        if s_std > eps(Float64)
            return (x .- mean(x)) ./ s_std
        else
            return zeros(Float64, n)
        end
    end
end

"""
    _zscore_for_remaining(corr_vector::AbstractVector{<:Real}, remaining_indices::AbstractVector{Int})

Z-score a correlation vector using robust statistics on only the values at `remaining_indices`.
"""
function _zscore_for_remaining(corr_vector::AbstractVector{<:Real}, remaining_indices::AbstractVector{Int})
    z = fill(NaN, length(corr_vector))
    if length(remaining_indices) >= 2
        remaining_vals = corr_vector[remaining_indices]
        z[remaining_indices] = _robust_z(remaining_vals)
    elseif length(remaining_indices) == 1
        z[remaining_indices] .= 0.0
    end
    return z
end

"""
    _calculate_eog_correlations(components_t::Matrix{Float64}, eog_signal::AbstractVector{<:Real})

Compute absolute correlation of each ICA component with the given EOG signal.
`components_t` has dimensions `(n_samples, n_components)`.
"""
function _calculate_eog_correlations(components_t::Matrix{Float64}, eog_signal::AbstractVector{<:Real})
    n_samples, n_components = size(components_t)
    corrs = zeros(Float64, n_components)
    for comp_idx = 1:n_components
        comp_signal = view(components_t, :, comp_idx)
        corrs[comp_idx] = abs(cor(comp_signal, eog_signal))
    end
    return corrs
end

"""
    _calculate_spatial_correlations(remaining_components::Vector{Int}, primary_components::Vector{Int}, mixing_matrix::Matrix{Float64}, n_components::Int)

Compute maximum spatial (topography) correlation between each remaining component and primary EOG components.
"""
function _calculate_spatial_correlations(
    remaining_components::Vector{Int},
    primary_components::Vector{Int},
    mixing_matrix::Matrix{Float64},
    n_components::Int,
)
    spatial_corrs = fill(NaN, n_components)
    if isempty(remaining_components) || isempty(primary_components)
        return spatial_corrs
    end

    for comp_idx in remaining_components
        comp_topography = view(mixing_matrix, :, comp_idx)
        max_corr = 0.0
        for primary_comp in primary_components
            primary_topography = view(mixing_matrix, :, primary_comp)
            c = abs(cor(comp_topography, primary_topography))
            if c > max_corr
                max_corr = c
            end
        end
        spatial_corrs[comp_idx] = max_corr
    end

    return spatial_corrs
end

"""
    _max_lagged_corr(x_norm::AbstractVector{Float64}, y_norm::AbstractVector{Float64}, lags)

Compute maximum absolute cross-correlation using BLAS dot products over pre-normalized signals.
"""
function _max_lagged_corr(x_norm::AbstractVector{Float64}, y_norm::AbstractVector{Float64}, lags)
    n = length(x_norm)
    best = 0.0
    for lag in lags
        a = max(1, 1 - lag)
        b = min(n, n - lag)
        a > b && continue
        c = abs(dot(view(x_norm, a:b), view(y_norm, (a+lag):(b+lag))))
        if c > best
            best = c
        end
    end
    return best
end

"""
    _calculate_lagged_correlations(remaining_components::Vector{Int}, primary_components::Vector{Int}, components_t::Matrix{Float64}, max_lag_samples::Int, lag_step::Int)

Compute maximum absolute lagged cross-correlation between remaining and primary components via pre-normalized BLAS dot products.
"""
function _calculate_lagged_correlations(
    remaining_components::Vector{Int},
    primary_components::Vector{Int},
    components_t::Matrix{Float64},
    max_lag_samples::Int,
    lag_step::Int,
)
    n_samples, n_components = size(components_t)
    lagged_corrs = fill(NaN, n_components)
    if isempty(remaining_components) || isempty(primary_components)
        return lagged_corrs
    end

    lags = (-max_lag_samples):lag_step:max_lag_samples

    # Pre-normalize primary components once
    primary_norms = Vector{Vector{Float64}}(undef, length(primary_components))
    for (i, p_idx) in enumerate(primary_components)
        sig = view(components_t, :, p_idx)
        denom = std(sig, corrected = false) * sqrt(n_samples)
        if denom > eps(Float64)
            primary_norms[i] = (sig .- mean(sig)) ./ denom
        else
            primary_norms[i] = zeros(Float64, n_samples)
        end
    end

    # Pre-normalize remaining components and compute max lagged correlation
    Threads.@threads for comp_idx in remaining_components
        sig = view(components_t, :, comp_idx)
        denom = std(sig, corrected = false) * sqrt(n_samples)
        if denom > eps(Float64)
            x_norm = (sig .- mean(sig)) ./ denom
            max_c = 0.0
            for y_norm in primary_norms
                c = _max_lagged_corr(x_norm, y_norm, lags)
                if c > max_c
                    max_c = c
                end
            end
            lagged_corrs[comp_idx] = max_c
        else
            lagged_corrs[comp_idx] = 0.0
        end
    end

    return lagged_corrs
end

"""
    _identify_secondary_components(remaining_components::Vector{Int}, corr_vector::Vector{Float64}, corr_z_vector::Vector{Float64}, min_corr_threshold::Real, z_threshold::Real)

Return indices of remaining components exceeding both minimum correlation and positive z-score thresholds.
"""
function _identify_secondary_components(
    remaining_components::Vector{Int},
    corr_vector::Vector{Float64},
    corr_z_vector::Vector{Float64},
    min_corr_threshold::Real,
    z_threshold::Real,
)
    isempty(remaining_components) && return Int[]
    corr_mask = corr_vector[remaining_components] .> min_corr_threshold
    # Only positive z-score is meaningful for non-negative correlations
    z_mask = corr_z_vector[remaining_components] .> z_threshold
    combined_mask = corr_mask .& z_mask
    return remaining_components[combined_mask]
end

# === Typed Empty DataFrame Constructors ===

function _empty_eog_metrics()
    return DataFrame(
        :Component => Int[],
        :vEOG_corr => Float64[],
        :vEOG_zscore => Float64[],
        :vEOG_spatial_corr => Float64[],
        :vEOG_spatial_corr_z => Float64[],
        :vEOG_temporal_corr => Float64[],
        :vEOG_temporal_corr_z => Float64[],
        :hEOG_corr => Float64[],
        :hEOG_zscore => Float64[],
        :hEOG_spatial_corr => Float64[],
        :hEOG_spatial_corr_z => Float64[],
        :hEOG_temporal_corr => Float64[],
        :hEOG_temporal_corr_z => Float64[],
    )
end

function _empty_ecg_metrics()
    return DataFrame(
        :Component => Int[],
        :num_peaks => Int[],
        :num_valid_ibis => Int[],
        :mean_ibi_s => Float64[],
        :std_ibi_s => Float64[],
        :peak_ratio => Float64[],
        :heart_rate_bpm => Float64[],
    )
end

function _empty_line_noise_metrics()
    return DataFrame(
        :Component => Int[],
        :line_power => Float64[],
        :flanking_power => Float64[],
        :power_ratio => Float64[],
        :harmonic_ratio => Float64[],
        :power_ratio_zscore => Float64[],
    )
end

# === EOG Identification ===

"""
    identify_eog_components(dat::ContinuousData, ica::InfoIca;
                          vEOG_channel::Symbol=:vEOG,
                          hEOG_channel::Symbol=:hEOG,
                          z_threshold::Real=3.0,
                          min_correlation::Real=0.5,
                          two_step::Bool=true,
                          sample_selection::Function = samples(),
                          interval_selection::Interval = times())

Identify ICA components potentially related to eye movements based on robust z-scored correlation with EOG channels.

Uses a two-step approach by default:
1. **Step 1**: Calculate correlations between all ICA components and EOG channels, compute robust z-scores (MAD-scaled), 
   and identify primary EOG components (those exceeding both the z-score and correlation thresholds).
2. **Step 2**: If `two_step=true`, calculate spatial correlations (topographies) and lagged temporal cross-correlations 
   between remaining components and primary EOG components. Identify secondary EOG components based on both spatial 
   or temporal correlation and z-score thresholds.

# Arguments
- `dat::ContinuousData`: The continuous data containing EOG channels.
- `ica::InfoIca`: The ICA result object.

# Keyword Arguments
- `vEOG_channel::Symbol`: Name of the vertical EOG channel (default: `:vEOG`).
- `hEOG_channel::Symbol`: Name of the horizontal EOG channel (default: `:hEOG`).
- `z_threshold::Real`: Robust Z-score threshold for step 1 identification (default: `3.0`). Primary components must exceed 
  this z-score threshold to be identified in step 1.
- `min_correlation::Real`: Minimum absolute correlation threshold with EOG channels (default: `0.5`).
- `two_step::Bool`: If true (default), perform step 2 to identify secondary components sharing spatial topographies or 
  lagged temporal profiles with primary components. If false, only step 1 primary components are returned.
- `sample_selection::Function`: Function to select samples from `dat.data` (default: `samples()`).
- `interval_selection::Interval`: Time interval filter (default: `times()`).

# Returns
- `Dict{Symbol, Vector{Int}}`: Dictionary containing:
  - `:vEOG`: Vector of component indices identified for vertical eye movements.
  - `:hEOG`: Vector of component indices identified for horizontal eye movements.
- `DataFrame`: DataFrame containing detailed metrics per component.
"""
function identify_eog_components(
    dat::ContinuousData,
    ica::InfoIca;
    vEOG_channel::Symbol = :vEOG,
    hEOG_channel::Symbol = :hEOG,
    sample_selection = samples(),
    interval_selection::Interval = times(),
    z_threshold::Real = 3.0,
    min_correlation::Real = 0.5,
    two_step::Bool = true,
    _precomputed_components::Union{Nothing,Tuple{Matrix{Float64},Int}} = nothing,
    _precomputed_samples::Union{Nothing,Vector{Int}} = nothing,
)
    n_components = size(ica.unmixing, 1)

    has_vEOG = vEOG_channel in propertynames(dat.data)
    has_hEOG = hEOG_channel in propertynames(dat.data)

    if !has_vEOG && !has_hEOG
        @minimal_warning "Neither $vEOG_channel nor $hEOG_channel found in data. Skipping EOG component identification."
        return Dict(:vEOG => Int[], :hEOG => Int[]), _empty_eog_metrics()
    end
    if !has_vEOG
        @minimal_warning "Vertical EOG channel :$vEOG_channel not found in data; skipping vertical EOG detection."
    end
    if !has_hEOG
        @minimal_warning "Horizontal EOG channel :$hEOG_channel not found in data; skipping horizontal EOG detection."
    end

    # Get samples to use (skip if pre-computed)
    selected_samples = if !isnothing(_precomputed_samples)
        _precomputed_samples
    else
        combined_sel = _combine_interval_sample(interval_selection, sample_selection)
        get_selected_samples(dat, combined_sel)
    end
    if isempty(selected_samples)
        @minimal_warning "No samples remaining after applying exclude criteria. Cannot identify eye components."
        return Dict(:vEOG => Int[], :hEOG => Int[]), _empty_eog_metrics()
    end

    if n_components <= 10
        @info "Small number of components ($n_components <= 10); z_threshold = $z_threshold may be conservative."
    end

    # Calculate components for valid samples (skip if pre-computed)
    components, n_comps_calc = if !isnothing(_precomputed_components)
        _precomputed_components
    else
        _prepare_ica_data_matrix(dat, ica, selected_samples)
    end

    # Store components as (n_samples, n_components) for contiguous column access
    components_t = Matrix{Float64}(components')

    # Step 1: EOG correlations and robust z-scores
    vEOG_corrs = fill(NaN, n_components)
    vEOG_corr_z = fill(NaN, n_components)
    primary_vEOG = Int[]
    if has_vEOG
        vEOG = dat.data[selected_samples, vEOG_channel]
        vEOG_corrs = _calculate_eog_correlations(components_t, vEOG)
        vEOG_corr_z = _robust_z(vEOG_corrs)
        primary_vEOG = findall((vEOG_corr_z .> z_threshold) .& (vEOG_corrs .> min_correlation))
    end

    hEOG_corrs = fill(NaN, n_components)
    hEOG_corr_z = fill(NaN, n_components)
    primary_hEOG = Int[]
    if has_hEOG
        hEOG = dat.data[selected_samples, hEOG_channel]
        hEOG_corrs = _calculate_eog_correlations(components_t, hEOG)
        hEOG_corr_z = _robust_z(hEOG_corrs)
        primary_hEOG = findall((hEOG_corr_z .> z_threshold) .& (hEOG_corrs .> min_correlation))
    end

    vEOG_spatial_corr = fill(NaN, n_components)
    vEOG_spatial_corr_z = fill(NaN, n_components)
    vEOG_temporal_corr = fill(NaN, n_components)
    vEOG_temporal_corr_z = fill(NaN, n_components)

    hEOG_spatial_corr = fill(NaN, n_components)
    hEOG_spatial_corr_z = fill(NaN, n_components)
    hEOG_temporal_corr = fill(NaN, n_components)
    hEOG_temporal_corr_z = fill(NaN, n_components)

    secondary_vEOG = Int[]
    secondary_hEOG = Int[]

    # Step 2: Compute spatial and lagged temporal correlations ONLY when two_step is true
    if two_step
        max_lag_samples = round(Int, 100.0 * dat.sample_rate / 1000.0) # 100 ms lag
        lag_step = max(1, round(Int, 5.0 * dat.sample_rate / 1000.0)) # 5 ms step

        if has_vEOG && !isempty(primary_vEOG)
            remaining_vEOG = setdiff(1:n_components, primary_vEOG)
            vEOG_spatial_corr = _calculate_spatial_correlations(remaining_vEOG, primary_vEOG, ica.mixing, n_components)
            vEOG_temporal_corr = _calculate_lagged_correlations(remaining_vEOG, primary_vEOG, components_t, max_lag_samples, lag_step)
            vEOG_spatial_corr_z = _zscore_for_remaining(vEOG_spatial_corr, remaining_vEOG)
            vEOG_temporal_corr_z = _zscore_for_remaining(vEOG_temporal_corr, remaining_vEOG)

            sec_v_spat =
                _identify_secondary_components(remaining_vEOG, vEOG_spatial_corr, vEOG_spatial_corr_z, min_correlation, z_threshold)
            sec_v_temp =
                _identify_secondary_components(remaining_vEOG, vEOG_temporal_corr, vEOG_temporal_corr_z, min_correlation, z_threshold)
            secondary_vEOG = union(sec_v_spat, sec_v_temp)
        end

        if has_hEOG && !isempty(primary_hEOG)
            remaining_hEOG = setdiff(1:n_components, primary_hEOG)
            hEOG_spatial_corr = _calculate_spatial_correlations(remaining_hEOG, primary_hEOG, ica.mixing, n_components)
            hEOG_temporal_corr = _calculate_lagged_correlations(remaining_hEOG, primary_hEOG, components_t, max_lag_samples, lag_step)
            hEOG_spatial_corr_z = _zscore_for_remaining(hEOG_spatial_corr, remaining_hEOG)
            hEOG_temporal_corr_z = _zscore_for_remaining(hEOG_temporal_corr, remaining_hEOG)

            sec_h_spat =
                _identify_secondary_components(remaining_hEOG, hEOG_spatial_corr, hEOG_spatial_corr_z, min_correlation, z_threshold)
            sec_h_temp =
                _identify_secondary_components(remaining_hEOG, hEOG_temporal_corr, hEOG_temporal_corr_z, min_correlation, z_threshold)
            secondary_hEOG = union(sec_h_spat, sec_h_temp)
        end
    end

    identified_vEOG = two_step ? union(primary_vEOG, secondary_vEOG) : primary_vEOG
    identified_hEOG = two_step ? union(primary_hEOG, secondary_hEOG) : primary_hEOG

    sort!(identified_vEOG)
    sort!(identified_hEOG)

    metrics_df = DataFrame(
        :Component => 1:n_components,
        :vEOG_corr => vEOG_corrs,
        :vEOG_zscore => vEOG_corr_z,
        :vEOG_spatial_corr => vEOG_spatial_corr,
        :vEOG_spatial_corr_z => vEOG_spatial_corr_z,
        :vEOG_temporal_corr => vEOG_temporal_corr,
        :vEOG_temporal_corr_z => vEOG_temporal_corr_z,
        :hEOG_corr => hEOG_corrs,
        :hEOG_zscore => hEOG_corr_z,
        :hEOG_spatial_corr => hEOG_spatial_corr,
        :hEOG_spatial_corr_z => hEOG_spatial_corr_z,
        :hEOG_temporal_corr => hEOG_temporal_corr,
        :hEOG_temporal_corr_z => hEOG_temporal_corr_z,
    )

    result_dict = Dict{Symbol,Vector{Int}}(:vEOG => identified_vEOG, :hEOG => identified_hEOG)

    return result_dict, metrics_df
end

# === ECG Identification ===

"""
    _find_unidirectional_peaks(data::AbstractVector; min_prominence_std::Real = 3.0, window_size::Int = 1)

Find local peaks in a single direction exceeding `mean + min_prominence_std * std`.
"""
function _find_unidirectional_peaks(data::AbstractVector; min_prominence_std::Real = 3.0, window_size::Int = 1)
    n = length(data)
    n < 2 * window_size + 1 && return Int[]
    peaks = Int[]
    m = mean(data)
    s = std(data)
    thresh = (s ≈ 0) ? m : m + min_prominence_std * s
    @inbounds for i = (window_size+1):(n-window_size)
        val = data[i]
        if val > thresh
            is_p = true
            for j = (i-window_size):(i-1)
                if val <= data[j]
                    is_p = false
                    break
                end
            end
            if is_p
                for j = (i+1):(i+window_size)
                    if val <= data[j]
                        is_p = false
                        break
                    end
                end
            end
            is_p && push!(peaks, i)
        end
    end
    return peaks
end

"""
    _calculate_ibi_metrics(peak_indices::Vector{Int}, sample_rate::Real, min_ibi_s::Real, max_ibi_s::Real)

Calculate inter-beat interval (IBI) metrics from peak indices.
"""
function _calculate_ibi_metrics(peak_indices::Vector{Int}, sample_rate::Real, min_ibi_s::Real, max_ibi_s::Real)
    num_peaks = length(peak_indices)
    num_peaks < 2 && return 0, NaN, NaN, 0.0

    ibis_s = diff(peak_indices) ./ sample_rate

    valid_ibi_mask = (ibis_s .>= min_ibi_s) .& (ibis_s .<= max_ibi_s)
    valid_ibis = ibis_s[valid_ibi_mask]
    num_valid_ibis = length(valid_ibis)

    peak_ratio = num_valid_ibis / (num_peaks - 1)

    if num_valid_ibis > 1
        mean_ibi = mean(valid_ibis)
        std_ibi = std(valid_ibis)
    elseif num_valid_ibis == 1
        mean_ibi = valid_ibis[1]
        std_ibi = 0.0
    else
        mean_ibi = NaN
        std_ibi = NaN
    end

    return num_valid_ibis, mean_ibi, std_ibi, peak_ratio
end

"""
    identify_ecg_components(dat::ContinuousData, ica::InfoIca;
                          min_bpm::Real=40, max_bpm::Real=120,
                          min_prominence_std::Real=4,
                          min_peaks::Int=100,
                          max_ibi_std_s::Real=0.2,
                          min_peak_ratio::Real=0.5,
                          sample_selection::Function = samples(),
                          interval_selection::Interval = times())

Identify ICA components potentially related to ECG artifacts based on peak detection
and interval regularity, evaluating both positive and negative R-wave polarities.
"""
function identify_ecg_components(
    dat::ContinuousData,
    ica::InfoIca;
    min_bpm::Real = 40,
    max_bpm::Real = 120,
    min_prominence_std::Real = 4,
    min_peaks::Int = 100,
    max_ibi_std_s::Real = 0.2,
    min_peak_ratio::Real = 0.5,
    sample_selection = samples(),
    interval_selection::Interval = times(),
    _precomputed_components::Union{Nothing,Tuple{Matrix{Float64},Int}} = nothing,
    _precomputed_samples::Union{Nothing,Vector{Int}} = nothing,
)
    selected_samples = if !isnothing(_precomputed_samples)
        _precomputed_samples
    else
        combined_sel = _combine_interval_sample(interval_selection, sample_selection)
        get_selected_samples(dat, combined_sel)
    end
    if isempty(selected_samples)
        @minimal_warning "No samples remaining after applying exclude criteria."
        return Int[], _empty_ecg_metrics()
    end

    components_subset, n_components = if !isnothing(_precomputed_components)
        _precomputed_components
    else
        _prepare_ica_data_matrix(dat, ica, selected_samples)
    end

    if n_components == 0
        return Int[], _empty_ecg_metrics()
    end

    min_ibi_s = 60.0 / max_bpm
    max_ibi_s = 60.0 / min_bpm

    # Adapt min_peaks for short recordings to prevent rejecting true ECG components
    duration_s = (maximum(selected_samples) - minimum(selected_samples) + 1) / Float64(dat.sample_rate)
    expected_peaks = (duration_s * min_bpm) / 60.0
    effective_min_peaks = min(min_peaks, max(3, round(Int, expected_peaks * 0.5)))
    if effective_min_peaks < min_peaks
        @minimal_warning "Selected recording duration ($(round(duration_s, digits=1)) s) is short for default min_peaks = $min_peaks; using adapted threshold min_peaks = $effective_min_peaks."
    end

    metrics = NamedTuple{
        (:Component, :num_peaks, :num_valid_ibis, :mean_ibi_s, :std_ibi_s, :peak_ratio, :heart_rate_bpm),
        Tuple{Int,Int,Int,Float64,Float64,Float64,Float64},
    }[]
    identified_ecg = Int[]

    n_samples = size(components_subset, 2)
    ts_zscored = Vector{Float64}(undef, n_samples)

    for comp_idx = 1:n_components
        component_ts = @view components_subset[comp_idx, :]

        m = mean(component_ts)
        s = std(component_ts)
        ts_zscored .= (component_ts .- m) ./ (s > eps(Float64) ? s : 1.0)

        # Detect peaks in both polarities to handle arbitrary ICA signs
        pk_pos = _find_unidirectional_peaks(ts_zscored; min_prominence_std = min_prominence_std)
        pk_neg = _find_unidirectional_peaks(-ts_zscored; min_prominence_std = min_prominence_std)

        # Map back to original sample coordinates to prevent sewing across non-contiguous gaps
        pos_orig = selected_samples[pk_pos]
        neg_orig = selected_samples[pk_neg]

        n_valid_pos, mean_pos, std_pos, ratio_pos = _calculate_ibi_metrics(pos_orig, dat.sample_rate, min_ibi_s, max_ibi_s)
        n_valid_neg, mean_neg, std_neg, ratio_neg = _calculate_ibi_metrics(neg_orig, dat.sample_rate, min_ibi_s, max_ibi_s)

        # Choose the polarity yielding more valid IBIs
        if n_valid_pos >= n_valid_neg
            num_peaks = length(pk_pos)
            num_valid_ibis, mean_ibi, std_ibi, peak_ratio = n_valid_pos, mean_pos, std_pos, ratio_pos
        else
            num_peaks = length(pk_neg)
            num_valid_ibis, mean_ibi, std_ibi, peak_ratio = n_valid_neg, mean_neg, std_neg, ratio_neg
        end

        has_sufficient_ibis = num_valid_ibis >= (effective_min_peaks - 1)
        has_low_ibi_variability = !isnan(std_ibi) && std_ibi <= max_ibi_std_s
        has_good_peak_ratio = peak_ratio >= min_peak_ratio

        is_ecg = has_sufficient_ibis && has_low_ibi_variability && has_good_peak_ratio
        is_ecg && push!(identified_ecg, comp_idx)

        heart_rate_bpm = (!isnan(mean_ibi) && mean_ibi > 0) ? (60.0 / mean_ibi) : NaN

        push!(
            metrics,
            (
                Component = comp_idx,
                num_peaks = num_peaks,
                num_valid_ibis = num_valid_ibis,
                mean_ibi_s = mean_ibi,
                std_ibi_s = std_ibi,
                peak_ratio = peak_ratio,
                heart_rate_bpm = heart_rate_bpm,
            ),
        )
    end

    metrics_df = DataFrame(metrics)
    sort!(identified_ecg)

    return identified_ecg, metrics_df
end

# === Spatial Kurtosis (Channel Noise) Identification ===

"""
    identify_spatial_kurtosis_components(ica::InfoIca; z_threshold::Real = 3.0)

Identify ICA components with high spatial kurtosis (localized, spot-like channel noise activity).
Uses robust MAD-based z-scoring across components.
"""
function identify_spatial_kurtosis_components(ica::InfoIca; z_threshold::Real = 3.0)
    n_components = size(ica.mixing, 2)
    spatial_kurtosis = Float64[]
    for i = 1:n_components
        weights = view(ica.mixing, :, i)
        k = kurtosis(weights)
        push!(spatial_kurtosis, k)
    end

    z_spatial_kurtosis = _robust_z(spatial_kurtosis)
    high_kurtosis_comps = findall(z_spatial_kurtosis .> z_threshold)
    sort!(high_kurtosis_comps)

    metrics_df = DataFrame(:Component => 1:n_components, :spatial_kurtosis => spatial_kurtosis, :z_spatial_kurtosis => z_spatial_kurtosis)

    return high_kurtosis_comps, metrics_df
end

# === Line Noise Identification ===

"""
    identify_line_noise_components(dat::ContinuousData, ica::InfoIca;
                                 line_freq::Real=50.0,
                                 freq_bandwidth::Real=1.0,
                                 z_threshold::Real=3.0,
                                 min_harmonic_power::Real=1.5,
                                 sample_selection::Function = samples(),
                                 interval_selection::Interval = times())

Identify ICA components with strong line noise characteristics using Welch periodograms
and spectral peakiness ratios relative to flanking baseline frequencies.
"""
function identify_line_noise_components(
    dat::ContinuousData,
    ica::InfoIca;
    sample_selection = samples(),
    interval_selection::Interval = times(),
    line_freq::Real = 50.0,
    freq_bandwidth::Real = 1.0,
    z_threshold::Real = 3.0,
    min_harmonic_power::Real = 1.5,
    _precomputed_components::Union{Nothing,Tuple{Matrix{Float64},Int}} = nothing,
    _precomputed_samples::Union{Nothing,Vector{Int}} = nothing,
)
    nyquist = dat.sample_rate / 2.0
    if line_freq >= nyquist
        @minimal_warning "Line frequency ($line_freq Hz) is at or above the Nyquist frequency ($nyquist Hz). Skipping line noise detection."
        return Int[], _empty_line_noise_metrics()
    end

    samples_to_use = if !isnothing(_precomputed_samples)
        _precomputed_samples
    else
        combined_sel = _combine_interval_sample(interval_selection, sample_selection)
        get_selected_samples(dat, combined_sel)
    end
    if isempty(samples_to_use)
        @minimal_warning "No samples remaining after applying exclude criteria. Cannot identify line noise components."
        return Int[], _empty_line_noise_metrics()
    end

    components, n_components = if !isnothing(_precomputed_components)
        _precomputed_components
    else
        _prepare_ica_data_matrix(dat, ica, samples_to_use)
    end

    n_samples = size(components, 2)
    components_t = Matrix{Float64}(components')

    # Welch periodogram with Hann window and 50% overlap using all samples
    n_per_seg = min(n_samples, nextpow(2, max(64, round(Int, 2.0 * dat.sample_rate))))
    n_overlap = div(n_per_seg, 2)

    first_pgram = DSP.welch_pgram(view(components_t, :, 1), n_per_seg, n_overlap; fs = dat.sample_rate, window = DSP.hanning)
    freqs = DSP.freq(first_pgram)
    n_freqs = length(freqs)

    psd = Matrix{Float64}(undef, n_freqs, n_components)
    psd[:, 1] .= DSP.power(first_pgram)

    for i = 2:n_components
        pgram = DSP.welch_pgram(view(components_t, :, i), n_per_seg, n_overlap; fs = dat.sample_rate, window = DSP.hanning)
        psd[:, i] .= DSP.power(pgram)
    end

    flank_bandwidth = 5.0
    line_band = findall(abs.(freqs .- line_freq) .<= freq_bandwidth)
    line_flank_band =
        findall((abs.(freqs .- line_freq) .> freq_bandwidth) .& (abs.(freqs .- line_freq) .<= freq_bandwidth + flank_bandwidth))

    harmonic2_freq = line_freq * 2.0
    harmonic2_band = findall(abs.(freqs .- harmonic2_freq) .<= freq_bandwidth)
    harmonic2_flank_band =
        findall((abs.(freqs .- harmonic2_freq) .> freq_bandwidth) .& (abs.(freqs .- harmonic2_freq) .<= freq_bandwidth + flank_bandwidth))

    metrics =
        NamedTuple{(:Component, :line_power, :flanking_power, :power_ratio, :harmonic_ratio),Tuple{Int,Float64,Float64,Float64,Float64}}[]

    for i = 1:n_components
        line_power = !isempty(line_band) ? mean(view(psd, line_band, i)) : 0.0
        flanking_power = !isempty(line_flank_band) ? mean(view(psd, line_flank_band, i)) : mean(view(psd, :, i))
        power_ratio = line_power / (flanking_power + eps())

        harmonic_ratio = if !isempty(harmonic2_band)
            h2_power = mean(view(psd, harmonic2_band, i))
            h2_flank = !isempty(harmonic2_flank_band) ? mean(view(psd, harmonic2_flank_band, i)) : flanking_power
            h2_power / (h2_flank + eps())
        else
            NaN
        end

        push!(
            metrics,
            (
                Component = i,
                line_power = line_power,
                flanking_power = flanking_power,
                power_ratio = power_ratio,
                harmonic_ratio = harmonic_ratio,
            ),
        )
    end

    metrics_df = DataFrame(metrics)
    power_ratio_z = _robust_z(metrics_df.power_ratio)
    metrics_df[!, :power_ratio_zscore] = power_ratio_z

    high_power_ratio_mask = power_ratio_z .> z_threshold
    # Only require 2nd harmonic if it is observable below Nyquist
    has_harmonic_mask = ifelse.(isnan.(metrics_df.harmonic_ratio), true, metrics_df.harmonic_ratio .> min_harmonic_power)

    line_noise_comps = findall(high_power_ratio_mask .& has_harmonic_mask)
    sort!(line_noise_comps)

    return line_noise_comps, metrics_df
end

# === ArtifactComponents Structure and Display ===

"""
    ArtifactComponents

A structure holding all identified artifact components from ICA analysis.

# Fields
- `eog::Dict{Symbol, Vector{Int}}`: Dictionary with :vEOG and :hEOG keys containing identified eye movement components
- `ecg::Vector{Int}`: Vector of identified ECG/heartbeat components  
- `line_noise::Vector{Int}`: Vector of identified line noise components
- `channel_noise::Vector{Int}`: Vector of identified high spatial kurtosis (channel noise) components
"""
struct ArtifactComponents
    eog::Dict{Symbol,Vector{Int}}
    ecg::Vector{Int}
    line_noise::Vector{Int}
    channel_noise::Vector{Int}
end

"""
    combine_artifact_components(eog_comps, ecg_comps, line_noise_comps, channel_noise_comps)

Combine all identified artifact components into a single `ArtifactComponents` structure.
"""
function combine_artifact_components(
    eog_comps::Dict{Symbol,Vector{Int}},
    ecg_comps::Vector{Int},
    line_noise_comps::Vector{Int},
    channel_noise_comps::Vector{Int},
)
    return ArtifactComponents(eog_comps, ecg_comps, line_noise_comps, channel_noise_comps)
end

"""
    get_all_ica_components(artifacts::ArtifactComponents)
    get_all_components(artifacts::ArtifactComponents)

Get all unique artifact component indices as a single sorted vector.
"""
function get_all_ica_components(artifacts::ArtifactComponents)
    all_comps = Set{Int}()
    for (_, comps) in artifacts.eog
        union!(all_comps, comps)
    end
    union!(all_comps, artifacts.ecg)
    union!(all_comps, artifacts.line_noise)
    union!(all_comps, artifacts.channel_noise)
    return sort(collect(all_comps))
end

const get_all_components = get_all_ica_components

function Base.show(io::IO, artifacts::ArtifactComponents)
    all_comps = get_all_ica_components(artifacts)
    print(io, "ArtifactComponents($(length(all_comps)) components: $(all_comps))")
end

function Base.show(io::IO, ::MIME"text/plain", artifacts::ArtifactComponents)
    all_comps = get_all_ica_components(artifacts)
    println(io, "ICA Artifact Components")
    println(io, "  vEOG: $(get(artifacts.eog, :vEOG, Int[]))")
    println(io, "  hEOG: $(get(artifacts.eog, :hEOG, Int[]))")
    println(io, "  ECG: $(artifacts.ecg)")
    println(io, "  Line Noise: $(artifacts.line_noise)")
    println(io, "  Channel Noise: $(artifacts.channel_noise)")
    print(io, "  All: $all_comps")
end

# === Unified Public API ===

"""
    identify_components(dat::ContinuousData, ica::InfoIca; method = :correlation,
                        sample_selection::Function = samples(), interval_selection::Interval = times(),
                        eog::NamedTuple = (;), ecg::NamedTuple = (;),
                        line_noise::NamedTuple = (;), channel_noise::NamedTuple = (;),
                        kwargs...)

Identify all types of artifact components in one unified call.
Supports both top-level general keyword arguments (e.g. `z_threshold = 2.5`) and per-detector
groups (`eog = (; z_threshold = 2.0)`, `ecg = (; min_bpm = 50)`).

# Arguments
- `dat::ContinuousData`: The continuous data
- `ica::InfoIca`: The ICA result object
- `method::Symbol`: Component identification method (default: `:correlation`)
- `sample_selection::Function`: Sample selection filter (default: `samples()`)
- `interval_selection::Interval`: Time interval filter (default: `times()`)
- `eog::NamedTuple`: Parameters specifically forwarded to `identify_eog_components`
- `ecg::NamedTuple`: Parameters specifically forwarded to `identify_ecg_components`
- `line_noise::NamedTuple`: Parameters specifically forwarded to `identify_line_noise_components`
- `channel_noise::NamedTuple`: Parameters specifically forwarded to `identify_spatial_kurtosis_components`
- `kwargs...`: Top-level keyword arguments automatically routed to applicable detectors.

# Returns
- `ArtifactComponents`: Combined structure containing all identified artifact components
- `Dict{Symbol, DataFrame}`: Dictionary containing metrics DataFrames for each component type
"""
function identify_components(
    dat::ContinuousData,
    ica::InfoIca;
    method::Symbol = :correlation,
    sample_selection = samples(),
    interval_selection::Interval = times(),
    eog::NamedTuple = (;),
    ecg::NamedTuple = (;),
    line_noise::NamedTuple = (;),
    channel_noise::NamedTuple = (;),
    kwargs...,
)
    if method == :correlation
        return _identify_components_correlation(
            dat,
            ica;
            sample_selection = sample_selection,
            interval_selection = interval_selection,
            eog = eog,
            ecg = ecg,
            line_noise = line_noise,
            channel_noise = channel_noise,
            kwargs...,
        )
    else
        @minimal_error("Unknown component identification method: :$method. Supported methods: :correlation")
    end
end

function _identify_components_correlation(
    dat::ContinuousData,
    ica::InfoIca;
    sample_selection = samples(),
    interval_selection::Interval = times(),
    eog::NamedTuple = (;),
    ecg::NamedTuple = (;),
    line_noise::NamedTuple = (;),
    channel_noise::NamedTuple = (;),
    kwargs...,
)
    # Pre-compute shared data ONCE
    combined_sel = _combine_interval_sample(interval_selection, sample_selection)
    selected_samples = get_selected_samples(dat, combined_sel)
    precomputed = _prepare_ica_data_matrix(dat, ica, selected_samples)

    # Route top-level keyword arguments to relevant detectors
    eog_keys = (:vEOG_channel, :hEOG_channel, :z_threshold, :min_correlation, :two_step)
    ecg_keys = (:min_bpm, :max_bpm, :min_prominence_std, :min_peaks, :max_ibi_std_s, :min_peak_ratio)
    line_noise_keys = (:line_freq, :freq_bandwidth, :z_threshold, :min_harmonic_power)
    channel_noise_keys = (:z_threshold,)

    all_known = union(eog_keys, ecg_keys, line_noise_keys, channel_noise_keys)
    unknown = [k for k in keys(kwargs) if !(k in all_known)]
    if !isempty(unknown)
        @minimal_warning "Unrecognized keyword arguments passed to identify_components: $(join(unknown, ", "))."
    end

    eog_merged = merge(NamedTuple(k => v for (k, v) in pairs(kwargs) if k in eog_keys), eog)
    ecg_merged = merge(NamedTuple(k => v for (k, v) in pairs(kwargs) if k in ecg_keys), ecg)
    line_noise_merged = merge(NamedTuple(k => v for (k, v) in pairs(kwargs) if k in line_noise_keys), line_noise)
    channel_noise_merged = merge(NamedTuple(k => v for (k, v) in pairs(kwargs) if k in channel_noise_keys), channel_noise)

    # Identify EOG components
    eog_comps, eog_metrics_df = identify_eog_components(
        dat,
        ica;
        sample_selection = sample_selection,
        interval_selection = interval_selection,
        _precomputed_components = precomputed,
        _precomputed_samples = selected_samples,
        eog_merged...,
    )

    # Identify ECG components
    ecg_comps, ecg_metrics_df = identify_ecg_components(
        dat,
        ica;
        sample_selection = sample_selection,
        interval_selection = interval_selection,
        _precomputed_components = precomputed,
        _precomputed_samples = selected_samples,
        ecg_merged...,
    )

    # Identify line noise components
    line_noise_comps, line_noise_metrics_df = identify_line_noise_components(
        dat,
        ica;
        sample_selection = sample_selection,
        interval_selection = interval_selection,
        _precomputed_components = precomputed,
        _precomputed_samples = selected_samples,
        line_noise_merged...,
    )

    # Identify channel noise components (spatial kurtosis)
    channel_noise_comps, channel_noise_metrics_df = identify_spatial_kurtosis_components(ica; channel_noise_merged...)

    # Combine all components
    artifacts = combine_artifact_components(eog_comps, ecg_comps, line_noise_comps, channel_noise_comps)

    # Combine all metrics
    metrics = Dict{Symbol,DataFrame}(
        :eog_metrics => eog_metrics_df,
        :ecg_metrics => ecg_metrics_df,
        :line_noise_metrics => line_noise_metrics_df,
        :channel_noise_metrics => channel_noise_metrics_df,
    )

    return artifacts, metrics
end
