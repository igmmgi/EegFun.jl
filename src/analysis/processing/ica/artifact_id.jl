

"""
    identify_eog_components(dat::ContinuousData, ica::InfoIca;
                          vEOG_channel::Symbol=:vEOG,
                          hEOG_channel::Symbol=:hEOG,
                          z_threshold::Float64=3.0,
                          two_step::Bool=true,
                          sample_selection::Function = samples())

Identify ICA components potentially related to eye movements based on z-scored correlation with EOG channels.

Uses a two-step approach by default:
1. **Step 1**: Calculate correlations between all ICA components and EOG channels, compute standard z-scores, 
   and identify primary EOG components (those exceeding both the z-score and correlation thresholds).
2. **Step 2**: Calculate spatial correlations (topographies) between remaining components and primary EOG components. 
   Identify secondary EOG components based on spatial correlation threshold. This approach validates that secondary 
   components share similar spatial patterns with primary components, which is more aligned with ICA principles.

This two-step approach helps detect secondary EOG components that share spatial topographies with primary components, 
even if their temporal correlations with EOG channels are not as strong.

# Arguments
- `dat::ContinuousData`: The continuous data containing EOG channels.
- `ica::InfoIca`: The ICA result object.

# Keyword Arguments
- `vEOG_channel::Symbol`: Name of the vertical EOG channel (default: :vEOG).
- `hEOG_channel::Symbol`: Name of the horizontal EOG channel (default: :hEOG).
- `z_threshold::Float64`: Absolute Z-score threshold for step1 identification (default: 3.0). Primary components must exceed 
  this z-score threshold to be identified in step1.
- `min_correlation::Float64`: For step1, minimum absolute correlation threshold with EOG channels (default: 0.6). Primary 
  components must exceed this correlation threshold in addition to the z-score threshold. For step2, this is the minimum 
  spatial correlation threshold with primary component topographies required to identify secondary components.
- `two_step::Bool`: If true (default), use two-step identification: step1 identifies primary components based on z-score and 
  EOG correlation, step2 identifies secondary components based on spatial correlation with primary components. If false, 
  use single-step identification with standard z-scores and correlation thresholds on all components simultaneously.
- `sample_selection::Function`: Function to select samples from `dat.data`. Defaults to `samples()`. Only 
  selected samples are used for correlation calculation.

# Returns
- `Dict{Symbol, Vector{Int}}`: Dictionary containing:
  - `:vEOG`: Vector of component indices identified for vertical eye movements.
  - `:hEOG`: Vector of component indices identified for horizontal eye movements.
- `DataFrame`: DataFrame containing detailed metrics per component:
  - `:Component`: Component index (1 to n_components)
  - `:vEOG_corr`: Absolute correlation with vertical EOG channel
  - `:vEOG_zscore` or `:vEOG_zscore_step1`: Z-score of correlation with vertical EOG channel.
    When `two_step=true`, only step1 z-scores are provided.
  - `:vEOG_spatial_corr`: (Only when `two_step=true`) Maximum spatial correlation (topography) with any primary vEOG 
    component. `NaN` for primary vEOG components identified in step1.
  - `:vEOG_spatial_corr_z`: Z-score of `:vEOG_spatial_corr`.
  - `:vEOG_temporal_corr`: (Only when `two_step=true`) Maximum lagged temporal correlation with any primary vEOG 
    component. `NaN` for primary vEOG components identified in step1.
  - `:vEOG_temporal_corr_z`: Z-score of `:vEOG_temporal_corr`.
  - `:hEOG_corr`: Absolute correlation with horizontal EOG channel
  - `:hEOG_zscore` or `:hEOG_zscore_step1`: Z-score of correlation with horizontal EOG channel.
    When `two_step=true`, only step1 z-scores are provided.
  - `:hEOG_spatial_corr`: (Only when `two_step=true`) Maximum spatial correlation (topography) with any primary hEOG 
    component. `NaN` for primary hEOG components identified in step1.
  - `:hEOG_spatial_corr_z`: Z-score of `:hEOG_spatial_corr`.
  - `:hEOG_temporal_corr`: (Only when `two_step=true`) Maximum lagged temporal correlation with any primary hEOG 
    component. `NaN` for primary hEOG components identified in step1.
  - `:hEOG_temporal_corr_z`: Z-score of `:hEOG_temporal_corr`.

# Examples
```julia
# Basic usage with default two-step approach
eog_comps, metrics = identify_eog_components(dat, ica_result)

# Use single-step identification
eog_comps, metrics = identify_eog_components(dat, ica_result, two_step=false)

# Custom threshold and sample selection
eog_comps, metrics = identify_eog_components(dat, ica_result, 
    z_threshold=2.5,
    sample_selection=samples_not(:is_extreme_value_100))
```
"""
function identify_eog_components(
    dat::ContinuousData,
    ica::InfoIca;
    vEOG_channel::Symbol = :vEOG,
    hEOG_channel::Symbol = :hEOG,
    sample_selection::Function = samples(),
    interval_selection::Interval = times(),
    z_threshold::Real = 3.0,
    min_correlation::Float64 = 0.5,
    two_step::Bool = true,
    _precomputed_components::Union{Nothing,Tuple{Matrix{Float64},Int}} = nothing,
    _precomputed_samples::Union{Nothing,Vector{Int}} = nothing,
)

    # Check basic inputs - return nothing as first value if EOG channels are missing.
    # NOTE: Currently requires both vEOG and hEOG to be present.
    if !(vEOG_channel in propertynames(dat.data)) || !(hEOG_channel in propertynames(dat.data))
        return nothing, DataFrame()
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
        return Dict(:vEOG => Int[], :hEOG => Int[]), DataFrame()
    end

    # Function to calculate initial EOG-component correlations 
    """Compute absolute correlation of each ICA component with the given EOG signal."""
    function calculate_correlations(eog_signal)
        corrs = zeros(n_components)
        comp_buf = Vector{Float64}(undef, size(components, 2))
        for comp_idx = 1:n_components
            comp_buf .= @view components[comp_idx, :]
            corrs[comp_idx] = abs(cor(comp_buf, eog_signal))
        end
        return corrs
    end

    # Function to calculate spatial correlations between identified vEOG/hEOG components and remaining components (used in Step2)
    """Compute spatial (topography) correlations between remaining and primary EOG components."""
    function calculate_spatial_correlations(remaining_components, primary_components, mixing_matrix, min_corr_threshold)

        spatial_corrs = fill(NaN, n_components)
        secondary_components = Int[]

        if !isempty(remaining_components) && !isempty(primary_components)
            # Calculate spatial correlations (topographies) between remaining components and each primary component taking max correlation
            spatial_corrs_remaining = zeros(length(remaining_components))
            for (idx, comp_idx) in enumerate(remaining_components)
                comp_topography = mixing_matrix[:, comp_idx]
                max_corr = 0.0
                for primary_comp in primary_components
                    primary_topography = mixing_matrix[:, primary_comp]
                    corr_val = abs(cor(comp_topography, primary_topography))
                    max_corr = max(max_corr, corr_val)
                end
                spatial_corrs_remaining[idx] = max_corr
            end

            # Map spatial correlations back to original component indices
            spatial_corrs[remaining_components] = spatial_corrs_remaining

            # Find secondary components (spatial correlation > min_correlation)
            secondary_idx = findall(spatial_corrs_remaining .> min_corr_threshold)
            secondary_components = remaining_components[secondary_idx]
        end

        return spatial_corrs, secondary_components
    end

    # Function to calculate lagged correlations between components without huge crosscor allocations
    """Compute maximum absolute lagged cross-correlation between remaining and primary components."""
    function calculate_lagged_correlations(remaining_components, primary_components, comp_matrix, max_lag_samples, lag_step)
        lagged_corrs = fill(NaN, n_components)
        if isempty(remaining_components) || isempty(primary_components)
            return lagged_corrs
        end

        lags = (-max_lag_samples):lag_step:max_lag_samples
        n_s = size(comp_matrix, 2)

        # Pre-allocate contiguous buffers to eliminate cache misses from row-major stride
        comp_buf = Vector{Float64}(undef, n_s)
        eog_buf = Vector{Float64}(undef, n_s)

        for comp_idx in remaining_components
            comp_buf .= @view comp_matrix[comp_idx, :]
            m_c = mean(comp_buf)
            s_c = std(comp_buf, corrected = false)

            max_corr = 0.0
            for eog_comp_idx in primary_components
                eog_buf .= @view comp_matrix[eog_comp_idx, :]
                m_e = mean(eog_buf)
                s_e = std(eog_buf, corrected = false)
                denom = n_s * s_c * s_e

                # Manual zero-allocation cross correlation for the few specific lags
                max_corr_val = 0.0
                for lag in lags
                    sum_xy = 0.0
                    @inbounds @simd for i = max(1, 1-lag):min(n_s, n_s-lag)
                        sum_xy += (eog_buf[i] - m_e) * (comp_buf[i+lag] - m_c)
                    end
                    c = abs(sum_xy / denom)
                    if c > max_corr_val
                        max_corr_val = c
                    end
                end

                if max_corr_val > max_corr
                    max_corr = max_corr_val
                end
            end

            lagged_corrs[comp_idx] = max_corr
        end

        return lagged_corrs
    end

    # Get EOG signals for valid samples only
    vEOG = dat.data[selected_samples, vEOG_channel]
    hEOG = dat.data[selected_samples, hEOG_channel]

    # Calculate components for valid samples (skip if pre-computed)
    components, n_components = if !isnothing(_precomputed_components)
        _precomputed_components
    else
        _prepare_ica_data_matrix(dat, ica, selected_samples)
    end

    # Step 1: Always calculate correlations and z-scores (used in both single-step and two-step modes)
    vEOG_corrs = calculate_correlations(vEOG)
    hEOG_corrs = calculate_correlations(hEOG)

    vEOG_corr_z = StatsBase.zscore(vEOG_corrs)
    hEOG_corr_z = StatsBase.zscore(hEOG_corrs)

    # Identify primary components (components meeting both z-score and correlation thresholds)
    primary_vEOG = findall((abs.(vEOG_corr_z) .> z_threshold) .& (abs.(vEOG_corrs) .> min_correlation))
    primary_hEOG = findall((abs.(hEOG_corr_z) .> z_threshold) .& (abs.(hEOG_corrs) .> min_correlation))

    # Step 2: Always calculate spatial and temporal correlations between identified components and remaining
    remaining_vEOG = setdiff(1:n_components, primary_vEOG)
    remaining_hEOG = setdiff(1:n_components, primary_hEOG)

    # Calculate vEOG/hEOG step2: spatial correlations with primary vEOG/hEOG components
    vEOG_spatial_corr, _ = calculate_spatial_correlations(remaining_vEOG, primary_vEOG, ica.mixing, min_correlation)
    hEOG_spatial_corr, _ = calculate_spatial_correlations(remaining_hEOG, primary_hEOG, ica.mixing, min_correlation)

    # Lag range: +- 100ms  
    # Convert to samples based on sample rate
    max_lag_samples = round(Int, 100.0 * dat.sample_rate / 1000.0) # 100 ms lag
    lag_step = max(1, round(Int, 5.0 * dat.sample_rate / 1000.0)) # 5 ms step

    # Calculate lagged correlations (zero-allocation, uses pre-computed component means)
    vEOG_temporal_corr = calculate_lagged_correlations(remaining_vEOG, primary_vEOG, components, max_lag_samples, lag_step)
    hEOG_temporal_corr = calculate_lagged_correlations(remaining_hEOG, primary_hEOG, components, max_lag_samples, lag_step)

    # Helper function to compute z-scores only for remaining components
    """Z-score a correlation vector using only the values at `remaining_indices`."""
    function zscore_for_remaining(corr_vector, remaining_indices)
        z = fill(NaN, length(corr_vector))
        if length(remaining_indices) >= 2
            remaining_vals = corr_vector[remaining_indices]
            if std(remaining_vals) > 0.0
                z[remaining_indices] = StatsBase.zscore(remaining_vals)
            end
        end
        return z
    end

    # Calculate z-scores for spatial and temporal correlations
    vEOG_spatial_corr_z = zscore_for_remaining(vEOG_spatial_corr, remaining_vEOG)
    vEOG_temporal_corr_z = zscore_for_remaining(vEOG_temporal_corr, remaining_vEOG)
    hEOG_spatial_corr_z = zscore_for_remaining(hEOG_spatial_corr, remaining_hEOG)
    hEOG_temporal_corr_z = zscore_for_remaining(hEOG_temporal_corr, remaining_hEOG)

    # Identify secondary components: require BOTH correlation > min_correlation AND z-score > z_threshold
    """Return indices of remaining components exceeding both correlation and z-score thresholds."""
    function identify_secondary_components(remaining_components, corr_vector, corr_z_vector, min_corr_threshold, z_threshold)
        isempty(remaining_components) && return Int[]
        corr_mask = corr_vector[remaining_components] .> min_corr_threshold
        z_mask = abs.(corr_z_vector[remaining_components]) .> z_threshold
        combined_mask = corr_mask .& z_mask
        return remaining_components[combined_mask]
    end

    # Identify secondary components for vEOG: spatial OR temporal (each requires both corr and z-score)
    secondary_vEOG_spatial =
        identify_secondary_components(remaining_vEOG, vEOG_spatial_corr, vEOG_spatial_corr_z, min_correlation, z_threshold)
    secondary_vEOG_temporal =
        identify_secondary_components(remaining_vEOG, vEOG_temporal_corr, vEOG_temporal_corr_z, min_correlation, z_threshold)

    # Identify secondary components for hEOG: spatial OR temporal (each requires both corr and z-score)
    secondary_hEOG_spatial =
        identify_secondary_components(remaining_hEOG, hEOG_spatial_corr, hEOG_spatial_corr_z, min_correlation, z_threshold)
    secondary_hEOG_temporal =
        identify_secondary_components(remaining_hEOG, hEOG_temporal_corr, hEOG_temporal_corr_z, min_correlation, z_threshold)

    # Combine spatial and temporal secondary components
    secondary_vEOG = union(secondary_vEOG_spatial, secondary_vEOG_temporal)
    secondary_hEOG = union(secondary_hEOG_spatial, secondary_hEOG_temporal)

    # Combine primary and secondary components (secondary only added if two_step is true)
    identified_vEOG = two_step ? union(primary_vEOG, secondary_vEOG) : primary_vEOG
    identified_hEOG = two_step ? union(primary_hEOG, secondary_hEOG) : primary_hEOG

    # DataFrame metrics
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

    sort!(identified_vEOG)
    sort!(identified_hEOG)

    result_dict = Dict{Symbol,Vector{Int}}(:vEOG => identified_vEOG, :hEOG => identified_hEOG)

    return result_dict, metrics_df
end

"""
    identify_ecg_components(dat::ContinuousData, ica::InfoIca;
                              min_bpm::Real=50, max_bpm::Real=110,
                              min_prominence_std::Real=4,
                              min_peaks::Int=100,
                              max_ibi_std_s::Real=0.2,
                              min_peak_ratio::Real=0.5,
                              sample_selection::Function = samples(),
    interval_selection::Interval = times(),
)

Identify ICA components potentially related to ECG artifacts based on peak detection
and interval regularity, using only samples consistent with ICA calculation.

# Arguments
- `dat::ContinuousData`: The continuous data (needed for sampling rate and sample selection columns).
- `ica::InfoIca`: The ICA result object.

# Keyword Arguments
- `min_bpm::Real`: Minimum plausible heart rate in beats per minute (default: 40).
- `max_bpm::Real`: Maximum plausible heart rate in beats per minute (default: 120).
- `min_prominence_std::Real`: Minimum peak prominence in standard deviations above mean for z-scored time series (default: 4).
- `min_peaks::Int`: Minimum number of valid inter-beat intervals required (default: 100). 
  Note: Since `num_valid_ibis` is the number of valid intervals between peaks (which is `num_peaks - 1` if all are valid),
  the check uses `num_valid_ibis >= (min_peaks - 1)` to account for this relationship.
- `max_ibi_std_s::Real`: Maximum standard deviation of the inter-beat intervals (in seconds) for component to be flagged (default: 0.2).
- `min_peak_ratio::Real`: Minimum ratio of valid inter-beat intervals to total inter-peak intervals (default: 0.5). 
  This ensures that a sufficient proportion of detected peaks fall within the plausible heart rate range.
- `sample_selection::Function`: Function to select samples from `dat.data`. Defaults to `samples()`. Only 
  selected samples are used for component calculation and peak detection.

# Returns
- `Vector{Int}`: Sorted vector of indices identified as potential ECG components.
- `DataFrame`: DataFrame containing metrics for each component (calculated on the selected samples):
  - `:Component`: Component index (1 to n).
  - `:num_peaks`: Number of detected prominent peaks.
  - `:num_valid_ibis`: Number of inter-beat intervals within the plausible BPM range.
  - `:mean_ibi_s`: Mean inter-beat interval in seconds (if num_valid_ibis > 0).
  - `:std_ibi_s`: Standard deviation of inter-beat intervals in seconds (if num_valid_ibis > 1).
  - `:peak_ratio`: Ratio of valid inter-beat intervals to total inter-peak intervals.
  - `:heart_rate_bpm`: Estimated heart rate in beats per minute (if mean_ibi_s is valid).
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
    sample_selection::Function = samples(),
    interval_selection::Interval = times(),
    _precomputed_components::Union{Nothing,Tuple{Matrix{Float64},Int}} = nothing,
    _precomputed_samples::Union{Nothing,Vector{Int}} = nothing,
)

    # Data Preparation (skip if pre-computed)
    selected_samples = if !isnothing(_precomputed_samples)
        _precomputed_samples
    else
        combined_sel = _combine_interval_sample(interval_selection, sample_selection)
        get_selected_samples(dat, combined_sel)
    end
    if isempty(selected_samples)
        @minimal_warning "No samples remaining after applying exclude criteria."
        return Int[], DataFrame()
    end

    # Calculate components for valid samples (skip if pre-computed)
    components_subset, n_components = if !isnothing(_precomputed_components)
        _precomputed_components
    else
        _prepare_ica_data_matrix(dat, ica, selected_samples)
    end

    # Early return if no components
    if n_components == 0
        return Int[],
        DataFrame(
            Component = Int[],
            num_peaks = Int[],
            num_valid_ibis = Int[],
            mean_ibi_s = Float64[],
            std_ibi_s = Float64[],
            peak_ratio = Float64[],
            heart_rate_bpm = Float64[],
        )
    end

    # Convert BPM to plausible IBI range
    min_ibi_s = 60.0 / max_bpm
    max_ibi_s = 60.0 / min_bpm

    metrics = []
    identified_ecg = Int[]

    # Pre-allocate buffer for z-scoring
    n_samples = size(components_subset, 2)
    ts_zscored = Vector{Float64}(undef, n_samples)

    for comp_idx = 1:n_components
        component_ts = @view components_subset[comp_idx, :]

        # Z-score the time series for consistent peak detection across components
        m = mean(component_ts)
        s = std(component_ts)
        ts_zscored .= (component_ts .- m) ./ s

        # Find prominent peaks 
        peak_indices = _find_peaks(ts_zscored; min_prominence_std = min_prominence_std)
        num_peaks = length(peak_indices)

        # Calculate IBI metrics
        num_valid_ibis, mean_ibi, std_ibi, peak_ratio = _calculate_ibi_metrics(peak_indices, dat.sample_rate, min_ibi_s, max_ibi_s)

        # Check if component meets ECG criteria
        has_sufficient_ibis = num_valid_ibis >= (min_peaks - 1)
        has_low_ibi_variability = std_ibi <= max_ibi_std_s
        has_good_peak_ratio = peak_ratio >= min_peak_ratio

        is_ecg = has_sufficient_ibis && has_low_ibi_variability && has_good_peak_ratio
        is_ecg && push!(identified_ecg, comp_idx)

        # Calculate heart rate if we have valid IBI
        if isnan(mean_ibi) || mean_ibi <= 0
            heart_rate_bpm = NaN
        else
            heart_rate_bpm = 60.0 / mean_ibi
        end

        # Store metrics
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

    # Build metrics DataFrame
    metrics_df = DataFrame(metrics)

    sort!(identified_ecg)

    return identified_ecg, metrics_df

end

"""
    _calculate_ibi_metrics(peak_indices::Vector{Int}, sample_rate::Real, min_ibi_s::Real, max_ibi_s::Real)

Calculate inter-beat interval (IBI) metrics from peak indices.

# Arguments
- `peak_indices::Vector{Int}`: Indices of detected peaks
- `sample_rate::Real`: Sampling rate in Hz
- `min_ibi_s::Real`: Minimum valid IBI in seconds
- `max_ibi_s::Real`: Maximum valid IBI in seconds

# Returns
- `num_valid_ibis::Int`: Number of valid IBIs
- `mean_ibi::Float64`: Mean IBI in seconds (NaN if insufficient data)
- `std_ibi::Float64`: Standard deviation of IBIs in seconds (NaN if insufficient data)
- `peak_ratio::Float64`: Ratio of valid IBIs to total inter-peak intervals
"""
function _calculate_ibi_metrics(peak_indices::Vector{Int}, sample_rate::Real, min_ibi_s::Real, max_ibi_s::Real)

    num_peaks = length(peak_indices)
    num_peaks < 2 && return 0, NaN, NaN, 0.0

    # Calculate IBIs
    ibis_s = diff(peak_indices) ./ sample_rate

    # Filter valid IBIs
    valid_ibi_mask = (ibis_s .>= min_ibi_s) .& (ibis_s .<= max_ibi_s)
    valid_ibis = ibis_s[valid_ibi_mask]
    num_valid_ibis = length(valid_ibis)

    # Calculate peak ratio
    peak_ratio = num_valid_ibis / (num_peaks - 1)

    # Calculate statistics
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
    identify_spatial_kurtosis_components(ica::InfoIca; z_threshold::Float64 = 3.0)

Identify ICA components with high spatial kurtosis (localized, spot-like activity).

Spatial kurtosis measures how localized the component's topography is. High spatial kurtosis
indicates that the component's activity is concentrated in a small number of channels (spot-like),
which is characteristic of channel noise or artifacts.

# Arguments
- `ica::InfoIca`: The ICA result object.

# Keyword Arguments
- `z_threshold::Float64`: Z-score threshold for identifying high spatial kurtosis components (default: 3.0).

# Returns
- `Vector{Int}`: Indices of components with high spatial kurtosis.
- `DataFrame`: DataFrame containing spatial kurtosis values and z-scores for all components.
"""
function identify_spatial_kurtosis_components(ica::InfoIca; z_threshold::Real = 3.0)

    # Calculate spatial kurtosis for each component's weights
    n_components = size(ica.mixing, 2)
    spatial_kurtosis = Float64[]
    for i = 1:n_components
        weights = ica.mixing[:, i]
        k = kurtosis(weights)
        push!(spatial_kurtosis, k)
    end

    # Calculate z-scores of spatial kurtosis values
    z_spatial_kurtosis = StatsBase.zscore(spatial_kurtosis)

    # Identify components with high spatial kurtosis (using z-scores)
    high_kurtosis_comps = findall(z_spatial_kurtosis .> z_threshold)  # Only positive deviations (localized activity)
    sort!(high_kurtosis_comps)

    # Create metrics DataFrame
    metrics_df = DataFrame(:Component => 1:n_components, :spatial_kurtosis => spatial_kurtosis, :z_spatial_kurtosis => z_spatial_kurtosis)

    return high_kurtosis_comps, metrics_df
end


"""
    identify_line_noise_components(dat::ContinuousData, ica::InfoIca;
                                 line_freq::Real=50.0,
                                 freq_bandwidth::Real=1.0,
                                 z_threshold::Float64=3.0,
                                 min_harmonic_power::Real=0.5)

Identify ICA components with strong line noise characteristics using spectral peakiness.

Spectral peakiness compares power at the line frequency to nearby flanking frequencies
(±5 Hz excluding the line band), rather than to the entire broadband spectrum. This avoids
false positives from components with generally elevated power that happens to include the
line frequency range.

# Arguments
- `dat::ContinuousData`: The continuous data.
- `ica::InfoIca`: The ICA result object.

# Keyword Arguments
- `sample_selection::Function`: Function to select samples (default: `samples()`).
- `line_freq::Real`: Line frequency in Hz (default: 50.0 for European power).
- `freq_bandwidth::Real`: Bandwidth around line frequency to consider (default: 1.0 Hz).
- `z_threshold::Float64`: Z-score threshold for identifying line noise components (default: 3.0).
- `min_harmonic_power::Real`: Minimum peakiness ratio of 2nd harmonic (100 Hz vs. its flanking frequencies) (default: 1.5).

# Returns
- `Vector{Int}`: Indices of components with strong line noise characteristics.
- `DataFrame`: DataFrame containing spectral metrics for all components:
  - `:Component`: Component index
  - `:line_power`: Mean power at the line frequency band
  - `:flanking_power`: Mean power at flanking frequencies (used as baseline for peakiness)
  - `:power_ratio`: Spectral peakiness ratio (line_power / flanking_power)
  - `:harmonic_ratio`: Spectral peakiness ratio at the 2nd harmonic
  - `:power_ratio_zscore`: Z-score of power_ratio across components
"""
function identify_line_noise_components(
    dat::ContinuousData,
    ica::InfoIca;
    sample_selection::Function = samples(),
    interval_selection::Interval = times(),
    line_freq::Real = 50.0,
    freq_bandwidth::Real = 1.0,
    z_threshold::Real = 3.0,
    min_harmonic_power::Real = 1.5,
    _precomputed_components::Union{Nothing,Tuple{Matrix{Float64},Int}} = nothing,
    _precomputed_samples::Union{Nothing,Vector{Int}} = nothing,
)

    # Get samples to use (skip if pre-computed)
    samples_to_use = if !isnothing(_precomputed_samples)
        _precomputed_samples
    else
        combined_sel = _combine_interval_sample(interval_selection, sample_selection)
        get_selected_samples(dat, combined_sel)
    end
    if isempty(samples_to_use)
        @minimal_warning "No samples remaining after applying exclude criteria. Cannot identify line noise components."
        return Int[], DataFrame()
    end

    # Calculate components for valid samples (skip if pre-computed)
    components, n_components = if !isnothing(_precomputed_components)
        _precomputed_components
    else
        _prepare_ica_data_matrix(dat, ica, samples_to_use)
    end

    # Calculate power spectrum for each component
    # Use a reasonable FFT size (power of 2, capped at 2^16 points for performance)
    n_samples = size(components, 2)
    nfft = min(nextpow(2, n_samples), 65536)  # Properly cap at 2^16 points
    freqs = FFTW.rfftfreq(nfft, dat.sample_rate)
    n_freqs = length(freqs)
    psd = zeros(n_freqs, n_components)

    # Pre-allocate FFT buffer and plan FFT
    signal_fft = zeros(Float64, nfft)
    fft_plan = FFTW.plan_rfft(signal_fft)

    for i = 1:n_components
        signal = @view components[i, :]

        # Prepare signal for FFT: truncate or zero-pad to nfft points
        if n_samples > nfft
            signal_fft .= @view signal[1:nfft]
        elseif n_samples < nfft
            signal_fft[1:n_samples] .= signal
            signal_fft[(n_samples+1):end] .= 0.0
        else
            signal_fft .= signal
        end

        # Calculate power spectral density
        psd[:, i] .= abs2.(fft_plan * signal_fft)
    end

    # Define frequency bands
    # Flanking bandwidth for spectral peakiness (Hz on each side of line freq, excluding the line band itself)
    flank_bandwidth = 5.0

    # Line frequency band (narrow band around 50 Hz)
    line_band = findall(abs.(freqs .- line_freq) .<= freq_bandwidth)

    # Flanking bands: frequencies near the line frequency but excluding the line band
    # e.g., for 50 Hz with freq_bandwidth=1 and flank_bandwidth=5: uses 45-49 Hz and 51-55 Hz
    line_flank_band =
        findall((abs.(freqs .- line_freq) .> freq_bandwidth) .& (abs.(freqs .- line_freq) .<= freq_bandwidth + flank_bandwidth))

    # 2nd harmonic bands
    harmonic2_freq = line_freq * 2
    harmonic2_band = findall(abs.(freqs .- harmonic2_freq) .<= freq_bandwidth)
    harmonic2_flank_band =
        findall((abs.(freqs .- harmonic2_freq) .> freq_bandwidth) .& (abs.(freqs .- harmonic2_freq) .<= freq_bandwidth + flank_bandwidth))

    # Calculate metrics for each component
    metrics = []
    for i = 1:n_components
        # Get power at line frequency band
        line_power = mean(psd[line_band, i])

        # Spectral peakiness: compare line frequency power to nearby flanking frequencies
        flanking_power = !isempty(line_flank_band) ? mean(psd[line_flank_band, i]) : mean(psd[:, i])

        # Power ratio: how much does 50 Hz stand out above its immediate neighbours
        power_ratio = line_power / (flanking_power + eps())

        # Harmonic peakiness: does 100 Hz also stand out above its flanking frequencies
        harmonic_ratio = if !isempty(harmonic2_band)
            harmonic2_power = mean(psd[harmonic2_band, i])
            harmonic2_flanking = !isempty(harmonic2_flank_band) ? mean(psd[harmonic2_flank_band, i]) : flanking_power
            harmonic2_power / (harmonic2_flanking + eps())
        else
            NaN  # 2nd harmonic frequency above Nyquist
        end

        # Store metrics
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

    # Create metrics DataFrame
    metrics_df = DataFrame(metrics)

    # Calculate z-scores of power ratios
    power_ratio_z = StatsBase.zscore(metrics_df.power_ratio)
    metrics_df[!, :power_ratio_zscore] = power_ratio_z

    # Identify components with strong line noise characteristics
    # Step 1: High power ratio (z-score > threshold)
    high_power_ratio_mask = power_ratio_z .> z_threshold

    # Step 2: Must have significant 2nd harmonic (> min_harmonic_power)
    has_harmonic_mask = metrics_df.harmonic_ratio .> min_harmonic_power

    # Combine both criteria
    line_noise_comps = findall(high_power_ratio_mask .& has_harmonic_mask)
    sort!(line_noise_comps)

    return line_noise_comps, metrics_df

end

"""
    ArtifactComponents

A structure to hold all identified artifact components from ICA analysis.

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

Combine all identified artifact components into a single ArtifactComponents structure.

# Arguments
- `eog_comps::Dict{Symbol, Vector{Int}}`: EOG components dictionary
- `ecg_comps::Vector{Int}`: ECG components vector
- `line_noise_comps::Vector{Int}`: Line noise components vector
- `channel_noise_comps::Vector{Int}`: Channel noise components vector

# Returns
- `ArtifactComponents`: Combined structure containing all artifact components
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
    identify_components(dat::ContinuousData, ica::InfoIca; method = :correlation,
    sample_selection::Function = samples(), interval_selection::Interval = times(), kwargs...)

Identify all types of artifact components in one unified call.

# Arguments
- `dat::ContinuousData`: The continuous data
- `ica::InfoIca`: The ICA result object
- `method::Symbol`: Component identification method:
  - `:correlation` (default): Correlation-based detection (EOG, ECG, line noise, spatial kurtosis).
    The original EegFun approach using z-scored correlations with EOG/ECG channels,
    spectral analysis for line noise, and spatial kurtosis for channel noise.
- `sample_selection::Function`: Sample selection function (default: samples())
- `kwargs...`: Additional keyword arguments passed to individual identification functions

# Returns
- `ArtifactComponents`: Combined structure containing all identified artifact components
- `Dict{Symbol, DataFrame}`: Dictionary containing metrics DataFrames for each component type:
  - `:eog_metrics`: EOG component metrics
  - `:ecg_metrics`: ECG component metrics  
  - `:line_noise_metrics`: Line noise component metrics
  - `:channel_noise_metrics`: Channel noise component metrics

# Examples
```julia
# Default correlation-based identification
artifacts, metrics = identify_components(dat, ica_result)

# With custom sample selection
artifacts, metrics = identify_components(dat, ica_result, 
    sample_selection = samples_not(:is_extreme_value_100))
```
"""
function identify_components(
    dat::ContinuousData,
    ica::InfoIca;
    method::Symbol = :correlation,
    sample_selection::Function = samples(),
    interval_selection::Interval = times(),
    kwargs...,
)
    if method == :correlation
        return _identify_components_correlation(
            dat,
            ica;
            sample_selection = sample_selection,
            interval_selection = interval_selection,
            kwargs...,
        )
    else
        @minimal_error("Unknown component identification method: :$method. " * "Supported methods: :correlation")
    end
end

"""
    _identify_components_correlation(dat, ica; kwargs...) -> (ArtifactComponents, Dict)

The original correlation-based component identification method.
This is the internal implementation dispatched by `identify_components(method = :correlation)`.
"""
function _identify_components_correlation(
    dat::ContinuousData,
    ica::InfoIca;
    sample_selection::Function = samples(),
    interval_selection::Interval = times(),
    kwargs...,
)
    # Pre-compute shared data ONCE (avoids 3x redundant matrix operations)
    combined_sel = _combine_interval_sample(interval_selection, sample_selection)
    selected_samples = get_selected_samples(dat, combined_sel)
    precomputed = _prepare_ica_data_matrix(dat, ica, selected_samples)

    # Identify EOG components (pass through use_robust_zscore if provided)
    eog_comps, eog_metrics_df = identify_eog_components(
        dat,
        ica;
        sample_selection = sample_selection,
        _precomputed_components = precomputed,
        _precomputed_samples = selected_samples,
        kwargs...,
    )

    # Identify ECG components
    ecg_comps, ecg_metrics_df = identify_ecg_components(
        dat,
        ica;
        sample_selection = sample_selection,
        _precomputed_components = precomputed,
        _precomputed_samples = selected_samples,
        kwargs...,
    )

    # Identify line noise components
    line_noise_comps, line_noise_metrics_df = identify_line_noise_components(
        dat,
        ica;
        sample_selection = sample_selection,
        _precomputed_components = precomputed,
        _precomputed_samples = selected_samples,
        kwargs...,
    )

    # Identify channel noise components (spatial kurtosis - uses only ica.mixing, no data matrix needed)
    channel_noise_comps, channel_noise_metrics_df = identify_spatial_kurtosis_components(ica; kwargs...)

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

"""
    get_all_components(artifacts::ArtifactComponents)

Get all unique artifact component indices as a single vector.

# Arguments
- `artifacts::ArtifactComponents`: The artifact components structure

# Returns
- `Vector{Int}`: Sorted vector of all unique artifact component indices
"""
function get_all_ica_components(artifacts::ArtifactComponents)
    all_comps = Set{Int}()

    # Add EOG components
    for (_, comps) in artifacts.eog
        union!(all_comps, comps)
    end

    # Add other components
    union!(all_comps, artifacts.ecg)
    union!(all_comps, artifacts.line_noise)
    union!(all_comps, artifacts.channel_noise)

    return sort(collect(all_comps))
end



function Base.show(io::IO, artifacts::ArtifactComponents)
    # Get all components for and print short summary
    all_comps = get_all_ica_components(artifacts)
    println(io, "ICA Artifact Components")
    println(io, "vEOG: $(artifacts.eog[:vEOG])")
    println(io, "hEOG: $(artifacts.eog[:hEOG])")
    println(io, "ECG: $(artifacts.ecg)")
    println(io, "Line Noise: $(artifacts.line_noise)")
    println(io, "Channel Noise: $(artifacts.channel_noise)")
    println(io, "All: $all_comps")
end
