function _robust_z(x::AbstractVector{<:Real})
    n = length(x)
    n == 0 && return Float64[]
    n == 1 && return [0.0]
    med = median(x)
    mad = median(abs.(x .- med))
    s = 1.482602218505602 * mad
    if s > eps(Float64)
        return (x .- med) ./ s
    end
    s_std = std(x)
    return s_std > eps(Float64) ? (x .- mean(x)) ./ s_std : zeros(Float64, n)
end

function _zscore_for_remaining(corr_vector::AbstractVector{<:Real}, remaining_indices::AbstractVector{Int})
    z = fill(NaN, length(corr_vector))
    if length(remaining_indices) >= 2
        z[remaining_indices] = _robust_z(corr_vector[remaining_indices])
    elseif length(remaining_indices) == 1
        z[remaining_indices] .= 0.0
    end
    return z
end

function _calculate_eog_correlations(components_t::Matrix{Float64}, eog_signal::AbstractVector{<:Real})
    n_components = size(components_t, 2)
    return [abs(cor(view(components_t, :, comp_idx), eog_signal)) for comp_idx = 1:n_components]
end

function _calculate_spatial_correlations(
    remaining_components::Vector{Int},
    primary_components::Vector{Int},
    mixing_matrix::Matrix{Float64},
    n_components::Int,
)
    spatial_corrs = fill(NaN, n_components)
    (isempty(remaining_components) || isempty(primary_components)) && return spatial_corrs

    for comp_idx in remaining_components
        comp_topo = view(mixing_matrix, :, comp_idx)
        max_corr = 0.0
        for primary_comp in primary_components
            c = abs(cor(comp_topo, view(mixing_matrix, :, primary_comp)))
            c > max_corr && (max_corr = c)
        end
        spatial_corrs[comp_idx] = max_corr
    end
    return spatial_corrs
end

function _max_lagged_corr(x_norm::AbstractVector{Float64}, y_norm::AbstractVector{Float64}, lags)
    n = length(x_norm)
    best = 0.0
    for lag in lags
        a = max(1, 1 - lag)
        b = min(n, n - lag)
        a > b && continue
        c = abs(dot(view(x_norm, a:b), view(y_norm, (a+lag):(b+lag))))
        c > best && (best = c)
    end
    return best
end

function _normalize_unit_l2(sig::AbstractVector{Float64})
    n = length(sig)
    denom = std(sig, corrected = false) * sqrt(n)
    return denom > eps(Float64) ? (sig .- mean(sig)) ./ denom : zeros(Float64, n)
end

function _calculate_lagged_correlations(
    remaining_components::Vector{Int},
    primary_components::Vector{Int},
    components_t::Matrix{Float64},
    max_lag_samples::Int,
    lag_step::Int,
)
    n_samples, n_components = size(components_t)
    lagged_corrs = fill(NaN, n_components)
    (isempty(remaining_components) || isempty(primary_components)) && return lagged_corrs

    lags = (-max_lag_samples):lag_step:max_lag_samples
    primary_norms = [_normalize_unit_l2(view(components_t, :, p_idx)) for p_idx in primary_components]

    Threads.@threads for comp_idx in remaining_components
        sig = view(components_t, :, comp_idx)
        denom = std(sig, corrected = false) * sqrt(n_samples)
        if denom > eps(Float64)
            x_norm = (sig .- mean(sig)) ./ denom
            max_c = 0.0
            for y_norm in primary_norms
                c = _max_lagged_corr(x_norm, y_norm, lags)
                c > max_c && (max_c = c)
            end
            lagged_corrs[comp_idx] = max_c
        else
            lagged_corrs[comp_idx] = 0.0
        end
    end
    return lagged_corrs
end

function _identify_secondary_components(
    remaining_components::Vector{Int},
    corr_vector::Vector{Float64},
    corr_z_vector::Vector{Float64},
    min_corr_threshold::Real,
    z_threshold::Real,
)
    return [c for c in remaining_components if corr_vector[c] > min_corr_threshold && corr_z_vector[c] > z_threshold]
end

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

function _resolve_artifact_samples(
    dat::ContinuousData,
    sample_selection,
    interval_selection::Interval,
    precomputed_samples::Union{Nothing,Vector{Int}},
)
    !isnothing(precomputed_samples) && return precomputed_samples
    return get_selected_samples(dat, _combine_interval_sample(interval_selection, sample_selection))
end

function _resolve_artifact_components(
    dat::ContinuousData,
    ica::InfoIca,
    selected_samples::Vector{Int},
    precomputed_components::Union{Nothing,Tuple{Matrix{Float64},Int}},
)
    return !isnothing(precomputed_components) ? precomputed_components : _prepare_ica_data_matrix(dat, ica, selected_samples)
end

function _evaluate_eog_channel(
    dat::ContinuousData,
    selected_samples::Vector{Int},
    eog_channel::Symbol,
    has_eog::Bool,
    components_t::Matrix{Float64},
    mixing::Matrix{Float64},
    n_components::Int,
    min_correlation::Real,
    z_threshold::Real,
    two_step::Bool,
    max_lag_samples::Int,
    lag_step::Int,
)
    corrs = fill(NaN, n_components)
    corr_z = fill(NaN, n_components)
    spatial_corr = fill(NaN, n_components)
    spatial_corr_z = fill(NaN, n_components)
    temporal_corr = fill(NaN, n_components)
    temporal_corr_z = fill(NaN, n_components)

    has_eog || return Int[], corrs, corr_z, spatial_corr, spatial_corr_z, temporal_corr, temporal_corr_z

    eog_sig = dat.data[selected_samples, eog_channel]
    corrs = _calculate_eog_correlations(components_t, eog_sig)
    corr_z = _robust_z(corrs)
    primary = findall((corr_z .> z_threshold) .& (corrs .> min_correlation))

    identified = primary
    if two_step && !isempty(primary)
        remaining = setdiff(1:n_components, primary)
        spatial_corr = _calculate_spatial_correlations(remaining, primary, mixing, n_components)
        temporal_corr = _calculate_lagged_correlations(remaining, primary, components_t, max_lag_samples, lag_step)
        spatial_corr_z = _zscore_for_remaining(spatial_corr, remaining)
        temporal_corr_z = _zscore_for_remaining(temporal_corr, remaining)

        sec_spat = _identify_secondary_components(remaining, spatial_corr, spatial_corr_z, min_correlation, z_threshold)
        sec_temp = _identify_secondary_components(remaining, temporal_corr, temporal_corr_z, min_correlation, z_threshold)
        identified = union(primary, sec_spat, sec_temp)
    end

    sort!(identified)
    return identified, corrs, corr_z, spatial_corr, spatial_corr_z, temporal_corr, temporal_corr_z
end

"""
    identify_eog_components(dat::ContinuousData, ica::InfoIca;
                            vEOG_channel=:vEOG, hEOG_channel=:hEOG,
                            z_threshold=3.0, min_correlation=0.5, two_step=true,
                            sample_selection=samples(), interval_selection=times())

Identify vertical and horizontal EOG ICA components via MAD-scaled robust z-scored correlation
(plus optional `two_step` spatial topography and lagged temporal cross-correlation matching).
Returns `(Dict(:vEOG => ..., :hEOG => ...), metrics_df)`.
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
    !has_vEOG && @minimal_warning "Vertical EOG channel :$vEOG_channel not found in data; skipping vertical EOG detection."
    !has_hEOG && @minimal_warning "Horizontal EOG channel :$hEOG_channel not found in data; skipping horizontal EOG detection."

    selected_samples = _resolve_artifact_samples(dat, sample_selection, interval_selection, _precomputed_samples)
    if isempty(selected_samples)
        @minimal_warning "No samples remaining after applying exclude criteria. Cannot identify eye components."
        return Dict(:vEOG => Int[], :hEOG => Int[]), _empty_eog_metrics()
    end

    n_components <= 10 && @info "Small number of components ($n_components <= 10); z_threshold = $z_threshold may be conservative."

    components, _ = _resolve_artifact_components(dat, ica, selected_samples, _precomputed_components)
    components_t = Matrix{Float64}(components')

    max_lag_samples = round(Int, 100.0 * dat.sample_rate / 1000.0)
    lag_step = max(1, round(Int, 5.0 * dat.sample_rate / 1000.0))

    identified_vEOG, vEOG_corrs, vEOG_corr_z, vEOG_spatial_corr, vEOG_spatial_corr_z, vEOG_temporal_corr, vEOG_temporal_corr_z =
        _evaluate_eog_channel(
            dat,
            selected_samples,
            vEOG_channel,
            has_vEOG,
            components_t,
            ica.mixing,
            n_components,
            min_correlation,
            z_threshold,
            two_step,
            max_lag_samples,
            lag_step,
        )

    identified_hEOG, hEOG_corrs, hEOG_corr_z, hEOG_spatial_corr, hEOG_spatial_corr_z, hEOG_temporal_corr, hEOG_temporal_corr_z =
        _evaluate_eog_channel(
            dat,
            selected_samples,
            hEOG_channel,
            has_hEOG,
            components_t,
            ica.mixing,
            n_components,
            min_correlation,
            z_threshold,
            two_step,
            max_lag_samples,
            lag_step,
        )

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

    return Dict{Symbol,Vector{Int}}(:vEOG => identified_vEOG, :hEOG => identified_hEOG), metrics_df
end

function _find_unidirectional_peaks(
    data::AbstractVector;
    min_prominence_std::Real = 3.0,
    window_size::Int = 1,
    min_distance::Int = 1,
    sample_coords::Union{Nothing,AbstractVector{<:Integer}} = nothing,
)
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

    if min_distance > 1 && length(peaks) > 1
        coords = isnothing(sample_coords) ? peaks : [Int(sample_coords[p]) for p in peaks]
        order = sortperm(peaks; by = p -> data[p], rev = true)
        suppressed = falses(length(peaks))
        for k in order
            suppressed[k] && continue
            c_k = coords[k]
            lo = searchsortedfirst(coords, c_k - min_distance + 1)
            hi = searchsortedlast(coords, c_k + min_distance - 1)
            for j = lo:hi
                j != k && (suppressed[j] = true)
            end
        end
        peaks = [peaks[j] for j in eachindex(peaks) if !suppressed[j]]
    end

    return peaks
end

function _calculate_ibi_metrics(peak_indices::Vector{Int}, sample_rate::Real, min_ibi_s::Real, max_ibi_s::Real)
    num_peaks = length(peak_indices)
    num_peaks < 2 && return 0, NaN, NaN, 0.0

    inv_fs = 1.0 / Float64(sample_rate)
    valid_ibis = Float64[]
    @inbounds for i = 1:(num_peaks-1)
        ibi = (peak_indices[i+1] - peak_indices[i]) * inv_fs
        (min_ibi_s <= ibi <= max_ibi_s) && push!(valid_ibis, ibi)
    end

    num_valid_ibis = length(valid_ibis)
    peak_ratio = num_valid_ibis / (num_peaks - 1)

    if num_valid_ibis > 1
        mean_ibi, std_ibi = mean(valid_ibis), std(valid_ibis)
    elseif num_valid_ibis == 1
        mean_ibi, std_ibi = valid_ibis[1], 0.0
    else
        mean_ibi, std_ibi = NaN, NaN
    end

    return num_valid_ibis, mean_ibi, std_ibi, peak_ratio
end

"""
    identify_ecg_components(dat::ContinuousData, ica::InfoIca;
                            min_bpm=40, max_bpm=120, min_prominence_std=4,
                            min_peaks=100, max_ibi_std_s=0.2, min_peak_ratio=0.5,
                            sample_selection=samples(), interval_selection=times())

Identify cardiac (ECG) ICA components via dual-polarity QRS peak detection and inter-beat interval
regularity. Returns `(ecg_indices, metrics_df)`.
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
    selected_samples = _resolve_artifact_samples(dat, sample_selection, interval_selection, _precomputed_samples)
    if isempty(selected_samples)
        @minimal_warning "No samples remaining after applying exclude criteria."
        return Int[], _empty_ecg_metrics()
    end

    components_subset, n_components = _resolve_artifact_components(dat, ica, selected_samples, _precomputed_components)
    n_components == 0 && return Int[], _empty_ecg_metrics()

    min_ibi_s = 60.0 / max_bpm
    max_ibi_s = 60.0 / min_bpm

    min_s, max_s = extrema(selected_samples)
    duration_s = (max_s - min_s + 1) / Float64(dat.sample_rate)
    expected_peaks = (duration_s * min_bpm) / 60.0
    effective_min_peaks = min(min_peaks, max(3, round(Int, expected_peaks * 0.5)))
    if effective_min_peaks < min_peaks
        @minimal_warning "Selected recording duration ($(round(duration_s, digits=1)) s) is short for default min_peaks = $min_peaks; using adapted threshold min_peaks = $effective_min_peaks."
    end

    num_peaks_col = Vector{Int}(undef, n_components)
    num_valid_col = Vector{Int}(undef, n_components)
    mean_ibi_col = Vector{Float64}(undef, n_components)
    std_ibi_col = Vector{Float64}(undef, n_components)
    peak_ratio_col = Vector{Float64}(undef, n_components)
    bpm_col = Vector{Float64}(undef, n_components)
    identified_ecg = Int[]

    n_samples = size(components_subset, 2)
    ts_zscored = Vector{Float64}(undef, n_samples)
    refractory_samples = max(1, round(Int, 0.2 * dat.sample_rate))

    for comp_idx = 1:n_components
        component_ts = @view components_subset[comp_idx, :]
        m = mean(component_ts)
        s = std(component_ts)
        inv_s = s > eps(Float64) ? (1.0 / s) : 1.0
        @inbounds @simd for j = 1:n_samples
            ts_zscored[j] = (component_ts[j] - m) * inv_s
        end

        pk_pos = _find_unidirectional_peaks(
            ts_zscored;
            min_prominence_std = min_prominence_std,
            min_distance = refractory_samples,
            sample_coords = selected_samples,
        )
        @. ts_zscored = -ts_zscored
        pk_neg = _find_unidirectional_peaks(
            ts_zscored;
            min_prominence_std = min_prominence_std,
            min_distance = refractory_samples,
            sample_coords = selected_samples,
        )

        n_valid_pos, mean_pos, std_pos, ratio_pos = _calculate_ibi_metrics(selected_samples[pk_pos], dat.sample_rate, min_ibi_s, max_ibi_s)
        n_valid_neg, mean_neg, std_neg, ratio_neg = _calculate_ibi_metrics(selected_samples[pk_neg], dat.sample_rate, min_ibi_s, max_ibi_s)

        if n_valid_pos >= n_valid_neg
            num_peaks = length(pk_pos)
            num_valid_ibis, mean_ibi, std_ibi, peak_ratio = n_valid_pos, mean_pos, std_pos, ratio_pos
        else
            num_peaks = length(pk_neg)
            num_valid_ibis, mean_ibi, std_ibi, peak_ratio = n_valid_neg, mean_neg, std_neg, ratio_neg
        end

        if num_valid_ibis >= (effective_min_peaks - 1) && !isnan(std_ibi) && std_ibi <= max_ibi_std_s && peak_ratio >= min_peak_ratio
            push!(identified_ecg, comp_idx)
        end

        num_peaks_col[comp_idx] = num_peaks
        num_valid_col[comp_idx] = num_valid_ibis
        mean_ibi_col[comp_idx] = mean_ibi
        std_ibi_col[comp_idx] = std_ibi
        peak_ratio_col[comp_idx] = peak_ratio
        bpm_col[comp_idx] = (!isnan(mean_ibi) && mean_ibi > 0) ? (60.0 / mean_ibi) : NaN
    end

    metrics_df = DataFrame(
        :Component => 1:n_components,
        :num_peaks => num_peaks_col,
        :num_valid_ibis => num_valid_col,
        :mean_ibi_s => mean_ibi_col,
        :std_ibi_s => std_ibi_col,
        :peak_ratio => peak_ratio_col,
        :heart_rate_bpm => bpm_col,
    )
    sort!(identified_ecg)
    return identified_ecg, metrics_df
end

"""
    identify_spatial_kurtosis_components(ica::InfoIca; z_threshold = 3.0)

Identify focal channel-noise ICA components via MAD-scaled robust z-scores of mixing-weight kurtosis.
"""
function identify_spatial_kurtosis_components(ica::InfoIca; z_threshold::Real = 3.0)
    n_components = size(ica.mixing, 2)
    spatial_kurtosis = [kurtosis(view(ica.mixing, :, i)) for i = 1:n_components]
    z_spatial_kurtosis = _robust_z(spatial_kurtosis)
    high_kurtosis_comps = sort!(findall(z_spatial_kurtosis .> z_threshold))

    metrics_df = DataFrame(:Component => 1:n_components, :spatial_kurtosis => spatial_kurtosis, :z_spatial_kurtosis => z_spatial_kurtosis)
    return high_kurtosis_comps, metrics_df
end

"""
    identify_line_noise_components(dat::ContinuousData, ica::InfoIca;
                                   line_freq=50.0, freq_bandwidth=1.0,
                                   z_threshold=3.0, min_harmonic_power=1.5,
                                   sample_selection=samples(), interval_selection=times())

Identify line-noise ICA components via Welch PSD peak-to-flank ratios at `line_freq` (and its 2nd harmonic
when below Nyquist). Returns `(line_noise_indices, metrics_df)`.
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

    samples_to_use = _resolve_artifact_samples(dat, sample_selection, interval_selection, _precomputed_samples)
    if isempty(samples_to_use)
        @minimal_warning "No samples remaining after applying exclude criteria. Cannot identify line noise components."
        return Int[], _empty_line_noise_metrics()
    end

    components, n_components = _resolve_artifact_components(dat, ica, samples_to_use, _precomputed_components)
    n_samples = size(components, 2)
    (n_components == 0 || n_samples == 0) && return Int[], _empty_line_noise_metrics()

    n_per_seg = min(n_samples, nextpow(2, max(64, round(Int, 2.0 * dat.sample_rate))))
    n_overlap = div(n_per_seg, 2)
    wcfg = DSP.Periodograms.WelchConfig(n_samples, Float64; n = n_per_seg, noverlap = n_overlap, fs = dat.sample_rate, window = DSP.hanning)

    sig_buf = Vector{Float64}(undef, n_samples)
    @inbounds @simd for s = 1:n_samples
        sig_buf[s] = components[1, s]
    end
    first_pgram = DSP.welch_pgram(sig_buf, wcfg)
    freqs = DSP.freq(first_pgram)
    n_freqs = length(freqs)

    psd = Matrix{Float64}(undef, n_freqs, n_components)
    psd[:, 1] .= DSP.power(first_pgram)
    for i = 2:n_components
        @inbounds @simd for s = 1:n_samples
            sig_buf[s] = components[i, s]
        end
        psd[:, i] .= DSP.power(DSP.welch_pgram(sig_buf, wcfg))
    end

    flank_bandwidth = 5.0
    line_band = findall(f -> abs(f - line_freq) <= freq_bandwidth, freqs)
    line_flank_band = findall(f -> freq_bandwidth < abs(f - line_freq) <= freq_bandwidth + flank_bandwidth, freqs)

    harmonic2_freq = line_freq * 2.0
    harmonic2_band = findall(f -> abs(f - harmonic2_freq) <= freq_bandwidth, freqs)
    harmonic2_flank_band = findall(f -> freq_bandwidth < abs(f - harmonic2_freq) <= freq_bandwidth + flank_bandwidth, freqs)

    line_power_col = Vector{Float64}(undef, n_components)
    flank_power_col = Vector{Float64}(undef, n_components)
    power_ratio_col = Vector{Float64}(undef, n_components)
    harmonic_ratio_col = Vector{Float64}(undef, n_components)

    for i = 1:n_components
        lp = !isempty(line_band) ? mean(view(psd, line_band, i)) : 0.0
        fp = !isempty(line_flank_band) ? mean(view(psd, line_flank_band, i)) : mean(view(psd, :, i))
        pr = lp / (fp + eps())

        hr = if !isempty(harmonic2_band)
            h2_power = mean(view(psd, harmonic2_band, i))
            h2_flank = !isempty(harmonic2_flank_band) ? mean(view(psd, harmonic2_flank_band, i)) : fp
            h2_power / (h2_flank + eps())
        else
            NaN
        end

        line_power_col[i] = lp
        flank_power_col[i] = fp
        power_ratio_col[i] = pr
        harmonic_ratio_col[i] = hr
    end

    power_ratio_z = _robust_z(power_ratio_col)
    metrics_df = DataFrame(
        :Component => 1:n_components,
        :line_power => line_power_col,
        :flanking_power => flank_power_col,
        :power_ratio => power_ratio_col,
        :harmonic_ratio => harmonic_ratio_col,
        :power_ratio_zscore => power_ratio_z,
    )

    line_noise_comps = sort!(
        findall(
            i -> power_ratio_z[i] > z_threshold && (isnan(harmonic_ratio_col[i]) || harmonic_ratio_col[i] > min_harmonic_power),
            1:n_components,
        ),
    )
    return line_noise_comps, metrics_df
end

"""
    ArtifactComponents

Container of identified ICA artifact component indices (`eog`, `ecg`, `line_noise`, `channel_noise`, `manual`).
"""
struct ArtifactComponents
    eog::Dict{Symbol,Vector{Int}}
    ecg::Vector{Int}
    line_noise::Vector{Int}
    channel_noise::Vector{Int}
    manual::Vector{Int}
end

ArtifactComponents(eog::Dict{Symbol,Vector{Int}}, ecg::Vector{Int}, line_noise::Vector{Int}, channel_noise::Vector{Int}) =
    ArtifactComponents(eog, ecg, line_noise, channel_noise, Int[])

"""
    combine_artifact_components(eog_comps, ecg_comps, line_noise_comps, channel_noise_comps) -> ArtifactComponents

Combine per-category ICA artifact indices into an [`ArtifactComponents`](@ref) struct.
"""
combine_artifact_components(
    eog_comps::Dict{Symbol,Vector{Int}},
    ecg_comps::Vector{Int},
    line_noise_comps::Vector{Int},
    channel_noise_comps::Vector{Int},
) = ArtifactComponents(eog_comps, ecg_comps, line_noise_comps, channel_noise_comps)

"""
    get_all_ica_components(artifacts::ArtifactComponents) -> Vector{Int}

Return all unique flagged ICA component indices in `artifacts` as a sorted vector.
"""
function get_all_ica_components(artifacts::ArtifactComponents)
    all_comps = Set{Int}()
    for comps in values(artifacts.eog)
        union!(all_comps, comps)
    end
    union!(all_comps, artifacts.ecg, artifacts.line_noise, artifacts.channel_noise, artifacts.manual)
    return sort!(collect(all_comps))
end

"""
    get_all_components(artifacts::ArtifactComponents) -> Vector{Int}

Alias for [`get_all_ica_components`](@ref).
"""
get_all_components(artifacts::ArtifactComponents) = get_all_ica_components(artifacts)

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
    !isempty(artifacts.manual) && println(io, "  Manual: $(artifacts.manual)")
    print(io, "  All: $all_comps")
end

"""
    identify_components(dat::ContinuousData, ica::InfoIca; method = :correlation,
                        sample_selection = samples(), interval_selection = times(),
                        eog = (;), ecg = (;), line_noise = (;), channel_noise = (;),
                        kwargs...) -> (ArtifactComponents, Dict{Symbol, DataFrame})

Identify EOG, ECG, line-noise, and focal channel-noise ICA components in a single pass,
sharing precomputed ICA activations across detectors.
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
    method === :correlation || @minimal_error("Unknown component identification method: :$method. Supported methods: :correlation")

    selected_samples = get_selected_samples(dat, _combine_interval_sample(interval_selection, sample_selection))
    precomputed = _prepare_ica_data_matrix(dat, ica, selected_samples)

    eog_keys = (:vEOG_channel, :hEOG_channel, :z_threshold, :min_correlation, :two_step)
    ecg_keys = (:min_bpm, :max_bpm, :min_prominence_std, :min_peaks, :max_ibi_std_s, :min_peak_ratio)
    line_noise_keys = (:line_freq, :freq_bandwidth, :z_threshold, :min_harmonic_power)
    channel_noise_keys = (:z_threshold,)

    all_known = (eog_keys..., ecg_keys..., line_noise_keys..., channel_noise_keys...)
    unknown = [k for k in keys(kwargs) if k ∉ all_known]
    !isempty(unknown) && @minimal_warning "Unrecognized keyword arguments passed to identify_components: $(join(unknown, ", "))."

    eog_merged = merge(NamedTuple(k => v for (k, v) in pairs(kwargs) if k in eog_keys), eog)
    ecg_merged = merge(NamedTuple(k => v for (k, v) in pairs(kwargs) if k in ecg_keys), ecg)
    line_noise_merged = merge(NamedTuple(k => v for (k, v) in pairs(kwargs) if k in line_noise_keys), line_noise)
    channel_noise_merged = merge(NamedTuple(k => v for (k, v) in pairs(kwargs) if k in channel_noise_keys), channel_noise)

    eog_comps, eog_metrics_df = identify_eog_components(
        dat,
        ica;
        sample_selection = sample_selection,
        interval_selection = interval_selection,
        _precomputed_components = precomputed,
        _precomputed_samples = selected_samples,
        eog_merged...,
    )
    ecg_comps, ecg_metrics_df = identify_ecg_components(
        dat,
        ica;
        sample_selection = sample_selection,
        interval_selection = interval_selection,
        _precomputed_components = precomputed,
        _precomputed_samples = selected_samples,
        ecg_merged...,
    )
    line_noise_comps, line_noise_metrics_df = identify_line_noise_components(
        dat,
        ica;
        sample_selection = sample_selection,
        interval_selection = interval_selection,
        _precomputed_components = precomputed,
        _precomputed_samples = selected_samples,
        line_noise_merged...,
    )
    channel_noise_comps, channel_noise_metrics_df = identify_spatial_kurtosis_components(ica; channel_noise_merged...)

    artifacts = combine_artifact_components(eog_comps, ecg_comps, line_noise_comps, channel_noise_comps)
    metrics = Dict{Symbol,DataFrame}(
        :eog_metrics => eog_metrics_df,
        :ecg_metrics => ecg_metrics_df,
        :line_noise_metrics => line_noise_metrics_df,
        :channel_noise_metrics => channel_noise_metrics_df,
    )
    return artifacts, metrics
end
