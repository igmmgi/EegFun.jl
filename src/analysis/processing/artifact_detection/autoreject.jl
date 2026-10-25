# === AUTOREJECT: ADAPTIVE ARTIFACT REJECTION AND REPAIR ===
#
# Implementation of the Autoreject algorithm based on:
# Jas, M., Engemann, D. A., Bekhti, Y., Raimondo, F., & Gramfort, A. (2017).
# "Autoreject: Automated artifact rejection for MEG and EEG data".
# NeuroImage, 159, 417-429. doi:10.1016/j.neuroimage.2017.06.030

"""
    AutoRejectResult

Result of adaptive artifact rejection and channel repair (Autoreject approach).

# Fields
- `thresholds::Dict{Symbol,Float64}`: Learned peak-to-peak amplitude threshold per channel (μV)
- `mode::Symbol`: Rejection mode (`:local` or `:global`)
- `consensus_perc::Float64`: Consensus rejection threshold κ (fraction of bad channels before epoch is rejected)
- `max_repair_perc::Float64`: Maximum repair threshold ρ (fraction of bad channels allowed to be repaired)
- `channels::Vector{Symbol}`: Analyzed channels
- `ptp::Matrix{Float64}`: Peak-to-peak amplitude matrix (n_epochs × n_channels)
- `bad_matrix::BitMatrix`: Binary matrix where `true` indicates a channel exceeded its threshold in an epoch
- `epochs_clean::Vector{Int}`: Global indices of clean epochs (zero bad channels)
- `epochs_repaired::Vector{Int}`: Global indices of epochs that had channels repaired
- `epochs_rejected::Vector{Int}`: Global indices of epochs rejected
- `repaired_channels::OrderedDict{Int,Vector{Symbol}}`: Mapping of global epoch index => channels repaired
- `rejection_info::Vector{EpochRejectionInfo}`: Standard EpochRejectionInfo per condition for pipeline compatibility
- `summary_table::DataFrame`: Summary DataFrame of repair/rejection metrics
- `channel_table::DataFrame`: Per-channel threshold and artifact occurrence table
"""
struct AutoRejectResult
    thresholds::Dict{Symbol,Float64}
    mode::Symbol
    consensus_perc::Float64
    max_repair_perc::Float64
    channels::Vector{Symbol}
    ptp::Matrix{Float64}
    bad_matrix::BitMatrix
    epochs_clean::Vector{Int}
    epochs_repaired::Vector{Int}
    epochs_rejected::Vector{Int}
    repaired_channels::OrderedDict{Int,Vector{Symbol}}
    rejection_info::Vector{EpochRejectionInfo}
    summary_table::DataFrame
    channel_table::DataFrame
end

function Base.show(io::IO, res::AutoRejectResult)
    n_tot = length(res.epochs_clean) + length(res.epochs_repaired) + length(res.epochs_rejected)
    n_kept = length(res.epochs_clean) + length(res.epochs_repaired)
    retention = n_tot > 0 ? round(100 * n_kept / n_tot, digits = 1) : 0.0
    th_vals = collect(values(res.thresholds))
    min_th = isempty(th_vals) ? 0.0 : minimum(th_vals)
    med_th = isempty(th_vals) ? 0.0 : median(th_vals)
    max_th = isempty(th_vals) ? 0.0 : maximum(th_vals)

    print(io, "AutoRejectResult(")
    print(io, "mode=:$(res.mode), ")
    print(io, "epochs=$(n_kept)/$(n_tot) kept ($(retention)%), ")
    print(io, "repaired=$(length(res.epochs_repaired)), ")
    print(io, "rejected=$(length(res.epochs_rejected)), ")
    print(io, "thresholds=[min=$(round(min_th, digits=1)), med=$(round(med_th, digits=1)), max=$(round(max_th, digits=1))] μV")
    print(io, ")")
end

# === CORE AUTOREJECT COMPUTATION ===

"""
    _kfold_splits(n_total::Int, n_folds::Int, seed::Union{Int,Nothing}) -> Vector{Vector{Int}}

Generate k-fold cross-validation partition of indices `1:n_total`.
"""
function _kfold_splits(n_total::Int, n_folds::Int, seed::Union{Int,Nothing})
    rng = isnothing(seed) ? Random.default_rng() : Random.MersenneTwister(seed)
    perm = randperm(rng, n_total)
    folds = Vector{Vector{Int}}(undef, n_folds)
    fold_size = div(n_total, n_folds)
    rem_samples = rem(n_total, n_folds)

    start_idx = 1
    for f = 1:n_folds
        cur_size = fold_size + (f <= rem_samples ? 1 : 0)
        end_idx = start_idx + cur_size - 1
        folds[f] = perm[start_idx:end_idx]
        start_idx = end_idx + 1
    end
    return folds
end

"""
    _candidate_thresholds(ptp_values::AbstractVector{<:Real}, n_candidates::Int, quant_range::Tuple{Float64,Float64}) -> Vector{Float64}

Generate sorted, unique candidate thresholds from quantiles of peak-to-peak amplitudes.
"""
function _candidate_thresholds(ptp_values::AbstractVector{<:Real}, n_candidates::Int, quant_range::Tuple{Float64,Float64} = (0.4, 0.99))
    isempty(ptp_values) && return [100.0]

    min_p, max_p = extrema(ptp_values)
    if min_p ≈ max_p
        return [max_p > 0 ? Float64(max_p * 1.05) : 100.0]
    end

    q_levels = range(quant_range[1], quant_range[2], length = max(5, n_candidates))
    candidates = Float64.(unique(quantile(ptp_values, q_levels)))
    filter!(x -> x > 0 && isfinite(x), candidates)

    if isempty(candidates)
        candidates = [Float64(max(100.0, max_p))]
    end

    sort!(candidates)
    return candidates
end

"""
    _cv_channel_threshold(ch_matrix::Matrix{Float64}, ptp_col::Vector{Float64}, candidates::Vector{Float64}, folds::Vector{Vector{Int}}) -> Float64

Estimate optimal threshold for a single channel via k-fold cross-validation.
`ch_matrix` has dimensions `(n_timepoints, n_epochs)`.
"""
function _cv_channel_threshold(
    ch_matrix::Matrix{Float64},
    ptp_col::Vector{Float64},
    candidates::Vector{Float64},
    folds::Vector{Vector{Int}},
)
    n_timepoints, n_epochs = size(ch_matrix)
    n_folds = length(folds)
    n_candidates = length(candidates)

    if n_candidates == 1
        return candidates[1]
    end

    total_errors = zeros(Float64, n_candidates)

    for fold_idx = 1:n_folds
        test_indices = folds[fold_idx]
        isempty(test_indices) && continue

        train_indices = setdiff(1:n_epochs, test_indices)
        isempty(train_indices) && continue

        # Compute robust median signal of the test fold at each time point
        test_median = zeros(Float64, n_timepoints)
        for t = 1:n_timepoints
            test_median[t] = median(@view ch_matrix[t, test_indices])
        end

        # Pre-filter training epochs by threshold candidates
        for (c_idx, th) in enumerate(candidates)
            good_train = [i for i in train_indices if ptp_col[i] <= th]

            if isempty(good_train)
                total_errors[c_idx] += Inf
                continue
            end

            # Mean training signal across good training epochs
            err_sq = 0.0
            n_good = length(good_train)
            for t = 1:n_timepoints
                s_sum = 0.0
                for ep_i in good_train
                    s_sum += ch_matrix[t, ep_i]
                end
                train_mean_t = s_sum / n_good
                diff_t = train_mean_t - test_median[t]
                err_sq += diff_t * diff_t
            end

            total_errors[c_idx] += sqrt(err_sq / n_timepoints)
        end
    end

    # Return candidate with minimum cross-validated error
    best_idx = argmin(total_errors)
    return candidates[best_idx]
end

"""
    _repair_epochs_cached!(epochs_data::Vector{DataFrame}, epochs_to_repair::Vector{Int},
                           bad_matrix::BitMatrix, channels::Vector{Symbol}, layout::Layout;
                           repair_method::Symbol = :spherical_spline, m::Int = 4, lambda::Float64 = 1e-5)
                           -> OrderedDict{Int, Vector{Symbol}}

Interpolate bad channels in epochs marked for repair, caching the interpolation operator
for epochs that share the identical pattern of bad channels.
"""
function _repair_epochs_cached!(
    epochs_data::Vector{DataFrame},
    epochs_to_repair::Vector{Int},
    bad_matrix::BitMatrix,
    channels::Vector{Symbol},
    layout::Layout;
    repair_method::Symbol = :spherical_spline,
    m::Int = 4,
    lambda::Float64 = 1e-5,
)
    repaired_channels = OrderedDict{Int,Vector{Symbol}}()
    isempty(epochs_to_repair) && return repaired_channels

    # Determine 3D coordinates availability for spherical splines
    has_spline_coords = false
    if repair_method == :spherical_spline
        _ensure_coordinates_3d!(layout)
        has_spline_coords = has_3d_coords(layout)
    end

    # Channel symbol to index in layout/channels
    ch_indices = Dict(ch => i for (i, ch) in enumerate(channels))
    n_channels = length(channels)

    # Cache for spherical spline interpolation operators: pattern => (good_indices, bad_indices, operator_matrix)
    spline_operator_cache = Dict{Vector{Symbol},Tuple{Vector{Int},Vector{Int},Matrix{Float64}}}()

    # Pre-extract coordinates for all channels if using spherical spline
    coords = zeros(Float64, n_channels, 3)
    if has_spline_coords
        for (i, ch) in enumerate(channels)
            idx = findfirst(x -> x == ch, layout.data.label)
            if !isnothing(idx) && hasproperty(layout.data, :x3)
                coords[i, 1] = Float64(layout.data.x3[idx])
                coords[i, 2] = Float64(layout.data.y3[idx])
                coords[i, 3] = Float64(layout.data.z3[idx])
            end
        end
    end

    for ep_idx in epochs_to_repair
        bad_mask = @view bad_matrix[ep_idx, :]
        bad_chs = channels[bad_mask]
        isempty(bad_chs) && continue

        # Sort pattern to allow reuse across epochs with identical bad channels
        pattern = sort(bad_chs)
        epoch_df = epochs_data[ep_idx]

        if repair_method == :spherical_spline && has_spline_coords
            bad_idx = [ch_indices[ch] for ch in pattern if haskey(ch_indices, ch)]
            good_idx = setdiff(1:n_channels, bad_idx)

            if length(good_idx) >= 4 && !isempty(bad_idx)
                if !haskey(spline_operator_cache, pattern)
                    pos_good = coords[good_idx, :]
                    pos_bad = coords[bad_idx, :]
                    interp_mat = _make_interpolation_matrix(pos_good, pos_bad; m = m, lambda = lambda)
                    spline_operator_cache[pattern] = (good_idx, bad_idx, interp_mat)
                end

                g_idx, b_idx, op_mat = spline_operator_cache[pattern]

                # Extract data matrix (channels × timepoints)
                epoch_mat = Matrix(epoch_df[:, channels])'
                # Apply interpolation operator
                epoch_mat[b_idx, :] = op_mat * epoch_mat[g_idx, :]
                # Store back
                epoch_df[:, channels] = epoch_mat'

                repaired_channels[ep_idx] = pattern
                continue
            end
        end

        # Fallback to neighbor interpolation if spline coordinates not available or failed
        if !isnothing(layout.neighbours)
            epoch_mat = Matrix(epoch_df[:, channels])'
            _repair_channels_neighbor!(epoch_mat, channels, pattern, layout.neighbours)
            epoch_df[:, channels] = epoch_mat'
            repaired_channels[ep_idx] = pattern
        end
    end

    return repaired_channels
end

# === PUBLIC AUTOREJECT API ===

"""
    autoreject!(dat::Vector{EpochData}; kwargs...) -> AutoRejectResult
    autoreject!(dat::EpochData; kwargs...) -> AutoRejectResult

Adaptive artifact rejection and channel repair (Autoreject approach) mutating `dat` in-place.

When `reject = true` (default), rejected epochs are removed from `dat`. When `reject = false`,
repairable bad channels are interpolated in-place but bad epochs are not pruned (useful when
further interactive review is desired).

# Arguments
- `dat`: `EpochData` or `Vector{EpochData}` to clean.

# Keyword Arguments
- `mode::Symbol`: Rejection mode (`:local` or `:global`, default: `:local`).
  - `:local`: Learns per-channel thresholds via cross-validation, repairs epochs with few bad channels,
    and rejects epochs with excessive bad channels.
  - `:global`: Learns a single global threshold across all channels, reject-only.
- `channel_selection::Function`: Channel selector (default: `channels()`).
- `interval_selection::Interval`: Time window for peak-to-peak evaluation (default: `times()`).
- `n_folds::Int`: Number of cross-validation folds (default: 10).
- `n_candidates::Int`: Number of candidate thresholds to evaluate per channel (default: 30).
- `candidate_quantiles::Tuple{Float64,Float64}`: Quantile range for candidate thresholds (default: `(0.4, 0.99)`).
- `consensus_perc::Float64`: Maximum fraction of bad channels in an epoch before the entire epoch is rejected (κ, default: 0.2).
- `max_repair_perc::Float64`: Maximum fraction of bad channels in an epoch that can be repaired (ρ, default: 0.1).
- `repair_method::Symbol`: Repair method (`:spherical_spline` or `:neighbor_interpolation`, default: `:spherical_spline`).
- `m::Int`: Legendre polynomial order for spherical splines (default: 4).
- `lambda::Float64`: Regularization parameter for spherical splines (default: 1e-5).
- `seed::Union{Int,Nothing}`: Random seed for CV splits (default: 42).
- `reject::Bool`: Whether to remove rejected epochs from `dat` (default: true).
- `verbose::Bool`: Whether to log progress (default: true).

# Returns
- `AutoRejectResult`: Complete result with learned thresholds, repair/rejection counts, and summary tables.

# References
Jas, M., Engemann, D. A., Bekhti, Y., Raimondo, F., & Gramfort, A. (2017). Autoreject: Automated artifact rejection
for MEG and EEG data. NeuroImage, 159, 417-429.
"""
function autoreject!(
    dat::Vector{EpochData};
    mode::Symbol = :local,
    channel_selection::Function = channels(),
    interval_selection::Interval = times(),
    n_folds::Int = 10,
    n_candidates::Int = 30,
    candidate_quantiles::Tuple{Float64,Float64} = (0.4, 0.99),
    consensus_perc::Float64 = 0.2,
    max_repair_perc::Float64 = 0.1,
    repair_method::Symbol = :spherical_spline,
    m::Int = 4,
    lambda::Float64 = 1e-5,
    seed::Union{Int,Nothing} = 42,
    reject::Bool = true,
    verbose::Bool = true,
)
    isempty(dat) && @minimal_error("Cannot run autoreject on empty epoch vector")

    mode ∉ [:local, :global] && @minimal_error("Unknown mode: $mode. Use :local or :global.")
    (consensus_perc < 0 || consensus_perc > 1) && @minimal_error("consensus_perc must be between 0.0 and 1.0")
    (max_repair_perc < 0 || max_repair_perc > 1) && @minimal_error("max_repair_perc must be between 0.0 and 1.0")

    # Layout from first condition
    layout = dat[1].layout
    selected_channels = get_selected_channels(dat[1], channel_selection; include_meta = false, include_extra = false)
    isempty(selected_channels) && @minimal_error("No channels selected for autoreject")

    # Flatten all epochs across conditions to pool thresholds (prevents condition-specific bias)
    all_epoch_dfs = DataFrame[]
    epoch_mapping = Tuple{Int,Int}[] # (condition_idx, local_epoch_idx)

    for (cond_idx, ep_dat) in enumerate(dat)
        for loc_idx = 1:n_epochs(ep_dat)
            push!(all_epoch_dfs, ep_dat.data[loc_idx])
            push!(epoch_mapping, (cond_idx, loc_idx))
        end
    end

    n_total_epochs = length(all_epoch_dfs)
    n_channels = length(selected_channels)

    if n_total_epochs == 0
        @minimal_warning("No epochs found to run autoreject.")
        empty_res = AutoRejectResult(
            Dict{Symbol,Float64}(),
            mode,
            consensus_perc,
            max_repair_perc,
            selected_channels,
            zeros(Float64, 0, n_channels),
            falses(0, n_channels),
            Int[],
            Int[],
            Int[],
            OrderedDict{Int,Vector{Symbol}}(),
            EpochRejectionInfo[],
            DataFrame(),
            DataFrame(),
        )
        return empty_res
    end

    verbose && @info "Running Autoreject (mode = :$mode) on $n_total_epochs epochs across $n_channels channels..."

    # Determine sample window
    first_epoch = all_epoch_dfs[1]
    combined_sel = _combine_interval_sample(interval_selection, samples())
    sample_indices = get_selected_samples(first_epoch, combined_sel)
    isempty(sample_indices) && @minimal_error("Selected interval contains no samples")
    n_timepoints = length(sample_indices)

    # 1. Compute Peak-to-Peak (PTP) amplitude matrix: (n_total_epochs × n_channels)
    ptp_matrix = zeros(Float64, n_total_epochs, n_channels)
    for (ep_idx, ep_df) in enumerate(all_epoch_dfs)
        for (ch_idx, ch) in enumerate(selected_channels)
            col = ep_df[!, ch]
            s_view = @view col[sample_indices]
            ptp_matrix[ep_idx, ch_idx] = Float64(maximum(s_view) - minimum(s_view))
        end
    end

    # 2. Threshold Estimation via Cross-Validation
    thresholds = Dict{Symbol,Float64}()

    # Determine number of CV folds based on dataset size
    eff_folds = n_folds
    if n_total_epochs < 4
        @minimal_warning("Very few epochs ($n_total_epochs); falling back to 90th percentile threshold without cross-validation.")
        eff_folds = 1
    elseif n_total_epochs < eff_folds
        eff_folds = max(2, div(n_total_epochs, 2))
        verbose && @info "Adjusted cross-validation folds to $eff_folds due to trial count ($n_total_epochs)"
    end

    folds = eff_folds > 1 ? _kfold_splits(n_total_epochs, eff_folds, seed) : [collect(1:n_total_epochs)]

    if mode == :global
        # Global mode: single threshold across all channels
        all_ptps = vec(ptp_matrix)
        global_candidates = _candidate_thresholds(all_ptps, n_candidates, candidate_quantiles)

        if eff_folds > 1
            # Evaluate across all channels simultaneously
            candidate_errors = zeros(Float64, length(global_candidates))
            for f_indices in folds
                train_idx = setdiff(1:n_total_epochs, f_indices)
                (isempty(f_indices) || isempty(train_idx)) && continue

                for (c_idx, th) in enumerate(global_candidates)
                    total_ch_err = 0.0
                    for (ch_i, ch) in enumerate(selected_channels)
                        # Extract test median
                        ch_mat = zeros(Float64, n_timepoints, n_total_epochs)
                        for (ep_i, ep_df) in enumerate(all_epoch_dfs)
                            ch_mat[:, ep_i] = ep_df[sample_indices, ch]
                        end
                        test_med = [median(@view ch_mat[t, f_indices]) for t = 1:n_timepoints]

                        good_train = [i for i in train_idx if ptp_matrix[i, ch_i] <= th]
                        if isempty(good_train)
                            total_ch_err += Inf
                        else
                            diff_norm = 0.0
                            n_g = length(good_train)
                            for t = 1:n_timepoints
                                tr_mean = sum(@view ch_mat[t, good_train]) / n_g
                                diff_norm += (tr_mean - test_med[t])^2
                            end
                            total_ch_err += sqrt(diff_norm / n_timepoints)
                        end
                    end
                    candidate_errors[c_idx] += total_ch_err
                end
            end
            best_global_th = global_candidates[argmin(candidate_errors)]
        else
            best_global_th = quantile(all_ptps, 0.9)
        end

        for ch in selected_channels
            thresholds[ch] = best_global_th
        end

    else
        # Local mode: optimize threshold per channel
        # Pre-extract channel data arrays to maximize loop efficiency
        ch_arrays = [zeros(Float64, n_timepoints, n_total_epochs) for _ = 1:n_channels]
        for (ep_i, ep_df) in enumerate(all_epoch_dfs)
            for (ch_i, ch) in enumerate(selected_channels)
                ch_arrays[ch_i][:, ep_i] = ep_df[sample_indices, ch]
            end
        end

        threshold_results = Vector{Float64}(undef, n_channels)

        Threads.@threads for ch_idx = 1:n_channels
            ptp_col = ptp_matrix[:, ch_idx]
            candidates = _candidate_thresholds(ptp_col, n_candidates, candidate_quantiles)

            if eff_folds > 1
                best_th = _cv_channel_threshold(ch_arrays[ch_idx], ptp_col, candidates, folds)
            else
                best_th = quantile(ptp_col, 0.9)
            end
            threshold_results[ch_idx] = best_th
        end

        for (ch_idx, ch) in enumerate(selected_channels)
            thresholds[ch] = threshold_results[ch_idx]
        end
    end

    # 3. Decision Matrix: Bad Channels per Epoch
    bad_matrix = falses(n_total_epochs, n_channels)
    for ep_idx = 1:n_total_epochs
        for ch_idx = 1:n_channels
            ch = selected_channels[ch_idx]
            if ptp_matrix[ep_idx, ch_idx] > thresholds[ch]
                bad_matrix[ep_idx, ch_idx] = true
            end
        end
    end

    # 4. Classify Epochs: Clean, Repairable, Rejected
    epochs_clean = Int[]
    epochs_repaired = Int[]
    epochs_rejected = Int[]

    # In global mode: pure rejection (no repair)
    effective_max_repair = mode == :global ? 0.0 : max_repair_perc
    k_max = round(Int, consensus_perc * n_channels)
    r_max = round(Int, effective_max_repair * n_channels)
    r_max = min(r_max, k_max)

    for ep_idx = 1:n_total_epochs
        n_bad = sum(@view bad_matrix[ep_idx, :])
        if n_bad == 0
            push!(epochs_clean, ep_idx)
        elseif n_bad <= r_max && mode != :global
            push!(epochs_repaired, ep_idx)
        else
            push!(epochs_rejected, ep_idx)
        end
    end

    # 5. Repair Channels for Repairable Epochs
    repaired_channels = OrderedDict{Int,Vector{Symbol}}()
    if !isempty(epochs_repaired) && mode != :global
        verbose && @info "Repairing bad channels across $(length(epochs_repaired)) epochs using $repair_method..."
        repaired_channels = _repair_epochs_cached!(
            all_epoch_dfs,
            epochs_repaired,
            bad_matrix,
            selected_channels,
            layout;
            repair_method = repair_method,
            m = m,
            lambda = lambda,
        )
    end

    # 6. Map back to conditions and construct EpochRejectionInfo
    rejection_info_per_condition = EpochRejectionInfo[]

    for (c_idx, ep_dat) in enumerate(dat)
        local_total = n_epochs(ep_dat)
        cond_num = ep_dat.condition
        cond_name = ep_dat.condition_name

        c_rejections = Rejection[]
        c_repaired = OrderedDict{Int,Vector{Symbol}}()
        c_skipped = OrderedDict{Int,Vector{Symbol}}()

        for (g_idx, (m_c_idx, m_loc_idx)) in enumerate(epoch_mapping)
            m_c_idx != c_idx && continue

            # If rejected
            if g_idx in epochs_rejected
                for ch_i in findall(@view bad_matrix[g_idx, :])
                    push!(c_rejections, Rejection(selected_channels[ch_i], m_loc_idx))
                end
                # In case no channel exceeded but was marked rejected
                if isempty(c_rejections) || last(c_rejections).epoch != m_loc_idx
                    push!(c_rejections, Rejection(selected_channels[1], m_loc_idx))
                end
            end

            # If repaired
            if haskey(repaired_channels, g_idx)
                c_repaired[m_loc_idx] = repaired_channels[g_idx]
            end
        end

        info = EpochInfo(cond_num, cond_name, local_total)
        median_th = median(collect(values(thresholds)))

        rej_info = EpochRejectionInfo(
            "autoreject",
            info,
            length(c_rejections),
            median_th,
            c_rejections,
            0.0,
            nothing,
            unique_rejections(c_rejections),
            c_repaired,
            c_skipped,
        )
        push!(rejection_info_per_condition, rej_info)
    end

    # 7. Apply epoch rejection if requested
    if reject
        for (c_idx, ep_dat) in enumerate(dat)
            r_info = rejection_info_per_condition[c_idx]
            reject_epochs!(ep_dat, r_info)
        end
    end

    # 8. Create Summary and Channel Reporting Tables
    n_kept = length(epochs_clean) + length(epochs_repaired)
    retention_pct = n_total_epochs > 0 ? round(100 * n_kept / n_total_epochs, digits = 1) : 0.0

    all_repaired_chs = Symbol[]
    for ch_list in values(repaired_channels)
        append!(all_repaired_chs, ch_list)
    end
    unique_repaired_count = length(unique(all_repaired_chs))

    th_vals = [thresholds[ch] for ch in selected_channels]

    summary_df = DataFrame(
        Metric = [
            "Total Epochs",
            "Clean Epochs (no artifacts)",
            "Repaired Epochs",
            "Rejected Epochs",
            "Retained Epochs (clean + repaired)",
            "Retention Rate (%)",
            "Unique Channels Repaired",
            "Mean Threshold (μV)",
            "Min Threshold (μV)",
            "Max Threshold (μV)",
            "Consensus Threshold (κ)",
            "Max Repair Threshold (ρ)",
        ],
        Value = [
            string(n_total_epochs),
            string(length(epochs_clean)),
            string(length(epochs_repaired)),
            string(length(epochs_rejected)),
            string(n_kept),
            Printf.@sprintf("%.1f%%", retention_pct),
            string(unique_repaired_count),
            Printf.@sprintf("%.2f", mean(th_vals)),
            Printf.@sprintf("%.2f", minimum(th_vals)),
            Printf.@sprintf("%.2f", maximum(th_vals)),
            Printf.@sprintf("%.2f (%d chs)", consensus_perc, k_max),
            Printf.@sprintf("%.2f (%d chs)", effective_max_repair, r_max),
        ],
    )

    channel_df = DataFrame(
        channel = selected_channels,
        threshold_uV = [round(thresholds[ch], digits = 2) for ch in selected_channels],
        bad_epochs_count = [sum(bad_matrix[:, c]) for c = 1:n_channels],
        bad_epochs_pct = [round(100 * sum(bad_matrix[:, c]) / max(1, n_total_epochs), digits = 1) for c = 1:n_channels],
    )

    result = AutoRejectResult(
        thresholds,
        mode,
        consensus_perc,
        effective_max_repair,
        selected_channels,
        ptp_matrix,
        bad_matrix,
        epochs_clean,
        epochs_repaired,
        epochs_rejected,
        repaired_channels,
        rejection_info_per_condition,
        summary_df,
        channel_df,
    )

    verbose && @info "Autoreject completed: $n_kept/$n_total_epochs epochs retained ($retention_pct%)."
    return result
end

function autoreject!(dat::EpochData; kwargs...)
    vec_res = autoreject!([dat]; kwargs...)
    return vec_res
end

"""
    autoreject(dat::Vector{EpochData}; kwargs...) -> Tuple{Vector{EpochData}, AutoRejectResult}
    autoreject(dat::EpochData; kwargs...) -> Tuple{EpochData, AutoRejectResult}

Non-mutating version of `autoreject!`. Returns a deep copy of `dat` and the `AutoRejectResult`.
"""
function autoreject(dat::Vector{EpochData}; kwargs...)
    dat_copy = [copy(d) for d in dat]
    res = autoreject!(dat_copy; kwargs...)
    return dat_copy, res
end

function autoreject(dat::EpochData; kwargs...)
    dat_copy = copy(dat)
    res = autoreject!(dat_copy; kwargs...)
    return dat_copy, res
end
