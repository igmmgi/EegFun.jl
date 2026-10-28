# === AUTOREJECT: ADAPTIVE ARTIFACT REJECTION AND REPAIR ===
#
# Adaptation of the Autoreject algorithm based on:
# Jas, M., Engemann, D. A., Bekhti, Y., Raimondo, F., & Gramfort, A. (2017).
# "Autoreject: Automated artifact rejection for MEG and EEG data".
# NeuroImage, 159, 417-429. doi:10.1016/j.neuroimage.2017.06.030
#
# Adaptations and implementation details in EegFun.jl:
# 1. Channel thresholds (τ_c) are learned via K-fold cross-validation on peak-to-peak amplitude.
# 2. Fast candidate evaluation evaluates all candidate thresholds simultaneously per fold via BLAS matrix multiplication.
# 3. Consensus threshold (κ) and maximum repair threshold (ρ) determine epoch classification:
#    - Epochs with bad channels > κ * N are rejected.
#    - Epochs with bad channels <= κ * N are retained:
#      - If bad channels <= ρ * N, all bad channels are interpolated.
#      - If ρ * N < bad channels <= κ * N, the worst ρ * N channels (ranked by ptp / threshold)
#        are interpolated, while the remaining bad channels are retained with tracking in skipped_channels.
# 4. Spherical spline interpolation uses unit-sphere geometry (Perrin et al., 1989) with
#    operator caching across epochs that share identical bad-channel patterns.
# 5. Global mode learns a single threshold across all channels by summing per-channel cross-validation RMSEs.

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
- `skipped_channels::OrderedDict{Int,Vector{Symbol}}`: Mapping of global epoch index => channels left uninterpolated (in middle band)
- `rejection_info::Vector{EpochRejectionInfo}`: Standard EpochRejectionInfo per condition for pipeline compatibility
- `summary_table::DataFrame`: Summary DataFrame of repair/rejection metrics
- `channel_table::DataFrame`: Per-channel threshold and artifact occurrence table
- `n_folds::Int`: Number of cross-validation folds evaluated
- `n_candidates::Int`: Number of candidate thresholds evaluated
- `seed::Union{Int,Nothing}`: Random seed used for CV splits
- `repair_method::Symbol`: Repair method used (`:spherical_spline` or `:neighbor_interpolation`)
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
    skipped_channels::OrderedDict{Int,Vector{Symbol}}
    rejection_info::Vector{EpochRejectionInfo}
    summary_table::DataFrame
    channel_table::DataFrame
    n_folds::Int
    n_candidates::Int
    seed::Union{Int,Nothing}
    repair_method::Symbol
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
    print(io, "folds=$(res.n_folds), ")
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
    _candidate_thresholds(ptp_values::AbstractVector{<:Real}, n_candidates::Int,
                          quant_range::Tuple{Float64,Float64} = (0.40, 1.0)) -> Vector{Float64}

Generate sorted, unique candidate thresholds from quantiles of peak-to-peak amplitudes.
Always includes a candidate slightly above the maximum peak-to-peak amplitude (1.01 * max)
so that completely clean channels have an option to reject 0 epochs.
"""
function _candidate_thresholds(ptp_values::AbstractVector{<:Real}, n_candidates::Int, quant_range::Tuple{Float64,Float64} = (0.40, 1.0))
    isempty(ptp_values) && return [100.0]

    min_p, max_p = extrema(ptp_values)
    if min_p ≈ max_p || n_candidates <= 1
        return [Float64(max_p > 0 ? max_p * 1.01 : 100.0)]
    end

    q_levels = range(quant_range[1], quant_range[2], length = max(5, n_candidates - 1))
    candidates = Float64.(unique(quantile(ptp_values, q_levels)))
    push!(candidates, Float64(max_p * 1.01))
    filter!(x -> x > 0 && isfinite(x), candidates)

    if isempty(candidates)
        candidates = [Float64(max(100.0, max_p * 1.01))]
    end

    sort!(candidates)
    unique!(candidates)
    return candidates
end

"""
    _cv_candidate_errors(X::Matrix{Float64}, ptp_col::Vector{Float64},
                         candidates::Vector{Float64}, folds::Vector{Vector{Int}}) -> Vector{Float64}

Evaluate cross-validation RMSE across all candidate thresholds simultaneously using BLAS matrix multiplication.
`X` has dimensions `(n_timepoints, n_epochs)`.
"""
function _cv_candidate_errors(X::Matrix{Float64}, ptp_col::Vector{Float64}, candidates::Vector{Float64}, folds::Vector{Vector{Int}})
    n_timepoints, n_epochs = size(X)
    n_cands = length(candidates)
    errors = zeros(Float64, n_cands)
    in_test = falses(n_epochs)

    for test_idx in folds
        isempty(test_idx) && continue
        fill!(in_test, false)
        in_test[test_idx] .= true

        # Robust test median across test epochs for this fold: (n_timepoints,)
        test_med = vec(median(view(X, :, test_idx); dims = 2))

        # Weight matrix: W is (n_epochs × n_cands)
        # Column c contains weights for training epochs where ptp <= candidates[c]
        W = zeros(Float64, n_epochs, n_cands)
        @inbounds for c = 1:n_cands
            th = candidates[c]
            for ep = 1:n_epochs
                if !in_test[ep] && ptp_col[ep] <= th
                    W[ep, c] = 1.0
                end
            end
        end

        counts = vec(sum(W; dims = 1))
        for c = 1:n_cands
            if counts[c] == 0
                errors[c] = Inf
            else
                W[:, c] ./= counts[c]
            end
        end

        # M is (n_timepoints × n_cands) containing training means for all candidates in one BLAS GEMM call
        M = X * W

        for c = 1:n_cands
            if isfinite(errors[c])
                diff_norm_sq = sum(abs2, view(M, :, c) .- test_med)
                errors[c] += sqrt(diff_norm_sq / n_timepoints)
            end
        end
    end

    return errors
end

"""
    _best_candidate_index(errors::Vector{Float64}, tol::Float64 = 1e-7) -> Int

Select the index of the minimum cross-validation error among candidate thresholds.
If multiple candidates achieve the minimum within numerical tolerance `tol` (e.g. due to
identical training partitions or floating-point noise from multi-threaded BLAS GEMM),
consistently selects the first (lowest) candidate threshold for deterministic results across platforms.
"""
function _best_candidate_index(errors::Vector{Float64}, tol::Float64 = 1e-7)
    all(isinf, errors) && return length(errors)
    best_idx = 1
    min_e = errors[1]
    @inbounds for i = 2:length(errors)
        if errors[i] < min_e - tol
            min_e = errors[i]
            best_idx = i
        end
    end
    return best_idx
end

"""
    _repair_epochs_cached!(epochs_data::Vector{DataFrame}, epochs_to_repair::Vector{Int},
                           repair_matrix::BitMatrix, target_channels::Vector{Symbol}, layout::Layout;
                           repair_method::Symbol = :spherical_spline, m::Int = 4, lambda::Float64 = 1e-5,
                           sample_indices = Colon())
                           -> OrderedDict{Int, Vector{Symbol}}

Interpolate bad channels in epochs marked for repair, caching the interpolation operator
for epochs that share the identical pattern of bad channels. Writes only to repaired channel columns.
"""
function _repair_epochs_cached!(
    epochs_data::Vector{DataFrame},
    epochs_to_repair::Vector{Int},
    repair_matrix::BitMatrix,
    target_channels::Vector{Symbol},
    layout::Layout;
    repair_method::Symbol = :spherical_spline,
    m::Int = 4,
    lambda::Real = 1e-5,
    sample_indices = Colon(),
)
    repaired_channels = OrderedDict{Int,Vector{Symbol}}()
    isempty(epochs_to_repair) && return repaired_channels

    n_channels = length(target_channels)
    has_spline_coords = false
    local U_all::Matrix{Float64}
    local f_legendre::Vector{Float64}

    if repair_method == :spherical_spline
        ch_layout_indices = Int[]
        all_found = true
        for ch in target_channels
            idx = findfirst(==(ch), layout.data.label)
            if isnothing(idx)
                all_found = false
                break
            end
            push!(ch_layout_indices, idx)
        end

        if all_found && (hasproperty(layout.data, :inc) || hasproperty(layout.data, :x2))
            try
                U_all, _ = _spline_geometry(layout.data, ch_layout_indices)
                f_legendre = _g_factors(m, 50)
                has_spline_coords = true
            catch e
                @minimal_warning "Failed to construct spherical geometry from layout: $(sprint(showerror, e)). Falling back to neighbour interpolation if available."
                has_spline_coords = false
            end
        else
            @minimal_warning "Layout is missing electrode coordinates for one or more target channels. Falling back to neighbour interpolation if available."
            has_spline_coords = false
        end
    end

    # Cache for spherical spline interpolation operators: pattern => (good_indices, bad_indices, operator_matrix)
    spline_operator_cache = Dict{Vector{Symbol},Tuple{Vector{Int},Vector{Int},Matrix{Float64}}}()

    for ep_idx in epochs_to_repair
        bad_mask = @view repair_matrix[ep_idx, :]
        bad_idx = findall(bad_mask)
        isempty(bad_idx) && continue

        bad_chs = target_channels[bad_idx]
        pattern = sort(bad_chs)
        epoch_df = epochs_data[ep_idx]
        n_timepoints = nrow(epoch_df)

        if repair_method == :spherical_spline && has_spline_coords
            good_idx = setdiff(1:n_channels, bad_idx)

            if length(good_idx) >= 4
                if !haskey(spline_operator_cache, pattern)
                    try
                        U_good = U_all[good_idx, :]
                        U_bad = U_all[bad_idx, :]
                        op_mat = _spline_good_to_bad(U_good, U_bad, f_legendre, lambda)
                        spline_operator_cache[pattern] = (good_idx, bad_idx, op_mat)
                    catch e
                        @minimal_warning "Spherical spline operator generation failed for bad channels $pattern: $e. Falling back to neighbour interpolation."
                    end
                end

                if haskey(spline_operator_cache, pattern)
                    g_idx, b_idx, op_mat = spline_operator_cache[pattern]

                    # Extract only good channels
                    good_data = Matrix{Float64}(undef, length(g_idx), n_timepoints)
                    for (row_i, g_i) in enumerate(g_idx)
                        good_data[row_i, :] = epoch_df[sample_indices, target_channels[g_i]]
                    end

                    # Compute interpolated potentials for bad channels
                    repaired_data = op_mat * good_data

                    # Write only to bad channel columns in epoch DataFrame
                    for (row_i, b_i) in enumerate(b_idx)
                        epoch_df[sample_indices, target_channels[b_i]] = repaired_data[row_i, :]
                    end

                    repaired_channels[ep_idx] = pattern
                    continue
                end
            end
        end

        # Fallback to neighbour interpolation if spline coordinates not available or failed
        if !isnothing(layout.neighbours)
            epoch_mat = Matrix(epoch_df[:, target_channels])'
            actually_repaired = Symbol[]
            for bad_ch in pattern
                if haskey(layout.neighbours, bad_ch)
                    good_neighs = [n for n in layout.neighbours[bad_ch].channels if n ∉ pattern && n in target_channels]
                    if !isempty(good_neighs)
                        push!(actually_repaired, bad_ch)
                    end
                end
            end

            if !isempty(actually_repaired)
                _repair_channels_neighbor!(epoch_mat, target_channels, actually_repaired, layout.neighbours)
                for bad_ch in actually_repaired
                    ch_idx = findfirst(==(bad_ch), target_channels)
                    epoch_df[sample_indices, bad_ch] = epoch_mat[ch_idx, :]
                end
                repaired_channels[ep_idx] = actually_repaired
            end
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
    and rejects epochs with excessive bad channels (> consensus_perc).
  - `:global`: Learns a single global threshold across all channels by summing CV errors, reject-only.
- `channel_selection::Function`: Channel selector (default: `channels()`).
- `interval_selection::Interval`: Time window for peak-to-peak evaluation (default: `times()`).
- `n_folds::Int`: Number of cross-validation folds (default: 10, minimum: 2).
- `n_candidates::Int`: Number of candidate thresholds to evaluate per channel (default: 30).
- `candidate_quantiles::Tuple{Float64,Float64}`: Quantile range for candidate thresholds (default: `(0.40, 1.0)`).
- `consensus_perc::Float64`: Maximum fraction of bad channels in an epoch before the entire epoch is rejected (κ, default: 0.2).
- `max_repair_perc::Float64`: Maximum fraction of bad channels in an epoch that can be repaired (ρ, default: 0.1).
- `repair_method::Symbol`: Repair method (`:spherical_spline` or `:neighbor_interpolation`, default: `:spherical_spline`).
- `m::Int`: Legendre polynomial order for spherical splines (default: 4).
- `lambda::Float64`: Regularization parameter for spherical splines (default: 1e-5).
- `seed::Union{Int,Nothing}`: Random seed for CV splits (default: 42).
- `reject::Bool`: Whether to remove rejected epochs from `dat` (default: true).
- `verbose::Bool`: Whether to log progress (default: true).

# Performance
Local mode utilizes `Threads.@threads` for multi-threaded channel threshold estimation when Julia
is started with multiple threads (`julia -t auto` or `JULIA_NUM_THREADS`).

# Returns
- `AutoRejectResult`: Complete result with learned thresholds, repair/rejection counts, and summary tables.

# References
Jas, M., Engemann, D. A., Bekhti, Y., Raimondo, F., & Gramfort, A. (2017). Autoreject: Automated artifact rejection
for MEG and EEG data. NeuroImage, 159, 417-429.
"""
function autoreject!(
    dat::AbstractVector{<:EpochData};
    mode::Symbol = :local,
    channel_selection = channels(),
    interval_selection::Interval = times(),
    n_folds::Int = 10,
    n_candidates::Int = 30,
    candidate_quantiles::Tuple{Real,Real} = (0.40, 1.0),
    consensus_perc::Real = 0.2,
    max_repair_perc::Real = 0.1,
    repair_method::Symbol = :spherical_spline,
    m::Int = 4,
    lambda::Real = 1e-5,
    seed::Union{Integer,Nothing} = 42,
    reject::Bool = true,
    verbose::Bool = true,
)
    isempty(dat) && @minimal_error("Cannot run autoreject on empty epoch vector")

    mode ∉ [:local, :global] && @minimal_error("Unknown mode: $mode. Use :local or :global.")
    repair_method ∉ [:spherical_spline, :neighbor_interpolation] &&
        @minimal_error("Unknown repair method: $repair_method. Allowed: :spherical_spline, :neighbor_interpolation")
    n_folds < 2 && @minimal_error("n_folds must be at least 2, got $n_folds")
    (consensus_perc < 0.0 || consensus_perc > 1.0) && @minimal_error("consensus_perc must be between 0.0 and 1.0")
    (max_repair_perc < 0.0 || max_repair_perc > 1.0) && @minimal_error("max_repair_perc must be between 0.0 and 1.0")

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
            OrderedDict{Int,Vector{Symbol}}(),
            EpochRejectionInfo[],
            DataFrame(),
            DataFrame(),
            n_folds,
            n_candidates,
            seed,
            repair_method,
        )
        return empty_res
    end

    # Determine sample window
    first_epoch = all_epoch_dfs[1]
    combined_sel = _combine_interval_sample(interval_selection, samples())
    sample_indices = get_selected_samples(first_epoch, combined_sel)
    isempty(sample_indices) && @minimal_error("Selected interval contains no samples")
    n_timepoints = length(sample_indices)

    # 1. Compute Peak-to-Peak (PTP) amplitude matrix in a single pass using extrema
    ptp_matrix = zeros(Float64, n_total_epochs, n_channels)
    for (ep_idx, ep_df) in enumerate(all_epoch_dfs)
        for (ch_idx, ch) in enumerate(selected_channels)
            col = ep_df[!, ch]
            s_view = @view col[sample_indices]
            min_v, max_v = extrema(s_view)
            ptp_matrix[ep_idx, ch_idx] = Float64(max_v - min_v)
        end
    end

    # Fallback for datasets with very few epochs (< 4):
    # Cross-validation cannot reliably estimate thresholds on fewer than 4 trials.
    # Rather than flagging arbitrary percentiles, skip rejection with a warning.
    if n_total_epochs < 4
        @minimal_warning "Too few epochs ($n_total_epochs < 4) for cross-validated Autoreject threshold estimation. Skipping adaptive artifact rejection."

        rejection_info_per_condition = EpochRejectionInfo[]
        for ep_dat in dat
            cond_num = ep_dat.condition
            cond_name = ep_dat.condition_name
            local_total = n_epochs(ep_dat)
            info = EpochInfo(cond_num, cond_name, local_total)
            push!(
                rejection_info_per_condition,
                EpochRejectionInfo(
                    "autoreject",
                    info,
                    0,
                    0.0,
                    Rejection[],
                    0.0,
                    nothing,
                    Rejection[],
                    OrderedDict{Int,Vector{Symbol}}(),
                    OrderedDict{Int,Vector{Symbol}}(),
                ),
            )
        end

        summary_df = DataFrame(Metric = ["Total Epochs", "Status"], Value = [string(n_total_epochs), "Skipped (< 4 epochs)"])
        channel_df = DataFrame(
            channel = selected_channels,
            threshold_uV = fill(Inf, n_channels),
            bad_epochs_count = zeros(Int, n_channels),
            bad_epochs_pct = zeros(Float64, n_channels),
        )

        return AutoRejectResult(
            Dict(ch => Inf for ch in selected_channels),
            mode,
            consensus_perc,
            max_repair_perc,
            selected_channels,
            ptp_matrix,
            falses(n_total_epochs, n_channels),
            collect(1:n_total_epochs),
            Int[],
            Int[],
            OrderedDict{Int,Vector{Symbol}}(),
            OrderedDict{Int,Vector{Symbol}}(),
            rejection_info_per_condition,
            summary_df,
            channel_df,
            n_folds,
            n_candidates,
            seed,
            repair_method,
        )
    end

    verbose && @info "Running Autoreject (mode = :$mode) on $n_total_epochs epochs across $n_channels channels..."

    # 2. Threshold Estimation via Cross-Validation
    thresholds = Dict{Symbol,Float64}()

    eff_folds = n_folds
    if n_total_epochs < eff_folds
        eff_folds = max(2, div(n_total_epochs, 2))
        verbose && @info "Adjusted cross-validation folds to $eff_folds due to trial count ($n_total_epochs)"
    end

    folds = _kfold_splits(n_total_epochs, eff_folds, seed)

    if mode == :global
        # Global mode: single threshold across all channels evaluated simultaneously
        all_ptps = vec(ptp_matrix)
        global_candidates = _candidate_thresholds(all_ptps, n_candidates, candidate_quantiles)
        total_candidate_errors = zeros(Float64, length(global_candidates))

        # Sum per-channel CV errors across global candidates without repeated matrix allocations
        for ch_idx = 1:n_channels
            ch = selected_channels[ch_idx]
            X = Matrix{Float64}(undef, n_timepoints, n_total_epochs)
            for (ep_i, ep_df) in enumerate(all_epoch_dfs)
                X[:, ep_i] = ep_df[sample_indices, ch]
            end
            err = _cv_candidate_errors(X, ptp_matrix[:, ch_idx], global_candidates, folds)
            total_candidate_errors .+= err
        end

        best_idx = _best_candidate_index(total_candidate_errors)
        best_global_th = global_candidates[best_idx]

        for ch in selected_channels
            thresholds[ch] = best_global_th
        end

    else
        # Local mode: optimize threshold per channel via multi-threading
        # Each thread allocates its own (n_timepoints × n_total_epochs) matrix, keeping peak memory bounded
        threshold_results = Vector{Float64}(undef, n_channels)

        Threads.@threads for ch_idx = 1:n_channels
            ch = selected_channels[ch_idx]
            ptp_col = ptp_matrix[:, ch_idx]
            candidates = _candidate_thresholds(ptp_col, n_candidates, candidate_quantiles)

            X = Matrix{Float64}(undef, n_timepoints, n_total_epochs)
            for (ep_i, ep_df) in enumerate(all_epoch_dfs)
                X[:, ep_i] = ep_df[sample_indices, ch]
            end

            err = _cv_candidate_errors(X, ptp_col, candidates, folds)
            best_idx = _best_candidate_index(err)
            threshold_results[ch_idx] = candidates[best_idx]
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

    # Warn if any channel has artifacts in a very large percentage of epochs (>= 30%)
    for ch_idx = 1:n_channels
        ch = selected_channels[ch_idx]
        bad_count = sum(@view bad_matrix[:, ch_idx])
        bad_pct = 100.0 * bad_count / max(1, n_total_epochs)
        if bad_pct >= 30.0
            @minimal_warning "Channel :$ch is corrupted in $(round(bad_pct, digits=1))% of epochs ($bad_count/$n_total_epochs). Consider continuous bad-channel repair or rejection before epoching."
        end
    end

    # 4. Classify Epochs: Clean, Repairable, Rejected (Middle band handling between ρ and κ)
    epochs_clean = Int[]
    epochs_repaired = Int[]
    epochs_rejected = Int[]
    skipped_channels = OrderedDict{Int,Vector{Symbol}}()

    # Consensus rejection count (k_max): reject if bad channels > k_max
    k_max = round(Int, consensus_perc * n_channels)
    # Maximum repairable channel count (r_max)
    r_max = mode == :global ? 0 : round(Int, max_repair_perc * n_channels)

    # repair_matrix tracks which channels will actually be repaired
    repair_matrix = copy(bad_matrix)

    for ep_idx = 1:n_total_epochs
        bad_idx = findall(@view bad_matrix[ep_idx, :])
        n_bad = length(bad_idx)

        if n_bad == 0
            push!(epochs_clean, ep_idx)
        elseif mode == :global || n_bad > k_max
            # Exceeds consensus threshold: epoch is rejected
            push!(epochs_rejected, ep_idx)
        else
            # 0 < n_bad <= k_max: epoch is retained!
            push!(epochs_repaired, ep_idx)

            if n_bad <= r_max
                # All bad channels will be repaired
            else
                # Middle band (r_max < n_bad <= k_max):
                # Interpolate the worst r_max channels by relative peak-to-peak amplitude (ptp / threshold)
                bad_severity = [(ch_i, ptp_matrix[ep_idx, ch_i] / thresholds[selected_channels[ch_i]]) for ch_i in bad_idx]
                sort!(bad_severity, by = x -> x[2], rev = true)

                channels_to_skip_idx = [x[1] for x in bad_severity[(r_max+1):end]]

                # Clear repair flag for skipped channels
                for ch_i in channels_to_skip_idx
                    repair_matrix[ep_idx, ch_i] = false
                end

                skipped_channels[ep_idx] = [selected_channels[ch_i] for ch_i in channels_to_skip_idx]
            end
        end
    end

    # 5. Repair Channels for Repairable Epochs
    repaired_channels = OrderedDict{Int,Vector{Symbol}}()
    if !isempty(epochs_repaired) && mode != :global
        verbose && @info "Repairing bad channels across $(length(epochs_repaired)) epochs using $repair_method..."
        repaired_channels = _repair_epochs_cached!(
            all_epoch_dfs,
            epochs_repaired,
            repair_matrix,
            selected_channels,
            layout;
            repair_method = repair_method,
            m = m,
            lambda = lambda,
            sample_indices = sample_indices,
        )

        # Handle failed repairs: any epoch where repair could not be performed is moved to epochs_rejected
        failed_repairs = [ep for ep in epochs_repaired if !haskey(repaired_channels, ep)]
        if !isempty(failed_repairs)
            filter!(ep -> haskey(repaired_channels, ep), epochs_repaired)
            append!(epochs_rejected, failed_repairs)
            unique!(epochs_rejected)
            verbose && @minimal_warning "$(length(failed_repairs)) epochs could not be repaired and were moved to rejected epochs."
        end
    end

    # 6. Map back to conditions and construct EpochRejectionInfo
    rejection_info_per_condition = EpochRejectionInfo[]
    rejected_set = Set(epochs_rejected)

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
            if g_idx in rejected_set
                bad_chs = findall(@view bad_matrix[g_idx, :])
                if !isempty(bad_chs)
                    for ch_i in bad_chs
                        push!(c_rejections, Rejection(selected_channels[ch_i], m_loc_idx))
                    end
                else
                    push!(c_rejections, Rejection(selected_channels[1], m_loc_idx))
                end
            end

            # If repaired
            if haskey(repaired_channels, g_idx)
                c_repaired[m_loc_idx] = repaired_channels[g_idx]
            end

            # If skipped channels (middle band)
            if haskey(skipped_channels, g_idx)
                c_skipped[m_loc_idx] = skipped_channels[g_idx]
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
            "CV Folds",
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
            Printf.@sprintf("%.2f (%d chs)", max_repair_perc, r_max),
            string(eff_folds),
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
        max_repair_perc,
        selected_channels,
        ptp_matrix,
        bad_matrix,
        epochs_clean,
        epochs_repaired,
        epochs_rejected,
        repaired_channels,
        skipped_channels,
        rejection_info_per_condition,
        summary_df,
        channel_df,
        eff_folds,
        n_candidates,
        seed,
        repair_method,
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

Non-mutating version of `autoreject!`. Deep-copies `dat` and returns `(cleaned_data, result)`.
Leaves input `dat` completely unmodified.
"""
function autoreject(dat::AbstractVector{<:EpochData}; kwargs...)
    dat_copy = [copy(d) for d in dat]
    res = autoreject!(dat_copy; kwargs...)
    return dat_copy, res
end

function autoreject(dat::EpochData; kwargs...)
    dat_copy = copy(dat)
    res = autoreject!(dat_copy; kwargs...)
    return dat_copy, res
end
