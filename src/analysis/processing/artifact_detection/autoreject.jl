# Adaptive artifact rejection and channel repair (Jas et al., 2017; NeuroImage, 159, 417-429).

"""
    AutoRejectResult

Result of adaptive artifact rejection and channel repair ([`autoreject`](@ref) / [`autoreject!`](@ref)).

# Fields
- `thresholds::Dict{Symbol,Float64}`: Learned peak-to-peak amplitude threshold per channel (μV).
- `mode::Symbol`: Rejection mode (`:local` or `:global`).
- `consensus_perc::Float64`: Consensus rejection fraction κ (`bad > κ * N` rejects epoch).
- `max_repair_perc::Float64`: Maximum repair fraction ρ (`bad <= ρ * N` interpolates all; middle band repairs worst `ρ * N`).
- `channels::Vector{Symbol}`: Analyzed channels.
- `ptp::Matrix{Float64}`: Peak-to-peak amplitude matrix (`n_epochs × n_channels`).
- `bad_matrix::BitMatrix`: Channel × epoch threshold violations.
- `epochs_clean::Vector{Int}`, `epochs_repaired::Vector{Int}`, `epochs_rejected::Vector{Int}`: Global epoch indices.
- `repaired_channels::OrderedDict{Int,Vector{Symbol}}`, `skipped_channels::OrderedDict{Int,Vector{Symbol}}`: Per-epoch channel tracking.
- `rejection_info::Vector{EpochRejectionInfo}`: Per-condition `EpochRejectionInfo`.
- `summary_table::DataFrame`, `channel_table::DataFrame`: Summary and per-channel tables.
- `n_folds::Int`, `n_candidates::Int`, `seed::Union{Int,Nothing}`, `repair_method::Symbol`: Run configuration.
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
    n_kept = length(res.epochs_clean) + length(res.epochs_repaired)
    n_tot = n_kept + length(res.epochs_rejected)
    retention = n_tot > 0 ? round(100 * n_kept / n_tot; digits = 1) : 0.0
    th_vals = collect(values(res.thresholds))
    min_th, med_th, max_th = isempty(th_vals) ? (0.0, 0.0, 0.0) : (minimum(th_vals), median(th_vals), maximum(th_vals))

    print(
        io,
        "AutoRejectResult(mode=:$(res.mode), folds=$(res.n_folds), epochs=$(n_kept)/$(n_tot) kept ($(retention)%), ",
        "repaired=$(length(res.epochs_repaired)), rejected=$(length(res.epochs_rejected)), ",
        "thresholds=[min=$(round(min_th, digits=1)), med=$(round(med_th, digits=1)), max=$(round(max_th, digits=1))] μV)",
    )
end

function _kfold_splits(n_total::Int, n_folds::Int, seed::Union{Int,Nothing})
    rng = isnothing(seed) ? Random.default_rng() : Random.MersenneTwister(seed)
    perm = randperm(rng, n_total)
    folds = Vector{Vector{Int}}(undef, n_folds)
    fold_size, rem_samples = divrem(n_total, n_folds)
    start_idx = 1
    for f = 1:n_folds
        end_idx = start_idx + fold_size + (f <= rem_samples ? 1 : 0) - 1
        folds[f] = perm[start_idx:end_idx]
        start_idx = end_idx + 1
    end
    return folds
end

function _candidate_thresholds(ptp_values::AbstractVector{<:Real}, n_candidates::Int, quant_range::Tuple{Real,Real} = (0.40, 1.0))
    isempty(ptp_values) && return [100.0]
    min_p, max_p = extrema(ptp_values)
    (min_p ≈ max_p || n_candidates <= 1) && return [Float64(max_p > 0 ? max_p * 1.01 : 100.0)]

    q_levels = range(Float64(quant_range[1]), Float64(quant_range[2]); length = max(5, n_candidates - 1))
    candidates = Float64.(unique(quantile(ptp_values, q_levels)))
    push!(candidates, Float64(max_p * 1.01))
    filter!(x -> x > 0 && isfinite(x), candidates)
    isempty(candidates) && return [Float64(max(100.0, max_p * 1.01))]
    sort!(candidates)
    unique!(candidates)
    return candidates
end

function _cv_candidate_errors(X::Matrix{Float64}, ptp_col::AbstractVector{Float64}, candidates::Vector{Float64}, folds::Vector{Vector{Int}})
    n_timepoints, n_epochs = size(X)
    n_cands = length(candidates)
    errors = zeros(Float64, n_cands)
    in_test = falses(n_epochs)
    W = zeros(Float64, n_epochs, n_cands)
    M = zeros(Float64, n_timepoints, n_cands)
    test_med = Vector{Float64}(undef, n_timepoints)
    med_buf = Vector{Float64}(undef, maximum(length, folds; init = 0))
    inv_nt = 1.0 / n_timepoints

    for test_idx in folds
        k_len = length(test_idx)
        k_len == 0 && continue
        fill!(in_test, false)
        in_test[test_idx] .= true

        buf = view(med_buf, 1:k_len)
        @inbounds for t = 1:n_timepoints
            for j = 1:k_len
                buf[j] = X[t, test_idx[j]]
            end
            test_med[t] = median!(buf)
        end

        @inbounds for c = 1:n_cands
            th = candidates[c]
            cnt = 0
            for ep = 1:n_epochs
                keep = !in_test[ep] && ptp_col[ep] <= th
                W[ep, c] = ifelse(keep, 1.0, 0.0)
                cnt += keep
            end
            if cnt == 0
                errors[c] = Inf
            else
                inv_cnt = 1.0 / cnt
                @simd for ep = 1:n_epochs
                    W[ep, c] *= inv_cnt
                end
            end
        end

        mul!(M, X, W)

        @inbounds for c = 1:n_cands
            if isfinite(errors[c])
                s = 0.0
                @simd for t = 1:n_timepoints
                    d = M[t, c] - test_med[t]
                    s += d * d
                end
                errors[c] += sqrt(s * inv_nt)
            end
        end
    end
    return errors
end

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

@inline function _ptp_range(col::AbstractVector{<:Real}, sample_indices)
    mn, mx = extrema(@view col[sample_indices])
    return Float64(mx - mn)
end

@inline _copy_sample_slice!(dst::AbstractVector{Float64}, src::AbstractVector{<:Real}, sample_indices) =
    copyto!(dst, @view src[sample_indices])

function _fill_channel_epochs!(X::Matrix{Float64}, dfs::Vector{DataFrame}, col_idx::Int, sample_indices)
    for (ep_i, ep_df) in enumerate(dfs)
        _copy_sample_slice!(view(X, :, ep_i), ep_df[!, col_idx], sample_indices)
    end
    return X
end

function _channel_cv_threshold(
    all_epoch_dfs::Vector{DataFrame},
    col_idx::Int,
    sample_indices,
    ptp_col::AbstractVector{Float64},
    n_candidates::Int,
    candidate_quantiles::Tuple{Real,Real},
    folds::Vector{Vector{Int}},
)
    candidates = _candidate_thresholds(ptp_col, n_candidates, candidate_quantiles)
    X = Matrix{Float64}(undef, length(sample_indices), length(all_epoch_dfs))
    _fill_channel_epochs!(X, all_epoch_dfs, col_idx, sample_indices)
    err = _cv_candidate_errors(X, ptp_col, candidates, folds)
    return candidates[_best_candidate_index(err)]
end

function _repair_epochs_cached!(
    epochs_data::Vector{DataFrame},
    epochs_to_repair::Vector{Int},
    repair_matrix::BitMatrix,
    target_channels::Vector{Symbol},
    layout::Layout;
    repair_method::Symbol = :spherical_spline,
    m::Int = 4,
    lambda::Real = 1e-5,
)
    repaired_channels = OrderedDict{Int,Vector{Symbol}}()
    isempty(epochs_to_repair) && return repaired_channels

    n_channels = length(target_channels)
    has_spline_coords = false
    local U_all::Matrix{Float64}
    local f_legendre::Vector{Float64}

    if repair_method == :spherical_spline
        ch_layout_indices = Int[]
        for ch in target_channels
            idx = findfirst(==(ch), layout.data.label)
            isnothing(idx) && break
            push!(ch_layout_indices, idx)
        end

        if length(ch_layout_indices) == n_channels && (hasproperty(layout.data, :inc) || hasproperty(layout.data, :x2))
            try
                U_all, _ = _spline_geometry(layout.data, ch_layout_indices)
                f_legendre = _g_factors(m, 50)
                has_spline_coords = true
            catch e
                @minimal_warning "Failed to construct spherical geometry from layout: $(sprint(showerror, e)). Falling back to neighbor interpolation if available."
            end
        else
            @minimal_warning "Layout is missing electrode coordinates for one or more target channels. Falling back to neighbor interpolation if available."
        end
    end

    spline_cache = Dict{Vector{Symbol},Tuple{Vector{Int},Vector{Int},Matrix{Float64}}}()

    for ep_idx in epochs_to_repair
        bad_idx = findall(@view repair_matrix[ep_idx, :])
        isempty(bad_idx) && continue

        pattern = sort(target_channels[bad_idx])
        epoch_df = epochs_data[ep_idx]
        n_tp = nrow(epoch_df)

        if repair_method == :spherical_spline && has_spline_coords
            good_idx = setdiff(1:n_channels, bad_idx)
            if length(good_idx) >= 4
                if !haskey(spline_cache, pattern)
                    try
                        op_T = collect(transpose(_spline_good_to_bad(U_all[good_idx, :], U_all[bad_idx, :], f_legendre, lambda)))
                        spline_cache[pattern] = (good_idx, bad_idx, op_T)
                    catch e
                        @minimal_warning "Spherical spline operator generation failed for bad channels $pattern: $e. Falling back to neighbor interpolation."
                    end
                end

                if haskey(spline_cache, pattern)
                    g_idx, b_idx, op_T = spline_cache[pattern]
                    good_data = Matrix{Float64}(undef, n_tp, length(g_idx))
                    for (col_i, g_i) in enumerate(g_idx)
                        copyto!(view(good_data, :, col_i), epoch_df[!, target_channels[g_i]])
                    end
                    repaired_data = good_data * op_T
                    for (col_i, b_i) in enumerate(b_idx)
                        epoch_df[!, target_channels[b_i]] = repaired_data[:, col_i]
                    end
                    repaired_channels[ep_idx] = pattern
                    continue
                end
            end
        end

        if !isnothing(layout.neighbors)
            plans = _resolve_neighbor_plans(target_channels, pattern, layout.neighbors; bad_channels = pattern, log_repairs = false)
            if !isempty(plans)
                _apply_neighbor_plans!(epoch_df, plans)
                repaired_channels[ep_idx] = [bad_ch for (_, bad_ch, _, _, _) in plans]
            end
        end
    end
    return repaired_channels
end

"""
    autoreject!(dat::AbstractVector{<:EpochData}; kwargs...) -> AutoRejectResult
    autoreject!(dat::EpochData; kwargs...) -> AutoRejectResult

Adaptive cross-validated artifact rejection and channel repair (Jas et al., 2017), mutating `dat` in-place.

# Keyword Arguments
- `mode::Symbol = :local`: `:local` (per-channel CV thresholds + repair/reject) or `:global` (single global threshold, reject-only).
- `channel_selection = channels()`: Channel selector.
- `interval_selection::Interval = times()`: Time window for peak-to-peak evaluation.
- `n_folds::Int = 10`, `n_candidates::Int = 30`, `candidate_quantiles = (0.40, 1.0)`: CV search parameters.
- `consensus_perc::Real = 0.2`, `max_repair_perc::Real = 0.1`: Consensus rejection (κ) and max repair (ρ) channel fractions.
- `repair_method::Symbol = :spherical_spline`: `:spherical_spline` (`m = 4`, `lambda = 1e-5`) or `:neighbor_interpolation`.
- `seed::Union{Integer,Nothing} = 42`, `reject::Bool = true`, `verbose::Bool = true`.
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
    mode in (:local, :global) || @minimal_error("Unknown mode: $mode. Use :local or :global.")
    repair_method in (:spherical_spline, :neighbor_interpolation) ||
        @minimal_error("Unknown repair method: $repair_method. Allowed: :spherical_spline, :neighbor_interpolation")
    n_folds < 2 && @minimal_error("n_folds must be at least 2, got $n_folds")
    (0.0 <= consensus_perc <= 1.0) || @minimal_error("consensus_perc must be between 0.0 and 1.0")
    (0.0 <= max_repair_perc <= 1.0) || @minimal_error("max_repair_perc must be between 0.0 and 1.0")

    layout = dat[1].layout
    selected_channels = get_selected_channels(dat[1], channel_selection; include_meta = false, include_extra = false)
    isempty(selected_channels) && @minimal_error("No channels selected for autoreject")
    n_channels = length(selected_channels)

    all_epoch_dfs = DataFrame[]
    epoch_mapping = Tuple{Int,Int}[]
    for (cond_idx, ep_dat) in enumerate(dat)
        for loc_idx = 1:n_epochs(ep_dat)
            push!(all_epoch_dfs, ep_dat.data[loc_idx])
            push!(epoch_mapping, (cond_idx, loc_idx))
        end
    end

    n_total_epochs = length(all_epoch_dfs)
    if n_total_epochs == 0
        @minimal_warning("No epochs found to run autoreject.")
        return AutoRejectResult(
            Dict{Symbol,Float64}(),
            mode,
            Float64(consensus_perc),
            Float64(max_repair_perc),
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
    end

    first_epoch = all_epoch_dfs[1]
    sample_indices = get_selected_samples(first_epoch, _combine_interval_sample(interval_selection, samples()))
    isempty(sample_indices) && @minimal_error("Selected interval contains no samples")
    n_timepoints = length(sample_indices)
    col_indices = [columnindex(first_epoch, ch) for ch in selected_channels]

    ptp_matrix = Matrix{Float64}(undef, n_total_epochs, n_channels)
    for (ep_idx, ep_df) in enumerate(all_epoch_dfs)
        for ch_idx = 1:n_channels
            ptp_matrix[ep_idx, ch_idx] = _ptp_range(ep_df[!, col_indices[ch_idx]], sample_indices)
        end
    end

    if n_total_epochs < 4
        @minimal_warning "Too few epochs ($n_total_epochs < 4) for cross-validated Autoreject threshold estimation. Skipping adaptive artifact rejection."
        rej_infos = [
            EpochRejectionInfo(
                "autoreject",
                EpochInfo(d.condition, d.condition_name, n_epochs(d)),
                0,
                0.0,
                Rejection[],
                0.0,
                nothing,
                Rejection[],
                OrderedDict{Int,Vector{Symbol}}(),
                OrderedDict{Int,Vector{Symbol}}(),
            ) for d in dat
        ]
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
            Float64(consensus_perc),
            Float64(max_repair_perc),
            selected_channels,
            ptp_matrix,
            falses(n_total_epochs, n_channels),
            collect(1:n_total_epochs),
            Int[],
            Int[],
            OrderedDict{Int,Vector{Symbol}}(),
            OrderedDict{Int,Vector{Symbol}}(),
            rej_infos,
            summary_df,
            channel_df,
            n_folds,
            n_candidates,
            seed,
            repair_method,
        )
    end

    verbose && @info "Running Autoreject (mode = :$mode) on $n_total_epochs epochs across $n_channels channels..."

    eff_folds = n_total_epochs < n_folds ? max(2, div(n_total_epochs, 2)) : n_folds
    eff_folds != n_folds && verbose && @info "Adjusted cross-validation folds to $eff_folds due to trial count ($n_total_epochs)"
    folds = _kfold_splits(n_total_epochs, eff_folds, seed)

    th_vals = Vector{Float64}(undef, n_channels)
    if mode == :global
        global_candidates = _candidate_thresholds(vec(ptp_matrix), n_candidates, candidate_quantiles)
        total_errors = zeros(Float64, length(global_candidates))
        X_global = Matrix{Float64}(undef, n_timepoints, n_total_epochs)
        for ch_idx = 1:n_channels
            _fill_channel_epochs!(X_global, all_epoch_dfs, col_indices[ch_idx], sample_indices)
            total_errors .+= _cv_candidate_errors(X_global, view(ptp_matrix, :, ch_idx), global_candidates, folds)
        end
        fill!(th_vals, global_candidates[_best_candidate_index(total_errors)])
    else
        Threads.@threads for ch_idx = 1:n_channels
            th_vals[ch_idx] = _channel_cv_threshold(
                all_epoch_dfs,
                col_indices[ch_idx],
                sample_indices,
                view(ptp_matrix, :, ch_idx),
                n_candidates,
                candidate_quantiles,
                folds,
            )
        end
    end
    thresholds = Dict{Symbol,Float64}(selected_channels[i] => th_vals[i] for i = 1:n_channels)

    bad_matrix = ptp_matrix .> reshape(th_vals, 1, :)
    for ch_idx = 1:n_channels
        bad_count = count(@view bad_matrix[:, ch_idx])
        bad_pct = 100.0 * bad_count / n_total_epochs
        if bad_pct >= 30.0
            @minimal_warning "Channel :$(selected_channels[ch_idx]) is corrupted in $(round(bad_pct, digits=1))% of epochs ($bad_count/$n_total_epochs). Consider continuous bad-channel repair or rejection before epoching."
        end
    end

    epochs_clean = Int[]
    epochs_repaired = Int[]
    epochs_rejected = Int[]
    skipped_channels = OrderedDict{Int,Vector{Symbol}}()
    k_max = round(Int, consensus_perc * n_channels)
    r_max = mode == :global ? 0 : round(Int, max_repair_perc * n_channels)
    repair_matrix = copy(bad_matrix)

    for ep_idx = 1:n_total_epochs
        bad_idx = findall(@view bad_matrix[ep_idx, :])
        n_bad = length(bad_idx)
        if n_bad == 0
            push!(epochs_clean, ep_idx)
        elseif mode == :global || n_bad > k_max
            push!(epochs_rejected, ep_idx)
        else
            push!(epochs_repaired, ep_idx)
            if n_bad > r_max
                bad_severity = [(ch_i, ptp_matrix[ep_idx, ch_i] / th_vals[ch_i]) for ch_i in bad_idx]
                sort!(bad_severity; by = last, rev = true)
                skip_idx = [ch_i for (ch_i, _) in @view bad_severity[(r_max+1):end]]
                repair_matrix[ep_idx, skip_idx] .= false
                skipped_channels[ep_idx] = selected_channels[skip_idx]
            end
        end
    end

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
        )
        failed_repairs = [ep for ep in epochs_repaired if !haskey(repaired_channels, ep)]
        if !isempty(failed_repairs)
            filter!(ep -> haskey(repaired_channels, ep), epochs_repaired)
            append!(epochs_rejected, failed_repairs)
            unique!(epochs_rejected)
            verbose && @minimal_warning "$(length(failed_repairs)) epochs could not be repaired and were moved to rejected epochs."
        end
    end

    rejection_info_per_condition = Vector{EpochRejectionInfo}(undef, length(dat))
    rejected_set = Set(epochs_rejected)
    median_th = median(th_vals)

    for (c_idx, ep_dat) in enumerate(dat)
        c_rejections = Rejection[]
        c_repaired = OrderedDict{Int,Vector{Symbol}}()
        c_skipped = OrderedDict{Int,Vector{Symbol}}()

        for (g_idx, (m_c_idx, m_loc_idx)) in enumerate(epoch_mapping)
            m_c_idx == c_idx || continue
            if g_idx in rejected_set
                bad_chs = findall(@view bad_matrix[g_idx, :])
                if isempty(bad_chs)
                    push!(c_rejections, Rejection(selected_channels[1], m_loc_idx))
                else
                    for ch_i in bad_chs
                        push!(c_rejections, Rejection(selected_channels[ch_i], m_loc_idx))
                    end
                end
            end
            haskey(repaired_channels, g_idx) && (c_repaired[m_loc_idx] = repaired_channels[g_idx])
            haskey(skipped_channels, g_idx) && (c_skipped[m_loc_idx] = skipped_channels[g_idx])
        end

        rejection_info_per_condition[c_idx] = EpochRejectionInfo(
            "autoreject",
            EpochInfo(ep_dat.condition, ep_dat.condition_name, n_epochs(ep_dat)),
            length(c_rejections),
            median_th,
            c_rejections,
            0.0,
            nothing,
            unique_rejections(c_rejections),
            c_repaired,
            c_skipped,
        )
    end

    if reject
        for (c_idx, ep_dat) in enumerate(dat)
            reject_epochs!(ep_dat, rejection_info_per_condition[c_idx])
        end
    end

    n_kept = length(epochs_clean) + length(epochs_repaired)
    retention_pct = round(100 * n_kept / n_total_epochs; digits = 1)
    unique_repaired_count = length(unique(Iterators.flatten(values(repaired_channels))))

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

    bad_counts = [count(@view bad_matrix[:, c]) for c = 1:n_channels]
    channel_df = DataFrame(
        channel = selected_channels,
        threshold_uV = round.(th_vals; digits = 2),
        bad_epochs_count = bad_counts,
        bad_epochs_pct = round.(100 .* bad_counts ./ n_total_epochs; digits = 1),
    )

    verbose && @info "Autoreject completed: $n_kept/$n_total_epochs epochs retained ($retention_pct%)."
    return AutoRejectResult(
        thresholds,
        mode,
        Float64(consensus_perc),
        Float64(max_repair_perc),
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
end

autoreject!(dat::EpochData; kwargs...) = autoreject!([dat]; kwargs...)

"""
    autoreject(dat::AbstractVector{<:EpochData}; kwargs...) -> Tuple{Vector{EpochData}, AutoRejectResult}
    autoreject(dat::EpochData; kwargs...) -> Tuple{EpochData, AutoRejectResult}

Non-mutating version of [`autoreject!`](@ref). Returns `(cleaned_data, result)` without modifying `dat`.
"""
function autoreject(dat::AbstractVector{<:EpochData}; kwargs...)
    dat_copy = [copy(d) for d in dat]
    return dat_copy, autoreject!(dat_copy; kwargs...)
end

function autoreject(dat::EpochData; kwargs...)
    dat_copy = copy(dat)
    return dat_copy, autoreject!(dat_copy; kwargs...)
end
