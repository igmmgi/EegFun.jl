const _EPOCH_Z_MEASURES = (:variance, :max, :min, :abs, :range, :kurtosis)

"""
    detect_bad_epochs_automatic(dat::EpochData; z_criterion::Real=3, abs_criterion::Real=100,
                                channel_selection=channels(), interval_selection::Interval=times(),
                                z_measures::AbstractVector{Symbol}=[:variance, :max, :min, :abs, :range, :kurtosis],
                                name::String="rejection_info") -> EpochRejectionInfo

Detect bad channel-epoch pairs using across-epoch z-score thresholds (`z_criterion`, set `0` to disable)
and/or an absolute voltage threshold in μV (`abs_criterion`, set `0` to disable).

At least one of `z_criterion` or `abs_criterion` must be positive. Supported `z_measures` are
`:variance`, `:max`, `:min`, `:abs`, `:range`, and `:kurtosis`.
"""
function detect_bad_epochs_automatic(
    dat::EpochData;
    z_criterion::Real = 3,
    abs_criterion::Real = 100,
    channel_selection = channels(),
    interval_selection::Interval = times(),
    z_measures::AbstractVector{Symbol} = [:variance, :max, :min, :abs, :range, :kurtosis],
    name::String = "rejection_info",
)::EpochRejectionInfo

    @info "--------------------------------"
    @info "Condition: $(dat.condition) ($(dat.condition_name)) - Detecting bad epochs"

    z_criterion < 0 && @minimal_error("Z-criterion must be non-negative")
    abs_criterion < 0 && @minimal_error("Absolute criterion must be non-negative")
    z_criterion == 0 && abs_criterion == 0 && @minimal_error("One of z_criterion or abs_criterion must be > 0")

    invalid = [m for m in z_measures if m ∉ _EPOCH_Z_MEASURES]
    isempty(invalid) || @minimal_error("Invalid measures: $(invalid).")

    selected_channels = get_selected_channels(dat, channel_selection, include_meta = false, include_extra = false)
    @info "Selected channels: $(_print_vector(selected_channels))"
    isempty(selected_channels) && @minimal_error("No channels selected for epoch rejection")

    selected_samples = get_selected_samples(dat, interval_selection)
    @info "Selected interval: $(length(selected_samples)) samples"
    isempty(selected_samples) && @minimal_error("No samples selected for epoch rejection")

    z_viol, abs_viol = _calculate_epoch_metrics(dat, selected_channels, selected_samples, Float64(z_criterion), Float64(abs_criterion))

    rejected_info = Rejection[]
    z_lists = ntuple(_ -> Rejection[], 6)

    if z_criterion > 0
        for measure in z_measures
            mi = findfirst(==(measure), _EPOCH_Z_MEASURES)::Int
            list = z_lists[mi]
            for (ci, ch) in enumerate(selected_channels)
                for epoch_idx in z_viol[mi, ci]
                    r = Rejection(ch, epoch_idx)
                    push!(rejected_info, r)
                    push!(list, r)
                end
            end
        end
    end

    abs_rejections = if abs_criterion > 0
        abs_list = Rejection[]
        for (ci, ch) in enumerate(selected_channels)
            for epoch_idx in abs_viol[ci]
                r = Rejection(ch, epoch_idx)
                push!(rejected_info, r)
                push!(abs_list, r)
            end
        end
        abs_list
    else
        nothing
    end

    z_rejections = z_criterion > 0 ? ZScoreRejectionInfo(collect(Symbol, z_measures), z_lists...) : nothing
    info = EpochInfo(dat.condition, dat.condition_name, n_epochs(dat))

    return EpochRejectionInfo(
        name,
        info,
        length(rejected_info),
        abs_criterion,
        abs_rejections,
        z_criterion,
        z_rejections,
        unique_rejections(rejected_info),
        nothing,
        nothing,
    )
end

detect_bad_epochs_automatic(dat::AbstractVector{<:EpochData}; kwargs...) = detect_bad_epochs_automatic.(dat; kwargs...)

get_rejected(info::EpochRejectionInfo) = info.rejected
get_rejected(info::AbstractVector{<:EpochRejectionInfo}) = get_rejected.(info)

@inline function _epoch_sample_stats(col::AbstractVector{<:Real}, samples_idx, compute_z::Bool)
    v = @views col[samples_idx]
    mn, mx = extrema(v)
    ab = max(abs(mn), abs(mx))
    compute_z || return mn, mx, ab, 0.0, 0.0
    return mn, mx, ab, var(v), kurtosis(v)
end

function _calculate_epoch_metrics(
    dat::EpochData,
    selected_channels::AbstractVector{Symbol},
    selected_samples::AbstractVector{<:Integer},
    z_criterion::Float64,
    abs_criterion::Float64,
)
    n_ch = length(selected_channels)
    z_viol = [Int[] for _ = 1:6, _ = 1:n_ch]
    abs_viol = [Int[] for _ = 1:n_ch]
    n_ep = length(dat.data)
    n_ep == 0 && return z_viol, abs_viol

    first_epoch = first(dat.data)
    col_indices = [columnindex(first_epoch, ch) for ch in selected_channels]

    samples_idx = if selected_samples isa AbstractUnitRange
        selected_samples
    elseif length(selected_samples) > 1 &&
           issorted(selected_samples) &&
           (selected_samples[end] - selected_samples[1] + 1 == length(selected_samples))
        selected_samples[1]:selected_samples[end]
    else
        selected_samples
    end

    compute_z = z_criterion > 0
    check_abs = abs_criterion > 0

    Threads.@threads for ci = 1:n_ch
        col_idx = col_indices[ci]
        variances = Vector{Float64}(undef, compute_z ? n_ep : 0)
        max_vals = Vector{Float64}(undef, compute_z ? n_ep : 0)
        min_vals = Vector{Float64}(undef, compute_z ? n_ep : 0)
        abs_vals = Vector{Float64}(undef, compute_z ? n_ep : 0)
        ranges = Vector{Float64}(undef, compute_z ? n_ep : 0)
        kurtoses = Vector{Float64}(undef, compute_z ? n_ep : 0)

        for e = 1:n_ep
            mn, mx, ab, vr, kt = _epoch_sample_stats(dat.data[e][!, col_idx], samples_idx, compute_z)
            check_abs && ab > abs_criterion && push!(abs_viol[ci], e)
            if compute_z
                min_vals[e] = mn
                max_vals[e] = mx
                abs_vals[e] = ab
                ranges[e] = mx - mn
                variances[e] = vr
                kurtoses[e] = kt
            end
        end

        if compute_z
            for (mi, stat_vec) in enumerate((variances, max_vals, min_vals, abs_vals, ranges, kurtoses))
                zs = zscore(stat_vec)
                out = z_viol[mi, ci]
                @inbounds for e = 1:n_ep
                    abs(zs[e]) > z_criterion && push!(out, e)
                end
            end
        end
    end

    return z_viol, abs_viol
end
