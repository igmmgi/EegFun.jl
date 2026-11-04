@inline function _contiguous_sample_range(selected_samples::Vector{Int})
    n = length(selected_samples)
    n == 0 && return 1:0
    first_s = @inbounds selected_samples[1]
    last_s = @inbounds selected_samples[n]
    last_s - first_s + 1 == n || return nothing
    @inbounds for i = 2:n
        selected_samples[i] == first_s + i - 1 || return nothing
    end
    return first_s:last_s
end

function _extract_channel_sample_matrix(dat::DataFrame, selected_samples::Vector{Int}, selected_channels::Vector{Symbol})::Matrix{Float64}
    n_s = length(selected_samples)
    n_ch = length(selected_channels)
    X = Matrix{Float64}(undef, n_s, n_ch)
    contig = _contiguous_sample_range(selected_samples)

    for (c, ch) in enumerate(selected_channels)
        col = dat[!, ch]
        if col isa Vector{Float64}
            col_f = col::Vector{Float64}
            if !isnothing(contig)
                @inbounds copyto!(X, (c - 1) * n_s + 1, col_f, first(contig), n_s)
            else
                @inbounds @simd for s = 1:n_s
                    X[s, c] = col_f[selected_samples[s]]
                end
            end
        elseif !isnothing(contig)
            s0 = first(contig) - 1
            @inbounds for s = 1:n_s
                X[s, c] = Float64(col[s0+s])
            end
        else
            @inbounds for s = 1:n_s
                X[s, c] = Float64(col[selected_samples[s]])
            end
        end
    end
    return X
end

"""
    correlation_matrix(dat::SingleDataFrameEeg; sample_selection = samples(),
                       interval_selection::Interval = times(), channel_selection = channels(),
                       include_extra::Bool = false) -> DataFrame

Compute the Pearson correlation matrix across selected EEG channels, returned as a `DataFrame`
with a leading `:row` column of channel symbols.
"""
function correlation_matrix(
    dat::SingleDataFrameEeg;
    sample_selection = samples(),
    interval_selection::Interval = times(),
    channel_selection = channels(),
    include_extra::Bool = false,
)::DataFrame
    selected_channels = get_selected_channels(dat, channel_selection; include_meta = false, include_extra = include_extra)
    isempty(selected_channels) && @minimal_error("No channels selected for correlation matrix")
    selected_samples = get_selected_samples(dat, _combine_interval_sample(interval_selection, sample_selection))
    return _correlation_matrix(dat.data, selected_samples, selected_channels)
end

function _extract_unit_centered_matrix(dat::DataFrame, selected_samples::Vector{Int}, selected_channels::Vector{Symbol})::Matrix{Float64}
    n_s = length(selected_samples)
    n_ch = length(selected_channels)
    Z = Matrix{Float64}(undef, n_s, n_ch)
    contig = _contiguous_sample_range(selected_samples)
    inv_n = 1.0 / n_s

    for (c, ch) in enumerate(selected_channels)
        col = dat[!, ch]
        if col isa Vector{Float64} && !isnothing(contig)
            col_f = col::Vector{Float64}
            s0 = first(contig) - 1
            sum_x = 0.0
            @inbounds @simd for s = 1:n_s
                sum_x += col_f[s0+s]
            end
            μ = sum_x * inv_n
            ss = 0.0
            @inbounds @simd for s = 1:n_s
                d = col_f[s0+s] - μ
                Z[s, c] = d
                ss = muladd(d, d, ss)
            end
            if ss > 0.0
                inv_norm = 1.0 / sqrt(ss)
                @inbounds @simd for s = 1:n_s
                    Z[s, c] *= inv_norm
                end
            else
                @inbounds fill!(@view(Z[:, c]), NaN)
            end
        else
            sum_x = 0.0
            if col isa Vector{Float64}
                col_f = col::Vector{Float64}
                @inbounds @simd for s = 1:n_s
                    sum_x += col_f[selected_samples[s]]
                end
                μ = sum_x * inv_n
                ss = 0.0
                @inbounds @simd for s = 1:n_s
                    d = col_f[selected_samples[s]] - μ
                    Z[s, c] = d
                    ss = muladd(d, d, ss)
                end
            elseif !isnothing(contig)
                s0 = first(contig) - 1
                @inbounds for s = 1:n_s
                    sum_x += Float64(col[s0+s])
                end
                μ = sum_x * inv_n
                ss = 0.0
                @inbounds for s = 1:n_s
                    d = Float64(col[s0+s]) - μ
                    Z[s, c] = d
                    ss = muladd(d, d, ss)
                end
            else
                @inbounds for s = 1:n_s
                    sum_x += Float64(col[selected_samples[s]])
                end
                μ = sum_x * inv_n
                ss = 0.0
                @inbounds for s = 1:n_s
                    d = Float64(col[selected_samples[s]]) - μ
                    Z[s, c] = d
                    ss = muladd(d, d, ss)
                end
            end
            if ss > 0.0
                inv_norm = 1.0 / sqrt(ss)
                @inbounds @simd for s = 1:n_s
                    Z[s, c] *= inv_norm
                end
            else
                @inbounds fill!(@view(Z[:, c]), NaN)
            end
        end
    end
    return Z
end

function _cor_to_dataframe(R::Matrix{Float64}, row_labels::Vector{Symbol}, col_labels::Vector{Symbol})::DataFrame
    n_cols = length(col_labels)
    out_names = Vector{Symbol}(undef, n_cols + 1)
    out_cols = Vector{AbstractVector}(undef, n_cols + 1)
    out_names[1] = :row
    out_cols[1] = copy(row_labels)
    for j = 1:n_cols
        out_names[j+1] = col_labels[j]
        out_cols[j+1] = R[:, j]
    end
    return DataFrame(out_cols, out_names; copycols = false)
end

function _correlation_matrix(dat::DataFrame, selected_samples::Vector{Int}, selected_channels::Vector{Symbol})::DataFrame
    if length(selected_samples) <= 1
        X = _extract_channel_sample_matrix(dat, selected_samples, selected_channels)
        return _cor_to_dataframe(cor(X), selected_channels, selected_channels)
    end
    Z = _extract_unit_centered_matrix(dat, selected_samples, selected_channels)
    R = Z' * Z
    n_ch = length(selected_channels)
    @inbounds for j = 1:n_ch
        for i = 1:n_ch
            R[i, j] = i == j && !isnan(R[i, i]) ? 1.0 : clamp(R[i, j], -1.0, 1.0)
        end
    end
    return _cor_to_dataframe(R, selected_channels, selected_channels)
end

"""
    correlation_matrix_dual_selection(dat::SingleDataFrameEeg; sample_selection = samples(),
                                      interval_selection::Interval = times(),
                                      channel_selection1 = channels(), channel_selection2 = channels(),
                                      include_extra_selection1::Bool = false,
                                      include_extra_selection2::Bool = true) -> DataFrame

Compute the cross-correlation matrix between two channel sets (`channel_selection1` as rows,
`channel_selection2` as columns) over a shared sample selection.
"""
function correlation_matrix_dual_selection(
    dat::SingleDataFrameEeg;
    sample_selection = samples(),
    interval_selection::Interval = times(),
    channel_selection1 = channels(),
    channel_selection2 = channels(),
    include_extra_selection1::Bool = false,
    include_extra_selection2::Bool = true,
)::DataFrame
    selected_channels1 = get_selected_channels(dat, channel_selection1; include_meta = false, include_extra = include_extra_selection1)
    selected_channels2 = get_selected_channels(dat, channel_selection2; include_meta = false, include_extra = include_extra_selection2)
    isempty(selected_channels1) && @minimal_error("No channels selected for first channel set")
    isempty(selected_channels2) && @minimal_error("No channels selected for second channel set")

    selected_samples = get_selected_samples(dat, _combine_interval_sample(interval_selection, sample_selection))
    return _correlation_matrix_dual_selection(dat.data, selected_samples, selected_channels1, selected_channels2)
end

function _correlation_matrix_dual_selection(
    dat::DataFrame,
    selected_samples::Vector{Int},
    selected_channels1::Vector{Symbol},
    selected_channels2::Vector{Symbol},
)::DataFrame
    if length(selected_samples) <= 1
        matrix1 = _extract_channel_sample_matrix(dat, selected_samples, selected_channels1)
        matrix2 = _extract_channel_sample_matrix(dat, selected_samples, selected_channels2)
        return _cor_to_dataframe(cor(matrix1, matrix2), selected_channels1, selected_channels2)
    end
    Z1 = _extract_unit_centered_matrix(dat, selected_samples, selected_channels1)
    Z2 = _extract_unit_centered_matrix(dat, selected_samples, selected_channels2)
    R = Z1' * Z2
    @inbounds for k in eachindex(R)
        R[k] = clamp(R[k], -1.0, 1.0)
    end
    return _cor_to_dataframe(R, selected_channels1, selected_channels2)
end

"""
    channel_joint_probability(dat::SingleDataFrameEeg; sample_selection = samples(),
                              interval_selection::Interval = times(), channel_selection = channels(),
                              include_extra::Bool = false, threshold::Real = 5.0,
                              normalize::Int = 2, discret::Int = 1000) -> DataFrame

Compute per-channel joint probability (negative log-likelihood) and threshold-based rejection flags.
"""
function channel_joint_probability(
    dat::SingleDataFrameEeg;
    sample_selection = samples(),
    interval_selection::Interval = times(),
    channel_selection = channels(),
    include_extra::Bool = false,
    threshold::Real = 5.0,
    normalize::Int = 2,
    discret::Int = 1000,
)::DataFrame
    selected_channels = get_selected_channels(dat, channel_selection; include_meta = false, include_extra = include_extra)
    isempty(selected_channels) && @minimal_error("No channels selected for joint probability calculation")
    selected_samples = get_selected_samples(dat, _combine_interval_sample(interval_selection, sample_selection))
    return _channel_joint_probability(
        dat.data,
        selected_samples,
        selected_channels;
        threshold = threshold,
        normval = normalize,
        discret = discret,
    )
end

function _vec_joint_probability_value(data::AbstractVector{<:Real}, bins::Int, sortbox::Vector{Int})::Float64
    n = length(data)
    n == 0 && @minimal_error("Cannot compute probability for empty data vector")

    if bins > 0
        min_v = Float64(@inbounds data[1])
        max_v = min_v
        @inbounds for i = 2:n
            v = Float64(data[i])
            v < min_v && (min_v = v)
            v > max_v && (max_v = v)
        end
        range_val = max_v - min_v
        range_val == 0.0 && return n * log(Float64(n))

        fill!(sortbox, 0)
        bm1 = bins - 1
        @inbounds for i = 1:n
            bin = clamp(floor(Int, (Float64(data[i]) - min_v) / range_val * bm1) + 1, 1, bins)
            sortbox[bin] += 1
        end

        inv_n = 1.0 / n
        jp_val = 0.0
        @inbounds for b = 1:bins
            c = sortbox[b]
            c > 0 && (jp_val -= c * log(c * inv_n))
        end
        return jp_val
    else
        μ, σ = mean(data), std(data)
        σ == 0 && return n * log(Float64(n))
        inv_sqrt2pi = 1.0 / sqrt(2π)
        sum_p = 0.0
        sum_z2 = 0.0
        @inbounds for i = 1:n
            z = (Float64(data[i]) - μ) / σ
            z2 = z * z
            sum_z2 += z2
            sum_p += exp(-0.5 * z2) * inv_sqrt2pi
        end
        return 0.5 * sum_z2 - n * log(inv_sqrt2pi / sum_p)
    end
end

function _normalize_and_threshold_jp!(jp::Vector{Float64}, threshold::Real, normalize::Int)
    if normalize != 0
        tmpjp = normalize == 2 ? _trim_extremes(jp) : jp
        s = isempty(tmpjp) ? 0.0 : std(tmpjp)
        if !isfinite(s) || s == 0
            fill!(jp, 0.0)
        else
            μ = mean(tmpjp)
            jp .= (jp .- μ) ./ s
        end
    end
    rej = threshold != 0 ? (jp .> threshold) : falses(length(jp))
    return jp, rej
end

function _channel_joint_probability(
    dat::DataFrame,
    selected_samples::Vector{Int},
    selected_channels::Vector{Symbol};
    threshold::Real = 5.0,
    normval::Int = 2,
    discret::Int = 1000,
)
    @info "channel_joint_probability: Computing probability for channels $(_print_vector(selected_channels))"
    n_ch = length(selected_channels)
    n_s = length(selected_samples)
    contig = _contiguous_sample_range(selected_samples)
    jp = Vector{Float64}(undef, n_ch)

    _foreach_chunk(n_ch) do ch_range
        sortbox = discret > 0 ? Vector{Int}(undef, discret) : Int[]
        sample_buf = isnothing(contig) ? Vector{Float64}(undef, n_s) : Float64[]
        for c in ch_range
            col = dat[!, selected_channels[c]]
            if !isnothing(contig)
                jp[c] = _vec_joint_probability_value(@view(col[contig]), discret, sortbox)
            else
                if col isa Vector{Float64}
                    col_f = col::Vector{Float64}
                    @inbounds @simd for i = 1:n_s
                        sample_buf[i] = col_f[selected_samples[i]]
                    end
                else
                    @inbounds for i = 1:n_s
                        sample_buf[i] = Float64(col[selected_samples[i]])
                    end
                end
                jp[c] = _vec_joint_probability_value(sample_buf, discret, sortbox)
            end
        end
    end

    jp, indelec = _normalize_and_threshold_jp!(jp, threshold, normval)
    return DataFrame(channel = copy(selected_channels), jp = jp, rejection = indelec; copycols = false)
end

function _joint_probability(signal::AbstractMatrix{<:Real}, threshold::Real, normalize::Int, discret::Int = 1000)
    nbchan = size(signal, 1)
    jp = Vector{Float64}(undef, nbchan)
    sortbox = discret > 0 ? Vector{Int}(undef, discret) : Int[]

    @inbounds for rc = 1:nbchan
        jp[rc] = _vec_joint_probability_value(@view(signal[rc, :]), discret, sortbox)
    end
    return _normalize_and_threshold_jp!(jp, threshold, normalize)
end

function compute_probability!(probaMap::Vector{Float64}, data::AbstractVector{<:Real}, bins::Int)::Vector{Float64}
    isempty(data) && @minimal_error("Cannot compute probability for empty data vector")
    n = length(data)

    if bins > 0
        min_val, max_val = extrema(data)
        range_val = max_val - min_val
        if range_val == 0
            fill!(probaMap, 1.0 / n)
            return probaMap
        end

        sortbox = zeros(Int, bins)
        bm1 = bins - 1
        @inbounds for x in data
            bin = clamp(floor(Int, (x - min_val) / range_val * bm1) + 1, 1, bins)
            sortbox[bin] += 1
        end
        inv_n = 1.0 / n
        @inbounds for (i, x) in enumerate(data)
            bin = clamp(floor(Int, (x - min_val) / range_val * bm1) + 1, 1, bins)
            probaMap[i] = sortbox[bin] * inv_n
        end
    else
        μ, σ = mean(data), std(data)
        if σ == 0
            fill!(probaMap, 1.0 / n)
            return probaMap
        end
        inv_sqrt2pi = 1.0 / sqrt(2π)
        sum_p = 0.0
        @inbounds for (i, x) in enumerate(data)
            z = (x - μ) / σ
            p = exp(-0.5 * z * z) * inv_sqrt2pi
            probaMap[i] = p
            sum_p += p
        end
        probaMap ./= sum_p
    end
    return probaMap
end

function _trim_extremes(x::Vector{Float64})
    n = length(x)
    trim = round(Int, n * 0.1)
    sorted = sort(x)
    lo = trim + 1
    hi = n - trim
    return lo <= hi ? view(sorted, lo:hi) : view(sorted, 1:n)
end

"""
    get_eog_channels(eog_cfg::EogConfig) -> Vector{Symbol}

Return `[vEOG_output_symbol, hEOG_output_symbol]` from `eog_cfg`.
"""
get_eog_channels(eog_cfg::EogConfig) = [Symbol(eog_cfg.vEOG_channels[3][1]), Symbol(eog_cfg.hEOG_channels[3][1])]

function add_zscore_columns!(df::DataFrame, exclude_columns::Vector{Symbol} = [:row])
    numeric_columns = [col for col in propertynames(df) if eltype(df[!, col]) <: Number && !(col in exclude_columns)]
    isempty(numeric_columns) && @minimal_error("No numeric columns found for z-score calculation")

    for col in numeric_columns
        values = df[!, col]
        s = std(values)
        df[!, Symbol("z_$(col)")] = s == 0 ? zeros(length(values)) : (values .- mean(values)) ./ s
    end
    return df
end

add_zscore_columns(df::DataFrame, exclude_columns::Vector{Symbol} = [:row]) = add_zscore_columns!(copy(df), exclude_columns)

"""
    correlation_matrix_eog(dat::SingleDataFrameEeg, eog_cfg::EogConfig; sample_selection = samples(),
                           interval_selection::Interval = times(), channel_selection = channels(),
                           include_extra::Bool = false) -> DataFrame

Compute the correlation matrix between EEG channels (`channel_selection`) and the EOG channels configured in `eog_cfg`.
"""
function correlation_matrix_eog(
    dat::SingleDataFrameEeg,
    eog_cfg::EogConfig;
    sample_selection = samples(),
    interval_selection::Interval = times(),
    channel_selection = channels(),
    include_extra::Bool = false,
)::DataFrame
    return correlation_matrix_dual_selection(
        dat;
        sample_selection = sample_selection,
        interval_selection = interval_selection,
        channel_selection1 = channel_selection,
        channel_selection2 = channels(get_eog_channels(eog_cfg)),
        include_extra_selection1 = include_extra,
        include_extra_selection2 = true,
    )
end

"""
    identify_bad_channels(summary_df::DataFrame, joint_prob_df::DataFrame; zvar_criterion::Real = 3.0) -> Vector{Symbol}

Return the union of channels flagged by `abs(zvar) > zvar_criterion` in `summary_df` or `rejection == true` in `joint_prob_df`.
"""
function identify_bad_channels(summary_df::DataFrame, joint_prob_df::DataFrame; zvar_criterion::Real = 3.0)::Vector{Symbol}
    bad_by_zvar = summary_df[abs.(summary_df.zvar) .> zvar_criterion, :channel]
    bad_by_jp = joint_prob_df[joint_prob_df.rejection, :channel]
    return unique(vcat(bad_by_zvar, bad_by_jp))
end

"""
    partition_channels_by_eog_correlation(bad_channels::Vector{Symbol}, eog_correlation_df::DataFrame;
                                          eog_channels::Vector{Symbol} = [:hEOG, :vEOG],
                                          threshold::Real = 0.3, use_z::Bool = false) -> Tuple{Vector{Symbol},Vector{Symbol}}

Split `bad_channels` into `(non_eog_bad, eog_related)` according to whether any correlation in `eog_correlation_df`
exceeds `threshold` in magnitude.
"""
function partition_channels_by_eog_correlation(
    bad_channels::Vector{Symbol},
    eog_correlation_df::DataFrame;
    eog_channels::Vector{Symbol} = [:hEOG, :vEOG],
    threshold::Real = 0.3,
    use_z::Bool = false,
)::Tuple{Vector{Symbol},Vector{Symbol}}
    isempty(bad_channels) && return bad_channels, Symbol[]
    if isempty(eog_channels)
        @minimal_warning "No EOG channels provided to partition_channels_by_eog_correlation."
        return bad_channels, Symbol[]
    end

    cols_to_use = [c for c in (use_z ? [Symbol("z_$(ch)") for ch in eog_channels] : eog_channels) if hasproperty(eog_correlation_df, c)]
    row_col = eog_correlation_df.row
    eog_related = Symbol[]
    for ch in bad_channels
        r_idx = findfirst(==(ch), row_col)
        if !isnothing(r_idx) && any(c -> abs(eog_correlation_df[r_idx, c]) > threshold, cols_to_use)
            push!(eog_related, ch)
        end
    end
    return setdiff(bad_channels, eog_related), eog_related
end

"""
    check_channel_neighbors(bad_channels::Vector{Symbol}, layout::Layout) -> Vector{Symbol}

Return channels in `bad_channels` that have at least 2 spatial neighbors in `layout` and none of their neighbors in `bad_channels`.
"""
function check_channel_neighbors(bad_channels::Vector{Symbol}, layout::Layout)::Vector{Symbol}
    (isempty(bad_channels) || isnothing(layout.neighbors)) && return Symbol[]
    bad_set = Set(bad_channels)
    repairable = Symbol[]
    for bad_ch in bad_channels
        if haskey(layout.neighbors, bad_ch)
            nbrs = layout.neighbors[bad_ch].channels
            if length(nbrs) >= 2 && all(ch -> !(ch in bad_set), nbrs)
                push!(repairable, bad_ch)
            end
        end
    end
    return repairable
end
