# =============================================================================
# CORRELATION ANALYSIS
# =============================================================================

"""
    correlation_matrix(dat::SingleDataFrameEeg; channel_selection::Function = channels())

Calculate the correlation matrix between EEG channels.

# Arguments
- `dat::SingleDataFrameEeg`: The EEG data object
- `channel_selection::Function`: Function that returns boolean vector for channel filtering (default: channels() - all channels)

# Returns
- `DataFrame`: Correlation matrix with channel names as both row and column names

# Examples
```julia
# Calculate correlation matrix for all channels
cm = correlation_matrix(dat)

# Calculate correlation matrix for specific channels
cm = correlation_matrix(dat, channel_selection = channels([:Fp1, :Fp2, :F3, :F4]))
```
"""
function correlation_matrix(
    dat::SingleDataFrameEeg;
    sample_selection::Function = samples(),
    channel_selection::Function = channels(),
    include_extra::Bool = false,
)::DataFrame
    selected_channels =
        eegfun.get_selected_channels(dat, channel_selection; include_meta = false, include_extra = include_extra)
    isempty(selected_channels) && @minimal_error_throw "No channels selected for correlation matrix"
    selected_samples = eegfun.get_selected_samples(dat, sample_selection)
    return _correlation_matrix(dat.data, selected_samples, selected_channels)
end

"""
    _correlation_matrix(df::DataFrame, channels::Vector{Symbol})

Internal function to calculate correlation matrix for specified channels.

# Arguments
- `df::DataFrame`: DataFrame containing the data
- `channels::Vector{Symbol}`: Vector of channel symbols to include

# Returns
- `DataFrame`: Correlation matrix with channel names as both row and column names
"""
function _correlation_matrix(
    dat::DataFrame,
    selected_samples::Vector{Int},
    selected_channels::Vector{Symbol},
)::DataFrame
    selected_data = select(dat, selected_channels)[selected_samples, :]
    df = DataFrame(cor(Matrix(selected_data)), selected_channels)
    insertcols!(df, 1, :row => selected_channels)
    return df
end

"""
    correlation_matrix_dual_selection(dat::SingleDataFrameEeg; 
                                     sample_selection::Function = samples(),
                                     channel_selection1::Function = channels(),
                                     channel_selection2::Function = channels(),
                                     include_extra_selection1::Bool = false,
                                     include_extra_selection2::Bool = true)

Calculate correlation matrix between two sets of channels using a single sample selection.

This function is particularly useful for calculating correlations between all EEG channels 
and EOG channels using the same time points.

# Arguments
- `dat::SingleDataFrameEeg`: The EEG data object
- `sample_selection::Function`: Function for sample selection (default: samples())
- `channel_selection1::Function`: Function for channel selection for first set (default: channels())
- `channel_selection2::Function`: Function for channel selection for second set (default: channels())
- `include_extra_selection1::Bool`: Whether to include extra channels for first set (default: false)
- `include_extra_selection2::Bool`: Whether to include extra channels for second set (default: true)

# Returns
- `DataFrame`: Correlation matrix with first channel set as rows and second channel set as columns

# Examples
```julia
# Calculate correlations between all channels and EOG channels
cm = correlation_matrix_dual_selection(dat, 
    sample_selection = samples(),  # All samples
    channel_selection1 = channels(),  # All EEG channels
    channel_selection2 = channels([:vEOG, :hEOG])  # EOG channels
)

# Calculate correlations between frontal channels and EOG channels
cm = correlation_matrix_dual_selection(dat,
    sample_selection = samples_not(:is_extreme_value_100),
    channel_selection1 = channels([:Fp1, :Fp2, :F3, :F4, :F5, :F6, :F7, :F8]),
    channel_selection2 = channels([:vEOG, :hEOG])
)

# Include extra channels for EOG set but not EEG set
cm = correlation_matrix_dual_selection(dat,
    sample_selection = samples(),
    channel_selection1 = channels(),  # EEG channels only
    channel_selection2 = channels([:vEOG, :hEOG]),  # EOG channels
    include_extra_selection1 = false,  # No extra channels for EEG
    include_extra_selection2 = true    # Include extra channels for EOG
)
```
"""
function correlation_matrix_dual_selection(
    dat::SingleDataFrameEeg;
    sample_selection::Function = samples(),
    channel_selection1::Function = channels(),
    channel_selection2::Function = channels(),
    include_extra_selection1::Bool = false,
    include_extra_selection2::Bool = true,
)::DataFrame
    # Get selected channels for both sets with separate include_extra settings
    selected_channels1 =
        get_selected_channels(dat, channel_selection1; include_meta = false, include_extra = include_extra_selection1)
    selected_channels2 =
        get_selected_channels(dat, channel_selection2; include_meta = false, include_extra = include_extra_selection2)

    isempty(selected_channels1) && @minimal_error_throw "No channels selected for first channel set"
    isempty(selected_channels2) && @minimal_error_throw "No channels selected for second channel set"

    # Get selected samples (same for both channel sets)
    selected_samples = get_selected_samples(dat, sample_selection)

    return _correlation_matrix_dual_selection(dat.data, selected_samples, selected_channels1, selected_channels2)
end

"""
    _correlation_matrix_dual_selection(df::DataFrame, samples::Vector{Int}, channels1::Vector{Symbol}, 
                                      channels2::Vector{Symbol})

Internal function to calculate correlation matrix between two sets of channels with the same sample selection.

# Arguments
- `df::DataFrame`: DataFrame containing the data
- `samples::Vector{Int}`: Sample indices (same for both channel sets)
- `channels1::Vector{Symbol}`: First set of channels (rows)
- `channels2::Vector{Symbol}`: Second set of channels (columns)

# Returns
- `DataFrame`: Correlation matrix with first channel set as rows and second channel set as columns
"""
function _correlation_matrix_dual_selection(
    dat::DataFrame,
    selected_samples::Vector{Int},
    selected_channels1::Vector{Symbol},
    selected_channels2::Vector{Symbol},
)::DataFrame
    # Get data for both channel sets using the same sample selection
    data1 = select(dat, selected_channels1)[selected_samples, :]
    data2 = select(dat, selected_channels2)[selected_samples, :]

    # Convert to matrices and calculate correlations
    matrix1 = Matrix(data1)  # channels1 × samples
    matrix2 = Matrix(data2)  # channels2 × samples

    # Calculate correlation matrix: channels1 × channels2
    corr_matrix = cor(matrix1, matrix2)

    # Create DataFrame with proper column and row names
    df = DataFrame(corr_matrix, selected_channels2)
    insertcols!(df, 1, :row => selected_channels1)

    return df
end

"""
    channel_joint_probability(dat::SingleDataFrameEeg; channel_selection::Function = channels(), threshold::Real = 5.0, normalize::Int = 2, discret::Int = 1000)

Calculate joint probability of extreme values across EEG channels.

# Arguments
- `dat::SingleDataFrameEeg`: The EEG data object
- `channel_selection::Function`: Function that returns boolean vector for channel filtering (default: channels() - all channels)
- `threshold::Real`: Threshold for extreme value detection (default: 5.0)
- `normalize::Int`: Normalization method (default: 2)
- `discret::Int`: Number of discretization bins (default: 1000)

# Returns
- `DataFrame`: Joint probability data with channel names and probability values

# Examples
```julia
# Calculate joint probability for all channels
jp = channel_joint_probability(dat)

# Calculate joint probability with custom threshold
jp = channel_joint_probability(dat, threshold = 0.3)
```
"""
function channel_joint_probability(
    dat::SingleDataFrameEeg;
    sample_selection::Function = samples(),
    channel_selection::Function = channels(),
    include_extra::Bool = false,
    threshold::Real = 5.0,
    normalize::Int = 2,
    discret::Int = 1000,
)::DataFrame
    selected_channels =
        get_selected_channels(dat, channel_selection; include_meta = false, include_extra = include_extra)
    isempty(selected_channels) && @minimal_error_throw "No channels selected for joint probability calculation"
    selected_samples = get_selected_samples(dat, sample_selection)

    return _channel_joint_probability(
        dat.data,
        selected_samples,
        selected_channels;
        threshold = threshold,
        normval = normalize,
        discret = discret,
    )
end

"""
    _channel_joint_probability(df::DataFrame, channels::Vector{Symbol}, threshold::Float64, normalize::Int, discret::Int)

Internal function to calculate joint probability for specified channels.

# Arguments
- `df::DataFrame`: DataFrame containing the data
- `channels::Vector{Symbol}`: Vector of channel symbols to include
- `threshold::Float64`: Threshold for extreme value detection
- `normalize::Int`: Normalization method
- `discret::Int`: Number of discretization bins

# Returns
- `DataFrame`: Joint probability data with channel names and probability values
"""
function _channel_joint_probability(
    dat::DataFrame,
    selected_samples::Vector{Int},
    selected_channels::Vector{Symbol};
    threshold::Float64 = 5.0,
    normval::Int = 2,
    discret::Int = 1000,
)::DataFrame
    @info "channel_joint_probability: Computing probability for channels $(print_vector(selected_channels))"

    # Select the specified channels and filter by samples
    data = select(dat[selected_samples, :], selected_channels)

    # Convert to matrix and compute joint probability
    jp, indelec = _joint_probability(Matrix(data)', threshold, normval)
    return DataFrame(channel = selected_channels, jp = jp, rejection = indelec)
end

"""
    _joint_probability(signal::AbstractMatrix{Float64}, threshold::Float64, normalize::Int, discret::Int = 1000)

Calculate joint probability of extreme values in a signal matrix.

# Arguments
- `signal::AbstractMatrix{Float64}`: Signal data matrix (channels × samples)
- `threshold::Float64`: Threshold for extreme value detection
- `normalize::Int`: Normalization method
- `discret::Int`: Number of discretization bins

# Returns
- `Tuple{Vector{Float64}, Vector{Bool}}`: Joint probability values for each channel and rejection flags
"""
function _joint_probability(signal::AbstractMatrix{Float64}, threshold::Float64, normalize::Int, discret::Int = 1000)
    nbchan = size(signal, 1)
    jp = zeros(nbchan)
    dataProba = Vector{Float64}(undef, size(signal, 2)) # Pre-allocate

    @inbounds for rc = 1:nbchan
        compute_probability!(dataProba, view(signal, rc, :), discret)
        jp[rc] = -sum(log, dataProba)
    end

    # Normalize the joint probability
    if normalize != 0
        tmpjp = normalize == 2 ? _trim_extremes(jp) : jp
        jp .= (jp .- mean(tmpjp)) ./ std(tmpjp)
    end

    rej = threshold != 0 ? abs.(jp) .> threshold : falses(nbchan)
    return jp, rej
end

"""
    compute_probability!(probaMap::Vector{Float64}, data::AbstractVector{Float64}, bins::Int)::Vector{Float64}

Compute probability distribution for a data vector.

# Arguments
- `probaMap::Vector{Float64}`: Pre-allocated vector to store probability values
- `data::AbstractVector{Float64}`: Input data vector
- `bins::Int`: Number of bins for discretization

# Returns
- `Vector{Float64}`: Probability distribution (same as probaMap)
"""
function compute_probability!(probaMap::Vector{Float64}, data::AbstractVector{Float64}, bins::Int)::Vector{Float64}
    isempty(data) && @minimal_error_throw "Cannot compute probability for empty data vector"

    if bins > 0
        min_val, max_val = extrema(data)
        range_val = max_val - min_val

        # Handle case where all values are the same (range_val == 0)
        if range_val == 0
            # All values are identical - uniform probability
            fill!(probaMap, 1.0 / length(data))
            return probaMap
        end

        sortbox = zeros(Int, bins)

        # Single-pass binning and counting
        @inbounds for x in data
            bin = clamp(floor(Int, (x - min_val) / range_val * (bins - 1)) + 1, 1, bins)
            sortbox[bin] += 1
        end

        # Compute probabilities
        n = length(data)
        @inbounds for (i, x) in enumerate(data)
            bin = clamp(floor(Int, (x - min_val) / range_val * (bins - 1)) + 1, 1, bins)
            probaMap[i] = sortbox[bin] / n
        end
    else
        # Gaussian approximation
        μ, σ = mean(data), std(data)

        # Handle case where all values are the same (σ == 0)
        if σ == 0
            # All values are identical - uniform probability
            fill!(probaMap, 1.0 / length(data))
            return probaMap
        end

        inv_sqrt2pi = 1 / (√(2π))
        @inbounds for (i, x) in enumerate(data)
            z = (x - μ) / σ
            probaMap[i] = exp(-0.5 * z * z) * inv_sqrt2pi
        end
        sum_p = sum(probaMap)
        probaMap ./= sum_p
    end

    return probaMap

end


"""
    _trim_extremes(x::Vector{Float64})

Trim extreme values from a data vector using robust statistics.

# Arguments
- `x::Vector{Float64}`: Input data vector

# Returns
- `Vector{Float64}`: Data vector with extreme values trimmed
"""
function _trim_extremes(x::Vector{Float64})
    n = length(x)
    trim = round(Int, n * 0.1)
    sorted = sort(x)
    return view(sorted, (trim+1):(n-trim))
end

# =============================================================================
# CORRELATION MATRIX UTILITIES
# =============================================================================

"""
    get_eog_channels(eog_cfg::EogConfig)

Extract EOG channel symbols from EogConfig.

# Arguments
- `eog_cfg::EogConfig`: The EOG configuration object

# Returns
- `Vector{Symbol}`: Vector of EOG channel symbols (e.g., [:vEOG, :hEOG])

# Examples
```julia
eog_cfg = EogConfig(50.0, 30.0, [["Fp1", "Fp2"], ["IO1", "IO2"], ["vEOG"]], [["F9"], ["F10"], ["hEOG"]])
eog_channels = get_eog_channels(eog_cfg)
# Returns [:vEOG, :hEOG]
```
"""
function get_eog_channels(eog_cfg::EogConfig)
    vEOG_channel = Symbol(eog_cfg.vEOG_channels[3][1])
    hEOG_channel = Symbol(eog_cfg.hEOG_channels[3][1])
    return [vEOG_channel, hEOG_channel]
end

"""
    add_zscore_columns!(df::DataFrame, exclude_columns::Vector{Symbol} = [:row])

Add z-score columns to a correlation matrix DataFrame.

This function adds z-score columns for each numeric column in the DataFrame,
excluding specified columns (default: :row).

# Arguments
- `df::DataFrame`: The correlation matrix DataFrame to modify
- `exclude_columns::Vector{Symbol}`: Columns to exclude from z-score calculation (default: [:row])

# Returns
- `DataFrame`: The modified DataFrame with z-score columns added

# Examples
```julia
# Add z-score columns to correlation matrix
cm = correlation_matrix_dual_selection(dat, 
    channel_selection1 = channels(),
    channel_selection2 = channels([:vEOG, :hEOG])
)
add_zscore_columns!(cm)  # Adds z_vEOG, z_hEOG columns

# Exclude additional columns from z-score calculation
add_zscore_columns!(cm, [:row, :channel_name])
```
"""
function add_zscore_columns!(df::DataFrame, exclude_columns::Vector{Symbol} = [:row])
    # Get numeric columns (excluding specified columns)
    numeric_columns = filter(col -> eltype(df[!, col]) <: Number && !(col in exclude_columns), names(df))

    if isempty(numeric_columns)
        @minimal_error "No numeric columns found for z-score calculation"
        return df
    end

    # Add z-score columns for each numeric column
    for col in numeric_columns
        z_col_name = Symbol("z_$(col)")
        values = df[!, col]

        # Calculate z-scores (handle case where all values are the same)
        if std(values) == 0
            df[!, z_col_name] = zeros(length(values))
        else
            df[!, z_col_name] = (values .- mean(values)) ./ std(values)
        end
    end

    return df
end

"""
    add_zscore_columns(df::DataFrame, exclude_columns::Vector{Symbol} = [:row])

Non-mutating version of add_zscore_columns!.

# Arguments
- `df::DataFrame`: The correlation matrix DataFrame
- `exclude_columns::Vector{Symbol}`: Columns to exclude from z-score calculation (default: [:row])

# Returns
- `DataFrame`: A new DataFrame with z-score columns added

# Examples
```julia
# Create new DataFrame with z-score columns
cm_with_z = add_zscore_columns(correlation_matrix_dual_selection(dat, 
    channel_selection1 = channels(),
    channel_selection2 = channels([:vEOG, :hEOG])
))
```
"""
function add_zscore_columns(df::DataFrame, exclude_columns::Vector{Symbol} = [:row])
    df_copy = copy(df)
    add_zscore_columns!(df_copy, exclude_columns)
    return df_copy
end

"""
    correlation_matrix_eog(dat::SingleDataFrameEeg, eog_cfg::EogConfig; 
                          sample_selection::Function = samples(),
                          channel_selection::Function = channels(),
                          include_extra::Bool = false)

Calculate correlation matrix between EEG channels and EOG channels using EogConfig.

This is a convenience function that extracts EOG channel names from the EogConfig
and calls correlation_matrix_dual_selection internally.

# Arguments
- `dat::SingleDataFrameEeg`: The EEG data object
- `eog_cfg::EogConfig`: The EOG configuration object
- `sample_selection::Function`: Function for sample selection (default: samples())
- `channel_selection::Function`: Function for channel selection for EEG channels (default: channels())
- `include_extra::Bool`: Whether to include extra channels (default: false)

# Returns
- `DataFrame`: Correlation matrix with EEG channels as rows and EOG channels as columns

# Examples
```julia
# Calculate correlations between all EEG channels and EOG channels
cm = correlation_matrix_eog(dat, eog_cfg)

# Calculate correlations with specific sample and channel selection
cm = correlation_matrix_eog(dat, eog_cfg,
    sample_selection = samples_not(:is_extreme_value_100),
    channel_selection = channels([:Fp1, :Fp2, :F3, :F4])
)

# Add z-score columns
cm = correlation_matrix_eog(dat, eog_cfg)
add_zscore_columns!(cm)
```
"""
function correlation_matrix_eog(
    dat::SingleDataFrameEeg,
    eog_cfg::EogConfig;
    sample_selection::Function = samples(),
    channel_selection::Function = channels(),
    include_extra::Bool = false,
)::DataFrame
    # Call the dual selection function
    return correlation_matrix_dual_selection(
        dat;
        sample_selection = sample_selection,
        channel_selection1 = channel_selection,
        channel_selection2 = channels(get_eog_channels(eog_cfg)),
        include_extra_selection1 = include_extra,
        include_extra_selection2 = true,  # EOG channels are typically extra channels
    )
end

# =============================================================================
# BAD CHANNEL IDENTIFICATION
# =============================================================================

"""
    identify_bad_channels(summary_df::DataFrame, joint_prob_df::DataFrame; 
                         zvar_criterion::Real = 3.0)::Vector{Symbol}

Identify bad channels based on channel summary z-variance and joint probability criteria.

# Arguments
- `summary_df::DataFrame`: Channel summary DataFrame (output from channel_summary)
- `joint_prob_df::DataFrame`: Joint probability DataFrame (output from channel_joint_probability)
- `zvar_criterion::Real`: Z-variance threshold for bad channel identification (default: 3.0)

# Returns
- `Vector{Symbol}`: Vector of channel names identified as bad channels

# Examples
```julia
# Get channel summary and joint probability
summary_df = channel_summary(dat)
joint_prob_df = channel_joint_probability(dat)

# Identify bad channels with default criteria
bad_channels = identify_bad_channels(summary_df, joint_prob_df)

# Use custom z-variance criterion
bad_channels = identify_bad_channels(summary_df, joint_prob_df, zvar_criterion = 2.5)
```
"""
function identify_bad_channels(
    summary_df::DataFrame,
    joint_prob_df::DataFrame;
    zvar_criterion::Real = 3.0,
)::Vector{Symbol}

    # Identify bad channels based on z-variance criterion
    bad_by_zvar = summary_df[abs.(summary_df.zvar).>zvar_criterion, :channel]

    # Identify bad channels based on joint probability criterion
    bad_by_jp = joint_prob_df[joint_prob_df.rejection, :channel]

    # Combine both criteria (union of bad channels)
    all_bad_channels = unique(vcat(bad_by_zvar, bad_by_jp))

    return all_bad_channels
end

"""
    partition_channels_by_eog_correlation(bad_channels::Vector{Symbol}, eog_correlation_df::DataFrame;
                                          eog_channels::Vector{Symbol} = [:hEOG, :vEOG],
                                          threshold::Real = 0.3,
                                          use_z::Bool = false)::Tuple{Vector{Symbol},Vector{Symbol}}

Partition bad channels into two groups based on correlation with EOG columns:
- First element: bad channels NOT highly correlated with EOG (retain for non-EOG handling)
- Second element: bad channels correlated with EOG (prefer ICA handling)

# Arguments
- `bad_channels::Vector{Symbol}`: Vector of bad channel names
- `eog_correlation_df::DataFrame`: EOG correlation matrix DataFrame
- `eog_channels::Vector{Symbol}`: EOG columns to use (default: `[:hEOG, :vEOG]`)
- `threshold::Real`: Correlation threshold (default: 0.3)
- `use_z::Bool`: If true, use z-scored equivalents (e.g., `:z_hEOG`, `:z_vEOG`)
"""
function partition_channels_by_eog_correlation(
    bad_channels::Vector{Symbol},
    eog_correlation_df::DataFrame;
    eog_channels::Vector{Symbol} = [:hEOG, :vEOG],
    threshold::Real = 0.3,
    use_z::Bool = false,
)::Tuple{Vector{Symbol},Vector{Symbol}}

    if isempty(bad_channels)
        return bad_channels, Symbol[]
    end

    if isempty(eog_channels)
        @minimal_warning "No EOG channels provided to partition_channels_by_eog_correlation."
        return bad_channels, Symbol[]
    end

    cols_to_use = use_z ? [Symbol("z_$(ch)") for ch in eog_channels] : eog_channels

    eog_related = Symbol[]
    for ch in bad_channels
        rows = eog_correlation_df[eog_correlation_df.row.==ch, :]
        if nrow(rows) > 0
            if any(hasproperty(rows, c) && abs(rows[1, c]) > threshold for c in cols_to_use)
                push!(eog_related, ch)
            end
        end
    end

    non_eog_bad = setdiff(bad_channels, eog_related)
    return non_eog_bad, eog_related
end

"""
    check_channel_neighbors(bad_channels::Vector{Symbol}, layout::Layout)::Vector{Symbol}

Check which bad channels can be repaired using neighbor interpolation.
A channel is repairable only if:
- ALL of its neighbors are good (not in the bad_channels list)
- There are at least 2 good neighbors (required for meaningful interpolation)
This prevents bad channels from being used to repair other bad channels.

# Arguments
- `bad_channels::Vector{Symbol}`: Vector of bad channel names
- `layout::Layout`: Layout object containing neighbor information

# Returns
- `Vector{Symbol}`: Bad channels that have ALL good neighbors (can be repaired)

# Examples
```julia
# Check which bad channels can be repaired
repairable_channels = check_channel_neighbors(bad_channels, layout)
```
"""
function check_channel_neighbors(bad_channels::Vector{Symbol}, layout::Layout)::Vector{Symbol}

    if isempty(bad_channels)
        return bad_channels
    end

    repairable_channels = Symbol[]

    for bad_ch in bad_channels
        # Check if this channel has neighbors defined in the layout
        if !isnothing(layout.neighbours) && haskey(layout.neighbours, bad_ch)
            neighbors = layout.neighbours[bad_ch]
            # Only repair if ALL neighbors are good (not in bad_channels list)
            # AND there are at least 2 good neighbors (need at least 2 for meaningful interpolation)
            good_neighbors = setdiff(neighbors.channels, bad_channels)
            if length(good_neighbors) == length(neighbors.channels) && length(good_neighbors) >= 2
                push!(repairable_channels, bad_ch)
            end
        end
    end

    return repairable_channels
end
