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
                                     include_extra1::Bool = false,
                                     include_extra2::Bool = false)

Calculate correlation matrix between two sets of channels using a single sample selection.

This function is particularly useful for calculating correlations between all EEG channels 
and EOG channels using the same time points.

# Arguments
- `dat::SingleDataFrameEeg`: The EEG data object
- `sample_selection::Function`: Function for sample selection (default: samples())
- `channel_selection1::Function`: Function for channel selection for first set (default: channels())
- `channel_selection2::Function`: Function for channel selection for second set (default: channels())
- `include_extra1::Bool`: Whether to include extra channels for first set (default: false)
- `include_extra2::Bool`: Whether to include extra channels for second set (default: false)

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
    include_extra1 = false,  # No extra channels for EEG
    include_extra2 = true    # Include extra channels for EOG
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
    selected_channels1 = get_selected_channels(dat, channel_selection1; include_meta = false, include_extra = include_extra_selection1)
    selected_channels2 = get_selected_channels(dat, channel_selection2; include_meta = false, include_extra = include_extra_selection2)
    
    isempty(selected_channels1) && @minimal_error_throw "No channels selected for first channel set"
    isempty(selected_channels2) && @minimal_error_throw "No channels selected for second channel set"
    
    # Get selected samples (same for both channel sets)
    selected_samples = get_selected_samples(dat, sample_selection)
    
    return _correlation_matrix_dual_selection(
        dat.data, 
        selected_samples, selected_channels1, selected_channels2
    )
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
    channel_joint_probability(dat::SingleDataFrameEeg; channel_selection::Function = channels(), threshold::Real = 0.5, normalize::Int = 1, discret::Int = 1000)

Calculate joint probability of extreme values across EEG channels.

# Arguments
- `dat::SingleDataFrameEeg`: The EEG data object
- `channel_selection::Function`: Function that returns boolean vector for channel filtering (default: channels() - all channels)
- `threshold::Real`: Threshold for extreme value detection (default: 0.5)
- `normalize::Int`: Normalization method (default: 1)
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
- `Vector{Float64}`: Joint probability values for each channel
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

    if bins > 0
        min_val, max_val = extrema(data)
        range_val = max_val - min_val
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
    numeric_columns = Base.filter(col -> eltype(df[!, col]) <: Number && !(col in exclude_columns), names(df))
    
    if isempty(numeric_columns)
        @warn "No numeric columns found for z-score calculation"
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
    return correlation_matrix_dual_selection(dat;
        sample_selection = sample_selection,
        channel_selection1 = channel_selection,
        channel_selection2 = channels(get_eog_channels(eog_cfg)),
        include_extra_selection1 = include_extra,
        include_extra_selection2 = true  # EOG channels are typically extra channels
    )
end
