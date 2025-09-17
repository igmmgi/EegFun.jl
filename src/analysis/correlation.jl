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
function correlation_matrix(dat::SingleDataFrameEeg; channel_selection::Function = channels())
    selected_channels = get_selected_channels(dat, channel_selection, include_meta = false, include_extra = false)
    if isempty(selected_channels)
        @minimal_error("No channels selected for correlation analysis")
    end
    
    return _correlation_matrix(dat.data, selected_channels)
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
function _correlation_matrix(df::DataFrame, channels::Vector{Symbol})
    # Extract data for selected channels
    data_matrix = Matrix(df[:, channels])
    
    # Calculate correlation matrix
    corr_matrix = cor(data_matrix)
    
    # Create DataFrame with channel names
    corr_df = DataFrame(corr_matrix, channels)
    corr_df[!, :channel] = channels
    select!(corr_df, :channel, channels...)
    
    return corr_df
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
function channel_joint_probability(dat::SingleDataFrameEeg; channel_selection::Function = channels(), threshold::Real = 0.5, normalize::Int = 1, discret::Int = 1000)
    selected_channels = get_selected_channels(dat, channel_selection, include_meta = false, include_extra = false)
    if isempty(selected_channels)
        @minimal_error("No channels selected for joint probability analysis")
    end
    
    return _channel_joint_probability(dat.data, selected_channels, threshold, normalize, discret)
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
function _channel_joint_probability(df::DataFrame, channels::Vector{Symbol}, threshold::Float64, normalize::Int, discret::Int)
    # Extract data for selected channels (transpose to get channels × samples)
    data_matrix = Matrix(df[:, channels])'
    
    # Calculate joint probability
    jp_values = _joint_probability(data_matrix, threshold, normalize, discret)
    
    # Create DataFrame with results
    jp_df = DataFrame(
        channel = channels,
        jp = jp_values
    )
    
    return jp_df
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
    n_channels, n_samples = size(signal)
    jp_values = zeros(Float64, n_channels)
    
    for ch in 1:n_channels
        # Get signal for this channel
        channel_signal = signal[ch, :]
        
        # Calculate probability map
        proba_map = zeros(Float64, discret)
        compute_probability!(proba_map, channel_signal, length(proba_map))
        
        # Calculate joint probability
        if normalize == 1
            jp_values[ch] = sum(proba_map .> threshold) / discret
        else
            jp_values[ch] = sum(proba_map .> threshold)
        end
    end
    
    return jp_values
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
    # Trim extreme values
    trimmed_data = _trim_extremes(data)
    
    # Calculate histogram with explicit bin edges to ensure correct number of bins
    min_val, max_val = extrema(trimmed_data)
    if min_val == max_val
        # Handle case where all values are the same
        probaMap .= 1.0 / length(probaMap)
        return probaMap
    end
    
    bin_edges = range(min_val, max_val, length = bins + 1)
    hist = fit(Histogram, trimmed_data, bin_edges)
    
    # Normalize to get probabilities
    probaMap .= hist.weights ./ sum(hist.weights)
    
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
    # Calculate robust statistics
    q25, q75 = quantile(x, [0.25, 0.75])
    iqr = q75 - q25
    
    # Define outlier bounds
    lower_bound = q25 - 1.5 * iqr
    upper_bound = q75 + 1.5 * iqr
    
    # Trim extreme values
    trimmed = copy(x)
    trimmed[trimmed .< lower_bound] .= lower_bound
    trimmed[trimmed .> upper_bound] .= upper_bound
    
    return trimmed
end
