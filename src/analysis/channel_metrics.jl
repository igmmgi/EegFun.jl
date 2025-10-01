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
    threshold::Real = 3.0,
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
