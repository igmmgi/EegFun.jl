



"""
    create_eeg_dataframe(data::BioSemiBDF.BioSemiData)::DataFrame

Creates a DataFrame containing EEG data from a BioSemiBDF object.

# Arguments
- `data::BioSemiBDF.BioSemiData`: The BioSemi data structure containing time, triggers, and channel data.

# Returns
A DataFrame with columns: `time`, `sample`, `triggers` (cleaned), and channel data from the BioSemiBDF.

# Note
The trigger data is automatically cleaned to detect only onset events, converting sustained trigger 
signals into single onset events. For example: [0, 1, 1, 0, 0, 2, 2, 2, 0, 0] becomes [0, 1, 0, 0, 0, 2, 0, 0, 0, 0].
"""
function create_eeg_dataframe(dat::BioSemiBDF.BioSemiData)::DataFrame
    @info "create_eeg_dataframe: Creating EEG DataFrame"
    df = hcat(
        DataFrame(file = filename(dat), time = dat.time, sample = 1:length(dat.time), triggers = _clean_triggers(dat.triggers.raw)),
        DataFrame(Float64.(dat.data), Symbol.(dat.header.channel_labels[1:(end-1)])),  # assumes last channel is trigger
    )
    return df
end


"""
    create_eeg_dataframe(dat::BioSemiBDF.BioSemiData, layout::DataFrame)::ContinuousData

Creates a ContinuousData object from a BioSemiBDF data structure and a layout DataFrame.

# Arguments
- `dat::BioSemiBDF.BioSemiData`: The BioSemi data structure containing EEG data.
- `layout::DataFrame`: The DataFrame containing layout information.

# Returns
A ContinuousData object containing the EEG data and layout information.

"""
function create_eeg_dataframe(dat::BioSemiBDF.BioSemiData, layout::Layout)::ContinuousData
    return ContinuousData(create_eeg_dataframe(dat), layout, dat.header.sample_rate[1], AnalysisInfo())
end






"""
    correlation_matrix(dat::ContinuousData; sample_selection::Function = samples(), channel_selection::Function = channels(), include_extra::Bool = false)::Matrix{Float64}

Calculates the correlation matrix for the EEG data.

# Arguments
- `dat::ContinuousData`: The ContinuousData object containing EEG data.
- `sample_selection::Function`: Function that returns boolean vector for sample filtering (default: include all samples).
- `channel_selection::Function`: Function that returns boolean vector for channel filtering (default: include all channels).
- `include_extra::Bool`: Whether to include additional channels (default: false).

# Returns
A matrix containing the correlation values between the specified channels.

# Examples

## Basic Usage
```julia
# Correlation matrix for layout channels only (default)
corr_matrix = correlation_matrix(dat)

# Correlation matrix for specific layout channels
corr_matrix = correlation_matrix(dat, channel_selection = channels([:Fp1, :Fp2, :F3, :F4]))

# Correlation matrix excluding reference channels from layout
corr_matrix = correlation_matrix(dat, channel_selection = channels_not([:M1, :M2]))
```

## Including Additional Channels
```julia
# Correlation matrix for additional channels (EOG, extreme value flags, etc.)
corr_matrix = correlation_matrix(dat, channel_selection = channels([:Fp1, :Fp2, :vEOG, :hEOG]))
```

## Sample Selection
```julia
# Correlation matrix only for good samples
corr_matrix = correlation_matrix(dat, sample_selection = samples_not(:is_extreme_value_100))

# Correlation matrix only within epoch windows
corr_matrix = correlation_matrix(dat, sample_selection = samples(:epoch_window))

# Correlation matrix for good samples within epochs
corr_matrix = correlation_matrix(dat, 
    sample_selection = samples_and([
        :epoch_window,
        samples_not(:is_extreme_value_100)
    ])
)
```

## Channel Selection
```julia
# Frontal channels only
corr_matrix = correlation_matrix(dat, channel_selection = channels(1:10))

# Channels starting with "F" (frontal)
corr_matrix = correlation_matrix(dat, channel_selection = x -> startswith.(string.(x), "F"))

# Parietal channels only
corr_matrix = correlation_matrix(dat, channel_selection = x -> startswith.(string.(x), "P"))

# Mix of layout and additional channels
corr_matrix = correlation_matrix(dat, channel_selection = channels([:Fp1, :Fp2, :vEOG, :hEOG]))
```

## Combined Selection
```julia
# Exclude reference channels and bad samples
corr_matrix = correlation_matrix(dat, 
    channel_selection = channels_not([:M1, :M2]),
    sample_selection = samples_not(:is_extreme_value_100)
)

# Only frontal channels, good samples
corr_matrix = correlation_matrix(dat, 
    channel_selection = channels(1:10),
    sample_selection = samples_and([
        :epoch_window, 
        samples_not(:is_extreme_value_100)
    ])
)

# Complex filtering: frontal channels, good samples, within epochs
corr_matrix = correlation_matrix(dat, 
    channel_selection = channels([:Fp1, :Fp2, :F3, :F4, :F5, :F6, :F7, :F8]),
    sample_selection = samples_and([
        :epoch_window, 
        samples_not(:is_extreme_value_100),
        samples_not(:is_vEOG),
        samples_not(:is_hEOG)
    ])
)
```

## Quality Control Applications
```julia
# 1. Detect extreme values
is_extreme_value!(dat, 100, channel_out = :is_extreme_100)

# 2. Get correlation matrix for good data
good_corr = correlation_matrix(dat, sample_selection = samples_not(:is_extreme_100))
```
"""
function correlation_matrix(
    dat::ContinuousData;
    sample_selection::Function = samples(),
    channel_selection::Function = channels(),
    include_extra::Bool = false,
)::DataFrame
    selected_channels = get_selected_channels(dat, channel_selection; include_extra = include_extra)
    selected_samples = get_selected_samples(dat, sample_selection)
    return _correlation_matrix(dat.data, selected_samples, selected_channels)
end

# Internal function for plain DataFrames with explicit channel specification
function _correlation_matrix(
    dat::DataFrame,
    selected_samples::Vector{Int},
    selected_channels::Vector{Symbol},
)::DataFrame
    # Select the specified channels and samples
    selected_data = select(dat, selected_channels)[selected_samples, :]
    # Compute correlation matrix
    df = DataFrame(cor(Matrix(selected_data)), selected_channels)
    insertcols!(df, 1, :row => selected_channels)
    return df
end


"""
    detect_eog_onsets!(dat::ContinuousData, criterion::Float64, channel_in::Symbol, channel_out::Symbol)

Detects EOG (electrooculogram) onsets in the EEG data based on a specified criterion.

# Arguments
- `dat::ContinuousData`: The ContinuousData object containing EEG data.
- `criterion::Real`: The threshold for detecting EOG onsets.
- `channel_in::Symbol`: The channel from which to detect EOG onsets.
- `channel_out::Symbol`: The channel where the detected EOG onsets will be recorded as boolean values, indicating the presence of an EOG event.

# Returns
Nothing. The function modifies the input data in place.

# Examples

## Basic Usage
```julia
# Detect vertical EOG onsets
detect_eog_onsets!(dat, 50.0, :vEOG, :is_vEOG)

# Detect horizontal EOG onsets
detect_eog_onsets!(dat, 30.0, :hEOG, :is_hEOG)

# Combine for any EOG artifact
combine_boolean_columns!(dat, [:is_vEOG, :is_hEOG], :or, output_column = :is_any_EOG)

# Create quality flags
combine_boolean_columns!(dat, [:is_extreme_value, :is_any_EOG], :nor, output_column = :is_good_data)

# Complex quality control (good samples = not extreme AND not any EOG)
combine_boolean_columns!(dat, [:is_extreme_value, :is_vEOG, :is_hEOG], :nor, output_column = :is_clean_data)

## Multiple EOG Channels
```julia
# Detect both vertical and horizontal EOG
detect_eog_onsets!(dat, 50.0, :vEOG, :is_vEOG)
detect_eog_onsets!(dat, 30.0, :hEOG, :is_hEOG)

# Combine for any EOG artifact
# dat.data[!, :is_any_EOG] = dat.data[!, :is_vEOG] .| dat.data[!, :is_hEOG]
```
"""
function detect_eog_onsets!(dat::ContinuousData, criterion::Real, channel_in::Symbol, channel_out::Symbol)
    @info "detect_eog_onsets!: Detecting EOG onsets in channel $(channel_in) with stepsize criterion $(criterion)"
    step_size = div(dat.sample_rate, 20)
    eog_signal = dat.data[1:step_size:end, channel_in]
    eog_diff = diff(eog_signal)
    eog_idx = findall(x -> abs(x) >= criterion, eog_diff)
    eog_idx = [idx for (i, idx) in enumerate(eog_idx) if i == 1 || (idx - eog_idx[i-1] > 2)] * step_size
    dat.data[!, channel_out] .= false
    dat.data[eog_idx, channel_out] .= true
    return nothing
end

# Internal function for plain DataFrames with explicit channel specification
function _is_extreme_value(
    dat::DataFrame,
    criterion::Number,
    selected_channels::Vector{Symbol},
    selected_samples::Vector{Int},
)::Vector{Bool}
    # Initialize result vector with false for all samples
    result = fill(false, nrow(dat))

    # Check for extreme values only in selected samples
    if !isempty(selected_samples)
        data_subset = select(dat[selected_samples, :], selected_channels)
        extreme_mask = any(x -> abs.(x) >= criterion, Matrix(data_subset), dims = 2)[:]
        result[selected_samples] = extreme_mask
    end

    return result
end

# Internal function for plain DataFrames with explicit channel specification
function _is_extreme_value!(
    dat::DataFrame,
    criterion::Number,
    selected_channels::Vector{Symbol};
    channel_out::Symbol = :is_extreme_value,
)
    dat[!, channel_out] .= any(x -> abs.(x) >= criterion, Matrix(select(dat, selected_channels)), dims = 2)[:]
end

"""
    is_extreme_value(dat::ContinuousData, criterion::Number; sample_selection::Function = samples(), channel_selection::Function = channels(), include_additional_channels::Bool = false)::Vector{Bool}

Checks if any values in the specified channels exceed a given criterion.

# Arguments
- `dat::ContinuousData`: The ContinuousData object containing EEG data.
- `criterion::Number`: The threshold for determining extreme values.
- `sample_selection::Function`: Function that returns boolean vector for sample filtering (default: include all samples).
- `channel_selection::Function`: Function that returns boolean vector for channel filtering (default: include all channels).
- `include_additional_channels::Bool`: Whether to include additional channels (default: false).

# Returns
A Boolean vector indicating whether any extreme values were found for each row. Only samples selected by sample_selection are checked for extreme values.

# Examples
```julia
# Check extreme values in layout channels only (default)
is_extreme_value(dat, 100)

# Check extreme values in specific layout channels
is_extreme_value(dat, 100, channel_selection = channels([:Fp1, :Fp2]))

# Check extreme values only in good samples (exclude already detected extreme values)
is_extreme_value(dat, 100, sample_selection = samples_not(:is_extreme_value_100))

# Check extreme values only within epoch windows
is_extreme_value(dat, 100, sample_selection = samples(:epoch_window))

# Check extreme values in additional channels (automatically detected)
is_extreme_value(dat, 100, channel_selection = channels([:Fp1, :vEOG, :is_extreme_value500]))

# Exclude reference channels from layout
is_extreme_value(dat, 100, channel_selection = channels_not([:M1, :M2]))

# Combined filtering: check specific channels only in good samples
is_extreme_value(dat, 100, 
    sample_selection = samples_and([:epoch_window, samples_not(:is_vEOG)]),
    channel_selection = channels([:Fp1, :Fp2, :F3, :F4])
)
```
"""
function is_extreme_value(
    dat::ContinuousData,
    criterion::Number;
    sample_selection::Function = samples(),
    channel_selection::Function = channels(),
    include_additional_channels::Bool = false,
)::Vector{Bool}

    selected_channels =
        get_selected_channels(dat, channel_selection; include_additional_channels = include_additional_channels)
    selected_samples = get_selected_samples(dat, sample_selection)

    @info "is_extreme_value: Checking for extreme values in channel $(print_vector(selected_channels)) with criterion $(criterion)"
    return _is_extreme_value(dat.data, criterion, selected_channels, selected_samples)
end

"""
    is_extreme_value!(dat::ContinuousData, criterion::Number; sample_selection::Function = samples(), channel_selection::Function = channels(), include_additional_channels::Bool = false, channel_out::Symbol = :is_extreme_value)

Checks if any values in the specified channels exceed a given criterion and adds the result as a new column.

# Arguments
- `dat::ContinuousData`: The ContinuousData object containing EEG data.
- `criterion::Number`: The threshold for determining extreme values.
- `sample_selection::Function`: Function that returns boolean vector for sample filtering (default: include all samples).
- `channel_selection::Function`: Function that returns boolean vector for channel filtering (default: include all channels).
- `include_additional_channels::Bool`: Whether to include additional channels (default: false).
- `channel_out::Symbol`: Name of the output column (default: :is_extreme_value).

# Returns
Nothing. The function modifies the input data in place.

# Examples
```julia
# Check extreme values in layout channels only (default)
is_extreme_value!(dat, 100)

# Check extreme values in specific layout channels
is_extreme_value!(dat, 100, channel_selection = channels([:Fp1, :Fp2]))

# Check extreme values only in good samples (exclude already detected extreme values)
is_extreme_value!(dat, 100, sample_selection = samples_not(:is_extreme_value_100))

# Check extreme values only within epoch windows
is_extreme_value!(dat, 100, sample_selection = samples(:epoch_window))

# Check extreme values in additional channels (automatically detected)
is_extreme_value!(dat, 100, channel_selection = channels([:Fp1, :vEOG, :is_extreme_value500]))

# Exclude reference channels from layout
is_extreme_value!(dat, 100, channel_selection = channels_not([:M1, :M2]))

# Combined filtering: check specific channels only in good samples
is_extreme_value!(dat, 100, 
    sample_selection = samples_and([:epoch_window, samples_not(:is_vEOG)]),
    channel_selection = channels([:Fp1, :Fp2, :F3, :F4])
)
```
"""
function is_extreme_value!(
    dat::ContinuousData,
    criterion::Number;
    sample_selection::Function = samples(),
    channel_selection::Function = channels(),
    include_meta::Bool = false,
    include_extra::Bool = false,
    channel_out::Symbol = :is_extreme_value,
)
    selected_channels =
        get_selected_channels(dat, channel_selection; include_meta = include_meta, include_extra = include_extra)
    selected_samples = get_selected_samples(dat, sample_selection)

    @info "is_extreme_value!: Checking for extreme values in channel $(print_vector(selected_channels)) with criterion $(criterion)"
    dat.data[!, channel_out] = _is_extreme_value(dat.data, criterion, selected_channels, selected_samples)
end



"""
    n_extreme_value(dat::ContinuousData, criterion::Number; sample_selection::Function = samples(), channel_selection::Function = channels(), include_additional_channels::Bool = false)::Int

Counts the number of extreme values in the specified channels.

# Arguments
- `dat::ContinuousData`: The ContinuousData object containing EEG data.
- `criterion::Number`: The threshold for determining extreme values.
- `sample_selection::Function`: Function that returns boolean vector for sample filtering (default: include all samples).
- `channel_selection::Function`: Function that returns boolean vector for channel filtering (default: include all channels).
- `include_additional_channels::Bool`: Whether to include additional channels (default: false).

# Returns
An integer count of the number of extreme values found in the selected samples and channels.

# Examples
```julia
# Count extreme values in layout channels only (default)
n_extreme_value(dat, 100)

# Count extreme values in specific layout channels
n_extreme_value(dat, 100, channel_selection = channels([:Fp1, :Fp2]))

# Count extreme values only in good samples (exclude already detected extreme values)
n_extreme_value(dat, 100, sample_selection = samples_not(:is_extreme_value_100))

# Count extreme values only within epoch windows
n_extreme_value(dat, 100, sample_selection = samples(:epoch_window))

# Count extreme values excluding reference channels from layout
n_extreme_value(dat, 100, channel_selection = channels_not([:M1, :M2]))

# Count extreme values in additional channels (EOG, extreme value flags, etc.)
n_extreme_value(dat, 100, channel_selection = channels([:Fp1, :vEOG, :is_extreme_value500]))

# Combined filtering: count extreme values in specific channels only in good samples
n_extreme_value(dat, 100, 
    sample_selection = samples_and([:epoch_window, samples_not(:is_vEOG)]),
    channel_selection = channels([:Fp1, :Fp2, :F3, :F4])
)
```
"""
function n_extreme_value(
    dat::ContinuousData,
    criterion::Number;
    sample_selection::Function = samples(),
    channel_selection::Function = channels(),
    include_additional_channels::Bool = false,
)::Int
    selected_channels =
        get_selected_channels(dat, channel_selection; include_additional_channels = include_additional_channels)
    selected_samples = get_selected_samples(dat, sample_selection)

    return _n_extreme_value(dat.data, criterion, selected_channels, selected_samples)
end

# Internal function for plain DataFrames with explicit channel specification
function _n_extreme_value(
    dat::DataFrame,
    criterion::Number,
    selected_channels::Vector{Symbol},
    selected_samples::Vector{Int},
)::Int
    @info "n_extreme_value: Counting extreme values in channel $(_print_vector(selected_channels)) with criterion $(criterion)"

    # Only count extreme values in selected samples
    if isempty(selected_samples)
        return 0
    end

    data_subset = select(dat[selected_samples, :], selected_channels)
    return sum(sum.(eachcol(abs.(data_subset) .>= criterion)))
end

# Internal function for plain DataFrames with explicit channel specification
function _channel_joint_probability(
    dat::DataFrame,
    selected_samples::Vector{Int},
    selected_channels::Vector{Symbol};
    threshold::Float64 = 5.0,
    normval::Int = 2,
)::DataFrame
    @info "channel_joint_probability: Computing probability for channels $(_print_vector(selected_channels))"

    # Select the specified channels and filter by samples
    data = select(dat[selected_samples, :], selected_channels)

    # Convert to matrix and compute joint probability
    jp, indelec = _joint_probability(Matrix(data)', threshold, normval)
    return DataFrame(channel = selected_channels, jp = jp, rejection = indelec)
end

"""
    channel_joint_probability(dat::ContinuousData; sample_selection::Function = samples(), channel_selection::Function = channels(), include_additional_channels::Bool = false, threshold::Real = 3.0, normval::Real = 2)::DataFrame

Computes joint probability for EEG channels.

# Arguments
- `dat::ContinuousData`: The ContinuousData object containing EEG data.
- `sample_selection::Function`: Function that returns boolean vector for sample filtering (default: include all samples).
- `channel_selection::Function`: Function that returns boolean vector for channel filtering (default: include all channels).
- `include_additional_channels::Bool`: Whether to include additional channels (default: false).
- `threshold::Real`: Threshold for joint probability (default: 3.0).
- `normval::Real`: Normalization value (default: 2).

# Returns
A DataFrame containing joint probability values for each channel.

# Examples
```julia
# Basic joint probability for layout channels only (default)
channel_joint_probability(dat)

# Filter samples where epoch_window is true
channel_joint_probability(dat, sample_selection = samples(:epoch_window))

# Filter to specific layout channels
channel_joint_probability(dat, channel_selection = channels([:Fp1, :Fp2]))

# Filter to additional channels (EOG, extreme value flags, etc.)
channel_joint_probability(dat, channel_selection = channels([:Fp1, :vEOG, :is_extreme_value500]))

# Combine both filters
channel_joint_probability(dat, 
    sample_selection = samples(:epoch_window),
    channel_selection = channels_not([:M1, :M2])
)
```
"""
function channel_joint_probability(
    dat::ContinuousData;
    sample_selection::Function = samples(),
    channel_selection::Function = channels(),
    include_additional_channels::Bool = false,
    threshold::Real = 3.0,
    normval::Real = 2,
)::DataFrame
    selected_channels =
        get_selected_channels(dat, channel_selection; include_additional_channels = include_additional_channels)
    selected_samples = get_selected_samples(dat, sample_selection)

    return _channel_joint_probability(
        dat.data,
        selected_samples,
        selected_channels;
        threshold = threshold,
        normval = normval,
    )
end


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

Computes the probability of each value in the data vector.

# Arguments
- `probaMap::Vector{Float64}`: The vector to store the computed probabilities.
- `data::AbstractVector{Float64}`: The data vector to compute the probabilities for.
- `bins::Int`: The number of bins to use for the probability computation.

# Returns
A vector of probabilities.

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


function _trim_extremes(x::Vector{Float64})
    n = length(x)
    trim = round(Int, n * 0.1)
    sorted = sort(x)
    return view(sorted, (trim+1):(n-trim))
end



"""
    get_mean_amplitude(erp_data::ErpData, time_window::Tuple{<:Real, <:Real})

Calculates the mean amplitude for each electrode within a specified time window.

# Arguments
- `erp_data::ErpData`: ERP data structure
- `time_window::Tuple{<:Real, <:Real}`: Time window as (start_time, end_time) in seconds

# Returns
- `DataFrame`: A DataFrame with electrode labels as column names and corresponding mean amplitudes
"""
function get_mean_amplitude(erp_data::ErpData, time_window::Tuple{<:Real,<:Real})
    # Find time indices within the window
    time_indices = findall(x -> x >= time_window[1] && x <= time_window[2], erp_data.time)

    if isempty(time_indices)
        error("No data points found within the specified time window")
    end

    # Calculate mean amplitude for each electrode
    mean_amplitudes = Dict{Symbol,Float64}()
    for electrode in erp_data.layout.label
        if haskey(erp_data.data, electrode)
            mean_amplitudes[electrode] = mean(erp_data.data[time_indices, electrode])
        end
    end

    return DataFrame(mean_amplitudes)
end













"""
    combine_boolean_columns!(dat::ContinuousData, columns::Vector{Symbol}, operation::Symbol; output_column::Symbol = :combined_flags)

Combines multiple boolean columns using a specified logical operation and stores the result in a new column.

# Arguments
- `dat::ContinuousData`: The ContinuousData object containing the boolean columns
- `columns::Vector{Symbol}`: Vector of column names to combine
- `operation::Symbol`: Logical operation to apply (:and, :or, :xor, :nand, :nor, :xnor)
- `output_column::Symbol`: Name of the output column (default: :combined_flags)

# Returns
Nothing. The function modifies the input data in place.

# Examples
```julia
# Combine EOG artifacts (any EOG = vertical OR horizontal)
combine_boolean_columns!(dat, [:is_vEOG, :is_hEOG], :or, output_column = :is_any_EOG)

# Combine quality flags (good data = NOT extreme AND NOT EOG)
combine_boolean_columns!(dat, [:is_extreme_value, :is_any_EOG], :nor, output_column = :is_good_data)

# Combine multiple artifact types (any artifact)
combine_boolean_columns!(dat, [:is_extreme_value, :is_vEOG, :is_hEOG], :or, output_column = :is_any_artifact)

# Create clean data flag (no artifacts)
combine_boolean_columns!(dat, [:is_extreme_value, :is_vEOG, :is_hEOG], :nor, output_column = :is_clean_data)
```

# Available Operations
- `:and` - All columns must be true (logical AND)
- `:or` - At least one column must be true (logical OR)
- `:nand` - Not all columns are true (logical NAND)
- `:nor` - No columns are true (logical NOR)
"""
function combine_boolean_columns!(
    dat::ContinuousData,
    columns::Vector{Symbol},
    operation::Symbol;
    output_column::Symbol = :combined_flags,
)
    # Input validation
    @assert !isempty(columns) "Must specify at least one column to combine"
    @assert all(col -> hasproperty(dat.data, col), columns) "All specified columns must exist in the data"
    @assert operation in [:and, :or, :nand, :nor] "Invalid operation. Must be one of: :and, :or, :nand, :nor"

    # Get the boolean columns
    bool_columns = [dat.data[!, col] for col in columns]

    # Apply the logical operation
    result = if operation == :and
        all.(zip(bool_columns...))
    elseif operation == :or
        any.(zip(bool_columns...))
    elseif operation == :nand
        .!(all.(zip(bool_columns...)))
    elseif operation == :nor
        .!(any.(zip(bool_columns...)))
    end

    # Store the result
    dat.data[!, output_column] = result

    @info "combine_boolean_columns!: Combined $(length(columns)) columns using :$operation operation into column :$output_column"
end
