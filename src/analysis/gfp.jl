"""
Global Field Power (GFP) and Global Dissimilarity calculation.

Global Field Power is a reference-independent measure of response strength computed
as the standard deviation across all channels at each time point. Global Dissimilarity
measures the rate of topographic change across time.

References:
- Lehmann, D. & Skrandies, W. (1980). Reference-free identification of components
  of checkerboard-evoked multichannel potential fields. Electroencephalography and
  Clinical Neurophysiology, 48, 609-621.
- Lehmann, D. & Skrandies, W. (1984). Spatial analysis of evoked potentials in
  man--a review. Progress in Neurobiology, 23, 227-250.
- Skrandies, W. (1990). Global Field Power and Topographic Similarity.
  Brain Topography, 3(1), 137-141.
"""

#=============================================================================
    CORE GFP CALCULATION FUNCTIONS
=============================================================================#

"""
    gfp(dat::ErpData; 
        channel_selection::Function = channels(),
        normalize::Bool = false)::DataFrame

Calculate Global Field Power (GFP) for ERP data.

Global Field Power is computed as the standard deviation across channels at each
time point, providing a reference-independent measure of global response strength.

# Arguments
- `dat::ErpData`: ERP data structure
- `channel_selection::Function`: Channel predicate for selecting channels to include (default: all channels)
- `normalize::Bool`: If true, normalize GFP to 0-100% range (default: false)

# Returns
- `DataFrame`: Contains columns `:time`, `:gfp`, and metadata columns (`:condition`, `:condition_name`, etc.)

# Examples
```julia
using eegfun, JLD2

# Load ERP data
erp_data = load("participant_1_erps.jld2", "erps")[1]

# Calculate GFP using all channels
gfp_result = gfp(erp_data)

# Calculate GFP using specific channels
gfp_result = gfp(erp_data, channel_selection = channels([:C3, :C4, :Cz, :CP3, :CP4, :CPz]))

# Calculate normalized GFP (0-100%)
gfp_result = gfp(erp_data, normalize = true)

# Access the values
time_vector = gfp_result.time
gfp_values = gfp_result.gfp
```

# Notes
- GFP is reference-independent and reflects the overall strength of the electric field
- High GFP values indicate strong, synchronized activity across channels
- Low GFP values indicate weak or desynchronized activity
- Normalization to 0-100% is useful for comparing across different datasets or conditions
"""
function gfp(dat::ErpData; channel_selection::Function = channels(), normalize::Bool = false)::DataFrame

    @info "Calculating Global Field Power (GFP)"

    # Get selected channels (exclude metadata)
    selected_channels = get_selected_channels(dat, channel_selection, include_meta = false, include_extra = false)

    if isempty(selected_channels)
        @minimal_error_throw("No channels selected for GFP calculation")
    end

    @info "Using $(length(selected_channels)) channels for GFP calculation"

    # Extract channel data as matrix: n_timepoints × n_channels
    channel_matrix = Matrix(dat.data[!, selected_channels])

    # Calculate GFP as standard deviation across channels at each time point
    # GFP(t) = std(channels at time t)
    gfp_values = vec(std(channel_matrix, dims = 2))

    # Normalize if requested
    if normalize
        gfp_min = minimum(gfp_values)
        gfp_max = maximum(gfp_values)
        if gfp_max - gfp_min > eps()  # Avoid division by zero
            gfp_values = ((gfp_values .- gfp_min) ./ (gfp_max - gfp_min)) .* 100
            @info "GFP normalized to 0-100% range"
        else
            @minimal_warning "GFP range is zero, normalization skipped"
        end
    end

    # Create output DataFrame with metadata
    result_df = DataFrame()

    # Copy metadata columns
    meta_cols = meta_labels(dat)
    for col in meta_cols
        result_df[!, col] = copy(dat.data[!, col])
    end

    # Add GFP column
    result_df[!, :gfp] = gfp_values

    @info "GFP calculation complete"
    return result_df
end


"""
    gfp(dat::Vector{ErpData}; 
        channel_selection::Function = channels(),
        normalize::Bool = false)::Vector{DataFrame}

Calculate Global Field Power for multiple ERP datasets (e.g., multiple conditions).

# Arguments
- `dat::Vector{ErpData}`: Vector of ERP data structures
- `channel_selection::Function`: Channel predicate for selecting channels (default: all channels)
- `normalize::Bool`: If true, normalize GFP to 0-100% range (default: false)

# Returns
- `Vector{DataFrame}`: Vector of DataFrames, each containing GFP for one condition

# Examples
```julia
# Load ERP data for multiple conditions
erps = load("participant_1_erps.jld2", "erps")

# Calculate GFP for all conditions
gfp_results = gfp(erps)

# Calculate normalized GFP
gfp_results = gfp(erps, normalize = true)

# Plot GFP for each condition
for (i, gfp_data) in enumerate(gfp_results)
    plot(gfp_data.time, gfp_data.gfp, label = "Condition \$i")
end
```
"""
function gfp(dat::Vector{ErpData}; channel_selection::Function = channels(), normalize::Bool = false)::Vector{DataFrame}

    @info "Calculating GFP for $(length(dat)) dataset(s)"

    results = DataFrame[]
    for (i, erp_data) in enumerate(dat)
        @info "Processing dataset $i/$(length(dat))"
        gfp_result = gfp(erp_data; channel_selection = channel_selection, normalize = normalize)
        push!(results, gfp_result)
    end

    return results
end


"""
    global_dissimilarity(dat::ErpData;
                         channel_selection::Function = channels(),
                         normalize::Bool = false)::DataFrame

Calculate Global Dissimilarity (GD) for ERP data.

Global Dissimilarity measures the rate of topographic change by quantifying how much
the normalized potential distribution changes from one time point to the next.

# Arguments
- `dat::ErpData`: ERP data structure
- `channel_selection::Function`: Channel predicate for selecting channels (default: all channels)
- `normalize::Bool`: If true, normalize GD to 0-100% range (default: false)

# Returns
- `DataFrame`: Contains columns `:time`, `:dissimilarity`, and metadata columns

# Examples
```julia
# Calculate global dissimilarity
gd_result = global_dissimilarity(erp_data)

# With specific channels and normalization
gd_result = global_dissimilarity(erp_data, 
                                 channel_selection = channels([:C3, :C4, :Cz]),
                                 normalize = true)
```

# Notes
- Global dissimilarity peaks indicate moments of rapid topographic change
- These peaks may indicate transitions between different brain states or ERP components
- Normalization to 0-100% facilitates comparison across datasets
"""
function global_dissimilarity(
    dat::ErpData;
    channel_selection::Function = channels(),
    normalize::Bool = false,
)::DataFrame

    @info "Calculating Global Dissimilarity (GD)"

    # Get selected channels
    selected_channels = get_selected_channels(dat, channel_selection, include_meta = false, include_extra = false)

    if isempty(selected_channels)
        @minimal_error_throw("No channels selected for dissimilarity calculation")
    end

    @info "Using $(length(selected_channels)) channels for dissimilarity calculation"

    # Extract channel data as matrix: n_timepoints × n_channels
    channel_matrix = Matrix(dat.data[!, selected_channels])
    n_timepoints, n_channels = size(channel_matrix)

    # Calculate GFP for normalization
    gfp_values = vec(std(channel_matrix, dims = 2))

    # Normalize channel data by GFP at each time point
    # This creates a matrix where each row has unit variance
    normalized_map = channel_matrix ./ gfp_values

    # Calculate differences between consecutive time points
    # diff(normalized_map, dims=1) gives (t+1) - t
    map_diff = diff(normalized_map, dims = 1)

    # Global dissimilarity is the mean absolute difference across channels
    gd_values = vec(mean(abs.(map_diff), dims = 2))

    # Replicate first value to maintain same length as time vector
    # (diff reduces length by 1)
    gd_values = vcat(gd_values[1], gd_values)

    # Normalize if requested
    if normalize
        gd_min = minimum(gd_values)
        gd_max = maximum(gd_values)
        if gd_max - gd_min > eps()
            gd_values = ((gd_values .- gd_min) ./ (gd_max - gd_min)) .* 100
            @info "Global Dissimilarity normalized to 0-100% range"
        else
            @minimal_warning "GD range is zero, normalization skipped"
        end
    end

    # Create output DataFrame with metadata
    result_df = DataFrame()

    # Copy metadata columns
    meta_cols = meta_labels(dat)
    for col in meta_cols
        result_df[!, col] = copy(dat.data[!, col])
    end

    # Add dissimilarity column
    result_df[!, :dissimilarity] = gd_values

    @info "Global Dissimilarity calculation complete"
    return result_df
end


"""
    global_dissimilarity(dat::Vector{ErpData};
                         channel_selection::Function = channels(),
                         normalize::Bool = false)::Vector{DataFrame}

Calculate Global Dissimilarity for multiple ERP datasets.

See single-dataset version for details.
"""
function global_dissimilarity(
    dat::Vector{ErpData};
    channel_selection::Function = channels(),
    normalize::Bool = false,
)::Vector{DataFrame}

    @info "Calculating Global Dissimilarity for $(length(dat)) dataset(s)"

    results = DataFrame[]
    for (i, erp_data) in enumerate(dat)
        @info "Processing dataset $i/$(length(dat))"
        gd_result = global_dissimilarity(erp_data; channel_selection = channel_selection, normalize = normalize)
        push!(results, gd_result)
    end

    return results
end


"""
    gfp_and_dissimilarity(dat::ErpData;
                          channel_selection::Function = channels(),
                          normalize::Bool = false)::DataFrame

Calculate both Global Field Power and Global Dissimilarity in one call.

This is a convenience function that computes both metrics efficiently.

# Arguments
- `dat::ErpData`: ERP data structure
- `channel_selection::Function`: Channel predicate for selecting channels (default: all channels)
- `normalize::Bool`: If true, normalize both metrics to 0-100% range (default: false)

# Returns
- `DataFrame`: Contains columns `:time`, `:gfp`, `:dissimilarity`, and metadata columns

# Examples
```julia
# Calculate both metrics
result = gfp_and_dissimilarity(erp_data)

# Access the values
time_vector = result.time
gfp_values = result.gfp
dissimilarity_values = result.dissimilarity

# With normalization
result = gfp_and_dissimilarity(erp_data, normalize = true)
```
"""
function gfp_and_dissimilarity(
    dat::ErpData;
    channel_selection::Function = channels(),
    normalize::Bool = false,
)::DataFrame

    @info "Calculating Global Field Power and Global Dissimilarity"

    # Get selected channels
    selected_channels = get_selected_channels(dat, channel_selection, include_meta = false, include_extra = false)

    if isempty(selected_channels)
        @minimal_error_throw("No channels selected for GFP/dissimilarity calculation")
    end

    @info "Using $(length(selected_channels)) channels"

    # Extract channel data as matrix
    channel_matrix = Matrix(dat.data[!, selected_channels])
    n_timepoints, n_channels = size(channel_matrix)

    # Calculate GFP
    gfp_values = vec(std(channel_matrix, dims = 2))

    # Normalize GFP if requested
    if normalize
        gfp_min = minimum(gfp_values)
        gfp_max = maximum(gfp_values)
        if gfp_max - gfp_min > eps()
            gfp_normalized = ((gfp_values .- gfp_min) ./ (gfp_max - gfp_min)) .* 100
        else
            gfp_normalized = gfp_values
        end
    else
        gfp_normalized = gfp_values
    end

    # Calculate Global Dissimilarity using GFP for normalization
    normalized_map = channel_matrix ./ gfp_values
    map_diff = diff(normalized_map, dims = 1)
    gd_values = vec(mean(abs.(map_diff), dims = 2))
    gd_values = vcat(gd_values[1], gd_values)

    # Normalize dissimilarity if requested
    if normalize
        gd_min = minimum(gd_values)
        gd_max = maximum(gd_values)
        if gd_max - gd_min > eps()
            gd_normalized = ((gd_values .- gd_min) ./ (gd_max - gd_min)) .* 100
        else
            gd_normalized = gd_values
        end
    else
        gd_normalized = gd_values
    end

    # Create output DataFrame
    result_df = DataFrame()

    # Copy metadata columns
    meta_cols = meta_labels(dat)
    for col in meta_cols
        result_df[!, col] = copy(dat.data[!, col])
    end

    # Add computed columns
    result_df[!, :gfp] = gfp_normalized
    result_df[!, :dissimilarity] = gd_normalized

    @info "GFP and Global Dissimilarity calculation complete"
    return result_df
end


"""
    gfp_and_dissimilarity(dat::Vector{ErpData};
                          channel_selection::Function = channels(),
                          normalize::Bool = false)::Vector{DataFrame}

Calculate both GFP and Global Dissimilarity for multiple ERP datasets.

See single-dataset version for details.
"""
function gfp_and_dissimilarity(
    dat::Vector{ErpData};
    channel_selection::Function = channels(),
    normalize::Bool = false,
)::Vector{DataFrame}

    @info "Calculating GFP and Global Dissimilarity for $(length(dat)) dataset(s)"

    results = DataFrame[]
    for (i, erp_data) in enumerate(dat)
        @info "Processing dataset $i/$(length(dat))"
        result = gfp_and_dissimilarity(erp_data; channel_selection = channel_selection, normalize = normalize)
        push!(results, result)
    end

    return results
end
