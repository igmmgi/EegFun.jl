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

# ==============================================================================
#   INTERNAL HELPERS
# ==============================================================================

"""Normalize values to 0-100% range. Returns original values if range is zero."""
function _normalize_to_percent(values::Vector{Float64})
    v_min = minimum(values)
    v_max = maximum(values)
    if v_max - v_min > eps()
        return ((values .- v_min) ./ (v_max - v_min)) .* 100
    else
        return values
    end
end

"""Build a result DataFrame with metadata columns from the source data, plus named result columns."""
function _build_gfp_result_df(dat::ErpData, result_pairs::Pair{Symbol,Vector{Float64}}...)
    result_df = DataFrame()

    # Copy metadata columns
    for col in meta_labels(dat)
        result_df[!, col] = copy(dat.data[!, col])
    end

    # Add result columns
    for (name, values) in result_pairs
        result_df[!, name] = values
    end

    return result_df
end

"""Get channel matrix and selected channels, with validation."""
function _get_channel_matrix(dat::ErpData, channel_selection::Function, operation::String)
    selected_channels = get_selected_channels(dat, channel_selection, include_meta = false, include_extra = false)

    if isempty(selected_channels)
        @minimal_error("No channels selected for $operation")
    end

    @info "Using $(length(selected_channels)) channels for $operation"
    channel_matrix = Matrix(dat.data[!, selected_channels])

    return channel_matrix
end

"""Apply a single-ErpData function across a vector of ErpData, with condition filtering."""
function _apply_to_conditions(
    fn::Function,
    dat::Vector{ErpData};
    condition_selection::Function = conditions(),
    operation::String = "calculation",
    kwargs...,
)::Vector{DataFrame}
    dat_filtered = dat[get_selected_conditions(dat, condition_selection)]

    @info "$operation for $(length(dat_filtered)) dataset(s)"

    return [fn(erp_data; kwargs...) for erp_data in dat_filtered]
end

# ==============================================================================
#   GFP
# ==============================================================================

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
using EegFun

# Load ERP data
erp_data = EegFun.read_data("participant_1_erps.jld2")

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

    channel_matrix = _get_channel_matrix(dat, channel_selection, "GFP calculation")

    # GFP(t) = std(channels at time t) - population std (corrected=false)
    gfp_values = vec(std(channel_matrix, dims = 2, corrected = false))

    if normalize
        gfp_values = _normalize_to_percent(gfp_values)
        @info "GFP normalized to 0-100% range"
    end

    @info "GFP calculation complete"
    return _build_gfp_result_df(dat, :gfp => gfp_values)
end


"""
    gfp(dat::Vector{ErpData}; 
        condition_selection::Function = conditions(),
        channel_selection::Function = channels(),
        normalize::Bool = false)::Vector{DataFrame}

Calculate Global Field Power for multiple ERP datasets (e.g., multiple conditions).

# Arguments
- `dat::Vector{ErpData}`: Vector of ERP data structures
- `condition_selection::Function`: Condition predicate for selecting conditions (default: all conditions)
- `channel_selection::Function`: Channel predicate for selecting channels (default: all channels)
- `normalize::Bool`: If true, normalize GFP to 0-100% range (default: false)

# Returns
- `Vector{DataFrame}`: Vector of DataFrames, each containing GFP for one condition

# Examples
```julia
# Load ERP data for multiple conditions
erps = EegFun.read_data("participant_1_erps.jld2")

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
function gfp(
    dat::Vector{ErpData};
    condition_selection::Function = conditions(),
    channel_selection::Function = channels(),
    normalize::Bool = false,
)::Vector{DataFrame}
    return _apply_to_conditions(gfp, dat; condition_selection, operation = "Calculating GFP", channel_selection, normalize)
end


# ==============================================================================
#   GLOBAL DISSIMILARITY
# ==============================================================================

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
function global_dissimilarity(dat::ErpData; channel_selection::Function = channels(), normalize::Bool = false)::DataFrame

    @info "Calculating Global Dissimilarity (GD)"

    channel_matrix = _get_channel_matrix(dat, channel_selection, "dissimilarity calculation")

    # Normalize channel data by GFP at each time point
    gfp_values = vec(std(channel_matrix, dims = 2, corrected = false))
    normalized_map = channel_matrix ./ gfp_values

    # GD(t) = sqrt(mean((u_i(t) - u_i(t-1))^2))
    map_diff = diff(normalized_map, dims = 1)
    gd_values = vec(sqrt.(mean(map_diff .^ 2, dims = 2)))

    # Replicate first value to maintain same length as time vector
    gd_values = vcat(gd_values[1], gd_values)

    if normalize
        gd_values = _normalize_to_percent(gd_values)
        @info "Global Dissimilarity normalized to 0-100% range"
    end

    @info "Global Dissimilarity calculation complete"
    return _build_gfp_result_df(dat, :dissimilarity => gd_values)
end


"""
    global_dissimilarity(dat::Vector{ErpData};
                         condition_selection::Function = conditions(),
                         channel_selection::Function = channels(),
                         normalize::Bool = false)::Vector{DataFrame}

Calculate Global Dissimilarity for multiple ERP datasets.

See single-dataset version for details.
"""
function global_dissimilarity(
    dat::Vector{ErpData};
    condition_selection::Function = conditions(),
    channel_selection::Function = channels(),
    normalize::Bool = false,
)::Vector{DataFrame}
    return _apply_to_conditions(
        global_dissimilarity,
        dat;
        condition_selection,
        operation = "Calculating Global Dissimilarity",
        channel_selection,
        normalize,
    )
end


# ==============================================================================
#   GFP AND DISSIMILARITY (combined)
# ==============================================================================

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
function gfp_and_dissimilarity(dat::ErpData; channel_selection::Function = channels(), normalize::Bool = false)::DataFrame

    @info "Calculating Global Field Power and Global Dissimilarity"

    channel_matrix = _get_channel_matrix(dat, channel_selection, "GFP/dissimilarity calculation")

    # Calculate GFP (population std)
    gfp_values = vec(std(channel_matrix, dims = 2, corrected = false))

    # Calculate Global Dissimilarity using GFP for normalization
    normalized_map = channel_matrix ./ gfp_values
    map_diff = diff(normalized_map, dims = 1)
    gd_values = vec(sqrt.(mean(map_diff .^ 2, dims = 2)))
    gd_values = vcat(gd_values[1], gd_values)

    if normalize
        gfp_values = _normalize_to_percent(gfp_values)
        gd_values = _normalize_to_percent(gd_values)
    end

    @info "GFP and Global Dissimilarity calculation complete"
    return _build_gfp_result_df(dat, :gfp => gfp_values, :dissimilarity => gd_values)
end


"""
    gfp_and_dissimilarity(dat::Vector{ErpData};
                          condition_selection::Function = conditions(),
                          channel_selection::Function = channels(),
                          normalize::Bool = false)::Vector{DataFrame}

Calculate both GFP and Global Dissimilarity for multiple ERP datasets.

# Arguments
- `dat::Vector{ErpData}`: Vector of ERP data structures
- `condition_selection::Function`: Condition predicate for selecting conditions (default: all conditions)
- `channel_selection::Function`: Channel predicate for selecting channels (default: all channels)
- `normalize::Bool`: If true, normalize both metrics to 0-100% range (default: false)

# Returns
- `Vector{DataFrame}`: Vector of DataFrames, each containing GFP and dissimilarity for one condition

See single-dataset version for additional details.
"""
function gfp_and_dissimilarity(
    dat::Vector{ErpData};
    condition_selection::Function = conditions(),
    channel_selection::Function = channels(),
    normalize::Bool = false,
)::Vector{DataFrame}
    return _apply_to_conditions(
        gfp_and_dissimilarity,
        dat;
        condition_selection,
        operation = "Calculating GFP and Global Dissimilarity",
        channel_selection,
        normalize,
    )
end
