function run_ica(
    dat::ContinuousData;
    n_components::Union{Nothing,Int} = nothing,
    exclude_channels::Vector{Symbol} = Symbol[],
    exclude_samples::Vector{Symbol} = Symbol[],
    include_samples::Vector{Symbol} = Symbol[],
    hp_filter::Bool = true,
    lp_filter::Bool = false,
    hp_freq::Float64 = 1.0,
    lp_freq::Float64 = 30.0,
    params::IcaPrms = IcaPrms(),
)
    # Create a copy of the data to avoid modifying the original
    dat_ica = deepcopy(dat)
    
    # Apply filters if requested
    if hp_filter
        @info "Applying high-pass filter"
        dat_ica = filter_data(dat_ica, "hp", "iir", hp_freq, order = 1)
    end
    if lp_filter
        @info "Applying low-pass filter"
        dat_ica = filter_data(dat_ica, "lp", "iir", lp_freq, order = 3)
    end

    # Get channels to use
    channels = setdiff(dat_ica.layout.label, exclude_channels)
    if isempty(channels)
        error("No channels available after excluding specified channels")
    end

    # Get samples to use
    samples = trues(size(dat_ica.data, 1))  # Start with all samples
    
    # Apply exclude_samples filters
    for col in exclude_samples
        if hasproperty(dat_ica.data, col)
            samples .&= .!dat_ica.data[!, col]
        end
    end
    
    # Apply include_samples filters
    for col in include_samples
        if hasproperty(dat_ica.data, col)
            samples .&= dat_ica.data[!, col]
        end
    end
    
    # Convert to indices
    sample_indices = findall(samples)
    if isempty(sample_indices)
        error("No samples available after applying sample filters")
    end

    # Set n_components if not specified
    if isnothing(n_components)
        n_components = length(channels) - 1
    elseif n_components > length(channels)
        @warn "Requested $n_components components but only $(length(channels)) channels available. Using $(length(channels) - 1) components instead."
        n_components = length(channels) - 1
    end

    @info "\nRunning ICA with $(length(channels)) channels and $(length(sample_indices)) samples, $(n_components) components"

    # Create data matrix and run ICA
    dat_for_ica = create_ica_data_matrix(dat_ica.data, channels, sample_indices)
    ica_result = infomax_ica(dat_for_ica, channels, n_components = n_components, params = params)

    return ica_result
end


function IcaPrms(;
    l_rate = 0.001,
    max_iter = 512,
    w_change = 1e-6,
    anneal_deg = 60.0,
    anneal_step = 0.9,
    blowup = 1e15,
    blowup_fac = 0.8,
    max_weight = 1e8,
    restart_factor = 0.9,
    degconst = 180.0 / π,
    default_stop = 1e-6,
)
    IcaPrms(
        l_rate,
        max_iter,
        w_change,
        anneal_deg,
        anneal_step,
        blowup,
        blowup_fac,
        max_weight,
        restart_factor,
        degconst,
        default_stop,
    )
end


function create_ica_data_matrix(dat::DataFrame, channels, samples)
    dat = dat[samples, :]
    # need to make sure we have unique samples as with longer epochs there is potential for overlap
    dat = unique(dat, :sample)
    # select only the channels we want
    dat = dat[!, intersect(propertynames(dat), channels)]
    return permutedims(Matrix(dat))
end


mutable struct WorkArrays
    weights::Matrix{Float64}
    BI::Matrix{Float64}
    u::Matrix{Float64}
    y::Matrix{Float64}
    data_block::Matrix{Float64}
    oldweights::Matrix{Float64}
    startweights::Matrix{Float64}
    weights_temp::Matrix{Float64}
    y_temp::Matrix{Float64}
    bi_weights::Matrix{Float64}
    wu_term::Matrix{Float64}
    delta::Matrix{Float64}
    olddelta::Matrix{Float64}
end

function create_work_arrays(n_components::Int, block_size::Int)
    weights = Matrix{Float64}(I, n_components, n_components)  # Initialize as identity matrix
    return WorkArrays(
        weights,  # weights - start as identity matrix
        block_size * Matrix{Float64}(I, n_components, n_components),  # BI
        zeros(n_components, block_size),  # u
        zeros(n_components, block_size),  # y
        zeros(n_components, block_size),  # data_block
        copy(weights),  # oldweights - copy of identity matrix
        copy(weights),  # startweights - copy of identity matrix
        zeros(n_components, n_components),  # weights_temp
        zeros(n_components, block_size),  # y_temp
        zeros(n_components, n_components),  # bi_weights
        zeros(n_components, n_components),  # wu_term
        zeros(1, n_components^2),  # delta
        zeros(1, n_components^2),   # olddelta
    )
end


function infomax_ica(
    dat_ica::Matrix{Float64},
    data_labels;
    n_components::Union{Nothing,Int} = nothing,
    params::IcaPrms = IcaPrms(),
)

    # Store original mean before removing it
    original_mean = vec(mean(dat_ica, dims = 2))
    
    # Center and scale data
    dat_ica .-= original_mean
    scale = sqrt(norm((dat_ica * dat_ica') / size(dat_ica, 2)))
    dat_ica ./= scale

    # PCA reduction
    if !isnothing(n_components)
        F = svd(dat_ica)
        # Ensure we're using the correct dimensions for PCA components
        pca_components = F.U[1:size(dat_ica, 1), 1:n_components]  # Explicitly specify both dimensions
        dat_ica = pca_components' * dat_ica
    else
        n_components = size(dat_ica, 1)
        pca_components = nothing
    end

    # sphering 
    sphere = inv(sqrt(cov(dat_ica, dims = 2)))
    dat_ica = sphere * dat_ica

    # initialize
    n_channels = size(dat_ica, 1)
    n_samples = size(dat_ica, 2)
    block = min(Int(floor(sqrt(n_samples / 3.0))), 512)
    work = create_work_arrays(n_channels, block)
    lastt = block * div(n_samples, block)
    step = 0
    wts_blowup = false
    change = 0.0
    oldchange = 0.0
    angledelta = 0.0

    # pre-allocate permutation vector
    permute_indices = Vector{Int}(undef, n_samples)

    @inbounds while step < params.max_iter
        randperm!(permute_indices)

        for t = 1:block:lastt
            block_end = min(t + block - 1, n_samples)
            block_size = block_end - t + 1

            # extract data block
            copyto!(view(work.data_block, :, 1:block_size), view(dat_ica, :, view(permute_indices, t:block_end)))

            # forward pass
            mul!(work.u, work.weights, work.data_block)
            @. work.y = 1 / (1 + exp(-work.u))
            @. work.y_temp = 1 - 2 * work.y

            # update weights 
            mul!(work.wu_term, work.y_temp, transpose(work.u))
            work.bi_weights .= work.BI .+ work.wu_term
            mul!(work.weights_temp, work.bi_weights, work.weights)
            @. work.weights += params.l_rate * work.weights_temp

            # boom?
            if maximum(abs, work.weights) > params.max_weight
                wts_blowup = true
                change = NaN
                break
            end
        end

        if !wts_blowup
            work.oldweights .-= work.weights
            step += 1
            work.delta .= reshape(work.oldweights, 1, :)
            change = dot(work.delta, work.delta)
        end

        if wts_blowup || isnan(change) || isinf(change)
            step = 0
            change = NaN
            wts_blowup = false
            params.l_rate *= params.restart_factor
            work.weights .= work.startweights
            work.oldweights .= work.startweights
            continue
        end

        if step > 2
            angledelta = acos(clamp(dot(work.delta, work.olddelta) / sqrt(change * oldchange), -1, 1))
            if params.degconst * angledelta > params.anneal_deg
                params.l_rate *= params.anneal_step
                work.olddelta .= work.delta
                oldchange = change
            end
        elseif step == 1
            work.olddelta .= work.delta
            oldchange = change
        end

        work.oldweights .= work.weights

        if step > 2 && change < params.w_change
            break
        elseif change > params.blowup
            params.l_rate *= params.blowup_fac
        end

        @info "Step $step, change = $change, lrate = $(params.l_rate), angle = $((params.degconst) * angledelta)"

    end

    # Final calculations
    if !isnothing(pca_components)
        work.weights = work.weights * sphere * pca_components'
    else
        work.weights = work.weights * sphere
    end

    mixing = pinv(work.weights)

    # calculate total variance explained and order
    meanvar = vec(sum(abs2, mixing, dims = 1) .* sum(abs2, dat_ica, dims = 2)' ./ (n_components * n_samples - 1))
    meanvar_normalized = meanvar ./ sum(meanvar)
    order = sortperm(meanvar_normalized, rev = true)

    return InfoIca(
        work.weights[order, :],
        mixing[:, order],
        sphere,
        meanvar_normalized[order],
        scale,
        original_mean,
        [Symbol("IC$i") for i = 1:size(work.weights, 1)],
        data_labels,
    )

end


function remove_ica_components(dat::DataFrame, ica::InfoIca, components_to_remove::Vector{Int})

    n_components = size(ica.unmixing, 1)
    if !all(1 .<= components_to_remove .<= n_components)
        throw(ArgumentError("Components must be between 1 and $n_components"))
    end
    
    dat_out = deepcopy(dat)
    
    # Get data dimensions
    n_channels = length(ica.data_label)
    
    # Get data and scale it
    data = permutedims(Matrix(dat_out[!, ica.data_label]))
    data .-= ica.mean
    data ./= ica.scale
    
    # Get removed activations before transformation
    removed_activations = view(ica.unmixing, components_to_remove, :) * data
    
    # Pre-compute the transformation matrix
    tra = Matrix(I, n_channels, n_channels) - 
          view(ica.mixing, :, components_to_remove) * 
          view(ica.unmixing, components_to_remove, :)
    
    # Apply transformation and restore scaling
    cleaned_data = tra * data
    cleaned_data .*= ica.scale
    cleaned_data .+= ica.mean
    
    # Create output DataFrame and assign result
    dat_out[!, ica.data_label] .= permutedims(cleaned_data)
    
    return dat_out, removed_activations

end

function restore_original_data(dat::DataFrame, ica::InfoIca, components_removed::Vector{Int}, removed_activations::Matrix{Float64})
    

    n_components = size(ica.unmixing, 1)
    if !all(1 .<= components_removed .<= n_components)
        throw(ArgumentError("Components must be between 1 and $n_components"))
    end
    
    dat_out = deepcopy(dat)
    
    # Get data and scale it
    data = permutedims(Matrix(dat_out[!, ica.data_label]))
    data .-= ica.mean
    data ./= ica.scale
    
    # Get current activations
    activations = ica.unmixing * data
    
    # Restore removed components
    activations[components_removed, :] .= removed_activations
    
    # Back to channel space and restore scaling
    restored_data = ica.mixing * activations
    restored_data .*= ica.scale
    restored_data .+= ica.mean
    
    # Create output and assign result
    dat_out[!, ica.data_label] .= permutedims(restored_data)
    
    return dat_out

end

function ica(dat::ContinuousData; kwargs...) 
    run_ica(dat; kwargs...)
end

"""
    identify_eye_components(ica_result::InfoIca, dat::ContinuousData;
                          vEOG_channel::Symbol=:vEOG,
                          hEOG_channel::Symbol=:hEOG,
                          z_threshold::Float64=3.0,
                          exclude_samples::Union{Nothing,Vector{Symbol}} = [:is_extreme_value])

Identify ICA components potentially related to eye movements based on z-scored correlation.

# Arguments
- `ica_result::InfoIca`: The ICA result object.
- `dat::ContinuousData`: The continuous data containing EOG channels.

# Keyword Arguments
- `vEOG_channel::Symbol`: Name of the vertical EOG channel (default: :vEOG).
- `hEOG_channel::Symbol`: Name of the horizontal EOG channel (default: :hEOG).
- `z_threshold::Float64`: Absolute Z-score threshold for identification (default: 3.0).
- `exclude_samples::Union{Nothing,Vector{Symbol}}`: Optional vector of Bool columns in `dat.data` marking samples to exclude. Defaults to `[:is_extreme_value]`.

# Returns
- `Dict{Symbol, Vector{Int}}`: Dictionary containing:
  - `:vEOG`: Vector of indices identified for vertical eye movements.
  - `:hEOG`: Vector of indices identified for horizontal eye movements.
- `DataFrame`: DataFrame containing detailed correlation metrics per component.
"""
function identify_eye_components(
    ica_result::InfoIca,
    dat::ContinuousData;
    vEOG_channel::Symbol = :vEOG,
    hEOG_channel::Symbol = :hEOG,
    z_threshold::Float64 = 3.0,
    exclude_samples::Union{Nothing, Vector{Symbol}} = nothing
)

    # Check basic inputs
    if !(vEOG_channel in propertynames(dat.data))
        error("Vertical EOG channel $vEOG_channel not found in data")
    end
    if !(hEOG_channel in propertynames(dat.data))
        error("Horizontal EOG channel $hEOG_channel not found in data")
    end

    # Get samples to use
    samples_to_use = _get_samples_to_use(dat, nothing, exclude_samples)
    if isempty(samples_to_use)
        @warn "No samples remaining after applying exclude criteria. Cannot identify eye components."
        return Dict(:vEOG => Int[], :hEOG => Int[]), DataFrame()
    end

    # Get EOG signals for valid samples only
    vEOG = dat.data[samples_to_use, vEOG_channel]
    hEOG = dat.data[samples_to_use, hEOG_channel]

    # Prepare data matrix for valid samples
    relevant_cols = vcat(ica_result.data_label)
    data_subset_df = dat.data[samples_to_use, relevant_cols]
    dat_matrix = permutedims(Matrix(data_subset_df))
    dat_matrix .-= mean(dat_matrix, dims=2)
    dat_matrix ./= ica_result.scale

    # Calculate components for valid samples
    components = ica_result.unmixing * dat_matrix
    n_components = size(components, 1)

    # Function to calculate correlations for all components
    function calculate_correlations(eog_signal)
        corrs = zeros(n_components)
        for comp_idx = 1:n_components 
             corrs[comp_idx] = abs(cor(components[comp_idx, :], eog_signal))
        end
        return corrs
    end

    identified_vEOG = Int[]
    identified_hEOG = Int[]
    vEOG_corr_z = Float64[]
    hEOG_corr_z = Float64[]
    vEOG_corrs = Float64[]
    hEOG_corrs = Float64[]

    # vEOG
    vEOG_corrs = calculate_correlations(vEOG)
    vEOG_corr_z = StatsBase.zscore(vEOG_corrs)
    identified_vEOG = findall(abs.(vEOG_corr_z) .> z_threshold)

    # hEOG
    hEOG_corrs = calculate_correlations(hEOG)
    hEOG_corr_z = StatsBase.zscore(hEOG_corrs)
    identified_hEOG = findall(abs.(hEOG_corr_z) .> z_threshold)

    sort!(identified_vEOG)
    sort!(identified_hEOG)

    result_dict = Dict{Symbol, Vector{Int}}( 
        :vEOG => identified_vEOG,
        :hEOG => identified_hEOG,
    )
    
    metrics_df = DataFrame(
        :Component => 1:n_components,
        :vEOG_corr => vEOG_corrs,
        :vEOG_zscore => vEOG_corr_z,
        :hEOG_corr => hEOG_corrs,
        :hEOG_zscore => hEOG_corr_z
    )

    return result_dict, metrics_df 
end


"""
    identify_ekg_components(ica_result::InfoIca, dat::ContinuousData;
                              min_bpm::Real=40, max_bpm::Real=120,
                              min_prominence_std::Real=2.5,
                              min_peaks::Int=10,
                              max_ibi_std_s::Real=0.05,
                              include_samples::Union{Nothing,Vector{Symbol}} = nothing,
                              exclude_samples::Union{Nothing,Vector{Symbol}} = nothing

Identify ICA components potentially related to EKG artifacts based on peak detection
and interval regularity, using only samples consistent with ICA calculation.

# Arguments
- `ica_result::InfoIca`: The ICA result object.
- `dat::ContinuousData`: The continuous data (needed for sampling rate `fs` and sample selection columns).

# Keyword Arguments
- `min_bpm::Real`: Minimum plausible heart rate in beats per minute (default: 40).
- `max_bpm::Real`: Maximum plausible heart rate in beats per minute (default: 120).
- `min_prominence_std::Real`: Minimum peak prominence in standard deviations above mean (default: 2.5).
- `min_peaks::Int`: Minimum number of prominent peaks required within plausible heart rate range (default: 10).
- `max_ibi_std_s::Real`: Maximum standard deviation of the inter-beat intervals (in seconds) for component to be flagged (default: 0.05).
- `include_samples::Union{Nothing,Vector{Symbol}}`: Optional vector of Bool columns in `dat.data` marking samples to *include*. Defaults to `nothing` (include all unless excluded).
- `exclude_samples::Union{Nothing,Vector{Symbol}}`: Optional vector of Bool columns in `dat.data` marking samples to *exclude*. Defaults to `[:is_extreme_value]`.

# Returns
- `Vector{Int}`: Sorted vector of indices identified as potential EKG components.
- `DataFrame`: DataFrame containing metrics for each component (calculated on the included samples):
  - `:Component`: Component index (1 to n).
  - `:num_peaks`: Number of detected prominent peaks.
  - `:num_valid_ibis`: Number of inter-beat intervals within the plausible BPM range.
  - `:mean_ibi_s`: Mean inter-beat interval in seconds (if num_valid_ibis > 0).
  - `:std_ibi_s`: Standard deviation of inter-beat intervals in seconds (if num_valid_ibis > 1).
  - `:is_ekg_artifact`: Boolean flag indicating if component met the criteria.
"""
function identify_ecg_components(
    ica_result::InfoIca,
    dat::ContinuousData;
    min_bpm::Real=40,
    max_bpm::Real=120,
    min_prominence_std::Real=3,      
    min_peaks::Int=10,                
    max_ibi_std_s::Real=0.2,         
    min_peak_ratio::Real=0.7,         
    include_samples::Union{Nothing,Vector{Symbol}} = nothing,
    exclude_samples::Union{Nothing,Vector{Symbol}} = nothing, 
)

    # Data Preparation 
    samples_to_use = _get_samples_to_use(dat, include_samples, exclude_samples)
    if isempty(samples_to_use)
        @warn "No samples remaining after applying exclude criteria."
        return Int[], DataFrame(Component=Int[], num_peaks=Int[], num_valid_ibis=Int[], mean_ibi_s=Float64[], std_ibi_s=Float64[], peak_ratio=Float64[], is_ecg_artifact=Bool[])
    end

    # Process data
    relevant_cols = ica_result.data_label
    data_subset_df = dat.data[samples_to_use, relevant_cols]
    dat_matrix_subset = permutedims(Matrix(data_subset_df)) 
    dat_matrix_subset .-= mean(dat_matrix_subset, dims=2)
    dat_matrix_subset ./= ica_result.scale
    
    components_subset = ica_result.unmixing * dat_matrix_subset 
    n_components = size(ica_result.unmixing, 1)

    # Convert BPM to plausible IBI range
    min_ibi_s = 60.0 / max_bpm
    max_ibi_s = 60.0 / min_bpm

    # Store results
    metrics = []
    identified_ecg = Int[]

    # Loop through components
    for comp_idx in 1:n_components
        ts = components_subset[comp_idx, :]

        # Find prominent peaks - apply stricter criteria
        # 1. Use a higher threshold
        peak_indices = _findpeaks((ts .* -1); min_prominence_std=min_prominence_std)

        num_peaks = length(peak_indices)
        mean_ibi = NaN
        std_ibi = NaN
        num_valid_ibis = 0
        peak_ratio = 0.0  
        is_ecg = false

        if num_peaks >= 2
            # Calculate IBIs
            ibis_s = diff(peak_indices) ./ dat.sample_rate

            # Filter valid IBIs
            valid_ibi_mask = (ibis_s .>= min_ibi_s) .& (ibis_s .<= max_ibi_s)
            valid_ibis = ibis_s[valid_ibi_mask]
            num_valid_ibis = length(valid_ibis)
            
            # Calculate peak ratio (new metric)
            peak_ratio = num_valid_ibis / (num_peaks - 1)

            if num_valid_ibis > 1 
                mean_ibi = mean(valid_ibis)
                std_ibi = std(valid_ibis)
                
                # Apply stricter criteria
                if num_valid_ibis >= (min_peaks - 1) && 
                   std_ibi <= max_ibi_std_s && 
                   peak_ratio >= min_peak_ratio  # Add peak ratio criterion
                    is_ecg = true
                    push!(identified_ecg, comp_idx)
                end
            elseif num_valid_ibis == 1 
                mean_ibi = valid_ibis[1]
                std_ibi = 0.0
            end
        end

        # Store metrics 
        push!(metrics, (
            Component=comp_idx,
            num_peaks=num_peaks,
            num_valid_ibis=num_valid_ibis,
            mean_ibi_s=mean_ibi,
            std_ibi_s=std_ibi,
            peak_ratio=peak_ratio,
            is_ecg_artifact=is_ecg
        ))
    end

    # Finalize results
    metrics_df = DataFrame(metrics)
    if isempty(metrics)
         metrics_df = DataFrame(Component=1:n_components, num_peaks=0, num_valid_ibis=0, 
                                mean_ibi_s=NaN, std_ibi_s=NaN, peak_ratio=0.0, is_ecg_artifact=false)
    end
    sort!(identified_ecg)

    return identified_ecg, metrics_df
end

"""
    identify_spatial_kurtosis_components(ica_result::InfoIca, dat::ContinuousData;
                                      exclude_samples::Union{Nothing,Vector{Symbol}} = [:is_extreme_value],
                                      z_threshold::Float64 = 3.0)

Identify ICA components with high spatial kurtosis (localized, spot-like activity).

# Arguments
- `ica_result::InfoIca`: The ICA result object.
- `dat::ContinuousData`: The continuous data.

# Keyword Arguments
- `exclude_samples::Union{Nothing,Vector{Symbol}}`: Optional vector of Bool columns in `dat.data` marking samples to exclude. Defaults to `[:is_extreme_value]`.
- `z_threshold::Float64`: Z-score threshold for identifying high spatial kurtosis components (default: 3.0).

# Returns
- `Vector{Int}`: Indices of components with high spatial kurtosis.
- `DataFrame`: DataFrame containing spatial kurtosis values and z-scores for all components.
"""
function identify_spatial_kurtosis_components(
    ica_result::InfoIca,
    dat::ContinuousData;
    exclude_samples::Union{Nothing,Vector{Symbol}} = [:is_extreme_value],
    z_threshold::Float64 = 3.0
)
    # Calculate spatial kurtosis for each component's weights
    n_components = size(ica_result.mixing, 2)
    spatial_kurtosis = Float64[]
    
    for i in 1:n_components
        # Get component weights
        weights = ica_result.mixing[:, i]
        # Calculate kurtosis of the weights
        k = kurtosis(weights)
        push!(spatial_kurtosis, k)
    end

    # Calculate z-scores of spatial kurtosis values
    spatial_kurtosis_z = StatsBase.zscore(spatial_kurtosis)

    # Identify components with high spatial kurtosis (using z-scores)
    high_kurtosis_comps = findall(spatial_kurtosis_z .> z_threshold)  # Only positive deviations (localized activity)
    sort!(high_kurtosis_comps)

    # Create metrics DataFrame
    metrics_df = DataFrame(
        :Component => 1:n_components,
        :SpatialKurtosis => spatial_kurtosis,
        :SpatialKurtosisZScore => spatial_kurtosis_z
    )

    return high_kurtosis_comps, metrics_df
end


"""
    identify_line_noise_components(ica_result::InfoIca, dat::ContinuousData;
                                 exclude_samples::Union{Nothing,Vector{Symbol}} = [:is_extreme_value],
                                 line_freq::Real=50.0,
                                 freq_bandwidth::Real=1.0,
                                 z_threshold::Float64=3.0,
                                 min_harmonic_power::Real=0.5)

Identify ICA components with strong line noise characteristics.

# Arguments
- `ica_result::InfoIca`: The ICA result object.
- `dat::ContinuousData`: The continuous data.

# Keyword Arguments
- `exclude_samples::Union{Nothing,Vector{Symbol}}`: Optional vector of Bool columns in `dat.data` marking samples to exclude. Defaults to `[:is_extreme_value]`.
- `line_freq::Real`: Line frequency in Hz (default: 50.0 for European power).
- `freq_bandwidth::Real`: Bandwidth around line frequency to consider (default: 1.0 Hz).
- `z_threshold::Float64`: Z-score threshold for identifying line noise components (default: 3.0).
- `min_harmonic_power::Real`: Minimum power ratio of harmonics relative to fundamental (default: 0.5).

# Returns
- `Vector{Int}`: Indices of components with strong line noise characteristics.
- `DataFrame`: DataFrame containing spectral metrics for all components.
"""
function identify_line_noise_components(
    ica_result::InfoIca,
    dat::ContinuousData;
    exclude_samples::Union{Nothing,Vector{Symbol}} = [:is_extreme_value],
    line_freq::Real=50.0,
    freq_bandwidth::Real=1.0,
    z_threshold::Float64=3.0,
    min_harmonic_power::Real=0.5
)
    # Get samples to use
    samples_to_use = _get_samples_to_use(dat, nothing, exclude_samples)
    if isempty(samples_to_use)
        @warn "No samples remaining after applying exclude criteria. Cannot identify line noise components."
        return Int[], DataFrame()
    end

    # Prepare data matrix for valid samples
    relevant_cols = vcat(ica_result.data_label)
    data_subset_df = dat.data[samples_to_use, relevant_cols]
    dat_matrix = permutedims(Matrix(data_subset_df))
    dat_matrix .-= mean(dat_matrix, dims=2)
    dat_matrix ./= ica_result.scale

    # Calculate components for valid samples
    components = ica_result.unmixing * dat_matrix
    n_components = size(components, 1)
    fs = dat.sample_rate

    # Calculate power spectrum for each component
    # Use a reasonable FFT size (power of 2, but not too large)
    nfft = min(nextpow(2, size(components, 2)), 2^16)  # Cap at 2^16 points
    freqs = FFTW.rfftfreq(nfft, fs)
    psd = zeros(length(freqs), n_components)
    
    for i in 1:n_components
        # Zero-pad or truncate to nfft points
        signal = components[i, :]
        if length(signal) > nfft
            signal = signal[1:nfft]
        elseif length(signal) < nfft
            signal = [signal; zeros(nfft - length(signal))]
        end
        psd[:, i] = abs2.(FFTW.rfft(signal))
    end

    # Find indices for line frequency and harmonics
    line_idx = findmin(abs.(freqs .- line_freq))[2]
    line_band = findall(abs.(freqs .- line_freq) .<= freq_bandwidth)
    
    # Calculate metrics for each component
    metrics = []
    for i in 1:n_components
        # Get power at line frequency and surrounding band
        line_power = mean(psd[line_band, i])
        
        # Calculate power in surrounding bands (excluding line frequency)
        surrounding_bands = setdiff(1:length(freqs), line_band)
        surrounding_power = mean(psd[surrounding_bands, i])
        
        # Calculate power ratio
        power_ratio = line_power / (surrounding_power + eps())
        
        # Check for harmonics (2x and 3x line frequency)
        harmonic_powers = Float64[]
        for h in 2:3
            harmonic_freq = line_freq * h
            harmonic_idx = findmin(abs.(freqs .- harmonic_freq))[2]
            harmonic_band = findall(abs.(freqs .- harmonic_freq) .<= freq_bandwidth)
            harmonic_power = mean(psd[harmonic_band, i])
            push!(harmonic_powers, harmonic_power / line_power)
        end
        
        # Store metrics
        push!(metrics, (
            Component=i,
            LinePower=line_power,
            SurroundingPower=surrounding_power,
            PowerRatio=power_ratio,
            Harmonic2Ratio=harmonic_powers[1],
            Harmonic3Ratio=harmonic_powers[2]
        ))
    end

    # Create metrics DataFrame
    metrics_df = DataFrame(metrics)
    
    # Calculate z-scores of power ratios
    power_ratio_z = StatsBase.zscore(metrics_df.PowerRatio)
    metrics_df[!, :PowerRatioZScore] = power_ratio_z

    # Identify components with strong line noise characteristics
    line_noise_comps = findall(power_ratio_z .> z_threshold)
    
    # Additional check for harmonics
    if !isempty(line_noise_comps)
        harmonic_mask = (metrics_df.Harmonic2Ratio .> min_harmonic_power) .| 
                       (metrics_df.Harmonic3Ratio .> min_harmonic_power)
        line_noise_comps = intersect(line_noise_comps, findall(harmonic_mask))
    end
    
    sort!(line_noise_comps)

    return line_noise_comps, metrics_df
end


"""
    identify_line_noise_components(ica_result::InfoIca, dat::ContinuousData;
                                 exclude_samples::Union{Nothing,Vector{Symbol}} = [:is_extreme_value],
                                 line_freq::Real=50.0,
                                 freq_bandwidth::Real=1.0,
                                 z_threshold::Float64=3.0,
                                 min_harmonic_power::Real=0.5)

Identify ICA components with strong line noise characteristics.

# Arguments
- `ica_result::InfoIca`: The ICA result object.
- `dat::ContinuousData`: The continuous data.

# Keyword Arguments
- `exclude_samples::Union{Nothing,Vector{Symbol}}`: Optional vector of Bool columns in `dat.data` marking samples to exclude. Defaults to `[:is_extreme_value]`.
- `line_freq::Real`: Line frequency in Hz (default: 50.0 for European power).
- `freq_bandwidth::Real`: Bandwidth around line frequency to consider (default: 1.0 Hz).
- `z_threshold::Float64`: Z-score threshold for identifying line noise components (default: 3.0).
- `min_harmonic_power::Real`: Minimum power ratio of harmonics relative to fundamental (default: 0.5).

# Returns
- `Vector{Int}`: Indices of components with strong line noise characteristics.
- `DataFrame`: DataFrame containing spectral metrics for all components.
"""
function identify_line_noise_components(
    ica_result::InfoIca,
    dat::ContinuousData;
    exclude_samples::Union{Nothing,Vector{Symbol}} = [:is_extreme_value],
    line_freq::Real=50.0,
    freq_bandwidth::Real=1.0,
    z_threshold::Float64=3.0,
    min_harmonic_power::Real=0.5
)
    # Get samples to use
    samples_to_use = _get_samples_to_use(dat, nothing, exclude_samples)
    if isempty(samples_to_use)
        @warn "No samples remaining after applying exclude criteria. Cannot identify line noise components."
        return Int[], DataFrame()
    end

    # Prepare data matrix for valid samples
    relevant_cols = vcat(ica_result.data_label)
    data_subset_df = dat.data[samples_to_use, relevant_cols]
    dat_matrix = permutedims(Matrix(data_subset_df))
    dat_matrix .-= mean(dat_matrix, dims=2)
    dat_matrix ./= ica_result.scale

    # Calculate components for valid samples
    components = ica_result.unmixing * dat_matrix
    n_components = size(components, 1)
    fs = dat.sample_rate

    # Calculate power spectrum for each component
    # Use a reasonable FFT size (power of 2, but not too large)
    nfft = min(nextpow(2, size(components, 2)), 2^16)  # Cap at 2^16 points
    freqs = FFTW.rfftfreq(nfft, fs)
    psd = zeros(length(freqs), n_components)
    
    for i in 1:n_components
        # Zero-pad or truncate to nfft points
        signal = components[i, :]
        if length(signal) > nfft
            signal = signal[1:nfft]
        elseif length(signal) < nfft
            signal = [signal; zeros(nfft - length(signal))]
        end
        psd[:, i] = abs2.(FFTW.rfft(signal))
    end

    # Find indices for line frequency and harmonics
    line_idx = findmin(abs.(freqs .- line_freq))[2]
    line_band = findall(abs.(freqs .- line_freq) .<= freq_bandwidth)
    
    # Calculate metrics for each component
    metrics = []
    for i in 1:n_components
        # Get power at line frequency and surrounding band
        line_power = mean(psd[line_band, i])
        
        # Calculate power in surrounding bands (excluding line frequency)
        surrounding_bands = setdiff(1:length(freqs), line_band)
        surrounding_power = mean(psd[surrounding_bands, i])
        
        # Calculate power ratio
        power_ratio = line_power / (surrounding_power + eps())
        
        # Check for harmonics (2x and 3x line frequency)
        harmonic_powers = Float64[]
        for h in 2:3
            harmonic_freq = line_freq * h
            harmonic_idx = findmin(abs.(freqs .- harmonic_freq))[2]
            harmonic_band = findall(abs.(freqs .- harmonic_freq) .<= freq_bandwidth)
            harmonic_power = mean(psd[harmonic_band, i])
            push!(harmonic_powers, harmonic_power / line_power)
        end
        
        # Store metrics
        push!(metrics, (
            Component=i,
            LinePower=line_power,
            SurroundingPower=surrounding_power,
            PowerRatio=power_ratio,
            Harmonic2Ratio=harmonic_powers[1],
            Harmonic3Ratio=harmonic_powers[2]
        ))
    end

    # Create metrics DataFrame
    metrics_df = DataFrame(metrics)
    
    # Calculate z-scores of power ratios
    power_ratio_z = StatsBase.zscore(metrics_df.PowerRatio)
    metrics_df[!, :PowerRatioZScore] = power_ratio_z

    # Identify components with strong line noise characteristics
    line_noise_comps = findall(power_ratio_z .> z_threshold)
    
    # Additional check for harmonics
    if !isempty(line_noise_comps)
        harmonic_mask = (metrics_df.Harmonic2Ratio .> min_harmonic_power) .| 
                       (metrics_df.Harmonic3Ratio .> min_harmonic_power)
        line_noise_comps = intersect(line_noise_comps, findall(harmonic_mask))
    end
    
    sort!(line_noise_comps)

    return line_noise_comps, metrics_df
end