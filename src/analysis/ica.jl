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


function remove_ica_components(dat::ContinuousData, ica::InfoIca, components_to_remove::Vector{Int})
    remove_ica_components(dat.data, ica, components_to_remove)
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

function restore_original_data(dat::ContinuousData, ica::InfoIca, components_removed::Vector{Int}, removed_activations::Matrix{Float64})
    restore_original_data(dat.data, ica, components_removed, removed_activations)
end

