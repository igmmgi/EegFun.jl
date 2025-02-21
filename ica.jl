using BioSemiBDF
using CSV
using DataFrames
using DSP
using GLMakie
using LibGEOS
using LinearAlgebra
using OrderedCollections
using Random
using Statistics
using StatsBase: kurtosis
using FFTW
using Printf


"""
    IcaPrms

Structure containing ICA algorithm parameters.

# Fields
- `l_rate::Float64`: Initial learning rate (default: 0.001)
- `max_iter::Int`: Maximum number of iterations (default: 512)
- `w_change::Float64`: Change threshold for stopping (default: 1e-6)
- `anneal_deg::Float64`: Angle for learning rate reduction (default: 60.0)
- `anneal_step::Float64`: Learning rate reduction factor (default: 0.9)
- `blowup::Float64`: Maximum weight change allowed (default: 1e15)
- `blowup_fac::Float64`: Learning rate reduction on blowup (default: 0.8)
- `max_weight::Float64`: Maximum weight magnitude (default: 1e8)
- `restart_factor::Float64`: Learning rate factor on restart (default: 0.9)
- `degconst::Float64`: Degrees to radians conversion (default: 180.0 / π)
- `default_stop::Float64`: Default weight stop criterion (default: 1e-6)
"""
mutable struct IcaPrms
    l_rate::Float64
    max_iter::Int
    w_change::Float64
    anneal_deg::Float64
    anneal_step::Float64
    blowup::Float64
    blowup_fac::Float64
    max_weight::Float64
    restart_factor::Float64
    degconst::Float64
    default_stop::Float64
end

# Outer constructor with default values
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

"""
    InfoIca

Structure containing ICA decomposition results.

# Fields
- `unmixing::Matrix{Float64}`: Transforms data to ICs
- `mixing::Matrix{Float64}`: Transforms ICs back to data
- `activation::Matrix{Float64}`: Activation matrix
- `variance::Vector{Float64}`: Variance explained by each component
- `scale::Float64`: Data scaling factor
- `mean::Vector{Float64}`: Mean of the data
- `ica_label::Vector{String}`: Component labels
- `data_label::Vector{String}`: Original data channel labels
"""
struct InfoIca
    unmixing::Matrix{Float64}
    mixing::Matrix{Float64}
    activation::Matrix{Float64}
    variance::Vector{Float64}
    scale::Float64
    mean::Vector{Float64}
    ica_label::Vector{String}
    data_label::Vector{String}
end



"""
    infomax_ica(data::Matrix{Float64}; params::IcaPrms = IcaPrms(), extended::Bool = false, ext_params::Union{Nothing, ExtendedIcaPrms} = nothing, n_components::Union{Nothing,Int} = nothing)

Perform Infomax ICA decomposition on EEG data using the algorithm from EEGLAB's runica.

# Arguments
- `data::Matrix{Float64}`: Data matrix (channels × samples)
- `params::IcaPrms`: ICA parameters
- `extended::Bool`: Whether to use extended Infomax
- `ext_params::ExtendedIcaPrms`: Extended ICA parameters
- `n_components::Union{Nothing,Int}`: Number of components for PCA reduction (default: no reduction)

# Returns
- `InfoICA`: Structure containing:
  - `topo`: Topography matrix (mixing matrix)
  - `unmixing`: Unmixing matrix (weights)
  - `label`: Component labels
"""

function to_data_frame(dat::EpochData)
    return vcat(dat.data...)
end

function to_data_frame(dat::Vector{EpochData})
    return vcat([vcat(dat[idx].data[:]...) for idx in eachindex(dat)]...)
end

function create_ica_data_matrix(dat::DataFrame, channels; samples_to_include = nothing)
    if isnothing(samples_to_include)
        samples_to_include = dat.sample
    end
    dat = dat[samples_to_include, :]
    # need to make sure we have unique samples as with longer epochs there is potential for overlap
    dat = unique(dat, :sample)
    # select only the channels we want
    dat = dat[!, intersect(names(dat), channels)]
    return permutedims(Matrix(dat))
end




# function infomax_ica(
#     dat::ContinuousData,
#     data_labels;
#     n_components::Union{Nothing,Int} = nothing,
#     params::IcaPrms = IcaPrms(),
# )
#     # select actual eeg data columns
#     dat = permutedims(Matrix(dat.data[!, intersect(names(dat.data), data_labels)]))
#     return infomax_ica(dat, data_labels, params = params, n_components = n_components)
# end


# Pre-allocate all arrays in a single struct for better cache locality
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
    sample_rate::Int = 256,
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

    # Calculate total variance explained
    meanvar = vec(sum(abs2, mixing, dims = 1) .* sum(abs2, dat_ica, dims = 2)' ./ (n_components * n_samples - 1))
    meanvar_normalized = meanvar ./ sum(meanvar)
    order = sortperm(meanvar_normalized, rev = true)

    return InfoIca(
        work.weights[order, :],
        mixing[:, order],
        work.weights,
        meanvar_normalized[order],
        scale,
        original_mean,
        ["IC$i" for i = 1:size(work.weights, 1)],
        data_labels,
    )
end



function remove_ica_components(dat::ContinuousData, ica_result::InfoIca, components_to_remove::Vector{Int})
    dat_out = deepcopy(dat)

    ica_channels = ica_result.data_label

    # Create transformation matrix for ICA channels
    tra = Matrix(I, length(ica_channels), length(ica_channels)) - ica_result.mixing[:, components_to_remove] * ica_result.unmixing[components_to_remove, :]

    # Apply transformation to ICA channels
    cleaned_ica_data = tra * Matrix(dat_out.data[!, ica_channels])'

    dat_out.data[!, ica_channels] = cleaned_ica_data'

    return dat_out

end

function remove_ica_components(dat::DataFrame, ica_result::InfoIca, components_to_remove::Vector{Int})
    dat_out = deepcopy(dat)

    ica_channels = ica_result.data_label

    # Create transformation matrix for ICA channels
    tra = Matrix(I, length(ica_channels), length(ica_channels)) - ica_result.mixing[:, components_to_remove] * ica_result.unmixing[components_to_remove, :]

    # Apply transformation to ICA channels
    cleaned_ica_data = tra * Matrix(dat[!, ica_channels])'

    dat_out[!, ica_channels] = cleaned_ica_data'

    return dat_out

end

function restore_original_data(dat::ContinuousData, ica_result::InfoIca, components_removed::Vector{Int})
    println("Components being restored: ", components_removed)
    println("Scale factor: ", ica_result.scale)
    println("Mean values: ", ica_result.mean[1:5])  # First 5 means
    
    dat_out = deepcopy(dat)
    
    # Reconstruct the full ICA decomposition
    ica_components = ica_result.mixing[:, components_removed] * ica_result.unmixing[components_removed, :] * Matrix(dat_out.data[!, ica_result.data_label])'
    
    # Create new data with components added
    new_data = Matrix(dat_out.data[!, ica_result.data_label]) .+ ica_components' * ica_result.scale
    
    # Add mean
    new_data .+= ica_result.mean'
    
    # Replace the columns in the DataFrame
    dat_out.data[!, ica_result.data_label] .= new_data
    
    return dat_out
end


function restore_original_data(dat::DataFrame, ica_result::InfoIca, components_removed::Vector{Int})
    dat_out = deepcopy(dat)
    # Reconstruct the full ICA decomposition
    ica_components = ica_result.mixing[:, components_removed] * ica_result.unmixing[components_removed, :] * Matrix(dat_out[!, ica_result.data_label])'

    # Restore the original data by adding the ICA components and scaling
    dat_out[!, ica_result.data_label] .+= ica_components' * ica_result.scale .+ ica_result.mean'

    return dat_out
end












