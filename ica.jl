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
    InfoICA

Structure containing ICA decomposition results.

# Fields
- `topo::Matrix`: Topography matrix (mixing matrix)
- `unmixing::Matrix`: Unmixing matrix (weights)
- `label::Vector{String}`: Component labels
"""
struct InfoIca
    weights::Matrix{Float64}
    sphere::Matrix{Float64}
    mixing::Matrix{Float64}
    unmixing::Matrix{Float64}
    scale::Float64
    ica_label::Vector{String}
    data_label::Vector{String}
end


"""
    log_progress(step::Int, max_iter::Int, change::Float64, l_rate::Float64, angledelta::Float64)
Log progress of ICA computation.
"""
function log_progress(step::Int, change::Float64, l_rate::Float64, angledelta::Float64)
    @info "Step $step, change = $change, lrate = $l_rate, angle = $angledelta"
end

function compute_sphere_matrix(data::Matrix{Float64})
    # Compute the covariance matrix
    cov_matrix = cov(data, dims = 2)
    # Compute the sphering matrix
    sphere = 2.0 * inv(sqrt(cov_matrix))
    # Compute scale factors (diagonal elements of the inverse sphering matrix)
    scale_factors = diag(inv(sphere))
    return sphere, scale_factors
end

function apply_sphering!(data::Matrix{Float64}, sphere::Matrix{Float64})
    data .= sphere * data
    return data
end

function demean(data::Matrix{Float64})
    return data .- mean(data, dims = 2)
end


function scale_data(data::Matrix{Float64})
    scale = sqrt(norm((data * data') / size(data, 2)))
    scaled_data = data / scale
    println("using scaled data")
    return scaled_data, scale
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


function infomax_ica(
    dat::Matrix{Float64},
    data_labels;
    n_components::Union{Nothing,Int} = nothing,
    params::IcaPrms = IcaPrms(),
)
    dat_ica = demean(copy(dat))
    scale = sqrt(norm((dat_ica * dat_ica') / size(dat_ica, 2)))
    dat_ica ./= scale

    # PCA reduction
    if !isnothing(n_components)
        F = svd(dat_ica)
        dat_ica = F.U[:, 1:n_components]' * dat_ica
        pca_components = view(F.U, :, 1:n_components)
    else
        n_components = size(dat_ica, 1)
        pca_components = nothing
    end

    # Sphering 
    sphere = 2.0 * inv(sqrt(cov(dat_ica, dims = 2)))
    dat_ica .= sphere * dat_ica

    # Initialize
    n_channels = size(dat_ica, 1)
    n_samples = size(dat_ica, 2)
    block = min(Int(floor(sqrt(n_samples / 3.0))), 512)
    lastt = block * div(n_samples, block)

    # Pre-allocate arrays
    weights = Matrix{Float64}(I, n_components, n_components)
    oldweights = copy(weights)
    startweights = copy(weights)
    u = zeros(n_channels, block)
    y = zeros(n_channels, block)
    y_temp = zeros(n_channels, block)
    data_block = zeros(n_channels, block)
    delta = zeros(1, n_channels^2)
    olddelta = zeros(1, n_channels^2)
    BI = block * Matrix{Float64}(I, n_channels, n_channels)
    weights_temp = similar(weights)
    bi_weights = similar(weights)
    wu_term = similar(weights)

    # Main loop
    step = 0
    wts_blowup = false
    change = 0.0
    oldchange = 0.0
    angledelta = 0.0

    @inbounds while step < params.max_iter
        permute_indices = randperm(n_samples)

        for t = 1:block:lastt
            block_end = min(t + block - 1, n_samples)
            block_size = block_end - t + 1

            # Extract data block
            data_block[:, 1:block_size] .= view(dat_ica, :, permute_indices[t:block_end])

            # Forward pass
            mul!(u, weights, data_block)
            y .= 1 ./ (1 .+ exp.(-u))
            y_temp .= 1 .- 2 .* y

            # Weight update
            mul!(wu_term, y_temp, u')
            bi_weights .= BI .+ wu_term
            mul!(weights_temp, bi_weights, weights)
            weights .+= params.l_rate .* weights_temp

            if maximum(abs, weights) > params.max_weight
                wts_blowup = true
                change = NaN
                break
            end
        end

        if !wts_blowup
            oldweights .-= weights
            step += 1
            delta .= reshape(oldweights, 1, :)
            change = dot(delta, delta)
        end

        if wts_blowup || isnan(change) || isinf(change)
            step = 0
            change = NaN
            wts_blowup = false
            params.l_rate *= params.restart_factor
            weights .= startweights
            oldweights .= startweights
            continue
        end

        if step > 2
            angledelta = acos(clamp(dot(delta, olddelta) / sqrt(change * oldchange), -1, 1))
            if params.degconst * angledelta > params.anneal_deg
                params.l_rate *= params.anneal_step
                olddelta .= delta
                oldchange = change
            end
        elseif step == 1
            olddelta .= delta
            oldchange = change
        end

        oldweights .= weights

        if step > 2 && change < params.w_change
            break
        elseif change > params.blowup
            params.l_rate *= params.blowup_fac
        end

        log_progress(step, change, params.l_rate, params.degconst * angledelta)
    end

    # Final calculations
    if !isnothing(pca_components)
        weights = weights * sphere * pca_components'
    else
        weights = weights * sphere
    end

    sphere_final = Matrix{Float64}(I, size(weights, 2), size(weights, 2))
    unmixing = weights * sphere_final
    mixing = pinv(unmixing)

    meanvar = vec(sum(abs2, mixing, dims = 1) .* sum(abs2, dat_ica, dims = 2)' ./ (n_components * n_samples - 1))
    order = sortperm(meanvar, rev = true)

    return InfoIca(
        weights[order, :],
        sphere_final,
        mixing[:, order],
        unmixing[order, :],
        scale,
        ["IC$i" for i = 1:size(weights, 1)],
        data_labels,
    )
end

