# using CairoMakie
using BioSemiBDF
using CSV
using DataFrames
using DSP
using GLMakie
using LibGEOS
using LinearAlgebra
using MAT
using OrderedCollections
using Random
using Statistics
using StatsBase: kurtosis


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
        new(
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
function infomax_ica(
    dat::Matrix{Float64},
    data_labels;
    params::IcaPrms = IcaPrms(),
    n_components::Union{Nothing,Int} = nothing,
)

    dat_ica = copy(dat)
    dat_ica = demean(dat_ica)
    dat_ica, scale = scale_data(dat_ica)

    # Apply PCA reduction if n_components is specified
    if !isnothing(n_components)
        # Compute SVD directly on centered data for efficiency
        F = svd(dat_ica)
        # Take only the first n_components
        dat_ica = F.U[:, 1:n_components]' * dat_ica
        # Store PCA components for later reconstruction
        pca_components = F.U[:, 1:n_components] # eigenvectors
    else
        pca_components = nothing
        n_components = size(dat_ica, 1)
    end

    # Sphering 
    sphere, scale_factors = compute_sphere_matrix(dat_ica)
    dat_ica = apply_sphering!(dat_ica, sphere)

    # Pre-allocate all matrices and temporaries
    n_channels = size(dat_ica, 1)
    n_samples = size(dat_ica, 2)
    n_channels_square = n_channels^2
    block = Int(floor(sqrt(n_samples / 3.0)))
    nblock = div(n_samples, block)
    lastt = (nblock - 1) * block + 1

    # Pre-allocate workspace arrays
    weights = Matrix{Float64}(I, n_components, n_components)
    BI = block * Matrix{Float64}(I, n_channels, n_channels)
    u = zeros(n_channels, block)
    y = similar(u)
    data_block = zeros(n_channels, block)
    oldweights = copy(weights)
    startweights = copy(weights)
    permute_indices = Vector{Int}(undef, n_samples)

    # Initialize training variables
    step = 0
    blockno = 1
    wts_blowup = false
    change = 0.0
    oldchange = 0.0
    degconst = 180.0 / π
    angledelta = 0.0
    olddelta = zeros(n_channels_square)
    delta = zeros(n_channels_square)

    # TODO: min/max learning rate?
    while step < params.max_iter

        randperm!(permute_indices)

        for t = 1:block:lastt

            block_indices = permute_indices[t:min(t + block - 1, end)]
            @views data_block .= dat_ica[:, block_indices]
            mul!(u, weights, data_block)
            @. y = 1 / (1 + exp(-u))
            weights = weights + params.l_rate * (BI + (1 .- 2 * y) * u') * weights

            if maximum(abs.(weights)) > params.max_weight
                wts_blowup = true
                change = NaN
                break
            end

            blockno += 1

        end

        if !wts_blowup

            oldwtchange = weights - oldweights
            step += 1

            # Compute and print weight and update angle changes
            angledelta = 0.0
            delta = reshape(oldwtchange, 1, n_channels * n_channels)
            change = dot(delta, delta)
        end

        # Handle weight blowup
        if wts_blowup || isnan(change) || isinf(change)
            @warn "Weight blowup detected - reducing learning rate and restarting"
            step = 0
            change = NaN
            wts_blowup = false
            blockno = 1
            params.l_rate *= params.restart_factor
            weights = copy(startweights)
            oldweights = copy(startweights)
            continue
        end

        # Compute angle change (matching MATLAB)
        if step > 2
            angledelta = acos(dot(delta, olddelta) / sqrt(change * oldchange))
        end
        log_progress(step, change, params.l_rate, degconst * angledelta)

        # Compute weight changes matching MATLAB implementation
        oldweights = copy(weights)

        if step > 2
            if degconst * angledelta > params.anneal_deg
                params.l_rate *= params.anneal_step
                olddelta = copy(delta)
                oldchange = change
            end
        elseif step == 1
            olddelta = copy(delta)
            oldchange = change
        end

        if step > 2 && change < params.w_change
            @info "Weights stabilized - stopping at step $step"
            break
        elseif change > params.blowup      # if weights blow up,
            params.l_rate = params.l_rate * DEFAULT_BLOWUP_FAC    # keep trying
        end

    end

    if !isnothing(pca_components)
        weights = weights * sphere * pca_components[:, 1:n_components]'
    else
        weights = weights * sphere
    end

    # Calculate mixing matrix (topography)
    sphere = Matrix{Float64}(I, size(weights)[2], size(weights)[2])
    winv = pinv(weights * sphere)

    meanvar = sum(winv .^ 2, dims = 1) .* sum((dat_ica') .^ 2, dims = 1) ./ ((n_components * size(dat_ica)[2]) - 1)
    order = sortperm(vec(meanvar), rev = true)  # Get indices in descending order of meanvar

    # Reorder matrices and labels
    weights = weights[order, :]
    unmixing = weights * sphere
    mixing = pinv(unmixing)

    # Generate labels based on actual number of components
    labels = ["IC$i" for i = 1:size(weights, 1)]

    return InfoIca(weights, sphere, mixing, unmixing, scale, labels, data_labels)
end




# function read_mat_file(filename)
#     file = matopen(filename)
#     dat = read(file)
#     close(file)
#     return dat
# end
# dat = read_mat_file("dat.mat")["dat"]
# dat = read_mat_file("dat1.mat")["dat"]
# # @time output = infomax_ica(dat)
# @time output = infomax_ica(dat, n_components = 68)
# @time output = infomax_ica(dat, n_components = 10)








