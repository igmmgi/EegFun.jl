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
- `use_bias::Bool`: Whether to use bias term (default: true)
- `anneal_deg::Float64`: Angle for learning rate reduction (default: 60.0)
- `anneal_step::Float64`: Learning rate reduction factor (default: 0.9)
- `blowup::Float64`: Maximum weight change allowed (default: 1e15)
- `blowup_fac::Float64`: Learning rate reduction on blowup (default: 0.8)
- `n_small_angle::Int`: Number of angles below threshold (default: 20)
- `max_weight::Float64`: Maximum weight magnitude (default: 1e8)
- `restart_factor::Float64`: Learning rate factor on restart (default: 0.9)
- `degconst::Float64`: Degrees to radians conversion (default: 180.0 / π)
- `default_stop::Float64`: Default weight stop criterion (default: 1e-6)
- `momentum::Float64`: Momentum parameter (default: 0.0)
"""
mutable struct IcaPrms
    l_rate::Float64
    max_iter::Int
    w_change::Float64
    use_bias::Bool
    anneal_deg::Float64
    anneal_step::Float64
    blowup::Float64
    blowup_fac::Float64
    n_small_angle::Int
    max_weight::Float64
    restart_factor::Float64
    degconst::Float64
    default_stop::Float64
    momentum::Float64

    function IcaPrms(;
        l_rate = 0.001,
        max_iter = 512,
        w_change = 1e-6,
        use_bias = true,
        anneal_deg = 60.0,
        anneal_step = 0.9,
        blowup = 1e15,
        blowup_fac = 0.8,
        n_small_angle = 20,
        max_weight = 1e8,
        restart_factor = 0.9,
        degconst = 180.0 / π,
        default_stop = 1e-6,
        momentum = 0.0,
    )
        new(
            l_rate,
            max_iter,
            w_change,
            use_bias,
            anneal_deg,
            anneal_step,
            blowup,
            blowup_fac,
            n_small_angle,
            max_weight,
            restart_factor,
            degconst,
            default_stop,
            momentum,
        )
    end
end

"""
    ExtendedIcaPrms

Structure containing extended ICA algorithm parameters.

# Fields
- `ext_blocks::Int`: Number of blocks between kurtosis estimation (default: 1)
- `n_subgauss::Int`: Number of sub-gaussian components (default: 1)
- `extmomentum::Float64`: Momentum for kurtosis averaging (default: 0.5)
- `signsbias::Float64`: Bias for sign changes (default: 0.02)
- `signcount_threshold::Int`: Threshold for sign change detection (default: 25)
- `signcount_step::Int`: Step size for sign change adaptation (default: 2)
- `kurt_size::Int`: Number of samples for kurtosis estimation (default: min(6000, n_samples))
"""
mutable struct ExtendedIcaPrms
    ext_blocks::Int
    n_subgauss::Int
    extmomentum::Float64
    signsbias::Float64
    signcount_threshold::Int
    signcount_step::Int
    kurt_size::Int

    # default parameters matching MATLAB
    function ExtendedIcaPrms(;
        ext_blocks = 1,
        n_subgauss = 1,
        extmomentum = 0.5,
        signsbias = 0.02,
        signcount_threshold = 25,
        signcount_step = 2,
        kurt_size = 6000,
    )
        new(ext_blocks, n_subgauss, extmomentum, signsbias, signcount_threshold, signcount_step, kurt_size)
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
    unmixing::Matrix{Float64}
    topo::Matrix{Float64}
    label::Vector{String}
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
    # Spher the data
    data .= sphere * data
    return data
end

"""
    infomax_ica(data::Matrix{Float64}; params::IcaPrms = IcaPrms(), extended::Bool = false, ext_params::Union{Nothing, ExtendedIcaPrms} = nothing, n_components::Union{Nothing,Int} = nothing, do_sphering::Bool = true)

Perform Infomax ICA decomposition on EEG data using the algorithm from EEGLAB's runica.

# Arguments
- `data::Matrix{Float64}`: Data matrix (channels × samples)
- `params::IcaPrms`: ICA parameters
- `extended::Bool`: Whether to use extended Infomax
- `ext_params::ExtendedIcaPrms`: Extended ICA parameters
- `n_components::Union{Nothing,Int}`: Number of components for PCA reduction (default: no reduction)
- `do_sphering::Bool`: Whether to perform sphering (default: true)

# Returns
- `InfoICA`: Structure containing:
  - `topo`: Topography matrix (mixing matrix)
  - `unmixing`: Unmixing matrix (weights)
  - `label`: Component labels
"""
function infomax_ica(
    dat::Matrix{Float64};
    extended::Bool = true,
    params::IcaPrms = IcaPrms(),
    ext_params::ExtendedIcaPrms = ExtendedIcaPrms(),
    n_components::Union{Nothing,Int} = nothing,
    do_sphering::Bool = true,
)

    dat_ica = copy(dat)
    # Subtract row means from data
    dat_ica .-= mean(dat_ica, dims = 2)

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
    end

    #weights = Matrix{Float64}(I, n_components, n_channels)
    # Apply sphering if requested
    if do_sphering
        sphere, scale_factors = compute_sphere_matrix(dat_ica)
        dat_ica = apply_sphering!(dat_ica, sphere)
    else
        sphere = Matrix{Float64}(I, size(dat_ica, 1), size(dat_ica, 1))
        scale_factors = ones(size(dat, 1))
    end

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
    bias = zeros(n_channels, 1)
    onesrow = ones(1, block)
    u = zeros(n_channels, block)
    y = similar(u)
    delta_weights = similar(weights)
    temp_matrix = similar(weights)
    yu = similar(weights)
    uu = similar(weights)
    data_block = zeros(n_channels, block)
    oldweights = copy(weights)
    startweights = copy(weights)
    permute_indices = Vector{Int}(undef, n_samples)

    # For extended ICA - pre-allocate more arrays to avoid allocations
    if extended
        signs = ones(n_channels)
        signs[1:ext_params.n_subgauss] .= -1
        signs = Diagonal(signs)
        old_kurt = zeros(n_channels)
        oldsigns = Diagonal(copy(diag(signs)))  # Store as diagonal matrix
        signcount = 0
        partact = zeros(n_channels, min(ext_params.kurt_size, n_samples))
        kk = zeros(n_channels)
        m2 = zeros(n_channels)
        m4 = zeros(n_channels)
        new_signs = zeros(n_channels)
    end

    # Initialize training variables
    step = 0
    blockno = 1
    wts_blowup = false
    change = 0.0
    oldchange = 0.0
    degconst = 180.0 / π
    angledelta = 0.0
    lrates = zeros(params.max_iter)
    olddelta = zeros(n_channels_square)
    delta = zeros(n_channels_square)

    # Constants for learning rate management
    # MIN_LRATE = 1e-6  # Minimum learning rate
    # MAX_LRATE = 0.1   # Maximum learning rate
    DEFAULT_RESTART_FAC = 0.9  # Factor for reducing learning rate on restart

    while step < params.max_iter
        # step += 1
        # wts_blowup = false

        randperm!(permute_indices)

        for t = 1:block:lastt

            block_indices = permute_indices[t:min(t + block - 1, end)]
            @views data_block .= dat_ica[:, block_indices]
            mul!(u, weights, data_block)
            #u .+= bias .* onesrow

            if extended
                # Extended ICA update with stability checks
                @views @. y = tanh(u)
                mul!(yu, y, u')
                mul!(uu, u, u')

                # Modified weight update for extended ICA
                @views @. temp_matrix = BI - signs * yu - uu
                mul!(delta_weights, temp_matrix, weights)

                if maximum(abs.(delta_weights)) > params.blowup
                    wts_blowup = true
                    break
                end

                @views @. weights += params.l_rate * delta_weights

                # Update signs periodically using kurtosis
                if blockno % ext_params.ext_blocks == 0
                    if ext_params.kurt_size < n_samples
                        rp = view(permute_indices, 1:ext_params.kurt_size)
                        @views mul!(partact, weights, dat[:, rp])
                    else
                        @views mul!(partact, weights, dat)
                    end

                    # Compute kurtosis more efficiently
                    @views for i = 1:n_channels
                        m2[i] = mean(view(partact, i, :) .^ 2)^2
                        m4[i] = mean(view(partact, i, :) .^ 4)
                        kk[i] = (m4[i] / m2[i]) - 3.0
                    end

                    # Apply momentum to kurtosis estimates
                    if ext_params.extmomentum > 0
                        @. kk = ext_params.extmomentum * old_kurt + (1.0 - ext_params.extmomentum) * kk
                        copyto!(old_kurt, kk)
                    end

                    # Update signs and adjust update frequency
                    @. new_signs = sign(kk + ext_params.signsbias)
                    new_signs_diag = Diagonal(new_signs)

                    # Compare with oldsigns to track stability
                    if new_signs_diag == oldsigns
                        signcount += 1
                        if signcount >= ext_params.signcount_threshold
                            ext_params.ext_blocks *= ext_params.signcount_step
                            signcount = 0
                            @debug "Extending blocks to $(ext_params.ext_blocks)"
                        end
                    else
                        signcount = 0
                    end

                    oldsigns = signs  # Store current signs before updating
                    signs = new_signs_diag
                end
            else
                # Standard ICA update matching MATLAB implementation

                # Standard ICA update matching MATLAB implementation
                @. y = 1 / (1 + exp(-u))
                weights = weights .+ params.l_rate * (BI + (1 .- 2*y) * u') .* weights


                # @. y = 1 / (1 + exp(-u))
                # y .*= -2.0
                # y .+= 1.0
                # mul!(yu, y, u')
                # mul!(weights, (BI + yu), weights)

                # Add momentum if enabled
                #if params.momentum > 0
                #    weights = weights + params.momentum * delta_weights
                #end


                # Standard weight update
                # weights = weights + params.l_rate * (BI + (1 .- 2 .* y) * u') * weights

                if maximum(abs.(weights)) > params.max_weight
                    wts_blowup = true
                    change = NaN
                    break
                end
            end

            blockno += 1
        end

        if !wts_blowup
            oldwtchange = weights - oldweights
            step += 1

            # Compute and print weight and update angle changes
            lrates[step] = params.l_rate
            angledelta = 0.0
            delta = reshape(oldwtchange, 1, n_channels * n_channels)
            change = dot(delta, delta)
        end

        # Handle weight blowup
        if wts_blowup || isnan(change) || isinf(change)  # if weights blow up
            @warn "Weight blowup detected - reducing learning rate and restarting"
            step = 0
            change = params.default_stop
            wts_blowup = false
            blockno = 1
            params.l_rate *= params.restart_factor
            weights = copy(startweights)
            oldweights = copy(startweights)
            prevweights = copy(startweights)
            prevwtchange = zeros(size(weights))
            bias = zeros(size(weights, 1), 1)
            continue
        end

        # Compute weight changes matching MATLAB implementation
        oldwtchange = weights - oldweights
        delta = reshape(oldwtchange, :)
        change = dot(delta, delta)

        # if !isfinite(change)
        #     @warn "Weight change not finite at step $step"
        #     weights = copy(startweights)  # Reset to starting weights
        #     params.l_rate *= DEFAULT_RESTART_FAC

        #     # if params.l_rate < MIN_LRATE
        #     #     @warn "Learning rate too small ($params.l_rate) - terminating"
        #     #     break
        #     # end

        #     step = 0
        #     change = params.default_stop
        #     continue
        # end

        # Compute angle change (matching MATLAB)
        if step > 2
            angledelta = acos(dot(delta, olddelta) / sqrt(change * oldchange))
        end

        # Update learning rate with bounds checking
        lrates[step] = params.l_rate

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

        # Check convergence using default_stop
        if step > 2 && change < params.default_stop
            @info "Convergence reached at step $step with change = $change"
            break
        elseif change > params.blowup
            params.l_rate *= params.blowup_fac
            # params.l_rate = max(params.l_rate, MIN_LRATE)
        end

        copyto!(oldweights, weights)

        log_progress(step, change, params.l_rate, degconst * angledelta)

    end

    weights = weights * sphere * pca_components[:, 1:n_components]'
    sphere = Matrix{Float64}(I, size(weights)[2], size(weights)[2])


    # Calculate mixing matrix (topography)
    topo = pinv(weights * sphere)

    # Standardize signs based on maximum absolute value in mixing matrix
    for i = 1:size(topo, 2)
        max_abs_idx = argmax(abs.(topo[:, i]))
        if topo[max_abs_idx, i] < 0
            topo[:, i] .*= -1
            weights[i, :] .*= -1
        end
    end



    # Rescale mixing matrix if sphering was applied
    if do_sphering
        # If PCA was applied, expand scale_factors to match original data dimensions
        # if !isnothing(n_components)
        #     expanded_scale_factors = ones(size(dat, 1))
        #     expanded_scale_factors[1:n_components] = scale_factors
        # else
        #     expanded_scale_factors = scale_factors
        # end
        topo .*= scale_factors'
        weights ./= scale_factors
    end

    meanvar = sum(topo .^ 2, dims = 1) .* sum((dat_ica') .^ 2, dims = 1) ./ ((n_components * size(dat_ica)[2]) - 1)
    order = sortperm(vec(meanvar), rev = true)  # Get indices in descending order of meanvar

    # Reorder matrices and labels
    topo = topo[:, order]
    weights = weights[order, :]

    unmixing = weights * sphere
    topo = pinv(unmixing)

    # Generate labels based on actual number of components
    n_out = size(weights, 1)
    labels = ["IC$i" for i = 1:n_out]

    return InfoIca(unmixing, topo, labels)
end

"""
    pre_whiten(dat; return_scale=false)

Whiten the data by standardizing each channel.

# Arguments
- `dat::Matrix{Float64}`: Data matrix (channels × samples)
- `return_scale::Bool=false`: Whether to return scaling factors

# Returns
- If return_scale=false: `Matrix{Float64}`: Whitened data
- If return_scale=true: `Tuple{Matrix{Float64}, Vector{Float64}}`: (Whitened data, scaling factors)
"""
function pre_whiten(dat::Matrix{Float64}; return_scale::Bool = false)
    scale_factors = std(dat, dims = 2)
    whitened = dat ./ scale_factors

    return return_scale ? (whitened, scale_factors) : whitened
end

"""
    rescale_data(dat::Matrix{Float64}, scale_factors::Vector{Float64})

Rescale the data back to original scale after ICA.

# Arguments
- `dat::Matrix{Float64}`: Data matrix to rescale (channels × samples)
- `scale_factors::Vector{Float64}`: Original scaling factors from pre_whiten

# Returns
- `Matrix{Float64}`: Rescaled data
"""
function rescale_data(dat::Matrix{Float64}, scale_factors::Vector{Float64})
    return dat .* scale_factors
end




function read_mat_file(filename)
    file = matopen(filename)
    dat = read(file)
    close(file)
    return dat
end
dat = read_mat_file("dat.mat")["dat"]

# @time output = infomax_ica(dat, extended = false, n_components=68)
@time output = infomax_ica(dat, extended = false, n_components = 68)

# @time output = infomax_ica(dat, extended = false, n_components = 10)
#@time output = infomax_ica(dat, extended = true)
#@time output = infomax_ica(dat, extended = true, n_components = 10)







