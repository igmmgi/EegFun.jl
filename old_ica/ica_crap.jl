# using CairoMakie
using Base.Threads
using BioSemiBDF
using CSV
using DSP
using DataFrames
using GLMakie
using GLMakie
using JLD2
using LibGEOS
using LinearAlgebra
using LinearAlgebra
using LinearAlgebra
using MAT
using OrderedCollections
using Printf
using Random
using Random
using Random
using ScatteredInterpolation
using Statistics
using StatsBase
using StatsBase: kurtosis


"""
    IcaPrms

Structure containing ICA algorithm parameters.

# Fields
- `l_rate::Float64`: Initial learning rate
- `max_iter::Int`: Maximum number of iterations
- `w_change::Float64`: Change threshold for stopping
- `use_bias::Bool`: Whether to use bias term
- `anneal_deg::Float64`: Angle for learning rate reduction
- `anneal_step::Float64`: Learning rate reduction factor
- `blowup::Float64`: Maximum weight change allowed
- `blowup_fac::Float64`: Learning rate reduction on blowup
- `n_small_angle::Int`: Number of angles below threshold
- `max_weight::Float64`: Maximum weight magnitude
- `restart_factor::Float64`: Learning rate factor on restart
- `degconst::Float64`: Degrees to radians conversion
- `default_stop::Float64`: Default weight stop criterion
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

    # default parameters
    function IcaPrms(;
        l_rate = 0.001,
        max_iter = 512,
        w_change = 0.000001,
        use_bias = true,
        anneal_deg = 60.0,
        anneal_step = 0.9,
        blowup = 1e4,
        blowup_fac = 0.8,
        n_small_angle = 20,
        max_weight = 1e8,
        restart_factor = 0.9,
        degconst = 180.0 / π,
        default_stop = 0.000001,
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
    validate_inputs(data::Matrix{Float64}, params::IcaPrms, ext_params::Union{Nothing, ExtendedIcaPrms})

Validate inputs and parameters for ICA computation.

# Arguments
- `data`: Input data matrix (channels × samples)
- `params`: ICA parameters
- `ext_params`: Extended ICA parameters (if applicable)

# Throws
- `ArgumentError` if inputs are invalid
"""
function validate_inputs(data::Matrix{Float64}, params::IcaPrms, ext_params::Union{Nothing, ExtendedIcaPrms})
    # Check data dimensions
    n_channels, n_samples = size(data)
    if n_channels < 2
        throw(ArgumentError("Data must have at least 2 channels"))
    end
    if n_samples < n_channels
        throw(ArgumentError("Number of samples must be greater than or equal to number of channels"))
    end

    # Check for NaN/Inf values
    if any(!isfinite, data)
        throw(ArgumentError("Data contains NaN or Inf values"))
    end

    # Validate ICA parameters
    if params.l_rate <= 0
        throw(ArgumentError("Learning rate must be positive"))
    end
    if params.max_iter < 1
        throw(ArgumentError("Maximum iterations must be at least 1"))
    end
    if params.default_stop <= 0
        throw(ArgumentError("Convergence threshold must be positive"))
    end

    # Validate extended ICA parameters
    if ext_params !== nothing
        if ext_params.n_subgauss < 0 || ext_params.n_subgauss > n_channels
            throw(ArgumentError("Number of sub-Gaussian components must be between 0 and $n_channels"))
        end
        if ext_params.extmomentum < 0 || ext_params.extmomentum > 1
            throw(ArgumentError("Momentum must be between 0 and 1"))
        end
    end
end






"""
    log_progress(step::Int, max_iter::Int, change::Float64, l_rate::Float64, angledelta::Float64)
Log progress of ICA computation.
"""
function log_progress(step::Int, max_iter::Int, change::Float64, l_rate::Float64, angledelta::Float64)
    if step % 10 == 0 || step == 1 || step == max_iter
        progress = 100 * step / max_iter
        @info "Step $step ($(round(progress, digits=2))%): change = $change, lrate = $l_rate, angle = $angledelta"
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
    topo::Matrix{Float64}
    unmixing::Matrix{Float64}
    label::Vector{String}
end

"""
    IcaWorkspace

Structure containing pre-allocated arrays for ICA computation.

# Fields
- `weights`: Current weight matrix
- `oldweights`: Previous weight matrix
- `startweights`: Initial weight matrix
- `bias`: Bias vector
- `data_block`: Current data block
- `u`: Activations
- `y`: Nonlinearity output
- `delta_weights`: Weight updates
- `temp_matrix`: Temporary matrix for computations
- `yu`: y * u'
- `uu`: u * u'
- `BI`: Block size * identity matrix
- `onesrow`: Row vector of ones
- `permute_indices`: Permutation indices for data blocks
- `olddelta`: Previous weight change vector
- `degconst`: Degrees to radians conversion factor
- `signs`: Signs for extended ICA (Diagonal matrix)
- `old_kurt`: Previous kurtosis values for extended ICA
- `oldsigns`: Previous signs for extended ICA
- `partact`: Partial activations for kurtosis estimation
- `kk`: Current kurtosis values
- `m2`: Second moments
- `m4`: Fourth moments
- `new_signs`: New signs for extended ICA
"""
struct IcaWorkspace
    weights::Matrix{Float64}
    oldweights::Matrix{Float64}
    startweights::Matrix{Float64}
    bias::Vector{Float64}
    data_block::Matrix{Float64}
    u::Matrix{Float64}
    y::Matrix{Float64}
    delta_weights::Matrix{Float64}
    temp_matrix::Matrix{Float64}
    yu::Matrix{Float64}
    uu::Matrix{Float64}
    BI::Matrix{Float64}
    onesrow::Matrix{Float64}
    permute_indices::Vector{Int}
    olddelta::Vector{Float64}
    degconst::Float64
    
    # Extended ICA specific
    signs::Union{Nothing, Diagonal{Float64}}
    old_kurt::Union{Nothing, Vector{Float64}}
    oldsigns::Union{Nothing, Diagonal{Float64}}
    partact::Union{Nothing, Matrix{Float64}}
    kk::Union{Nothing, Vector{Float64}}
    m2::Union{Nothing, Vector{Float64}}
    m4::Union{Nothing, Vector{Float64}}
    new_signs::Union{Nothing, Vector{Float64}}

    function IcaWorkspace(n_channels::Int, n_samples::Int, block::Int, extended::Bool, ext_params::Union{Nothing, ExtendedIcaPrms})
        # Initialize weights with small random values
        weights = 0.01 * randn(n_channels, n_channels)
        weights ./= norm(weights)

        # Initialize other arrays
        oldweights = copy(weights)
        startweights = copy(weights)
        bias = zeros(n_channels)
        data_block = zeros(n_channels, block)
        u = zeros(n_channels, block)
        y = similar(u)
        delta_weights = similar(weights)
        temp_matrix = similar(weights)
        yu = similar(weights)
        uu = similar(weights)
        BI = block * Matrix{Float64}(I, n_channels, n_channels)
        onesrow = ones(1, block)
        permute_indices = Vector{Int}(undef, n_samples)

        # Initialize extended ICA arrays if needed
        if extended
            signs = Diagonal(vcat(-ones(ext_params.n_subgauss), ones(n_channels - ext_params.n_subgauss)))
            old_kurt = zeros(n_channels)
            oldsigns = copy(signs)
            partact = zeros(n_channels, min(ext_params.kurt_size, n_samples))
            kk = zeros(n_channels)
            m2 = zeros(n_channels)
            m4 = zeros(n_channels)
            new_signs = zeros(n_channels)
        else
            signs = nothing
            old_kurt = nothing
            oldsigns = nothing
            partact = nothing
            kk = nothing
            m2 = nothing
            m4 = nothing
            new_signs = nothing
        end

        return IcaWorkspace(weights, oldweights, startweights, bias, data_block, u, y, delta_weights, temp_matrix, yu, uu, BI, onesrow, permute_indices, signs, old_kurt, oldsigns, partact, kk, m2, m4, new_signs)
    end
end

"""
    initialize_extended_ica!(ws::IcaWorkspace, n_channels::Int, n_samples::Int, ext_params::ExtendedIcaPrms)

Initialize extended ICA-specific arrays.
"""
function initialize_extended_ica!(ws::IcaWorkspace, n_channels::Int, n_samples::Int, ext_params::ExtendedIcaPrms)
    ws.signs = Diagonal(vcat(-ones(ext_params.n_subgauss), ones(n_channels - ext_params.n_subgauss)))
    ws.old_kurt = zeros(n_channels)
    ws.oldsigns = copy(ws.signs)
    ws.partact = zeros(n_channels, min(ext_params.kurt_size, n_samples))
    ws.kk = zeros(n_channels)
    ws.m2 = zeros(n_channels)
    ws.m4 = zeros(n_channels)
    ws.new_signs = zeros(n_channels)
end

"""
    initialize_workspace(data, extended, ext_params)

Initialize workspace arrays for ICA computation.
"""
function initialize_workspace(data::Matrix{Float64}, extended::Bool, ext_params::Union{Nothing, ExtendedIcaPrms})
    n_channels, n_samples = size(data)
    block = Int(floor(sqrt(n_samples / 3.0)))
    return IcaWorkspace(n_channels, n_samples, block, extended, ext_params)
end

"""
    infomax_ica(dat; kwargs...)

Perform Infomax ICA decomposition on EEG data using the algorithm from EEGLAB's runica.

# Arguments
- `dat::Matrix{Float64}`: Data matrix (channels × samples)
- `extended::Bool=true`: Whether to use extended Infomax
- `params::ICAParams`: Algorithm parameters
- `ext_params::ExtendedICAParams`: Extended algorithm parameters
- `n_components::Union{Nothing,Int}=nothing`: Number of components for PCA reduction (default: no reduction)
- `do_whitening::Bool=true`: Whether to perform whitening

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
    do_whitening::Bool = true,
)
    # Validate inputs
    validate_inputs(dat, params, ext_params)

    # Store original dimensions
    
    original_n_channels = size(dat, 1)

    # Scale data to unit variance
    data_scaled = dat ./ std(dat, dims=2)

    # Apply whitening if requested
    if do_whitening
        data_whitened, scale_factors = pre_whiten(dat, return_scale = true)
    else
        data_whitened = dat
        scale_factors = ones(size(dat, 1))
    end

    # Apply PCA reduction if n_components is specified
    if !isnothing(n_components)
        # Ensure n_components doesn't exceed data dimensions
        # n_components = min(n_components, minimum(size(data_whitened)))

        # Center the data
        data_centered = data_whitened .- mean(data_whitened, dims = 2)

        # Compute SVD directly on centered data for efficiency
        F = svd(data_centered)

        # Take only the first n_components
        dat_reduced = F.U[:, 1:n_components]' * data_centered

        # Store PCA components for later reconstruction
        pca_components = F.U[:, 1:n_components]
    else
        dat_reduced = data_whitened
        pca_components = nothing
    end

    # Pre-allocate all matrices and temporaries
    n_channels = size(dat_reduced, 1)
    n_samples = size(dat_reduced, 2)
    n_channels_square = n_channels^2
    block = Int(floor(sqrt(n_samples / 3.0)))
    nblock = div(n_samples, block)
    lastt = (nblock - 1) * block + 1

    # Initialize workspace
    workspace = IcaWorkspace(n_channels, n_samples, block, extended, ext_params)
    weights = workspace.weights
    oldweights = workspace.oldweights
    startweights = workspace.startweights
    bias = workspace.bias
    data_block = workspace.data_block
    u = workspace.u
    y = workspace.y
    delta_weights = workspace.delta_weights
    temp_matrix = workspace.temp_matrix
    yu = workspace.yu
    uu = workspace.uu
    BI = workspace.BI
    onesrow = workspace.onesrow
    permute_indices = workspace.permute_indices

    # Extended ICA variables
    if extended
        signs = workspace.signs
        old_kurt = workspace.old_kurt
        oldsigns = workspace.oldsigns
        partact = workspace.partact
        kk = workspace.kk
        m2 = workspace.m2
        m4 = workspace.m4
        new_signs = workspace.new_signs
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
    MIN_LRATE = 1e-6  # Minimum learning rate
    MAX_LRATE = 0.1   # Maximum learning rate
    DEFAULT_RESTART_FAC = 0.9  # Factor for reducing learning rate on restart

    # Set initial learning rate based on data size
    params.l_rate = min(0.001, 1.0 / sqrt(n_samples))

    while step < params.max_iter
        step += 1

        # Log progress
        log_progress(step, params.max_iter, change, params.l_rate, degconst * angledelta)

        wts_blowup = false
        
        randperm!(permute_indices)
        
        for t = 1:block:lastt
            process_block!(workspace, dat_reduced, t, block, params, extended ? ext_params : nothing)
            
            if maximum(abs.(workspace.delta_weights)) > params.blowup
                wts_blowup = true
                break
            end
            
            blockno += 1
        end

        # Handle weight blowup
        if wts_blowup
            @warn "Weight change too large at step $step - reducing learning rate"
            weights = copy(startweights)  # Reset to starting weights
            params.l_rate *= DEFAULT_RESTART_FAC
            
            if params.l_rate < MIN_LRATE
                @warn "Learning rate too small ($params.l_rate) - terminating"
                break
            end
            
            step = 0  # Reset step count
            change = params.default_stop
            wts_blowup = false
            blockno = 1
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
        #     
        #     # if params.l_rate < MIN_LRATE
        #     #     @warn "Learning rate too small ($params.l_rate) - terminating"
        #     # end
        #     
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
                # params.l_rate = max(params.l_rate, MIN_LRATE)  # Don't go below minimum
                # params.l_rate = min(params.l_rate, MAX_LRATE)  # Don't go above maximum
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
            params.l_rate = max(params.l_rate, MIN_LRATE)
        end
        
        copyto!(oldweights, weights)
        
        @info "Step $step: change = $change, lrate = $(params.l_rate), angle = $(degconst * angledelta)"
    end

    # Post-process results
    topo, weights, labels = postprocess_results(weights, dat_reduced, pca_components, do_whitening ? scale_factors : nothing)

    return InfoIca(topo, weights, labels)
end

"""
    postprocess_results(weights::Matrix{Float64}, data::Matrix{Float64}, pca_components::Union{Nothing, Matrix{Float64}}, scale_factors::Union{Nothing, Vector{Float64}})

Post-process ICA results including sign standardization, component reordering, and matrix adjustments.

# Arguments
- `weights`: Unmixing matrix
- `data`: Processed data matrix
- `pca_components`: PCA components (if dimensionality reduction was applied)
- `scale_factors`: Whitening scale factors (if whitening was applied)

# Returns
- `topo`: Mixing matrix (channels × components)
- `weights`: Unmixing matrix (components × channels)
- `labels`: Component labels
"""
function postprocess_results(weights::Matrix{Float64}, data::Matrix{Float64}, pca_components::Union{Nothing, Matrix{Float64}}, scale_factors::Union{Nothing, Vector{Float64}})
    # Calculate mixing matrix (topography)
    topo = pinv(weights)

    # Adjust for PCA if applied
    if !isnothing(pca_components)
        topo = pca_components * topo
        weights = weights * pca_components'
    end

    # Standardize signs based on maximum absolute value in mixing matrix
    for i = 1:size(topo, 2)
        max_abs_idx = argmax(abs.(topo[:, i]))
        if topo[max_abs_idx, i] < 0
            topo[:, i] .*= -1
            weights[i, :] .*= -1
        end
    end

    # Adjust for whitening if applied
    if !isnothing(scale_factors)
        topo .*= scale_factors
        weights ./= scale_factors'
    end

    # Reorder components based on explained variance
    proj = weights * data
    vars = vec(var(proj, dims=2))
    order = sortperm(vars, rev=true)

    # Reorder matrices
    topo = topo[:, order]
    weights = weights[order, :]

    # Generate labels
    n_out = size(weights, 1)
    labels = ["IC$i" for i = 1:n_out]

    return topo, weights, labels
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

function process_block!(workspace::IcaWorkspace, data::Matrix{Float64}, t::Int, block::Int, params::IcaPrms, ext_params::Union{Nothing, ExtendedIcaPrms})
    @views block_indices = workspace.permute_indices[t:min(t+block-1, end)]
    @views copyto!(workspace.data_block, data[:, block_indices])
    
    # Compute activations
    @inbounds mul!(workspace.u, workspace.weights, workspace.data_block)
    @inbounds workspace.u .+= workspace.bias .* workspace.onesrow

    # Check for NaN/Inf in activations
    if any(!isfinite, workspace.u)
        @warn "Activations contain NaN/Inf values"
    end

    if ext_params !== nothing
        # Extended ICA update
        @inbounds @. workspace.y = tanh(workspace.u)
        @inbounds mul!(workspace.yu, workspace.y, workspace.u')
        @inbounds mul!(workspace.uu, workspace.u, workspace.u')
        
        @inbounds @. workspace.temp_matrix = workspace.BI - workspace.signs * workspace.yu - workspace.uu
        @inbounds mul!(workspace.delta_weights, workspace.temp_matrix, workspace.weights)
    else
        # Standard ICA update
        @inbounds @. workspace.y = 1 / (1 + exp(-workspace.u))
        @inbounds @. workspace.y = workspace.y * -2.0 + 1.0
        @inbounds mul!(workspace.yu, workspace.y, workspace.u')
        @inbounds mul!(workspace.delta_weights, workspace.BI + workspace.yu, workspace.weights)
    end

    # Check for NaN/Inf in weight updates
    if any(!isfinite, workspace.delta_weights)
        @warn "Weight updates contain NaN/Inf values"
    end

    # Update weights (in-place)
    @inbounds @. workspace.weights += params.l_rate * workspace.delta_weights

    # Update bias if needed
    if params.use_bias
        @views workspace.bias .+= sum(workspace.y, dims=2) .* (-2.0 * params.l_rate / block)
    end
end

# Generate synthetic EEG-like data
# function generate_test_data(;
#     n_channels = 32,     # Number of EEG channels
#     n_samples = 10000,   # Number of time points
#     n_sources = 4,       # Number of source signals
#     sampling_rate = 250,  # Hz
# )
#     # Time vector
#     t = range(0, length = n_samples, step = 1 / sampling_rate)
# 
#     # Generate source signals
#     sources = zeros(n_sources, n_samples)
# 
#     # Source 1: Alpha wave (10 Hz oscillation)
#     sources[1, :] = sin.(2π * 10 * t)
# 
#     # Source 2: Beta wave (20 Hz oscillation)
#     sources[2, :] = 0.5 * sin.(2π * 20 * t)
# 
#     # Source 3: Slow drift
#     sources[3, :] = 0.3 * sin.(2π * 0.1 * t)
# 
#     # Source 4: Random spikes (simulating muscle artifacts)
#     spike_locations = rand(1:n_samples, 50)
#     sources[4, spike_locations] .= randn(length(spike_locations)) * 2
# 
#     # Create random mixing matrix
#     Random.seed!(42)  # for reproducibility
#     mixing_matrix = randn(n_channels, n_sources)
# 
#     # Mix sources
#     mixed_data = mixing_matrix * sources
# 
#     # Add some noise
#     mixed_data .+= 0.1 * randn(size(mixed_data))
# 
#     return mixed_data, sources, mixing_matrix
# end

function read_mat_file(filename)
    file = matopen(filename)
    dat = read(file)
    close(file)
    return dat
end
dat = read_mat_file("dat.mat")["dat"]

@time output = infomax_ica(dat, extended = false)
@time output = infomax_ica(dat, extended = false, n_components = 10)
@time output = infomax_ica(dat, extended = true)
@time output = infomax_ica(dat, extended = true, n_components = 10)


"""
    benchmark_ica(data::Matrix{Float64}; kwargs...)

Benchmark the ICA computation and return timing results.
"""
function benchmark_ica(data::Matrix{Float64}; kwargs...)
    @info "Running ICA benchmark..."
    results = @timed infomax_ica(data; kwargs...)
    @info "Benchmark completed: time = $(results.time) seconds, memory = $(results.bytes / 1e6) MB"
    return results
end

"""
    test_ica()

Run unit tests for ICA functionality.
"""
function test_ica()
    @info "Running ICA unit tests..."
    
    # Generate synthetic data
    data = rand(32, 1000)  # 32 channels, 1000 time points
    
    # Test standard ICA
    result = infomax_ica(data, extended=false)
    @assert size(result.topo) == (32, 32) "Mixing matrix size incorrect"
    @assert size(result.unmixing) == (32, 32) "Unmixing matrix size incorrect"
    
    # Test extended ICA
    result = infomax_ica(data, extended=true)
    @assert size(result.topo) == (32, 32) "Mixing matrix size incorrect"
    
    @info "All tests passed!"
end

function check_convergence(workspace::IcaWorkspace, params::IcaPrms, step::Int, change::Float64, oldchange::Float64)
    # Compute weight change
    oldwtchange = workspace.weights - workspace.oldweights
    delta = vec(oldwtchange)
    
    # Check for NaN/Inf in delta
    if any(!isfinite, delta)
        @warn "Weight change contains NaN/Inf values - resetting weights"
        copyto!(workspace.weights, workspace.startweights)
        params.l_rate *= 0.9
        return false, change  # Indicate reset
    end
    
    # Compute change
    change = dot(delta, delta)
    
    # Check for NaN/Inf in change
    if !isfinite(change)
        @warn "Change value is NaN/Inf - resetting weights"
        copyto!(workspace.weights, workspace.startweights)
        params.l_rate *= 0.9
        return false, change  # Indicate reset
    end
    
    # Check convergence
    if step > 2 && change < params.default_stop
        @info "Convergence reached at step $step with change = $change"
        return true, change  # Indicate convergence
    end
    
    return false, change  # Continue training
end






