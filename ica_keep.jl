using ProgressMeter
using Statistics
using Random
using LinearAlgebra
using Base.Threads
using StatsBase: kurtosis

"""
    ICAParams

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
"""
mutable struct ICAParams
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
   
    # default parameters
    function ICAParams(;
        l_rate=0.001,
        max_iter=200,
        w_change=1e-12,
        use_bias=true,
        anneal_deg=60.0,
        anneal_step=0.9,
        blowup=1e4,
        blowup_fac=0.5,
        n_small_angle=20,
        max_weight=1e8,
        restart_factor=0.9,
        degconst=180.0/π
    )
        new(l_rate, max_iter, w_change, use_bias, anneal_deg, anneal_step,
            blowup, blowup_fac, n_small_angle, max_weight, restart_factor, degconst)
    end
end

"""
    ExtendedICAParams

Structure containing extended ICA algorithm parameters.

# Fields
- `ext_blocks::Int`: Number of blocks between kurtosis estimation
- `n_subgauss::Int`: Number of sub-gaussian components
- `extmomentum::Float64`: Momentum for kurtosis averaging
- `signsbias::Float64`: Bias for sign changes
- `signcount_threshold::Int`: Threshold for sign change detection
- `signcount_step::Int`: Step size for sign change adaptation
- `kurt_size::Int`: Number of samples for kurtosis estimation
"""
mutable struct ExtendedICAParams
    ext_blocks::Int
    n_subgauss::Int
    extmomentum::Float64
    signsbias::Float64
    signcount_threshold::Int
    signcount_step::Int
    kurt_size::Int
    
    # default parameters
    function ExtendedICAParams(;
        ext_blocks=1,
        n_subgauss=1,
        extmomentum=0.5,
        signsbias=0.02,
        signcount_threshold=25,
        signcount_step=2,
        kurt_size=6000
    )
        new(ext_blocks, n_subgauss, extmomentum, signsbias,
            signcount_threshold, signcount_step, kurt_size)
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
struct InfoICA
    topo::Matrix{Float64}
    unmixing::Matrix{Float64}
    label::Vector{String}
end

"""
    infomax_ica(dat; kwargs...)

Perform Infomax ICA decomposition on EEG data using the algorithm from EEGLAB's runica.

# Arguments
- `dat::Matrix{Float64}`: Data matrix (channels × samples)
- `extended::Bool=true`: Whether to use extended Infomax
- `params::ICAParams`: Algorithm parameters
- `ext_params::ExtendedICAParams`: Extended algorithm parameters

# Returns
- `InfoICA`: Structure containing:
  - `topo`: Topography matrix (mixing matrix)
  - `unmixing`: Unmixing matrix (weights)
  - `label`: Component labels
"""
function infomax_ica(
    dat::Matrix{Float64};
    extended::Bool=true,
    params::ICAParams=ICAParams(),
    ext_params::ExtendedICAParams=ExtendedICAParams(),
)
    # Pre-allocate reused matrices and vectors
    n_channels, n_samples = size(dat)
    n_channels_square = n_channels^2
    block = Int(floor(sqrt(n_samples / 3.0)))
    nblock = div(n_samples, block)
    lastt = (nblock - 1) * block + 1

    # Pre-allocate matrices used in the loop
    weights = Matrix{Float64}(I, n_channels, n_channels)
    BI = block * Matrix{Float64}(I, n_channels, n_channels)
    bias = zeros(n_channels, 1)
    onesrow = ones(1, block)
    u = zeros(n_channels, block)
    y = similar(u)
    delta_weights = similar(weights)
    
    # For extended ICA
    if extended
        signs = ones(n_channels)
        signs[1:ext_params.n_subgauss] .= -1
        old_kurt = zeros(n_channels)
        oldsigns = copy(signs)
        tpartact = zeros(n_channels, ext_params.kurt_size)
        kurt = zeros(n_channels)
    end

    # Save initial weights
    startweights = copy(weights)
    oldweights = copy(startweights)
    
    # Initialize training variables
    step = 0
    count_small_angle = 0
    wts_blowup = false
    blockno = 0
    signcount = 0
    initial_ext_blocks = ext_params.ext_blocks
    olddelta = zeros(n_channels_square)
    oldchange = 0.0
    
    # Pre-allocate permutation vector
    permute = Vector{Int}(undef, n_samples)

    while step < params.max_iter
        # Reuse permutation vector
        randperm!(permute)

        for t in 1:block:lastt
            # Reuse pre-allocated matrices instead of creating new ones
            mul!(u, weights, view(dat, :, permute[t:t+block-1]))
            u .+= bias .* onesrow

            if extended
                # Extended ICA update
                y .= tanh.(u)
                mul!(delta_weights, weights, BI .- signs[:] .* (y * u') - (u * u'))
                weights .+= params.l_rate .* delta_weights
                
                if params.use_bias
                    bias .+= params.l_rate .* sum(y, dims=2) .* -2.0
                end
            else
                # Standard ICA update
                y .= 1 ./ (1 .+ exp.(-u))
                mul!(delta_weights, (BI + (1.0 .- 2.0 .* y) * u'), weights)
                weights .+= params.l_rate .* delta_weights
                
                if params.use_bias
                    bias .+= params.l_rate .* sum(1.0 .- 2.0 .* y, dims=2)
                end
            end

            # Check for weight blowup
            if maximum(abs.(weights)) > params.max_weight
                wts_blowup = true
                break
            end

            blockno += 1

            # Extended ICA kurtosis estimation
            if extended && ext_params.ext_blocks > 0 && 
               blockno % ext_params.ext_blocks == 0
                
                if ext_params.kurt_size < n_samples
                    rp = rand(1:n_samples, ext_params.kurt_size)
                    mul!(tpartact, weights, view(dat, :, rp))
                else
                    mul!(tpartact, weights, dat)
                end

                # Calculate kurtosis in-place
                for i in 1:n_channels
                    kurt[i] = kurtosis(view(tpartact, i, :))
                end
                
                if ext_params.extmomentum != 0
                    kurt .= ext_params.extmomentum .* old_kurt .+ 
                           (1.0 .- ext_params.extmomentum) .* kurt
                    old_kurt .= kurt
                end
                
                # Update signs in-place
                signs .= sign.(kurt .+ ext_params.signsbias)
                ndiff = sum(signs .!= oldsigns)
                
                if ndiff == 0
                    signcount += 1
                else
                    signcount = 0
                end
                
                oldsigns = signs

                if signcount >= ext_params.signcount_threshold
                    ext_params.ext_blocks = floor(
                        Int, 
                        ext_params.ext_blocks * ext_params.signcount_step
                    )
                    signcount = 0
                end
            end
        end

        if !wts_blowup
            oldwtchange = weights .- oldweights
            step += 1
            angledelta = 0.0
            delta = vec(oldwtchange)
            change = sum(delta .* delta)

            if step > 2
                angledelta = acos(sum(delta .* olddelta) / sqrt(change * oldchange))
                angledelta *= params.degconst
            end

            @info "Step $step: lrate $(params.l_rate), wchange $change, angledelta $angledelta"

            # Anneal learning rate
            oldweights = copy(weights)
            
            if angledelta > params.anneal_deg
                params.l_rate *= params.anneal_step
                olddelta = delta
                oldchange = change
                count_small_angle = 0
            else
                if step == 1
                    olddelta = delta
                    oldchange = change
                end
                
                if !isnothing(params.n_small_angle)
                    count_small_angle += 1
                    if count_small_angle > params.n_small_angle
                        params.max_iter = step
                    end
                end
            end

            # Apply stopping criteria
            if step > 2 && change < params.w_change
                step = params.max_iter
            elseif change > params.blowup
                params.l_rate *= params.blowup_fac
            end
        else
            # Restart if weights blow up
            step = 0
            wts_blowup = false
            blockno = 1
            params.l_rate *= params.restart_factor
            weights = copy(startweights)
            oldweights = copy(startweights)
            olddelta = zeros(1, n_channels_square)
            bias = zeros(n_channels, 1)
            ext_params.ext_blocks = initial_ext_blocks

            if extended
                signs = ones(n_channels)
                signs[1:ext_params.n_subgauss] .= -1
                oldsigns = zeros(n_channels)
            end

            @info "Lowering learning rate to $(params.l_rate) and restarting..."
        end
    end

    # Calculate mixing matrix (topography)
    topo = pinv(weights)
    
    # Generate default labels
    n_components = size(weights, 1)
    labels = ["IC$i" for i in 1:n_components]
    
    return InfoICA(topo, weights, labels)
end

"""
    pre_whiten(dat)

Whiten the data by standardizing each channel.

# Arguments
- `dat::Matrix{Float64}`: Data matrix (channels × samples)

# Returns
- `Matrix{Float64}`: Whitened data
"""
function pre_whiten(dat::Matrix{Float64})
    return dat ./ std(dat, dims=2)
end

"""
    pca_reduction(dat)

Perform PCA dimensionality reduction.

# Arguments
- `dat::Matrix{Float64}`: Data matrix (samples × features)

# Returns
- `Tuple{Matrix{Float64}, Vector{Float64}, Matrix{Float64}}`: 
    (components, explained variance, rotation matrix)
"""
function pca_reduction(dat::Matrix{Float64})
    # Center the data
    dat = copy(dat)
    dat .-= mean(dat, dims=1)

    # Perform SVD
    U, S, V = svd(dat, full=false)
    
    # Ensure consistent signs
    max_abs_cols = argmax(abs.(U), dims=1)
    signs = sign.(U[max_abs_cols])
    U .*= signs
    V .*= signs

    # Calculate explained variance
    explained_variance = (S .^ 2) ./ (size(dat, 1) - 1)
    
    # Scale components
    U .*= sqrt(size(dat, 1) - 1)

    return U, explained_variance, V
end

using GLMakie
using Random
using Statistics
using LinearAlgebra

# Generate synthetic EEG-like data
function generate_test_data(; 
    n_channels=32,     # Number of EEG channels
    n_samples=10000,   # Number of time points
    n_sources=4,       # Number of source signals
    sampling_rate=250  # Hz
)
    # Time vector
    t = range(0, length=n_samples, step=1/sampling_rate)
    
    # Generate source signals
    sources = zeros(n_sources, n_samples)
    
    # Source 1: Alpha wave (10 Hz oscillation)
    sources[1, :] = sin.(2π * 10 * t)
    
    # Source 2: Beta wave (20 Hz oscillation)
    sources[2, :] = 0.5 * sin.(2π * 20 * t)
    
    # Source 3: Slow drift
    sources[3, :] = 0.3 * sin.(2π * 0.1 * t)
    
    # Source 4: Random spikes (simulating muscle artifacts)
    spike_locations = rand(1:n_samples, 50)
    sources[4, spike_locations] .= randn(length(spike_locations)) * 2
    
    # Create random mixing matrix
    Random.seed!(42)  # for reproducibility
    mixing_matrix = randn(n_channels, n_sources)
    
    # Mix sources
    mixed_data = mixing_matrix * sources
    
    # Add some noise
    mixed_data .+= 0.1 * randn(size(mixed_data))
    
    return mixed_data, sources, mixing_matrix
end

# # Generate test data
# n_channels = 32
# n_samples = 10000
# data, true_sources, true_mixing = generate_test_data(
#     n_channels=n_channels,
#     n_samples=n_samples
# )
# 
# # Preprocess data
data_whitened = pre_whiten(Float64.(Matrix(test_data)))
# # Run ICA
@btime output = infomax_ica(permutedims(data_whitened), extended=false)
# 
# # Get independent components
# components = output.unmixing * data_whitened

