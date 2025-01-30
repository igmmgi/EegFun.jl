
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
    # Store original dimensions
    
    original_n_channels = size(dat, 1)

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

    # Pre-allocate workspace arrays
    weights = Matrix{Float64}(I, n_channels, n_channels)
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
    MIN_LRATE = 1e-6  # Minimum learning rate
    MAX_LRATE = 0.1   # Maximum learning rate
    DEFAULT_RESTART_FAC = 0.9  # Factor for reducing learning rate on restart

    while step < params.max_iter
        step += 1
        wts_blowup = false
        
        randperm!(permute_indices)
        
        for t = 1:block:lastt
            block_indices = permute_indices[t:min(t+block-1, end)]
            @views data_block .= dat_reduced[:, block_indices]
            
            mul!(u, weights, data_block)
            u .+= bias .* onesrow

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
                        @views mul!(partact, weights, dat_reduced[:, rp])
                    else
                        @views mul!(partact, weights, dat_reduced)
                    end
                    
                    # Compute kurtosis more efficiently
                    @views for i in 1:n_channels
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
                @. y = 1 / (1 + exp(-u))
                y .*= -2.0
                y .+= 1.0
                mul!(yu, y, u')
                mul!(delta_weights, (BI + yu), weights)
                
                if maximum(abs.(delta_weights)) > params.blowup
                    wts_blowup = true
                    break
                end
                
                @. weights += params.l_rate * delta_weights
            end

            if params.use_bias
                @views bias .+= sum(y, dims=2) .* (-2.0 * params.l_rate / block)
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

        if !isfinite(change)
            @warn "Weight change not finite at step $step"
            weights = copy(startweights)  # Reset to starting weights
            params.l_rate *= DEFAULT_RESTART_FAC
            
            # if params.l_rate < MIN_LRATE
            #     @warn "Learning rate too small ($params.l_rate) - terminating"
            #     break
            # end
            
            step = 0
            change = params.default_stop
            continue
        end

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

    # Calculate mixing matrix (topography)
    topo = pinv(weights)

    # If PCA was applied, adjust the mixing/unmixing matrices
    if !isnothing(n_components)
        # Adjust mixing matrix (original_channels × n_components)
        topo = pca_components * topo

        # Adjust unmixing matrix (n_components × original_channels)
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

    # Rescale mixing matrix if whitening was applied
    if do_whitening
        topo .*= scale_factors
        weights ./= scale_factors'
    end

    # Reorder components based on explained variance
    proj = weights * dat  # Get component activations
    vars = vec(var(proj, dims = 2))  # Calculate variance of each component
    order = sortperm(vars, rev = true)  # Get indices in descending order

    # Reorder matrices and labels
    topo = topo[:, order]
    weights = weights[order, :]

    # Generate labels based on actual number of components
    n_out = size(weights, 1)
    labels = ["IC$i" for i = 1:n_out]

    return InfoIca(topo, weights, labels)
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

@time output = infomax_ica(dat, extended = false, n_components = 5)

# Get independent components
components = output.unmixing * data_whitened





# # Generate test data
# n_channels = 32
# n_samples = 10000
# data, true_sources, true_mixing = generate_test_data(
#     n_channels=n_channels,
#     n_samples=n_samples
# )
# 
# # Preprocess data
# data_whitened = pre_whiten(Float64.(Matrix(test_data)))
# # # Run ICA
# output = infomax_ica(permutedims(data_whitened), extended=false)
# 
# # Get independent components
# components = output.unmixing * data_whitened
# 
# """
#     remove_components(data::Matrix{Float64}, ica::InfoIca, components_to_remove::Vector{Int})
# 
# Remove specified ICA components from the data and reconstruct the signal.
# 
# # Arguments
# - `data::Matrix{Float64}`: Original data matrix (channels × samples)
# - `ica::InfoIca`: ICA decomposition result
# - `components_to_remove::Vector{Int}`: Indices of components to remove (e.g., [1, 4] to remove first and fourth component)
# 
# # Returns
# - `Matrix{Float64}`: Reconstructed data with specified components removed
# """
# function remove_components(data::Matrix{Float64}, ica::InfoIca, components_to_remove::Vector{Int})
#     # Get ICA activation matrix (components × samples)
#     activation = ica.unmixing * data
# 
#     # Zero out the components we want to remove
#     activation[components_to_remove, :] .= 0
# 
#     # Reconstruct the data using the modified components
#     reconstructed_data = ica.topo * activation
# 
#     return reconstructed_data
# end
# 
# """
#     get_components(data::Matrix{Float64}, ica::InfoIca)
# 
# Get the ICA component activations.
# 
# # Arguments
# - `data::Matrix{Float64}`: Original data matrix (channels × samples)
# - `ica::InfoIca`: ICA decomposition result
# 
# # Returns
# - `Matrix{Float64}`: Component activation matrix (components × samples)
# """
# function get_components(data::Matrix{Float64}, ica::InfoIca)
#     return ica.unmixing * data
# end
# 
# """
#     plot_components(data::Matrix{Float64}, ica::InfoIca; 
#                    n_components=nothing, time_window=nothing, 
#                    sampling_rate=250)
# 
# Create plots of ICA components including time series and topographies.
# 
# # Arguments
# - `data::Matrix{Float64}`: Original data matrix (channels × samples)
# - `ica::InfoIca`: ICA decomposition result
# - `n_components::Union{Nothing,Int}=nothing`: Number of components to plot (default: all)
# - `time_window::Union{Nothing,Tuple{Int,Int}}=nothing`: Time range to plot (samples)
# - `sampling_rate::Int=250`: Sampling rate in Hz for time axis
# 
# # Returns
# - `Figure`: Makie figure containing component plots
# """
# function plot_components(
#     data::Matrix{Float64},
#     ica::InfoIca;
#     n_components = nothing,
#     time_window = nothing,
#     sampling_rate = 256,
# )
#     # Calculate component activations (n_components × n_timepoints)
#     components = get_components(data, ica)  # Using our existing get_components function
#     n_total = size(components, 1)
#     n_show = isnothing(n_components) ? n_total : min(n_components, n_total)
# 
#     # Set time window
#     if isnothing(time_window)
#         time_window = (1, size(data, 2))
#     end
#     t_start, t_end = time_window
#     time = range(t_start / sampling_rate, t_end / sampling_rate, length = t_end - t_start + 1)
# 
#     # Create figure
#     fig = Figure(resolution = (1000, 100 * n_show))
# 
#     for i = 1:n_show
#         # Time series plot
#         ax_time = Axis(fig[i, 1], xlabel = "Time (s)", ylabel = "Amplitude", title = "Component $(i)")
# 
#         lines!(ax_time, time, components[i, t_start:t_end])
# 
#         # Topography plot
#         ax_topo = Axis(fig[i, 2], aspect = 1, title = "Topography")
# 
#         # Plot topography as a simple heatmap
#         heatmap!(ax_topo, reshape(ica.topo[:, i], :, 1))
#         hidedecorations!(ax_topo)
#     end
# 
#     return fig
# end
# 
# """
#     plot_component_comparison(original::Matrix{Float64}, cleaned::Matrix{Float64};
#                             channel::Int=1, time_window=nothing, sampling_rate=250)
# 
# Create a comparison plot of original vs cleaned data for a single channel.
# 
# # Arguments
# - `original::Matrix{Float64}`: Original data matrix (channels × samples)
# - `cleaned::Matrix{Float64}`: Cleaned data matrix (channels × samples)
# - `channel::Int=1`: Channel to plot
# - `time_window::Union{Nothing,Tuple{Int,Int}}=nothing`: Time range to plot (samples)
# - `sampling_rate::Int=250`: Sampling rate in Hz for time axis
# 
# # Returns
# - `Figure`: Makie figure containing comparison plot
# """
# function plot_component_comparison(
#     original::Matrix{Float64},
#     cleaned::Matrix{Float64};
#     channel::Int = 1,
#     time_window = nothing,
#     sampling_rate = 256,
# )
#     # Set time window
#     if isnothing(time_window)
#         time_window = (1, size(original, 2))
#     end
#     t_start, t_end = time_window
#     time = range(t_start / sampling_rate, t_end / sampling_rate, length = t_end - t_start + 1)
# 
#     # Create figure
#     fig = Figure(resolution = (1000, 400))
# 
#     ax = Axis(fig[1, 1], xlabel = "Time (s)", ylabel = "Amplitude", title = "Channel $channel: Original vs Cleaned")
# 
#     lines!(ax, time, original[channel, t_start:t_end], label = "Original")
#     lines!(ax, time, cleaned[channel, t_start:t_end], label = "Cleaned")
#     axislegend()
# 
#     return fig
# end
# 
# # Plot all components
# # fig1 = plot_components(permutedims(data_whitened), output, n_components=5)
# 
# # # Plot first 5 components with a specific time window
# # fig2 = plot_components(data_whitened, output, 
# #     n_components=5, 
# #     time_window=(1, 1000),  # first 4 seconds at 250 Hz
# #     sampling_rate=250
# # )
# 
# # Compare original vs cleaned data for frontal channel
# # fig3 = plot_component_comparison(data_whitened, cleaned_whitened,
# #     channel=1,  # typically Fp1 or similar for blink artifacts
# #     time_window=(1, 1000),
# #     sampling_rate=250
# # )
# 
# 
# 
# # Preprocess with whitening, keeping track of scaling factors
# # data_whitened, scale_factors = pre_whiten(data, return_scale=true)
# # 
# # # Run ICA
# # ica_result = infomax_ica(data_whitened)
# # 
# # # Remove components
# # cleaned_whitened = remove_components(data_whitened, ica_result, [1, 3])
# # 
# # # Rescale back to original units
# # cleaned_data = rescale_data(cleaned_whitened, scale_factors)







