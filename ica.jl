using Statistics
using Random
using LinearAlgebra
using Base.Threads
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
   
    # default parameters
    function IcaPrms(;
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
    ExtendedIcaPrms

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
mutable struct ExtendedIcaPrms
    ext_blocks::Int
    n_subgauss::Int
    extmomentum::Float64
    signsbias::Float64
    signcount_threshold::Int
    signcount_step::Int
    kurt_size::Int
    
    # default parameters
    function ExtendedIcaPrms(;
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

# Returns
- `InfoICA`: Structure containing:
  - `topo`: Topography matrix (mixing matrix)
  - `unmixing`: Unmixing matrix (weights)
  - `label`: Component labels
"""
function infomax_ica(
    dat::Matrix{Float64};
    extended::Bool=true,
    params::IcaPrms=IcaPrms(),
    ext_params::ExtendedIcaPrms=ExtendedIcaPrms(),
)
    # Pre-allocate all matrices and temporaries
    n_channels, n_samples = size(dat)
    n_channels_square = n_channels^2
    block = Int(floor(sqrt(n_samples / 3.0)))
    nblock = div(n_samples, block)
    lastt = (nblock - 1) * block + 1

    # Pre-allocate all matrices used in the loop
    weights = Matrix{Float64}(I, n_channels, n_channels)
    BI = block * Matrix{Float64}(I, n_channels, n_channels)
    bias = zeros(n_channels, 1)
    onesrow = ones(1, block)
    u = zeros(n_channels, block)
    y = similar(u)
    delta_weights = similar(weights)
    temp_matrix = similar(weights)  # Temporary matrix for intermediate calculations
    yu = similar(weights)  # Pre-allocate y * u' matrix
    uu = similar(weights)  # Pre-allocate u * u' matrix
    
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
    
    # Pre-allocate permutation vector and reuse it
    permute = Vector{Int}(undef, n_samples)
    
    # Pre-allocate vectors for kurtosis calculation
    if extended
        rp = Vector{Int}(undef, ext_params.kurt_size)
    end

    while step < params.max_iter
        randperm!(permute)

        for t in 1:block:lastt
            dat_view = view(dat, :, permute[t:t+block-1])
            mul!(u, weights, dat_view)
            u .+= bias .* onesrow

            if extended
                # Extended ICA update - minimize temporary allocations
                y .= tanh.(u)
                mul!(yu, y, u')
                mul!(uu, u, u')
                mul!(temp_matrix, weights, BI)
                temp_matrix .-= signs[:] .* yu
                temp_matrix .-= uu
                mul!(delta_weights, weights, temp_matrix)
                @. weights += params.l_rate * delta_weights
                
                if params.use_bias
                    @views bias .+= params.l_rate .* sum(y, dims=2) .* -2.0
                end
            else
                # Standard ICA update - minimize temporary allocations
                y .= 1 ./ (1 .+ exp.(-u))
                mul!(yu, (1.0 .- 2.0 .* y), u')
                mul!(temp_matrix, BI + yu, weights)
                @. weights += params.l_rate * temp_matrix
                
                if params.use_bias
                    @views bias .+= params.l_rate .* sum(1.0 .- 2.0 .* y, dims=2)
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
function pre_whiten(dat::Matrix{Float64}; return_scale::Bool=false)
    scale_factors = std(dat, dims=2)
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
# data_whitened = pre_whiten(Float64.(Matrix(test_data)))
# # # Run ICA
# output = infomax_ica(permutedims(data_whitened), extended=false)
# 
# # Get independent components
# components = output.unmixing * data_whitened

"""
    remove_components(data::Matrix{Float64}, ica::InfoIca, components_to_remove::Vector{Int})

Remove specified ICA components from the data and reconstruct the signal.

# Arguments
- `data::Matrix{Float64}`: Original data matrix (channels × samples)
- `ica::InfoIca`: ICA decomposition result
- `components_to_remove::Vector{Int}`: Indices of components to remove (e.g., [1, 4] to remove first and fourth component)

# Returns
- `Matrix{Float64}`: Reconstructed data with specified components removed
"""
function remove_components(
    data::Matrix{Float64}, 
    ica::InfoIca, 
    components_to_remove::Vector{Int}
)
    # Get ICA activation matrix (components × samples)
    activation = ica.unmixing * data
    
    # Zero out the components we want to remove
    activation[components_to_remove, :] .= 0
    
    # Reconstruct the data using the modified components
    reconstructed_data = ica.topo * activation
    
    return reconstructed_data
end

"""
    get_components(data::Matrix{Float64}, ica::InfoIca)

Get the ICA component activations.

# Arguments
- `data::Matrix{Float64}`: Original data matrix (channels × samples)
- `ica::InfoIca`: ICA decomposition result

# Returns
- `Matrix{Float64}`: Component activation matrix (components × samples)
"""
function get_components(data::Matrix{Float64}, ica::InfoIca)
    return ica.unmixing * data
end

"""
    plot_components(data::Matrix{Float64}, ica::InfoIca; 
                   n_components=nothing, time_window=nothing, 
                   sampling_rate=250)

Create plots of ICA components including time series and topographies.

# Arguments
- `data::Matrix{Float64}`: Original data matrix (channels × samples)
- `ica::InfoIca`: ICA decomposition result
- `n_components::Union{Nothing,Int}=nothing`: Number of components to plot (default: all)
- `time_window::Union{Nothing,Tuple{Int,Int}}=nothing`: Time range to plot (samples)
- `sampling_rate::Int=250`: Sampling rate in Hz for time axis

# Returns
- `Figure`: Makie figure containing component plots
"""
function plot_components(
    data::Matrix{Float64}, 
    ica::InfoIca;
    n_components=nothing,
    time_window=nothing,
    sampling_rate=256
)
    # Calculate component activations (n_components × n_timepoints)
    components = get_components(data, ica)  # Using our existing get_components function
    n_total = size(components, 1)
    n_show = isnothing(n_components) ? n_total : min(n_components, n_total)
    
    # Set time window
    if isnothing(time_window)
        time_window = (1, size(data, 2))
    end
    t_start, t_end = time_window
    time = range(t_start/sampling_rate, t_end/sampling_rate, length=t_end-t_start+1)
    
    # Create figure
    fig = Figure(resolution=(1000, 100 * n_show))
    
    for i in 1:n_show
        # Time series plot
        ax_time = Axis(fig[i, 1],
            xlabel = "Time (s)",
            ylabel = "Amplitude",
            title = "Component $(i)")
        
        lines!(ax_time, time, components[i, t_start:t_end])
        
        # Topography plot
        ax_topo = Axis(fig[i, 2],
            aspect=1,
            title = "Topography")
        
        # Plot topography as a simple heatmap
        heatmap!(ax_topo, reshape(ica.topo[:, i], :, 1))
        hidedecorations!(ax_topo)
    end
    
    return fig
end

"""
    plot_component_comparison(original::Matrix{Float64}, cleaned::Matrix{Float64};
                            channel::Int=1, time_window=nothing, sampling_rate=250)

Create a comparison plot of original vs cleaned data for a single channel.

# Arguments
- `original::Matrix{Float64}`: Original data matrix (channels × samples)
- `cleaned::Matrix{Float64}`: Cleaned data matrix (channels × samples)
- `channel::Int=1`: Channel to plot
- `time_window::Union{Nothing,Tuple{Int,Int}}=nothing`: Time range to plot (samples)
- `sampling_rate::Int=250`: Sampling rate in Hz for time axis

# Returns
- `Figure`: Makie figure containing comparison plot
"""
function plot_component_comparison(
    original::Matrix{Float64},
    cleaned::Matrix{Float64};
    channel::Int=1,
    time_window=nothing,
    sampling_rate=256
)
    # Set time window
    if isnothing(time_window)
        time_window = (1, size(original, 2))
    end
    t_start, t_end = time_window
    time = range(t_start/sampling_rate, t_end/sampling_rate, length=t_end-t_start+1)
    
    # Create figure
    fig = Figure(resolution=(1000, 400))
    
    ax = Axis(fig[1, 1],
        xlabel = "Time (s)",
        ylabel = "Amplitude",
        title = "Channel $channel: Original vs Cleaned")
    
    lines!(ax, time, original[channel, t_start:t_end], label="Original")
    lines!(ax, time, cleaned[channel, t_start:t_end], label="Cleaned")
    axislegend()
    
    return fig
end

# Plot all components
# fig1 = plot_components(permutedims(data_whitened), output, n_components=5)

# # Plot first 5 components with a specific time window
# fig2 = plot_components(data_whitened, output, 
#     n_components=5, 
#     time_window=(1, 1000),  # first 4 seconds at 250 Hz
#     sampling_rate=250
# )

# Compare original vs cleaned data for frontal channel
# fig3 = plot_component_comparison(data_whitened, cleaned_whitened,
#     channel=1,  # typically Fp1 or similar for blink artifacts
#     time_window=(1, 1000),
#     sampling_rate=250
# )



# Preprocess with whitening, keeping track of scaling factors
# data_whitened, scale_factors = pre_whiten(data, return_scale=true)
# 
# # Run ICA
# ica_result = infomax_ica(data_whitened)
# 
# # Remove components
# cleaned_whitened = remove_components(data_whitened, ica_result, [1, 3])
# 
# # Rescale back to original units
# cleaned_data = rescale_data(cleaned_whitened, scale_factors)