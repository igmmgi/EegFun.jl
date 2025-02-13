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
    default_stop = 1e-6
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
        default_stop
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




function infomax_ica(
    dat::ContinuousData,
    data_labels;
    n_components::Union{Nothing,Int} = nothing,
    params::IcaPrms = IcaPrms(),
)
    # select actual eeg data columns
    dat = permutedims(Matrix(dat.data[!, intersect(names(dat.data), data_labels)]))
    return infomax_ica(dat, data_labels, params = params, n_components = n_components)
end


function infomax_ica(
    dat::Matrix{Float64},
    data_labels;
    n_components::Union{Nothing,Int} = nothing,
    params::IcaPrms = IcaPrms(),
)
    # Preprocessing - minimize copies
    dat_ica = demean(dat)
    scale = sqrt(norm((dat_ica * dat_ica') / size(dat_ica, 2)))
    dat_ica ./= scale

    # PCA reduction - optimize memory usage
    if !isnothing(n_components)
        F = svd(dat_ica)
        dat_ica = F.U[:, 1:n_components]' * dat_ica
        pca_components = view(F.U, :, 1:n_components)
    else
        n_components = size(dat_ica, 1)
        pca_components = nothing
    end

    # Sphering - optimize computation
    sphere = 2.0 * inv(sqrt(cov(dat_ica, dims=2)))
    mul!(dat_ica, sphere, dat_ica)

    # Pre-allocate with exact sizes
    n_channels = size(dat_ica, 1)
    n_samples = size(dat_ica, 2)
    block = min(Int(floor(sqrt(n_samples / 3.0))), 512) # Cap block size
    lastt = block * div(n_samples, block)

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

    work = WorkArrays(
        Matrix{Float64}(I, n_components, n_components),  # weights
        block * Matrix{Float64}(I, n_channels, n_channels),  # BI
        zeros(n_channels, block),  # u
        zeros(n_channels, block),  # y
        zeros(n_channels, block),  # data_block
        copy(Matrix{Float64}(I, n_components, n_components)),  # oldweights
        copy(Matrix{Float64}(I, n_components, n_components)),  # startweights
        similar(Matrix{Float64}(I, n_components, n_components)),  # weights_temp
        zeros(n_channels, block),  # y_temp
        similar(Matrix{Float64}(I, n_components, n_components)),  # bi_weights
        similar(Matrix{Float64}(I, n_components, n_components)),  # wu_term
        zeros(1, n_channels^2),  # delta
        zeros(1, n_channels^2)   # olddelta
    )

    # Pre-allocate permutation vector
    permute_indices = Vector{Int}(undef, n_samples)

    # Training variables
    step = 0
    wts_blowup = false
    change = 0.0
    oldchange = 0.0
    angledelta = 0.0

    # Main loop with optimizations
    @inbounds while step < params.max_iter
        randperm!(permute_indices)

        @inbounds for t = 1:block:lastt
            block_end = min(t + block - 1, n_samples)
            block_size = block_end - t + 1
            
            # Use direct indexing instead of views where possible
            @simd for i in 1:n_channels
                @simd for j in 1:block_size
                    work.data_block[i,j] = dat_ica[i, permute_indices[t+j-1]]
                end
            end
            
            # Optimized matrix operations
            mul!(work.u, work.weights, work.data_block)
            
            # Vectorized operations
            @simd for i in eachindex(work.u)
                work.y[i] = 1 / (1 + exp(-work.u[i]))
                work.y_temp[i] = 1 - 2work.y[i]
            end

            mul!(work.wu_term, work.y_temp, work.u')
            
            @simd for i in eachindex(work.bi_weights)
                work.bi_weights[i] = work.BI[i] + work.wu_term[i]
            end
            
            mul!(work.weights_temp, work.bi_weights, work.weights)
            
            @simd for i in eachindex(work.weights)
                work.weights[i] += params.l_rate * work.weights_temp[i]
            end

            if maximum(abs, work.weights) > params.max_weight
                wts_blowup = true
                change = NaN
                break
            end
        end

        # Weight updates using pre-allocated arrays
        if !wts_blowup
            @. work.oldweights = work.weights - work.oldweights
            step += 1
            copyto!(work.delta, reshape(work.oldweights, 1, :))
            change = dot(work.delta, work.delta)
        end

        # Handle blowup with minimal allocations
        if wts_blowup || isnan(change) || isinf(change)
            step = 0
            change = NaN
            wts_blowup = false
            params.l_rate *= params.restart_factor
            copyto!(work.weights, work.startweights)
            copyto!(work.oldweights, work.startweights)
            continue
        end

        # Angle calculations with pre-allocated arrays
        if step > 2
            angledelta = acos(clamp(dot(work.delta, work.olddelta) / sqrt(change * oldchange), -1, 1))
            if params.degconst * angledelta > params.anneal_deg
                params.l_rate *= params.anneal_step
                copyto!(work.olddelta, work.delta)
                oldchange = change
            end
        elseif step == 1
            copyto!(work.olddelta, work.delta)
            oldchange = change
        end

        copyto!(work.oldweights, work.weights)

        # Convergence checks
        if step > 2 && change < params.w_change
            break
        elseif change > params.blowup
            params.l_rate *= params.blowup_fac
        end

        log_progress(step, change, params.l_rate, params.degconst * angledelta)
    end

    # Final optimized calculations
    if !isnothing(pca_components)
        work.weights = work.weights * sphere * pca_components'
    else
        work.weights = work.weights * sphere
    end

    sphere_final = Matrix{Float64}(I, size(work.weights, 2), size(work.weights, 2))
    unmixing = work.weights * sphere_final
    mixing = pinv(unmixing)

    # Optimized variance calculation
    meanvar = vec(sum(abs2, mixing, dims=1) .* sum(abs2, dat_ica, dims=2)' ./ (n_components * n_samples - 1))
    order = sortperm(meanvar, rev=true)

    return InfoIca(
        work.weights[order, :],
        sphere_final,
        mixing[:, order],
        unmixing[order, :],
        scale,
        ["IC$i" for i in 1:size(work.weights, 1)],
        data_labels
    )
end

"""
    tf_morlet(
        signal::Vector{Float64}, 
        fs::Real, 
        freqs::Vector{Float64}; 
        n_cycles::Union{Int,Vector{Int}}=7,
        fft_plan=nothing
    )

Compute time-frequency decomposition using Morlet wavelets with FFT plan optimization.

# Arguments
- `signal`: Input time series
- `fs`: Sampling frequency
- `freqs`: Frequencies to analyze
- `n_cycles`: Number of cycles for the Morlet wavelet
- `fft_plan`: Pre-computed FFT plan (optional)

# Returns
- Time-frequency power matrix
"""
function tf_morlet(
    signal::Vector{Float64}, 
    fs::Real, 
    freqs::Vector{Float64}; 
    n_cycles::Union{Int,Vector{Int}}=7,
    fft_plan=nothing
)
    n_times = length(signal)
    n_freqs = length(freqs)
    n_cycles = n_cycles isa Int ? fill(n_cycles, n_freqs) : n_cycles
    
    # Pre-allocate output
    tfr = zeros(Complex{Float64}, n_freqs, n_times)
    
    # Create or use FFT plan
    if isnothing(fft_plan)
        fft_plan = plan_fft(signal)
        signal_fft = fft_plan * signal
    else
        signal_fft = fft_plan * signal
    end
    
    # Frequency vector for FFT
    fft_freqs = FFTW.fftfreq(n_times, fs)
    
    # Pre-allocate wavelet
    wavelet = zeros(Complex{Float64}, n_times)
    
    # Compute wavelets and multiply with signal in frequency domain
    @inbounds for f_idx in 1:n_freqs
        # Create Morlet wavelet in frequency domain
        f = freqs[f_idx]
        σ = n_cycles[f_idx] / (2π * f)
        
        @simd for i in eachindex(fft_freqs)
            freq = fft_freqs[i]
            wavelet[i] = exp(-2im * π * freq * (-n_times/2/fs)) * 
                        exp(-((2π * freq - 2π * f)^2 * σ^2) / 2)
        end
        
        # Normalize
        wavelet ./= maximum(abs.(wavelet))
        
        # Multiply with signal in frequency domain and inverse FFT
        @views tfr[f_idx, :] .= ifft(signal_fft .* wavelet)
    end
    
    return abs2.(tfr)  # Return power
end

# Helper function to create and cache FFT plan
function create_fft_plan(signal::Vector{Float64})
    return plan_fft(signal; flags=FFTW.MEASURE)
end

"""
    generate_signal(
        freqs::Vector{Float64}, 
        amplitudes::Vector{Float64}, 
        duration::Float64, 
        fs::Float64; 
        noise_level::Float64=0.1
    ) -> Vector{Float64}

Generate a test signal composed of multiple sinusoids with optional noise.

# Arguments
- `freqs`: Vector of frequencies in Hz
- `amplitudes`: Vector of amplitudes for each frequency
- `duration`: Signal duration in seconds
- `fs`: Sampling frequency in Hz
- `noise_level`: Standard deviation of Gaussian noise to add (default: 0.1)

# Returns
- `Vector{Float64}`: Generated signal

# Example
```julia
# Generate 1s signal with 10Hz and 20Hz components
signal = generate_signal([10.0, 20.0], [1.0, 0.5], 1.0, 1000.0)
```
"""
function generate_signal(
    freqs::Vector{Float64}, 
    amplitudes::Vector{Float64}, 
    duration::Float64, 
    fs::Float64;
    noise_level::Float64=0.1
)
    @assert length(freqs) == length(amplitudes) "Number of frequencies must match number of amplitudes"
    
    # Pre-allocate time vector and signal
    n_samples = Int(round(duration * fs))
    t = range(0, duration, length=n_samples)
    signal = zeros(n_samples)
    
    # Add sinusoidal components
    @inbounds for (i, (freq, amp)) in enumerate(zip(freqs, amplitudes))
        @. signal += amp * sin(2π * freq * t)
    end
    
    # Add noise if requested
    if noise_level > 0
        signal .+= randn(n_samples) .* noise_level
    end
    
    return signal
end

# Convenience method for single frequency
function generate_signal(
    freq::Real, 
    amplitude::Real, 
    duration::Real, 
    fs::Real;
    noise_level::Real=0.1
)
    generate_signal([Float64(freq)], [Float64(amplitude)], 
                   Float64(duration), Float64(fs), 
                   noise_level=Float64(noise_level))
end

