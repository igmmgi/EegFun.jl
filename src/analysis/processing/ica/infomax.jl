

# === ICA (Infomax) ===
# === INFOMAX ICA IMPLEMENTATION ===

"""Standard Infomax ICA implementation (super-Gaussian sources only)."""
function infomax_ica(dat_ica::Matrix{Float64}, layout::Layout, filename::String; n_components::Int, params::IcaPrms = IcaPrms())

    # Store original mean before removing it
    original_mean = vec(mean(dat_ica, dims = 2))

    n_channels, n_samples = size(dat_ica)

    # Center and scale data
    dat_ica .-= original_mean
    scale = sqrt(norm((dat_ica * transpose(dat_ica)) / n_samples))
    dat_ica ./= scale

    # PCA reduction using high-precision SVD
    # This prevents smallest components from degrading into numerical noise
    F = svd(dat_ica)
    pca_components = F.U[:, 1:n_components]

    # Analytic sphering using exact singular values (skips a 260 MiB workspace allocation)
    eigenvalues = (F.S[1:n_components] .^ 2) ./ (n_samples - 1)
    sphere = diagm(1.0 ./ sqrt.(eigenvalues))

    # Apply transform and create final workspace matrix (n_components × n_samples)
    transform_matrix = sphere * transpose(pca_components)
    dat_ica_sphered = Matrix{Float64}(undef, n_components, n_samples)
    mul!(dat_ica_sphered, transform_matrix, dat_ica)
    dat_ica = dat_ica_sphered

    # initialize
    n_channels = size(dat_ica, 1)
    n_samples = size(dat_ica, 2)

    # Detect GPU hardware backend
    gpu_active = false
    if params.use_gpu
        if is_gpu_available()
            gpu_active = true
            @info "[GPU ACTIVATED] Running Infomax ICA on $(gpu_device_name())..."

            # Dynamically calculate GPU block size
            cpu_block = min(nextpow(2, max(Int(floor(sqrt(n_samples / 3.0))), 1)), 512)
            gpu_block = min(nextpow(2, cpu_block * 4), 1024)
            block = min(n_samples, max(gpu_block, 32))

            work = create_work_arrays(n_channels, block)
        else
            @minimal_warning "Requested GPU acceleration (use_gpu=true), but no functional GPU package (CUDA.jl, AMDGPU.jl, Metal.jl) has been loaded. Please run 'using CUDA' or 'using AMDGPU' before calling run_ica. Falling back to CPU."
            raw_block = Int(floor(sqrt(n_samples / 3.0)))
            block = min(nextpow(2, max(raw_block, 1)), 512)
            work = create_work_arrays(n_channels, block)
        end
    else
        raw_block = Int(floor(sqrt(n_samples / 3.0)))
        block = min(nextpow(2, max(raw_block, 1)), 512)
        work = create_work_arrays(n_channels, block)
    end

    if gpu_active
        _infomax_optimize_gpu!(work, dat_ica, params, block)
    else
        _infomax_optimize_cpu!(work, dat_ica, params, block)
    end

    # Final calculations
    work.weights = work.weights * sphere * pca_components'
    mixing = pinv(work.weights)

    # calculate total variance explained and order
    meanvar = vec(sum(abs2, mixing, dims = 1) .* sum(abs2, dat_ica, dims = 2)' ./ (n_components * n_samples - 1))
    meanvar_normalized = meanvar ./ sum(meanvar)
    order = sortperm(meanvar_normalized, rev = true)

    return InfoIca(
        filename,
        work.weights[order, :],
        mixing[:, order],
        sphere,
        meanvar_normalized[order],
        scale,
        original_mean,
        [Symbol("IC$i") for i = 1:size(work.weights, 1)],
        Dict{Int,Matrix{Float64}}(),
        layout,
        falses(size(work.weights, 1)),  # Regular Infomax: all super-Gaussian (all false)
    )

end

function _infomax_optimize_cpu!(work::WorkArrays, dat_ica::Matrix{Float64}, params::IcaPrms, block::Int)
    n_channels, n_samples = size(dat_ica)
    n_components = n_channels

    step = 0
    wts_blowup = false
    change = 0.0
    oldchange = 0.0
    angledelta = 0.0

    permute_indices = Vector{Int}(undef, n_samples)

    @inbounds while step < params.max_iter
        randperm!(permute_indices)

        for t = 1:block:n_samples
            block_end = min(t + block - 1, n_samples)
            block_size = block_end - t + 1

            # extract data block
            @inbounds for j = 1:block_size
                idx = permute_indices[t+j-1]
                @simd for i = 1:n_channels
                    work.data_block[i, j] = dat_ica[i, idx]
                end
            end

            # forward pass & weight update
            if block_size < block
                mul!(view(work.u, :, 1:block_size), work.weights, view(work.data_block, :, 1:block_size))
            else
                mul!(work.u, work.weights, work.data_block)
            end

            @fastmath @inbounds for i = 1:n_components
                @simd for j = 1:block_size
                    work.y[i, j] = 1.0 - 2.0 / (1.0 + exp(-work.u[i, j]))
                end
            end

            # update weights 
            if block_size < block
                mul!(work.wu_term, view(work.y, :, 1:block_size), transpose(view(work.u, :, 1:block_size)))
                mul!(work.weights_temp, work.wu_term, work.weights)
                @. work.weights += params.l_rate * (block_size * work.weights + work.weights_temp)
            else
                mul!(work.wu_term, work.y, transpose(work.u))
                mul!(work.weights_temp, work.wu_term, work.weights)
                @. work.weights += params.l_rate * (block_size * work.weights + work.weights_temp)
            end

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
            work.delta .= work.oldweights
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
            change < params.w_change && break
        elseif step == 1
            work.olddelta .= work.delta
            oldchange = change
        end

        work.oldweights .= work.weights
        change > params.blowup && (params.l_rate *= params.blowup_fac)

        if step == 1 || step % 10 == 0
            @info Printf.@sprintf(
                "Infomax step %d, change = %.7f, lrate = %.7f, angle = %.1f",
                step,
                change,
                params.l_rate,
                (params.degconst) * angledelta
            )
        end
    end
end


# === ICA extended (Infomax extended) ===
# === EXTENDED INFOMAX ICA IMPLEMENTATION ===

"""
    infomax_extended_ica(dat_ica::Matrix{Float64}, layout::Layout; n_components::Int, params::IcaPrms = IcaPrms())

Extended Infomax algorithm implementation that can separate both sub-Gaussian and super-Gaussian sources.

Extended Infomax is an enhancement of the standard Infomax algorithm that adapts to the statistical
properties of source signals. It uses kurtosis-based switching to handle both sub-Gaussian (kurtosis < 0)
and super-Gaussian (kurtosis > 0) sources, making it more versatile than standard Infomax.

# Arguments
- `dat_ica::Matrix{Float64}`: Data matrix (channels × samples), should be preprocessed
- `layout::Layout`: Layout information for channels
- `n_components::Int`: Number of ICA components to extract
- `params::IcaPrms`: ICA parameters (uses max_iter, w_change from params)

# Returns
`InfoIca` with unmixing, mixing, sphere, variance, and metadata.

# References
Lee, T. W., Girolami, M., & Sejnowski, T. J. (1999). Independent component analysis using an extended 
infomax algorithm for mixed subgaussian and supergaussian sources. Neural computation, 11(2), 417-441.
"""
function infomax_extended_ica(dat_ica::Matrix{Float64}, layout::Layout, filename::String; n_components::Int, params::IcaPrms = IcaPrms())

    # Store original mean before removing it
    original_mean = vec(mean(dat_ica, dims = 2))

    n_channels, n_samples = size(dat_ica)

    # Center and scale data
    dat_ica .-= original_mean
    scale = sqrt(norm((dat_ica * transpose(dat_ica)) / n_samples))
    dat_ica ./= scale

    # PCA reduction using high-precision SVD
    # This prevents smallest components from degrading into numerical noise
    F = svd(dat_ica)
    pca_components = F.U[:, 1:n_components]

    # Analytic sphering using exact singular values (skips a 260 MiB workspace allocation)
    eigenvalues = (F.S[1:n_components] .^ 2) ./ (n_samples - 1)
    sphere = diagm(1.0 ./ sqrt.(eigenvalues))

    # Apply transform and create final workspace matrix (n_components × n_samples)
    transform_matrix = sphere * transpose(pca_components)
    dat_ica_sphered = Matrix{Float64}(undef, n_components, n_samples)
    mul!(dat_ica_sphered, transform_matrix, dat_ica)
    dat_ica = dat_ica_sphered

    # initialize
    n_channels = size(dat_ica, 1)
    n_samples = size(dat_ica, 2)

    # Track kurtosis signs for each component (true = sub-Gaussian, false = super-Gaussian)
    is_sub_gaussian = falses(n_channels)

    # Detect GPU hardware backend
    gpu_active = false
    if params.use_gpu
        if is_gpu_available()
            gpu_active = true
            @info "[GPU ACTIVATED] Running Extended Infomax ICA on $(gpu_device_name())..."

            # Dynamically calculate GPU block size
            cpu_block = min(nextpow(2, max(Int(floor(sqrt(n_samples / 3.0))), 1)), 512)
            gpu_block = min(nextpow(2, cpu_block * 4), 1024)
            block = min(n_samples, max(gpu_block, 32))

            work = create_work_arrays(n_channels, block)
        else
            @minimal_warning "Requested GPU acceleration (use_gpu=true), but no functional GPU package (CUDA.jl, AMDGPU.jl, Metal.jl) has been loaded. Please run 'using CUDA' or 'using AMDGPU' before calling run_ica. Falling back to CPU."
            raw_block = Int(floor(sqrt(n_samples / 3.0)))
            block = min(nextpow(2, max(raw_block, 1)), 512)
            work = create_work_arrays(n_channels, block)
        end
    else
        raw_block = Int(floor(sqrt(n_samples / 3.0)))
        block = min(nextpow(2, max(raw_block, 1)), 512)
        work = create_work_arrays(n_channels, block)
    end

    if gpu_active
        _infomax_extended_optimize_gpu!(work, dat_ica, params, block, is_sub_gaussian)
    else
        _infomax_extended_optimize_cpu!(work, dat_ica, params, block, is_sub_gaussian)
    end

    # Final calculations
    work.weights = work.weights * sphere * pca_components'
    mixing = pinv(work.weights)

    # Log final kurtosis distribution
    n_sub_final = count(is_sub_gaussian)
    n_super_final = n_channels - n_sub_final
    @info "Extended Infomax completed: sup/sub-gauss ($n_super_final/$n_sub_final)"

    # calculate total variance explained and order
    meanvar = vec(sum(abs2, mixing, dims = 1) .* sum(abs2, dat_ica, dims = 2)' ./ (n_components * n_samples - 1))
    meanvar_normalized = meanvar ./ sum(meanvar)
    order = sortperm(meanvar_normalized, rev = true)

    return InfoIca(
        filename,
        work.weights[order, :],
        mixing[:, order],
        sphere,
        meanvar_normalized[order],
        scale,
        original_mean,
        [Symbol("IC$i") for i = 1:size(work.weights, 1)],
        Dict{Int,Matrix{Float64}}(),
        layout,
        is_sub_gaussian[order],
    )
end

function _infomax_extended_optimize_cpu!(
    work::WorkArrays,
    dat_ica::Matrix{Float64},
    params::IcaPrms,
    block::Int,
    is_sub_gaussian::BitVector,
)
    n_channels, n_samples = size(dat_ica)
    n_components = n_channels

    old_kurtosis = zeros(Float64, n_channels)
    extmomentum = 0.5

    step = 0
    wts_blowup = false
    change = 0.0
    oldchange = 0.0
    angledelta = 0.0

    permute_indices = Vector{Int}(undef, n_samples)

    n_super = n_channels
    n_sub = 0

    @inbounds while step < params.max_iter
        randperm!(permute_indices)

        for t = 1:block:n_samples
            block_end = min(t + block - 1, n_samples)
            block_size = block_end - t + 1

            # extract data block
            @inbounds for j = 1:block_size
                idx = permute_indices[t+j-1]
                @simd for i = 1:n_channels
                    work.data_block[i, j] = dat_ica[i, idx]
                end
            end

            # forward pass
            if block_size < block
                mul!(view(work.u, :, 1:block_size), work.weights, view(work.data_block, :, 1:block_size))
            else
                mul!(work.u, work.weights, work.data_block)
            end

            # Extended Infomax: use different nonlinearities based on kurtosis sign
            @fastmath @inbounds for i = 1:n_channels
                if !is_sub_gaussian[i] # Super-Gaussian
                    @simd for j = 1:block_size
                        work.y[i, j] = 1.0 - 2.0 / (1.0 + exp(-work.u[i, j]))
                    end
                else # Sub-Gaussian
                    @simd for j = 1:block_size
                        work.y[i, j] = -tanh(work.u[i, j])
                    end
                end
            end

            # update weights 
            if block_size < block
                mul!(work.wu_term, view(work.y, :, 1:block_size), transpose(view(work.u, :, 1:block_size)))
                mul!(work.weights_temp, work.wu_term, work.weights)
                @. work.weights += params.l_rate * (block_size * work.weights + work.weights_temp)
            else
                mul!(work.wu_term, work.y, transpose(work.u))
                mul!(work.weights_temp, work.wu_term, work.weights)
                @. work.weights += params.l_rate * (block_size * work.weights + work.weights_temp)
            end

            if maximum(abs, work.weights) > params.max_weight
                wts_blowup = true
                change = NaN
                break
            end
        end

        if !wts_blowup
            work.oldweights .-= work.weights
            step += 1
            work.delta .= work.oldweights
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

        # Extended Infomax: update kurtosis signs periodically
        if step > 10 && step % 10 == 0
            kurtsize = min(2000, n_samples)
            n_switched = 0
            kurtosis_values = Vector{Float64}(undef, n_channels)

            # CPU Kurtosis logic
            activations = Matrix{Float64}(undef, n_channels, kurtsize)
            if kurtsize < n_samples
                rp = randperm(n_samples)[1:kurtsize]
                @inbounds for j = 1:kurtsize
                    idx = rp[j]
                    for i = 1:n_channels
                        sum_val = 0.0
                        @simd for k = 1:n_channels
                            sum_val += work.weights[i, k] * dat_ica[k, idx]
                        end
                        activations[i, j] = sum_val
                    end
                end
            else
                mul!(activations, work.weights, dat_ica)
            end

            for i = 1:n_channels
                u2_sum = 0.0
                u4_sum = 0.0
                @inbounds @simd for j = 1:kurtsize
                    val = activations[i, j]
                    val2 = val * val
                    u2_sum += val2
                    u4_sum += val2 * val2
                end
                u2_mean_sq = (u2_sum / kurtsize)^2
                if u2_mean_sq > eps(Float64)
                    u4_mean = u4_sum / kurtsize
                    kurtosis_raw = (u4_mean / u2_mean_sq) - 3.0

                    if old_kurtosis[i] != 0.0
                        kurtosis = extmomentum * old_kurtosis[i] + (1.0 - extmomentum) * kurtosis_raw
                    else
                        kurtosis = kurtosis_raw
                    end
                    old_kurtosis[i] = kurtosis
                    kurtosis_values[i] = kurtosis

                    kurtosis_threshold = 0.05
                    old_is_sub = is_sub_gaussian[i]
                    is_sub_gaussian[i] = kurtosis < -kurtosis_threshold ? true : (kurtosis > kurtosis_threshold ? false : old_is_sub)

                    if old_is_sub != is_sub_gaussian[i]
                        n_switched += 1
                    end
                else
                    kurtosis_values[i] = 0.0
                end
            end

            # Update kurtosis counts
            n_sub = count(is_sub_gaussian)
            n_super = n_channels - n_sub
        end

        if step > 2
            angledelta = acos(clamp(dot(work.delta, work.olddelta) / sqrt(change * oldchange), -1, 1))
            if params.degconst * angledelta > params.anneal_deg
                params.l_rate *= params.anneal_step
                work.olddelta .= work.delta
                oldchange = change
            end
            change < params.w_change && break
        elseif step == 1
            work.olddelta .= work.delta
            oldchange = change
        end

        work.oldweights .= work.weights
        change > params.blowup && (params.l_rate *= params.blowup_fac)

        if step == 1 || step % 10 == 0
            @info Printf.@sprintf(
                "Extended-Infomax step %d, change = %.7f, lrate = %.7f, angle = %.1f, sup/sub-gauss: %d/%d",
                step,
                change,
                params.l_rate,
                (params.degconst) * angledelta,
                n_super,
                n_sub
            )
        end
    end
end


# === GPU versions of above (Infomax) ===
function _infomax_optimize_gpu!(work::WorkArrays, dat_ica::Matrix{Float64}, params::IcaPrms, block::Int)
    n_channels, n_samples = size(dat_ica)
    n_components = n_channels

    step = 0
    wts_blowup = false
    change = 0.0
    oldchange = 0.0
    angledelta = 0.0

    permute_indices = Vector{Int}(undef, n_samples)

    dat_ica_gpu = gpu_array(Float32.(dat_ica))
    weights_gpu = gpu_array(Float32.(work.weights))
    u_gpu = gpu_array(zeros(Float32, n_components, block))
    y_gpu = gpu_array(zeros(Float32, n_components, block))
    wu_gpu = gpu_array(zeros(Float32, n_components, n_components))
    wtemp_gpu = gpu_array(zeros(Float32, n_components, n_components))

    @inbounds while step < params.max_iter
        randperm!(permute_indices)

        perm_gpu = gpu_array(permute_indices)
        dat_ica_gpu_shuffled = dat_ica_gpu[:, perm_gpu]

        u_view_full = view(u_gpu, :, 1:block)
        y_view_full = view(y_gpu, :, 1:block)
        u_view_full_t = transpose(u_view_full)

        for t = 1:block:n_samples
            block_end = min(t + block - 1, n_samples)
            block_size = block_end - t + 1

            try
                data_block_gpu = view(dat_ica_gpu_shuffled, :, t:block_end)

                if block_size == block
                    mul!(u_view_full, weights_gpu, data_block_gpu)
                    @. y_view_full = 1.0f0 - 2.0f0 / (1.0f0 + exp(-u_view_full))
                    mul!(wu_gpu, y_view_full, u_view_full_t)
                else
                    u_view = view(u_gpu, :, 1:block_size)
                    y_view = view(y_gpu, :, 1:block_size)
                    mul!(u_view, weights_gpu, data_block_gpu)
                    @. y_view = 1.0f0 - 2.0f0 / (1.0f0 + exp(-u_view))
                    mul!(wu_gpu, y_view, transpose(u_view))
                end

                mul!(wtemp_gpu, wu_gpu, weights_gpu)
                @. weights_gpu += Float32(params.l_rate) * (Float32(block_size) * weights_gpu + wtemp_gpu)
            catch e
                @error "GPU processing failed in inner loop, throwing error." exception = (e, catch_backtrace())
                rethrow(e)
            end
        end

        # Check boom on GPU once per step
        if maximum(abs, weights_gpu) > params.max_weight
            wts_blowup = true
            change = NaN
        end

        if !wts_blowup
            work.weights .= Array(weights_gpu)
            work.oldweights .-= work.weights
            step += 1
            work.delta .= work.oldweights
            change = dot(work.delta, work.delta)
        end

        if wts_blowup || isnan(change) || isinf(change)
            step = 0
            change = NaN
            wts_blowup = false
            params.l_rate *= params.restart_factor
            work.weights .= work.startweights
            work.oldweights .= work.startweights
            weights_gpu .= gpu_array(Float32.(work.weights))
            continue
        end

        if step > 2
            angledelta = acos(clamp(dot(work.delta, work.olddelta) / sqrt(change * oldchange), -1, 1))
            if params.degconst * angledelta > params.anneal_deg
                params.l_rate *= params.anneal_step
                work.olddelta .= work.delta
                oldchange = change
            end
            change < params.w_change && break
        elseif step == 1
            work.olddelta .= work.delta
            oldchange = change
        end

        work.oldweights .= work.weights
        change > params.blowup && (params.l_rate *= params.blowup_fac)

        if step == 1 || step % 10 == 0
            @info Printf.@sprintf(
                "Infomax step %d, change = %.7f, lrate = %.7f, angle = %.1f",
                step,
                change,
                params.l_rate,
                (params.degconst) * angledelta
            )
        end
    end

    work.weights .= Array(weights_gpu)
end

function _infomax_extended_optimize_gpu!(
    work::WorkArrays,
    dat_ica::Matrix{Float64},
    params::IcaPrms,
    block::Int,
    is_sub_gaussian::BitVector,
)
    n_channels, n_samples = size(dat_ica)
    n_components = n_channels

    old_kurtosis = zeros(Float64, n_channels)
    extmomentum = 0.5

    step = 0
    wts_blowup = false
    change = 0.0
    oldchange = 0.0
    angledelta = 0.0

    permute_indices = Vector{Int}(undef, n_samples)
    n_super = n_channels
    n_sub = 0

    dat_ica_gpu = gpu_array(Float32.(dat_ica))
    weights_gpu = gpu_array(Float32.(work.weights))
    u_gpu = gpu_array(zeros(Float32, n_components, block))
    y_gpu = gpu_array(zeros(Float32, n_components, block))
    wu_gpu = gpu_array(zeros(Float32, n_components, n_components))
    wtemp_gpu = gpu_array(zeros(Float32, n_components, n_components))
    is_sub_gaussian_gpu = gpu_array(is_sub_gaussian)

    kurtsize = min(2000, n_samples)
    activations_gpu = gpu_array(zeros(Float32, n_channels, kurtsize))

    @inbounds while step < params.max_iter
        randperm!(permute_indices)

        perm_gpu = gpu_array(permute_indices)
        dat_ica_gpu_shuffled = dat_ica_gpu[:, perm_gpu]

        u_view_full = view(u_gpu, :, 1:block)
        y_view_full = view(y_gpu, :, 1:block)
        u_view_full_t = transpose(u_view_full)

        for t = 1:block:n_samples
            block_end = min(t + block - 1, n_samples)
            block_size = block_end - t + 1

            try
                data_block_gpu = view(dat_ica_gpu_shuffled, :, t:block_end)

                if block_size == block
                    mul!(u_view_full, weights_gpu, data_block_gpu)
                    @. y_view_full = ifelse(is_sub_gaussian_gpu, -tanh(u_view_full), 1.0f0 - 2.0f0 / (1.0f0 + exp(-u_view_full)))
                    mul!(wu_gpu, y_view_full, u_view_full_t)
                else
                    u_view = view(u_gpu, :, 1:block_size)
                    y_view = view(y_gpu, :, 1:block_size)
                    mul!(u_view, weights_gpu, data_block_gpu)
                    @. y_view = ifelse(is_sub_gaussian_gpu, -tanh(u_view), 1.0f0 - 2.0f0 / (1.0f0 + exp(-u_view)))
                    mul!(wu_gpu, y_view, transpose(u_view))
                end

                mul!(wtemp_gpu, wu_gpu, weights_gpu)
                @. weights_gpu += Float32(params.l_rate) * (Float32(block_size) * weights_gpu + wtemp_gpu)
            catch e
                @error "GPU processing failed in inner loop, throwing error." exception = (e, catch_backtrace())
                rethrow(e)
            end
        end

        # Check boom on GPU once per step
        if maximum(abs, weights_gpu) > params.max_weight
            wts_blowup = true
            change = NaN
        end

        if !wts_blowup
            work.weights .= Array(weights_gpu)
            work.oldweights .-= work.weights
            step += 1
            work.delta .= work.oldweights
            change = dot(work.delta, work.delta)
        end

        if wts_blowup || isnan(change) || isinf(change)
            step = 0
            change = NaN
            wts_blowup = false
            params.l_rate *= params.restart_factor
            work.weights .= work.startweights
            work.oldweights .= work.startweights
            weights_gpu .= gpu_array(Float32.(work.weights))
            continue
        end

        # Extended Infomax: update kurtosis signs periodically
        if step > 10 && step % 10 == 0
            rp_cpu = randperm(n_samples)[1:kurtsize]
            rp_gpu = gpu_array(rp_cpu)
            dat_subset = dat_ica_gpu[:, rp_gpu]
            mul!(activations_gpu, weights_gpu, dat_subset)

            # Compute raw kurtosis entirely on the GPU using mapreduce patterns or broadcasting
            act_sq = activations_gpu .^ 2
            u2_sum = sum(act_sq, dims = 2)
            u2_mean_sq = (u2_sum ./ kurtsize) .^ 2

            act_quad = act_sq .^ 2
            u4_mean = sum(act_quad, dims = 2) ./ kurtsize

            # Fetch back to CPU for momentum smoothing
            u2_mean_sq_cpu = Array(u2_mean_sq)
            u4_mean_cpu = Array(u4_mean)

            n_switched = 0

            for i = 1:n_channels
                if u2_mean_sq_cpu[i] > eps(Float64)
                    kurtosis_raw = (u4_mean_cpu[i] / u2_mean_sq_cpu[i]) - 3.0

                    if old_kurtosis[i] != 0.0
                        kurtosis = extmomentum * old_kurtosis[i] + (1.0 - extmomentum) * kurtosis_raw
                    else
                        kurtosis = kurtosis_raw
                    end
                    old_kurtosis[i] = kurtosis

                    kurtosis_threshold = 0.05
                    old_is_sub = is_sub_gaussian[i]
                    is_sub_gaussian[i] = kurtosis < -kurtosis_threshold ? true : (kurtosis > kurtosis_threshold ? false : old_is_sub)

                    if old_is_sub != is_sub_gaussian[i]
                        n_switched += 1
                    end
                end
            end

            # Push updated boolean vector to GPU
            is_sub_gaussian_gpu .= gpu_array(is_sub_gaussian)

            # Update kurtosis counts
            n_sub = count(is_sub_gaussian)
            n_super = n_channels - n_sub
        end

        if step > 2
            angledelta = acos(clamp(dot(work.delta, work.olddelta) / sqrt(change * oldchange), -1, 1))
            if params.degconst * angledelta > params.anneal_deg
                params.l_rate *= params.anneal_step
                work.olddelta .= work.delta
                oldchange = change
            end
            change < params.w_change && break
        elseif step == 1
            work.olddelta .= work.delta
            oldchange = change
        end

        work.oldweights .= work.weights
        change > params.blowup && (params.l_rate *= params.blowup_fac)

        if step == 1 || step % 10 == 0
            @info Printf.@sprintf(
                "Extended-Infomax step %d, change = %.7f, lrate = %.7f, angle = %.1f, sup/sub-gauss: %d/%d",
                step,
                change,
                params.l_rate,
                (params.degconst) * angledelta,
                n_super,
                n_sub
            )
        end
    end

    work.weights .= Array(weights_gpu)
end


