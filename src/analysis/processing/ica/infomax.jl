# === ICA (Infomax) ===
# === INFOMAX ICA IMPLEMENTATION ===

"""
    _ica_preprocess(dat_ica::Matrix{Float64}, n_components::Int, params::IcaPrms)

Shared preprocessing helper for Infomax and Extended Infomax ICA:
- Validates dimensions and sample size
- Detects rank deficiency (e.g. from average referencing)
- Centers and scales data without mutating the caller's matrix
- Copies params to avoid mutating caller parameters
- Computes PCA reduction and analytic sphering via SVD
- Configures CPU or GPU work buffers and block size
"""
function _ica_preprocess(
    dat_ica::Matrix{Float64},
    n_components::Int,
    params::IcaPrms,
)
    n_channels, n_samples = size(dat_ica)

    # 1. Guards on n_components
    if n_components <= 0
        @minimal_error "n_components must be positive, got $n_components"
    end
    if n_components > n_channels
        @minimal_error "Requested $n_components components but only $n_channels channels available."
    end

    # 2. Sample size check
    min_samples = 20 * n_components^2
    if n_samples < min_samples
        @minimal_warning "Sample size ($n_samples) is small for $n_components components (< 20 * n_components^2 = $min_samples). ICA decomposition may be unreliable."
    end

    # 3. Non-mutating copy of inputs
    dat = copy(dat_ica)
    prms = deepcopy(params)

    # 4. Center and scale data
    original_mean = vec(mean(dat, dims = 2))
    dat .-= original_mean
    scale = sqrt(norm((dat * transpose(dat)) / n_samples))
    dat ./= scale

    # 5. SVD and rank check
    F = svd!(dat)
    s_max = F.S[1]
    s_target = F.S[n_components]
    if s_max > 0 && (s_target / s_max < 1e-7 || s_target < 1e-10)
        eff_rank = count(s -> (s / s_max) > 1e-7, F.S)
        @minimal_error "Data is rank-deficient for requested $n_components components (singular value ratio = $(round(s_target / s_max, sigdigits=3))). Estimated rank is $eff_rank. Please set n_components <= $eff_rank (e.g. n_channels - 1 for average-referenced data)."
    end

    pca_components = F.U[:, 1:n_components]
    eigenvalues = (F.S[1:n_components] .^ 2) ./ (n_samples - 1)
    sphere = diagm(1.0 ./ sqrt.(eigenvalues))

    # Analytic sphering: sqrt(n_samples - 1) * V^T matches sphere * pca_components' * dat
    dat_sphered = Matrix{Float64}(undef, n_components, n_samples)
    dat_sphered .= sqrt(n_samples - 1) .* F.Vt[1:n_components, :]

    # 6. GPU backend detection and block sizing
    gpu_active = false
    if prms.use_gpu
        if is_gpu_available()
            gpu_active = true
            @info "[GPU ACTIVATED] Running ICA on $(gpu_device_name())..."
            cpu_block = min(nextpow(2, max(Int(floor(sqrt(n_samples / 3.0))), 1)), 512)
            gpu_block = min(nextpow(2, cpu_block * 4), 1024)
            block = min(n_samples, max(gpu_block, 32))
            work = create_work_arrays(n_components, block)
        else
            @minimal_warning "Requested GPU acceleration (use_gpu=true), but no functional GPU package (CUDA.jl, AMDGPU.jl, Metal.jl) has been loaded. Please run 'using CUDA' or 'using AMDGPU' before calling run_ica. Falling back to CPU."
            raw_block = Int(floor(sqrt(n_samples / 3.0)))
            block = min(nextpow(2, max(raw_block, 1)), 512)
            work = create_work_arrays(n_components, block)
        end
    else
        raw_block = Int(floor(sqrt(n_samples / 3.0)))
        block = min(nextpow(2, max(raw_block, 1)), 512)
        work = create_work_arrays(n_components, block)
    end

    return (
        dat_sphered = dat_sphered,
        sphere = sphere,
        pca_components = pca_components,
        scale = scale,
        original_mean = original_mean,
        prms = prms,
        block = block,
        gpu_active = gpu_active,
        work = work,
        n_samples = n_samples,
    )
end

"""
    _ica_postprocess(work, sphere, pca_components, scale, original_mean, n_samples, n_components, filename, layout, is_sub_gaussian)

Shared postprocessing helper for Infomax and Extended Infomax ICA:
- Computes unmixing matrix from sensor space
- Computes mixing matrix via pseudoinverse
- Calculates component variance explained from activation power and mixing norms (EEGLAB runica method)
- Orders components by descending variance explained
- Constructs and returns the InfoIca result
"""
function _ica_postprocess(
    W_ica::Matrix{Float64},
    sphere::Matrix{Float64},
    pca_components::Matrix{Float64},
    scale::Float64,
    original_mean::Vector{Float64},
    n_samples::Int,
    n_components::Int,
    filename::String,
    layout::Layout,
    is_sub_gaussian::BitVector,
)
    # Exact activation power on whitened data (EEGLAB runica method):
    # Activation variance = (n_samples - 1) * sum(abs2, W_ica, dims=2)
    act_power = (n_samples - 1) .* vec(sum(abs2, W_ica; dims = 2))

    # Full unmixing matrix from sensor space
    unmixing = W_ica * sphere * transpose(pca_components)
    mixing = pinv(unmixing)

    # Total variance explained per component and order descending
    meanvar = vec(sum(abs2, mixing; dims = 1)) .* act_power ./ (n_components * n_samples - 1)
    meanvar_normalized = meanvar ./ sum(meanvar)
    order = sortperm(meanvar_normalized, rev = true)

    return InfoIca(
        filename,
        unmixing[order, :],
        mixing[:, order],
        sphere,
        meanvar_normalized[order],
        scale,
        original_mean,
        [Symbol("IC$i") for i = 1:size(unmixing, 1)],
        Dict{Int,Matrix{Float64}}(),
        layout,
        is_sub_gaussian[order],
    )
end

function _ica_postprocess(
    work::WorkArrays,
    sphere::Matrix{Float64},
    pca_components::Matrix{Float64},
    scale::Float64,
    original_mean::Vector{Float64},
    n_samples::Int,
    n_components::Int,
    filename::String,
    layout::Layout,
    is_sub_gaussian::BitVector,
)
    return _ica_postprocess(
        work.weights,
        sphere,
        pca_components,
        scale,
        original_mean,
        n_samples,
        n_components,
        filename,
        layout,
        is_sub_gaussian,
    )
end

"""Standard Infomax ICA implementation (super-Gaussian sources only)."""
function infomax_ica(
    dat_ica::Matrix{Float64},
    layout::Layout,
    filename::String;
    n_components::Int,
    params::IcaPrms = IcaPrms(),
    rng::AbstractRNG = Random.default_rng(),
)
    prep = _ica_preprocess(dat_ica, n_components, params)

    if prep.gpu_active
        _infomax_optimize_gpu!(prep.work, prep.dat_sphered, prep.prms, prep.block, rng)
    else
        _infomax_optimize_cpu!(prep.work, prep.dat_sphered, prep.prms, prep.block, rng)
    end

    return _ica_postprocess(
        prep.work,
        prep.sphere,
        prep.pca_components,
        prep.scale,
        prep.original_mean,
        prep.n_samples,
        n_components,
        filename,
        layout,
        falses(n_components),
    )
end

function _infomax_optimize_cpu!(
    work::WorkArrays,
    dat_ica::Matrix{Float64},
    params::IcaPrms,
    block::Int,
    rng::AbstractRNG = Random.default_rng(),
)
    n_channels, n_samples = size(dat_ica)
    n_components = n_channels

    step = 0
    wts_blowup = false
    change = 0.0
    oldchange = 0.0
    angledelta = 0.0
    converged = false

    permute_indices = Vector{Int}(undef, n_samples)

    @inbounds while step < params.max_iter
        randperm!(rng, permute_indices)

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
            # Using -tanh(0.5 * u) avoids @fastmath overflow on exp(-u)
            if block_size < block
                mul!(view(work.u, :, 1:block_size), work.weights, view(work.data_block, :, 1:block_size))
                len = n_components * block_size
                @fastmath @inbounds @simd for k = 1:len
                    work.y[k] = -tanh(0.5 * work.u[k])
                end
                mul!(work.wu_term, view(work.y, :, 1:block_size), transpose(view(work.u, :, 1:block_size)))
            else
                mul!(work.u, work.weights, work.data_block)
                @fastmath @inbounds @simd for k in eachindex(work.y)
                    work.y[k] = -tanh(0.5 * work.u[k])
                end
                mul!(work.wu_term, work.y, transpose(work.u))
            end

            mul!(work.weights_temp, work.wu_term, work.weights)
            @. work.weights += params.l_rate * (block_size * work.weights + work.weights_temp)

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
            if params.l_rate < 1e-6
                @minimal_error "Infomax failed to converge: learning rate fell below minimum floor (1e-6) due to repeated weight blow-ups. Check data for rank deficiency or artifacts."
            end
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
            if change < params.w_change
                converged = true
                break
            end
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

    if !converged
        @minimal_warning "Infomax reached maximum iterations ($(params.max_iter)) without converging (final change = $change, target = $(params.w_change))."
    end
end


# === ICA extended (Infomax extended) ===
# === EXTENDED INFOMAX ICA IMPLEMENTATION ===

"""
    infomax_extended_ica(dat_ica::Matrix{Float64}, layout::Layout, filename::String; n_components::Int, params::IcaPrms = IcaPrms(), rng::AbstractRNG = Random.default_rng())

Extended Infomax algorithm implementation that can separate both sub-Gaussian and super-Gaussian sources.

Extended Infomax is an enhancement of the standard Infomax algorithm that adapts to the statistical
properties of source signals. It uses kurtosis-based switching to handle both sub-Gaussian (kurtosis < 0)
and super-Gaussian (kurtosis > 0) sources, following the natural gradient rule of Lee, Girolami, & Sejnowski (1999).

# Arguments
- `dat_ica::Matrix{Float64}`: Data matrix (channels × samples), should be preprocessed
- `layout::Layout`: Layout information for channels
- `filename::String`: Filename or label for the dataset
- `n_components::Int`: Number of ICA components to extract
- `params::IcaPrms`: ICA parameters (uses max_iter, w_change from params)
- `rng::AbstractRNG`: Random number generator for reproducible permutations

# Returns
`InfoIca` with unmixing, mixing, sphere, variance, and metadata.

# References
Lee, T. W., Girolami, M., & Sejnowski, T. J. (1999). Independent component analysis using an extended 
infomax algorithm for mixed subgaussian and supergaussian sources. Neural computation, 11(2), 417-441.
"""
function infomax_extended_ica(
    dat_ica::Matrix{Float64},
    layout::Layout,
    filename::String;
    n_components::Int,
    params::IcaPrms = IcaPrms(),
    rng::AbstractRNG = Random.default_rng(),
)
    prep = _ica_preprocess(dat_ica, n_components, params)
    is_sub_gaussian = falses(n_components)

    if prep.gpu_active
        _infomax_extended_optimize_gpu!(prep.work, prep.dat_sphered, prep.prms, prep.block, is_sub_gaussian, rng)
    else
        _infomax_extended_optimize_cpu!(prep.work, prep.dat_sphered, prep.prms, prep.block, is_sub_gaussian, rng)
    end

    n_sub_final = count(is_sub_gaussian)
    n_super_final = n_components - n_sub_final
    @info "Extended Infomax completed: sup/sub-gauss ($n_super_final/$n_sub_final)"

    return _ica_postprocess(
        prep.work,
        prep.sphere,
        prep.pca_components,
        prep.scale,
        prep.original_mean,
        prep.n_samples,
        n_components,
        filename,
        layout,
        is_sub_gaussian,
    )
end

function _infomax_extended_optimize_cpu!(
    work::WorkArrays,
    dat_ica::Matrix{Float64},
    params::IcaPrms,
    block::Int,
    is_sub_gaussian::BitVector,
    rng::AbstractRNG = Random.default_rng(),
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
    converged = false

    permute_indices = Vector{Int}(undef, n_samples)

    # Pre-allocate kurtosis estimation buffers
    kurtsize = min(2000, n_samples)
    kurt_data = Matrix{Float64}(undef, n_channels, kurtsize)
    activations = Matrix{Float64}(undef, n_channels, kurtsize)

    n_super = n_channels
    n_sub = 0

    @inbounds while step < params.max_iter
        randperm!(rng, permute_indices)

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

            # forward pass: u = W * x
            if block_size < block
                u_block = view(work.u, :, 1:block_size)
                y_block = view(work.y, :, 1:block_size)
                mul!(u_block, work.weights, view(work.data_block, :, 1:block_size))

                # Lee et al. (1999) update rule:
                # y = tanh(u)
                len = n_channels * block_size
                @fastmath @inbounds @simd for k = 1:len
                    work.y[k] = tanh(work.u[k])
                end

                # wu_term = y * u'
                mul!(work.wu_term, y_block, transpose(u_block))

                # Apply K = -1 for sub-Gaussian sources, K = +1 for super-Gaussian sources
                for i = 1:n_channels
                    if is_sub_gaussian[i]
                        work.wu_term[i, :] .*= -1.0
                    end
                end

                # Add u * u' term: wu_term += u * u'
                mul!(work.wu_term, u_block, transpose(u_block), 1.0, 1.0)
            else
                mul!(work.u, work.weights, work.data_block)

                @fastmath @inbounds @simd for k in eachindex(work.y)
                    work.y[k] = tanh(work.u[k])
                end

                mul!(work.wu_term, work.y, transpose(work.u))

                for i = 1:n_channels
                    if is_sub_gaussian[i]
                        work.wu_term[i, :] .*= -1.0
                    end
                end

                mul!(work.wu_term, work.u, transpose(work.u), 1.0, 1.0)
            end

            # weights_temp = wu_term * weights
            mul!(work.weights_temp, work.wu_term, work.weights)

            # weights += l_rate * (block_size * weights - weights_temp)
            # Note the minus sign according to Lee et al. (1999): ΔW ∝ (B*I - wu_term) * W
            @. work.weights += params.l_rate * (block_size * work.weights - work.weights_temp)

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
            if params.l_rate < 1e-6
                @minimal_error "Extended Infomax failed to converge: learning rate fell below minimum floor (1e-6) due to repeated weight blow-ups. Check data for rank deficiency or artifacts."
            end
            work.weights .= work.startweights
            work.oldweights .= work.startweights
            continue
        end

        # Extended Infomax: update kurtosis signs periodically
        if step > 10 && step % 10 == 0
            # Sample using already-shuffled permute_indices (zero extra allocations)
            if kurtsize < n_samples
                @inbounds for j = 1:kurtsize
                    idx = permute_indices[j]
                    @simd for k = 1:n_channels
                        kurt_data[k, j] = dat_ica[k, idx]
                    end
                end
                mul!(activations, work.weights, kurt_data)
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

                    kurtosis_threshold = 0.05
                    old_is_sub = is_sub_gaussian[i]
                    is_sub_gaussian[i] = kurtosis < -kurtosis_threshold ? true : (kurtosis > kurtosis_threshold ? false : old_is_sub)
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
            if change < params.w_change
                converged = true
                break
            end
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

    if !converged
        @minimal_warning "Extended Infomax reached maximum iterations ($(params.max_iter)) without converging (final change = $change, target = $(params.w_change))."
    end
end


# === GPU versions of above (Infomax) ===
function _infomax_optimize_gpu!(
    work::WorkArrays,
    dat_ica::Matrix{Float64},
    params::IcaPrms,
    block::Int,
    rng::AbstractRNG = Random.default_rng(),
)
    n_channels, n_samples = size(dat_ica)
    n_components = n_channels

    step = 0
    wts_blowup = false
    change = 0.0
    oldchange = 0.0
    angledelta = 0.0
    converged = false

    permute_indices = Vector{Int}(undef, n_samples)
    perm_gpu = gpu_array(zeros(Int, n_samples))

    dat_ica_gpu = gpu_array(Float32.(dat_ica))
    dat_ica_gpu_shuffled = gpu_array(zeros(Float32, n_channels, n_samples))
    weights_gpu = gpu_array(Float32.(work.weights))
    weights_f32_cpu = Matrix{Float32}(undef, n_components, n_components)
    startweights_f32 = Float32.(work.startweights)

    u_gpu = gpu_array(zeros(Float32, n_components, block))
    y_gpu = gpu_array(zeros(Float32, n_components, block))
    wu_gpu = gpu_array(zeros(Float32, n_components, n_components))
    wtemp_gpu = gpu_array(zeros(Float32, n_components, n_components))

    # Pre-construct reusable views outside the loop
    u_view_full = view(u_gpu, :, 1:block)
    y_view_full = view(y_gpu, :, 1:block)
    u_view_full_t = transpose(u_view_full)
    block_f32 = Float32(block)

    last_block_size = rem(n_samples, block)
    has_remainder = last_block_size > 0
    u_view_last = has_remainder ? view(u_gpu, :, 1:last_block_size) : u_view_full
    y_view_last = has_remainder ? view(y_gpu, :, 1:last_block_size) : y_view_full
    u_view_last_t = transpose(u_view_last)
    last_block_f32 = Float32(last_block_size)

    # Pre-construct data block views into preallocated dat_ica_gpu_shuffled
    block_starts = 1:block:n_samples
    data_block_views = [view(dat_ica_gpu_shuffled, :, t:min(t + block - 1, n_samples)) for t in block_starts]
    n_blocks = length(block_starts)

    @inbounds while step < params.max_iter
        randperm!(rng, permute_indices)

        copyto!(perm_gpu, permute_indices)
        dat_ica_gpu_shuffled .= view(dat_ica_gpu, :, perm_gpu)

        l_rate_f32 = Float32(params.l_rate)

        for b_idx = 1:n_blocks
            data_block_gpu = data_block_views[b_idx]
            block_size = size(data_block_gpu, 2)

            if block_size == block
                mul!(u_view_full, weights_gpu, data_block_gpu)
                @. y_view_full = -tanh(0.5f0 * u_view_full)
                mul!(wu_gpu, y_view_full, u_view_full_t)
                mul!(wtemp_gpu, wu_gpu, weights_gpu)
                @. weights_gpu += l_rate_f32 * (block_f32 * weights_gpu + wtemp_gpu)
            else
                mul!(u_view_last, weights_gpu, data_block_gpu)
                @. y_view_last = -tanh(0.5f0 * u_view_last)
                mul!(wu_gpu, y_view_last, u_view_last_t)
                mul!(wtemp_gpu, wu_gpu, weights_gpu)
                @. weights_gpu += l_rate_f32 * (last_block_f32 * weights_gpu + wtemp_gpu)
            end
        end

        # Check boom on GPU once per step
        if maximum(abs, weights_gpu) > params.max_weight
            wts_blowup = true
            change = NaN
        end

        if !wts_blowup
            copyto!(weights_f32_cpu, weights_gpu)
            work.weights .= weights_f32_cpu
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
            if params.l_rate < 1e-6
                @minimal_error "Infomax failed to converge on GPU: learning rate fell below minimum floor (1e-6) due to repeated weight blow-ups."
            end
            work.weights .= work.startweights
            work.oldweights .= work.startweights
            copyto!(weights_gpu, startweights_f32)
            continue
        end

        if step > 2
            angledelta = acos(clamp(dot(work.delta, work.olddelta) / sqrt(change * oldchange), -1, 1))
            if params.degconst * angledelta > params.anneal_deg
                params.l_rate *= params.anneal_step
                work.olddelta .= work.delta
                oldchange = change
            end
            if change < params.w_change
                converged = true
                break
            end
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

    if !converged
        @minimal_warning "Infomax reached maximum iterations ($(params.max_iter)) without converging (final change = $change, target = $(params.w_change))."
    end

    copyto!(weights_f32_cpu, weights_gpu)
    work.weights .= weights_f32_cpu
end

function _infomax_extended_optimize_gpu!(
    work::WorkArrays,
    dat_ica::Matrix{Float64},
    params::IcaPrms,
    block::Int,
    is_sub_gaussian::BitVector,
    rng::AbstractRNG = Random.default_rng(),
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
    converged = false

    permute_indices = Vector{Int}(undef, n_samples)
    n_super = n_channels
    n_sub = 0

    perm_gpu = gpu_array(zeros(Int, n_samples))
    dat_ica_gpu = gpu_array(Float32.(dat_ica))
    dat_ica_gpu_shuffled = gpu_array(zeros(Float32, n_channels, n_samples))

    weights_gpu = gpu_array(Float32.(work.weights))
    weights_f32_cpu = Matrix{Float32}(undef, n_components, n_components)
    startweights_f32 = Float32.(work.startweights)

    u_gpu = gpu_array(zeros(Float32, n_components, block))
    y_gpu = gpu_array(zeros(Float32, n_components, block))
    wu_gpu = gpu_array(zeros(Float32, n_components, n_components))
    wtemp_gpu = gpu_array(zeros(Float32, n_components, n_components))

    # K signs vector: +1 for super-Gaussian, -1 for sub-Gaussian
    signs_cpu = Float32[is_sub_gaussian[i] ? -1.0f0 : 1.0f0 for i in 1:n_channels]
    signs_gpu = gpu_array(reshape(signs_cpu, n_channels, 1))

    kurtsize = min(2000, n_samples)
    activations_gpu = gpu_array(zeros(Float32, n_channels, kurtsize))
    act_sq = gpu_array(zeros(Float32, n_channels, kurtsize))
    u2_sum = gpu_array(zeros(Float32, n_channels, 1))
    u4_sum = gpu_array(zeros(Float32, n_channels, 1))
    u2_mean_sq_cpu = Vector{Float32}(undef, n_channels)
    u4_mean_cpu = Vector{Float32}(undef, n_channels)
    inv_kurtsize = 1.0f0 / Float32(kurtsize)

    # Pre-construct reusable views outside the loop
    u_view_full = view(u_gpu, :, 1:block)
    y_view_full = view(y_gpu, :, 1:block)
    u_view_full_t = transpose(u_view_full)
    block_f32 = Float32(block)

    last_block_size = rem(n_samples, block)
    has_remainder = last_block_size > 0
    u_view_last = has_remainder ? view(u_gpu, :, 1:last_block_size) : u_view_full
    y_view_last = has_remainder ? view(y_gpu, :, 1:last_block_size) : y_view_full
    u_view_last_t = transpose(u_view_last)
    last_block_f32 = Float32(last_block_size)

    # Pre-construct data block views into preallocated dat_ica_gpu_shuffled
    block_starts = 1:block:n_samples
    data_block_views = [view(dat_ica_gpu_shuffled, :, t:min(t + block - 1, n_samples)) for t in block_starts]
    n_blocks = length(block_starts)

    @inbounds while step < params.max_iter
        randperm!(rng, permute_indices)

        copyto!(perm_gpu, permute_indices)
        dat_ica_gpu_shuffled .= view(dat_ica_gpu, :, perm_gpu)

        l_rate_f32 = Float32(params.l_rate)

        for b_idx = 1:n_blocks
            data_block_gpu = data_block_views[b_idx]
            block_size = size(data_block_gpu, 2)

            if block_size == block
                mul!(u_view_full, weights_gpu, data_block_gpu)
                @. y_view_full = tanh(u_view_full)
                mul!(wu_gpu, y_view_full, u_view_full_t)
                @. wu_gpu *= signs_gpu
                mul!(wu_gpu, u_view_full, u_view_full_t, 1.0f0, 1.0f0)
                mul!(wtemp_gpu, wu_gpu, weights_gpu)
                @. weights_gpu += l_rate_f32 * (block_f32 * weights_gpu - wtemp_gpu)
            else
                mul!(u_view_last, weights_gpu, data_block_gpu)
                @. y_view_last = tanh(u_view_last)
                mul!(wu_gpu, y_view_last, u_view_last_t)
                @. wu_gpu *= signs_gpu
                mul!(wu_gpu, u_view_last, u_view_last_t, 1.0f0, 1.0f0)
                mul!(wtemp_gpu, wu_gpu, weights_gpu)
                @. weights_gpu += l_rate_f32 * (last_block_f32 * weights_gpu - wtemp_gpu)
            end
        end

        # Check boom on GPU once per step
        if maximum(abs, weights_gpu) > params.max_weight
            wts_blowup = true
            change = NaN
        end

        if !wts_blowup
            copyto!(weights_f32_cpu, weights_gpu)
            work.weights .= weights_f32_cpu
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
            if params.l_rate < 1e-6
                @minimal_error "Extended Infomax failed to converge on GPU: learning rate fell below minimum floor (1e-6) due to repeated weight blow-ups."
            end
            work.weights .= work.startweights
            work.oldweights .= work.startweights
            copyto!(weights_gpu, startweights_f32)
            continue
        end

        # Extended Infomax: update kurtosis signs periodically
        if step > 10 && step % 10 == 0
            # dat_ica_gpu_shuffled is already randomized by perm_gpu
            kurt_view = view(dat_ica_gpu_shuffled, :, 1:kurtsize)
            mul!(activations_gpu, weights_gpu, kurt_view)

            # Compute raw kurtosis entirely on the GPU
            @. act_sq = activations_gpu^2
            sum!(u2_sum, act_sq)
            @. act_sq = act_sq^2
            sum!(u4_sum, act_sq)

            # Fetch back to CPU for momentum smoothing
            copyto!(u2_mean_sq_cpu, u2_sum)
            copyto!(u4_mean_cpu, u4_sum)

            @. u2_mean_sq_cpu = (u2_mean_sq_cpu * inv_kurtsize)^2
            @. u4_mean_cpu = u4_mean_cpu * inv_kurtsize

            signs_changed = false

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
                        signs_changed = true
                    end
                end
            end

            # Push updated signs vector to GPU if any changed
            if signs_changed
                signs_cpu .= Float32[is_sub_gaussian[i] ? -1.0f0 : 1.0f0 for i in 1:n_channels]
                copyto!(signs_gpu, reshape(signs_cpu, n_channels, 1))
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
            if change < params.w_change
                converged = true
                break
            end
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

    if !converged
        @minimal_warning "Extended Infomax reached maximum iterations ($(params.max_iter)) without converging (final change = $change, target = $(params.w_change))."
    end

    copyto!(weights_f32_cpu, weights_gpu)
    work.weights .= weights_f32_cpu
end
