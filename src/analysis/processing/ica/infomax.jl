function _ica_preprocess(dat_ica::Matrix{Float64}, n_components::Int, params::IcaPrms)
    n_channels, n_samples = size(dat_ica)

    n_components <= 0 && @minimal_error("n_components must be positive, got $n_components")
    n_components > n_channels && @minimal_error("Requested $n_components components but only $n_channels channels available.")

    min_samples = 20 * n_components^2
    if n_samples < min_samples
        @minimal_warning "Sample size ($n_samples) is small for $n_components components (< 20 * n_components^2 = $min_samples). ICA decomposition may be unreliable."
    end

    prms = deepcopy(params)
    original_mean = vec(mean(dat_ica, dims = 2))
    dat = dat_ica .- original_mean
    cov_mat = Matrix{Float64}(undef, n_channels, n_channels)
    mul!(cov_mat, dat, transpose(dat), 1.0 / n_samples, 0.0)
    scale = sqrt(norm(cov_mat))

    local pca_components::Matrix{Float64}
    local sphere::Matrix{Float64}
    local dat_sphered::Matrix{Float64}

    used_eigen = false
    if n_samples > 2 * n_channels && scale > 0.0
        cov_scaled = (n_samples / (scale * scale)) .* cov_mat
        ef = eigen!(Symmetric(cov_scaled))
        l_max = ef.values[end]
        l_target = ef.values[end-n_components+1]
        if l_max > 0.0 && l_target > 1e-8 * l_max && l_target > 1e-20
            idx = n_channels:-1:(n_channels-n_components+1)
            pca_components = ef.vectors[:, idx]
            @inbounds for j = 1:n_components
                col = @view pca_components[:, j]
                argmax(abs, col) < 0.0 && (col .*= -1.0)
            end
            eigenvalues = ef.values[idx] ./ (n_samples - 1)
            inv_sqrt_ev = 1.0 ./ sqrt.(eigenvalues)
            sphere = diagm(inv_sqrt_ev)
            whiten_mat = inv_sqrt_ev .* transpose(pca_components)
            dat_sphered = Matrix{Float64}(undef, n_components, n_samples)
            mul!(dat_sphered, whiten_mat, dat, 1.0 / scale, 0.0)
            used_eigen = true
        end
    end

    if !used_eigen
        dat ./= scale
        F = svd!(dat)
        s_max = F.S[1]
        s_target = F.S[n_components]
        if s_max > 0 && (s_target / s_max < 1e-7 || s_target < 1e-10)
            eff_rank = count(s -> (s / s_max) > 1e-7, F.S)
            @minimal_error "Data is rank-deficient for requested $n_components components (singular value ratio = $(round(s_target / s_max, sigdigits=3))). Estimated rank is $eff_rank. Please set n_components <= $eff_rank (e.g. n_channels - 1 for average-referenced data)."
        end

        pca_components = F.U[:, 1:n_components]
        eigenvalues = (@view(F.S[1:n_components]) .^ 2) ./ (n_samples - 1)
        sphere = diagm(1.0 ./ sqrt.(eigenvalues))

        dat_sphered = Matrix{Float64}(undef, n_components, n_samples)
        whiten_scale = sqrt(n_samples - 1)
        @views @. dat_sphered = whiten_scale * F.Vt[1:n_components, :]
    end

    cpu_block = min(nextpow(2, max(Int(floor(sqrt(n_samples / 3.0))), 1)), 512)
    gpu_active = false
    if prms.use_gpu
        if is_gpu_available()
            gpu_active = true
            @info "[GPU ACTIVATED] Running ICA on $(gpu_device_name())..."
            gpu_block = min(nextpow(2, cpu_block * 4), 1024)
            block = min(n_samples, max(gpu_block, 32))
        else
            @minimal_warning "Requested GPU acceleration (use_gpu=true), but no functional GPU package (CUDA.jl, AMDGPU.jl, Metal.jl) has been loaded. Please run 'using CUDA' or 'using AMDGPU' before calling run_ica. Falling back to CPU."
            block = cpu_block
        end
    else
        block = cpu_block
    end
    work = create_work_arrays(n_components, block)

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
    act_power = (n_samples - 1) .* vec(sum(abs2, W_ica; dims = 2))
    unmixing = W_ica * sphere * transpose(pca_components)
    mixing = pinv(unmixing)

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

_ica_postprocess(
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
) = _ica_postprocess(work.weights, sphere, pca_components, scale, original_mean, n_samples, n_components, filename, layout, is_sub_gaussian)

"""
    infomax_ica(dat_ica, layout, filename; n_components, params = IcaPrms(), rng = Random.default_rng()) -> InfoIca

Standard natural-gradient Infomax ICA decomposition (super-Gaussian sources).
"""
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
        _infomax_optimize_gpu!(prep.work, prep.dat_sphered, prep.prms, prep.block, nothing, rng)
    else
        _infomax_optimize_cpu!(prep.work, prep.dat_sphered, prep.prms, prep.block, nothing, rng)
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

"""
    infomax_extended_ica(dat_ica, layout, filename; n_components, params = IcaPrms(), rng = Random.default_rng()) -> InfoIca

Extended Infomax ICA decomposition separating mixed super-Gaussian and sub-Gaussian sources via
kurtosis-based sign switching (Lee, Girolami, & Sejnowski, 1999).
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
        _infomax_optimize_gpu!(prep.work, prep.dat_sphered, prep.prms, prep.block, is_sub_gaussian, rng)
    else
        _infomax_optimize_cpu!(prep.work, prep.dat_sphered, prep.prms, prep.block, is_sub_gaussian, rng)
    end

    n_sub_final = count(is_sub_gaussian)
    @info "Extended Infomax completed: sup/sub-gauss ($(n_components - n_sub_final)/$n_sub_final)"

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

@inline function _update_extended_kurtosis_signs!(
    is_sub_gaussian::BitVector,
    old_kurtosis::Vector{Float64},
    u2_mean_sq_fn,
    u4_mean_fn,
    n_channels::Int,
)
    extmomentum = 0.5
    kurtosis_threshold = 0.05
    signs_changed = false

    @inbounds for i = 1:n_channels
        u2_msq = u2_mean_sq_fn(i)
        if u2_msq > eps(Float64)
            kurtosis_raw = (u4_mean_fn(i) / u2_msq) - 3.0
            kurtosis = old_kurtosis[i] != 0.0 ? (extmomentum * old_kurtosis[i] + (1.0 - extmomentum) * kurtosis_raw) : kurtosis_raw
            old_kurtosis[i] = kurtosis

            old_is_sub = is_sub_gaussian[i]
            new_is_sub = kurtosis < -kurtosis_threshold ? true : (kurtosis > kurtosis_threshold ? false : old_is_sub)
            is_sub_gaussian[i] = new_is_sub
            signs_changed |= (old_is_sub != new_is_sub)
        end
    end
    return signs_changed
end

# SIMD-vectorizable Float64 tanh via range reduction (|z| <= 1.25), [7/6] Padé, and 4 angle-doublings (< 4.5e-16 max error)
@inline function _simd_tanh(x::Float64)
    x_c = ifelse(x > 20.0, 20.0, ifelse(x < -20.0, -20.0, x))
    z = x_c * 0.0625
    z2 = z * z
    num = z * evalpoly(z2, (135135.0, 17325.0, 378.0, 1.0))
    den = evalpoly(z2, (135135.0, 62370.0, 3150.0, 28.0))
    t = num / den
    t = (t + t) / (1.0 + t * t)
    t = (t + t) / (1.0 + t * t)
    t = (t + t) / (1.0 + t * t)
    return (t + t) / (1.0 + t * t)
end

@inline function _infomax_cpu_block_step!(
    work::WorkArrays,
    u_blk::M1,
    y_blk::M2,
    d_blk::M3,
    block_size::Int,
    l_rate::Float64,
    signs::S,
) where {M1<:AbstractMatrix{Float64},M2<:AbstractMatrix{Float64},M3<:AbstractMatrix{Float64},S}
    n_channels = size(work.weights, 1)
    mul!(u_blk, work.weights, d_blk)
    if !isnothing(signs)
        @inbounds for j = 1:block_size
            @simd ivdep for i = 1:n_channels
                u_ij = u_blk[i, j]
                y_blk[i, j] = u_ij + signs[i] * _simd_tanh(u_ij)
            end
        end
        mul!(work.wu_term, y_blk, transpose(u_blk))
        mul!(work.weights_temp, work.wu_term, work.weights)
        @inbounds @. work.weights += l_rate * (block_size * work.weights - work.weights_temp)
    else
        @inbounds for j = 1:block_size
            @simd ivdep for i = 1:n_channels
                y_blk[i, j] = -_simd_tanh(0.5 * u_blk[i, j])
            end
        end
        mul!(work.wu_term, y_blk, transpose(u_blk))
        mul!(work.weights_temp, work.wu_term, work.weights)
        @inbounds @. work.weights += l_rate * (block_size * work.weights + work.weights_temp)
    end
    return nothing
end

function _infomax_optimize_cpu!(
    work::WorkArrays,
    dat_ica::Matrix{Float64},
    params::IcaPrms,
    block::Int,
    is_sub_gaussian::Union{Nothing,BitVector},
    rng::AbstractRNG = Random.default_rng(),
)
    n_channels, n_samples = size(dat_ica)
    extended = !isnothing(is_sub_gaussian)
    alg_name = extended ? "Extended Infomax" : "Infomax"

    step = 0
    wts_blowup = false
    change = 0.0
    oldchange = 0.0
    angledelta = 0.0
    converged = false

    permute_indices = Vector{Int}(undef, n_samples)

    kurtsize = extended ? min(6000, n_samples) : 0
    old_kurtosis = extended ? zeros(Float64, n_channels) : Float64[]
    signs = extended ? Float64[(is_sub_gaussian::BitVector)[i] ? -1.0 : 1.0 for i = 1:n_channels] : nothing
    kurt_data = (extended && kurtsize < n_samples) ? Matrix{Float64}(undef, n_channels, kurtsize) : Matrix{Float64}(undef, 0, 0)
    activations = extended ? Matrix{Float64}(undef, n_channels, kurtsize) : Matrix{Float64}(undef, 0, 0)
    n_super = n_channels
    n_sub = 0

    @inbounds while step < params.max_iter
        randperm!(rng, permute_indices)

        for t = 1:block:n_samples
            block_end = min(t + block - 1, n_samples)
            block_size = block_end - t + 1

            for j = 1:block_size
                idx = permute_indices[t+j-1]
                @simd for i = 1:n_channels
                    work.data_block[i, j] = dat_ica[i, idx]
                end
            end

            if block_size < block
                _infomax_cpu_block_step!(
                    work,
                    view(work.u, :, 1:block_size),
                    view(work.y, :, 1:block_size),
                    view(work.data_block, :, 1:block_size),
                    block_size,
                    params.l_rate,
                    signs,
                )
            else
                _infomax_cpu_block_step!(work, work.u, work.y, work.data_block, block_size, params.l_rate, signs)
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
            if params.l_rate < 1e-6
                @minimal_error "$alg_name failed to converge: learning rate fell below minimum floor (1e-6) due to repeated weight blow-ups. Check data for rank deficiency or artifacts."
            end
            work.weights .= work.startweights
            work.oldweights .= work.startweights
            if extended
                fill!(is_sub_gaussian::BitVector, false)
                fill!(old_kurtosis, 0.0)
                fill!(signs::Vector{Float64}, 1.0)
                n_super = n_channels
                n_sub = 0
            end
            continue
        end

        if extended && step > 10 && step % 10 == 0
            sub_flags = is_sub_gaussian::BitVector
            if kurtsize < n_samples
                for j = 1:kurtsize
                    idx = permute_indices[j]
                    @simd for k = 1:n_channels
                        kurt_data[k, j] = dat_ica[k, idx]
                    end
                end
                mul!(activations, work.weights, kurt_data)
            else
                mul!(activations, work.weights, dat_ica)
            end

            inv_k = 1.0 / kurtsize
            signs_changed = _update_extended_kurtosis_signs!(
                sub_flags,
                old_kurtosis,
                i -> begin
                    s2 = 0.0
                    @inbounds @simd for j = 1:kurtsize
                        v = activations[i, j]
                        s2 += v * v
                    end
                    (s2 * inv_k)^2
                end,
                i -> begin
                    s4 = 0.0
                    @inbounds @simd for j = 1:kurtsize
                        v2 = activations[i, j]^2
                        s4 += v2 * v2
                    end
                    s4 * inv_k
                end,
                n_channels,
            )
            if signs_changed
                signs_vec = signs::Vector{Float64}
                for i = 1:n_channels
                    signs_vec[i] = sub_flags[i] ? -1.0 : 1.0
                end
            end
            n_sub = count(sub_flags)
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
            if extended
                @info Printf.@sprintf(
                    "Extended-Infomax step %d, change = %.7f, lrate = %.7f, angle = %.1f, sup/sub-gauss: %d/%d",
                    step,
                    change,
                    params.l_rate,
                    params.degconst * angledelta,
                    n_super,
                    n_sub,
                )
            else
                @info Printf.@sprintf(
                    "Infomax step %d, change = %.7f, lrate = %.7f, angle = %.1f",
                    step,
                    change,
                    params.l_rate,
                    params.degconst * angledelta,
                )
            end
        end
    end

    if !converged
        @minimal_warning "$alg_name reached maximum iterations ($(params.max_iter)) without converging (final change = $change, target = $(params.w_change))."
    end
    return nothing
end

function _infomax_optimize_gpu!(
    work::WorkArrays,
    dat_ica::Matrix{Float64},
    params::IcaPrms,
    block::Int,
    is_sub_gaussian::Union{Nothing,BitVector},
    rng::AbstractRNG = Random.default_rng(),
)
    n_channels, n_samples = size(dat_ica)
    extended = !isnothing(is_sub_gaussian)
    alg_name = extended ? "Extended Infomax" : "Infomax"

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
    weights_f32_cpu = Matrix{Float32}(undef, n_channels, n_channels)
    startweights_f32 = Float32.(work.startweights)

    u_gpu = gpu_array(zeros(Float32, n_channels, block))
    y_gpu = gpu_array(zeros(Float32, n_channels, block))
    wu_gpu = gpu_array(zeros(Float32, n_channels, n_channels))
    wtemp_gpu = gpu_array(zeros(Float32, n_channels, n_channels))

    old_kurtosis = extended ? zeros(Float64, n_channels) : Float64[]
    signs_cpu = extended ? Float32[(is_sub_gaussian::BitVector)[i] ? -1.0f0 : 1.0f0 for i = 1:n_channels] : Float32[]
    signs_gpu = extended ? gpu_array(reshape(signs_cpu, n_channels, 1)) : nothing
    kurtsize = extended ? min(6000, n_samples) : 0
    activations_gpu = extended ? gpu_array(zeros(Float32, n_channels, kurtsize)) : nothing
    act_sq = extended ? gpu_array(zeros(Float32, n_channels, kurtsize)) : nothing
    u2_sum = extended ? gpu_array(zeros(Float32, n_channels, 1)) : nothing
    u4_sum = extended ? gpu_array(zeros(Float32, n_channels, 1)) : nothing
    u2_mean_sq_cpu = extended ? Vector{Float32}(undef, n_channels) : Float32[]
    u4_mean_cpu = extended ? Vector{Float32}(undef, n_channels) : Float32[]
    inv_kurtsize = extended ? 1.0f0 / Float32(kurtsize) : 0.0f0
    n_super = n_channels
    n_sub = 0

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

    block_starts = 1:block:n_samples
    data_block_views = [view(dat_ica_gpu_shuffled, :, t:min(t+block-1, n_samples)) for t in block_starts]
    n_blocks = length(block_starts)

    @inbounds while step < params.max_iter
        randperm!(rng, permute_indices)
        copyto!(perm_gpu, permute_indices)
        dat_ica_gpu_shuffled .= view(dat_ica_gpu, :, perm_gpu)
        l_rate_f32 = Float32(params.l_rate)

        for b_idx = 1:n_blocks
            data_block_gpu = data_block_views[b_idx]
            is_full = size(data_block_gpu, 2) == block
            u_v = is_full ? u_view_full : u_view_last
            y_v = is_full ? y_view_full : y_view_last
            u_vt = is_full ? u_view_full_t : u_view_last_t
            bs_f32 = is_full ? block_f32 : last_block_f32

            mul!(u_v, weights_gpu, data_block_gpu)
            if extended
                @. y_v = u_v + signs_gpu * tanh(u_v)
                mul!(wu_gpu, y_v, u_vt)
                mul!(wtemp_gpu, wu_gpu, weights_gpu)
                @. weights_gpu += l_rate_f32 * (bs_f32 * weights_gpu - wtemp_gpu)
            else
                @. y_v = -tanh(0.5f0 * u_v)
                mul!(wu_gpu, y_v, u_vt)
                mul!(wtemp_gpu, wu_gpu, weights_gpu)
                @. weights_gpu += l_rate_f32 * (bs_f32 * weights_gpu + wtemp_gpu)
            end
        end

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
                @minimal_error "$alg_name failed to converge on GPU: learning rate fell below minimum floor (1e-6) due to repeated weight blow-ups."
            end
            work.weights .= work.startweights
            work.oldweights .= work.startweights
            copyto!(weights_gpu, startweights_f32)
            if extended
                fill!(is_sub_gaussian::BitVector, false)
                fill!(old_kurtosis, 0.0)
                fill!(signs_cpu, 1.0f0)
                copyto!(signs_gpu, reshape(signs_cpu, n_channels, 1))
                n_super = n_channels
                n_sub = 0
            end
            continue
        end

        if extended && step > 10 && step % 10 == 0
            sub_flags = is_sub_gaussian::BitVector
            mul!(activations_gpu, weights_gpu, view(dat_ica_gpu_shuffled, :, 1:kurtsize))
            @. act_sq = activations_gpu^2
            sum!(u2_sum, act_sq)
            @. act_sq = act_sq^2
            sum!(u4_sum, act_sq)

            copyto!(u2_mean_sq_cpu, u2_sum)
            copyto!(u4_mean_cpu, u4_sum)
            @. u2_mean_sq_cpu = (u2_mean_sq_cpu * inv_kurtsize)^2
            @. u4_mean_cpu = u4_mean_cpu * inv_kurtsize

            signs_changed = _update_extended_kurtosis_signs!(
                sub_flags,
                old_kurtosis,
                i -> Float64(u2_mean_sq_cpu[i]),
                i -> Float64(u4_mean_cpu[i]),
                n_channels,
            )
            if signs_changed
                for i = 1:n_channels
                    signs_cpu[i] = sub_flags[i] ? -1.0f0 : 1.0f0
                end
                copyto!(signs_gpu, reshape(signs_cpu, n_channels, 1))
            end
            n_sub = count(sub_flags)
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
            if extended
                @info Printf.@sprintf(
                    "Extended-Infomax step %d, change = %.7f, lrate = %.7f, angle = %.1f, sup/sub-gauss: %d/%d",
                    step,
                    change,
                    params.l_rate,
                    params.degconst * angledelta,
                    n_super,
                    n_sub,
                )
            else
                @info Printf.@sprintf(
                    "Infomax step %d, change = %.7f, lrate = %.7f, angle = %.1f",
                    step,
                    change,
                    params.l_rate,
                    params.degconst * angledelta,
                )
            end
        end
    end

    if !converged
        @minimal_warning "$alg_name reached maximum iterations ($(params.max_iter)) without converging (final change = $change, target = $(params.w_change))."
    end

    copyto!(weights_f32_cpu, weights_gpu)
    work.weights .= weights_f32_cpu
    return nothing
end
