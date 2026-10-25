# === Amica ===
# === AMICA ICA IMPLEMENTATION ===

# Define learning rate struct to avoid Amica.jl dependency
mutable struct AmicaLearningRate{T<:Real}
    lrate::T
    lrate0::T
    min_lrate::T
    max_lrate::T
    lratefact::T
    shapelrate::T
    shapelrate0::T
    shapelratefact::T
    numdecs::Int
    maxdecs::Int
    newtrate::T
    newt_ramp::Int
    minrho::T
    maxrho::T
end

function AmicaLearningRate(
    T::Type{<:Real} = Float64;
    lrate0 = 0.1,
    min_lrate = 1e-8,
    max_lrate = 1.0,
    lratefact = 0.5,
    shapelrate0 = 0.05,
    shapelratefact = 0.5,
    maxdecs = 3,
    newtrate = 1.0,
    newt_ramp = 10,
    minrho = 1.0,
    maxrho = 2.0,
)
    return AmicaLearningRate{T}(
        T(lrate0),
        T(lrate0),
        T(min_lrate),
        T(max_lrate),
        T(lratefact),
        T(shapelrate0),
        T(shapelrate0),
        T(shapelratefact),
        0,
        maxdecs,
        T(newtrate),
        newt_ramp,
        T(minrho),
        T(maxrho),
    )
end

# Native Amica Work Arrays (Zero-Allocation)

mutable struct AmicaWorkArrays{T<:Real}
    g_times_sources::Array{T,2}
    sum_z::Array{T,2}
    kp::Array{T,2}
    dmu_numer::Array{T,2}
    dmu_denom::Array{T,2}
    dbeta_denom::Array{T,2}
    dlambda_numer::Array{T,2}
    drho_numer::Array{T,2}
    newton_sigma2::Array{T,1}
    Lt_accum::Array{T,1}

    dA::Array{T,2}
    Lt::Array{T,1}
    newton_kappa::Array{T,1}
    newton_lambda::Array{T,1}

    # Preallocated block arrays (distinct buffers to avoid aliasing)
    source_signals::Array{T,2}
    y::Array{T,3}
    y_rho::Array{T,3}
    Q::Array{T,3}
    z::Array{T,3}
    fp::Array{T,3}
    scratch3::Array{T,3}
    Qmax::Array{T,2}
    logexp::Array{T,2}
    zsum::Array{T,2}
    g_block_sum::Array{T,2}
end

function AmicaWorkArrays(T::Type{<:Real}, N::Int, n::Int, m::Int, block_size::Int)
    return AmicaWorkArrays{T}(
        zeros(T, n, n),
        zeros(T, n, m),
        zeros(T, n, m),
        zeros(T, n, m),
        zeros(T, n, m),
        zeros(T, n, m),
        zeros(T, n, m),
        zeros(T, n, m),
        zeros(T, n),
        zeros(T, N),
        zeros(T, n, n),
        zeros(T, N),
        zeros(T, n),
        zeros(T, n),
        zeros(T, block_size, n),
        zeros(T, block_size, n, m),
        zeros(T, block_size, n, m),
        zeros(T, block_size, n, m),
        zeros(T, block_size, n, m),
        zeros(T, block_size, n, m),
        zeros(T, block_size, n, m),
        zeros(T, block_size, n),
        zeros(T, block_size, n),
        zeros(T, block_size, n),
        zeros(T, block_size, n),
    )
end

function reset_accumulators!(wa::AmicaWorkArrays{T}) where {T}
    fill!(wa.g_times_sources, zero(T))
    fill!(wa.sum_z, zero(T))
    fill!(wa.kp, zero(T))
    fill!(wa.dmu_numer, zero(T))
    fill!(wa.dmu_denom, zero(T))
    fill!(wa.dbeta_denom, zero(T))
    fill!(wa.dlambda_numer, zero(T))
    fill!(wa.drho_numer, zero(T))
    fill!(wa.newton_sigma2, zero(T))
    fill!(wa.Lt_accum, zero(T))
end

function process_amica_blocks!(
    wa_threads::Vector{AmicaWorkArrays{T}},
    data::AbstractMatrix{T},
    W::AbstractMatrix{T},
    scale::AbstractMatrix{T},
    location::AbstractMatrix{T},
    shape::AbstractMatrix{T},
    proportions::AbstractMatrix{T},
    newton_active::Bool,
    block_size::Int,
) where {T<:Real}
    N, n = size(data)
    m = size(scale, 2)
    num_blocks = cld(N, block_size)
    n_threads = length(wa_threads)

    # Avoid oversubscription between Julia threads and BLAS
    old_blas_threads = LinearAlgebra.BLAS.get_num_threads()
    LinearAlgebra.BLAS.set_num_threads(1)
    try
        Threads.@threads for tid = 1:n_threads
            wa = wa_threads[tid]
            for block = tid:n_threads:num_blocks
                lower = (block - 1) * block_size + 1
                upper = min(N, block * block_size)
                n_samples = upper - lower + 1
                rng = lower:upper

                source_signals = view(wa.source_signals, 1:n_samples, :)
                y = view(wa.y, 1:n_samples, :, :)
                y_rho = view(wa.y_rho, 1:n_samples, :, :)
                Q = view(wa.Q, 1:n_samples, :, :)
                z = view(wa.z, 1:n_samples, :, :)
                fp = view(wa.fp, 1:n_samples, :, :)
                scratch3 = view(wa.scratch3, 1:n_samples, :, :)
                Qmax = view(wa.Qmax, 1:n_samples, :)
                logexp = view(wa.logexp, 1:n_samples, :)
                zsum = view(wa.zsum, 1:n_samples, :)
                g_block_sum = view(wa.g_block_sum, 1:n_samples, :)

                mul!(source_signals, view(data, rng, :), transpose(W))

                # Pass 1: compute y, y_rho, and Q in a single fused loop
                for j = 1:m
                    for i = 1:n
                        scale_ij = scale[i, j]
                        loc_ij = location[i, j]
                        shape_ij = shape[i, j]
                        exp_pow = shape_ij - T(1.0)
                        qconst = -log(T(2)) - SpecialFunctions.loggamma(T(1) + T(1) / shape_ij) + log(proportions[i, j]) + log(scale_ij)
                        @fastmath @inbounds @simd for s = 1:n_samples
                            y_val = scale_ij * (source_signals[s, i] - loc_ij)
                            y[s, i, j] = y_val
                            abs_y = abs(y_val)
                            val = ifelse(abs_y == zero(T), T(1e-16), abs_y)
                            y_rho_val = exp(exp_pow * log(val))
                            y_rho[s, i, j] = y_rho_val
                            Q[s, i, j] = qconst - (y_rho_val * abs_y)
                        end
                    end
                end

                maximum!(view(Qmax, :, :, (1:1)), Q)
                @fastmath @inbounds @simd for I in CartesianIndices((1:n_samples, 1:n, 1:m))
                    scratch3[I] = exp(Q[I] - Qmax[I[1], I[2]])
                end
                sum!(view(logexp, :, :, (1:1)), scratch3)
                @fastmath @inbounds @simd for I in CartesianIndices((1:n_samples, 1:n))
                    logexp[I] = log(logexp[I]) + Qmax[I]
                end

                @fastmath @inbounds @simd for I in CartesianIndices((1:n_samples, 1:n, 1:m))
                    z[I] = exp(Q[I] - logexp[I[1], I[2]]) + T(1e-15)
                end

                for i = 1:n
                    @fastmath @inbounds @simd for s = 1:n_samples
                        wa.Lt_accum[lower+s-1] += logexp[s, i]
                    end
                end

                sum!(view(zsum, :, :, (1:1)), z)
                @fastmath @inbounds @simd for I in CartesianIndices((1:n_samples, 1:n, 1:m))
                    z[I] = z[I] / zsum[I[1], I[2]]
                end

                # Pass 2: compute fp, scratch3, drho, and accumulate into local scalar registers
                for j = 1:m
                    for i = 1:n
                        shape_ij = shape[i, j]
                        scale_ij = scale[i, j]
                        is_shape_le_2 = (shape_ij <= T(2))
                        acc_drho = zero(T)
                        acc_sum_z = zero(T)
                        acc_dmu_numer = zero(T)
                        acc_kp = zero(T)
                        acc_dbeta_denom = zero(T)
                        acc_dmu_denom = zero(T)
                        acc_dlambda_numer = zero(T)

                        @fastmath @inbounds @simd for s = 1:n_samples
                            y_val = y[s, i, j]
                            y_rho_val = y_rho[s, i, j]
                            z_val = z[s, i, j]

                            fp_val = y_rho_val * sign(y_val) * shape_ij
                            fp[s, i, j] = fp_val
                            scratch3[s, i, j] = scale_ij * z_val * fp_val

                            sc_val = y_rho_val * abs(y_val)
                            sc = ifelse(sc_val >= T(1.0e-16), z_val * log(sc_val) * sc_val, zero(T))
                            acc_drho += sc

                            acc_sum_z += z_val
                            acc_dmu_numer += fp_val * z_val
                            acc_kp += fp_val * z_val * fp_val
                            # dbeta_denom must be accumulated unconditionally for all shapes
                            acc_dbeta_denom += fp_val * z_val * y_val

                            if is_shape_le_2
                                y_safe = ifelse(y_val == zero(T), T(1e-16), y_val)
                                acc_dmu_denom += (z_val * fp_val / y_safe) * scale_ij
                            else
                                acc_dmu_denom += (z_val * fp_val * fp_val) * scale_ij
                            end
                            acc_dlambda_numer += z_val * (fp_val * y_val - one(T))^2
                        end

                        wa.drho_numer[i, j] += acc_drho
                        wa.sum_z[i, j] += acc_sum_z
                        wa.dmu_numer[i, j] += acc_dmu_numer
                        wa.kp[i, j] += acc_kp
                        wa.dbeta_denom[i, j] += acc_dbeta_denom
                        wa.dmu_denom[i, j] += acc_dmu_denom
                        wa.dlambda_numer[i, j] += acc_dlambda_numer
                    end
                end

                sum!(view(g_block_sum, :, :, (1:1)), scratch3)
                mul!(wa.g_times_sources, transpose(g_block_sum), source_signals, one(T), one(T))

                if newton_active
                    for i = 1:n
                        acc_sig = zero(T)
                        @fastmath @inbounds @simd for s = 1:n_samples
                            acc_sig += source_signals[s, i]^2
                        end
                        wa.newton_sigma2[i] += acc_sig / N
                    end
                end
            end
        end
    finally
        LinearAlgebra.BLAS.set_num_threads(old_blas_threads)
    end
end

function update_amica_parameters!(
    wa_threads::Vector{AmicaWorkArrays{T}},
    data::AbstractMatrix{T},
    A::AbstractMatrix{T},
    LLdetS::T,
    scale::AbstractMatrix{T},
    location::AbstractMatrix{T},
    shape::AbstractMatrix{T},
    proportions::AbstractMatrix{T},
    newton_active::Bool,
    m::Int,
    upd_shape::Bool,
    shapelrate::T,
    minrho::T,
    maxrho::T,
    block_size::Int,
) where {T<:Real}
    N, n = size(data)
    ldet = -logabsdet(A)[1]

    for wa in wa_threads
        reset_accumulators!(wa)
    end

    W = inv(A)
    process_amica_blocks!(wa_threads, data, W, scale, location, shape, proportions, newton_active, block_size)

    # Reduce thread results into wa_threads[1]
    wa = wa_threads[1]
    for tid = 2:length(wa_threads)
        wa_t = wa_threads[tid]
        wa.Lt_accum .+= wa_t.Lt_accum
        wa.drho_numer .+= wa_t.drho_numer
        wa.g_times_sources .+= wa_t.g_times_sources
        wa.newton_sigma2 .+= wa_t.newton_sigma2
        wa.sum_z .+= wa_t.sum_z
        wa.dmu_numer .+= wa_t.dmu_numer
        wa.kp .+= wa_t.kp
        wa.dmu_denom .+= wa_t.dmu_denom
        wa.dbeta_denom .+= wa_t.dbeta_denom
        wa.dlambda_numer .+= wa_t.dlambda_numer
    end

    wa.Lt .= ldet .+ LLdetS
    wa.Lt .+= wa.Lt_accum
    wa.dA .= I(n) - wa.g_times_sources ./ N
    LL_iter = sum(wa.Lt) / (N * n)

    if m > 1
        # Avoid collapsed mixture weights: apply floor and re-normalize across components
        min_prop = T(1e-4)
        raw_prop = max.(wa.sum_z ./ N, min_prop)
        proportions .= raw_prop ./ sum(raw_prop, dims = 2)
    end

    if newton_active
        dkap = (wa.kp ./ (proportions .* N)) .* scale .^ 2
        wa.newton_kappa .= dropdims(sum(proportions .* dkap, dims = 2), dims = 2)
        wa.newton_lambda .= dropdims(sum(proportions .* (wa.dlambda_numer ./ wa.sum_z .+ dkap .* location .^ 2), dims = 2), dims = 2)
    end

    if m > 1
        dmu_denom_safe = ifelse.(abs.(wa.dmu_denom) .< T(1e-15), sign.(wa.dmu_denom) .* T(1e-15), wa.dmu_denom)
        location .+= clamp.(wa.dmu_numer ./ dmu_denom_safe, T(-10.0), T(10.0))
    end

    scale_ratio = sqrt.(max.(wa.sum_z, T(1e-15)) ./ max.(wa.dbeta_denom, T(1e-15)))
    scale .*= clamp.(scale_ratio, T(0.1), T(10.0))

    if upd_shape
        shape .= clamp.(
            shape .+ (shapelrate .* (1 .- (shape ./ SpecialFunctions.digamma.(1 .+ 1 ./ shape)) .* wa.drho_numer ./ wa.sum_z)),
            minrho,
            maxrho,
        )
    end

    return LL_iter
end

function update_amica_mixing!(
    A::AbstractMatrix{T},
    dA::AbstractMatrix{T},
    newton_kappa::AbstractVector{T},
    newton_lambda::AbstractVector{T},
    newton_sigma2::AbstractVector{T},
    iter::Int,
    do_newton::Bool,
    no_newton::Bool,
    newt_start_iter::Int,
    lrate::AmicaLearningRate{T},
) where {T<:Real}
    n = size(A, 1)

    if (do_newton && !no_newton && iter >= newt_start_iter)
        if iter == newt_start_iter
            lrate.numdecs = 0
        end

        B = similar(dA)
        posdef = true
        for i = 1:n
            for k = 1:n
                if i == k
                    if isfinite(newton_lambda[i]) && abs(newton_lambda[i]) > eps(T)
                        B[i, k] = dA[i, k] / newton_lambda[i]
                    else
                        B[i, k] = zero(T)
                        posdef = false
                    end
                else
                    sk1 = newton_sigma2[i] * newton_kappa[k]
                    sk2 = newton_sigma2[k] * newton_kappa[i]
                    denom = sk1 * sk2 - one(T)
                    if isfinite(sk1) && isfinite(sk2) && isfinite(denom) && denom > eps(T)
                        B[i, k] = (sk1 * dA[i, k] - dA[k, i]) / denom
                    else
                        B[i, k] = zero(T)
                        posdef = false
                    end
                end
            end
        end

        if posdef
            lrate.lrate = min(lrate.newtrate, lrate.lrate + min(T(1.0) / lrate.newt_ramp, lrate.lrate))
            A .-= lrate.lrate .* A * B
            return no_newton
        else
            no_newton = true
            lrate.lrate = min(lrate.lrate0, lrate.lrate + min(T(1.0) / lrate.newt_ramp, lrate.lrate))
            A .-= lrate.lrate .* A * dA
            return no_newton
        end
    else
        lrate.lrate = min(lrate.lrate0, lrate.lrate + min(T(1.0) / lrate.newt_ramp, lrate.lrate))
        A .-= lrate.lrate .* A * dA
        return no_newton
    end
end

function reparameterize_amica!(A::AbstractMatrix{T}, location::AbstractMatrix{T}, scale::AbstractMatrix{T}) where {T<:Real}
    tau = dropdims(sqrt.(sum(A .^ 2, dims = 1)), dims = 1)
    mask = tau .> zero(T)
    A .= ifelse.(reshape(mask, 1, size(mask)...), A ./ tau', A)
    location .= ifelse.(mask, location .* tau, location)
    scale .= ifelse.(mask, scale ./ tau, scale)
end

function calculate_amica_lrate!(
    LL::AbstractVector,
    iter::Int,
    newt_start_iter::Int,
    do_newton::Bool,
    lrate::AmicaLearningRate{T},
) where {T<:Real}
    if LL[iter] < LL[iter-1]
        if lrate.lrate <= lrate.min_lrate
            return true
        else
            lrate.lrate *= lrate.lratefact
            lrate.shapelrate *= lrate.shapelratefact
            lrate.numdecs += 1
            if lrate.numdecs >= lrate.maxdecs
                lrate.lrate0 *= lrate.lratefact
                if iter > newt_start_iter
                    lrate.shapelrate0 *= lrate.shapelratefact
                end
                if do_newton && (iter > newt_start_iter)
                    lrate.newtrate *= lrate.lratefact
                end
                lrate.numdecs = 0
            end
        end
    end
    return false
end

"""
    amica_ica(dat_ica::Matrix{Float64}, layout::Layout, filename::String; n_components::Int, params::IcaPrms = IcaPrms(), rng::AbstractRNG = Random.default_rng())

Run the AMICA (Adaptive Mixture Independent Component Analysis) algorithm.
AMICA models source distributions using flexible mixtures of generalized Gaussians, adaptively
estimating component location, scale, shape (sub-/super-Gaussianity), and mixing proportions.

Includes second-order Newton acceleration, parameter snapshot rollback on non-finite log-likelihood,
SIMD register accumulation, and multi-threaded block processing.
"""
function amica_ica(
    dat_ica::Matrix{Float64},
    layout::Layout,
    filename::String;
    n_components::Int,
    params::IcaPrms = IcaPrms(),
    rng::AbstractRNG = Random.default_rng(),
)
    # Preprocessing: dimension checks, sample size warning, non-mutating copy, rank deficiency detection, analytic sphering
    prep = _ica_preprocess(dat_ica, n_components, params)

    n_samples = prep.n_samples

    # AMICA native math expects data as (N_samples, n_components)
    dat_ica_t = Matrix{Float64}(prep.dat_sphered')

    LLdetS = logabsdet(prep.sphere)[1]

    # AMICA parameters with dedicated fields from IcaPrms
    m = prep.prms.amica_m > 0 ? prep.prms.amica_m : 3
    maxiter = prep.prms.amica_max_iter > 0 ? prep.prms.amica_max_iter : 2000
    lrate_val = prep.prms.amica_lrate > 0 ? prep.prms.amica_lrate : 0.1
    lrate = AmicaLearningRate(Float64; lrate0 = lrate_val)

    # Dynamic block size based on N
    cpu_block = min(nextpow(2, max(Int(floor(sqrt(n_samples / 3.0))), 1)), 512)
    block_size = min(n_samples, max(cpu_block, 32))

    # Initialize A: small random ±0.005, diagonal = 1.0, then normalize columns
    Wtmp = rand(rng, Float64, n_components, n_components)
    A = 0.01 .* (0.5 .- Wtmp)
    for i = 1:n_components
        A[i, i] = 1.0
        A[:, i] = A[:, i] / norm(A[:, i])
    end

    proportions = fill(1.0 / m, n_components, m)
    shape = fill(1.5, n_components, m)

    # Initialize location centered around 0 + small noise ±0.05
    location = zeros(Float64, n_components, m)
    for j = 1:m
        location[:, j] .= Float64(j - 1 - (m - 1) / 2)
    end
    location .+= 0.05 .* (1.0 .- 2.0 .* rand(rng, Float64, n_components, m))

    # Initialize scale in range [0.95, 1.05]
    scale = ones(Float64, n_components, m) .+ 0.1 .* (0.5 .- rand(rng, Float64, n_components, m))

    n_threads = Threads.nthreads()
    wa_threads = [AmicaWorkArrays(Float64, n_samples, n_components, m, block_size) for _ = 1:n_threads]
    LL = fill(-Inf, maxiter)

    do_newton = true
    no_newton = false
    newt_start_iter = 50
    update_shape = true

    snap() = (copy(A), copy(scale), copy(location), copy(shape), copy(proportions))
    restore!(s) = begin
        copyto!(A, s[1])
        copyto!(scale, s[2])
        copyto!(location, s[3])
        copyto!(shape, s[4])
        copyto!(proportions, s[5])
    end

    best_ll = -Inf
    best_snap = snap()
    converged = false
    consecutive_conv = 0
    tol = prep.prms.default_stop > 0 ? prep.prms.default_stop : 1e-7
    last_iter = 0

    if prep.gpu_active
        A_final = _amica_optimize_gpu!(dat_ica_t, A, LLdetS, scale, location, shape, proportions, m, prep.prms, lrate)
        A .= A_final
    else
        for iter = 1:maxiter
            last_iter = iter
            cur_snap = snap()

            LL_iter = update_amica_parameters!(
                wa_threads,
                dat_ica_t,
                A,
                LLdetS,
                scale,
                location,
                shape,
                proportions,
                do_newton && iter >= newt_start_iter,
                m,
                update_shape,
                lrate.shapelrate,
                lrate.minrho,
                lrate.maxrho,
                block_size,
            )

            if !isfinite(LL_iter) || isnan(LL_iter)
                @minimal_warning "AMICA encountered non-finite log-likelihood ($LL_iter) at iteration $iter. Rolling back to best model."
                restore!(best_snap)
                break
            end

            if LL_iter > best_ll
                best_ll = LL_iter
                best_snap = cur_snap
            end

            LL[iter] = LL_iter

            wa = wa_threads[1]
            no_newton = update_amica_mixing!(
                A,
                wa.dA,
                wa.newton_kappa,
                wa.newton_lambda,
                wa.newton_sigma2,
                iter,
                do_newton,
                no_newton,
                newt_start_iter,
                lrate,
            )

            reparameterize_amica!(A, location, scale)

            if iter > 1
                if calculate_amica_lrate!(LL, iter, newt_start_iter, do_newton, lrate)
                    @info "AMICA learning rate reached minimum ($(lrate.min_lrate)) at iteration $iter. Stopping."
                    converged = true
                    break
                end

                if abs(LL[iter] - LL[iter-1]) < tol
                    consecutive_conv += 1
                    if consecutive_conv >= 3
                        @info "AMICA converged at iteration $iter (|ΔLL| < $tol for 3 consecutive steps)."
                        converged = true
                        break
                    end
                else
                    consecutive_conv = 0
                end
            end

            if iter == 1 || iter % 10 == 0
                @info Printf.@sprintf("AMICA iter %d, LL = %.7f", iter, LL[iter])
            end
        end

        # Ensure parameters correspond to the best model if not converged or last step degraded
        if best_ll > -Inf && (!converged || (last_iter > 0 && LL[last_iter] < best_ll))
            restore!(best_snap)
        end
    end

    W_final = inv(A)

    # Compute unmixed activations on whitened data to evaluate excess kurtosis for sub-Gaussian identification
    U = W_final * prep.dat_sphered
    is_sub_gaussian = falses(n_components)
    for i = 1:n_components
        u_i = view(U, i, :)
        m_i = mean(u_i)
        v_i = mean(abs2, u_i .- m_i)
        if v_i > 1e-12
            m4_i = mean(x -> (x - m_i)^4, u_i)
            kurt = (m4_i / (v_i^2)) - 3.0
            is_sub_gaussian[i] = (kurt < -0.05)
        end
    end

    return _ica_postprocess(
        W_final,
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

# === GPU versions of Amica ===

function _amica_optimize_gpu!(
    dat_ica_t::Matrix{Float64},
    A_cpu::Matrix{Float64},
    LLdetS::Float64,
    scale_cpu::Matrix{Float64},
    location_cpu::Matrix{Float64},
    shape_cpu::Matrix{Float64},
    proportions_cpu::Matrix{Float64},
    m::Int,
    params::IcaPrms,
    lrate::AmicaLearningRate{Float64},
)
    N, n = size(dat_ica_t)

    dat_gpu = gpu_array(Float32.(dat_ica_t))
    A = gpu_array(Float32.(A_cpu))

    # Preallocate GPU full-batch arrays to avoid 400MB+ allocations per iteration
    y = gpu_array(zeros(Float32, N, n, m))
    y_rho = gpu_array(zeros(Float32, N, n, m))
    Q = gpu_array(zeros(Float32, N, n, m))
    z = gpu_array(zeros(Float32, N, n, m))
    fp = gpu_array(zeros(Float32, N, n, m))
    scratch3 = gpu_array(zeros(Float32, N, n, m))
    sc = gpu_array(zeros(Float32, N, n, m))

    Qmax = gpu_array(zeros(Float32, N, n, 1))
    logexp = gpu_array(zeros(Float32, N, n, 1))
    zsum = gpu_array(zeros(Float32, N, n, 1))
    g_block_sum = gpu_array(zeros(Float32, N, n, 1))

    scale_z_fp = gpu_array(zeros(Float32, N, n, m))
    fp_z_y = gpu_array(zeros(Float32, N, n, m))
    dlambda_tmp = gpu_array(zeros(Float32, N, n, m))
    dmu_denom_1 = gpu_array(zeros(Float32, N, n, m))
    dmu_denom_2 = gpu_array(zeros(Float32, N, n, m))

    qconst_cpu = zeros(Float32, 1, n, m)
    source_signals = gpu_array(zeros(Float32, N, n))

    scale = gpu_array(Float32.(scale_cpu))
    location = gpu_array(Float32.(location_cpu))
    shape = gpu_array(Float32.(shape_cpu))
    proportions = gpu_array(Float32.(proportions_cpu))

    scale_3d = reshape(scale, 1, n, m)
    location_3d = reshape(location, 1, n, m)
    shape_3d = reshape(shape, 1, n, m)

    dA = gpu_array(zeros(Float32, n, n))

    maxiter = params.amica_max_iter > 0 ? params.amica_max_iter : 2000
    LL = fill(-Inf, maxiter)
    do_newton = true
    no_newton = false
    newt_start_iter = 50
    update_shape = true

    snap_gpu() = (copy(A_cpu), copy(scale_cpu), copy(location_cpu), copy(shape_cpu), copy(proportions_cpu))
    restore_gpu!(s) = begin
        copyto!(A_cpu, s[1])
        copyto!(scale_cpu, s[2])
        copyto!(location_cpu, s[3])
        copyto!(shape_cpu, s[4])
        copyto!(proportions_cpu, s[5])
    end

    best_ll = -Inf
    best_snap = snap_gpu()
    converged = false
    consecutive_conv = 0
    tol = params.default_stop > 0 ? params.default_stop : 1e-7
    last_iter = 0

    for iter = 1:maxiter
        last_iter = iter
        cur_snap = snap_gpu()

        W_c = inv(A_cpu)
        W = gpu_array(Float32.(W_c))

        mul!(source_signals, dat_gpu, transpose(W))

        for j = 1:m, i = 1:n
            qconst_cpu[1, i, j] =
                -log(2.0f0) - SpecialFunctions.loggamma(1.0f0 + 1.0f0 / Float32(shape_cpu[i, j])) +
                log(Float32(proportions_cpu[i, j])) +
                log(Float32(scale_cpu[i, j]))
        end
        qconst = gpu_array(qconst_cpu)

        ss_3d = reshape(source_signals, N, n, 1)

        @. y = scale_3d * (ss_3d - location_3d)
        @. y_rho = exp((shape_3d - 1.0f0) * log(max(abs(y), 1.0f-16)))
        @. Q = qconst - (y_rho * abs(y))

        maximum!(Qmax, Q)
        @. scratch3 = exp(Q - Qmax)
        sum!(logexp, scratch3)
        @. logexp = log(logexp) + Qmax

        @. z = exp(Q - logexp) + 1.0f-15
        sum!(zsum, z)
        @. z = z / zsum

        @. fp = y_rho * sign(y) * shape_3d

        @. sc = y_rho * abs(y)
        @. sc = ifelse(sc >= 1.0f-16, z * log(sc) * sc, 0.0f0)
        drho_numer = dropdims(sum(sc, dims = 1), dims = 1)

        @. scale_z_fp = scale_3d * z * fp
        sum!(g_block_sum, scale_z_fp)
        g_times_sources = transpose(dropdims(g_block_sum, dims = 3)) * source_signals

        newton_sigma2 = dropdims(sum(source_signals .^ 2, dims = 1), dims = 1) ./ Float32(N)

        sum_z = dropdims(sum(z, dims = 1), dims = 1)
        dmu_numer = dropdims(sum(fp .* z, dims = 1), dims = 1)
        kp = dropdims(sum(fp .* z .* fp, dims = 1), dims = 1)

        mask = shape_3d .<= 2.0f0
        y_safe = @. ifelse(y == 0.0f0, 1.0f-16, y)

        @. dmu_denom_1 = (z * fp / y_safe) * scale_3d
        @. dmu_denom_2 = (z * fp * fp) * scale_3d
        dmu_denom = dropdims(sum(ifelse.(mask, dmu_denom_1, dmu_denom_2), dims = 1), dims = 1)

        # Unconditional accumulation for dbeta_denom
        @. fp_z_y = fp * z * y
        dbeta_denom = dropdims(sum(fp_z_y, dims = 1), dims = 1)

        @. dlambda_tmp = z * (fp * y - 1.0f0)^2
        dlambda_numer = dropdims(sum(dlambda_tmp, dims = 1), dims = 1)

        ldet = -logabsdet(A_cpu)[1]
        # Sum columns on GPU, copy (1, n, 1) to CPU and sum in Float64 to eliminate Float32 summation noise
        sum_logexp_cols = sum(logexp, dims = 1)
        Lt_sum = sum(Float64, Array(vec(sum_logexp_cols)))
        LL_iter = (Float64(ldet + LLdetS) * Float64(N) + Lt_sum) / Float64(N * n)

        if !isfinite(LL_iter) || isnan(LL_iter)
            @minimal_warning "AMICA GPU encountered non-finite log-likelihood ($LL_iter) at iteration $iter. Rolling back to best model."
            restore_gpu!(best_snap)
            break
        end

        if LL_iter > best_ll
            best_ll = LL_iter
            best_snap = cur_snap
        end

        LL[iter] = LL_iter

        @. dA = -g_times_sources / Float32(N)
        dA[diagind(dA)] .+= 1.0f0

        # Convert accumulators to Float64
        sum_z_c = Float64.(Array(sum_z))
        dmu_numer_c = Float64.(Array(dmu_numer))
        dmu_denom_c = Float64.(Array(dmu_denom))
        dbeta_denom_c = Float64.(Array(dbeta_denom))
        drho_numer_c = Float64.(Array(drho_numer))
        kp_c = Float64.(Array(kp))
        dlambda_numer_c = Float64.(Array(dlambda_numer))
        newton_sigma2_c = Float64.(Array(newton_sigma2))

        if m > 1
            min_prop = 1.0e-4
            raw_prop = max.(sum_z_c ./ Float64(N), min_prop)
            proportions_cpu .= raw_prop ./ sum(raw_prop, dims = 2)
        end

        newton_active = do_newton && iter >= newt_start_iter
        if newton_active
            dkap = @. (kp_c / (proportions_cpu * N)) * scale_cpu^2
            newton_kappa_c = vec(sum(proportions_cpu .* dkap, dims = 2))
            newton_lambda_c = vec(sum(proportions_cpu .* (dlambda_numer_c ./ sum_z_c .+ dkap .* location_cpu .^ 2), dims = 2))
        else
            newton_kappa_c = zeros(Float64, n)
            newton_lambda_c = zeros(Float64, n)
        end

        if m > 1
            dmu_denom_safe = ifelse.(abs.(dmu_denom_c) .< 1e-15, sign.(dmu_denom_c) .* 1e-15, dmu_denom_c)
            location_cpu .+= clamp.(dmu_numer_c ./ dmu_denom_safe, -10.0, 10.0)
        end

        scale_ratio = sqrt.(max.(sum_z_c, 1e-15) ./ max.(dbeta_denom_c, 1e-15))
        scale_cpu .*= clamp.(scale_ratio, 0.1, 10.0)

        if update_shape
            @. shape_cpu = clamp(
                shape_cpu +
                (lrate.shapelrate * (1.0 - (shape_cpu / SpecialFunctions.digamma(1.0 + 1.0 / shape_cpu)) * drho_numer_c / sum_z_c)),
                lrate.minrho,
                lrate.maxrho,
            )
        end

        scale .= gpu_array(Float32.(scale_cpu))
        location .= gpu_array(Float32.(location_cpu))
        shape .= gpu_array(Float32.(shape_cpu))
        proportions .= gpu_array(Float32.(proportions_cpu))

        dA_c = Float64.(Array(dA))
        no_newton = update_amica_mixing!(
            A_cpu,
            dA_c,
            newton_kappa_c,
            newton_lambda_c,
            newton_sigma2_c,
            iter,
            do_newton,
            no_newton,
            newt_start_iter,
            lrate,
        )

        reparameterize_amica!(A_cpu, location_cpu, scale_cpu)

        A .= gpu_array(Float32.(A_cpu))
        location .= gpu_array(Float32.(location_cpu))
        scale .= gpu_array(Float32.(scale_cpu))

        if iter > 1
            if calculate_amica_lrate!(LL, iter, newt_start_iter, do_newton, lrate)
                @info "AMICA GPU learning rate reached minimum ($(lrate.min_lrate)) at iteration $iter. Stopping."
                converged = true
                break
            end
            if abs(LL[iter] - LL[iter-1]) < tol
                consecutive_conv += 1
                if consecutive_conv >= 3
                    @info "AMICA GPU converged at iteration $iter (|ΔLL| < $tol for 3 consecutive steps)."
                    converged = true
                    break
                end
            else
                consecutive_conv = 0
            end
        end

        if iter == 1 || iter % 10 == 0
            @info Printf.@sprintf("AMICA iter %d, LL = %.7f", iter, LL[iter])
        end
    end

    if best_ll > -Inf && (!converged || (last_iter > 0 && LL[last_iter] < best_ll))
        restore_gpu!(best_snap)
    end

    return A_cpu
end
