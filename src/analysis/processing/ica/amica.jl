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
    step_halved::Bool
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
        false,
    )
end

struct _AmicaSnapshot{T<:Real}
    A::Matrix{T}
    scale::Matrix{T}
    location::Matrix{T}
    shape::Matrix{T}
    proportions::Matrix{T}
end

_make_amica_snap(A::Matrix{T}, scale::Matrix{T}, location::Matrix{T}, shape::Matrix{T}, proportions::Matrix{T}) where {T<:Real} =
    _AmicaSnapshot{T}(copy(A), copy(scale), copy(location), copy(shape), copy(proportions))

@inline function _save_amica_snap!(
    snap::_AmicaSnapshot{T},
    A::Matrix{T},
    scale::Matrix{T},
    location::Matrix{T},
    shape::Matrix{T},
    proportions::Matrix{T},
) where {T<:Real}
    copyto!(snap.A, A)
    copyto!(snap.scale, scale)
    copyto!(snap.location, location)
    copyto!(snap.shape, shape)
    copyto!(snap.proportions, proportions)
    return nothing
end

@inline function _restore_amica_snap!(
    A::Matrix{T},
    scale::Matrix{T},
    location::Matrix{T},
    shape::Matrix{T},
    proportions::Matrix{T},
    snap::_AmicaSnapshot{T},
) where {T<:Real}
    copyto!(A, snap.A)
    copyto!(scale, snap.scale)
    copyto!(location, snap.location)
    copyto!(shape, snap.shape)
    copyto!(proportions, snap.proportions)
    return nothing
end

mutable struct AmicaWorkArrays{T<:Real}
    g_times_sources::Matrix{T}
    sum_z::Matrix{T}
    kp::Matrix{T}
    dmu_numer::Matrix{T}
    dmu_denom::Matrix{T}
    dbeta_denom::Matrix{T}
    dlambda_numer::Matrix{T}
    drho_numer::Matrix{T}
    newton_sigma2::Vector{T}
    ll_sum::T

    dA::Matrix{T}
    newton_kappa::Vector{T}
    newton_lambda::Vector{T}
    B_buf::Matrix{T}
    A_step::Matrix{T}

    source_signals::Matrix{T}
    y::Array{T,3}
    y_rho::Array{T,3}
    Q::Array{T,3}
    z::Array{T,3}
    fp::Array{T,3}
    scratch3::Array{T,3}
    Qmax::Matrix{T}
    logexp::Matrix{T}
    zsum::Matrix{T}
    g_block_sum::Matrix{T}
end

function AmicaWorkArrays(T::Type{<:Real}, n::Int, m::Int, block_size::Int)
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
        zero(T),
        zeros(T, n, n),
        zeros(T, n),
        zeros(T, n),
        zeros(T, n, n),
        zeros(T, n, n),
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
    wa.ll_sum = zero(T)
    return nothing
end

function _process_amica_worker_blocks!(
    wa::AmicaWorkArrays{T},
    tid::Int,
    n_workers::Int,
    num_blocks::Int,
    data::AbstractMatrix{T},
    W::AbstractMatrix{T},
    qconsts::Matrix{T},
    scale::AbstractMatrix{T},
    location::AbstractMatrix{T},
    shape::AbstractMatrix{T},
    newton_active::Bool,
    block_size::Int,
) where {T<:Real}
    N, n = size(data)
    m = size(scale, 2)

    for block = tid:n_workers:num_blocks
        lower = (block - 1) * block_size + 1
        upper = min(N, block * block_size)
        n_samples = upper - lower + 1

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

        mul!(source_signals, view(data, lower:upper, :), transpose(W))

        @inbounds for j = 1:m, i = 1:n
            scale_ij = scale[i, j]
            loc_ij = location[i, j]
            exp_pow = shape[i, j] - T(1.0)
            qconst = qconsts[i, j]
            @fastmath @simd for s = 1:n_samples
                y_val = scale_ij * (source_signals[s, i] - loc_ij)
                y[s, i, j] = y_val
                abs_y = abs(y_val)
                val = ifelse(abs_y == zero(T), T(1e-16), abs_y)
                y_rho_val = exp(exp_pow * log(val))
                y_rho[s, i, j] = y_rho_val
                Q[s, i, j] = qconst - (y_rho_val * abs_y)
            end
        end

        maximum!(view(Qmax, :, :, 1:1), Q)
        @fastmath @inbounds @simd for I in CartesianIndices((1:n_samples, 1:n, 1:m))
            scratch3[I] = exp(Q[I] - Qmax[I[1], I[2]])
        end
        sum!(view(logexp, :, :, 1:1), scratch3)

        acc_ll = zero(T)
        @fastmath @inbounds for i = 1:n
            @simd for s = 1:n_samples
                le = log(logexp[s, i]) + Qmax[s, i]
                logexp[s, i] = le
                acc_ll += le
            end
        end
        wa.ll_sum += acc_ll

        @fastmath @inbounds @simd for I in CartesianIndices((1:n_samples, 1:n, 1:m))
            z[I] = exp(Q[I] - logexp[I[1], I[2]]) + T(1e-15)
        end
        sum!(view(zsum, :, :, 1:1), z)
        @fastmath @inbounds @simd for I in CartesianIndices((1:n_samples, 1:n, 1:m))
            z[I] /= zsum[I[1], I[2]]
        end

        @inbounds for j = 1:m, i = 1:n
            shape_ij = shape[i, j]
            scale_ij = scale[i, j]
            is_shape_le_2 = shape_ij <= T(2)
            acc_drho = zero(T)
            acc_sum_z = zero(T)
            acc_dmu_numer = zero(T)
            acc_kp = zero(T)
            acc_dbeta_denom = zero(T)
            acc_dmu_denom = zero(T)
            acc_dlambda_numer = zero(T)

            @fastmath @simd for s = 1:n_samples
                y_val = y[s, i, j]
                y_rho_val = y_rho[s, i, j]
                z_val = z[s, i, j]

                fp_val = y_rho_val * sign(y_val) * shape_ij
                fp[s, i, j] = fp_val
                scratch3[s, i, j] = scale_ij * z_val * fp_val

                sc_val = y_rho_val * abs(y_val)
                acc_drho += ifelse(sc_val >= T(1.0e-16), z_val * log(sc_val) * sc_val, zero(T))
                acc_sum_z += z_val
                acc_dmu_numer += fp_val * z_val
                acc_kp += fp_val * z_val * fp_val
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

        sum!(view(g_block_sum, :, :, 1:1), scratch3)
        mul!(wa.g_times_sources, transpose(g_block_sum), source_signals, one(T), one(T))

        if newton_active
            @inbounds for i = 1:n
                acc_sig = zero(T)
                @fastmath @simd for s = 1:n_samples
                    acc_sig += source_signals[s, i]^2
                end
                wa.newton_sigma2[i] += acc_sig / N
            end
        end
    end
    return nothing
end

function process_amica_blocks!(
    wa_threads::Vector{AmicaWorkArrays{T}},
    data::AbstractMatrix{T},
    W::AbstractMatrix{T},
    qconsts::Matrix{T},
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
    n_workers = min(length(wa_threads), num_blocks)

    @inbounds for j = 1:m, i = 1:n
        qconsts[i, j] = -log(T(2)) - SpecialFunctions.loggamma(T(1) + T(1) / shape[i, j]) + log(proportions[i, j]) + log(scale[i, j])
    end

    if n_workers == 1
        _process_amica_worker_blocks!(wa_threads[1], 1, 1, num_blocks, data, W, qconsts, scale, location, shape, newton_active, block_size)
    else
        old_blas_threads = LinearAlgebra.BLAS.get_num_threads()
        LinearAlgebra.BLAS.set_num_threads(1)
        try
            Threads.@threads for tid = 1:n_workers
                _process_amica_worker_blocks!(
                    wa_threads[tid],
                    tid,
                    n_workers,
                    num_blocks,
                    data,
                    W,
                    qconsts,
                    scale,
                    location,
                    shape,
                    newton_active,
                    block_size,
                )
            end
        finally
            LinearAlgebra.BLAS.set_num_threads(old_blas_threads)
        end
    end
    return n_workers
end

function update_amica_parameters!(
    wa_threads::Vector{AmicaWorkArrays{T}},
    qconsts::Matrix{T},
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

    for wa_i in wa_threads
        reset_accumulators!(wa_i)
    end

    W = inv(A)
    n_workers = process_amica_blocks!(wa_threads, data, W, qconsts, scale, location, shape, proportions, newton_active, block_size)

    wa = wa_threads[1]
    for tid = 2:n_workers
        wa_t = wa_threads[tid]
        wa.ll_sum += wa_t.ll_sum
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

    inv_N = one(T) / T(N)
    @inbounds for j = 1:n, i = 1:n
        wa.dA[i, j] = (i == j ? one(T) : zero(T)) - wa.g_times_sources[i, j] * inv_N
    end
    LL_iter = (T(N) * (ldet + LLdetS) + wa.ll_sum) / T(N * n)

    if m > 1
        min_prop = T(1e-4)
        @inbounds for i = 1:n
            row_sum = zero(T)
            for j = 1:m
                p = max(wa.sum_z[i, j] * inv_N, min_prop)
                proportions[i, j] = p
                row_sum += p
            end
            inv_row = one(T) / row_sum
            for j = 1:m
                proportions[i, j] *= inv_row
            end
        end
    end

    if newton_active
        @inbounds for i = 1:n
            kap_i = zero(T)
            lam_i = zero(T)
            for j = 1:m
                prop_ij = proportions[i, j]
                dkap = (wa.kp[i, j] / (prop_ij * T(N))) * scale[i, j]^2
                kap_i += prop_ij * dkap
                lam_i += prop_ij * (wa.dlambda_numer[i, j] / wa.sum_z[i, j] + dkap * location[i, j]^2)
            end
            wa.newton_kappa[i] = kap_i
            wa.newton_lambda[i] = lam_i
        end
    end

    if m > 1
        @inbounds for j = 1:m, i = 1:n
            denom = wa.dmu_denom[i, j]
            denom_safe = abs(denom) < T(1e-15) ? sign(denom) * T(1e-15) : denom
            location[i, j] += clamp(wa.dmu_numer[i, j] / denom_safe, T(-10.0), T(10.0))
        end
    end

    @inbounds for j = 1:m, i = 1:n
        ratio = sqrt(max(wa.sum_z[i, j], T(1e-15)) / max(wa.dbeta_denom[i, j], T(1e-15)))
        scale[i, j] *= clamp(ratio, T(0.1), T(10.0))
    end

    if upd_shape
        @inbounds for j = 1:m, i = 1:n
            sh = shape[i, j]
            upd = shapelrate * (one(T) - (sh / SpecialFunctions.digamma(one(T) + one(T) / sh)) * wa.drho_numer[i, j] / wa.sum_z[i, j])
            shape[i, j] = clamp(sh + upd, minrho, maxrho)
        end
    end

    return LL_iter
end

@inline function _advance_amica_lrate!(lrate::AmicaLearningRate{T}, cap::T) where {T<:Real}
    if lrate.step_halved
        lrate.lrate = min(cap, lrate.lrate)
        lrate.step_halved = false
    else
        lrate.lrate = min(cap, lrate.lrate + min(one(T) / lrate.newt_ramp, lrate.lrate))
    end
    return nothing
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
    lrate::AmicaLearningRate{T};
    B_buf::Union{Nothing,AbstractMatrix{T}} = nothing,
    A_step::Union{Nothing,AbstractMatrix{T}} = nothing,
) where {T<:Real}
    n = size(A, 1)
    step_mat = isnothing(A_step) ? similar(A) : A_step

    if do_newton && !no_newton && iter >= newt_start_iter
        iter == newt_start_iter && (lrate.numdecs = 0)

        B = isnothing(B_buf) ? similar(dA) : B_buf
        posdef = true
        @inbounds for k = 1:n, i = 1:n
            if i == k
                nl_i = newton_lambda[i]
                if isfinite(nl_i) && abs(nl_i) > eps(T)
                    B[i, k] = dA[i, k] / nl_i
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

        if posdef
            _advance_amica_lrate!(lrate, lrate.newtrate)
            mul!(step_mat, A, B)
            @. A -= lrate.lrate * step_mat
            return no_newton
        end
        no_newton = true
    end

    _advance_amica_lrate!(lrate, lrate.lrate0)
    mul!(step_mat, A, dA)
    @. A -= lrate.lrate * step_mat
    return no_newton
end

function reparameterize_amica!(A::AbstractMatrix{T}, location::AbstractMatrix{T}, scale::AbstractMatrix{T}) where {T<:Real}
    n_rows, n_cols = size(A)
    m = size(location, 2)
    @inbounds for col = 1:n_cols
        ssq = zero(T)
        @simd for row = 1:n_rows
            ssq += A[row, col]^2
        end
        tau = sqrt(ssq)
        if tau > zero(T)
            inv_tau = one(T) / tau
            @simd for row = 1:n_rows
                A[row, col] *= inv_tau
            end
            for j = 1:m
                location[col, j] *= tau
                scale[col, j] *= inv_tau
            end
        end
    end
    return nothing
end

function calculate_amica_lrate!(
    LL::AbstractVector,
    iter::Int,
    newt_start_iter::Int,
    do_newton::Bool,
    lrate::AmicaLearningRate{T},
) where {T<:Real}
    if LL[iter] < LL[iter-1]
        lrate.lrate <= lrate.min_lrate && return true
        lrate.lrate *= lrate.lratefact
        lrate.shapelrate *= lrate.shapelratefact
        lrate.step_halved = true
        lrate.numdecs += 1
        if lrate.numdecs >= lrate.maxdecs
            lrate.lrate0 *= lrate.lratefact
            iter > newt_start_iter && (lrate.shapelrate0 *= lrate.shapelratefact)
            (do_newton && iter > newt_start_iter) && (lrate.newtrate *= lrate.lratefact)
            lrate.numdecs = 0
        end
    end
    return false
end

"""
    amica_ica(dat_ica, layout, filename; n_components, params = IcaPrms(), rng = Random.default_rng()) -> InfoIca

Run Adaptive Mixture Independent Component Analysis (AMICA), modeling source densities with
generalized Gaussian mixtures, second-order Newton updates, and log-likelihood rollback protection.
"""
function amica_ica(
    dat_ica::Matrix{Float64},
    layout::Layout,
    filename::String;
    n_components::Int,
    params::IcaPrms = IcaPrms(),
    rng::AbstractRNG = Random.default_rng(),
)
    prep = _ica_preprocess(dat_ica, n_components, params)
    n_samples = prep.n_samples
    dat_ica_t = Matrix{Float64}(prep.dat_sphered')
    LLdetS = logabsdet(prep.sphere)[1]

    m = prep.prms.amica_m > 0 ? prep.prms.amica_m : 3
    maxiter = prep.prms.amica_max_iter > 0 ? prep.prms.amica_max_iter : 2000
    lrate_val = prep.prms.amica_lrate > 0 ? prep.prms.amica_lrate : 0.1
    lrate = AmicaLearningRate(Float64; lrate0 = lrate_val)

    cpu_block = min(nextpow(2, max(Int(floor(sqrt(n_samples / 3.0))), 1)), 512)
    block_size = min(n_samples, max(cpu_block, 32))

    Wtmp = rand(rng, Float64, n_components, n_components)
    A = 0.01 .* (0.5 .- Wtmp)
    for i = 1:n_components
        A[i, i] = 1.0
        col = view(A, :, i)
        col ./= norm(col)
    end

    proportions = fill(1.0 / m, n_components, m)
    shape = fill(1.5, n_components, m)

    location = zeros(Float64, n_components, m)
    for j = 1:m
        location[:, j] .= Float64(j - 1 - (m - 1) / 2)
    end
    location .+= 0.05 .* (1.0 .- 2.0 .* rand(rng, Float64, n_components, m))
    scale = 1.0 .+ 0.1 .* (0.5 .- rand(rng, Float64, n_components, m))

    if prep.gpu_active
        A_final = _amica_optimize_gpu!(dat_ica_t, A, LLdetS, scale, location, shape, proportions, m, prep.prms, lrate)
        A .= A_final
    else
        n_workers = max(1, min(Threads.nthreads(), cld(n_samples, block_size)))
        wa_threads = [AmicaWorkArrays(Float64, n_components, m, block_size) for _ = 1:n_workers]
        qconsts = Matrix{Float64}(undef, n_components, m)
        LL = fill(-Inf, maxiter)

        do_newton = true
        no_newton = false
        newt_start_iter = 50
        update_shape = true

        best_ll = -Inf
        cur_snap = _make_amica_snap(A, scale, location, shape, proportions)
        best_snap = _make_amica_snap(A, scale, location, shape, proportions)
        converged = false
        aborted_nonfinite = false
        consecutive_conv = 0
        tol = prep.prms.default_stop > 0 ? prep.prms.default_stop : 1e-7
        last_iter = 0

        for iter = 1:maxiter
            last_iter = iter
            _save_amica_snap!(cur_snap, A, scale, location, shape, proportions)

            LL_iter = update_amica_parameters!(
                wa_threads,
                qconsts,
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

            if !isfinite(LL_iter)
                @minimal_warning "AMICA encountered non-finite log-likelihood ($LL_iter) at iteration $iter. Rolling back to best model."
                aborted_nonfinite = true
                _restore_amica_snap!(A, scale, location, shape, proportions, best_snap)
                break
            end

            if LL_iter > best_ll
                best_ll = LL_iter
                cur_snap, best_snap = best_snap, cur_snap
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
                lrate;
                B_buf = wa.B_buf,
                A_step = wa.A_step,
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

        if !converged && !aborted_nonfinite
            @minimal_warning "AMICA did not converge within $maxiter iterations (best LL = $(round(best_ll, digits=7))). Returning best model encountered."
        end

        if best_ll > -Inf && (!converged || (last_iter > 0 && LL[last_iter] < best_ll))
            _restore_amica_snap!(A, scale, location, shape, proportions, best_snap)
        end
    end

    W_final = inv(A)
    U = dat_ica_t * transpose(W_final)
    is_sub_gaussian = falses(n_components)
    inv_samples = 1.0 / n_samples
    @inbounds for i = 1:n_components
        s1 = 0.0
        @simd for s = 1:n_samples
            s1 += U[s, i]
        end
        m_i = s1 * inv_samples
        s2 = 0.0
        s4 = 0.0
        @simd for s = 1:n_samples
            d = U[s, i] - m_i
            d2 = d * d
            s2 += d2
            s4 += d2 * d2
        end
        v_i = s2 * inv_samples
        if v_i > 1e-12
            is_sub_gaussian[i] = ((s4 * inv_samples) / (v_i * v_i) - 3.0) < -0.05
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

    best_ll = -Inf
    cur_snap = _make_amica_snap(A_cpu, scale_cpu, location_cpu, shape_cpu, proportions_cpu)
    best_snap = _make_amica_snap(A_cpu, scale_cpu, location_cpu, shape_cpu, proportions_cpu)
    converged = false
    aborted_nonfinite = false
    consecutive_conv = 0
    tol = params.default_stop > 0 ? params.default_stop : 1e-7
    last_iter = 0

    for iter = 1:maxiter
        last_iter = iter
        _save_amica_snap!(cur_snap, A_cpu, scale_cpu, location_cpu, shape_cpu, proportions_cpu)

        W = gpu_array(Float32.(inv(A_cpu)))
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

        @. fp_z_y = fp * z * y
        dbeta_denom = dropdims(sum(fp_z_y, dims = 1), dims = 1)

        @. dlambda_tmp = z * (fp * y - 1.0f0)^2
        dlambda_numer = dropdims(sum(dlambda_tmp, dims = 1), dims = 1)

        ldet = -logabsdet(A_cpu)[1]
        Lt_sum = sum(Float64, Array(vec(sum(logexp, dims = 1))))
        LL_iter = (Float64(ldet + LLdetS) * Float64(N) + Lt_sum) / Float64(N * n)

        if !isfinite(LL_iter)
            @minimal_warning "AMICA GPU encountered non-finite log-likelihood ($LL_iter) at iteration $iter. Rolling back to best model."
            aborted_nonfinite = true
            _restore_amica_snap!(A_cpu, scale_cpu, location_cpu, shape_cpu, proportions_cpu, best_snap)
            break
        end

        if LL_iter > best_ll
            best_ll = LL_iter
            cur_snap, best_snap = best_snap, cur_snap
        end
        LL[iter] = LL_iter

        @. dA = -g_times_sources / Float32(N)
        dA[diagind(dA)] .+= 1.0f0

        sum_z_c = Float64.(Array(sum_z))
        dmu_numer_c = Float64.(Array(dmu_numer))
        dmu_denom_c = Float64.(Array(dmu_denom))
        dbeta_denom_c = Float64.(Array(dbeta_denom))
        drho_numer_c = Float64.(Array(drho_numer))
        kp_c = Float64.(Array(kp))
        dlambda_numer_c = Float64.(Array(dlambda_numer))
        newton_sigma2_c = Float64.(Array(newton_sigma2))

        if m > 1
            raw_prop = max.(sum_z_c ./ Float64(N), 1.0e-4)
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

        no_newton = update_amica_mixing!(
            A_cpu,
            Float64.(Array(dA)),
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

    if !converged && !aborted_nonfinite
        @minimal_warning "AMICA GPU did not converge within $maxiter iterations (best LL = $(round(best_ll, digits=7))). Returning best model encountered."
    end

    if best_ll > -Inf && (!converged || (last_iter > 0 && LL[last_iter] < best_ll))
        _restore_amica_snap!(A_cpu, scale_cpu, location_cpu, shape_cpu, proportions_cpu, best_snap)
    end

    return A_cpu
end
