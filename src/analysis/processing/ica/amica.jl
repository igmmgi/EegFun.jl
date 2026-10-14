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

    # Preallocated block arrays
    source_signals::Array{T,2}
    y::Array{T,3}
    y_rho::Array{T,3}
    Q::Array{T,3}
    z::Array{T,3}
    fp::Array{T,3}
    scratch3::Array{T,3}
    scratch_N_n::Array{T,2}
    scratch_N_n2::Array{T,2}
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

    Threads.@threads for tid = 1:n_threads
        wa = wa_threads[tid]
        for block = tid:n_threads:num_blocks
            lower = (block - 1) * block_size + 1
            upper = min(N, block * block_size)
            n_samples = upper - lower + 1
            rng = lower:upper

            source_signals = view(wa.source_signals, 1:n_samples, :)
            y = view(wa.y, (1:n_samples), :, :)
            y_rho = view(wa.y_rho, (1:n_samples), :, :)
            Q = view(wa.Q, (1:n_samples), :, :)
            z = view(wa.z, (1:n_samples), :, :)
            fp = view(wa.fp, (1:n_samples), :, :)
            scratch3 = view(wa.scratch3, (1:n_samples), :, :)
            Qmax = view(wa.scratch_N_n, 1:n_samples, :)
            logexp = view(wa.scratch_N_n2, 1:n_samples, :)
            zsum = view(wa.scratch_N_n, 1:n_samples, :)
            g_block_sum = view(wa.scratch_N_n, 1:n_samples, :)

            mul!(source_signals, view(data, rng, :), transpose(W))

            for j = 1:m
                for i = 1:n
                    @fastmath @inbounds @simd for s = 1:n_samples
                        y[s, i, j] = scale[i, j] * (source_signals[s, i] - location[i, j])
                    end
                end
            end

            for j = 1:m
                for i = 1:n
                    @fastmath @inbounds @simd for s = 1:n_samples
                        val = y[s, i, j]
                        val = ifelse(val == zero(T), T(1e-16), val)
                        y_rho[s, i, j] = exp((shape[i, j] - T(1.0)) * log(abs(val)))
                    end
                end
            end

            for j = 1:m
                for i = 1:n
                    qconst = -log(T(2)) - SpecialFunctions.loggamma(T(1) + T(1) / shape[i, j]) + log(proportions[i, j]) + log(scale[i, j])
                    @fastmath @inbounds @simd for s = 1:n_samples
                        Q[s, i, j] = qconst - (y_rho[s, i, j] * abs(y[s, i, j]))
                    end
                end
            end

            maximum!(view(Qmax, :, :, (1:1)), Q)
            @fastmath @inbounds @simd for I in eachindex(Q)
                scratch3[I] = exp(Q[I] - Qmax[I[1], I[2]])
            end
            sum!(view(logexp, :, :, (1:1)), scratch3)
            @fastmath @inbounds @simd for I in eachindex(logexp)
                logexp[I] = log(logexp[I]) + Qmax[I]
            end

            @fastmath @inbounds @simd for I in eachindex(Q)
                z[I] = exp(Q[I] - logexp[I[1], I[2]]) + T(1e-15)
            end

            for i = 1:n
                @fastmath @inbounds @simd for s = 1:n_samples
                    wa.Lt_accum[lower+s-1] += logexp[s, i]
                end
            end

            sum!(view(zsum, :, :, (1:1)), z)
            @fastmath @inbounds @simd for I in eachindex(z)
                z[I] = z[I] / zsum[I[1], I[2]]
            end

            for j = 1:m
                for i = 1:n
                    @fastmath @inbounds @simd for s = 1:n_samples
                        fp[s, i, j] = y_rho[s, i, j] * sign(y[s, i, j]) * shape[i, j]
                    end
                end
            end

            for j = 1:m
                for i = 1:n
                    @fastmath @inbounds @simd for s = 1:n_samples
                        sc = y_rho[s, i, j] * abs(y[s, i, j])
                        sc = ifelse(sc >= T(1.0e-16), z[s, i, j] * log(sc) * sc, T(0.0))
                        wa.drho_numer[i, j] += sc
                    end
                end
            end

            for j = 1:m
                for i = 1:n
                    @fastmath @inbounds @simd for s = 1:n_samples
                        scratch3[s, i, j] = scale[i, j] * z[s, i, j] * fp[s, i, j]
                    end
                end
            end
            sum!(view(g_block_sum, :, :, (1:1)), scratch3)
            mul!(wa.g_times_sources, transpose(g_block_sum), source_signals, one(T), one(T))

            if newton_active
                for i = 1:n
                    @fastmath @inbounds @simd for s = 1:n_samples
                        wa.newton_sigma2[i] += (source_signals[s, i]^2) / N
                    end
                end
            end

            for j = 1:m
                for i = 1:n
                    @fastmath @inbounds @simd for s = 1:n_samples
                        z_val = z[s, i, j]
                        fp_val = fp[s, i, j]
                        y_val = y[s, i, j]
                        wa.sum_z[i, j] += z_val
                        wa.dmu_numer[i, j] += fp_val * z_val
                        wa.kp[i, j] += fp_val * z_val * fp_val
                        if shape[i, j] <= T(2)
                            y_safe = ifelse(y_val == zero(T), T(1e-16), y_val)
                            wa.dmu_denom[i, j] += (z_val * fp_val / y_safe) * scale[i, j]
                            wa.dbeta_denom[i, j] += fp_val * z_val * y_val
                        else
                            wa.dmu_denom[i, j] += (z_val * fp_val * fp_val) * scale[i, j]
                        end
                        wa.dlambda_numer[i, j] += z_val * (fp_val * y_val - T(1.0))^2
                    end
                end
            end
        end
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
        proportions .= ifelse.(wa.sum_z .>= T(0), wa.sum_z ./ N, T(1) / N)
    end

    if newton_active
        dkap = (wa.kp ./ (proportions .* N)) .* scale .^ 2
        wa.newton_kappa .= dropdims(sum(proportions .* dkap, dims = 2), dims = 2)
        wa.newton_lambda .= dropdims(sum(proportions .* (wa.dlambda_numer ./ wa.sum_z .+ dkap .* location .^ 2), dims = 2), dims = 2)
    end

    if m > 1
        location .+= wa.dmu_numer ./ wa.dmu_denom
    end

    scale .*= sqrt.(wa.sum_z ./ wa.dbeta_denom)

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

function amica_ica(dat_ica::Matrix{Float64}, layout::Layout, filename::String; n_components::Int, params::IcaPrms = IcaPrms())

    n_channels, n_samples = size(dat_ica)

    # PCA / Sphering identically to infomax
    original_mean = vec(mean(dat_ica, dims = 2))
    dat_ica .-= original_mean

    data_scale = sqrt(norm((dat_ica * transpose(dat_ica)) / n_samples))
    dat_ica ./= data_scale

    # Apply SVD
    F = svd(dat_ica)
    pca_components = F.U[:, 1:n_components]
    eigenvalues = (F.S[1:n_components] .^ 2) ./ max(n_samples - 1, 1)
    sphere = diagm(1.0 ./ sqrt.(eigenvalues))

    transform_matrix = sphere * transpose(pca_components)
    dat_ica_sphered = Matrix{Float64}(undef, n_components, n_samples)
    mul!(dat_ica_sphered, transform_matrix, dat_ica)
    dat_ica = dat_ica_sphered

    # AMICA native math expects data as (N_samples, n_components)
    dat_ica_t = Matrix{Float64}(dat_ica')

    LLdetS = logabsdet(sphere)[1]

    # Amica parameters
    m = 3 # number of mixtures, could be parameterized
    maxiter = params.max_iter > 0 ? params.max_iter : 2000
    lrate_val = params.l_rate > 0 ? params.l_rate : 0.1
    lrate = AmicaLearningRate(Float64; lrate0 = lrate_val)

    # Dynamic block size based on N
    cpu_block = min(nextpow(2, max(Int(floor(sqrt(n_samples / 3.0))), 1)), 512)
    block_size = min(n_samples, max(cpu_block, 32))

    rng = Xoshiro(42)

    # Initialize A to match Fortran/Amica.jl: small random ±0.005, diagonal = 1.0, then normalize columns
    Wtmp = rand(rng, Float64, n_components, n_components)
    A = 0.01 .* (0.5 .- Wtmp)
    for i = 1:n_components
        A[i, i] = 1.0
        A[:, i] = A[:, i] / norm(A[:, i])
    end

    proportions = fill(1.0 / m, n_components, m)
    shape = fill(1.5, n_components, m)

    # Initialize location to match Fortran: centered around 0 (e.g., -1, 0, 1 for m=3) + small noise ±0.05
    location = zeros(Float64, n_components, m)
    for j = 1:m
        location[:, j] .= Float64(j - 1 - (m - 1) / 2)
    end
    location .+= 0.05 .* (1.0 .- 2.0 .* rand(rng, Float64, n_components, m))

    # Initialize scale to match Fortran: values in range [0.95, 1.05]
    scale = ones(Float64, n_components, m) .+ 0.1 .* (0.5 .- rand(rng, Float64, n_components, m))

    n_threads = Threads.nthreads()
    wa_threads = [AmicaWorkArrays(Float64, n_samples, n_components, m, block_size) for _ = 1:n_threads]
    LL = fill(-Inf, maxiter)

    do_newton = true
    no_newton = false
    newt_start_iter = 50
    update_shape = true

    # Detect GPU hardware backend
    gpu_active = false
    if params.use_gpu
        if is_gpu_available()
            gpu_active = true
            @info "[GPU ACTIVATED] Running AMICA ICA on $(gpu_device_name())..."
        else
            @minimal_warning "Requested GPU acceleration (use_gpu=true), but no functional GPU package (CUDA.jl, AMDGPU.jl, Metal.jl) has been loaded. Please run 'using CUDA' or 'using AMDGPU' before calling run_ica. Falling back to CPU."
        end
    end

    # Main loop
    if gpu_active
        LL_gpu = _amica_optimize_gpu!(dat_ica_t, A, LLdetS, scale, location, shape, proportions, m, params, lrate)
        LL .= LL_gpu
    else
        for iter = 1:maxiter
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
                if isnan(LL[iter])
                    @warn("Got NaN! Exiting ...")
                    break
                end
                if calculate_amica_lrate!(LL, iter, newt_start_iter, do_newton, lrate)
                    break
                end

                # Convergence check (from Amica)
                if abs(LL[iter] - LL[iter-1]) < 1e-7
                    break
                end
            end

            if iter == 1 || iter % 10 == 0
                @info "AMICA iter $iter, LL = $(LL[iter])"
            end
        end
    end

    # Final calculations
    W = inv(A)
    weights = W * sphere * pca_components'
    mixing = pinv(weights)

    # Calculate total variance explained and order
    # Since dat_ica was overwritten with the sphered data, its row variances are uniform.
    # The true variance of source i is proportional to the squared norm of the i-th row of W.
    source_vars = vec(sum(abs2, W, dims = 2))
    meanvar = vec(sum(abs2, mixing, dims = 1)) .* source_vars
    meanvar_normalized = meanvar ./ sum(meanvar)
    order = sortperm(meanvar_normalized, rev = true)

    return InfoIca(
        filename,
        weights[order, :],
        mixing[:, order],
        sphere,
        meanvar_normalized[order],
        data_scale,
        original_mean,
        [Symbol("IC$i") for i = 1:size(weights, 1)],
        Dict{Int,Matrix{Float64}}(),
        layout,
        falses(size(weights, 1)),
    )
end



# === GPU versions of Amica ===

# GPU versions of Amica
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

    # Preallocate GPU full-batch arrays
    y = gpu_array(zeros(Float32, N, n, m))
    y_rho = gpu_array(zeros(Float32, N, n, m))
    Q = gpu_array(zeros(Float32, N, n, m))
    z = gpu_array(zeros(Float32, N, n, m))
    fp = gpu_array(zeros(Float32, N, n, m))

    qconst_cpu = zeros(Float32, 1, n, m)
    source_signals = gpu_array(zeros(Float32, N, n))

    scale = gpu_array(Float32.(scale_cpu))
    location = gpu_array(Float32.(location_cpu))
    shape = gpu_array(Float32.(shape_cpu))
    proportions = gpu_array(Float32.(proportions_cpu))

    # We will reshape these to 1 x n x m during broadcasting
    scale_3d = reshape(scale, 1, n, m)
    location_3d = reshape(location, 1, n, m)
    shape_3d = reshape(shape, 1, n, m)
    proportions_3d = reshape(proportions, 1, n, m)

    dA = gpu_array(zeros(Float32, n, n))

    LL = fill(-Inf, params.max_iter > 0 ? params.max_iter : 2000)
    do_newton = true
    no_newton = false
    newt_start_iter = 50
    update_shape = true
    maxiter = params.max_iter > 0 ? params.max_iter : 2000

    for iter = 1:maxiter
        W_c = inv(A_cpu)
        W = gpu_array(Float32.(W_c))

        mul!(source_signals, dat_gpu, transpose(W))

        # We need qconst
        for j = 1:m, i = 1:n
            qconst_cpu[1, i, j] =
                -log(2.0f0) - SpecialFunctions.loggamma(1.0f0 + 1.0f0 / Float32(shape_cpu[i, j])) +
                log(Float32(proportions_cpu[i, j])) +
                log(Float32(scale_cpu[i, j]))
        end
        qconst = gpu_array(qconst_cpu)

        ss_3d = reshape(source_signals, N, n, 1)

        # y = scale * (source_signals - location)
        @. y = scale_3d * (ss_3d - location_3d)

        # y_rho = exp((shape - 1.0) * log(abs(y)))
        @. y_rho = exp((shape_3d - 1.0f0) * log(max(abs(y), 1.0f-16)))

        # Q = qconst - (y_rho * abs(y))
        @. Q = qconst - (y_rho * abs(y))

        Qmax = maximum(Q, dims = 3)

        scratch3 = exp.(Q .- Qmax)
        logexp = log.(sum(scratch3, dims = 3)) .+ Qmax

        @. z = exp(Q - logexp) + 1.0f-15

        zsum = sum(z, dims = 3)
        @. z = z / zsum

        @. fp = y_rho * sign(y) * shape_3d

        sc = @. y_rho * abs(y)
        sc = @. ifelse(sc >= 1.0f-16, z * log(sc) * sc, 0.0f0)
        drho_numer = dropdims(sum(sc, dims = 1), dims = 1) # n x m

        g_block_sum = sum(scale_3d .* z .* fp, dims = 3) # N x n x 1
        g_times_sources = transpose(dropdims(g_block_sum, dims = 3)) * source_signals # n x n

        newton_sigma2 = dropdims(sum(source_signals .^ 2, dims = 1), dims = 1) ./ Float32(N) # n

        sum_z = dropdims(sum(z, dims = 1), dims = 1) # n x m
        dmu_numer = dropdims(sum(fp .* z, dims = 1), dims = 1) # n x m
        kp = dropdims(sum(fp .* z .* fp, dims = 1), dims = 1) # n x m

        mask = shape_3d .<= 2.0f0
        y_safe = @. ifelse(y == 0.0f0, 1.0f-16, y)

        dmu_denom_1 = @. (z * fp / y_safe) * scale_3d
        dmu_denom_2 = @. (z * fp * fp) * scale_3d
        dmu_denom = dropdims(sum(ifelse.(mask, dmu_denom_1, dmu_denom_2), dims = 1), dims = 1)

        dbeta_denom = dropdims(sum(ifelse.(mask, fp .* z .* y, 0.0f0), dims = 1), dims = 1)
        dlambda_numer = dropdims(sum(z .* (fp .* y .- 1.0f0) .^ 2, dims = 1), dims = 1)

        ldet = -logabsdet(A_cpu)[1]
        Lt_sum = sum(logexp) # scalar
        Lt = Float32(ldet + LLdetS) * Float32(N) + Lt_sum
        LL_iter = Lt / (N * n)
        LL[iter] = Float64(LL_iter)

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
            @. proportions_cpu = ifelse(sum_z_c >= 0.0, sum_z_c / N, 1.0 / N)
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
            @. location_cpu += dmu_numer_c / dmu_denom_c
        end

        @. scale_cpu *= sqrt(sum_z_c / dbeta_denom_c)

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
            if isnan(LL[iter])
                @warn("Got NaN! Exiting ...")
                break
            end
            if calculate_amica_lrate!(LL, iter, newt_start_iter, do_newton, lrate)
                break
            end
            if abs(LL[iter] - LL[iter-1]) < 1e-7
                break
            end
        end

        if iter == 1 || iter % 10 == 0
            @info "AMICA iter $iter, LL = $(LL[iter])"
        end
    end

    return LL
end

"""AMICA ICA implementation natively integrated for EegFun.jl"""
