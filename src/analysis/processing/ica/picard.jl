const _PICARD_THREAD_ELEMS = 16_384

# Avoid oversubscribing hardware threads when OpenBLAS worker threads are also active in Picard
@inline _picard_num_tasks(T_len::Int) = clamp(Sys.CPU_THREADS - BLAS.get_num_threads(), 1, min(Threads.nthreads(), cld(T_len, 256)))

@inline function _log_cosh(y::T) where {T<:AbstractFloat}
    ay = abs(y)
    return ay > T(30.0) ? ay : ay + log1p(exp(T(-2.0) * ay))
end

function _picard_y_loss_threaded(
    Y::Matrix{T},
    signs::Vector{T},
    extended::Bool,
    thread_sums::Vector{Float64},
    DY::Matrix{T},
    alpha::T,
) where {T<:AbstractFloat}
    N_ch, T_len = size(Y)
    fused = alpha != zero(T)
    nth = min(length(thread_sums), _picard_num_tasks(T_len))
    chunk_size = cld(T_len, nth)
    fill!(thread_sums, 0.0)
    Threads.@threads :static for tid = 1:nth
        j_start = (tid - 1) * chunk_size + 1
        j_end = min(tid * chunk_size, T_len)
        local_sum = 0.0
        if extended
            @inbounds for j = j_start:j_end, i = 1:N_ch
                y = fused ? Y[i, j] + alpha * DY[i, j] : Y[i, j]
                local_sum += Float64(signs[i] * _log_cosh(y) + T(0.5) * y^2)
            end
        else
            @inbounds for j = j_start:j_end, i = 1:N_ch
                y = fused ? Y[i, j] + alpha * DY[i, j] : Y[i, j]
                local_sum += Float64(signs[i] * _log_cosh(y))
            end
        end
        thread_sums[tid] = local_sum
    end
    return T(sum(@view(thread_sums[1:nth])) / Float64(T_len))
end

function _picard_y_loss_serial(Y::Matrix{T}, signs::Vector{T}, extended::Bool, DY::Matrix{T}, alpha::T) where {T<:AbstractFloat}
    N_ch, T_len = size(Y)
    fused = alpha != zero(T)
    total = 0.0
    if extended
        @inbounds for j = 1:T_len, i = 1:N_ch
            y = fused ? Y[i, j] + alpha * DY[i, j] : Y[i, j]
            total += Float64(signs[i] * _log_cosh(y) + T(0.5) * y^2)
        end
    else
        @inbounds for j = 1:T_len, i = 1:N_ch
            y = fused ? Y[i, j] + alpha * DY[i, j] : Y[i, j]
            total += Float64(signs[i] * _log_cosh(y))
        end
    end
    return T(total / Float64(T_len))
end

@inline function _picard_y_loss(
    Y::Matrix{T},
    signs::Vector{T},
    extended::Bool,
    thread_sums::Vector{Float64};
    DY::Matrix{T} = Y,
    alpha::T = zero(T),
) where {T<:AbstractFloat}
    if length(Y) >= _PICARD_THREAD_ELEMS && Threads.nthreads() > 1
        return _picard_y_loss_threaded(Y, signs, extended, thread_sums, DY, alpha)
    end
    return _picard_y_loss_serial(Y, signs, extended, DY, alpha)
end

function _picard_loss(Y::Matrix{T}, W::Matrix{T}, signs::Vector{T}, extended::Bool, thread_sums::Vector{Float64}) where {T<:AbstractFloat}
    logdetW, _ = logabsdet(W)
    return T(-logdetW + _picard_y_loss(Y, signs, extended, thread_sums))
end

function _picard_solve_hessian!(Z::AbstractMatrix{T}, h::AbstractMatrix{T}, G::AbstractMatrix{T}) where {T<:AbstractFloat}
    N = size(Z, 1)
    ht = transpose(h)
    Gt = transpose(G)
    @. Z = (ht * G - Gt) / (h * ht - one(T))
    @inbounds for i = 1:N
        Z[i, i] = G[i, i] / (h[i, i] + one(T))
    end
    return nothing
end

function _picard_regularize_hessian!(h::Matrix{T}, lambda_min::T) where {T<:AbstractFloat}
    N = size(h, 1)
    @inbounds for j = 1:N, i = (j+1):N
        hij = h[i, j]
        hji = h[j, i]
        eig = T(0.5) * (hij + hji - sqrt((hij - hji)^2 + T(4.0)))
        if eig < lambda_min
            update = lambda_min - eig
            h[i, j] += update
            h[j, i] += update
        end
    end
    return nothing
end

mutable struct LBFGSBuffer{T<:AbstractFloat,MatT<:AbstractMatrix{T}}
    s_buf::Vector{MatT}
    y_buf::Vector{MatT}
    r_buf::Vector{T}
    a_buf::Vector{T}
    m::Int
    len::Int
    head::Int
end

function LBFGSBuffer(T::Type{<:AbstractFloat}, N::Int, m::Int)
    return LBFGSBuffer{T,Matrix{T}}([zeros(T, N, N) for _ = 1:m], [zeros(T, N, N) for _ = 1:m], zeros(T, m), zeros(T, m), m, 0, 1)
end

function LBFGSBuffer(sample_mat::MatT, m::Int) where {T<:AbstractFloat,MatT<:AbstractMatrix{T}}
    N = size(sample_mat, 1)
    s_buf = [similar(sample_mat, N, N) for _ = 1:m]
    y_buf = [similar(sample_mat, N, N) for _ = 1:m]
    for k = 1:m
        fill!(s_buf[k], zero(T))
        fill!(y_buf[k], zero(T))
    end
    return LBFGSBuffer{T,MatT}(s_buf, y_buf, zeros(T, m), zeros(T, m), m, 0, 1)
end

@inline function _lbfgs_reset!(buf::LBFGSBuffer)
    buf.len = 0
    buf.head = 1
    return nothing
end

function _lbfgs_push!(buf::LBFGSBuffer{T,MatT}, direction::MatT, y_diff::MatT, r_val::T) where {T,MatT<:AbstractMatrix{T}}
    idx = buf.head
    copyto!(buf.s_buf[idx], direction)
    copyto!(buf.y_buf[idx], y_diff)
    buf.r_buf[idx] = r_val
    buf.head = (idx % buf.m) + 1
    buf.len = min(buf.len + 1, buf.m)
    return nothing
end

@inline function _lbfgs_get(buf::LBFGSBuffer, k::Int)
    idx = mod1(buf.head - k, buf.m)
    return buf.s_buf[idx], buf.y_buf[idx], buf.r_buf[idx]
end

function _picard_l_bfgs_direction!(
    Z::MatT,
    q::MatT,
    G::MatT,
    h::MatT,
    buf::LBFGSBuffer{T,MatT},
) where {T<:AbstractFloat,MatT<:AbstractMatrix{T}}
    copyto!(q, G)
    m_len = buf.len

    for k = 1:m_len
        s, y, r = _lbfgs_get(buf, k)
        a = r * dot(s, q)
        buf.a_buf[k] = a
        @. q -= a * y
    end

    _picard_solve_hessian!(Z, h, q)

    for k = m_len:-1:1
        s, y, r = _lbfgs_get(buf, k)
        a = buf.a_buf[k]
        beta = r * dot(y, Z)
        @. Z += (a - beta) * s
    end

    @. Z = -Z
    return nothing
end

function _picard_line_search(
    Y::Matrix{T},
    W::Matrix{T},
    Y_new::Matrix{T},
    W_new::Matrix{T},
    DY::Matrix{T},
    DW::Matrix{T},
    transform::Matrix{T},
    direction::Matrix{T},
    signs::Vector{T},
    current_loss::T,
    logdetW::T,
    ls_tries::Int,
    extended::Bool,
    thread_sums::Vector{Float64},
) where {T<:AbstractFloat}
    N = size(Y, 1)
    alpha = one(T)

    mul!(DY, direction, Y)
    mul!(DW, direction, W)

    for _ = 1:ls_tries
        @. transform = alpha * direction
        @inbounds for i = 1:N
            transform[i, i] += one(T)
        end
        step_logdet, _ = logabsdet(transform)
        logdet_total = logdetW + T(step_logdet)

        new_loss = -logdet_total + _picard_y_loss(Y, signs, extended, thread_sums; DY = DY, alpha = alpha)
        if isfinite(new_loss) && new_loss < current_loss
            @. Y_new = Y + alpha * DY
            @. W_new = W + alpha * DW
            return true, new_loss, alpha, logdet_total
        end
        alpha /= T(2.0)
    end

    copyto!(Y_new, Y)
    copyto!(W_new, W)
    return false, current_loss, alpha, logdetW
end

function _picard_eval_nonlinear_threaded!(psiY::Matrix{T}, psidY::Matrix{T}, Y_square::Matrix{T}, Y::Matrix{T}) where {T<:AbstractFloat}
    N, T_len = size(Y)
    nth = _picard_num_tasks(T_len)
    chunk_size = cld(T_len, nth)
    Threads.@threads :static for tid = 1:nth
        j_start = (tid - 1) * chunk_size + 1
        j_end = min(tid * chunk_size, T_len)
        @inbounds for j = j_start:j_end
            @simd for i = 1:N
                yij = Y[i, j]
                @fastmath ty = tanh(yij)
                psiY[i, j] = ty
                psidY[i, j] = one(T) - ty * ty
                Y_square[i, j] = yij * yij
            end
        end
    end
    return nothing
end

function _picard_scale_psid_threaded!(psidY::Matrix{T}, signs::Vector{T}) where {T<:AbstractFloat}
    N, T_len = size(psidY)
    nth = _picard_num_tasks(T_len)
    chunk_size = cld(T_len, nth)
    Threads.@threads :static for tid = 1:nth
        j_start = (tid - 1) * chunk_size + 1
        j_end = min(tid * chunk_size, T_len)
        @inbounds for j = j_start:j_end
            @simd for i = 1:N
                psidY[i, j] = psidY[i, j] * signs[i] + one(T)
            end
        end
    end
    return nothing
end

function _picard_optimize_cpu!(W::Matrix{T}, dat_ica::Matrix{T}, params::IcaPrms, extended::Bool) where {T<:AbstractFloat}
    N, T_len = size(dat_ica)
    m = params.picard_m
    tol = T(params.default_stop)
    lambda_min = T(params.picard_lambda_min)
    ls_tries = params.picard_ls_tries
    max_iter = params.max_iter
    use_threads = (N * T_len >= _PICARD_THREAD_ELEMS) && (Threads.nthreads() > 1)
    inv_T = one(T) / T(T_len)

    Y = copy(dat_ica)
    Y_new = similar(Y)
    W_new = similar(W)
    DY = zeros(T, N, T_len)
    DW = zeros(T, N, N)
    transform = zeros(T, N, N)
    y_diff = zeros(T, N, N)

    psiY = zeros(T, N, T_len)
    psidY = zeros(T, N, T_len)
    Y_square = zeros(T, N, T_len)

    G = zeros(T, N, N)
    G_old = zeros(T, N, N)
    h = zeros(T, N, N)
    direction = zeros(T, N, N)
    q = zeros(T, N, N)

    signs = ones(T, N)
    old_signs = ones(T, N)
    lbfgs = LBFGSBuffer(T, N, m)
    thread_sums = zeros(Float64, Threads.nthreads())

    logdetW = T(logabsdet(W)[1])
    current_loss = _picard_loss(Y, W, signs, extended, thread_sums)
    gradient_norm = zero(T)

    C = extended ? zeros(T, N, N) : zeros(T, 0, 0)
    C_orig = extended ? zeros(T, N, N) : zeros(T, 0, 0)
    C_tmp = extended ? zeros(T, N, N) : zeros(T, 0, 0)
    mean_psid = extended ? zeros(T, N, 1) : zeros(T, 0, 0)

    if extended
        mul!(C_orig, dat_ica, transpose(dat_ica))
        C_orig .*= inv_T
        mul!(C_tmp, W, C_orig)
        mul!(C, C_tmp, transpose(W))
    end

    for n = 1:max_iter
        if use_threads
            _picard_eval_nonlinear_threaded!(psiY, psidY, Y_square, Y)
        else
            @fastmath @inbounds @simd for k in eachindex(Y)
                yij = Y[k]
                ty = tanh(yij)
                psiY[k] = ty
                psidY[k] = one(T) - ty * ty
                Y_square[k] = yij * yij
            end
        end

        mul!(G, psiY, transpose(Y))
        G .*= inv_T

        if extended
            sum!(mean_psid, psidY)
            sign_change = (n == 1)
            @inbounds for i = 1:N
                k_val = (mean_psid[i, 1] * inv_T) * C[i, i] - G[i, i]
                s_val = k_val >= zero(T) ? one(T) : -one(T)
                s_val != old_signs[i] && (sign_change = true)
                signs[i] = s_val
                old_signs[i] = s_val
            end

            @inbounds for j_col = 1:N, i = 1:N
                G[i, j_col] = G[i, j_col] * signs[i] + C[i, j_col]
            end

            if use_threads && length(psidY) >= 262_144
                _picard_scale_psid_threaded!(psidY, signs)
            else
                @inbounds for j = 1:T_len
                    @simd for i = 1:N
                        psidY[i, j] = psidY[i, j] * signs[i] + one(T)
                    end
                end
            end

            if sign_change
                current_loss = _picard_loss(Y, W, signs, extended, thread_sums)
                _lbfgs_reset!(lbfgs)
            end
        end

        mul!(h, psidY, transpose(Y_square))
        h .*= inv_T
        _picard_regularize_hessian!(h, lambda_min)

        @inbounds for i = 1:N
            G[i, i] -= one(T)
        end

        gradient_norm = maximum(abs, G)
        if gradient_norm < tol
            @info "Picard converged at step $n (gradient norm: $gradient_norm)"
            return true, n, signs, gradient_norm
        end

        if n > 1
            @. y_diff = G - G_old
            sy = dot(direction, y_diff)
            if sy > eps(T) * norm(direction) * norm(y_diff)
                _lbfgs_push!(lbfgs, direction, y_diff, one(T) / sy)
            else
                _lbfgs_reset!(lbfgs)
            end
        end
        copyto!(G_old, G)

        _picard_l_bfgs_direction!(direction, q, G, h, lbfgs)

        ls_ok, new_loss, alpha, logdetW = _picard_line_search(
            Y,
            W,
            Y_new,
            W_new,
            DY,
            DW,
            transform,
            direction,
            signs,
            current_loss,
            logdetW,
            ls_tries,
            extended,
            thread_sums,
        )

        if !ls_ok
            _picard_solve_hessian!(direction, h, G)
            @. direction = -direction
            _lbfgs_reset!(lbfgs)
            ls_ok, new_loss, alpha, logdetW = _picard_line_search(
                Y,
                W,
                Y_new,
                W_new,
                DY,
                DW,
                transform,
                direction,
                signs,
                current_loss,
                logdetW,
                ls_tries,
                extended,
                thread_sums,
            )

            if !ls_ok
                @info "Picard line search failed to find an improvement. Stopping early at step $n (gradient norm: $gradient_norm)."
                return false, n, signs, gradient_norm
            end
        end

        @. direction *= alpha
        copyto!(Y, Y_new)
        copyto!(W, W_new)
        current_loss = new_loss

        if extended
            mul!(C_tmp, W, C_orig)
            mul!(C, C_tmp, transpose(W))
        end

        if n == 1 || n % 10 == 0
            @info Printf.@sprintf("Picard step %d, gradient norm = %.7f, loss = %.7f", n, gradient_norm, current_loss)
        end
    end

    return false, max_iter, signs, gradient_norm
end

"""
    picard_ica(dat_ica, layout, filename; n_components, extended = false,
               params = IcaPrms(), rng = Random.default_rng()) -> InfoIca

Run Preconditioned ICA for Real Data (Picard; Ablin, Cardoso, & Gramfort, 2018), using an exact
2×2 block Hessian preconditioner with L-BFGS (`extended = true` separates mixed sub-/super-Gaussian sources).
"""
function picard_ica(
    dat_ica::Matrix{Float64},
    layout::Layout,
    filename::String;
    n_components::Int,
    extended::Bool = false,
    params::IcaPrms = IcaPrms(),
    rng::AbstractRNG = Random.default_rng(),
)
    ext_str = extended ? "Extended-" : ""
    @info "Running $(ext_str)Picard ICA (use_gpu=$(params.use_gpu)): $(size(dat_ica, 1)) channels x $(size(dat_ica, 2)) samples -> $n_components components"

    prep = _ica_preprocess(dat_ica, n_components, params)

    if prep.gpu_active
        W_init = Matrix{Float32}(I, n_components, n_components)
        converged, n_iter, signs_f32, grad_norm = _picard_optimize_gpu!(W_init, Float32.(prep.dat_sphered), prep.prms, extended)
        W_final = Float64.(W_init)
        signs = Float64.(signs_f32)
    else
        W_final = Matrix{Float64}(I, n_components, n_components)
        converged, n_iter, signs, grad_norm = _picard_optimize_cpu!(W_final, prep.dat_sphered, prep.prms, extended)
    end

    if !converged
        @minimal_warning "Picard ICA reached maximum iterations or stopped early at iteration $n_iter without meeting tolerance ($(prep.prms.default_stop)) (final gradient norm: $(round(grad_norm, sigdigits=4)))."
    end

    is_sub_gaussian = extended ? BitVector(signs .< 0.0) : falses(n_components)
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

function _picard_y_loss_gpu(
    Y_gpu::AbstractMatrix{T},
    signs_gpu::AbstractVector{T},
    extended::Bool,
    loss_tmp::Union{AbstractMatrix{T},Nothing} = nothing,
) where {T<:AbstractFloat}
    N, T_len = size(Y_gpu)
    s = reshape(signs_gpu, N, 1)

    loss_mat = if loss_tmp === nothing
        extended ? (@. s * (abs(Y_gpu) + log1p(exp(T(-2.0) * abs(Y_gpu)))) + T(0.5) * Y_gpu^2) :
        (@. s * (abs(Y_gpu) + log1p(exp(T(-2.0) * abs(Y_gpu)))))
    else
        if extended
            @. loss_tmp = s * (abs(Y_gpu) + log1p(exp(T(-2.0) * abs(Y_gpu)))) + T(0.5) * Y_gpu^2
        else
            @. loss_tmp = s * (abs(Y_gpu) + log1p(exp(T(-2.0) * abs(Y_gpu))))
        end
        loss_tmp
    end

    row_sums = sum(loss_mat, dims = 2)
    return T(sum(Float64, Array(vec(row_sums))) / Float64(T_len))
end

function _picard_line_search_gpu!(
    Y::AbstractMatrix{T},
    W::AbstractMatrix{T},
    Y_new::AbstractMatrix{T},
    W_new::AbstractMatrix{T},
    DY::AbstractMatrix{T},
    DW::AbstractMatrix{T},
    direction::AbstractMatrix{T},
    direction_cpu::Matrix{T},
    transform_cpu::Matrix{T},
    signs::AbstractVector{T},
    loss_tmp::AbstractMatrix{T},
    current_loss::T,
    logdetW::T,
    ls_tries::Int,
    extended::Bool,
) where {T<:AbstractFloat}
    N = size(Y, 1)
    alpha = one(T)

    copyto!(direction_cpu, direction)
    mul!(DY, direction, Y)
    mul!(DW, direction, W)

    for _ = 1:ls_tries
        @. transform_cpu = alpha * direction_cpu
        @inbounds for i = 1:N
            transform_cpu[i, i] += one(T)
        end
        logdet_total = logdetW + T(logabsdet(transform_cpu)[1])

        @. Y_new = Y + alpha * DY
        trial_loss = -logdet_total + _picard_y_loss_gpu(Y_new, signs, extended, loss_tmp)

        if isfinite(trial_loss) && trial_loss < current_loss
            @. W_new = W + alpha * DW
            return true, trial_loss, alpha, logdet_total
        end
        alpha /= T(2.0)
    end

    copyto!(Y_new, Y)
    copyto!(W_new, W)
    return false, current_loss, alpha, logdetW
end

function _picard_solve_hessian_gpu!(
    Z::AbstractMatrix{T},
    h::AbstractMatrix{T},
    G::AbstractMatrix{T},
    I_diag::AbstractMatrix{Bool},
) where {T<:AbstractFloat}
    ht = transpose(h)
    Gt = transpose(G)
    @. Z = ifelse(I_diag, G / (h + one(T)), (ht * G - Gt) / (h * ht - one(T)))
    return nothing
end

function _picard_optimize_gpu!(W_cpu::Matrix{Float32}, dat_ica_cpu::Matrix{Float32}, params::IcaPrms, extended::Bool)
    N, T_len = size(dat_ica_cpu)
    m = params.picard_m
    tol = Float32(params.default_stop)
    lambda_min = Float32(params.picard_lambda_min)
    ls_tries = params.picard_ls_tries
    max_iter = params.max_iter

    dat_ica_gpu = gpu_array(copy(dat_ica_cpu))
    Y = copy(dat_ica_gpu)
    W = gpu_array(copy(W_cpu))

    Y_new, W_new = similar(Y), similar(W)
    DY, DW = similar(Y), similar(W)
    loss_tmp = similar(Y)
    psiY, psidY, Y_square = similar(Y), similar(Y), similar(Y)
    G, G_old, h = similar(W), similar(W), similar(W)
    direction, q, y_diff = similar(W), similar(W), similar(W)

    signs = gpu_array(ones(Float32, N))
    I_diag = gpu_array(Matrix{Bool}(I, N, N))

    h_cpu = zeros(Float32, N, N)
    G_cpu = zeros(Float32, N, N)
    direction_cpu = zeros(Float32, N, N)
    transform_cpu = zeros(Float32, N, N)
    signs_cpu = extended ? ones(Float32, N) : Float32[]
    old_signs_cpu = extended ? ones(Float32, N) : Float32[]
    mean_psid_gpu = extended ? gpu_array(zeros(Float32, N, 1)) : nothing
    mean_psid_cpu = extended ? zeros(Float32, N) : Float32[]
    C_cpu = extended ? zeros(Float32, N, N) : Matrix{Float32}(undef, 0, 0)

    lbfgs = LBFGSBuffer(G, m)
    logdetW = Float32(logabsdet(W_cpu)[1])
    current_loss = Float32(-logdetW + _picard_y_loss_gpu(Y, signs, extended, loss_tmp))
    gradient_norm = 0.0f0

    C = extended ? similar(W) : gpu_array(zeros(Float32, 0, 0))
    C_tmp = extended ? similar(W) : gpu_array(zeros(Float32, 0, 0))
    C_orig_cpu = extended ? zeros(Float32, N, N) : zeros(Float32, 0, 0)
    if extended
        mul!(C_orig_cpu, dat_ica_cpu, transpose(dat_ica_cpu))
        C_orig_cpu ./= Float32(T_len)
    end
    C_orig = gpu_array(C_orig_cpu)
    if extended
        mul!(C_tmp, W, C_orig)
        mul!(C, C_tmp, transpose(W))
    end

    for n = 1:max_iter
        @. psiY = tanh(Y)
        @. psidY = 1.0f0 - psiY^2

        mul!(G, psiY, transpose(Y))
        G ./= Float32(T_len)

        if extended
            sum!(mean_psid_gpu, psidY)
            copyto!(mean_psid_cpu, mean_psid_gpu)
            copyto!(C_cpu, C)
            copyto!(G_cpu, G)

            sign_change = (n == 1)
            for i = 1:N
                k_val = (mean_psid_cpu[i] / Float32(T_len)) * C_cpu[i, i] - G_cpu[i, i]
                s_val = k_val >= 0.0f0 ? 1.0f0 : -1.0f0
                s_val != old_signs_cpu[i] && (sign_change = true)
                signs_cpu[i] = s_val
                old_signs_cpu[i] = s_val
            end

            copyto!(signs, signs_cpu)
            @. G = G * signs + C
            @. psidY = psidY * signs + 1.0f0

            if sign_change
                current_loss = Float32(-logdetW + _picard_y_loss_gpu(Y, signs, extended, loss_tmp))
                _lbfgs_reset!(lbfgs)
            end
        end

        @. Y_square = Y^2
        mul!(h, psidY, transpose(Y_square))
        h ./= Float32(T_len)

        copyto!(h_cpu, h)
        _picard_regularize_hessian!(h_cpu, lambda_min)
        copyto!(h, h_cpu)

        copyto!(G_cpu, G)
        for i = 1:N
            G_cpu[i, i] -= 1.0f0
        end
        copyto!(G, G_cpu)

        gradient_norm = maximum(abs, G_cpu)
        if gradient_norm < tol
            @info "Picard converged at step $n (gradient norm: $gradient_norm)"
            copyto!(W_cpu, W)
            return true, n, signs_cpu, gradient_norm
        end

        if n > 1
            @. y_diff = G - G_old
            sy = dot(direction, y_diff)
            if sy > eps(Float32) * norm(direction) * norm(y_diff)
                _lbfgs_push!(lbfgs, direction, y_diff, 1.0f0 / sy)
            else
                _lbfgs_reset!(lbfgs)
            end
        end
        copyto!(G_old, G)

        _picard_l_bfgs_direction!(direction, q, G, h, lbfgs)

        ls_ok, new_loss, alpha, logdetW = _picard_line_search_gpu!(
            Y,
            W,
            Y_new,
            W_new,
            DY,
            DW,
            direction,
            direction_cpu,
            transform_cpu,
            signs,
            loss_tmp,
            current_loss,
            logdetW,
            ls_tries,
            extended,
        )

        if !ls_ok
            _picard_solve_hessian_gpu!(direction, h, G, I_diag)
            @. direction = -direction
            _lbfgs_reset!(lbfgs)
            ls_ok, new_loss, alpha, logdetW = _picard_line_search_gpu!(
                Y,
                W,
                Y_new,
                W_new,
                DY,
                DW,
                direction,
                direction_cpu,
                transform_cpu,
                signs,
                loss_tmp,
                current_loss,
                logdetW,
                ls_tries,
                extended,
            )

            if !ls_ok
                @info "Picard line search failed to find an improvement. Stopping early at step $n (gradient norm: $gradient_norm)."
                copyto!(W_cpu, W)
                return false, n, signs_cpu, gradient_norm
            end
        end

        @. direction *= alpha
        copyto!(Y, Y_new)
        copyto!(W, W_new)
        current_loss = new_loss

        n % 20 == 0 && mul!(Y, W, dat_ica_gpu)
        if extended
            mul!(C_tmp, W, C_orig)
            mul!(C, C_tmp, transpose(W))
        end

        if n == 1 || n % 10 == 0
            @info Printf.@sprintf("Picard step %d, gradient norm = %.7f, loss = %.7f", n, gradient_norm, current_loss)
        end
    end

    mul!(Y, W, dat_ica_gpu)
    copyto!(W_cpu, W)
    return false, max_iter, signs_cpu, gradient_norm
end
