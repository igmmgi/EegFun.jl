# === Picard ===
# === PICARD ICA IMPLEMENTATION (Preconditioned ICA for Real Data) ===

# Fast SIMD-friendly approximation of log(1 + exp(-2|y|))
@inline function _fast_log1p_exp(ax::T) where {T<:AbstractFloat}
    if ax > T(6.5)
        return zero(T)
    elseif ax > T(3.0)
        return exp(T(-2.0) * ax)
    else
        return @evalpoly(
            ax,
            T(0.69310808),
            T(-0.99854606),
            T(0.48711652),
            T(0.04765150),
            T(-0.17336261),
            T(0.09241291),
            T(-0.02471741),
            T(0.00345703),
            T(-0.00020136)
        )
    end
end

# Compute the Y-dependent loss terms, optionally fused with a linear step (Y + α*DY).
# thread_sums is a pre-allocated buffer of length ≥ nthreads() to avoid per-call allocation.
function _picard_y_loss(
    Y::Matrix{T},
    signs::Vector{T},
    extended::Bool,
    thread_sums::Vector{Float64};
    DY::Matrix{T} = Y,
    alpha::T = zero(T),
) where {T<:AbstractFloat}
    N_ch, T_len = size(Y)
    fused = alpha != zero(T)
    nth = Threads.nthreads()
    chunk_size = div(T_len, nth)
    fill!(thread_sums, 0.0)

    Threads.@threads for tid = 1:nth
        j_start = (tid - 1) * chunk_size + 1
        j_end = (tid == nth) ? T_len : tid * chunk_size

        acc = 0.0
        if extended
            @inbounds for j = j_start:j_end
                for i = 1:N_ch
                    y = fused ? Y[i, j] + alpha * DY[i, j] : Y[i, j]
                    abs_y = abs(y)
                    term = _fast_log1p_exp(abs_y)
                    acc += Float64(signs[i] * (abs_y + term) + T(0.5) * y^2)
                end
            end
        else
            @inbounds for j = j_start:j_end
                for i = 1:N_ch
                    y = fused ? Y[i, j] + alpha * DY[i, j] : Y[i, j]
                    abs_y = abs(y)
                    term = _fast_log1p_exp(abs_y)
                    acc += Float64(signs[i] * (abs_y + term))
                end
            end
        end
        thread_sums[tid] = acc
    end

    return T(sum(thread_sums) / Float64(T_len))
end

# Full Picard loss: -log|det(W)| + mean Y-dependent terms
function _picard_loss(Y::Matrix{T}, W::Matrix{T}, signs::Vector{T}, extended::Bool, thread_sums::Vector{Float64}) where {T<:AbstractFloat}
    logdetW, _ = logabsdet(W)
    return T(-logdetW + _picard_y_loss(Y, signs, extended, thread_sums))
end


# Solve the 2×2 block Hessian system element-wise.
function _picard_solve_hessian!(Z::AbstractMatrix{T}, h::AbstractMatrix{T}, G::AbstractMatrix{T}) where {T<:AbstractFloat}
    ht = transpose(h)
    Gt = transpose(G)
    @. Z = (ht * G - Gt) / (h * ht - one(T))
end

# Regularize off-diagonal Hessian elements
function _picard_regularize_hessian!(h::Matrix{T}, lambda_min::T) where {T<:AbstractFloat}
    N = size(h, 1)
    for j = 1:N
        for i = (j+1):N
            hij = h[i, j]
            hji = h[j, i]
            discr = sqrt((hij - hji)^2 + T(4.0))
            eig = T(0.5) * (hij + hji - discr)
            if eig < lambda_min
                update = lambda_min - eig
                h[i, j] += update
                h[j, i] += update
            end
        end
    end
end

# Pre-allocated circular buffer for L-BFGS history (avoids push!/popfirst!/copy! allocations)
mutable struct LBFGSBuffer{T<:AbstractFloat, MatT<:AbstractMatrix{T}}
    s_buf::Vector{MatT}   # pre-allocated direction history
    y_buf::Vector{MatT}   # pre-allocated gradient diff history
    r_buf::Vector{T}      # pre-allocated curvature scalars
    a_buf::Vector{T}      # pre-allocated alpha coefficients
    m::Int                # max history length
    len::Int              # current number of stored entries
    head::Int             # next write position (1-indexed)
end

function LBFGSBuffer(T::Type{<:AbstractFloat}, N::Int, m::Int)
    s_buf = [zeros(T, N, N) for _ = 1:m]
    y_buf = [zeros(T, N, N) for _ = 1:m]
    r_buf = zeros(T, m)
    a_buf = zeros(T, m)
    return LBFGSBuffer{T, Matrix{T}}(s_buf, y_buf, r_buf, a_buf, m, 0, 1)
end

function LBFGSBuffer(sample_mat::MatT, m::Int) where {T<:AbstractFloat, MatT<:AbstractMatrix{T}}
    N = size(sample_mat, 1)
    s_buf = [similar(sample_mat, N, N) for _ = 1:m]
    y_buf = [similar(sample_mat, N, N) for _ = 1:m]
    for k = 1:m
        fill!(s_buf[k], zero(T))
        fill!(y_buf[k], zero(T))
    end
    r_buf = zeros(T, m)
    a_buf = zeros(T, m)
    return LBFGSBuffer{T, MatT}(s_buf, y_buf, r_buf, a_buf, m, 0, 1)
end

@inline function _lbfgs_reset!(buf::LBFGSBuffer)
    buf.len = 0
    buf.head = 1
end

function _lbfgs_push!(buf::LBFGSBuffer{T, MatT}, direction::MatT, y_diff::MatT, r_val::T) where {T, MatT<:AbstractMatrix{T}}
    idx = buf.head
    copyto!(buf.s_buf[idx], direction)
    copyto!(buf.y_buf[idx], y_diff)
    buf.r_buf[idx] = r_val
    buf.head = (idx % buf.m) + 1
    buf.len = min(buf.len + 1, buf.m)
end

# Get the k-th most recent entry (1 = most recent, len = oldest)
@inline function _lbfgs_get(buf::LBFGSBuffer, k::Int)
    # head-1 is the most recently written slot
    idx = mod1(buf.head - k, buf.m)
    return buf.s_buf[idx], buf.y_buf[idx], buf.r_buf[idx]
end

# L-BFGS two-loop recursion with Hessian preconditioning
function _picard_l_bfgs_direction!(
    Z::MatT,
    q::MatT,
    G::MatT,
    h::MatT,
    buf::LBFGSBuffer{T, MatT},
) where {T<:AbstractFloat, MatT<:AbstractMatrix{T}}
    copyto!(q, G)

    m_len = buf.len
    # Backward pass: most recent (k=1) to oldest (k=m_len)
    for k = 1:m_len
        s, y, r = _lbfgs_get(buf, k)
        a = r * dot(s, q)
        buf.a_buf[k] = a
        @. q -= a * y
    end

    _picard_solve_hessian!(Z, h, q)

    # Forward pass: oldest (k=m_len) to most recent (k=1)
    for k = m_len:-1:1
        s, y, r = _lbfgs_get(buf, k)
        a = buf.a_buf[k]
        beta = r * dot(y, Z)
        @. Z += (a - beta) * s
    end

    @. Z = -Z
end


# Line search with fused loss computation and incremental logdet
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
        for i = 1:N
            transform[i, i] += one(T)
        end
        step_logdet, _ = logabsdet(transform)
        logdet_total = logdetW + T(step_logdet)

        y_loss = _picard_y_loss(Y, signs, extended, thread_sums; DY = DY, alpha = alpha)
        new_loss = -logdet_total + y_loss

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

# _picard_reset_lbfgs! is now replaced by _lbfgs_reset! on the LBFGSBuffer

# === CPU Optimization Loop ===

function _picard_optimize_cpu!(W::Matrix{T}, dat_ica::Matrix{T}, params::IcaPrms, extended::Bool) where {T<:AbstractFloat}
    N, T_len = size(dat_ica)
    m = params.picard_m
    tol = T(params.default_stop)
    lambda_min = T(params.picard_lambda_min)
    ls_tries = params.picard_ls_tries
    max_iter = params.max_iter

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

    G     = zeros(T, N, N)
    G_old = zeros(T, N, N)
    h     = zeros(T, N, N)

    direction = zeros(T, N, N)
    q = zeros(T, N, N)
    Z = zeros(T, N, N)

    signs     = ones(T, N)
    old_signs = ones(T, N)

    # Pre-allocated L-BFGS circular buffer (avoids copy/push!/popfirst! allocations)
    lbfgs = LBFGSBuffer(T, N, m)

    # Pre-allocated thread_sums buffer for _picard_y_loss
    thread_sums = zeros(Float64, Threads.nthreads())

    logdetW_val, _ = logabsdet(W)
    logdetW = T(logdetW_val)
    current_loss = _picard_loss(Y, W, signs, extended, thread_sums)
    sign_change = false

    # Extended mode covariance: C = W * C_orig * W'
    # C_orig = (dat_ica * dat_ica') / T_len
    C = extended ? zeros(T, N, N) : zeros(T, 0, 0)
    C_orig = extended ? zeros(T, N, N) : zeros(T, 0, 0)
    K = extended ? zeros(T, N) : zeros(T, 0)

    # Temporary buffer for C = W * C_orig * W'
    C_tmp = extended ? zeros(T, N, N) : zeros(T, 0, 0)

    if extended
        mul!(C_orig, dat_ica, transpose(dat_ica))
        C_orig ./= T(T_len)

        mul!(C_tmp, W, C_orig)
        mul!(C, C_tmp, transpose(W))
    end

    for n = 1:max_iter
        Threads.@threads for j = 1:T_len
            @inbounds @simd for i = 1:N
                yij = Y[i, j]
                @fastmath ty = tanh(yij)
                psiY[i, j] = ty
                psidY[i, j] = one(T) - ty^2
                Y_square[i, j] = yij * yij
            end
        end

        mul!(G, psiY, transpose(Y))
        G ./= T(T_len)

        if extended
            for i = 1:N
                mean_psid = sum(view(psidY, i, :)) / T_len
                K[i] = T(mean_psid) * C[i, i] - G[i, i]
            end

            @. signs = sign(K)
            if n > 1
                sign_change = false
                @inbounds for i = 1:N
                    if signs[i] != old_signs[i]
                        sign_change = true
                        break
                    end
                end
            end
            copyto!(old_signs, signs)

            for i = 1:N
                s = signs[i]
                for j_col = 1:N
                    G[i, j_col] *= s
                end
            end
            @. G += C

            Threads.@threads for j = 1:T_len
                @inbounds for i = 1:N
                    psidY[i, j] = psidY[i, j] * signs[i] + one(T)
                end
            end
        end

        mul!(h, psidY, transpose(Y_square))
        h ./= T(T_len)

        _picard_regularize_hessian!(h, lambda_min)

        for i = 1:N
            G[i, i] -= one(T)
        end

        gradient_norm = maximum(abs, G)
        if gradient_norm < tol
            @info "Picard converged at step $n (gradient norm: $gradient_norm)"
            return n
        end

        if n > 1
            @. y_diff = G - G_old
            _lbfgs_push!(lbfgs, direction, y_diff, one(T) / dot(direction, y_diff))
        end
        copyto!(G_old, G)

        if extended && sign_change
            current_loss = _picard_loss(Y, W, signs, extended, thread_sums)
            _lbfgs_reset!(lbfgs)
        end

        _picard_l_bfgs_direction!(Z, q, G, h, lbfgs)
        copyto!(direction, Z)

        converged, new_loss, alpha, logdetW = _picard_line_search(
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

        if !converged
            @. direction = -G
            _lbfgs_reset!(lbfgs)
            converged, new_loss, alpha, logdetW = _picard_line_search(
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

            if !converged
                @info "Line search failed to find an improvement. Stopping early."
                return n
            end
        end

        # Scale direction by alpha to store the actual step taken in L-BFGS history
        @. direction *= alpha

        copyto!(Y, Y_new)
        copyto!(W, W_new)
        current_loss = new_loss

        if extended
            # Exact C update: C = W_new * C_orig * W_new'
            mul!(C_tmp, W, C_orig)
            mul!(C, C_tmp, transpose(W))
        end

        if n == 1 || n % 10 == 0
            @info Printf.@sprintf("Picard step %d, gradient norm = %.7f, loss = %.7f", n, gradient_norm, current_loss)
        end
    end
    return max_iter
end

# =============================================================================

# === Public API ===

function picard_ica(
    dat_ica::Matrix{Float64},
    layout::Layout,
    filename::String;
    n_components::Int,
    extended::Bool = false,
    params::IcaPrms = IcaPrms(),
)
    ext_str = extended ? "Extended-" : ""
    @info "Running $(ext_str)Picard ICA (use_gpu=$(params.use_gpu)): $(size(dat_ica,1)) channels x $(size(dat_ica,2)) samples -> $n_components components"

    n_samples = size(dat_ica, 2)

    original_mean = vec(mean(dat_ica, dims = 2))

    dat_ica .-= original_mean
    scale = sqrt(norm((dat_ica * transpose(dat_ica)) / n_samples))
    dat_ica ./= scale

    F = svd(dat_ica)
    pca_components = F.U[:, 1:n_components]

    eigenvalues = (F.S[1:n_components] .^ 2) ./ (n_samples - 1)
    sphere = diagm(1.0 ./ sqrt.(eigenvalues))

    transform_matrix = sphere * transpose(pca_components)
    dat_ica_sphered = Matrix{Float64}(undef, n_components, n_samples)
    mul!(dat_ica_sphered, transform_matrix, dat_ica)
    dat_ica = dat_ica_sphered

    gpu_active = false
    if params.use_gpu
        if is_gpu_available()
            gpu_active = true
            @info "[GPU ACTIVATED] Running Picard ICA on $(gpu_device_name())..."
        else
            @minimal_warning "Requested GPU acceleration (use_gpu=true), but no functional GPU package (CUDA.jl, AMDGPU.jl, Metal.jl) has been loaded. Please run 'using CUDA' or 'using AMDGPU' before calling run_ica. Falling back to CPU."
        end
    end

    if gpu_active
        # GPU path strictly uses Float32 for memory bandwidth and speed
        W_init = Matrix{Float32}(I, n_components, n_components)
        dat_ica_f32 = Float32.(dat_ica)
        _picard_optimize_gpu!(W_init, dat_ica_f32, params, extended)
        W_final = Float64.(W_init)
    else
        # CPU path uses Float64 to stay consistent with CPU Infomax behavior
        W_init = Matrix{Float64}(I, n_components, n_components)
        _picard_optimize_cpu!(W_init, dat_ica, params, extended)
        W_final = W_init
    end

    unmixing_matrix = W_final * sphere * pca_components'
    mixing_matrix = pinv(unmixing_matrix)

    meanvar = vec(sum(abs2, mixing_matrix, dims = 1) .* sum(abs2, dat_ica, dims = 2)' ./ (n_components * n_samples - 1))
    meanvar_normalized = meanvar ./ sum(meanvar)
    order = sortperm(meanvar_normalized, rev = true)

    return InfoIca(
        filename,
        unmixing_matrix[order, :],
        mixing_matrix[:, order],
        sphere,
        meanvar_normalized[order],
        scale,
        original_mean,
        [Symbol("IC$i") for i = 1:n_components],
        Dict{Int,Matrix{Float64}}(),
        layout,
        falses(n_components),
    )
end

# === GPU versions of Picard ===

# GPU versions of Picard
# GPU native implementation of Picard Y-dependent loss (supports optional preallocated loss_tmp buffer)
function _picard_y_loss_gpu(
    Y_gpu::AbstractMatrix{T},
    signs_gpu::AbstractVector{T},
    extended::Bool,
    loss_tmp::Union{AbstractMatrix{T}, Nothing} = nothing,
) where {T<:AbstractFloat}
    N, T_len = size(Y_gpu)
    s = reshape(signs_gpu, N, 1)

    if loss_tmp === nothing
        if extended
            val = sum(@. s * (abs(Y_gpu) + log1p(exp(T(-2.0) * abs(Y_gpu)))) + T(0.5) * Y_gpu^2)
        else
            val = sum(@. s * (abs(Y_gpu) + log1p(exp(T(-2.0) * abs(Y_gpu)))))
        end
    else
        if extended
            @. loss_tmp = s * (abs(Y_gpu) + log1p(exp(T(-2.0) * abs(Y_gpu)))) + T(0.5) * Y_gpu^2
        else
            @. loss_tmp = s * (abs(Y_gpu) + log1p(exp(T(-2.0) * abs(Y_gpu))))
        end
        val = sum(loss_tmp)
    end

    return val / T(T_len)
end

# GPU line search helper with pre-allocated CPU logdet calculation
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
        for i = 1:N
            transform_cpu[i, i] += one(T)
        end
        step_logdet, _ = logabsdet(transform_cpu)
        logdet_total = logdetW + T(step_logdet)

        @. Y_new = Y + alpha * DY
        y_loss = _picard_y_loss_gpu(Y_new, signs, extended, loss_tmp)
        trial_loss = -logdet_total + y_loss

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

# GPU Optimization Loop
# =============================================================================

function _picard_optimize_gpu!(W_cpu::Matrix{Float32}, dat_ica_cpu::Matrix{Float32}, params::IcaPrms, extended::Bool)
    N, T_len = size(dat_ica_cpu)
    m = params.picard_m
    tol = Float32(params.default_stop)
    lambda_min = Float32(params.picard_lambda_min)
    ls_tries = params.picard_ls_tries
    max_iter = params.max_iter

    Y = gpu_array(copy(dat_ica_cpu))
    W = gpu_array(copy(W_cpu))

    Y_new = similar(Y)
    W_new = similar(W)
    DY = similar(Y)
    DW = similar(W)
    loss_tmp = similar(Y)

    psiY = similar(Y)
    psidY = similar(Y)
    Y_square = similar(Y)

    G     = similar(W)
    G_old = similar(W)
    h     = similar(W)

    direction = similar(W)
    q = similar(W)
    y_diff = similar(W)

    signs     = gpu_array(ones(Float32, N))
    old_signs = gpu_array(ones(Float32, N))

    # Pre-allocated CPU buffers to eliminate per-iteration host/device allocations
    h_cpu = zeros(Float32, N, N)
    G_cpu = zeros(Float32, N, N)
    direction_cpu = zeros(Float32, N, N)
    transform_cpu = zeros(Float32, N, N)
    W_cpu_tmp = copy(W_cpu)
    signs_cpu = extended ? ones(Float32, N) : Float32[]
    old_signs_cpu = extended ? ones(Float32, N) : Float32[]
    mean_psid_gpu = extended ? gpu_array(zeros(Float32, N, 1)) : nothing
    mean_psid_cpu = extended ? zeros(Float32, N) : Float32[]
    C_cpu = extended ? zeros(Float32, N, N) : Matrix{Float32}(undef, 0, 0)

    # Pre-allocated L-BFGS circular buffer (avoids push!/popfirst!/copy! GPU allocations)
    lbfgs = LBFGSBuffer(G, m)

    logdetW_val, _ = logabsdet(W_cpu_tmp)
    logdetW = Float32(logdetW_val)
    current_loss = Float32(-logdetW + _picard_y_loss_gpu(Y, signs, extended, loss_tmp))
    sign_change = false

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

            sign_change = false
            for i = 1:N
                k_val = (mean_psid_cpu[i] / Float32(T_len)) * C_cpu[i, i] - G_cpu[i, i]
                s_val = sign(k_val)
                if n > 1 && s_val != old_signs_cpu[i]
                    sign_change = true
                end
                signs_cpu[i] = s_val
                old_signs_cpu[i] = s_val
            end

            copyto!(signs, signs_cpu)

            @. G = G * signs + C
            @. psidY = psidY * signs + 1.0f0
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
            break
        end

        if n > 1
            @. y_diff = G - G_old
            r_val = 1.0f0 / dot(direction, y_diff)
            _lbfgs_push!(lbfgs, direction, y_diff, r_val)
        end
        copyto!(G_old, G)

        if extended && sign_change
            current_loss = Float32(-logdetW + _picard_y_loss_gpu(Y, signs, extended, loss_tmp))
            _lbfgs_reset!(lbfgs)
        end

        _picard_l_bfgs_direction!(direction, q, G, h, lbfgs)

        converged, new_loss, alpha, logdetW = _picard_line_search_gpu!(
            Y, W, Y_new, W_new, DY, DW, direction,
            direction_cpu, transform_cpu, signs, loss_tmp,
            current_loss, logdetW, ls_tries, extended,
        )

        if !converged
            @. direction = -G
            _lbfgs_reset!(lbfgs)

            converged, new_loss, alpha, logdetW = _picard_line_search_gpu!(
                Y, W, Y_new, W_new, DY, DW, direction,
                direction_cpu, transform_cpu, signs, loss_tmp,
                current_loss, logdetW, ls_tries, extended,
            )

            if !converged
                @info "Line search failed to find an improvement. Stopping early."
                break
            end
        end

        # Scale direction by alpha to store the actual step taken in L-BFGS history
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

    copyto!(W_cpu, W)
end

