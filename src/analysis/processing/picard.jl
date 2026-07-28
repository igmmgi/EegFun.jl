# =============================================================================
# PICARD ICA IMPLEMENTATION (Preconditioned ICA for Real Data)
# =============================================================================

function _picard_loss(Y::Matrix{Float32}, W::Matrix{Float32}, signs::Vector{Float32}, extended::Bool)
    N, T = size(Y)
    
    logdetW, _ = logabsdet(W)
    
    chan_loss = zeros(Float32, N)
    
    if extended
        chan_sq = zeros(Float32, N)
        Threads.@threads for i in 1:N
            loss_i = 0.0f0
            sq_i = 0.0f0
            @inbounds @simd for j in 1:T
                @fastmath val = abs(Y[i, j])
                loss_i += val + log1p(exp(-2.0f0 * val))
                sq_i += Y[i, j]^2
            end
            chan_loss[i] = loss_i
            chan_sq[i] = sq_i
        end
        
        loss_y = 0.0f0
        @inbounds for i in 1:N
            loss_y += signs[i] * (chan_loss[i] / Float32(T)) + 0.5f0 * (chan_sq[i] / Float32(T))
        end
    else
        Threads.@threads for i in 1:N
            loss_i = 0.0f0
            @inbounds @simd for j in 1:T
                @fastmath val = abs(Y[i, j])
                loss_i += val + log1p(exp(-2.0f0 * val))
            end
            chan_loss[i] = loss_i
        end
        
        loss_y = 0.0f0
        @inbounds for i in 1:N
            loss_y += signs[i] * (chan_loss[i] / Float32(T))
        end
    end
    
    return Float32(-logdetW + loss_y)
end

function _picard_loss(Y::AbstractMatrix{Float32}, W::AbstractMatrix{Float32}, signs::AbstractVector{Float32}, extended::Bool)
    N, T = size(Y)
    
    W_cpu = Array(W)
    logdetW, _ = logabsdet(W_cpu)
    
    chan_loss = mapreduce(y -> abs(y) + log1p(exp(-2.0f0 * abs(y))), +, Y, dims=2) ./ Float32(T)
    total_y_loss = sum(signs .* vec(chan_loss))
    
    if extended
        # Extended adds 0.5 * Y^2 / T to the loss
        chan_sq = mapreduce(y -> y^2, +, Y, dims=2) ./ Float32(T)
        total_y_loss += 0.5f0 * sum(vec(chan_sq))
    end
    
    return Float32(-logdetW + total_y_loss)
end

function _picard_solve_hessian!(Z::AbstractMatrix{Float32}, h::AbstractMatrix{Float32}, h_off::AbstractMatrix{Float32}, G::AbstractMatrix{Float32})
    # det = h * h' - h_off * h_off'
    # return (h' * G - h_off * G') / det
    # Element-wise operations! h_off is just ones(N, 1) * ones(1, N) in our case, so it's a matrix of 1s.
    ht = transpose(h)
    Gt = transpose(G)
    @. Z = (ht * G - Gt) / (h * ht - 1.0f0)
end

function _picard_regularize_hessian!(h::AbstractMatrix{Float32}, lambda_min::Float32)
    N = size(h, 1)
    # discr = sqrt.((h - h').^2 .+ 4.0)
    # eigenvalues = 0.5 .* (h .+ h' .- discr)
    # regularize where eigenvalues < lambda_min (off-diagonals only)
    for j in 1:N
        for i in (j+1):N
            hij = h[i, j]
            hji = h[j, i]
            discr = sqrt((hij - hji)^2 + 4.0f0)
            eig = 0.5f0 * (hij + hji - discr)
            if eig < lambda_min
                update = lambda_min - eig
                h[i, j] += update
                h[j, i] += update
            end
        end
    end
end

function _picard_l_bfgs_direction!(
    Z::AbstractMatrix{Float32},
    q::AbstractMatrix{Float32}, 
    G::AbstractMatrix{Float32}, 
    h::AbstractMatrix{Float32}, 
    s_list::Vector{<:AbstractMatrix{Float32}}, 
    y_list::Vector{<:AbstractMatrix{Float32}}, 
    r_list::Vector{Float32}, 
    a_list::Vector{Float32}
)
    copyto!(q, G)
    empty!(a_list)
    
    m_len = length(s_list)
    for k in m_len:-1:1
        s = s_list[k]
        y = y_list[k]
        r = r_list[k]
        
        alpha = r * dot(s, q)
        push!(a_list, alpha)
        @. q -= alpha * y
    end
    
    # Solve Hessian
    _picard_solve_hessian!(Z, h, h, q) # h_off = 1 implicit
    
    # Reverse a_list logic manually since we pushed in reverse
    for k in 1:m_len
        s = s_list[k]
        y = y_list[k]
        r = r_list[k]
        alpha = a_list[m_len - k + 1]
        
        beta = r * dot(y, Z)
        @. Z += (alpha - beta) * s
    end
    
    @. Z = -Z
end

function _picard_line_search(
    Y::AbstractMatrix{Float32}, 
    W::AbstractMatrix{Float32}, 
    Y_new::AbstractMatrix{Float32},
    W_new::AbstractMatrix{Float32},
    transform::AbstractMatrix{Float32},
    direction::AbstractMatrix{Float32}, 
    signs::AbstractVector{Float32}, 
    current_loss::Float32, 
    ls_tries::Int, 
    extended::Bool
)
    N = size(W, 1)
    alpha = 1.0f0
    
    for _ in 1:ls_tries
        # transform = I + alpha * direction
        @. transform = alpha * direction
        transform[diagind(transform)] .+= 1.0f0
        
        mul!(Y_new, transform, Y)
        mul!(W_new, transform, W)
        
        new_loss = _picard_loss(Y_new, W_new, signs, extended)
        
        if isfinite(new_loss) && new_loss < current_loss
            return true, new_loss, alpha
        end
        alpha /= 2.0f0
    end
    
    # Revert Y_new and W_new to old values since line search failed
    copyto!(Y_new, Y)
    copyto!(W_new, W)
    return false, current_loss, alpha
end

function _picard_optimize_cpu!(
    W::Matrix{Float32}, 
    dat_ica::Matrix{Float32}, 
    params::IcaPrms, 
    extended::Bool
)
    N, T = size(dat_ica)
    m = params.picard_m
    tol = Float32(params.default_stop)
    lambda_min = Float32(params.picard_lambda_min)
    ls_tries = params.picard_ls_tries
    max_iter = params.max_iter
    
    Y = copy(dat_ica)
    Y_new = similar(Y)
    W_new = similar(W)
    transform = zeros(Float32, N, N)
    
    psiY = zeros(Float32, N, T)
    psidY = zeros(Float32, N, T)
    Y_square = zeros(Float32, N, T)
    
    G = zeros(Float32, N, N)
    G_old = zeros(Float32, N, N)
    h = zeros(Float32, N, N)
    
    direction = zeros(Float32, N, N)
    q = zeros(Float32, N, N)
    Z = zeros(Float32, N, N)
    
    signs = ones(Float32, N)
    old_signs = ones(Float32, N)
    
    s_list = Matrix{Float32}[]
    y_list = Matrix{Float32}[]
    r_list = Float32[]
    a_list = Float32[]
    
    current_loss = _picard_loss(Y, W, signs, extended)
    sign_change = false
    
    C = zeros(Float32, N, N)
    if extended
        mul!(C, dat_ica, transpose(dat_ica))
        C ./= Float32(T)
    end
    
    K = zeros(Float32, N)
    
    for n in 1:max_iter
        # Score and derivative
        Threads.@threads for j in 1:T
            @inbounds @simd for i in 1:N
                @fastmath ty = tanh(Y[i, j])
                psiY[i, j] = ty
                psidY[i, j] = 1.0f0 - ty^2
            end
        end
        
        # G = psiY * Y^T / T
        mul!(G, psiY, transpose(Y))
        G ./= Float32(T)
        
        if extended
            # K = mean(psidY, axis=1) * diag(C) - diag(G)
            for i in 1:N
                mean_psid = sum(view(psidY, i, :)) / T
                K[i] = mean_psid * C[i, i] - G[i, i]
            end
            
            @. signs = sign(K)
            if n > 1
                sign_change = any(signs .!= old_signs)
            end
            copyto!(old_signs, signs)
            
            # G *= signs[:, None], psidY *= signs[:, None]
            for i in 1:N
                s = signs[i]
                for j in 1:N
                    G[i, j] *= s
                end
                for j in 1:T
                    psidY[i, j] *= s
                end
            end
            
            @. G += C
            @. psidY += 1.0f0
        end
        
        # h = psidY * Y_square^T / T
        @. Y_square = Y^2
        mul!(h, psidY, transpose(Y_square))
        h ./= Float32(T)
        
        _picard_regularize_hessian!(h, lambda_min)
        
        # G -= eye(N)
        for i in 1:N
            G[i, i] -= 1.0f0
        end
        
        gradient_norm = maximum(abs.(G))
        if gradient_norm < tol
            @info "Picard converged at step $n (gradient norm: $gradient_norm)"
            break
        end
        
        # Memory update
        if n > 1
            push!(s_list, copy(direction))
            push!(y_list, G .- G_old)
            push!(r_list, 1.0f0 / dot(direction, G .- G_old))
            
            if length(s_list) > m
                popfirst!(s_list)
                popfirst!(y_list)
                popfirst!(r_list)
            end
        end
        copyto!(G_old, G)
        
        if extended && sign_change
            current_loss = _picard_loss(Y, W, signs, extended)
            empty!(s_list)
            empty!(y_list)
            empty!(r_list)
        end
        
        _picard_l_bfgs_direction!(Z, q, G, h, s_list, y_list, r_list, a_list)
        copyto!(direction, Z)
        
        converged, new_loss, alpha = _picard_line_search(Y, W, Y_new, W_new, transform, direction, signs, current_loss, ls_tries, extended)
        
        if !converged
            @. direction = -G
            empty!(s_list)
            empty!(y_list)
            empty!(r_list)
            converged, new_loss, alpha = _picard_line_search(Y, W, Y_new, W_new, transform, direction, signs, current_loss, ls_tries, extended)
            
            if !converged
                @info "Line search failed to find an improvement. Stopping early."
                break
            end
        end
        
        copyto!(Y, Y_new)
        copyto!(W, W_new)
        current_loss = new_loss
        
        if extended
            # C = W * covariance * W^T
            # covariance is original C, so we need to compute W * (dat * dat' / T) * W'
            # which is just (W * dat) * (W * dat)' / T = Y * Y^T / T
            mul!(C, Y, transpose(Y))
            C ./= Float32(T)
        end
        
        current_loss = new_loss
        
        if n == 1 || n % 10 == 0
            @info Printf.@sprintf(
                "Picard step %d, gradient norm = %.7f, loss = %.7f",
                n, gradient_norm, current_loss
            )
        end
    end
end

function _picard_optimize_gpu!(
    W_cpu::Matrix{Float32}, 
    dat_ica_cpu::Matrix{Float32}, 
    params::IcaPrms, 
    extended::Bool
)
    N, T = size(dat_ica_cpu)
    m = params.picard_m
    tol = Float32(params.default_stop)
    lambda_min = Float32(params.picard_lambda_min)
    ls_tries = params.picard_ls_tries
    max_iter = params.max_iter
    
    Y = gpu_array(dat_ica_cpu)
    W = gpu_array(W_cpu)
    
    Y_new = gpu_array(zeros(Float32, N, T))
    W_new = gpu_array(zeros(Float32, N, N))
    transform = gpu_array(zeros(Float32, N, N))
    
    psiY = gpu_array(zeros(Float32, N, T))
    psidY = gpu_array(zeros(Float32, N, T))
    Y_square = gpu_array(zeros(Float32, N, T))
    
    G = gpu_array(zeros(Float32, N, N))
    G_old = gpu_array(zeros(Float32, N, N))
    h = gpu_array(zeros(Float32, N, N))
    
    direction = gpu_array(zeros(Float32, N, N))
    q = gpu_array(zeros(Float32, N, N))
    Z = gpu_array(zeros(Float32, N, N))
    
    signs = gpu_array(ones(Float32, N))
    old_signs = gpu_array(ones(Float32, N))
    
    # L-BFGS history lives on GPU
    s_list = typeof(G)[]
    y_list = typeof(G)[]
    r_list = Float32[]
    a_list = Float32[]
    
    current_loss = _picard_loss(Y, W, signs, extended)
    sign_change = false
    
    C = gpu_array(zeros(Float32, N, N))
    if extended
        mul!(C, Y, transpose(Y)) # since W is identity originally, Y = dat_ica
        C ./= Float32(T)
    end
    
    K = gpu_array(zeros(Float32, N))
    
    for n in 1:max_iter
        # Score and derivative
        @. psiY = tanh(Y)
        @. psidY = 1.0f0 - psiY^2
        
        # G = psiY * Y^T / T
        mul!(G, psiY, transpose(Y))
        G ./= Float32(T)
        
        if extended
            # K = mean(psidY, axis=1) * diag(C) - diag(G)
            mean_psid = vec(sum(psidY, dims=2) ./ Float32(T))
            mean_psid_cpu = Array(mean_psid)
            C_cpu = Array(C)
            G_cpu = Array(G)
            K_cpu = zeros(Float32, N)
            
            for i in 1:N
                K_cpu[i] = mean_psid_cpu[i] * C_cpu[i, i] - G_cpu[i, i]
            end
            
            K .= gpu_array(K_cpu)
            @. signs = sign(K)
            
            if n > 1
                sign_change = any(Array(signs) .!= Array(old_signs))
            end
            copyto!(old_signs, signs)
            
            # G *= signs[:, None], psidY *= signs[:, None]
            @. G = G * signs
            @. psidY = psidY * signs
            
            @. G += C
            @. psidY += 1.0f0
        end
        
        # h = psidY * Y_square^T / T
        @. Y_square = Y^2
        mul!(h, psidY, transpose(Y_square))
        h ./= Float32(T)
        
        # Regularize on CPU because of off-diagonal eigenvalue logic
        h_cpu = Array(h)
        _picard_regularize_hessian!(h_cpu, lambda_min)
        h .= gpu_array(h_cpu)
        
        # G -= eye(N)
        G_cpu = Array(G)
        for i in 1:N
            G_cpu[i, i] -= 1.0f0
        end
        G .= gpu_array(G_cpu)
        
        gradient_norm = maximum(abs.(G_cpu))
        if gradient_norm < tol
            @info "Picard converged at step $n (gradient norm: $gradient_norm)"
            break
        end
        
        # Memory update
        if n > 1
            push!(s_list, copy(direction))
            push!(y_list, G .- G_old)
            push!(r_list, 1.0f0 / dot(direction, G .- G_old))
            
            if length(s_list) > m
                popfirst!(s_list)
                popfirst!(y_list)
                popfirst!(r_list)
            end
        end
        copyto!(G_old, G)
        
        if extended && sign_change
            current_loss = _picard_loss(Y, W, signs, extended)
            empty!(s_list)
            empty!(y_list)
            empty!(r_list)
        end
        
        _picard_l_bfgs_direction!(Z, q, G, h, s_list, y_list, r_list, a_list)
        copyto!(direction, Z)
        
        converged, new_loss, alpha = _picard_line_search(Y, W, Y_new, W_new, transform, direction, signs, current_loss, ls_tries, extended)
        
        if !converged
            @. direction = -G
            empty!(s_list)
            empty!(y_list)
            empty!(r_list)
            converged, new_loss, alpha = _picard_line_search(Y, W, Y_new, W_new, transform, direction, signs, current_loss, ls_tries, extended)
            
            if !converged
                @info "Line search failed to find an improvement. Stopping early."
                break
            end
        end
        
        copyto!(Y, Y_new)
        copyto!(W, W_new)
        current_loss = new_loss
        
        if extended
            mul!(C, Y, transpose(Y))
            C ./= Float32(T)
        end
        
        current_loss = new_loss
        
        if n == 1 || n % 10 == 0
            @info Printf.@sprintf(
                "Picard step %d, gradient norm = %.7f, loss = %.7f",
                n, gradient_norm, current_loss
            )
        end
    end
    
    copyto!(W_cpu, Array(W))
end

function picard_ica(
    dat_ica::Matrix{Float64},
    layout::Layout,
    filename::String;
    n_components::Int,
    extended::Bool = false,
    params::IcaPrms = IcaPrms(),
)
    # Print Info
    ext_str = extended ? "Extended-" : ""
    @info "Running $(ext_str)Picard ICA (use_gpu=$(params.use_gpu)): $(size(dat_ica,1)) channels x $(size(dat_ica,2)) samples -> $n_components components"

    n_channels = size(dat_ica, 1)
    n_samples = size(dat_ica, 2)

    # Store original mean before removing it
    original_mean = vec(mean(dat_ica, dims = 2))

    # Center and scale data
    dat_ica .-= original_mean
    scale = sqrt(norm((dat_ica * transpose(dat_ica)) / n_samples))
    dat_ica ./= scale

    # PCA reduction using high-precision SVD
    F = svd(dat_ica)
    pca_components = F.U[:, 1:n_components]

    # Analytic sphering using exact singular values
    eigenvalues = (F.S[1:n_components] .^ 2) ./ (n_samples - 1)
    sphere = diagm(1.0 ./ sqrt.(eigenvalues))

    transform_matrix = sphere * transpose(pca_components)
    dat_ica_sphered = Matrix{Float64}(undef, n_components, n_samples)
    mul!(dat_ica_sphered, transform_matrix, dat_ica)
    dat_ica = dat_ica_sphered

    # Detect GPU hardware backend
    gpu_active = false
    if params.use_gpu
        if is_gpu_available()
            gpu_active = true
            @info "[GPU ACTIVATED] Running Picard ICA on $(gpu_device_name())..."
        else
            @minimal_warning "Requested GPU acceleration (use_gpu=true), but no functional GPU package (CUDA.jl, AMDGPU.jl, Metal.jl) has been loaded. Please run 'using CUDA' or 'using AMDGPU' before calling run_ica. Falling back to CPU."
        end
    end

    # Initialize W as Identity (or decorrelated random). FastICA/Infomax start with random.
    # Picard usually starts with Identity since data is already sphered.
    W_init = Matrix{Float32}(I, n_components, n_components)
    dat_ica_f32 = Float32.(dat_ica)

    if gpu_active
        _picard_optimize_gpu!(W_init, dat_ica_f32, params, extended)
    else
        _picard_optimize_cpu!(W_init, dat_ica_f32, params, extended)
    end

    W_final = Float64.(W_init)

    # Compute final matrices
    unmixing_matrix = W_final * sphere * pca_components'
    mixing_matrix = pinv(unmixing_matrix)

    # calculate total variance explained and order
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
        falses(n_components)
    )
end
