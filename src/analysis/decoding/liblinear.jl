# Low-level LIBLINEAR C bindings and pre-allocated workspaces for MVPA decoding.
#
# Provides zero-allocation C calls directly to LIBLINEAR for high-throughput
# cross-validated SVM decoding across time points and folds.

"""
    libsvm_classifier(
        X_train::AbstractMatrix{Float64},
        y_train::AbstractVector{Int},
        X_test::AbstractMatrix{Float64};
        cost::Float64 = 1.0,
    ) -> Vector{Int}

LIBSVM classifier for decoding.

Uses LIBSVM's linear solver directly with fast in-place decision value evaluation.

# Arguments
- `X_train::AbstractMatrix{Float64}`: Training data [n_samples × n_features]
- `y_train::AbstractVector{Int}`: Training labels [n_samples]
- `X_test::AbstractMatrix{Float64}`: Test data [n_samples × n_features]

# Keyword Arguments
- `cost::Float64`: Regularization parameter (default: 1.0). Larger values = less regularization.

# Returns
- `Vector{Int}`: Predicted class labels

# Examples
```julia
y_pred = libsvm_classifier(X_train, y_train, X_test, cost=1.0)
```
"""
function libsvm_classifier(
    X_train::AbstractMatrix{Float64},
    y_train::AbstractVector{Int},
    X_test::AbstractMatrix{Float64};
    cost::Real = 1.0,
)

    n_train = size(X_train, 1)
    n_test = size(X_test, 1)

    n_train == 0 && return Int[]
    n_test == 0 && return Int[]

    LIBSVM.LIBLINEAR.set_print(false)
    model = LIBSVM.LIBLINEAR.linear_train(
        y_train,
        transpose(X_train);
        solver_type = Int32(LIBSVM.Linearsolver.L2R_L2LOSS_SVC),
        C = Float64(cost),
    )
    w = model.w
    labels = model.labels
    nr_c = Int(model.nr_class)
    n_features = size(X_test, 2)
    y_pred = Vector{Int}(undef, n_test)

    if nr_c == 2
        l1 = labels[1]
        l2 = labels[2]
        @inbounds for i = 1:n_test
            val = 0.0
            for f = 1:n_features
                val += w[f] * X_test[i, f]
            end
            y_pred[i] = val > 0.0 ? l1 : l2
        end
    else
        @inbounds for i = 1:n_test
            best_c = 1
            max_val = -Inf
            for c = 1:nr_c
                val = 0.0
                for f = 1:n_features
                    val += w[(f-1)*nr_c+c] * X_test[i, f]
                end
                if val > max_val
                    max_val = val
                    best_c = c
                end
            end
            y_pred[i] = labels[best_c]
        end
    end

    return y_pred
end


# === PRE-ALLOCATED LIBLINEAR WORKSPACE ===

"""
Pre-allocated workspace for LIBLINEAR in the inner decoding loop.
Eliminates all per-fold/timepoint memory allocations during training.
"""
struct LiblinearWorkspace
    nodes::Matrix{LIBSVM.LIBLINEAR.FeatureNode}
    nodeptrs::Vector{Ptr{LIBSVM.LIBLINEAR.FeatureNode}}
    y_double::Vector{Float64}
    param::Vector{LIBSVM.LIBLINEAR.Parameter}
    problem::Vector{LIBSVM.LIBLINEAR.Problem}
    labels_buf::Vector{Cint}
    w_buf::Vector{Float64}
end

function LiblinearWorkspace(max_features::Int, max_samples::Int, cost::Real, max_classes::Int = 2)
    nodes = Matrix{LIBSVM.LIBLINEAR.FeatureNode}(undef, max_features + 1, max_samples)
    for i = 1:max_samples
        for j = 1:max_features
            nodes[j, i] = LIBSVM.LIBLINEAR.FeatureNode(Cint(j), 0.0)
        end
        nodes[max_features+1, i] = LIBSVM.LIBLINEAR.FeatureNode(Cint(-1), NaN)
    end
    nodeptrs = [pointer(nodes, (i - 1) * (max_features + 1) + 1) for i = 1:max_samples]
    y_double = Vector{Float64}(undef, max_samples)

    solver_type = Int32(LIBSVM.Linearsolver.L2R_L2LOSS_SVC)
    eps = 0.01
    param = [
        LIBSVM.LIBLINEAR.Parameter(
            solver_type,
            eps,
            Float64(cost),
            Cint(0),
            convert(Ptr{Cint}, C_NULL),
            convert(Ptr{Float64}, C_NULL),
            0.1,
            0.5,
            convert(Ptr{Float64}, C_NULL),
            Cint(1),
        ),
    ]
    problem = [LIBSVM.LIBLINEAR.Problem(Cint(max_samples), Cint(max_features), pointer(y_double), pointer(nodeptrs), -1.0)]
    labels_buf = Vector{Cint}(undef, max_classes)
    w_dim = max_classes == 2 ? 1 : max_classes
    w_buf = Vector{Float64}(undef, max_features * w_dim)
    return LiblinearWorkspace(nodes, nodeptrs, y_double, param, problem, labels_buf, w_buf)
end

"""
Zero-allocation C training call using the pre-allocated workspace.
Directly loads model weights and label assignments without wrapping in Julia arrays.
Supports both binary and multi-class classification.
"""
function _train_prealloc_zero!(ws::LiblinearWorkspace, X_tr::AbstractMatrix{Float64}, y_tr::AbstractVector{Int})
    n_feat = size(X_tr, 1)
    n_samp = size(X_tr, 2)

    @inbounds for i = 1:n_samp
        for j = 1:n_feat
            ws.nodes[j, i] = LIBSVM.LIBLINEAR.FeatureNode(Cint(j), X_tr[j, i])
        end
        ws.y_double[i] = Float64(y_tr[i])
    end

    ws.problem[1] = LIBSVM.LIBLINEAR.Problem(Cint(n_samp), Cint(n_feat), pointer(ws.y_double), pointer(ws.nodeptrs), -1.0)

    ptr = ccall(
        (:train, LIBSVM.LIBLINEAR.liblinear),
        Ptr{LIBSVM.LIBLINEAR.Model},
        (Ptr{LIBSVM.LIBLINEAR.Problem}, Ptr{LIBSVM.LIBLINEAR.Parameter}),
        ws.problem,
        ws.param,
    )
    m = unsafe_load(ptr)

    nr_c = Int(m.nr_class)
    w_dim = Int(m.nr_feature)
    w_number = nr_c == 2 ? 1 : nr_c
    total_w = w_dim * w_number
    unsafe_copyto!(pointer(ws.w_buf), m.w, total_w)

    ws.labels_buf .= -1
    ccall((:get_labels, LIBSVM.LIBLINEAR.liblinear), Cvoid, (Ptr{LIBSVM.LIBLINEAR.Model}, Ptr{Cint}), ptr, pointer(ws.labels_buf))
    ccall((:free_model_content, LIBSVM.LIBLINEAR.liblinear), Cvoid, (Ptr{LIBSVM.LIBLINEAR.Model},), ptr)
    return ws.w_buf, ws.labels_buf, nr_c
end
