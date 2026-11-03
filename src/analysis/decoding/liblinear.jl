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
    cost <= 0 && @minimal_error("SVM cost parameter must be positive, got $cost")
    (any(!isfinite, X_train) || any(!isfinite, X_test)) && @minimal_error("Input features contain NaN or Inf values")

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
        bias = 1.0,
    )
    w = model.w
    labels = model.labels
    nr_c = Int(model.nr_class)
    n_features = size(X_test, 2)
    has_bias = model.bias >= 0
    y_pred = Vector{Int}(undef, n_test)

    if nr_c == 2
        l1 = labels[1]
        l2 = labels[2]
        b = has_bias ? w[n_features+1] * model.bias : 0.0
        @inbounds for i = 1:n_test
            val = b
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
                val = has_bias ? w[n_features*nr_c+c] * model.bias : 0.0
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
Pre-allocated workspace for LIBLINEAR in the inner decoding loop (with intercept `bias = 1.0`).
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
    # Extra row for bias node (index n_feat + 1, value 1.0) plus sentinel node (-1, NaN)
    n_node_rows = max_features + 2
    nodes = Matrix{LIBSVM.LIBLINEAR.FeatureNode}(undef, n_node_rows, max_samples)
    for i = 1:max_samples
        for j = 1:max_features
            nodes[j, i] = LIBSVM.LIBLINEAR.FeatureNode(Cint(j), 0.0)
        end
        nodes[max_features+1, i] = LIBSVM.LIBLINEAR.FeatureNode(Cint(max_features + 1), 1.0)
        nodes[max_features+2, i] = LIBSVM.LIBLINEAR.FeatureNode(Cint(-1), NaN)
    end
    nodeptrs = [pointer(nodes, (i - 1) * n_node_rows + 1) for i = 1:max_samples]
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
    problem = [LIBSVM.LIBLINEAR.Problem(Cint(max_samples), Cint(max_features + 1), pointer(y_double), pointer(nodeptrs), 1.0)]
    labels_buf = Vector{Cint}(undef, max_classes)
    w_dim = max_classes == 2 ? 1 : max_classes
    w_buf = Vector{Float64}(undef, (max_features + 1) * w_dim)
    return LiblinearWorkspace(nodes, nodeptrs, y_double, param, problem, labels_buf, w_buf)
end

"""
Zero-allocation C training call using the pre-allocated workspace (with `bias = 1.0`).
Directly loads model weights (including intercept at `n_feat + 1`) and label assignments,
and properly destroys the C `Model` struct via `free_and_destroy_model`.
Supports both binary and multi-class classification.
"""
function _train_prealloc_zero!(ws::LiblinearWorkspace, X_tr::AbstractMatrix{Float64}, y_tr::AbstractVector{Int})
    n_feat = size(X_tr, 1)
    n_samp = size(X_tr, 2)

    @inbounds for i = 1:n_samp
        for j = 1:n_feat
            ws.nodes[j, i] = LIBSVM.LIBLINEAR.FeatureNode(Cint(j), X_tr[j, i])
        end
        ws.nodes[n_feat+1, i] = LIBSVM.LIBLINEAR.FeatureNode(Cint(n_feat + 1), 1.0)
        ws.nodes[n_feat+2, i] = LIBSVM.LIBLINEAR.FeatureNode(Cint(-1), NaN)
        ws.y_double[i] = Float64(y_tr[i])
    end

    nr_c = 0
    GC.@preserve ws begin
        ws.problem[1] = LIBSVM.LIBLINEAR.Problem(Cint(n_samp), Cint(n_feat + 1), pointer(ws.y_double), pointer(ws.nodeptrs), 1.0)

        ptr = ccall(
            (:train, LIBSVM.LIBLINEAR.liblinear),
            Ptr{LIBSVM.LIBLINEAR.Model},
            (Ptr{LIBSVM.LIBLINEAR.Problem}, Ptr{LIBSVM.LIBLINEAR.Parameter}),
            ws.problem,
            ws.param,
        )
        m = unsafe_load(ptr)

        nr_c = Int(m.nr_class)
        w_dim = Int(m.nr_feature) + (m.bias >= 0 ? 1 : 0)
        w_number = nr_c == 2 ? 1 : nr_c
        total_w = w_dim * w_number
        unsafe_copyto!(pointer(ws.w_buf), m.w, total_w)

        ws.labels_buf .= -1
        ccall((:get_labels, LIBSVM.LIBLINEAR.liblinear), Cvoid, (Ptr{LIBSVM.LIBLINEAR.Model}, Ptr{Cint}), ptr, pointer(ws.labels_buf))
        ptr_ref = Ref{Ptr{LIBSVM.LIBLINEAR.Model}}(ptr)
        ccall((:free_and_destroy_model, LIBSVM.LIBLINEAR.liblinear), Cvoid, (Ptr{Ptr{LIBSVM.LIBLINEAR.Model}},), ptr_ref)
    end
    return ws.w_buf, ws.labels_buf, nr_c
end
