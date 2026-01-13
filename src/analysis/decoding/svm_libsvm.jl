"""
Direct LIBSVM implementation for decoding (no MLJ dependency).

This module provides a direct interface to LIBSVM.jl for SVM classification,
bypassing MLJ to reduce dependencies while maintaining the same interface.
"""

using LIBSVM
using Base: SubArray

"""
    libsvm_classifier(
        X_train::AbstractMatrix{Float64},
        y_train::Vector{Int},
        X_test::AbstractMatrix{Float64};
        C::Float64 = 1.0,
        kernel = LIBSVM.Kernel.Linear,
    ) -> Vector{Int}

Direct LIBSVM classifier matching the interface of `_mlj_classifier`.

This function uses LIBSVM.jl directly (no MLJ wrapper) for SVM classification.

# Arguments
- `X_train::AbstractMatrix{Float64}`: Training data [n_samples × n_features]
- `y_train::Vector{Int}`: Training labels [n_samples]
- `X_test::AbstractMatrix{Float64}`: Test data [n_samples × n_features]

# Keyword Arguments
- `C::Float64`: Regularization parameter (default: 1.0). Larger values = less regularization.
- `kernel`: Kernel type (default: `LIBSVM.Kernel.Linear`)

# Returns
- `Vector{Int}`: Predicted class labels

# Examples
```julia
y_pred = libsvm_classifier(X_train, y_train, X_test, C=1.0)
```
"""
function libsvm_classifier(
    X_train::AbstractMatrix{Float64},
    y_train::Vector{Int},
    X_test::AbstractMatrix{Float64};
    C::Float64 = 1.0,
    kernel = LIBSVM.Kernel.Linear,
)::Vector{Int}
    
    n_train = size(X_train, 1)
    n_test = size(X_test, 1)
    
    n_train == 0 && return Int[]
    n_test == 0 && return Int[]
    
    # Get unique classes
    classes = sort(unique(y_train))
    n_classes = length(classes)
    
    if n_classes < 2
        error("Need at least 2 classes, got $(n_classes)")
    end
    
    # Create class mapping (LIBSVM expects labels starting from 1)
    class_to_label = Dict(c => i for (i, c) in enumerate(classes))
    y_train_mapped = [class_to_label[y] for y in y_train]
    
    # Convert to SparseMatrixCSC (LIBSVM format)
    # LIBSVM expects features in columns, samples in rows (transpose)
    X_train_sparse = SparseArrays.sparse(X_train')
    X_test_sparse = SparseArrays.sparse(X_test')
    
    # Train model using svmtrain (LIBSVM.jl API)
    # Note: svmtype should be a Type (LIBSVM.SVC), and use 'cost' not 'C'
    model = LIBSVM.svmtrain(X_train_sparse, y_train_mapped; 
        svmtype = LIBSVM.SVC,
        kernel = kernel,
        cost = C,
    )
    
    # Predict
    y_pred_mapped, _ = LIBSVM.svmpredict(model, X_test_sparse)
    
    # Map predictions back to original class labels
    y_pred = [classes[pred] for pred in y_pred_mapped]
    
    return y_pred
end
