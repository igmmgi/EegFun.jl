"""
Pegasos (Primal Estimated sub-GrAdient SOlver for SVM) implementation.

This module provides a pure Julia implementation of linear SVM using the Pegasos algorithm,
which solves the primal SVM optimization problem using stochastic gradient descent.

The implementation is designed to be a drop-in replacement for MLJ/LIBSVM classifiers
in the decoding pipeline, accepting matrices and integer labels directly.
"""

using LinearAlgebra
using Random

"""
    PegasosSVM

Trained Pegasos linear SVM model.

# Fields
- `weights::Vector{Float64}`: Weight vector (w)
- `bias::Float64`: Bias term (b)
- `classes::Vector{Int}`: Unique class labels seen during training
- `n_features::Int`: Number of features
"""
struct PegasosSVM
    weights::Vector{Float64}
    bias::Float64
    classes::Vector{Int}
    n_features::Int
end

"""
    train_pegasos(
        X_train::AbstractMatrix{Float64},
        y_train::Vector{Int};
        C::Float64 = 1.0,
        max_iter::Int = 1000,
        tolerance::Float64 = 1e-6,
        rng::AbstractRNG = Random.GLOBAL_RNG,
    ) -> PegasosSVM

Train a linear SVM using the Pegasos algorithm.

# Arguments
- `X_train::AbstractMatrix{Float64}`: Training data [n_samples × n_features]
- `y_train::Vector{Int}`: Training labels [n_samples]

# Keyword Arguments
- `C::Float64`: Regularization parameter (default: 1.0). Larger values = less regularization.
- `max_iter::Int`: Maximum number of iterations (default: 500, may need more for high-dimensional data)
- `tolerance::Float64`: Convergence tolerance (default: 1e-6)
- `rng::AbstractRNG`: Random number generator for shuffling (default: GLOBAL_RNG)

# Returns
- `PegasosSVM`: Trained model

# Algorithm
Pegasos solves the primal SVM problem:
    min_w (λ/2)||w||² + (1/n) Σ max(0, 1 - y_i(w·x_i + b))

where λ = 1/C and the algorithm uses stochastic gradient descent.

# Examples
```julia
X = rand(100, 10)
y = [rand([-1, 1]) for _ in 1:100]
model = train_pegasos(X, y, C=1.0, max_iter=1000)
```
"""
function train_pegasos(
    X_train::AbstractMatrix{Float64},
    y_train::Vector{Int};
    C::Float64 = 1.0,
    max_iter::Int = 1000,
    tolerance::Float64 = 1e-6,
    rng::AbstractRNG = Random.GLOBAL_RNG,
)::PegasosSVM

    n_samples, n_features = size(X_train)
    length(y_train) == n_samples || error("X_train and y_train must have same number of samples")

    # Get unique classes and check if binary
    classes = unique(y_train)
    n_classes = length(classes)

    if n_classes < 2
        error("Need at least 2 classes, got $(n_classes)")
    end

    # For binary classification, convert labels to -1 and +1
    if n_classes == 2
        # Map classes to -1 and +1
        class_map = Dict(classes[1] => -1, classes[2] => 1)
        y_binary = [class_map[y] for y in y_train]
        
        # Train binary classifier
        w, b = _train_binary_pegasos(X_train, y_binary, C, max_iter, tolerance, rng)
        
        return PegasosSVM(w, b, classes, n_features)
    else
        # Multi-class: use one-vs-rest
        error("Multi-class classification not yet implemented. Use binary classification (2 classes).")
    end
end

"""
    _train_binary_pegasos(
        X::AbstractMatrix{Float64},
        y::Vector{Int},  # Must be -1 or +1
        C::Float64,
        max_iter::Int,
        tolerance::Float64,
        rng::AbstractRNG,
    ) -> (Vector{Float64}, Float64)

Internal function to train binary Pegasos SVM.

Returns (weights, bias).
"""
function _train_binary_pegasos(
    X::AbstractMatrix{Float64},
    y::Vector{Int},  # Must be -1 or +1
    C::Float64,
    max_iter::Int,
    tolerance::Float64,
    rng::AbstractRNG,
)::Tuple{Vector{Float64}, Float64}

    n_samples, n_features = size(X)
    λ = 1.0 / C  # Regularization parameter

    # Initialize weights and bias
    w = zeros(Float64, n_features)
    b = 0.0

    # Pre-allocate indices vector (reused each iteration to avoid allocation)
    indices = collect(1:n_samples)
    
    # Pre-compute λ for faster access
    inv_λ = 1.0 / λ  # Store 1/λ to avoid division

    # Main training loop - optimized for speed
    # Use per-epoch learning rate instead of per-sample to prevent it from becoming too small
    for iter in 1:max_iter
        # Learning rate per epoch: η = 1 / (λ * iter) 
        # This prevents learning rate from decaying too fast
        η = inv_λ / iter
        
        # Shuffle data each epoch (in-place to avoid allocation)
        shuffle!(rng, indices)

        for idx in indices
            x_i = @view X[idx, :]
            y_i = y[idx]

            # Compute margin: y_i * (w·x_i + b)
            margin = y_i * (dot(w, x_i) + b)

            # Update if margin < 1 (misclassified or within margin)
            if margin < 1.0
                # Gradient step: w = (1 - η * λ) * w + η * y_i * x_i
                # Use in-place operations with BLAS
                scale = 1.0 - η * λ
                w .*= scale  # In-place scaling (faster)
                LinearAlgebra.axpy!(η * y_i, x_i, w)  # w += η * y_i * x_i (BLAS, in-place)
                b += η * y_i
            else
                # Shrink weights (regularization step) - in-place
                scale = 1.0 - η * λ
                w .*= scale
            end
        end

        # Skip convergence checks for maximum speed - just use fixed iterations
        # (LIBSVM doesn't do expensive convergence checks either)
    end

    return w, b
end

"""
    predict(model::PegasosSVM, X_test::AbstractMatrix{Float64}) -> Vector{Int}

Predict class labels for test data.

# Arguments
- `model::PegasosSVM`: Trained Pegasos model
- `X_test::AbstractMatrix{Float64}`: Test data [n_samples × n_features]

# Returns
- `Vector{Int}`: Predicted class labels

# Examples
```julia
model = train_pegasos(X_train, y_train)
predictions = predict(model, X_test)
```
"""
function predict(model::PegasosSVM, X_test::AbstractMatrix{Float64})::Vector{Int}
    n_test = size(X_test, 1)
    n_test == 0 && return Int[]

    # Check feature dimension matches
    size(X_test, 2) == model.n_features || 
        error("X_test has $(size(X_test, 2)) features, but model expects $(model.n_features)")

    # For binary classification
    if length(model.classes) == 2
        predictions = Vector{Int}(undef, n_test)
        
        for i in 1:n_test
            x_i = @view X_test[i, :]
            # Compute decision function: w·x + b
            decision = dot(model.weights, x_i) + model.bias
            
            # Predict class based on sign
            if decision >= 0
                predictions[i] = model.classes[2]  # Positive class
            else
                predictions[i] = model.classes[1]  # Negative class
            end
        end
        
        return predictions
    else
        error("Multi-class prediction not yet implemented")
    end
end

"""
    pegasos_classifier(
        X_train::AbstractMatrix{Float64},
        y_train::Vector{Int},
        X_test::AbstractMatrix{Float64};
        C::Float64 = 1.0,
        max_iter::Int = 1000,
        tolerance::Float64 = 1e-6,
        rng::AbstractRNG = Random.GLOBAL_RNG,
    ) -> Vector{Int}

Convenience function matching the interface of `_mlj_classifier`.

This function trains a Pegasos SVM and returns predictions, matching the exact
interface expected by the decoding pipeline.

# Arguments
- `X_train::AbstractMatrix{Float64}`: Training data [n_samples × n_features]
- `y_train::Vector{Int}`: Training labels [n_samples]
- `X_test::AbstractMatrix{Float64}`: Test data [n_samples × n_features]

# Keyword Arguments
- `C::Float64`: Regularization parameter (default: 1.0)
- `max_iter::Int`: Maximum iterations (default: 20, optimized for speed)
- `tolerance::Float64`: Convergence tolerance (default: 1e-6)
- `rng::AbstractRNG`: Random number generator (default: GLOBAL_RNG)

# Returns
- `Vector{Int}`: Predicted class labels

# Examples
```julia
y_pred = pegasos_classifier(X_train, y_train, X_test, C=1.0)
```
"""
function pegasos_classifier(
    X_train::AbstractMatrix{Float64},
    y_train::Vector{Int},
    X_test::AbstractMatrix{Float64};
    C::Float64 = 1.0,
        max_iter::Int = 500,
    tolerance::Float64 = 1e-6,
    rng::AbstractRNG = Random.GLOBAL_RNG,
)::Vector{Int}
    model = train_pegasos(X_train, y_train; C=C, max_iter=max_iter, tolerance=tolerance, rng=rng)
    return predict(model, X_test)
end
