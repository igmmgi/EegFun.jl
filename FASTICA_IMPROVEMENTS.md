# FastICA Algorithm Improvement Suggestions

Based on comparison with MultiVariateStats.jl and NeuroAnalyzer.jl implementations.

## Key Improvements

### 1. Switch to Symmetric FastICA (Recommended)

**Current**: Deflationary approach (one component at a time)
**Better**: Symmetric approach (all components simultaneously)

**Benefits**:
- Faster execution (typically 2-3x)
- More numerically stable
- Better parallelization potential
- Matches industry-standard implementations (sklearn, MultiVariateStats)

### 2. Use Symmetric Decorrelation

**Current**: Gram-Schmidt orthogonalization (error accumulation)
**Better**: `W = W * (W'W)^{-1/2}` (symmetric decorrelation)

**Benefits**:
- More numerically stable
- Avoids error accumulation from sequential orthogonalization
- Standard in FastICA literature

### 3. Vectorized Contrast Function Computation

**Current**: Element-wise loops for each component
**Better**: Process all components in parallel using matrix operations

**Benefits**:
- Better cache utilization
- Leverages BLAS operations
- Faster execution

### 4. Improved Convergence Checking

**Current**: Checks each component individually
**Better**: Check all components simultaneously using `max(abs.(abs.(diag(W*Wp')) .- 1))`

## Implementation Notes

### Helper Function Needed

```julia
"""
    _invsqrtm!(C::AbstractMatrix{<:Real})

Compute inv(sqrtm(C)) through symmetric eigenvalue decomposition.
In-place version that modifies C.
"""
function _invsqrtm!(C::AbstractMatrix{<:Real})
    n = size(C, 1)
    size(C, 2) == n || error("C must be a square matrix.")
    E = eigen!(Symmetric(C))
    U = E.vectors
    evs = E.values
    for i = 1:n
        @inbounds evs[i] = 1.0 / sqrt(sqrt(evs[i]))
    end
    rmul!(U, Diagonal(evs))
    return U * transpose(U)
end
```

### Contrast Function Types (like MultiVariateStats)

Consider using type-based dispatch for contrast functions:

```julia
abstract type ICAGDeriv end

struct Tanh{T} <: ICAGDeriv
    a::T
end

struct Gaus <: ICAGDeriv end

function update!(f::Tanh{T}, U::AbstractMatrix{T}, E::AbstractVector{T}) where {T}
    n, k = size(U)
    a = f.a
    @inbounds for j in 1:k
        _s = zero(T)
        @fastmath for i in 1:n
            t = tanh(a * U[i,j])
            U[i,j] = t
            _s += a * (1 - t^2)
        end
        E[j] = _s / n
    end
end
```

## Performance Comparison

Expected improvements:
- **Speed**: 2-3x faster for typical EEG data (64 channels, 100k samples)
- **Stability**: Better convergence, especially for many components
- **Memory**: Similar or slightly better (pre-allocated matrices)

## Migration Path

1. Add `_invsqrtm!` helper function
2. Implement symmetric FastICA as alternative (keep deflationary as fallback)
3. Add contrast function types for better extensibility
4. Benchmark both approaches
5. Make symmetric the default if performance is better

## Code Structure (Symmetric FastICA)

```julia
# Pre-allocated storage (like MultiVariateStats)
Wp = similar(W)                # previous version of W
U  = Matrix{Float64}(undef, n_samples, n_components)  # w'x & g(w'x)
Y  = Matrix{Float64}(undef, n_channels, n_components)  # E{x g(w'x)}
E1 = Vector{Float64}(undef, n_components)              # E{g'(w'x)}

# Main loop
for iteration in 1:max_iter
    copyto!(Wp, W)
    
    # Apply W: U = X' * W (all components at once)
    mul!(U, transpose(dat_ica), W)
    
    # Compute g(w'x) and E{g'(w'x)} for all components
    update!(contrast_fun, U, E1)
    
    # Compute E{x g(w'x)} for all components
    rmul!(mul!(Y, dat_ica, U), 1 / n_samples)
    
    # Update all components: Y - E1 .* W
    for j in 1:n_components
        @. W[:, j] = Y[:, j] - E1[j] * W[:, j]
    end
    
    # Symmetric decorrelation
    copyto!(W, W * _invsqrtm!(W'W))
    
    # Check convergence (all components at once)
    chg = maximum(abs.(abs.(diag(W*Wp')) .- 1))
    converged = (chg < tol)
    
    if converged
        break
    end
end
```

