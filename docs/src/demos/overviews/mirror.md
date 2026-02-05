This demo demonstrates mirror padding, a technique for extending time-series data by reflecting it at the edges.

### What is Mirror Padding?

Mirror padding extends your data at the beginning and/or end by creating a reflected copy:

```

Original:     [A B C D E]
Pre-padding:  [E D C B A | A B C D E]
Post-padding: [A B C D E | E D C B A]
Both:         [E D C B A | A B C D E | E D C B A]
```

This creates smooth, continuous edges that avoid discontinuities.

### Why Use Mirror Padding?

**Edge Artifact Reduction**:

Filtering operations can introduce artifacts at the edges of your data because:

- Filters "see" discontinuities at boundaries
- Edge effects propagate inward
- Lost information at boundaries

Mirror padding solves this by:

- Creating smooth, continuous edges
- Allowing filters to process edges naturally
- Removing padding after filtering to discard corrupted edge regions

**Typical Workflow**:

1. **Mirror** your data (add padding)
2. **Filter** (artifacts now occur in padding region)
3. **Remove padding** (discard corrupted edges, keep clean data)

### Padding Options

The `mirror()` function accepts three padding modes:

| Mode | Description |
|------|-------------|
| **:pre** | Mirror before the start |
| **:post** | Mirror after the end |
| **:both** | Mirror at both ends |

### Use Cases

**High-pass filtering**:

High-pass filters are particularly sensitive to edge effects. Mirror padding helps preserve data near epoch boundaries.

**Time-frequency analysis**:

Wavelet transforms and spectral analysis benefit from smooth edges.

**Concatenated epochs**:

When processing multiple epochs sequentially, padding prevents edge artifacts between epochs.

### Workflow Summary

This demo shows:

1. **Load and preprocess** continuous data
2. **Extract epochs** (-200 to 1000 ms)
3. **Apply mirror padding** with different modes (:pre, :post, :both)
4. **Visualize padded epochs** to see the effect
5. **Average to ERPs** and apply padding
6. **Visualize padded ERPs**

The visual comparison clearly shows how padding extends the time window while maintaining smooth edges.
