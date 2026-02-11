This demo explores the core functionality of segmenting continuous data into epochs (or trials) based on event markers, and basic manipulation of these segments.

### What are Epochs?

Epoching is the process of extracting specific time intervals around events (e.g., stimuli or responses) from a continuous recording. This allows for epoch-based analysis and averaging to reveal Event-Related Potentials (ERPs).

**Key features of EegFun's epoching**:
- Time-relative segmentation
- Multi-condition definition via trigger sequences

### Capabilities

- **Flexible Extraction**: Define intervals relative to trigger onset (e.g., -200ms to +1000ms).
- **Condition Matching**: Match triggers or sequences of triggers to specific condition names.

## Workflow Summary

1. **Data Preparation**: Load continuous data, apply layout, and perform basic preprocessing (filtering, re-referencing).
2. **Define Conditions**: Use `EpochCondition` to specify which triggers belong to which experimental condition.
3. **Extraction**: Use `extract_epochs()` to create segments. This returns a collection of `EpochedData` objects.
