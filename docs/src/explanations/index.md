# Explanations

Deep-dive explanations of EegFun.jl's architecture, design decisions, and core concepts.

## Data Structures

- [Data Structures](data-structures.md) - Core data types and their relationships

Understanding EegFun.jl's data structures is essential for effective use of the package. The data structures guide explains:

- **ContinuousData** - Raw, unprocessed EEG recordings
- **EpochData** - Segmented trials around events
- **ErpData** - Averaged event-related potentials
- **TimeFreqData** - Time-frequency decompositions
- **InfoIca** - Independent Component Analysis results
- **Layout** - Electrode positioning information

Each type is designed to work seamlessly with EegFun.jl's processing pipeline while maintaining flexibility for custom workflows.
