# Epoch Selection

`EegFun.jl` provides a system for selecting epochs from continuous data. Whether your experimental paradigm uses simple single triggers or complex sequences of events with timing constraints, `EegFun.jl` can handle it.

This tutorial covers the various ways you can define epochs, from the simplest single-trigger approach to advanced pattern matching.

## The `EpochCondition` Structure

At the heart of epoch selection is the `EpochCondition` structure. It allows you to define exactly what constitutes an epoch in your study.

```julia
@kwdef struct EpochCondition
    name::String
    trigger_sequences::Vector{Vector{Union{Int,Symbol,UnitRange{Int}}}}
    reference_index::Int = 1
    timing_pairs::Union{Nothing,Vector{Tuple{Int,Int}}} = nothing
    min_interval::Union{Nothing,Float64} = nothing
    max_interval::Union{Nothing,Float64} = nothing
    after::Union{Nothing,Int} = nothing
    before::Union{Nothing,Int} = nothing
end
```

## Simple Epoching: Single Values

The most common case is extracting an epoch around a single trigger value.

```julia
# Define a condition for trigger 1
condition = EpochCondition(
    name = "Stimulus", 
    trigger_sequences = [[1]]
)

# Extract epochs (-200ms to 1000ms around the trigger)
epochs = extract_epochs(dat, [condition], (-0.2, 1.0))
```

To match any of multiple trigger values for a single condition (OR logic), provide multiple sequences:

```julia
# Match trigger 1 OR trigger 2
condition = EpochCondition(
    name = "AllStimuli", 
    trigger_sequences = [[1], [2]]
)
```

## Sequence Matching

Sometimes an "event" is actually a sequence of triggers. For example, a target stimulus (1) followed by a response (10).

```julia
# Match the sequence [1, 10]
condition = EpochCondition(
    name = "TargetResponse", 
    trigger_sequences = [[1, 10]]
)
```

### Wildcards and Ranges

You can use wildcards and ranges within your sequences:

- `:any`: Matches any trigger value.
- `UnitRange` (e.g., `1:10`): Matches any value within the range.

```julia
# Match trigger 1, then any trigger, then trigger 3
condition = EpochCondition(
    name = "Wildcard", 
    trigger_sequences = [[1, :any, 3]]
)

# Match any trigger from 1 to 5, then trigger 10
condition = EpochCondition(
    name = "Range", 
    trigger_sequences = [[1:5, 10]]
)
```

## Onset (t=0) Reference

By default, the first trigger in a sequence is considered $t=0$. You can change this using `reference_index`.

```julia
# Sequence: [Warning (1), Stimulus (2), Response (10)]
# We want t=0 to be the Stimulus (index 2)
condition = EpochCondition(
    name = "StimulusOnset",
    trigger_sequences = [[1, 2, 10]],
    reference_index = 2
)
```

## Timing Constraints

You can restrict matches to sequences where triggers occur within specific time intervals. This is useful for filtering out trials with late responses or accidental double-triggers.

```julia
# Match [1, 10] only if they occur between 200ms and 800ms apart
condition = EpochCondition(
    name = "ValidResponse",
    trigger_sequences = [[1, 10]],
    timing_pairs = [(1, 2)], # Calculate interval between 1st and 2nd trigger
    min_interval = 0.2,
    max_interval = 0.8
)
```

## Position Constraints

You can also filter sequences based on whether they occur before or after certain "marker" triggers.

```julia
# Only find sequences that occur AFTER trigger 99 (e.g., start of a block)
condition = EpochCondition(
    name = "Block2",
    trigger_sequences = [[1, 2]],
    after = 99
)

# Only find sequences that occur BEFORE trigger 88 (e.g., end of experimental phase)
condition = EpochCondition(
    name = "Phase1",
    trigger_sequences = [[1, 2]],
    before = 88
)
```

> [!IMPORTANT]
> You cannot specify both `after` and `before` in a single `EpochCondition`.

## External TOML Configuration

For complex studies, defining your epoch conditions in an external TOML file is often cleaner. This separates your analysis parameters from your code.

### The TOML Format

Create a file (e.g., `epochs.toml`):

```toml
[epochs]
[[epochs.conditions]]
name = "Standard"
trigger_sequences = [[1]]

[[epochs.conditions]]
name = "Target"
trigger_sequences = [[2, 10]]
reference_index = 1
timing_pairs = [[1, 2]]
min_interval = 0.1
max_interval = 1.0

[[epochs.conditions]]
name = "BlockStart"
trigger_sequences = [[1:5]]
after = 100
```

### Loading and Using the TOML

```julia
using TOML

# Load the configuration
config = TOML.parsefile("epochs.toml")

# Parse into EpochCondition objects
conditions = condition_parse_epoch(config)

# Extract epochs
epochs = extract_epochs(dat, conditions, (-0.2, 1.0))
```
