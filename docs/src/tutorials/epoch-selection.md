# Epoch Selection

`EegFun.jl` provides a system for selecting epochs from continuous data. This includes experimental paradigms that use simple single triggers or complex sequences of events with timing constraints.

This tutorial covers the various ways you can define epochs, from the simplest single-trigger approach to advanced pattern matching.

## Quick Reference

| Feature | Field(s) | Example |
|---------|----------|---------|
| Single trigger | `trigger_sequences = [[1]]` | Epoch around trigger 1 |
| Multiple triggers (OR) | `trigger_sequences = [[1], [2]]` | Either trigger 1 or 2 |
| Trigger sequence | `trigger_sequences = [[1, 10]]` | Trigger 1 followed by 10 |
| Wildcard | `trigger_sequences = [[1, :any, 3]]` | 1, then anything, then 3 |
| Range | `trigger_sequences = [[1:5, 10]]` | Any of 1–5, then 10 |
| t=0 reference | `reference_index = 2` | Align to 2nd trigger in sequence |
| Timing constraint | `timing_pairs`, `min_interval`, `max_interval` | Only if triggers 200–800ms apart |
| Position: after | `after = 99` | Only after marker trigger 99 |
| Position: before | `before = 88` | Only before marker trigger 88 |

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
condition = EegFun.EpochCondition(
    name = "Condition1", 
    trigger_sequences = [[1]]
)

# Extract epoch (-200ms to 1000ms around the trigger)
epochs = EegFun.extract_epochs(dat, condition, (-0.2, 1.0))
```

To match any of multiple trigger values for a single condition (OR logic), provide multiple sequences:

```julia
# Match trigger 1 OR trigger 2
condition = EegFun.EpochCondition(
    name = "Condition1", 
    trigger_sequences = [[1], [2]] # matches either trigger 1 or trigger 2
)

# Extract epoch (-200ms to 1000ms around the trigger)
epochs = EegFun.extract_epochs(dat, condition, (-0.2, 1.0))
```

## Multiple Conditions

Most experiments have more than one condition. Pass a vector of `EpochCondition` objects to `extract_epochs` to extract epochs for each condition separately:

```julia
# Define two separate conditions
conditions = [
    EegFun.EpochCondition(name = "Condition1", trigger_sequences = [[1]]),
    EegFun.EpochCondition(name = "Condition2", trigger_sequences = [[2]]),
]

# Extract epochs for both conditions
epochs = EegFun.extract_epochs(dat, conditions, (-0.2, 1.0))
```

Each condition is matched independently, and the resulting `EpochData` retains the condition labels for later analysis (e.g., averaging, comparison).

## Sequence Matching

Sometimes an "event" is actually a sequence of triggers. For example, a target stimulus (1) followed by a response (10).

```julia
# Match the sequence [1, 2]
condition = EegFun.EpochCondition(
    name = "TargetResponse", 
    trigger_sequences = [[1, 2]] # a 1 followed by a 2 with t=0 being the 1
)
```

### Wildcards and Ranges

You can use wildcards and ranges within your sequences:

- `:any`: Matches any trigger value.
- `UnitRange` (e.g., `1:10`): Matches any value within the range.

```julia
# Match trigger 1, then any trigger, then trigger 3
condition = EegFun.EpochCondition(
    name = "Wildcard", 
    trigger_sequences = [[1, :any, 3]]
)

# Match any trigger from 1 to 5, then trigger 10
condition = EegFun.EpochCondition(
    name = "Range", 
    trigger_sequences = [[1:5, 10]]
)
```

## Onset (t=0) Reference

By default, the first trigger in a sequence is considered $t=0$. You can change this using `reference_index`.

```julia
# Sequence: [Warning (1), Stimulus (2), Response (10)]
# We want t=0 to be the Stimulus (index 2)
condition = EegFun.EpochCondition(
    name = "StimulusOnset",
    trigger_sequences = [[1, 2, 10]], # 1 followed by 2 followed by 10 with t=0 being the 2
    reference_index = 2
)
```

::: tip
`reference_index` refers to the position within the trigger sequence, not the trigger value itself. The sequence matcher internally tracks the *actual sample position* for each element in the matched sequence, so `reference_index = 2` correctly resolves to the sample where trigger 2 occurred — even when zero-valued samples are skipped between triggers.
:::

## Timing Constraints

You can restrict matches to sequences where triggers occur within specific time intervals. This is useful for filtering out trials with late responses or accidental double-triggers.

```julia
# Match [1, 10] only if they occur between 200ms and 800ms apart
condition = EegFun.EpochCondition(
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
# Only find sequences that occur AFTER trigger 99 (e.g., start of a experimental block)
condition = EegFun.EpochCondition(
    name = "Block2",
    trigger_sequences = [[1, 2]],
    after = 99
)

# Only find sequences that occur BEFORE trigger 88 (e.g., end of first half)
condition = EegFun.EpochCondition(
    name = "Phase1",
    trigger_sequences = [[1, 2]],
    before = 88
)
```

> [!NOTE]
> You cannot specify both `after` and `before` on the same condition — use separate conditions if needed.

## External TOML Configuration

For complex studies, defining your epoch conditions in an external TOML file is often cleaner and is the recommended approach for the pipeline procedure.

### The TOML Format

Create a file (e.g., `epochs.toml`):

```toml
[epochs]

[[epochs.conditions]]
name = "condition_1"
trigger_sequences = [[1]]

[[epochs.conditions]]
name = "condition_2"
trigger_sequences = [[2]]

[[epochs.conditions]]
name = "condition_3" 
trigger_sequences = [[3]]

[[epochs.conditions]]
name = "condition_4"
trigger_sequences = [[4]]
```

Here is a more advanced example using sequence matching, timing constraints, and position constraints:

```toml
[epochs]

# Stimulus followed by a response, time-locked to the stimulus
[[epochs.conditions]]
name = "stimulus_response1"
trigger_sequences = [[1, 10]]
reference_index = 1
timing_pairs = [[1, 2]]
min_interval = 0.1
max_interval = 0.8

[[epochs.conditions]]
name = "stimulus_response2"
trigger_sequences = [[2, 10]]
reference_index = 1
timing_pairs = [[1, 2]]
min_interval = 0.1
max_interval = 0.8
```

### Loading and Using the TOML

```julia
using TOML

# Load the configuration
config = TOML.parsefile("epochs.toml")

# Parse into EpochCondition objects
conditions = EegFun.condition_parse_epoch(config)

# Extract epochs
epochs = EegFun.extract_epochs(dat, conditions, (-0.2, 1.0))
```

## Next Steps

- **Artifact handling** — [Artifact Handling](artifact-handling.md) for detecting and rejecting bad epochs after extraction
- **Batch pipelines** — [Batch Processing](batch-processing.md) uses epoch condition TOML files as part of its automated workflow
- **Selection predicates** — [Selection Patterns](selection-patterns.md) for filtering channels, samples, and time intervals in downstream analysis
