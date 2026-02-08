# Automated Batch Pipelines

Analyzing an entire cohort is the ultimate goal of most EEG experiments. While you can write your own processing loops, `EegFun.jl` is built around a powerful automated pipeline system that handles the complexity of file discovery, error isolation, and cohort-wide reporting for you.

## 1. Single-Line Execution

The core of the system is the `preprocess_v1` function. By passing it a single TOML configuration file, you trigger a multi-stage transformation across all your participants.

```julia
using EegFun

# Run the entire batch analysis
preprocess_v1("my_study_config.toml")
```

## 2. The Internal Lifecycle

When you run an automated pipeline, it doesn't just "apply a script" to files. It manages a sophisticated lifecycle consisting of several high-level phases:

### Phase 1: Environment Setup
The pipeline first loads your shared resources:
- Parses the **TOML configuration** and merges it with package defaults.
- Loads the **Electrode Layout** and calculates spatial neighbors for artifact repair.
- Discovers all **Raw Data Files** matching your selection patterns.

### Phase 2: Continuous-Level Processing
For each file, the pipeline performs initial "heavy lifting" on the continuous signal:
- **Rereferencing**: Standardizing the electrode voltages.
- **Filtering**: Removing slow drifts (High-pass) and high-frequency noise (Low-pass).
- **Artifact Detection**: Finding segments with extreme amplitudes or bad channel profiles.
- **ICA Component Removal**: Automatically identifying and removing eye blinks, muscle noise, and line interference.
- **Electrode Interpolation**: Repairing "bad" scalp channels using their spatial neighbors.

### Phase 3: Epoching & Cleanup
Once the continuous data is cleaned, the pipeline focuses on the specific experimental windows:
- **Epoch Extraction**: Cutting the data into segments based on triggers.
- **Baseline Correction**: Ensuring each segment starts from a zero-volt average.
- **Per-Epoch Repair**: Interpolating electrodes that only became "bad" during specific segments of the task.
- **Epoch Rejection**: Final removal of any segments that remain too noisy to analyze.

### Phase 4: Cohort Consolidation
After the last participant is processed, the pipeline automatically generates summaries of the entire study:
- **Grand Average Counts**: A report on how many trials survived for each participant/condition.
- **Electrode Reliability**: A summary identifying which electrodes were frequently repaired across the cohort.
- **ICA Statistics**: Reports on the average number and type of components removed per subject.

## 3. Resilience and Isolation

The pipeline is built to be "fail-safe." 

1.  **Logical Isolation**: If a specific participant's file is corrupted or missing, the pipeline logs the error but continues processing the rest of the cohort.
2.  **Log Granularity**: You receive a **Study-level log** summarizing the entire run, plus **Individual-level logs** for every participant, containing every table and warning generated during their specific processing.
3.  **Traceability**: Each output file is marked with metadata about the pipeline version and parameters used to create it.

---

> [!TIP]
> You can find the reference implementation for this entire flow in `src/pipelines/pipeline_v1.jl`.
