This demo walks through a complete manual preprocessing pipeline from raw BioSemi data to cleaned ERPs. It mirrors the steps and order used in `EegFun.preprocess_v1`, making it a useful reference for understanding what the automated pipeline does or for building a custom workflow.

### When to Use This

- **Learning** — understand each preprocessing step and why it is applied in a specific order
- **Custom workflows** — adapt the pipeline to your specific experimental design
- **Debugging** — step through the pipeline interactively to diagnose issues with your data

### Processing Phases

The workflow follows three phases:

### Phase 1: Setup and Initial Preprocessing

- Load raw data and configure electrode layout with coordinates and neighbours
- Mark epoch intervals for targeted analysis
- Rereference to average
- Apply initial filters (0.1 Hz high-pass for ERPs)
- Calculate and detect EOG (vEOG, hEOG)
- Extract "original" epochs for later comparison

### Phase 2: Artifact Detection and Cleaning (Continuous Level)

- Compute channel summary statistics and joint probability
- Detect extreme values at two thresholds (250 μV for ICA exclusion, 75 μV for artifact rejection)
- Identify and partition bad channels (EOG-related vs. non-EOG-related)
- Run ICA on a separate 1 Hz high-pass filtered copy, identify and remove artifact components
- Repair bad channels via neighbour interpolation
- Recalculate EOG channels after ICA and repair

### Phase 3: Epoch Extraction and Epoch-Level Processing

- Extract epochs from cleaned continuous data
- Baseline correction
- Two rounds of automatic epoch detection with channel repair between rounds
- Reject remaining bad epochs
- Average to ERPs and compare against originals

### Key Takeaway

The demo produces both "original" (lightly preprocessed) and "good" (fully cleaned) ERPs so you can visually assess the effect of the entire pipeline.
