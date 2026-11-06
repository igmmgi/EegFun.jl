# Artifact & Rejection GUIs

This demo walks through EegFun's interactive GUIs for visual epoch rejection and before/after artifact inspection.

## Included Interactive Tools

| Function | Purpose |
| --- | --- |
| `detect_bad_epochs_interactive` | Paginated grid of single-trial epochs with per-epoch green/red rejection toggles |
| `plot_artifact_detection` | Interactive epoch-by-epoch and artifact-by-artifact waveform browser |
| `plot_artifact_repair` | Before/after comparison of raw vs. interpolated epochs |
| `plot_artifact_rejection` | Before/after comparison highlighting rejected epochs with red axis spines |

## Typical Workflow

1. Run automatic detection (`detect_bad_epochs_automatic`) to pre-flag candidate trials, or start from scratch.
2. Review and toggle trials visually with `detect_bad_epochs_interactive`, then extract indices via `get_rejected_epochs` or drop flagged trials directly with `reject_epochs`.
3. Inspect repair (`repair_artifacts`) or rejection results interactively using the `◀ Previous Artifact` / `Next Artifact ▶` navigation buttons.



## Code Examples

::: details Show Code

```julia
# Demo: Interactive Artifact & Epoch Rejection GUIs
# Short walkthrough for testing EegFun's interactive epoch rejection grid
# and before/after artifact inspection viewers.

using EegFun
using GLMakie

#######################################################################
# PREPARE EXAMPLE EPOCHED DATA & AUTOMATIC DETECTIONS
#######################################################################

raw = EegFun.read_raw_data(EegFun.example_path("data/bdf/example1.bdf"))
layout = EegFun.read_layout(EegFun.example_path("layouts/biosemi/biosemi72.csv"))
EegFun.polar_to_cartesian_xy!(layout)
dat = EegFun.create_eegfun_data(raw, layout)

EegFun.rereference!(dat, :avg)
EegFun.highpass_filter!(dat, 0.1)

epoch_cfg = [
    EegFun.EpochCondition(name = "Trigger1", trigger_sequences = [[1]]),
    EegFun.EpochCondition(name = "Trigger2", trigger_sequences = [[2]]),
]
epochs = EegFun.extract_epochs(dat, epoch_cfg, (-0.2, 1.0))
EegFun.baseline!(epochs)

bad_epochs = EegFun.detect_bad_epochs_automatic(epochs, z_criterion = 3, abs_criterion = 100)

#######################################################################
# 1. INTERACTIVE EPOCH REJECTION GRID (detect_bad_epochs_interactive)
#######################################################################

# Controls to test:
# - Per-epoch green/red Toggle switches (top-left of each subplot)
# - "◀ Prev" / "Next ▶" page navigation buttons
# - Bad-channels filter toggle & global Y-limit toggle
state = EegFun.detect_bad_epochs_interactive(epochs[1], artifact_info = bad_epochs[1], dims = (4, 4))
# state = EegFun.detect_bad_epochs_interactive(epochs, artifact_info = bad_epochs, dims = (4, 4))

# Purely manual review (starts with zero epochs flagged):
state_manual = EegFun.detect_bad_epochs_interactive(epochs[1], dims = (4, 4))
# state_manual = EegFun.detect_bad_epochs_interactive(epochs, dims = (4, 4))

# Retrieve flagged epoch indices or apply interactive rejections to data:
rejected_idx = EegFun.get_rejected_epochs(state)
epochs_clean = EegFun.reject_epochs(epochs[1], state)

#######################################################################
# 2. ARTIFACT DETECTION VIEWER (plot_artifact_detection)
#######################################################################

# Controls to test:
# - "◀ Previous" / "Next ▶" (step epoch-by-epoch)
# - "◀ Previous Artifact" / "Next Artifact ▶" (jump directly between flagged epochs)
EegFun.plot_artifact_detection(epochs[1], bad_epochs[1])

#######################################################################
# 3. ARTIFACT REPAIR VIEWER (plot_artifact_repair)
#######################################################################

# Interactive before/after comparison of original (top) vs repaired (bottom) epochs
epochs_repaired = EegFun.repair_artifacts(epochs, bad_epochs, method = :neighbor_interpolation)
EegFun.plot_artifact_repair(epochs[1], epochs_repaired[1], bad_epochs[1])

#######################################################################
# 4. ARTIFACT REJECTION VIEWER (plot_artifact_rejection)
#######################################################################

# Interactive before/after comparison; rejected epochs show red axis spines in the lower panel
epochs_rejected = EegFun.reject_epochs(epochs, bad_epochs)
EegFun.plot_artifact_rejection(epochs[1], epochs_rejected[1], bad_epochs[1])
```

:::

## See Also

- [API Reference](../../reference/index.md)
