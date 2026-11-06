This demo walks through EegFun's interactive GUIs for visual epoch rejection and before/after artifact inspection.

### Included Interactive Tools

| Function | Purpose |
| --- | --- |
| `detect_bad_epochs_interactive` | Paginated grid of single-trial epochs with per-epoch green/red rejection toggles |
| `plot_artifact_detection` | Interactive epoch-by-epoch and artifact-by-artifact waveform browser |
| `plot_artifact_repair` | Before/after comparison of raw vs. interpolated epochs |
| `plot_artifact_rejection` | Before/after comparison highlighting rejected epochs with red axis spines |

### Typical Workflow

1. Run automatic detection (`detect_bad_epochs_automatic`) to pre-flag candidate trials, or start from scratch.
2. Review and toggle trials visually with `detect_bad_epochs_interactive`, then extract indices via `get_rejected_epochs` or drop flagged trials directly with `reject_epochs`.
3. Inspect repair (`repair_artifacts`) or rejection results interactively using the `◀ Previous Artifact` / `Next Artifact ▶` navigation buttons.

