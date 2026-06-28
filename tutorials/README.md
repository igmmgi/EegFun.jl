# EegFun.jl Tutorials

This directory contains Julia source scripts (`.jl`) for EegFun.jl tutorials and demos.

## How it works

These scripts serve two purposes:

1. **Direct use** — You can run any script directly in Julia to learn EegFun.jl functionality
2. **Documentation source** — The `docs/generate_demos.jl` script converts these `.jl` files into markdown documentation pages in `docs/src/tutorials/`

## Directory structure

| Folder | Contents |
|---|---|
| `artifacts/` | Artifact detection tutorials |
| `data/` | Data access and selection helpers |
| `erp/` | ERP measurements, GFP, grand average, LRP, etc. |
| `import/` | Importing data from various formats (BioSemi, BrainVision, EDF, etc.) |
| `plotting/` | All plotting tutorials and figure saving |
| `preprocessing/` | Filtering, rereferencing, epoching, ICA, etc. |
| `statistics/` | Statistical analysis, decoding, RSA |
| `time_frequency/` | Time-frequency analysis (Morlet, multitaper, STFT) |
| `workflows/` | Batch processing, pipeline templates, preprocessing workflows |

## Custom overviews

Each tutorial page can have a custom overview section. These are stored in `docs/overviews/` (matching the same subfolder structure) and are automatically injected by `generate_demos.jl` when building the documentation.

## Building documentation

To regenerate the tutorial documentation pages:

```bash
julia --project=docs docs/generate_demos.jl
```

Or use the interactive documentation manager:

```bash
julia --project=docs docs/doc_manager.jl
```
