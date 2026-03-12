# Saving Figures

This demo shows how to save Makie figures in different file formats with resolution and quality settings.

### Backends

| Backend | Purpose |
|---------|---------|
| **GLMakie** | Interactive work (windows, zooming, panning) |
| **CairoMakie** | High-quality export (raster and vector) |

You only need `using CairoMakie` — you do **not** need to activate it. Pass `backend = CairoMakie` in the `save()` call instead.

### File Formats

| Format | Type | Key Option |
|--------|------|------------|
| **PNG** | Raster | `px_per_unit` — scale resolution (e.g. `4` for print quality) |
| **SVG** | Vector | `pt_per_unit` — scale canvas size in points |
| **PDF** | Vector | `pt_per_unit` — scale canvas size in points |

### PNG (Raster)

```julia
save("figure.png", fig)                    # default resolution
save("figure.png", fig; px_per_unit = 4)   # high resolution
save("figure.png", fig; update = false)    # save current interactive view as-is
```

### SVG & PDF (Vector)

Always pass `backend = CairoMakie` when saving vector formats. Without it, `save()` uses the active backend (GLMakie), which embeds a rasterized bitmap — not true vector graphics.

```julia
save("figure.svg", fig; backend = CairoMakie)
save("figure.pdf", fig; backend = CairoMakie)
save("figure.pdf", fig; backend = CairoMakie, pt_per_unit = 1.5)  # larger
```

### Quick Reference

```julia
using GLMakie      # interactive work
using CairoMakie   # loaded but not activated

# ... create your plot ...
save("raster.png", fig; px_per_unit = 4)                          # high-res PNG
save("vector.svg", fig; backend = CairoMakie)                     # true vector SVG
save("vector.pdf", fig; backend = CairoMakie, pt_per_unit = 1.5)  # scaled PDF
```


## Code Examples

::: details Show Code

```julia
# Demo: Saving Figures
# Shows how to save Makie figures in different file formats (PNG, SVG, PDF)
# with resolution and quality settings.
using EegFun
using GLMakie
using CairoMakie  # loaded for backend = CairoMakie in save(), not activated

const DEMO_OUTPUT = "./demos/output/"
mkpath(DEMO_OUTPUT)

# ── Prepare some data ───────────────────────────────────────────────────────
dat = EegFun.read_raw_data("./resources/data/bdf/example1.bdf");
layout = EegFun.read_layout("./resources/layouts/biosemi/biosemi72.csv");
EegFun.polar_to_cartesian_xy!(layout)
dat = EegFun.create_eegfun_data(dat, layout);
EegFun.rereference!(dat, :avg)
EegFun.highpass_filter!(dat, 0.1)

epoch_cfg = [
    EegFun.EpochCondition(name = "ExampleEpoch1", trigger_sequences = [[1]]),
    EegFun.EpochCondition(name = "ExampleEpoch2", trigger_sequences = [[2]]),
]
epochs = EegFun.extract_epochs(dat, epoch_cfg, (-2, 4))
erps = EegFun.average_epochs(epochs)

# Create a figure (GLMakie is the active backend → interactive window)
GLMakie.activate!()
result = EegFun.plot_erp(erps, channel_selection = EegFun.channels([:Cz, :PO7, :PO8]), average_channels = true)
fig = result.fig


# ╔══════════════════════════════════════════════════════════════════════════╗
# ║  PNG (raster)                                                            ║
# ╚══════════════════════════════════════════════════════════════════════════╝

# Default PNG — figure size (in units) determines pixel dimensions
save(joinpath(DEMO_OUTPUT, "erp_default.png"), fig)

# Higher resolution — px_per_unit scales every unit by this factor.
# px_per_unit = 4 gives ~4× the pixel count → great for posters / print.
save(joinpath(DEMO_OUTPUT, "erp_hires.png"), fig; px_per_unit = 4)

# Save the current interactive view exactly as-is (skip layout recalculation).
# Useful after you have interactively zoomed / panned in the GLMakie window.
save(joinpath(DEMO_OUTPUT, "erp_current_view.png"), fig; update = false)


# ╔══════════════════════════════════════════════════════════════════════════╗
# ║  SVG & PDF (vector)                                                      ║
# ╚══════════════════════════════════════════════════════════════════════════╝
# IMPORTANT: always pass `backend = CairoMakie` when saving vector formats.
# Without it, save() uses the currently active backend (GLMakie), which
# embeds a rasterized bitmap inside the SVG/PDF — not true vector graphics!

# --- SVG ---
# Ideal for web or documents where you need infinitely scalable graphics.
save(joinpath(DEMO_OUTPUT, "erp.svg"), fig; backend = CairoMakie)

# pt_per_unit controls the overall size of the SVG/PDF canvas in points.
# Larger values → bigger document; smaller → more compact.
save(joinpath(DEMO_OUTPUT, "erp_large.svg"), fig; backend = CairoMakie, pt_per_unit = 1.5)

# --- PDF ---
# The standard for journal submissions and LaTeX documents.
save(joinpath(DEMO_OUTPUT, "erp.pdf"), fig; backend = CairoMakie)

# Scaled-up PDF — useful when the default appears too small in a document.
save(joinpath(DEMO_OUTPUT, "erp_large.pdf"), fig; backend = CairoMakie, pt_per_unit = 1.5)

GLMakie.closeall()
```

:::

## See Also

- [API Reference](../../reference/index.md)
