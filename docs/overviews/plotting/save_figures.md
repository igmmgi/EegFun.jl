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
