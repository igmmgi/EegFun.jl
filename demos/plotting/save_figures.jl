# Demo: Saving Figures
# Shows how to save Makie figures in different file formats (PNG, SVG, PDF)
# with resolution and quality settings.
using EegFun
# Note: EegFun.example_path() resolves bundled example data paths.
# When using your own data, simply pass the file path directly, e.g.:
# dat = EegFun.read_raw_data("/path/to/your/data.bdf")
using GLMakie
using CairoMakie  # loaded for backend = CairoMakie in save(), not activated

const DEMO_OUTPUT = "./demos/output/"
mkpath(DEMO_OUTPUT)

# ── Prepare some data ───────────────────────────────────────────────────────
dat = EegFun.read_raw_data(EegFun.example_path("data/bdf/example1.bdf"));
layout = EegFun.read_layout(EegFun.example_path("layouts/biosemi/biosemi72.csv"));
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
