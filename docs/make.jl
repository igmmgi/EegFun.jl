using Documenter
using DocumenterVitepress

# Add the parent directory to the load path so we can load the local package
push!(LOAD_PATH, dirname(@__DIR__))
using EegFun

pages = [
    "Home" => "index.md",
    "Tutorials" => [
        "Getting Started" => "tutorials/getting-started.md",
        "Basic Preprocessing" => "tutorials/basic-preprocessing.md",
        "ERP Analysis" => "tutorials/erp-analysis.md",
        "ICA Workflow" => "tutorials/ica-workflow.md",
    ],
    "How-To Guides" => [
        "Filter Data" => "how-to/filter-data.md",
        "Create Epochs" => "how-to/create-epochs.md",
        "Topographic Plots" => "how-to/topographic-plots.md",
    ],
    "Explanations" => ["Data Structures" => "explanations/data-structures.md", "Statistical Methods" => "explanations/statistics.md"],
    "Demos" => [
        "Artifacts" => "demos/artifacts.md",
        "Baseline" => "demos/baseline.md",
        "Channel Metrics" => "demos/channel_metrics.md",
        "Channel Repair" => "demos/channel_repair.md",
        "Channel Summary" => "demos/channel_summary.md",
        "Data" => "demos/data.md",
        "Decoding" => "demos/decoding.md",
        "ERP Measurements" => "demos/erp_measurements.md",
        "ICA" => "demos/ica.md",
        "Mirror" => "demos/mirror.md",
        "Plot Artifacts" => "demos/plot_artifacts.md",
        "Plot Channel Spectrum" => "demos/plot_channel_spectrum.md",
        "Plot Channel Summary" => "demos/plot_channel_summary.md",
        "Plot Correlation Heatmap" => "demos/plot_correlation_heatmap.md",
        "Plot Databrowser" => "demos/plot_databrowser.md",
        "Plot Epochs" => "demos/plot_epochs.md",
        "Plot ERP" => "demos/plot_erp.md",
        "Plot ERP Image" => "demos/plot_erp_image.md",
        "Plot Filter" => "demos/plot_filter.md",
        "Plot Joint Probability" => "demos/plot_joint_probability.md",
        "Plot Layout" => "demos/plot_layout.md",
        "Plot Topography" => "demos/plot_topography.md",
        "Plot Triggers" => "demos/plot_triggers.md",
        "Rereference" => "demos/rereference.md",
        "Resample" => "demos/resample.md",
        "RSA" => "demos/rsa.md",
        "Statistics" => "demos/statistics.md",
        "TF Morlet" => "demos/tf_morlet.md",
        "TF Multitaper" => "demos/tf_multitaper.md",
        "TF STFT" => "demos/tf_stft.md",
    ],
    "Reference" => [
        "Overview" => "reference/index.md",
        "Types" => "reference/types.md",
        # TODO: Add more reference pages as needed
    ],
]

makedocs(;
    modules = [EegFun],
    authors = "igmmgi",
    sitename = "EegFun",
    repo = "https://github.com/igmmgi/EegFun.jl",
    format = DocumenterVitepress.MarkdownVitepress(repo = "https://github.com/igmmgi/EegFun.jl", devbranch = "main", devurl = "dev"),
    warnonly = [:linkcheck, :cross_references, :missing_docs],
    draft = false,
    source = "src",
    build = "build",
    pages = pages,
)

# Deploy built VitePress site (DocumenterVitepress.deploydocs required since v0.2 for correct base paths)
DocumenterVitepress.deploydocs(repo = "github.com/igmmgi/EegFun.jl", devbranch = "main", push_preview = true)
