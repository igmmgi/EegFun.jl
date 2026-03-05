using Documenter
using DocumenterVitepress

# Add the parent directory to the load path so we can load the local package
push!(LOAD_PATH, dirname(@__DIR__))
using EegFun


pages = [
    "Home" => "index.md",
    "Tutorials" => [
        "Overview" => "tutorials/index.md",
        "Getting Started" => ["tutorials/getting-started.md", "IDE Workflows" => "tutorials/ide-workflows.md"],
        "Manual Preprocessing" => "tutorials/manual-preprocessing.md",
        "Epoch Selection" => "tutorials/epoch-selection.md",
        "Artifact Handling" => "tutorials/artifact-handling.md",
        "Layouts and Neighbors" => "tutorials/layouts-neighbors.md",
        "Batch Processing" => "tutorials/batch-processing.md",
        "Selection Patterns" => "tutorials/selection-patterns.md",
        "Plot GUI" => "tutorials/plot-gui.md",
        "Preprocessing Walkthrough" => "demos/workflows/preprocessing_workflow.md",
        "Worked Examples" => [
            "Visual Attention (Posner Cueing)" => "demos/experiments/visual-attention.md",
            "N170 (Face/Body)" => "demos/experiments/n170.md",
        ],
    ],
    "Explanations" => [
        "Overview" => "explanations/index.md",
        "Data Structures" => "explanations/data-structures.md",
        "Layouts" => "explanations/layouts.md",
    ],
    "Demos" => [
        "Import" => [
            "BioSemi Import" => "demos/import/biosemi_import.md",
            "BrainVision Import" => "demos/import/brainvision_import.md",
            "EEGLAB Import" => "demos/import/eeglab_import.md",
            "FieldTrip Import" => "demos/import/fieldtrip_import.md",
            "Data" => "demos/import/data.md",
            "Data Persistence (JLD2)" => "demos/import/jld2.md",
        ],
        "Preprocessing" => [
            "Filter" => "demos/preprocessing/filter.md",
            "Rereference" => "demos/preprocessing/rereference.md",
            "Resample" => "demos/preprocessing/resample.md",
            "Baseline" => "demos/preprocessing/baseline.md",
            "Mirror" => "demos/preprocessing/mirror.md",
            "Triggers" => "demos/preprocessing/triggers.md",
            "Epoch Extraction" => "demos/preprocessing/epochs.md",
            "Channel Operations" => "demos/preprocessing/channel_operations.md",
            "Channel Metrics" => "demos/preprocessing/channel_metrics.md",
            "Channel Repair" => "demos/preprocessing/channel_repair.md",
            "Channel Summary" => "demos/preprocessing/channel_summary.md",
            "ICA" => "demos/preprocessing/ica.md",
        ],
        "Artifacts" => ["Artifacts" => "demos/artifacts/artifacts.md", "Artifact Detection" => "demos/artifacts/artifact_detection.md"],
        "ERP Analysis" => [
            "ERP Measurements" => "demos/erp/erp_measurements.md",
            "Global Field Power" => "demos/erp/gfp.md",
            "Grand Average" => "demos/erp/grand_average.md",
            "Jackknife Average" => "demos/erp/jackknife_average.md",
            "Lateralised Readiness Potential" => "demos/erp/lrp.md",
            "Realignment" => "demos/erp/realign.md",
            "Condition Operations" => "demos/erp/condition_operations.md",
        ],
        "Plotting" => [
            "Plot Artifacts" => "demos/plotting/plot_artifacts.md",
            "Plot Channel Spectrum" => "demos/plotting/plot_channel_spectrum.md",
            "Plot Channel Summary" => "demos/plotting/plot_channel_summary.md",
            "Plot Correlation Heatmap" => "demos/plotting/plot_correlation_heatmap.md",
            "Plot Epochs" => "demos/plotting/plot_epochs.md",
            "Plot ERP" => "demos/plotting/plot_erp.md",
            "Plot ERP Image" => "demos/plotting/plot_erp_image.md",
            "Plot ERP Measurements" => "demos/plotting/plot_erp_measurements.md",
            "Plot Filter" => "demos/plotting/plot_filter.md",
            "Plot Frequency Spectrum" => "demos/plotting/plot_frequency_spectrum.md",
            "Plot GFP" => "demos/plotting/plot_gfp.md",
            "Plot ICA" => "demos/plotting/plot_ica.md",
            "Plot Joint Probability" => "demos/plotting/plot_joint_probability.md",
            "Plot Layout" => "demos/plotting/plot_layout.md",
            "Plot Topography" => "demos/plotting/plot_topography.md",
            "Plot Triggers" => "demos/plotting/plot_triggers.md",
        ],
        "Interactive / GUI" => [
            "Plot Databrowser" => "demos/plotting/plot_databrowser.md",
            "Plot ERP Filter GUI" => "demos/plotting/plot_erp_filter_gui.md",
            "Plot ERP Measurement GUI" => "demos/plotting/plot_erp_measurement_gui.md",
        ],
        "Specialized Plotting" => [
            "Plot Decoding" => "demos/plotting/plot_decoding.md",
            "Plot RSA" => "demos/plotting/plot_rsa.md",
            "Plot Statistics" => "demos/plotting/plot_statistics.md",
            "Plot Time-Frequency" => "demos/plotting/plot_tf.md",
        ],
        "Time-Frequency" => [
            "TF Morlet" => "demos/time_frequency/tf_morlet.md",
            "TF Multitaper" => "demos/time_frequency/tf_multitaper.md",
            "TF STFT" => "demos/time_frequency/tf_stft.md",
            "TF Operations" => "demos/time_frequency/tf_operations.md",
        ],
        "Statistics" => [
            "Statistics" => "demos/statistics/statistics.md",
            "TF Statistics" => "demos/statistics/tf_stats_test.md",
            "Decoding" => "demos/statistics/decoding.md",
            "RSA" => "demos/statistics/rsa.md",
        ],
    ],
    "Cheatsheet" => "cheatsheet.md",
    "Reference" => ["Overview" => "reference/index.md", "Types" => "reference/types.md"],
]

makedocs(;
    modules = [EegFun],
    authors = "igmmgi",
    sitename = "EegFun",
    repo = "https://github.com/igmmgi/EegFun.jl",
    format = DocumenterVitepress.MarkdownVitepress(repo = "https://github.com/igmmgi/EegFun.jl", devbranch = "main"),
    warnonly = [:linkcheck, :cross_references, :missing_docs],
    draft = false,
    source = "src",
    build = "build",
    pages = pages,
)

# Deploy built VitePress site (DocumenterVitepress.deploydocs required since v0.2 for correct base paths)
DocumenterVitepress.deploydocs(repo = "github.com/igmmgi/EegFun.jl", devbranch = "main", push_preview = true)
