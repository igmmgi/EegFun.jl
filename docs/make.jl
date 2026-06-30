using Documenter
using DocumenterVitepress

# Add the parent directory to the load path so we can load the local package
push!(LOAD_PATH, dirname(@__DIR__))
using EegFun


pages = [
    "Home" => "index.md",
    "Getting Started" => [
        "tutorials/getting-started.md",
        "IDE Workflows" => "tutorials/ide-workflows.md",
        "Why Julia?" => "tutorials/why-julia.md",
        "Julia Basics" => "tutorials/julia-basics.md",
    ],
    "Tutorials" => [
        "Overview" => "tutorials/index.md",
        "Data Structures" => "explanations/data-structures.md",
        "Layouts" => "explanations/layouts.md",
        "Manual Preprocessing" => "tutorials/manual-preprocessing.md",
        "Epoch Selection" => "tutorials/epoch-selection.md",
        "Artifact Handling" => "tutorials/artifact-handling.md",
        "Layouts and Neighbors" => "tutorials/layouts-neighbors.md",
        "Import" => [
            "Data Import" => "tutorials/import/data_import.md",
            "BioSemi Import" => "tutorials/import/biosemi_import.md",
            "BrainVision Import" => "tutorials/import/brainvision_import.md",
            "EDF Import" => "tutorials/import/edf_import.md",
            "FIF Import" => "tutorials/import/fif_import.md",
            "XDF Import" => "tutorials/import/xdf_import.md",
            "EEGLAB Import" => "tutorials/import/eeglab_import.md",
            "FieldTrip Import" => "tutorials/import/fieldtrip_import.md",
            "Data" => "tutorials/import/data.md",
            "Data Persistence (JLD2)" => "tutorials/import/jld2.md",
            "BIDS Export" => "tutorials/import/bids_export.md",
        ],
        "Data" => [
            "Data" => "tutorials/data/data.md",
            "Data Access" => "tutorials/data/data_access.md",
            "Selection Helpers" => "tutorials/data/selection_helpers.md",
        ],
        "Preprocessing" => [
            "Filter" => "tutorials/preprocessing/filter.md",
            "Rereference" => "tutorials/preprocessing/rereference.md",
            "Resample" => "tutorials/preprocessing/resample.md",
            "Baseline" => "tutorials/preprocessing/baseline.md",
            "Mirror" => "tutorials/preprocessing/mirror.md",
            "Triggers" => "tutorials/preprocessing/triggers.md",
            "Epoch Extraction" => "tutorials/preprocessing/epochs.md",
            "Channel Operations" => "tutorials/preprocessing/channel_operations.md",
            "Channel Metrics" => "tutorials/preprocessing/channel_metrics.md",
            "Channel Repair" => "tutorials/preprocessing/channel_repair.md",
            "Channel Summary" => "tutorials/preprocessing/channel_summary.md",
            "ICA" => "tutorials/preprocessing/ica.md",
            "Layouts & Neighbours" => "tutorials/preprocessing/layouts.md",
            "Analysis Settings" => "tutorials/preprocessing/analysis_settings.md",
        ],
        "Artifacts" => [
            "Artifacts" => "tutorials/artifacts/artifacts.md",
            "Artifact Detection" => "tutorials/artifacts/artifact_detection.md",
        ],
        "ERP Analysis" => [
            "ERP Measurements" => "tutorials/erp/erp_measurements.md",
            "Global Field Power" => "tutorials/erp/gfp.md",
            "Grand Average" => "tutorials/erp/grand_average.md",
            "Jackknife Average" => "tutorials/erp/jackknife_average.md",
            "Lateralised Readiness Potential" => "tutorials/erp/lrp.md",
            "Realignment" => "tutorials/erp/realign.md",
            "Condition Operations" => "tutorials/erp/condition_operations.md",
        ],
        "Time-Frequency" => [
            "TF Analysis" => "tutorials/time_frequency/tf_analysis.md",
            "TF Morlet" => "tutorials/time_frequency/tf_morlet.md",
            "TF Multitaper" => "tutorials/time_frequency/tf_multitaper.md",
            "TF STFT" => "tutorials/time_frequency/tf_stft.md",
            "TF Operations" => "tutorials/time_frequency/tf_operations.md",
        ],
        "Statistics" => [
            "Statistics" => "tutorials/statistics/statistics.md",
            "TF Statistics" => "tutorials/statistics/tf_stats_test.md",
            "Decoding" => "tutorials/statistics/decoding.md",
            "RSA" => "tutorials/statistics/rsa.md",
        ],
        "Plotting" => [
            "Plot Artifacts" => "tutorials/plotting/plot_artifacts.md",
            "Plot Channel Spectrum" => "tutorials/plotting/plot_channel_spectrum.md",
            "Plot Channel Summary" => "tutorials/plotting/plot_channel_summary.md",
            "Plot Correlation Heatmap" => "tutorials/plotting/plot_correlation_heatmap.md",
            "Plot Epochs" => "tutorials/plotting/plot_epochs.md",
            "Plot ERP" => "tutorials/plotting/plot_erp.md",
            "Plot ERP Image" => "tutorials/plotting/plot_erp_image.md",
            "Plot ERP Measurements" => "tutorials/plotting/plot_erp_measurements.md",
            "Plot Filter" => "tutorials/plotting/plot_filter.md",
            "Plot Frequency Spectrum" => "tutorials/plotting/plot_frequency_spectrum.md",
            "Plot GFP" => "tutorials/plotting/plot_gfp.md",
            "Plot ICA" => "tutorials/plotting/plot_ica.md",
            "Plot Joint Probability" => "tutorials/plotting/plot_joint_probability.md",
            "Plot Layout" => "tutorials/plotting/plot_layout.md",
            "Plot Topography" => "tutorials/plotting/plot_topography.md",
            "Plot Triggers" => "tutorials/plotting/plot_triggers.md",
            "Saving Figures" => "tutorials/plotting/save_figures.md",
        ],
        "Interactive / GUI" => [
            "Plot Databrowser" => "tutorials/plotting/plot_databrowser.md",
            "Plot ERP Filter GUI" => "tutorials/plotting/plot_erp_filter_gui.md",
            "Plot ERP Measurement GUI" => "tutorials/plotting/plot_erp_measurement_gui.md",
        ],
        "Specialized Plotting" => [
            "Plot Decoding" => "tutorials/plotting/plot_decoding.md",
            "Plot RSA" => "tutorials/plotting/plot_rsa.md",
            "Plot Statistics" => "tutorials/plotting/plot_statistics.md",
            "Plot Time-Frequency" => "tutorials/plotting/plot_tf.md",
        ],
    ],
    "How-to Guides" => [
        "Filter Data" => "how-to/filter-data.md",
        "Selection Patterns" => "tutorials/selection-patterns.md",
        "TOML Format" => "tutorials/toml-format.md",
        "Batch Processing" => "tutorials/batch-processing.md",
        "Pipeline Templates" => "tutorials/workflows/pipeline_templates.md",
        "Preprocessing Walkthrough" => "tutorials/workflows/preprocessing_workflow.md",
        "Plot GUI" => "tutorials/plot-gui.md",
        "Benchmark" => "tutorials/benchmark.md",
        "Worked Examples" => [
            "Visual Attention (Posner Cueing)" => "tutorials/experiments/visual-attention.md",
            "N170 (Face/Body)" => "tutorials/experiments/n170.md",
        ],
    ],
    "Explanations" => [
        "Data Structures" => "explanations/data-structures.md",
        "Layouts" => "explanations/layouts.md",
        "Artifact Handling" => "tutorials/artifact-handling.md",
        "Manual Preprocessing" => "tutorials/manual-preprocessing.md",
    ],
    "Teaching Demos" => [
        "Signal Processing" => [
            "Sampling" => "tutorials/teaching/signal_processing/signal_example_sampling.md",
            "Composition" => "tutorials/teaching/signal_processing/signal_example_composition.md",
            "Dot Product" => "tutorials/teaching/signal_processing/signal_example_dotproduct.md",
            "Convolution" => "tutorials/teaching/signal_processing/signal_example_convolution.md",
            "Time-Frequency" => "tutorials/teaching/signal_processing/signal_example_tf.md",
        ],
        "ICA" => [
            "Introduction" => "tutorials/teaching/ica/signal_example_ica_0.md",
            "Mixture" => "tutorials/teaching/ica/signal_example_ica_mixture.md",
            "Separation" => "tutorials/teaching/ica/signal_example_ica_separation.md",
            "Geometry" => "tutorials/teaching/ica/signal_example_ica_geometry.md",
            "Sphering" => "tutorials/teaching/ica/signal_example_ica_sphering.md",
            "Optimization" => "tutorials/teaching/ica/signal_example_ica_optimization.md",
            "CLT" => "tutorials/teaching/ica/signal_example_ica_clt.md",
            "Infomax" => "tutorials/teaching/ica/signal_example_ica_infomax.md",
        ],
        "Machine Learning" => [
            "Decoding" => "tutorials/teaching/machine_learning/signal_example_decoding.md",
        ],
        "ERP" => [
            "Simulate ERP" => "tutorials/teaching/erp/simulate_erp.md",
        ],
    ],
    "Cheatsheet" => "cheatsheet.md",
    "Reference" => [
        "Overview" => "reference/index.md",
        "Analysis" => "reference/analysis.md",
        "Plotting" => "reference/plotting.md",
        "Types" => "reference/types.md",
    ],
]


# Check if we are building the PDF or the VitePress site
build_pdf = get(ENV, "BUILD_PDF", "false") == "true"

if build_pdf
    using DocumenterTypst
    build_format = DocumenterTypst.Typst()
    build_dir = "build_pdf"
    
    # Create a temporary source directory for the PDF build
    src_dir = "src_pdf"
    ispath(joinpath(@__DIR__, src_dir)) && rm(joinpath(@__DIR__, src_dir), force=true, recursive=true)
    cp(joinpath(@__DIR__, "src"), joinpath(@__DIR__, src_dir))
    
    # Copy public assets into assets/public so Documenter natively copies them to the build directory
    mkpath(joinpath(@__DIR__, src_dir, "assets", "public"))
    for item in readdir(joinpath(@__DIR__, "src", "public"))
        cp(joinpath(@__DIR__, "src", "public", item), joinpath(@__DIR__, src_dir, "assets", "public", item), force=true)
    end
    
    # Fix absolute image paths for Typst
    for (root, dirs, files) in walkdir(joinpath(@__DIR__, src_dir))
        for file in files
            if endswith(file, ".md")
                path = joinpath(root, file)
                content = read(path, String)
                # Prevent Typst syntax error "expected function, found content" for *17*(4) by stripping italics
                content = replace(content, r"\*([0-9]+)\*\s*\(" => s"\1(")
                # Remove VitePress custom container markers like `::: details Show Code` or `:::`
                content = replace(content, r"^:::\s*.*?\n"m => "")
                # Remove https images because Typst cannot download them without network flags
                content = replace(content, r"!\[([^\]]*)\]\(https?://[^\)]+\)" => "")
                # Replace ![alt](/path) with ![alt](/assets/public/path)
                content = replace(content, r"!\[([^\]]*)\]\((/[^\)]+)\)" => s"![\1](/assets/public\2)")
                write(path, content)
            end
        end
    end
    
    # Filter out duplicate pages to prevent Typst label errors
    # Explanations contains duplicate files already in Tutorials
    filter!(p -> p[1] != "Explanations", pages)
    
    # Remove Layouts from Tutorials because it is hundreds of images
    for section in pages
        if section[1] == "Tutorials"
            filter!(p -> p[1] != "Layouts", section[2])
        end
    end
else
    build_format = DocumenterVitepress.MarkdownVitepress(repo = "https://github.com/igmmgi/EegFun.jl", devbranch = "main")
    build_dir = "build"
    src_dir = "src"
end

try
    makedocs(;
        modules = [EegFun],
        authors = "igmmgi",
        sitename = "EegFun",
        repo = "https://github.com/igmmgi/EegFun.jl",
        format = build_format,
        warnonly = [:linkcheck, :cross_references, :missing_docs],
        draft = false,
        source = src_dir,
        build = build_dir,
        pages = pages,
    )
finally
    if build_pdf && ispath(joinpath(@__DIR__, "src_pdf"))
        rm(joinpath(@__DIR__, "src_pdf"), force=true, recursive=true)
    end
end

# Deploy built VitePress site (DocumenterVitepress.deploydocs required since v0.2 for correct base paths)
if !build_pdf
    DocumenterVitepress.deploydocs(repo = "github.com/igmmgi/EegFun.jl", devbranch = "main", push_preview = true)
end
