# Script to automatically generate demo documentation from demo scripts
using Printf

# Demo metadata: (filename_without_ext, title, description)
demos = [
    (
        "artifacts",
        "Artifacts",
        "Comprehensive artifact detection workflow including EOG detection, extreme value detection, bad epoch detection, artifact repair, and rejection.",
    ),
    ("baseline", "Baseline", "Baseline correction methods for ERP data."),
    ("channel_metrics", "Channel Metrics", "Calculate and visualize channel quality metrics."),
    ("channel_repair", "Channel Repair", "Repair bad channels using interpolation methods."),
    ("channel_summary", "Channel Summary", "Generate summary statistics for channels."),
    ("data", "Data", "Data loading and basic data structures."),
    ("decoding", "Decoding", "Multivariate pattern analysis (MVPA) decoding workflow with SVM classification and statistical testing."),
    (
        "erp_measurements",
        "ERP Measurements",
        "Extract quantitative features from ERP data including amplitude, latency, area, and peak-to-peak measurements.",
    ),
    (
        "ica",
        "ICA",
        "Complete ICA workflow including component identification, artifact removal, and visualization for continuous and epoched data.",
    ),
    ("mirror", "Mirror", "Mirror electrode positions for lateralized analyses."),
    ("plot_artifacts", "Plot Artifacts", "Visualization of detected artifacts in epoched data."),
    ("plot_channel_spectrum", "Plot Channel Spectrum", "Visualize frequency spectra for channels."),
    ("plot_channel_summary", "Plot Channel Summary", "Plot summary statistics for all channels."),
    ("plot_correlation_heatmap", "Plot Correlation Heatmap", "Visualize channel correlation matrices."),
    ("plot_databrowser", "Plot Databrowser", "Interactive data browser for continuous, epoch, and ICA data."),
    ("plot_epochs", "Plot Epochs", "Visualizes individual epochs with channel selection."),
    ("plot_erp", "Plot ERP", "Plots event-related potentials with customizable options."),
    ("plot_erp_image", "Plot ERP Image", "Creates ERP image plots showing trial-by-trial variations."),
    ("plot_filter", "Plot Filter", "Visualizes filter frequency and phase responses."),
    ("plot_joint_probability", "Plot Joint Probability", "Visualizes channel joint probability for quality assessment."),
    ("plot_layout", "Plot Layout", "Displays electrode layout configurations."),
    ("plot_topography", "Plot Topography", "Creates topographic scalp maps at specific time points."),
    ("plot_triggers", "Plot Triggers", "Visualizes event markers and triggers in continuous data."),
    ("rereference", "Rereference", "Re-reference EEG data to different reference schemes."),
    ("resample", "Resample", "Resample data to different sampling rates."),
    ("rsa", "RSA", "RSA workflow for analyzing representational geometries across conditions and time."),
    (
        "statistics",
        "Statistics",
        "Statistical analysis options for ERP data including analytic t-tests and cluster-based permutation tests.",
    ),
    ("tf_morlet", "TF Morlet", "Time-frequency analysis using Morlet wavelets with various cycle configurations."),
    ("tf_multitaper", "TF Multitaper", "Time-frequency analysis using multitaper method."),
    ("tf_stft", "TF STFT", "Time-frequency analysis using STFT with configurable window parameters."),
]

# Create demo markdown files
for (filename, title, description) in demos
    source_file = "demos/$(filename).jl"
    output_file = "docs/src/demos/$(filename).md"

    if !isfile(source_file)
        @warn "Source file not found: $source_file"
        continue
    end

    # Read source code
    source_code = read(source_file, String)

    # Create markdown content
    markdown_content = """
    # $title

    $description

    ## Overview

    Demonstrates $description

    ## Source Code

    ::: details Show Code
    ```julia
    $source_code```
    :::

    ## See Also

    - [API Reference](../reference/index.md)
    """

    # Write output file
    write(output_file, markdown_content)
    @info "Created: $output_file"
end

@info "Demo documentation generation complete!"
