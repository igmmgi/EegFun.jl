# Script to automatically generate demo documentation from demo scripts
using Printf

# Demo metadata: (filename_without_ext, title, description)
demos = [
    # Artifact Detection
    ("plot_artifacts", "Artifact Detection Visualization", "Demonstrates visualization of detected artifacts in epoched data."),
    ("plot_joint_probability", "Joint Probability Plot", "Visualizes channel joint probability for quality assessment."),
    (
        "artifacts",
        "Artifact Detection",
        "Comprehensive artifact detection workflow including EOG detection, extreme value detection, bad epoch detection, artifact repair, and rejection.",
    ),

    # ICA
    (
        "ica",
        "ICA Analysis",
        "Complete ICA workflow including component identification, artifact removal, and visualization for continuous and epoched data.",
    ),

    # Statistics
    (
        "erp_measurements",
        "ERP Measurements",
        "Extract quantitative features from ERP data including amplitude, latency, area, and peak-to-peak measurements.",
    ),
    (
        "statistics",
        "Statistical Analysis",
        "Statistical analysis options for ERP data including analytic t-tests and cluster-based permutation tests.",
    ),
    (
        "decoding",
        "MVPA Decoding",
        "Multivariate pattern analysis (MVPA) decoding workflow with SVM classification and statistical testing.",
    ),
    ("rsa", "Representational Similarity Analysis", "RSA workflow for analyzing representational geometries across conditions and time."),

    # Time-Frequency
    ("tf_morlet", "Morlet Wavelet Analysis", "Time-frequency analysis using Morlet wavelets with various cycle configurations."),
    ("tf_multitaper", "Multitaper Spectrograms", "Time-frequency analysis using multitaper method."),
    ("tf_stft", "Short-Time Fourier Transform", "Time-frequency analysis using STFT with configurable window parameters."),

    # Visualization
    ("plot_channel_spectrum", "Channel Frequency Spectrum", "Visualize frequency spectra for channels."),
    ("plot_channel_summary", "Channel Summary Plot", "Plot summary statistics for all channels."),
    ("plot_correlation_heatmap", "Correlation Heatmap", "Visualize channel correlation matrices."),
    ("plot_databrowser", "Interactive Data Browser", "Interactive data browser for continuous, epoch, and ICA data."),
    ("plot_epochs", "Epoch Visualization", "Visualizes individual epochs with channel selection."),
    ("plot_erp", "ERP Waveform Plot", "Plots event-related potentials with customizable options."),
    ("plot_erp_image", "ERP Image Plot", "Creates ERP image plots showing trial-by-trial variations."),
    ("plot_filter", "Filter Response Visualization", "Visualizes filter frequency and phase responses."),
    ("plot_layout", "Electrode Layout", "Displays electrode layout configurations."),
    ("plot_topography", "Topographic Maps", "Creates topographic scalp maps at specific time points."),
    ("plot_triggers", "Trigger Marker Visualization", "Visualizes event markers and triggers in continuous data."),
]

# Create demo markdown files
for (filename, title, description) in demos
    source_file = "demos/$(filename).jl"
    output_file = "docs/src/demos/$(replace(filename, "_" => "-")).md"

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

    ```julia
    $source_code```

    ## See Also

    - [API Reference](../reference/index.md)
    """

    # Write output file
    write(output_file, markdown_content)
    @info "Created: $output_file"
end

@info "Demo documentation generation complete!"
