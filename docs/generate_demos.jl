# Script to automatically generate demo documentation from demo scripts
using Printf

# Demo metadata: (filename_without_ext, title)
demos = [
    ("artifacts", "Artifacts"),
    ("baseline", "Baseline"),
    ("biosemi_import", "BioSemi Import"),
    ("brainvision_import", "BrainVision Import"),
    ("channel_metrics", "Channel Metrics"),
    ("channel_repair", "Channel Repair"),
    ("channel_summary", "Channel Summary"),
    ("data", "Data"),
    ("decoding", "Decoding"),
    ("erp_measurements", "ERP Measurements"),
    ("ica", "ICA"),
    ("mirror", "Mirror"),
    ("plot_artifacts", "Plot Artifacts"),
    ("plot_channel_spectrum", "Plot Channel Spectrum"),
    ("plot_channel_summary", "Plot Channel Summary"),
    ("plot_correlation_heatmap", "Plot Correlation Heatmap"),
    ("plot_databrowser", "Plot Databrowser"),
    ("plot_epochs", "Plot Epochs"),
    ("plot_erp", "Plot ERP"),
    ("plot_erp_image", "Plot ERP Image"),
    ("plot_filter", "Plot Filter"),
    ("plot_joint_probability", "Plot Joint Probability"),
    ("plot_layout", "Plot Layout"),
    ("plot_topography", "Plot Topography"),
    ("plot_triggers", "Plot Triggers"),
    ("rereference", "Rereference"),
    ("resample", "Resample"),
    ("rsa", "RSA"),
    ("statistics", "Statistics"),
    ("tf_morlet", "TF Morlet"),
    ("tf_multitaper", "TF Multitaper"),
    ("tf_stft", "TF STFT"),
]

# Create demo markdown files
for demo_info in demos
    filename, title = demo_info

    source_file = "demos/$(filename).jl"
    output_file = "docs/src/demos/$(filename).md"
    overview_file = "docs/src/demos/overviews/$(filename).md"

    if !isfile(source_file)
        @warn "Source file not found: $source_file"
        continue
    end

    # Read source code
    source_code = read(source_file, String)

    # Check for custom overview file
    if isfile(overview_file)
        @info "Using custom overview from: $overview_file"
        overview_content = read(overview_file, String)

        # Extract first paragraph as description (up to first blank line)
        lines = split(overview_content, '\n')
        description_lines = []
        for line in lines
            if isempty(strip(line))
                break
            end
            push!(description_lines, line)
        end
        description = join(description_lines, ' ')

        # Use full overview content
        overview_section = overview_content
    else
        # Fallback: use default overview and description
        @warn "No custom overview found for $filename, using default"
        description = "Demonstrates $title functionality."
        overview_section = """
        ## Overview

        $description
        """
    end

    # Create markdown content
    markdown_content = """
    # $title

    $description

    $overview_section

    ## Code Examples

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
