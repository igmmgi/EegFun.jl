# Script to automatically generate demo documentation from demo scripts
using Printf

# Demo metadata: (subfolder, filename_without_ext, title)
demos = [
    # Import
    ("import", "biosemi_import", "BioSemi Import"),
    ("import", "brainvision_import", "BrainVision Import"),
    ("import", "eeglab_import", "EEGLAB Import"),
    ("import", "fieldtrip_import", "FieldTrip Import"),
    ("import", "jld2", "Data Persistence (JLD2)"),

    # Data
    ("data", "data", "Data"),
    ("data", "data_access", "Data Access"),
    ("data", "selection_helpers", "Selection Helpers"),

    # Preprocessing
    ("preprocessing", "filter", "Filter"),
    ("preprocessing", "rereference", "Rereference"),
    ("preprocessing", "resample", "Resample"),
    ("preprocessing", "baseline", "Baseline"),
    ("preprocessing", "mirror", "Mirror"),
    ("preprocessing", "triggers", "Triggers"),
    ("preprocessing", "epochs", "Epoch Extraction"),
    ("preprocessing", "channel_operations", "Channel Operations"),
    ("preprocessing", "channel_metrics", "Channel Metrics"),
    ("preprocessing", "channel_repair", "Channel Repair"),
    ("preprocessing", "channel_summary", "Channel Summary"),
    ("preprocessing", "ica", "ICA"),
    ("preprocessing", "layouts", "Layouts & Neighbours"),
    ("preprocessing", "analysis_settings", "Analysis Settings"),

    # Artifacts
    ("artifacts", "artifact_detection", "Artifact Detection"),

    # ERP Analysis
    ("erp", "erp_measurements", "ERP Measurements"),
    ("erp", "gfp", "Global Field Power"),
    ("erp", "grand_average", "Grand Average"),
    ("erp", "jackknife_average", "Jackknife Average"),
    ("erp", "lrp", "Lateralised Readiness Potential"),
    ("erp", "realign", "Realignment"),
    ("erp", "condition_operations", "Condition Operations"),

    # Statistics
    ("statistics", "statistics", "Statistics"),
    ("statistics", "tf_stats_test", "TF Statistics"),
    ("statistics", "decoding", "Decoding"),
    ("statistics", "rsa", "RSA"),

    # Time-Frequency
    ("time_frequency", "tf_morlet", "TF Morlet"),
    ("time_frequency", "tf_multitaper", "TF Multitaper"),
    ("time_frequency", "tf_stft", "TF STFT"),
    ("time_frequency", "tf_operations", "TF Operations"),

    # Plotting
    ("plotting", "plot_artifacts", "Plot Artifacts"),
    ("plotting", "plot_channel_spectrum", "Plot Channel Spectrum"),
    ("plotting", "plot_channel_summary", "Plot Channel Summary"),
    ("plotting", "plot_correlation_heatmap", "Plot Correlation Heatmap"),
    ("plotting", "plot_databrowser", "Plot Databrowser"),
    ("plotting", "plot_decoding", "Plot Decoding"),
    ("plotting", "plot_epochs", "Plot Epochs"),
    ("plotting", "plot_erp", "Plot ERP"),
    ("plotting", "plot_erp_filter_gui", "Plot ERP Filter GUI"),
    ("plotting", "plot_erp_image", "Plot ERP Image"),
    ("plotting", "plot_erp_measurements", "Plot ERP Measurements"),
    ("plotting", "plot_erp_measurement_gui", "Plot ERP Measurement GUI"),
    ("plotting", "plot_filter", "Plot Filter"),
    ("plotting", "plot_frequency_spectrum", "Plot Frequency Spectrum"),
    ("plotting", "plot_gfp", "Plot GFP"),
    ("plotting", "plot_ica", "Plot ICA"),
    ("plotting", "plot_joint_probability", "Plot Joint Probability"),
    ("plotting", "plot_layout", "Plot Layout"),
    ("plotting", "plot_rsa", "Plot RSA"),
    ("plotting", "plot_statistics", "Plot Statistics"),
    ("plotting", "plot_tf", "Plot Time-Frequency"),
    ("plotting", "plot_topography", "Plot Topography"),
    ("plotting", "plot_triggers", "Plot Triggers"),

    # Workflows
    ("workflows", "preprocessing_workflow", "Preprocessing Workflow"),
    ("workflows", "batch_processing", "Batch Processing"),
    ("workflows", "pipeline_templates", "Pipeline Templates"),
]

# Create demo markdown files
for demo_info in demos
    subfolder, filename, title = demo_info

    source_file = "demos/$(subfolder)/$(filename).jl"
    output_file = "docs/src/demos/$(subfolder)/$(filename).md"
    overview_file = "docs/overviews/$(subfolder)/$(filename).md"

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

    $overview_section

    ## Code Examples

    ::: details Show Code

    ```julia
    $source_code```

    :::

    ## See Also

    - [API Reference](../../reference/index.md)
    """

    # Write output file
    write(output_file, markdown_content)
    @info "Created: $output_file"
end

@info "Demo documentation generation complete!"
