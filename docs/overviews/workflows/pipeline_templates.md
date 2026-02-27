This demo shows how to generate scaffold files for custom preprocessing pipelines and configuration files.

### Pipeline Templates

`generate_pipeline_template` creates a ready-to-use Julia file with:

- Structured preprocessing loop over participants
- Built-in logging setup and error handling
- Numbered processing steps with placeholder comments
- Configuration loading from TOML files

This saves time when starting a new preprocessing project by providing a well-organised starting point.

### Configuration Templates

`generate_config_template` creates a TOML configuration file with all available preprocessing parameters:

- **Files**: Input/output directory paths
- **Filter**: High-pass, low-pass, and ICA-specific filter settings
- **EOG**: EOG channel names and correlation thresholds
- **EEG**: Artifact and extreme value detection thresholds
- **ICA**: Whether to apply ICA and data percentage settings
- **Layout**: Neighbour distance criterion
- **Epoch**: Timing parameters (start/end times)

### Using Templates with Pipelines

The generated pipeline and config files are designed to work together:

1. Generate both template files
2. Edit the config TOML with your study-specific parameters
3. Customise the pipeline Julia file with your processing steps
4. Run the pipeline

## Workflow Summary

This demo covers:

- Generating a pipeline template file
- Customising template options (number of steps, subsections)
- Generating a configuration template
- Understanding the generated config sections
