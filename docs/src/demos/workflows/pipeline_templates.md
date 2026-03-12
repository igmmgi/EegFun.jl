# Pipeline Templates

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


## Code Examples

::: details Show Code

```julia
# Demo: Pipeline Templates
# Shows how to generate scaffold files for custom preprocessing pipelines
# and configuration templates.

using EegFun

#######################################################################
# GENERATE A PIPELINE TEMPLATE
#######################################################################

# Generate a basic pipeline template with default settings
# This creates a .jl file with all the boilerplate: logging, config loading,
# error handling, and a structured loop for custom processing steps
EegFun.generate_pipeline_template("my_custom_pipeline.jl")

# Generate with a custom function name and more steps
EegFun.generate_pipeline_template(
    "my_pipeline.jl",
    "my_preprocess";
    options = EegFun.PipelineTemplateOptions(num_steps = 8, subsections_per_step = 3),
)


#######################################################################
# GENERATE A CONFIG TEMPLATE
#######################################################################

# Generate a TOML configuration template with all available settings
EegFun.generate_config_template(filename = "my_config.toml")

# The generated config includes sections for:
# - files (input/output paths)
# - filter (highpass, lowpass, ICA-specific filters)
# - eog (EOG channel configuration and thresholds)
# - eeg (artifact and extreme value thresholds)
# - ica (whether to apply, percentage of data)
# - layout (neighbour criterion)
# - epoch timing (start/end times)

# You can read and inspect a config file
# cfg = EegFun.read_config("my_config.toml")
# EegFun.print_config(cfg)


#######################################################################
# CLEANUP
#######################################################################

# Remove generated template files (for demo purposes)
isfile("my_custom_pipeline.jl") && rm("my_custom_pipeline.jl")
isfile("my_pipeline.jl") && rm("my_pipeline.jl")
isfile("my_config.toml") && rm("my_config.toml")
```

:::

## See Also

- [API Reference](../../reference/index.md)
