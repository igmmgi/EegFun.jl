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
EegFun.generate_config_template("my_config.toml")

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
