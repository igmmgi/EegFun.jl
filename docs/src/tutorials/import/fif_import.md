# FIF Import

## Overview

Demonstrates FIF Import functionality.


## Code Examples

::: details Show Code

```julia
# Demo: Loading and Processing FIF Files
# Shows how to load Functional Image Format (.fif) files,
# create EegFun data structures, apply basic preprocessing,
# and visualize the data.

# Note: EegFun.example_path() resolves bundled example data paths.
# When using your own data, simply pass the file path directly, e.g.:
# dat = EegFun.read_raw_data("/path/to/your/data.fif")

using EegFun
using GLMakie

# Load raw FIF data
raw_fif = EegFun.read_raw_data(EegFun.example_path("data/fif/test_raw.fif"))

# Using the fallback auto-generated layout (if you don't have a specific .csv layout)
dat = EegFun.create_eegfun_data(raw_fif)

EegFun.all_data(dat)
EegFun.meta_data(dat)
EegFun.all_labels(dat)
EegFun.channel_labels(dat)
EegFun.meta_labels(dat)
EegFun.extra_labels(dat) # empty


# Basic preprocessing
EegFun.highpass_filter!(dat, 1.0) # 1 Hz high-pass filter
EegFun.rereference!(dat, :avg)    # Average reference

println("\nTrigger summary:")
EegFun.trigger_count(dat)

# Open interactive databrowser
EegFun.plot_databrowser(dat);
```

:::

## See Also

- [API Reference](../../reference/index.md)
