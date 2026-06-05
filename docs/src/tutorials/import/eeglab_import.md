# EEGLAB Import

This demo demonstrates importing EEGLAB `.set` files into EegFun.jl. EEGLAB is a widely used open-source toolboxes for EEG processing (in MATLAB). 

### About EEGLAB .set Format

The `.set` format is a MATLAB-based file that contains a header structure with all recording parameters, trigger information, and often precomputed info like ICA weights or epoch definitions. The actual EEG data may be stored within the `.set` file or in a separate `.fdt` file.

**Key features**:

- Supports both continuous and epoched data
- Often includes ICA components and weights
- Comprehensive metadata storage

### Import Capabilities

**Data loading**:

- Automatic detection of continuous vs. epoched data
- Import of ICA weights and sphere matrices, if available
- Mapping of event/trigger labels (hashed to triggers)
- Support for external data files (.fdt)

### Data Mapping

**EegFun.read_eeglab** handles mapping of EEGLAB structures to native EegFun types:

- EEGLAB Dataset → `ContinuousData`, `EpochData`, or `ErpData`
- ICA info → `ICA` structure
- Event labels → Available in `:trigger_info` column

## Workflow Summary

1. **Load Data**: Depending on the file content, `read_eeglab()` can return just data or data plus ICA information.
2. **Check Triggers**: EEGLAB often uses string labels for events. EegFun hashes these for its trigger system while preserving the original labels.
3. **Verification**: Use `plot_databrowser()` to verify the imported time series or epochs.
4. **ICA Visualization**: If ICA information was imported, you can immediately plot components and activations.


## Code Examples

::: details Show Code

```julia
# Demo: Loading and Processing EEGLAB .set Files
# Note: trigger/event strings are hashed for :triggers, but available in :trigger_info.
# Once loaded, all EegFun functions work as expected.

# Note: EegFun.example_path() resolves bundled example data paths.
# When using your own data, simply pass the file path directly, e.g.:
# dat = EegFun.read_eeglab("/path/to/your/data.set")

using EegFun

# this seems to be a raw continuous data file without any ICA info
dat = EegFun.read_eeglab(EegFun.example_path("data/eeglab/eeglab_data.set"))
EegFun.plot_databrowser(dat)
EegFun.trigger_count(dat)

# this seems to be a epoched data file with ica info
dat, ica = EegFun.read_eeglab(EegFun.example_path("data/eeglab/epochs.set"))
EegFun.plot_databrowser(dat)

# We can plot the ICA activations 
EegFun.plot_topography(ica, component_selection = EegFun.components([1]))
```

:::

## See Also

- [API Reference](../../reference/index.md)
