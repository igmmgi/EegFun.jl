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

**EegFun.read_eeglab** handles the mapping EEGLAB structures to native EegFun types:
- EEGLAB Dataset → `ContinuousData` or `EpochedData`
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
"""
Demo: Loading and Processing EEGLAB .set Files

This demo shows how to:

TODO: This can be considered work in progress. I am not familiar with eeglab *.set/*.fdt files, so 
the code here is a bit of a guesswork. But it does seem to work with the two example datasets 
I found in eeglab/sample_data.

NB. trigger/event strings are hashed for the :triggers column, but are available in the :trigger_info column

Once the data is loaded, all EegFun functions should work as expected as 
read_eeglab converts to EegFun types.
"""

using EegFun

# this seems to be a raw continuous data file without and ICA info
dat = EegFun.read_eeglab("./resources/data/eeglab/eeglab_data.set")
EegFun.plot_databrowser(dat)
EegFun.trigger_count(dat)

# this seems to be a epoched data file with ica info
dat, ica = EegFun.read_eeglab("./resources/data/eeglab/epochs.set")
EegFun.plot_databrowser(dat)

# We can plot the ICA activations 
EegFun.plot_topography(ica, component_selection = EegFun.components([1]))
```

:::

## See Also

- [API Reference](../../reference/index.md)
