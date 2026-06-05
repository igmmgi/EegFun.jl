# FieldTrip Import

This demo demonstrates importing FieldTrip `.mat` files into EegFun.jl. FieldTrip is a MATLAB toolbox for EEG/MEG analysis, and its structures are a common way to exchange processed data.

### About FieldTrip .mat Format

FieldTrip stores data in MATLAB structures (generic `.mat` files). 
EegFun.jl uses `MAT.jl` to parse them and map them onto EegFun types.

**Key features**:

- Supports raw, epoched, and averaged (timelock) data

### Import Capabilities

**Data loading**:

- Import of continuous `raw` data
- Import of segmented `epochs` data
- Import of averaged `timelock` (ERP) data

### Data Mapping

**EegFun.read_fieldtrip** maps MATLAB structures to EegFun:

- FieldTrip `raw` / `epochs` → `ContinuousData` / `EpochData`
- FieldTrip `timelock` → `ErpData`
- Sample rates and time axes are automatically reconstructed

## Workflow Summary

1. **Load Layout**: Load a layout file and prepare it for plotting.
2. **Import Data**: Use `read_fieldtrip()` providing the path and the layout.
3. **Validation**: Visualize the imported structures using `plot_databrowser()` for continuous/epoched data or `plot_erp()` for averaged data.


## Code Examples

::: details Show Code

```julia
# FieldTrip Import Demo
#
# Demonstrates loading FieldTrip .mat files

# Note: EegFun.example_path() resolves bundled example data paths.
# When using your own data, simply pass the file path directly, e.g.:
# dat = EegFun.read_fieldtrip("/path/to/your/data.mat", layout)

using EegFun

# Load layout (FieldTrip doesn't store layout with data)
layout = EegFun.read_layout(EegFun.example_path("layouts/biosemi/biosemi72.csv"))
EegFun.polar_to_cartesian_xy!(layout)

# Load continuous data
println("Loading continuous data...")
continuous_data = EegFun.read_fieldtrip(EegFun.example_path("data/fieldtrip/continuous.mat"), layout)
EegFun.trigger_count(continuous_data)
EegFun.plot_databrowser(continuous_data)

# Load epoched data  
println("\nLoading epoched data...")
epoch_data = EegFun.read_fieldtrip(EegFun.example_path("data/fieldtrip/epochs.mat"), layout)
EegFun.plot_databrowser(epoch_data)

# Load ERP data
println("\nLoading ERP data...")
erp_data = EegFun.read_fieldtrip(EegFun.example_path("data/fieldtrip/erp.mat"), layout)

# Biologische Psychologie Labor Tübingen Custom mat files
# Essentially slightly stripped down FieldTrip structures
println("\nLoading ERP data...")
epoch_data = EegFun.read_fieldtrip(EegFun.example_path("data/fieldtrip/custom_epochs.mat"), layout)

# Load ERP data
println("\nLoading ERP data...")
erp_data = EegFun.read_fieldtrip(EegFun.example_path("data/fieldtrip/custom_erp.mat"), layout)
EegFun.plot_erp(erp_data, layout = :grid)



```

:::

## See Also

- [API Reference](../../reference/index.md)
