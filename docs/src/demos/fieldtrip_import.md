# FieldTrip Import

This demo demonstrates importing FieldTrip `.mat` files into EegFun.jl. FieldTrip is a leading MATLAB toolbox for EEG/MEG analysis, and its structures are a common way to exchange processed data.

### About FieldTrip .mat Format

FieldTrip stores data in MATLAB structures (typically named `data`, `data_all`, etc.). Since these are generic `.mat` files, EegFun.jl uses `MAT.jl` to parse them and map them onto EegFun types.

**Key features**:
- Supports raw, epoched, and averaged (timelock) data
- Can handle complex data structures with custom fields
- Widely used in academic research for advanced analysis (TF, connectivity)

### Import Capabilities

**Data loading**:
- Import of continuous `raw` data
- Import of segmented `epochs` data
- Import of averaged `timelock` (ERP) data
- Support for custom Tübingen lab variants of FieldTrip structures

### Data Mapping

**EegFun.read_fieldtrip** maps MATLAB structures to EegFun:
- FieldTrip `raw` / `epochs` → `ContinuousData` / `EpochedData`
- FieldTrip `timelock` → `ERPData`
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

using EegFun

# Load layout (FieldTrip doesn't store layout with data)
layout = EegFun.read_layout("./resources/layouts/biosemi/biosemi72.csv")
EegFun.polar_to_cartesian_xy!(layout)

# Load continuous data
println("Loading continuous data...")
continuous_data = EegFun.read_fieldtrip("./resources/data/fieldtrip/continuous.mat", layout)
EegFun.trigger_count(continuous_data)
EegFun.plot_databrowser(continuous_data)

# Load epoched data  
println("\nLoading epoched data...")
epoch_data = EegFun.read_fieldtrip("./resources/data/fieldtrip/epochs.mat", layout)
EegFun.plot_databrowser(epoch_data)

# Load ERP data
println("\nLoading ERP data...")
erp_data = EegFun.read_fieldtrip("./resources/data/fieldtrip/erp.mat", layout)

# Biologische Psychologie Labor Tübingen Custom mat files
# Essentially slightly strippeepd down FieldTrip structures
println("\nLoading ERP data...")
epoch_data = EegFun.read_fieldtrip("./resources/data/fieldtrip/custom_epochs.mat", layout)

# Load ERP data
println("\nLoading ERP data...")
erp_data = EegFun.read_fieldtrip("./resources/data/fieldtrip/custom_erp.mat", layout)
EegFun.plot_erp(erp_data, layout = :grid)



```

:::

## See Also

- [API Reference](../reference/index.md)
