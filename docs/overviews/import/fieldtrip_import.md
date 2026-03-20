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
