This demo demonstrates importing FieldTrip `.mat` files into EegFun.jl, including continuous, epoched, and ERP data.

### What is FieldTrip Format?

FieldTrip is a MATLAB toolbox for MEG/EEG analysis. The `.mat` format stores:

- EEG/MEG data structures
- Trial definitions and timing
- Channel information
- Time vectors
- Event markers

### Import Capabilities

**Continuous data**:

- Raw time series from continuous recording
- Trial structure with consistent time vectors
- Event markers and triggers

**Epoched data**:

- Segmented trials around events
- Trial-specific timing information
- Flexible trial structure (different trial lengths supported)

**ERP data**:

- Averaged event-related potentials
- Single time series per condition
- Smaller file sizes

### Data Mapping

**EegFun.read_fieldtrip** converts FieldTrip structures to native EegFun types:

- FieldTrip continuous → `ContinuousData`
- FieldTrip epoched → `EpochData`
- FieldTrip timelock (ERP) → `ErpData`

### Layout Requirement

FieldTrip files do not always include electrode layout information. You must provide a layout file separately when importing FieldTrip data.

### Custom FieldTrip Variants

This demo also shows support for custom FieldTrip-like formats (e.g., our own simplified lab-specific variants). These follow similar structures but may have reduced metadata.

## Important Notes

FieldTrip import is functional but considered work-in-progress. The implementation has only been tested with a few basic FieldTrip datatypes (continuous, epoched, and timelock/ERP structures). More complex FieldTrip features (e.g., frequency data, source reconstructions) are not currently supported. 

## Workflow Summary

This demo shows FieldTrip import workflows:

### 1. Load Layout

- FieldTrip doesn't store layout with data
- Load electrode layout separately
- Convert to 2D coordinates

### 2. Import Continuous Data

- Load continuous FieldTrip structure
- Check trigger counts
- Visualize in databrowser

### 3. Import Epoched Data

- Load segmented trial data
- Visualize epochs in databrowser

### 4. Import ERP Data

- Load averaged ERP data
- Plot ERPs with various layouts

### 5. Custom Variants

- Import lab-specific FieldTrip-like formats
- Same workflow as standard FieldTrip
