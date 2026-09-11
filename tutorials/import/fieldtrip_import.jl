# FieldTrip Import Demo
#
# Demonstrates loading FieldTrip .mat files

# Note: EegFun.example_path() resolves bundled example data paths.
# When using your own data, simply pass the file path directly, e.g.:
# dat = EegFun.read_fieldtrip("/path/to/your/data.mat", layout)

using EegFun
using GLMakie

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



