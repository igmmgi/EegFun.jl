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



