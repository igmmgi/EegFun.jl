using EegFun

# Read raw data
dat = EegFun.read_raw_data("./resources/data/example1.bdf");

# Read and prepare layout file
layout_file = EegFun.read_layout("./resources/layouts/biosemi/biosemi72.csv");
EegFun.polar_to_cartesian_xy!(layout_file)

# Create EegFun data structure (EegFun.ContinuousData)
dat = EegFun.create_eeg_dataframe(dat, layout_file);

# Minimal preprocessing
EegFun.rereference!(dat, :avg)
EegFun.highpass_filter!(dat, 0.1)
EegFun.baseline!(dat)

# Define epoch conditions
epoch_cfg = [
    EegFun.EpochCondition(name = "Trigger1", trigger_sequences = [[1]]),
    EegFun.EpochCondition(name = "Trigger2", trigger_sequences = [[2]]),
]

# Extract epochs (-200ms to 1000ms around triggers)
epochs = EegFun.extract_epochs(dat, epoch_cfg, (-0.2, 1.0));

# Plot individual epochs for first condition
EegFun.plot_epochs(epochs[1], channel_selection = EegFun.channels([:Fp1, :Cz, :Oz]))

# Plot different channel selections
EegFun.plot_epochs(epochs[1], channel_selection = EegFun.channels([:Fp1]))
EegFun.plot_epochs(epochs[2], channel_selection = EegFun.channels([:Cz]))

# Baseline correction to stimulus onset (t=0)
EegFun.baseline!(epochs, (-0.2, 0.0))
EegFun.plot_epochs(epochs[1], channel_selection = EegFun.channels([:Cz]))

# Compare both conditions on same plot
erps = EegFun.average_epochs(epochs)
EegFun.plot_erp(erps, channel_selection = EegFun.channels([:Cz]))

# Demonstrate epoch rejection based on amplitude
# (assuming reject_epochs function exists)
# epochs_clean = EegFun.reject_epochs(epochs, threshold = 100.0)

# Show trial counts
println("Condition 1: $(length(epochs[1].epochs)) trials")
println("Condition 2: $(length(epochs[2].epochs)) trials")
