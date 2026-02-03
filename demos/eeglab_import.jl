"""
Demo: Loading and Processing EEGLAB .set Files

This demo shows how to:
1. Load EEGLAB .set files into EegFun
2. Perform basic ERP analysis
3. Apply preprocessing
4. Visualize results
"""

using EegFun

# read raw data
dat = EegFun.load_eeglab("./resources/data/eeglab/epochs.set");
ica = EegFun.load_eeglab_ica("./resources/data/eeglab/epochs.set");

# TODO: coordinates look slightly off!

epochs = EegFun.load_eeglab("./resources/data/eeglab/epochs.set")
ica = EegFun.load_eeglab_ica("./resources/data/eeglab/epochs.set")

# Compare first channel
println("Epoch layout: ", epochs.layout.data[1, :])
println("ICA layout:   ", ica.layout.data[1, :])


EegFun.plot_topography(ica, component_selection = EegFun.components([1]))

EegFun.plot_databrowser(dat)
# EegFun.plot_epochs(dat, layout = :grid)




