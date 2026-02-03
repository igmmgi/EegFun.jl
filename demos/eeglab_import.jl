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

EegFun.plot_databrowser(dat)
EegFun.plot_epochs(dat, layout = :grid)




