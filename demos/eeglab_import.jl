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
