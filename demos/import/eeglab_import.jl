# Demo: Loading and Processing EEGLAB .set Files
# Note: trigger/event strings are hashed for :triggers, but available in :trigger_info.
# Once loaded, all EegFun functions work as expected.

# Note: EegFun.example_path() resolves bundled example data paths.
# When using your own data, simply pass the file path directly, e.g.:
# dat = EegFun.read_eeglab("/path/to/your/data.set")

using EegFun

# this seems to be a raw continuous data file without any ICA info
dat = EegFun.read_eeglab(EegFun.example_path("data/eeglab/eeglab_data.set"))
EegFun.plot_databrowser(dat)
EegFun.trigger_count(dat)

# this seems to be a epoched data file with ica info
dat, ica = EegFun.read_eeglab(EegFun.example_path("data/eeglab/epochs.set"))
EegFun.plot_databrowser(dat)

# We can plot the ICA activations 
EegFun.plot_topography(ica, component_selection = EegFun.components([1]))
