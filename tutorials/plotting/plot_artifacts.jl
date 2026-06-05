# Demo: Artifact Detection Visualization
# Shows plots for visualizing artifact detection and rejection.

# Note: EegFun.example_path() resolves bundled example data paths.
# When using your own data, simply pass the file path directly, e.g.:
# dat = EegFun.read_raw_data("/path/to/your/data.bdf")

using EegFun

# read raw data
dat = EegFun.read_raw_data(EegFun.example_path("data/bdf/example1.bdf"));
dat = EegFun.create_eegfun_data(dat)

# minimal preprocessing
EegFun.rereference!(dat, :avg)
EegFun.highpass_filter!(dat, 1)

# extract epochs
epoch_cfg = EegFun.EpochCondition(name = "ExampleEpoch1", trigger_sequences = [[1]])
epochs = EegFun.extract_epochs(dat, epoch_cfg, (-2, 4))

# Artifacts
artifacts = EegFun.detect_bad_epochs_automatic(epochs)

# Plot artifacts
EegFun.plot_artifact_detection(epochs, artifacts)

# Or with custom criteria
artifacts = EegFun.detect_bad_epochs_automatic(epochs, abs_criterion = 100, z_criterion = 0)
EegFun.plot_artifact_detection(epochs, artifacts)

artifacts = EegFun.detect_bad_epochs_automatic(epochs, abs_criterion = 0, z_criterion = 3)
EegFun.plot_artifact_detection(epochs, artifacts)

# extract multiple epochs
epoch_cfg = [
    EegFun.EpochCondition(name = "ExampleEpoch1", trigger_sequences = [[1]]),
    EegFun.EpochCondition(name = "ExampleEpoch2", trigger_sequences = [[2]]),
]
epochs = EegFun.extract_epochs(dat, epoch_cfg, (-2, 4))

# Artifacts
artifacts = EegFun.detect_bad_epochs_automatic(epochs)

# Plot artifacts
EegFun.plot_artifact_detection(epochs, artifacts)
EegFun.plot_artifact_detection(epochs[1], artifacts[1])


# Interactive artifact detection
EegFun.detect_bad_epochs_interactive(epochs, dims = (4, 4))

# Or combine with automatic artifact detection
artifacts = EegFun.detect_bad_epochs_automatic(epochs)
EegFun.detect_bad_epochs_interactive(epochs, artifact_info = artifacts, dims = (4, 4))

# Or combine with automatic artifact detection
EegFun.detect_bad_epochs_interactive(epochs, artifact_info = artifacts, dims = (4, 4), colormap = :GnBu)




