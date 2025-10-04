# manual testing of detect_bad_epochs
using eegfun
using GLMakie

# Get some basic data with initial preprocessing steps (high-pass filter, epoch)
data_file = joinpath(@__DIR__, "..", "..", "..",  "Flank_C_3.bdf")
layout_file = eegfun.read_layout("./data/layouts/biosemi/biosemi72.csv");

dat = eegfun.read_bdf(data_file);
dat = eegfun.create_eeg_dataframe(dat, layout_file);

eegfun.filter_data!(dat, "hp", 1)
epoch_cfg = [eegfun.EpochCondition(name = "ExampleEpoch1", trigger_sequences = [[1]])]
epochs = eegfun.EpochData[]
for (idx, epoch) in enumerate(epoch_cfg)
    push!(epochs, eegfun.extract_epochs(dat, idx, epoch, -2, 4))
end

bad_epochs_automatic = eegfun.detect_bad_epochs(epochs[1]);
bad_epochs_automatic = eegfun.detect_bad_epochs(epochs[1], abs_criterion = 200);
bad_epochs_automatic = eegfun.detect_bad_epochs(epochs[1], z_criterion = 3.0);
bad_epochs_automatic = eegfun.detect_bad_epochs(epochs[1], z_criterion = 3.0, abs_criterion = 200);

eegfun.unique_rejections(bad_epochs_automatic.rejected_epochs) # should just be the same
eegfun.unique_channels(bad_epochs_automatic.rejected_epochs)
eegfun.unique_epochs(bad_epochs_automatic.rejected_epochs)

bad_epochs_manual = eegfun.reject_epochs_interactive(epochs[1], dims = (4, 4))
bad_epochs_manual = eegfun.reject_epochs_interactive(epochs[1], dims = (4, 4), artifact_info = bad_epochs_automatic) 
bad_epochs_manual = eegfun.reject_epochs_interactive(epochs[1], dims = (4, 4), artifact_info = bad_epochs_automatic, colormap = :seaborn_colorblind) 
