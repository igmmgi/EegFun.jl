using eegfun
using GLMakie

# Load some test data and perform minimal preprocessing
data_file = joinpath(@__DIR__, "..", "..", "..", "Flank_C_3.bdf")
layout_file = eegfun.read_layout("./data/layouts/biosemi/biosemi72.csv")
eegfun.polar_to_cartesian_xy!(layout_file)
eegfun.polar_to_cartesian_xyz!(layout_file)
dat = eegfun.read_bdf(data_file)
dat = eegfun.create_eeg_dataframe(dat, layout_file)
eegfun.rereference!(dat, :avg)
eegfun.filter_data!(dat, "hp", 0.1)

# extreme value bool column
eegfun.is_extreme_value!(dat, 100)

# channel joint probability
channel_joint_probability = eegfun.channel_joint_probability(dat)
channel_joint_probability =
    eegfun.channel_joint_probability(dat, sample_selection = eegfun.samples_not(:is_extreme_value_100))


# Calculate EOG signals
eegfun.channel_difference!(
    dat,
    channel_selection1 = eegfun.channels([:Fp1, :Fp2]),
    channel_selection2 = eegfun.channels([:IO1, :IO2]),
    channel_out = :vEOG,
); # vertical EOG = mean(Fp1, Fp2) - mean(IO1, I02)
eegfun.channel_difference!(
    dat,
    channel_selection1 = eegfun.channels([:F9]),
    channel_selection2 = eegfun.channels([:F10]),
    channel_out = :hEOG,
); # vertical EOG = mean(Fp1, Fp2) - mean(IO1, I02)
eegfun.detect_eog_onsets!(dat, 50, :vEOG, :is_vEOG)
eegfun.detect_eog_onsets!(dat, 30, :hEOG, :is_hEOG)


# Calculate correlations between all channels and EOG channels
cm = eegfun.correlation_matrix_dual_selection(
    dat,
    sample_selection = eegfun.samples(),  # All samples
    channel_selection1 = eegfun.channels(),  # All EEG channels
    channel_selection2 = eegfun.channels([:vEOG, :hEOG]),  # EOG channels
)
eegfun.add_zscore_columns!(cm)

bad_channels = [:Fp1, :AF3]
non_eog_related, eog_related = eegfun.partition_channels_by_eog_correlation(bad_channels, cm)
