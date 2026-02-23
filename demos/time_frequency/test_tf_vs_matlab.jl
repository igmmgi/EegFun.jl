# Test: Load FieldTrip epoch data, run tf_morlet, and plot
# Compare visually with MATLAB's tfMyData output

using EegFun

# 1. Load FieldTrip .mat file
layout = EegFun.read_layout("./resources/layouts/biosemi/biosemi72.csv")
dat = EegFun.read_fieldtrip("/media/ian/SSD2/Projects/flankerConflict/ControlAttentionTask/EEG/Epoch1/3_1_pre.mat", layout)

# 2. TF morlet (matching MATLAB: freqs=1:1:40, cycles=5, mirror padding)
tf = EegFun.tf_morlet(dat; frequencies = 1:1:40, cycles = 5, pad = :both)

# 3. Plot — compare with MATLAB's 3_1_tf.mat
EegFun.plot_tf(
    tf;
    channel_selection = EegFun.channels(:CPz),
    baseline_interval = (-0.5, -0.3),
    baseline_method = :db,
    colorrange = (-3.0, 3.0),
    colormap = :jet,
    interpolate = false,
)
