# Demo: Filtering
# Shows how to apply high-pass and low-pass filters to EEG data,
# compare IIR and FIR methods, and visualize filter characteristics.

using EegFun

#######################################################################
# LOAD DATA
#######################################################################

dat = EegFun.read_raw_data("./resources/data/bdf/example1.bdf")
layout = EegFun.read_layout("./resources/layouts/biosemi/biosemi72.csv")
EegFun.polar_to_cartesian_xy!(layout)
dat = EegFun.create_test_eegfun_data(dat, layout)


#######################################################################
# HIGH-PASS FILTER
#######################################################################

# Remove slow drifts with a 0.1 Hz high-pass (standard for ERPs)
EegFun.highpass_filter!(dat, 0.1)

# Stronger high-pass for ICA preprocessing (1 Hz)
dat_ica = copy(dat)
EegFun.highpass_filter!(dat_ica, 1.0)


#######################################################################
# LOW-PASS FILTER
#######################################################################

# Remove high-frequency noise (30 Hz low-pass)
EegFun.lowpass_filter!(dat, 30.0)

# Anti-aliasing filter before downsampling (40 Hz)
EegFun.lowpass_filter!(dat, 40.0, order = 3)


#######################################################################
# FILTER SPECIFIC CHANNELS
#######################################################################

# Apply filter to only certain channels
EegFun.lowpass_filter!(dat, 30.0, channel_selection = EegFun.channels([:Cz, :Pz, :Fz]))


#######################################################################
# IIR vs FIR FILTERS
#######################################################################

# Default: IIR (Butterworth) — fast, good for most uses
EegFun.highpass_filter!(dat, 0.1, filter_method = "iir", order = 1)

# FIR (windowed sinc) — linear phase, better for time-domain analysis
EegFun.highpass_filter!(dat, 0.1, filter_method = "fir")


#######################################################################
# VISUALIZE FILTER RESPONSE
#######################################################################

# Create a filter object for inspection
filt = EegFun.create_highpass_filter(0.1, 2048.0)
EegFun.plot_filter_response(filt)

filt_lp = EegFun.create_lowpass_filter(30.0, 2048.0)
EegFun.plot_filter_response(filt_lp)

# Compare IIR vs FIR
filt_iir = EegFun.create_highpass_filter(0.1, 2048.0, filter_method = "iir")
filt_fir = EegFun.create_highpass_filter(0.1, 2048.0, filter_method = "fir")
EegFun.plot_filter_response(filt_iir)
EegFun.plot_filter_response(filt_fir)
