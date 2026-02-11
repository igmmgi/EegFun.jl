# Demo: Plotting Power Spectrum from SpectrumData
# Shows how to visualise frequency spectra computed by compute_spectrum,
# with options for channel selection, axis scaling, and dB units.

using EegFun

#######################################################################
# 1. LOAD DATA
#######################################################################

dat = EegFun.read_raw_data("./resources/data/bdf/example1.bdf")
layout = EegFun.read_layout("./resources/layouts/biosemi/biosemi72.csv")
EegFun.polar_to_cartesian_xy!(layout)
dat = EegFun.create_eegfun_data(dat, layout)

EegFun.highpass_filter!(dat, 0.5)

#######################################################################
# 2. COMPUTE SPECTRUM
#######################################################################

spectrum = EegFun.compute_spectrum(dat)

#######################################################################
# 3. BASIC SPECTRUM PLOT — ALL CHANNELS
#######################################################################

EegFun.plot_frequency_spectrum(spectrum)

#######################################################################
# 4. SINGLE CHANNEL
#######################################################################

EegFun.plot_frequency_spectrum(spectrum,
    channel_selection = channels(:Cz)
)

#######################################################################
# 5. LOG SCALE AXES
#######################################################################

EegFun.plot_frequency_spectrum(spectrum,
    channel_selection = channels(:Cz),
    x_scale = :log10,
    y_scale = :log10
)

#######################################################################
# 6. DECIBEL UNITS
#######################################################################

EegFun.plot_frequency_spectrum(spectrum,
    channel_selection = channels([:Cz, :Pz]),
    unit = :dB,
    max_freq = 100.0
)

#######################################################################
# 7. CUSTOM STYLING
#######################################################################

EegFun.plot_frequency_spectrum(spectrum,
    channel_selection = channels([:Cz, :Oz]),
    linewidth = 3,
    line_alpha = 0.6,
    title = "Power Spectrum"
)
