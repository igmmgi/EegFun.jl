# Demo: Plotting RSA Results
# Shows how to visualise Representational Similarity Analysis results:
# RDM heatmaps, dissimilarity timecourses, and model correlations.

# Note: EegFun.example_path() resolves bundled example data paths.
# When using your own data, simply pass the file path directly, e.g.:
# dat = EegFun.read_raw_data("/path/to/your/data.bdf")

using EegFun
using GLMakie

#######################################################################
# LOAD DATA AND COMPUTE RSA
#######################################################################

dat = EegFun.read_raw_data(EegFun.example_path("data/bdf/example1.bdf"))
layout = EegFun.read_layout(EegFun.example_path("layouts/biosemi/biosemi72.csv"))
EegFun.polar_to_cartesian_xy!(layout)
dat = EegFun.create_eegfun_data(dat, layout)

EegFun.highpass_filter!(dat, 0.1)
EegFun.lowpass_filter!(dat, 30.0)

epoch_cfg = [
    EegFun.EpochCondition(name = "Condition1", trigger_sequences = [[1]]),
    EegFun.EpochCondition(name = "Condition2", trigger_sequences = [[2]]),
]
epochs = EegFun.extract_epochs(dat, epoch_cfg, (-0.2, 0.8))
EegFun.baseline!(epochs, (-0.2, 0.0))

# compute RSA
rsa_result = EegFun.rsa(epochs)

#######################################################################
# RDM HEATMAP — AVERAGE ACROSS TIME
#######################################################################

# visualise overall representational structure
EegFun.plot_rdm_heatmap(rsa_result)

#######################################################################
# RDM HEATMAP — SPECIFIC TIME POINT
#######################################################################

# RDM at 300 ms post-stimulus
EegFun.plot_rdm_heatmap(rsa_result, time_point = 0.3)

# RDM at time index 50
EegFun.plot_rdm_heatmap(rsa_result, time_point = 50)

#######################################################################
# DISSIMILARITY TIMECOURSE
#######################################################################

# plot dissimilarity over time for all condition pairs
EegFun.plot_rdm_timecourse(rsa_result)

# only specific pairs
EegFun.plot_rdm_timecourse(rsa_result, condition_pairs = [(1, 2)])

#######################################################################
# MODEL CORRELATIONS
#######################################################################

# compare RSA results against theoretical models
# (requires having run compare_models beforehand)
# EegFun.plot_model_correlations(rsa_result)
